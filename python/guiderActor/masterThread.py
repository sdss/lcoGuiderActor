"""master thread for guiderActor."""
import Queue, threading
import math, numpy, re
import time
import subprocess
import pyfits
import os.path
import scipy.interpolate

from guiderActor import Msg, GuiderState, MASTER, GCAMERA
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
import opscore.utility.tback as tback
import opscore.utility.YPF as YPF
import RO

from gimg.guiderImage import GuiderImageAnalysis
from gimg import GuiderExceptions

def adiff(a1, a2):
    """ return a1-a2, all in degrees. """

    dd = a1-a2
    if dd >= 180.0:
        dd -= 360.0
    elif dd <= -180.0:
        dd += 360.0

    return dd

class FakeCommand(object):
    def _respond(self, tag, text):
        print "%s %s" % (tag, text)
    def warn(self, text):
        self._respond('w', text)
    def respond(self, text):
        self._respond('i', text)
    def inform(self, text):
        self._respond('i', text)
    def diag(self, text):
        self._respond('d', text)
    def finish(self, text):
        self._respond(':', text)
    def fail(self, text):
        self._respond('f', text)

def processOneFile(gState, guiderFile):
    queues = dict(MASTER=Queue.Queue())

    guideStep(None, queues, gState.cmd, guiderFile, True)

def processOneProcFile(gState, guiderFile, cartFile, plateFile, actor=None, queues=None, cmd=None, guideCmd=None):
    gState.setGuideMode('axes', False)
    gState.setGuideMode('focus', False)
    gState.setGuideMode('scale', False)

    if not cmd: cmd = FakeCommand()
    if not guideCmd: guideCmd = FakeCommand()
    if not queues: queues = dict(MASTER=Queue.Queue())

    gState.cmd = guideCmd 
    guideStep(None, queues, cmd, gState, guiderFile, True)

def send_decenter_status(cmd, gState, frameNo):
    """Output the decenter status keywords, including the most recent frame number."""
    cmd.inform("decenter=%d, %s, %7.3f, %7.3f, %7.3f, %7.3f, %7.3f" % (
                frameNo, ("enabled" if gState.decenter else "disabled"),
                gState.decenterRA, gState.decenterDec,
                gState.decenterRot, gState.decenterFocus, gState.decenterScale))
    cmd.inform("mangaDither=%s"%(gState.mangaDither))

def scale_pid_with_alt(cmd, gState, actorState):
    """Change the PID coefficients with altitude, to deal with high-alt tracking."""
    alt = actorState.models["tcc"].keyVarDict['axePos'][1]
    if gState.scale_pid_with_alt(alt):
        gState.output_pid()

def _check_fiber(fiber, gState, guideCmd):
    """Check whether the current fiber should currently be enabled."""
    # necessary?
    if fiber.gProbe is None:
        guideCmd.warn('text="Gprobe %d was not listed in plugmap info"' % fiber.fiberid)
        return False

    # Center up on acquisition fibers only.
    if gState.centerUp and fiber.gProbe.fiber_type != "ACQUIRE":
        guideCmd.diag('text="Gprobe %d is disabled during Center Up."' % fiber.fiberid)
        return False
    else:
        if fiber.gProbe.disabled:
            guideCmd.diag('text="Gprobe %d is not enabled."' % fiber.fiberid)
        return fiber.gProbe.enabled
#...

def _do_one_fiber(fiber, gState, cmd, frameInfo, haLimWarn):
    """
    Process one single fiber, computing various scales and corrections.
    """
    gProbe = fiber.gProbe
    
    # dx, dy are the offsets on the ALTA guider image
    fiber.dx = frameInfo.guideCameraScale*(fiber.xs - fiber.xcen) + (gProbe.xFerruleOffset / 1000.)
    fiber.dy = frameInfo.guideCameraScale*(fiber.ys - fiber.ycen) + (gProbe.yFerruleOffset / 1000.)
    poserr = fiber.xyserr

    # theta is the angle to rotate (x, y) on the ALTA to (ra, alt)
    # phi is the orientation of the alignment hole measured clockwise from N
    # rotation is the anticlockwise rotation from x on the ALTA to the pin
    theta = 90                   # allow for 90 deg rot of camera view, should be -90
    theta += gProbe.rotation # allow for intrinsic fibre rotation
    try:
        theta -= gProbe.phi      # allow for orientation of alignment hole
    except Exception as e:
        cmd.warn('text="skipping phi-less probe %s"' % (fiber.fiberid))
        return
    
    gProbe.rotStar2Sky = theta # Squirrel the real angle away.

    #FIXME PH -- We should ignore gprobes not present on plate/pointing (MARVELS dual pointing)
    #               and ignore fibers not found in flat.
    #            However we probably want to record values of disabled fibers for diagnosis
    if numpy.isnan(fiber.dx) or numpy.isnan(fiber.dy) or numpy.isnan(poserr):
        cmd.warn("text=%s" %
                      qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber measured=(%g, %g), nominal=(%g,%g)" % (
                          fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, gProbe.xCenter, gProbe.yCenter)))
        return

    if fiber.flux < frameInfo.minStarFlux:
        cmd.warn("text=%s" %
                      qstr("Star in gprobe %d too faint for guiding flux %g < %g minimum flux" % (
                          fiber.fiberid, fiber.flux, frameInfo.minStarFlux)))
        gProbe.tooFaint = True
    else:
        gProbe.tooFaint = False

    if poserr == 0:
        cmd.warn("text=%s" %
                      qstr("position error is 0 for gprobe %d star=(%g, %g) fiber=(%g, %g) nominal=(%g,%g)" % (
                          fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, gProbe.xCenter, gProbe.yCenter)))
        return

    theta = math.radians(theta)
    ct, st = math.cos(theta), math.sin(theta)
    # error in guide star position; n.b. still in mm here
    dRA   =  fiber.dx*ct + fiber.dy*st
    dDec  = -fiber.dx*st + fiber.dy*ct
    dDec *= -1

    #FIXME PH -- calc dAlt and dAz for guiding diagnostics,(output as part of fiber?)
    
    # Apply refraction correction
    xRefractCorr = 0.0
    yRefractCorr = 0.0
    haTime = 0.0
    try:
        if gState.refractionBalance > 0:
            if frameInfo.wavelength in gProbe.haOffsetTimes:
                haTimes = gProbe.haOffsetTimes[frameInfo.wavelength]
                if frameInfo.dHA < haTimes[0]:
                    if not haLimWarn:
                        cmd.warn('text="dHA (%0.1f) is below interpolation table; using limit (%0.1f)"' % (frameInfo.dHA, haTimes[0]))
                        haLimWarn = True
                    haTime = haTimes[0]
                elif frameInfo.dHA > haTimes[-1]:
                    if not haLimWarn:
                        cmd.warn('text="dHA (%0.1f) is above interpolation table; using limit (%0.1f)"' % (frameInfo.dHA, haTimes[-1]))
                        haLimWarn = True
                    haTime = haTimes[-1]
                else:
                    haTime = frameInfo.dHA
    
                # I'm now assuming 0...offset, but it should be offset1...offset2
                xInterp = scipy.interpolate.interp1d(haTimes,
                                                     gProbe.haXOffsets[frameInfo.wavelength])
                xRefractCorr = gState.refractionBalance * xInterp(haTime)
                yInterp = scipy.interpolate.interp1d(haTimes,
                                                     gProbe.haYOffsets[frameInfo.wavelength])
                yRefractCorr = gState.refractionBalance * yInterp(haTime)
            else:
                # JKP: TODO: these warnings might be excessive?
                cmd.warn('text="No HA Offset Time available for probe %d at wavelength %d. No refraction offset calculated."'%(gProbe.id,frameInfo.wavelength))
        else:
            # Don't do anything if the refraction balance is 0.
            pass
    except Exception as e:
        cmd.diag('text="failed to calc refraction offsets for %s: %s"' % (frameInfo.wavelength, e))
        pass

    cmd.inform('refractionOffset=%d,%d,%0.1f,%0.4f,%0.6f,%0.6f' % (frameInfo.frameNo, fiber.fiberid,
                                                                        gState.refractionBalance,
                                                                        haTime,
                                                                        xRefractCorr*frameInfo.arcsecPerMM,
                                                                        yRefractCorr*frameInfo.arcsecPerMM))
    dRA -= xRefractCorr
    dDec -= yRefractCorr
    
    # Apply RA & Dec user guiding offsets to mimic different xy fibers centers
    # The guiderRMS will be calculated around the new effective fiber centers
    if gState.decenter:
        # apply decenter offset so that telescope moves (not the star)
        dRA  += gState.decenterRA/frameInfo.arcsecPerMM
        dDec += gState.decenterDec/frameInfo.arcsecPerMM
        #decenterRot applied after guide solution
    
    fiber.dRA = dRA
    fiber.dDec = dDec
    raCenter  = gProbe.xFocal
    decCenter = gProbe.yFocal
    
    cmd.inform("probe=%d,%2d,0x%02x, %7.2f,%7.2f, %7.3f,%4.0f, %7.2f,%6.2f,%6.2f, %7.2f,%6.2f" % (
        frameInfo.frameNo, fiber.fiberid, gProbe.gprobebits,
        fiber.dRA*frameInfo.arcsecPerMM, fiber.dDec*frameInfo.arcsecPerMM,
        fiber.fwhm, gProbe.focusOffset,
        fiber.flux, fiber.mag, gProbe.ref_mag, fiber.sky, fiber.skymag))
    
    if gProbe.tooFaint:
        return
    
    #Collect fwhms for good in focus stars
    if gProbe.atFocus and gProbe.good:
        frameInfo.inFocusFwhm.append(fiber.fwhm)
    
    #accumulate guiding errors for good stars used in fit
    frameInfo.guideRMS += fiber.dx**2 + fiber.dy**2
    frameInfo.guideXRMS += fiber.dx**2
    frameInfo.guideYRMS += fiber.dy**2
    frameInfo.nguideRMS += 1
    frameInfo.guideRaRMS += dRA**2
    frameInfo.guideDecRMS += dDec**2
    #guideAzRMS += fiber.dAz**2
    #guideAltRMS += fiber.dAlt**2

    frameInfo.b[0] += dRA
    frameInfo.b[1] += dDec
    frameInfo.b[2] += raCenter*dDec - decCenter*dRA

    frameInfo.A[0, 0] += 1
    frameInfo.A[0, 1] += 0
    frameInfo.A[0, 2] += -decCenter

    frameInfo.A[1, 0] += 0
    frameInfo.A[1, 1] += 1
    frameInfo.A[1, 2] += raCenter

    frameInfo.A[2, 2] += raCenter*raCenter + decCenter*decCenter

    # Now scale.  We don't actually solve for scale and axis updates
    # simultanously, and we don't allow for the axis update when
    # estimating the scale.
    frameInfo.b3 += raCenter*dRA + decCenter*dDec
#...

def _find_focus_one_fiber(fiber,gState,frameInfo,C,A,b):
    """Accumulate the focus for one fiber into A and b."""
    # required?
    if fiber.gProbe is None:
        return
    gProbe = gState.gprobes[fiber.fiberid]
    if not gProbe.enabled:
        return

    # FIXME -- do we want to include ACQUISITION fibers?
    # PH -- currently all valid enabled fibers are used so OK.
    rms = fiber.fwhm / frameInfo.sigmaToFWHM
    if numpy.isnan(rms):
        return

    rms *= frameInfo.micronsPerArcsec # in microns
    rmsErr = 1

    d = gProbe.focusOffset
    x = rms*rms - C*d*d
    xErr = 2*rms*rmsErr

    try:
        ivar = 1/(xErr*xErr)
    except ZeroDivisionError:
        ivar = 1e-5

    b[0] += x*ivar
    b[1] += x*d*ivar

    A[0, 0] += ivar
    A[0, 1] += d*ivar

    A[1, 1] += d*d*ivar
#...

def apply_radecrot(cmd, gState, actor, actorState, offsetRa, offsetDec, offsetRot):
    """Finish the calculation for the ra/dec/rot corrections and apply them."""
    if gState.guideAxes:
        cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                 cmdStr="offset arc %f, %f"%(-offsetRa, -offsetDec))
        if cmdVar.didFail:
            cmd.warn('text="Failed to issue offset"')

        if offsetRot:
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                     cmdStr="offset guide %f, %f, %g"%(0.0, 0.0, -offsetRot))
        if cmdVar.didFail:
            cmd.warn('text="Failed to issue offset in rotator"')

    elif gState.centerUp:
        # If we are in the middle of an fk5InFiber (or other TCC track/pterr),
        # adjust the calibration offsets
        doCalibOffset = actorState.models["tcc"].keyVarDict["objName"][0] == "position reference star"
        if doCalibOffset:
            cmd.warn('text="using arc offsets at a pointing star"')
            
        cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                 cmdStr="offset arc %f, %f" % \
                                 (-offsetRa, -offsetDec))
        if cmdVar.didFail:
            cmd.warn('text="Failed to issue centering offset"')
        else:
            if not doCalibOffset:
                gState.setGuideMode('axes', True)


def guideStep(actor, queues, cmd, gState, inFile, oneExposure,
              guiderImageAnalysis, output_verify='warn', camera='gcamera'):
    """ One step of the guide loop, based on the given guider file.

    Args:
        actor: the current actor instance: lets us send direct tcc commands.
        queues: the queue list, so we can put commands on it.
        cmd: the currently active command, for message passing.
        gState: an instance of GuiderState, holding information about the gprobes, etc.
        inFile: the name of the current raw gcamera file.
        oneExposure: True if we are only handling a single exposure.
        guiderImageAnalysis: an instance of that class, to process the raw image.
        output_verify: passed on to the fits writer. See the pyfits docs for more.
        camera: set to 'ecamera' to not search for fibers and skip fiber-related processing.
    """
    # Setup to solve for the axis and maybe scale offsets.  We work consistently
    # in mm on the focal plane, only converting to angles to command the TCC.
    guideCameraScale = gState.gcameraMagnification*gState.gcameraPixelSize*1e-3 # mm/pixel
    arcsecPerMM = 3600./gState.plugPlateScale   #arcsec per mm
    frameNo = int(re.search(r"([0-9]+)\.fits.*$", inFile).group(1))

    # Object to gather all per-frame guiding info into.
    frameInfo = GuiderState.FrameInfo(frameNo,arcsecPerMM,guideCameraScale,gState.plugPlateScale)

    actorState = guiderActor.myGlobals.actorState
    guideCmd = gState.cmd
    guideCmd.respond("processing=%s" % inFile)
    
    h = pyfits.getheader(inFile)
    flatfile = h.get('FLATFILE', None)
    flatcart = h.get('FLATCART', None)
    darkfile = h.get('DARKFILE', None)
    if not flatfile:
        guideCmd.fail('guideState="failed"; text=%s' % qstr("No flat image available"))
        gState.cmd = None
        return frameInfo
    
    if not darkfile:
        guideCmd.fail('guideState="failed"; text=%s' % qstr("No dark image available"))
        gState.cmd = None
        return frameInfo

    if flatcart != gState.cartridge:
        if False:
            guideCmd.fail('guideState="failed"; text=%s' % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                            flatcart, gState.cartridge)))
            gState.cmd = None
            return frameInfo
        else:
            guideCmd.warn("text=%s" % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                flatcart, gState.cartridge)))

    try:
        setPoint = actorState.models[camera].keyVarDict["cooler"][0]
        guideCmd.inform('text="guideStep GuiderImageAnalysis.findStars()..."')
        fibers = guiderImageAnalysis(cmd,inFile,gState.gprobes,setPoint=setPoint,bypassDark=actorState.bypassDark,camera=camera)
        guideCmd.inform("text='GuiderImageAnalysis.findStars() got %i fibers'" % len(fibers))
    except GuiderExceptions.BadReadError as e:
        guideCmd.warn('text=%s' %qstr("Skipping badly formatted image."))
        return frameInfo
    except GuiderExceptions.FlatError as e:
        guideCmd.fail('guideState="failed"; text=%s' %qstr("Error reading/processing %s flat: %s"%(camera,e)))
        gState.cmd = None
        return frameInfo
    except GuiderExceptions.GuiderError as e:
        guideCmd.fail('guideState="failed"; text=%s' %qstr("Error processing %s image: %s"%(camera,e)))
        gState.cmd = None
        return frameInfo
    except Exception as e:
        guideCmd.fail('guideState="failed"; text=%s' % qstr("Unknown error in processing guide images: %s" % e))
        gState.cmd = None
        tback.tback("GuideTest", e)
        return frameInfo

    # Don't need to do anything else with ecam images.
    if camera == 'ecamera':
        # No more processing needed, so just write the file.
        guiderImageAnalysis.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes, output_verify=output_verify)
        return frameInfo

    #
    # N.B. fiber.xFocal and fiber.yFocal are the offsets of the stars
    # wrt the center of the plate in mm; fiber.xcen/star.xs are in pixels,
    # so we need a scale for the guide camera.  Nominally the guide camera
    # has the same scale as the plug plate itself, but maybe it doesn't,
    # so we'll include a possible magnification
    #

    # Grab some times for refraction correction
    longitude = -105.82045
    UTC = RO.Astro.Tm.utcFromPySec(time.time() +
                                   actorState.models["tcc"].keyVarDict["utc_TAI"][0])
    LST = RO.Astro.Tm.lastFromUT1(UTC, longitude)
    
    try:
        # NEWTCC: is there a better keyword than this to get the current RA?
        RAkey = actorState.models["tcc"].keyVarDict["objNetPos"][0]
        RA = RAkey.getPos()
        HA = adiff(LST, RA)     # The corrections are indexed by degrees, happily.
        frameInfo.dHA = adiff(HA, gState.design_ha)
    except:
        guideCmd.error('text="Could not determine current RA from TCC objNetPos. Please issue: tcc show object /full"')
        guideCmd.warn('text="WARNING: refraction corrections to guiding will not work until this is dealt with."')
        RA = numpy.nan
        HA = numpy.nan
        frameInfo.dHA = 0
    guideCmd.diag('text="LST=%0.4f RA=%0.4f HA=%0.4f desHA=%0.4f dHA=%0.4f"' %
                  (LST, RA, HA, gState.design_ha,frameInfo.dHA))

    # scale the PID values to function better at high alt.
    scale_pid_with_alt(guideCmd, gState, actorState)

    # Set the decenter parameters in frameInfo, so they don't get lost if decentering
    # changes during processing.
    if gState.decenter:
        frameInfo.setDecenter(gState)
    else:
        frameInfo.setDecenter(None)

    # At present, only APOGEE uses refractionOffsets.
    # So, only this wavelength will have haOffsetTime specified for each gprobe.
    frameInfo.wavelength = 16600
    frameInfo.refractionBalance = gState.refractionBalance

    haLimWarn = False # so we only warn once about passing the HA limit for refraction balance
    for fiber in fibers:
        if _check_fiber(fiber, gState, guideCmd):
            _do_one_fiber(fiber, gState, guideCmd, frameInfo, haLimWarn)

    nStar = frameInfo.A[0, 0]
    if nStar == 0 or gState.inMotion:
        if nStar == 0:
            guideCmd.warn('text="No stars are available for guiding."')
        else:
            guideCmd.warn('text="Telescope moved during exposure -- skipping this image."')

        guiderImageAnalysis.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes, output_verify=output_verify)

        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.cmd = None
            return frameInfo

        #if guidingIsOK(cmd, actorState):
        #    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
        #                            expTime=gState.expTime))
        return frameInfo
        
    frameInfo.A[2, 0] = frameInfo.A[0, 2]
    frameInfo.A[2, 1] = frameInfo.A[1, 2]
    try:
        if nStar == 1:
            guideCmd.warn('text="Only one star is usable"')
            x = frameInfo.b
            x[2, 0] = 0 # no rotation
        else:
            x = numpy.linalg.solve(frameInfo.A, frameInfo.b)

        # convert from mm to degrees
        dRA = x[0, 0]/gState.plugPlateScale
        dDec = x[1, 0]/gState.plugPlateScale
        dRot = -math.degrees(x[2, 0]) # and from radians to degrees
        
        # TBD: we are not applying any rotation decentering at present.
        #PH Kludge add the decenter guiding rotation offset here for now (in degrees)
        #if gState.decenter:
        #    dRot += gState.decenterRot/3600.0
        
        frameInfo.dRA  = dRA
        frameInfo.dDec = dDec
        frameInfo.dRot = dRot
        
        # directly apply a shift for centerUp and decentering.
        # otherwise, apply the shift via the usual pid.
        dt = gState.update_pid_time('raDec', time.time())
        gState.update_pid_time('rot', time.time())
        if gState.centerUp or gState.decenterCmd:
            offsetRa = -dRA
            offsetDec = -dDec
            offsetRot = 0
        else:
            offsetRa  = -gState.pid["raDec"].update(dRA, dt=dt)
            offsetDec = -gState.pid["raDec"].update(dDec, dt=dt)
            offsetRot = -gState.pid["rot"].update(dRot, dt=dt) if nStar > 1 else 0 # don't update I

        frameInfo.filtRA  = offsetRa
        frameInfo.filtDec = offsetDec
        frameInfo.filtRot = offsetRot

        frameInfo.offsetRA  = offsetRa  if (gState.guideAxes or gState.centerUp) else 0.0
        frameInfo.offsetDec = offsetDec if (gState.guideAxes or gState.centerUp) else 0.0
        frameInfo.offsetRot = offsetRot if (gState.guideAxes or gState.centerUp) else 0.0

        guideCmd.respond("axisError=%g, %g, %g" % (3600*dRA, 3600*dDec, 3600*dRot))
        guideCmd.respond("axisChange=%g, %g, %g, %s" % (-3600*offsetRa, -3600*offsetDec, -3600*offsetRot,
                                                        "enabled" if gState.guideAxes else "disabled"))
        #calc FWHM with trimmed mean for 8 in focus fibers
        nFwhm = len(frameInfo.inFocusFwhm)
        trimLo = 1 if nFwhm > 4 else 0
        trimHi= nFwhm - trimLo
        nKept = nFwhm - 2*trimLo
        nReject = nFwhm - nKept
        meanFwhm = (sum(frameInfo.inFocusFwhm))/nFwhm if nFwhm>0 else numpy.nan
        tMeanFwhm = (sum(sorted(frameInfo.inFocusFwhm)[trimLo:trimHi]))/nKept if nKept>0 else numpy.nan
        #loKept = frameInfo.inFocusFwhm[trimLo]
        #hiKept = frameInfo.inFocusFwhm[(trimHi-1)]
        infoString = "fwhm=%d, %7.2f, %d, %d, %7.2f" % (frameNo, tMeanFwhm, nKept, nReject, meanFwhm)
        guideCmd.inform(infoString)
        
        frameInfo.meanFwhm = meanFwhm
        frameInfo.tMeanFwhm = tMeanFwhm

        #rms position error prior to this frame's correction
        try:
            frameInfo.guideRMS  = math.sqrt(frameInfo.guideRMS/frameInfo.nguideRMS) *arcsecPerMM
            frameInfo.guideXRMS = math.sqrt(frameInfo.guideXRMS/frameInfo.nguideRMS) *arcsecPerMM
            frameInfo.guideYRMS = math.sqrt(frameInfo.guideYRMS/frameInfo.nguideRMS) *arcsecPerMM
            frameInfo.guideRaRMS = math.sqrt(frameInfo.guideRaRMS/frameInfo.nguideRMS) *arcsecPerMM
            frameInfo.guideDecRMS = math.sqrt(frameInfo.guideDecRMS/frameInfo.nguideRMS) *arcsecPerMM
        except:
            frameInfo.guideRMS = numpy.nan
            frameInfo.guideXRMS = numpy.nan
            frameInfo.guideYRMS = numpy.nan
            frameInfo.guideRaRMS = numpy.nan
            frameInfo.guideDecRMS = numpy.nan

        #FIXME PH ---Need to calculate Az and Alt RMS in arcsec
        # guideAzRMS  = numpy.nan
        # guideAltRMS = numpy.nan
        #frameInfo.guideAzRMS = guideAzRMS
        #frameInfo.guideAltRMS = guideAltRMS

        apply_radecrot(cmd, gState, actor, actorState, offsetRa, offsetDec, offsetRot)

    except numpy.linalg.LinAlgError:
        guideCmd.warn("text=%s" % qstr("Unable to solve for axis offsets"))

    # don't bother with focus/scale!
    if nStar <= 1 or gState.centerUp:
        guiderImageAnalysis.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes, output_verify=output_verify)
        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.cmd = None
        return frameInfo

    #
    # Scale
    #
    dScale = frameInfo.b3/frameInfo.A[2, 2]
    dt = gState.update_pid_time('scale', time.time())
    offsetScale = -gState.pid["scale"].update(dScale, dt=dt)

    frameInfo.dScale = dScale
    frameInfo.filtScale = offsetScale
    frameInfo.offsetScale = offsetScale if gState.guideScale else 0.0

    guideCmd.respond("scaleError=%g" % (dScale))
    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                             "enabled" if gState.guideScale else "disabled"))
    # the below is used by the observers to track the scale deltas.
    guideCmd.inform('text="delta percentage scale correction = %g"' % (-dScale*100.))
    
    # There is (not terribly surprisingly) evidence of crosstalk between scale and focus adjustements.
    # So for now defer focus changes if we apply a scale change.
    blockFocusMove = False
        
    if gState.guideScale:
        # This should be a tiny bit bigger than one full M1 axial step.
        if abs(offsetScale) < 3.4e-7:
            cmd.diag('text="skipping small scale change=%0.8f"' % (offsetScale))
        else:
            # Clip to the motion we think is too big to apply at once.
            offsetScale = 1 + max(min(offsetScale, 2e-6), -2e-6)

            # Last chance to bailout.
            if offsetScale < 0.9995 or offsetScale > 1.0005:
                cmd.warn('text="NOT setting scarily large scale=%0.8f"' % (offsetScale))
            else:
                # blockFocusMove = True
                cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                         cmdStr="set scale=%.9f /mult" % (offsetScale))
                if cmdVar.didFail:
                    guideCmd.warn('text="Failed to issue scale change"')

    #Evaluate RMS on fit over fibers used in fits here
    #FIXME--PH not calculated yet
    guideFitRMS = numpy.nan
    nguideFitRMS = 0
    nguideRejectFitRMS = 0

    # RMS guiding error output has to be after scale estimation so the full fit residual can be reported
    guideCmd.inform("guideRMS=%5d,%4.3f,%4d,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4d,%4d,%4.3f,%4.3f" % (
        frameInfo.frameNo, frameInfo.guideRMS, frameInfo.nguideRMS,
        frameInfo.guideAzRMS, frameInfo.guideAltRMS,
        frameInfo.guideXRMS, frameInfo.guideYRMS, guideFitRMS, nguideFitRMS, nguideRejectFitRMS,
        frameInfo.guideRaRMS, frameInfo.guideDecRMS))

    #
    # Now focus. If the ith star is d_i out of focus, and the RMS of an
    # in-focus star would be r0, and we are Delta out of focus, we measure
    # an RMS size R_i
    #   R_i^2 = r0^2 + C (d_i + Delta)^2
    # i.e.
    #   R_i^2 - C d_i^2 = (r0^2 + C Delta^2) + 2 C Delta d_i
    # which is a linear equation for x == R_i^2 - C d_i^2
    #
    # If the secondary is half the diameter of the primary, the
    # RMS^2 size of an image of radius r is 5/8 r^2.  The f ratio
    # is f, so if the image is formed a distance d from focus, the
    # radius of the doughnut is d/(2 f) so
    # RMS^2 = 5/(32 f^2) d^2, i.e. C = 5/(32 f^2)
    #
    focalRatio = 5.0
    C = 5/(32.0*focalRatio*focalRatio)

    A = numpy.matrix(numpy.zeros(2*2).reshape([2,2]))
    b = numpy.matrix(numpy.zeros(2).reshape([2,1]))

    for fiber in fibers:
        _find_focus_one_fiber(fiber,gState,frameInfo,C,A,b)

    A[1, 0] = A[0, 1]
    try:
        x = numpy.linalg.solve(A, b)

        Delta = x[1, 0]/(2*C)
        try:
            rms0 = math.sqrt(x[0, 0] - C*Delta*Delta)/frameInfo.micronsPerArcsec
        except ValueError as e:
            rms0 = float("NaN")

        # Note sign change here.
        dFocus = -Delta*gState.dSecondary_dmm # mm to move the secondary
        dt = gState.update_pid_time('focus', time.time())
        offsetFocus = -gState.pid["focus"].update(dFocus, dt=dt)

        frameInfo.dFocus = dFocus
        frameInfo.filtFocus = offsetFocus
        frameInfo.offsetFocus = offsetFocus if gState.guideFocus else 0.0
        frameInfo.seeing = rms0*frameInfo.sigmaToFWHM   #in arc sec

        guideCmd.respond("seeing=%g" % (rms0*frameInfo.sigmaToFWHM))
        guideCmd.respond("focusError=%g" % (dFocus))
        guideCmd.respond("focusChange=%g, %s" % (offsetFocus, "enabled" if (gState.guideFocus and not blockFocusMove) else "disabled"))
        if gState.guideFocus and not blockFocusMove:
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="set focus=%f/incremental" % (offsetFocus),
                                     timeLim=20)

            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue focus offset"')
    except numpy.linalg.LinAlgError:
        guideCmd.respond("focusError=%g" % (numpy.nan))
        guideCmd.respond("focusChange=%g, %s" % (numpy.nan, "enabled" if (gState.guideFocus and not blockFocusMove) else "disabled"))
        guideCmd.warn("text=%s" % qstr("Unable to solve for focus offset"))
        x = None

    # Write output fits file for TUI
    guiderImageAnalysis.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes, output_verify=output_verify)
    return frameInfo
#...

def loadAllProbes(cmd, gState):
    """
    Read in information about the current guide probes from the platedb.
    
    The contents of the plPlugMap table are documented in the data model here:
    http://data.sdss3.org/datamodel/files/PLATELIST_DIR/runs/PLATERUN/plPlugMap.html
    """
    gState.allProbes = None
    try:
        path = 'catPlPlugMapM'
        cmd1 = "%s -c %s -m %s -p %s -f %s %s" % (path,
                                                  gState.cartridge, gState.fscanMJD,
                                                  gState.pointing, gState.fscanID,
                                                  gState.plate)
        try:
            cmd.diag('text=%s' % (qstr('running: %s' % (cmd1))))
            ret = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
            plugmapBlob, errText = ret.communicate()
        except subprocess.CalledProcessError as e:
            cmd.warn('text="failed to load plugmap file: %s"' % (e))
            return

        ypm = YPF.YPF(fromString=plugmapBlob)
        pm = ypm.structs['PLUGMAPOBJ'].asArray()
        
        # output information about the science program for this plate.
        # jkp TBD: may need to conver these to strings in some way
        # to catch ["marvels","apogee"], but I need to see how lists are
        # handled in the YPF, which requires an example...
        #instruments = ypm['instruments'].value
        #platetype = ypm['platetype'].value
        #cmd.info('scienceProgram=%s,%s'%(instruments,platetype))
        
        # It is useful to keep the object information as well,
        # so that we can put "any star down any hole". This is potentially
        # very useful for testing.
        # TBD: we'll probably need a new type here for MaNGA.
        keep = pm[numpy.where(((pm.holeType == "GUIDE") & (pm.objType == "NA"))
                              | (pm.holeType == "OBJECT"))]
        cmd.diag('text="kept %d probes"' % (len(keep)))
        gState.allProbes = keep
    except Exception as e:
        cmd.warn('text=%s' % (qstr("could not load all probe info: %s" % (e))))
    
def loadTccBlock(cmd, actorState, gState):
    """
    This is used for fk5InFiber (exclusively, I think).
    We should only call it immediately prior to using fk5InFiber, via:
       guider prepFk5InFiber

    !!!!!!!!!!!!!!!!!!11
    jkp TBD: Why do we do this separately from loadAllProbes?
    They issue the same catPlPlugMapM command, right?
    I think we can merge them together, and/or get rid of loadAllProbes and
    replace it with specific platedbActor calls to load just the bits we want.
    !!!!!!!!!!!!!!!!!!11
    """
    try:
        cmd1 = "catPlPlugMapM -c %s -m %s -p %s -f %s %s" % (gState.cartridge, gState.fscanMJD,
                                                             gState.pointing, gState.fscanID,
                                                             gState.plate)
        #plate10k
        scratchFile = '/tmp/v_ca1_%04d.dat' % (gState.plate)
        cmd1 += " | %s /dev/stdin > %s" % (os.path.join(os.environ['GUIDERACTOR_DIR'],
                                                        'bin',
                                                        'convertPlPlugMap.py'),
                                           scratchFile)
        cmd.diag('text=%s' % (qstr('running: %s' % (cmd1))))
        
        cmd2 = " %s %s" % (os.path.join(os.environ['GUIDERACTOR_DIR'],
                                        'bin',
                                        'xferBlock.py'),
                           scratchFile)
        cmd.diag('text=%s' % (qstr('running: %s' % (cmd2))))
        ret = subprocess.call(cmd1, shell=True)
        if ret < 0:
            raise RuntimeError("cat and convert job failed with %s" % (-ret))
        ret = subprocess.call(cmd2, shell=True)
        if ret < 0:
            raise RuntimeError("xfer job failed with %s" % (-ret))

    except Exception as e:
        cmd.warn('text=%s' % (qstr("could not load a per-cartridge instrument block: %s" % (e))))
#...

def make_movie(actorState, cmd, start):
    """Make a movie from guider frames, from start to the most recent."""
    if start is None or start <= 0:
        cmd.diag("text='No start frame defined for make_movie.'")
        # noone bothered to define the start, so we don't have anything to work with.
        return False
    
    # Ask gcamera for the next frame, which is the end frame+1
    endFrame = actorState.models['gcamera'].keyVarDict['nextSeqno'][0]
    # use the simulator frame number, if it's been configured.
    simulating = actorState.models['gcamera'].keyVarDict['simulating']
    if simulating[0]:
        endFrame = simulating[2]
    endFrame -= 1
    filename = actorState.models['guider'].keyVarDict['file'][1]
    endSearch = re.search(r"([0-9]+)\.fits*", filename)
    if not endSearch:
        return False
    else:
        end = int(endSearch.group(1))
    # Don't bother making the movie if we haven't been operating the guider for very long.
    if end - start < 5:
        cmd.diag("text='Too few exposures to bother making a movie out of.'")
        return False
    
    # callCommand produces a movie command that is not connected to the current command.
    # This will prevent "This command has already finished" complaints.
    actorState.actor.callCommand("makeMovie start=%d end=%d"%(start,end))
    return True
#...

def cal_finished(msg, name, guiderImageAnalysis, actorState, gState):
    """Generic handling of finished dark/flat frame."""
    cmd = msg.cmd
    cmd.respond("processing=%s" % msg.filename)
    # need ".*" in the regex, because we may or may not have gzipped files.
    # frameNo = int(re.search(r"([0-9]+)\.fits.*", msg.filename).group(1))
    
    header = pyfits.getheader(msg.filename)
    exptype = header.get('IMAGETYP')
    if exptype != name:
        cmd.fail('text="%s image processing ignoring a %s image!!"' % (name,exptype))
        return
        
    cmd.diag('text="cal_finished guiderImageAnalysis.analyze%s()..."'%name)
    try:
        func = ''.join(("analyze",name[0].upper(),name[1:]))
        # Always read the setPoint, so that it is as up-to-date as possible.
        camera = 'ecamera' if gState.plateType == 'ecamera' else 'gcamera'
        guiderImageAnalysis.camera = camera
        setPoint = actorState.models[camera].keyVarDict["cooler"][0]
        if name == 'flat':
            guiderImageAnalysis.analyzeFlat(msg.filename,gState.gprobes,cmd,setPoint)
        elif name == 'dark':
            guiderImageAnalysis.analyzeDark(msg.filename,cmd,setPoint)
        else:
            raise ValueError("Don't know how to finish a %s guider cal."%name)
    except GuiderExceptions.GuiderError as e:
        cmd.error("text=%s"%qstr(e))
        cmd.fail('guideState="failed"; text=%s' %qstr("%s failed. Error reading/processing guider %s."%(func,name)))
        gState.cmd = None
        return
    except Exception as e:
        tback.tback("cal_finished", e)
        cmd.fail('text="%s failed for an unknown reason: %s' % (func,e))
        return
    
    try:
        outname = guiderImageAnalysis.getProcessedOutputName(msg.filename)
        dirname, filename = os.path.split(outname)
        cmd.inform('file=%s/,%s' % (dirname, filename))
        cmd.finish('text="%s image processing done"'%name)
    except Exception as e:
        tback.tback("cal_finished", e)
        cmd.fail('text="failed to save flat: %s"' % (e))
#...

def dark_finished(msg, guiderImageAnalysis, actorState, gState):
    """Process a finished dark frame."""
    cal_finished(msg,'dark',guiderImageAnalysis,actorState,gState)

def flat_finished(msg, guiderImageAnalysis, actorState, gState):
    """Process a finished flat frame."""
    header = pyfits.getheader(msg.filename)
    darkfile = header.get('DARKFILE', None)
    if not darkfile:
        msg.cmd.fail("text=%s" % qstr("No dark image listed in flat header!!"))
        return
    cal_finished(msg,'flat',guiderImageAnalysis,actorState,gState)
#...

#
# Actual guider commands, and sub-commands.
#
def load_cartridge(msg, queues, gState, actorState):
    """
    Load cartridge information into the appropriate systems.
    1. Process the information from a loadCartridge command into gState.
    2. load the instrument block into the TCC.
    3. load all guide probes, for "any star->any fiber" tricks.
    """
    
    gState.deleteAllGprobes()

    gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
    gState.fscanMJD, gState.fscanID = msg.fscanMJD, msg.fscanID
    gState.boresight_ra, gState.boresight_dec = msg.boresight_ra, msg.boresight_dec
    gState.design_ha = msg.design_ha
    gState.plateType = msg.survey
    gState.surveyMode = msg.surveyMode
    for id, gProbe in msg.gprobes.items():
        gState.gprobes[id] = gProbe
            
    # Build and install an instrument block for this cartridge info
    # NOTE: TBD: we don't actually need to do loadTccBlock unless we want fk5infiber.
    # See ticket #2229 for how we should deal with this: the blocks end up litering
    # a directory on the tcc, and we can't upload them to the new tcc via scp anyway.
    # loadTccBlock(msg.cmd, actorState, gState)

    # newtcc NOTE: We still need to "set inst=spectro", but the rest of it we probably don't need.
    #cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=msg.cmd,
    #                                    cmdStr="set inst=spectro")#/gcview=%s/keep=(scaleFac)" % (gState.plate))
    #if cmdVar.didFail:
    #    msg.cmd.fail('text="Failed to set inst!"')
    
    loadAllProbes(msg.cmd, gState)
    for id,gProbe in gState.gprobes.items():
        test = (gState.allProbes.fiberId == id) & (gState.allProbes.holeType == 'GUIDE')
        if test.any(): # should only be one
            gProbe.ugriz = gState.allProbes.mag[test][0]
    
    # TBD: SDSS4: We may have to twiddle with this for coobserved plates.
    # What to do with APOGEEMANGA? Also use the surveyMode?
    gState.setRefractionBalance(gState.plateType, gState.surveyMode)
    
    # Report the cartridge status
    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
#...

def set_decenter(cmd, decenters, gState, enable):
    """Enable/disable decentered guiding, and set offset coordinates."""
    # if we didn't get an enable message, don't allow the position to be changed.
    if enable is None:
        if not gState.decenter:
            failMsg = "Decentered guiding must be enabled before location can be specified."
            cmd.fail("text='%s'"%failMsg)
            return
    
    gState.setDecenter(decenters,cmd,enable)
#...

def set_refraction(cmd, gState, corrRatio, plateType, surveyMode):
    """Set refraction balance to either a specific value or based on plateType/surveyMode."""
    if corrRatio is not None:
        gState.refractionBalance = corrRatio
    elif plateType is not None:
        gState.setRefractionBalance(plateType,surveyMode)

def start_guider(cmd, gState, actorState, queues, camera='gcamera', stack=1,
                 expTime=5, force=False):
    """Start taking and processing exposures with either guider or engineering camera."""

    if gState.cmd:
        errMsg = "The guider appears to already be running"
        if force:
            cmd.warn('text="%s; restarting"' % (errMsg))
        else:
            cmd.fail('text="%s"' % (errMsg))
            return

    if gState.cartridge <= 0:
        failMsg = "No cart/plate information: please load cartridge and try again."
        cmd.fail('guideState=failed; text="%s"' % failMsg)
        return
    if not guidingIsOK(cmd, actorState, force=force):
        failMsg = "Not ok to guide in current state."
        cmd.fail('guideState=failed; text="%s"' % failMsg)
        return

    if (expTime is not None and gState.expTime != expTime):
        gState.expTime = expTime
    if (stack is not None and gState.stack != stack):
        gState.stack = stack

    queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=False))

    gState.cmd = cmd

    # if a start frame was already defined, don't re-define it
    # e.g., make_movie hadn't run successfully.
    if not gState.startFrame:
        # Keep track of the first exposure number for generating movies.
        # Take nextSeqNo+1 because the current value may still be the one
        # issued from the gcamera flat command, which we don't want for this.
        try:
            gState.startFrame = actorState.models['gcamera'].keyVarDict['nextSeqno'][0]+1
            # If 'nextSeqno' hasn't been seen yet (e.g., guider was started after gcamera),
            # we need to get gcamera status first.
        except TypeError:
            cmdVar = actorState.actor.cmdr.call(actor="gcamera", forUserCmd=cmd, cmdStr="status")
            if cmdVar.didFail:
                failMsg = "Cannot get gcamera status to determine nextSeqNo!"
                cmd.fail('guideState=failed; text="%s"' % failMsg)
                return
            # now we can do this safely.
            gState.startFrame = actorState.models['gcamera'].keyVarDict['nextSeqno'][0]+1
        # if we're in simulation mode, use that number instead.
        simulating = actorState.models['gcamera'].keyVarDict['simulating']
        if simulating[0]:
            gState.startFrame = simulating[2]

    gState.reset_pid_terms()
    gState.cmd.respond("guideState=on")
    queues[GCAMERA].put(Msg(Msg.EXPOSE, gState.cmd, replyQueue=queues[MASTER],
                            expTime=gState.expTime, stack=gState.stack, camera=camera))


def stop_guider(cmd, gState, actorState, queues, frameNo, success):
    """Stop current guider exposure and stop taking new exposures."""

    # Try to generate a movie out of the recent guider frames.
    if make_movie(actorState,cmd,gState.startFrame):
        gState.startFrame = None

    if not gState.cmd:
        cmd.fail('text="The guider is already off"')
        return
    
    if success:
        cmd.respond("guideState=stopping")
        # cleanup any pending decenter commands (e.g. decenter off)
        send_decenter_status(gState.cmd,gState,frameNo)
        gState.finish_decenter()
        
        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, cmd, quiet=True, priority=Msg.MEDIUM))
        if gState.cmd != cmd:
            cmd.finish("guideState=off")
        
        gState.cmd.finish("guideState=off")
        gState.cmd = None
    else:
        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, cmd, quiet=True, priority=Msg.MEDIUM))
        if gState.cmd != cmd:
            cmd.fail("guideState=failed")

        gState.cmd.fail("guideState=failed")
        gState.cmd = None

def main(actor, queues):
    """Main loop for master thread"""

    threadName = "master"

    actorState = guiderActor.myGlobals.actorState
    timeout = actorState.timeout
    force = False                       # guide even if the petals are closed
    oneExposure = False                 # just take a single exposure
    frameInfo = None                    # to catch guideStep's return.
    gState = actorState.gState
    # need to wait a couple seconds to let the models sync up.
    time.sleep(3)
    setPoint = actorState.models["gcamera"].keyVarDict["cooler"][0]
    print 'Initial gcamera setPoint:',setPoint
    guiderImageAnalysis = GuiderImageAnalysis(setPoint)
    
    while True:
        try:
            msg = queues[MASTER].get(timeout=timeout)

            qlen = queues[MASTER].qsize()
            if qlen > 0 and msg.cmd:
                msg.cmd.diag("text=master thread has %d items after a .get()" % (qlen))
                
            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))

                return

            elif msg.type == Msg.CENTERUP:
                # Arrange for the next exposure to do a centerUp.
                if not gState.cmd:
                    msg.cmd.fail('text="The guider must be running in order to centerUp"')
                    continue
                else:
                    gState.centerUp = msg.cmd # Provide some way for this command to be finished.

                continue

            elif msg.type == Msg.STOP_GUIDING:
                success = getattr(msg,'success',True) # Succeed, unless told otherwise
                frameNo = getattr(frameInfo,'frameNo',None)
                stop_guider(msg.cmd, gState, actorState, queues, frameNo, success)

            elif msg.type == Msg.START_GUIDING:
                expTime = getattr(msg,'expTime',None)
                stack = getattr(msg,'stack',None)
                force = getattr(msg,'force',False)
                camera = getattr(msg,'camera','gcamera')
                oneExposure = getattr(msg,'oneExposure',False)
                start_guider(msg.cmd, gState, actorState, queues, camera=camera, stack=stack, expTime=expTime, force=force)

            elif msg.type == Msg.REPROCESS_FILE:
                processOneProcFile(gState, msg.filename, actor, queues, cmd=msg.cmd)
                msg.cmd.finish('text="I do hope that succeeded."')
                
            elif msg.type == Msg.READ_PLATE_FILES:
                processOneProcFile(gState, msg.filename, actor, queues, cmd=msg.cmd)
                msg.cmd.finish('text="I do so hope that succeeded."')
                
            elif msg.type == Msg.EXPOSURE_FINISHED:
                if not gState.cmd:    # exposure already finished
                    gState.inMotion = False
                    continue
                
                # TBD: need to check whether the telescope moved here, and
                # ignore this frame if a "tcc offset" was issued.
                # This requires something that monitors tccModel.moveItems[4:]
                # changing to 'Y' so we can flag it, and clear the flag only
                # after we get to this point and have checked it.

                if not msg.success:
                    gState.inMotion = False
                    queues[MASTER].put(Msg(Msg.STOP_GUIDING, gState.cmd, success=False))
                    continue

                camera = getattr(msg,'camera','gcamera')

                frameInfo = guideStep(actor, queues, msg.cmd, gState, msg.filename, oneExposure, guiderImageAnalysis, camera=camera)
                gState.inMotion = False
                
                # output the keywords after the decenter changes, and finish the command.
                if gState.decenterCmd:
                    send_decenter_status(gState.cmd,gState,frameInfo.frameNo)
                    gState.finish_decenter()
                    gState.reset_pid_terms()
                
                # Declare the centerUp to be finished.
                if gState.centerUp:
                    gState.centerUp.finish()
                    gState.centerUp = False
                    gState.reset_pid_terms()
                    # Stuff has changed; tell STUI.
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))
                    
                if not gState.cmd:    # something fatal happened in guideStep
                    continue

                #
                # Is there anything to indicate that we shouldn't be guiding?
                #
                if not guidingIsOK(msg.cmd, actorState, force=force):
                    queues[MASTER].put(Msg(Msg.STOP_GUIDING, gState.cmd))
                    continue
                #
                # Start the next exposure
                #
                if oneExposure:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                    gState.cmd = None
                else:
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, gState.cmd, replyQueue=queues[MASTER],
                                            expTime=gState.expTime, stack=gState.stack, camera=camera))
                
            elif msg.type == Msg.TAKE_FLAT:
                if gState.cartridge <= 0:
                    msg.cmd.fail('text="no valid cartridge is loaded"')
                    continue
                camera = 'ecamera' if gState.plateType == 'ecamera' else 'gcamera'
                queues[GCAMERA].put(Msg(Msg.EXPOSE, msg.cmd, replyQueue=queues[MASTER],
                                        expType="flat", expTime=msg.expTime,
                                        cartridge=gState.cartridge, camera=camera))
            
            elif msg.type == Msg.TAKE_DARK:
                camera = 'ecamera' if gState.plateType == 'ecamera' else 'gcamera'
                queues[GCAMERA].put(Msg(Msg.EXPOSE, msg.cmd, replyQueue=queues[MASTER],
                                        expType="dark", expTime=msg.expTime,
                                        stack=msg.stack, camera=camera))
            
            elif msg.type == Msg.DARK_FINISHED:
                if not msg.success:
                    msg.cmd.fail('text="something went wrong when taking the dark"')
                    continue
                dark_finished(msg,guiderImageAnalysis,actorState,gState)
            
            elif msg.type == Msg.FLAT_FINISHED:
                if not msg.success:
                    msg.cmd.fail('text="something went wrong when taking the flat"')
                    continue
                flat_finished(msg,guiderImageAnalysis,actorState,gState)

            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.startFrame = None # clear the start frame: don't need it any more!
                load_cartridge(msg, queues, gState, actorState)
                    
            elif msg.type == Msg.SET_PID:
                gState.pid[msg.axis].setPID(Kp=msg.Kp, Ti=msg.Ti, Td=msg.Td, Imax=msg.Imax, nfilt=msg.nfilt)

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_REFRACTION:
                set_refraction(msg.cmd, gState, msg.corrRatio, msg.plateType, msg.surveyMode)
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.STAR_IN_FIBER:
                if gState.allProbes == None:
                    msg.cmd.fail('text="the probes for this plate are not available"')
                    continue

                w = None
                if msg.probe:
                    w = numpy.where((gState.allProbes.spectrographId == 2) &
                                    (gState.allProbes.fiberId == msg.probe) &
                                    (gState.allProbes.holeType == 'OBJECT'))
                    w = w[0]
                elif msg.gprobe:
                    w = numpy.where((gState.allProbes.fiberId == msg.gprobe) &
                                    (gState.allProbes.holeType == 'GUIDE'))
                    w = w[0]
                if w == None or len(w) != 1:
                    msg.cmd.fail('text="no unique destination probe was specified"')
                    continue
                dstProbe = gState.allProbes[w]
                dstX = dstProbe.xFocal
                dstY = dstProbe.yFocal
                
                w = None
                if msg.fromProbe:
                    w = numpy.where((gState.allProbes.spectrographId == 2) &
                                    (gState.allProbes.fiberId == msg.fromProbe) &
                                    (gState.allProbes.holeType == 'OBJECT'))
                    w = w[0]
                    if len(w) != 1:
                        msg.cmd.fail('text="no unique source probe was specified"')
                        continue
                elif msg.fromGprobe:
                    w = numpy.where((gState.allProbes.fiberId == msg.fromGprobe) &
                                    (gState.allProbes.holeType == 'GUIDE'))
                    w = w[0]
                    if len(w) != 1:
                        msg.cmd.fail('text="no unique source probe was specified"')
                        continue
                if w != None:
                    srcProbe = gState.allProbes[w]
                    srcX = srcProbe.xFocal
                    srcY = srcProbe.yFocal
                else:
                    srcProbe = None
                    srcX, srcY = 0.0, 0.0

                dx = (dstX - srcX) / gState.plugPlateScale
                dy = (dstY - srcY) / gState.plugPlateScale

                # OK. In all cases disable corrections.
                # For the gprobe case turn on the guide loop.
                for axis in ('axes', 'focus', 'scale'):
                    gState.setGuideMode(axis, False)
                if msg.gprobe:
                    actorState.queues[MASTER].put(Msg(Msg.START_GUIDING, cmd=msg.cmd, force=True))
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))

                msg.cmd.warn('text="offsetting by dy, dx = %g,%g (%g, %g,%g"' %
                             (dy, dx, gState.plugPlateScale,
                              (dstY - srcY), (dstX - srcX)))
                if True:
                    cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=msg.cmd,
                                                        cmdStr="offset bore %g,%g /pabs/computed" % (dx, dy))
                    if cmdVar.didFail:
                        if guidingIsOK(msg.cmd, actorState):
                            msg.cmd.warn('text="Failed to offset, but axes are bypassed"')
                        else:
                            gState.inMotion = False
                            msg.cmd.fail('text="Failed to offset"')
                            continue
                
                msg.cmd.finish()
                
            elif msg.type == Msg.SET_GUIDE_MODE:
                gState.setGuideMode(msg.what, msg.enable)
                #
                # Report the cartridge status
                #
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.ENABLE_FIBER:
                if gState.plate < 0:
                    msg.cmd.fail("text=\"no plate is loaded\"")
                    continue
                try:
                    gState.setGprobeState(msg.fiber, enable=msg.enable)
                except KeyError:
                    msg.cmd.fail('text="Unknown fiber id or fiber type: %s."'%str(msg.fiber))

            elif msg.type == Msg.CHANGE_SCALE:
                """ Change telescope scale by a factor of (1 + 0.01*delta), or to scale
                    We want to do this here, in the guider, so that we can readily track
                    when to ignore new exposures.
                """

                cmd = msg.cmd

                scale = actorState.models["tcc"].keyVarDict["scaleFac"][0]
                if "delta" in cmd.cmd.keywords:
                    delta = float(cmd.cmd.keywords["delta"].values[0])

                    newScale = (1 + 0.01*delta)*scale
                else:
                    newScale = float(cmd.cmd.keywords["scale"].values[0])

                gState.inMotion = True  # Alert the end of exposure processing to skip one.
                cmd.inform('text="currentScale=%g  newScale=%g"' % (scale, newScale))
                cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                                    cmdStr="set scale=%.8f" % (newScale))
                if cmdVar.didFail:
                    gState.inMotion = False
                    cmd.fail('text="Failed to set scale"')
                else:
                    gState.reset_pid_terms(['focus','scale'])
                    cmd.finish('text="scale change completed"')

            elif msg.type == Msg.SET_SCALE:
                gState.setScales(plugPlateScale=msg.plugPlateScale,
                                 gcameraMagnification=msg.gcameraMagnification,
                                 gcameraPixelSize=msg.gcameraPixelSize,
                                 dSecondary_dmm=msg.dSecondary_dmm)
                
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_TIME:
                gState.expTime = msg.expTime
                gState.stack = getattr(msg,'stack',1)
                # this should be set once at thread start in guiderActor_main()
                if hasattr(msg,'readTime'):
                    gState.readTime = msg.readTime

                for k in gState.pid.keys():
                    # camera read time happens for each frame in a stack.
                    gState.pid[k].setPID(dt=((gState.expTime+gState.readTime)*gState.stack + 5)) # "+ 5" to allow for overhead

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                
            elif msg.type == Msg.DECENTER:
                enable = getattr(msg,'enable',None)
                decenters = getattr(msg,'decenters',{})
                set_decenter(msg.cmd, decenters, gState, enable)
            
            elif msg.type == Msg.STATUS:
                # Try to generate status even after we have failed.
                cmd = msg.cmd if msg.cmd.alive else actor.bcast
                
                cmd.respond("cartridgeLoaded=%d, %d, %s, %d, %d" % (
                    gState.cartridge, gState.plate, gState.pointing, gState.fscanMJD, gState.fscanID))
                cmd.respond("survey=%s, %s"%(qstr(gState.plateType), qstr(gState.surveyMode)))
                
                try:
                    if not msg.full:
                        if msg.finish:
                            msg.cmd.finish()
                        continue
                except AttributeError:
                    pass

                cmd.respond("guideState=%s" % ("on" if gState.cmd else "off"))
                cmd.inform('text="The guider is %s"' % ("running" if gState.cmd else "off"))
                cmd.inform('text="Decentering is %s"' % ("off" if not gState.decenter else "on"))
                if gState.decenter:
                    frameNo = getattr(frameInfo,'frameNo',-2)
                    send_decenter_status(cmd,gState,frameNo)

                # Some fiber IDs may be absent from gprobeBits.keys(), so start them all with UNKNOWN
                liveProbes = gState.gprobes.keys()
                if liveProbes:
                    gprobeBits = [GuiderState.UNKNOWN,]*(1 + max(liveProbes))
                    for gProbe in gState.gprobes.values():
                        if gProbe:
                            gprobeBits[gProbe.id] = "0x%02x"%gProbe.gprobebits
                    cmd.respond("gprobeBits=%s" % ", ".join(gprobeBits[1:]))
                
                cmd.respond("guideEnable=%s, %s, %s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                cmd.respond("expTime=%g" % (gState.expTime))
                cmd.respond("stack=%g" % (gState.stack))
                cmd.respond("scales=%g, %g, %g, %g" % (gState.plugPlateScale,
                                                       gState.gcameraMagnification, gState.gcameraPixelSize,
                                                       gState.dSecondary_dmm,))
                gState.output_pid(cmd)

                if gState.refractionBalance != 0.0:
                    cmd.warn('refractionBalance=%0.1f' % (gState.refractionBalance))
                else:
                    cmd.respond('refractionBalance=%0.1f' % (gState.refractionBalance))
                cmd.diag('text="design_ha=%0.1f"' % (gState.design_ha))

                if msg.finish:
                    cmd.finish()
            else:
                raise ValueError("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag('text="%s alive"' % threadName)
        except Exception as e:
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            if gState.cmd:
                gState.cmd.error('text="%s"' % errMsg)
            else:
                actor.bcast.error('text="%s"' % errMsg)
            gState.cmd = None
            # jkp NOTE: this gives a RunTimeError/max recurision depth if
            # guiderActor isn't built correctly (missing lib/libguide.so)
            tback.tback(errMsg, e)

            #import pdb; pdb.set_trace()
            try:
                print "\n".join(tback.tback(errMsg, e)[0]) # old versions of tback return None
            except:
                pass

            try:
                msg.replyQueue.put(Msg.EXIT, cmd=msg.cmd, success=False)
            except Exception as e:
                pass

def guidingIsOK(cmd, actorState, force=False):
    """Is it OK to be guiding?"""

    if force:
        return True

    bypassedNames = actorState.models["sop"].keyVarDict["bypassedNames"]

    ffsStatus = actorState.models["mcp"].keyVarDict["ffsStatus"]
    open, closed = 0, 0
    for s in ffsStatus:
        if s == None:
            cmd.warn('text="Failed to get state of flat field screen from MCP"')
            break

        open += int(s[0])
        closed += int(s[1])

    if open != 8:
        msg = "FF petals aren\'t all open"
        if 'ffs' in bypassedNames:
            cmd.warn('text="%s; guidingIsOk failed, but ffs is bypassed in sop"' % msg)
        else:
            cmd.warn('text="%s; aborting guiding"' % msg)
            return False
    
    # This lets guiderImageAnalysis know to ignore dark frames.
    actorState.bypassDark = 'guider_dark' in bypassedNames
    
#   should we allow guiding with lamps on if axes are disabled
#   check if lamps are actually ON
    ffLamp = actorState.models["mcp"].keyVarDict["ffLamp"]
    hgCdLamp = actorState.models["mcp"].keyVarDict["hgCdLamp"]
    neLamp = actorState.models["mcp"].keyVarDict["neLamp"]
    if (any(ffLamp) and 'lamp_ff' not in bypassedNames) or \
       (any(hgCdLamp) and 'lamp_hgcd' not in bypassedNames) or \
       (any(neLamp) and 'lamp_ne' not in bypassedNames):
        cmd.warn('text="Calibration lamp on; aborting guiding"')
        return False

#   check if non sensed lamps are commanded ON
    uvLamp = actorState.models["mcp"].keyVarDict["uvLampCommandedOn"]
    whtLamp = actorState.models["mcp"].keyVarDict["whtLampCommandedOn"]
    if uvLamp.getValue() or whtLamp.getValue():
        cmd.warn('text="Calibration lamp commanded on; aborting guiding"')
        return False

    tccModel = actorState.models['tcc']
    axisCmdState = tccModel.keyVarDict['axisCmdState']
    if any(x.lower() != 'tracking' for x in axisCmdState):
        if 'axes' in bypassedNames:
            cmd.warn('text="TCC motion failed, but axis motions are bypassed in sop"')
        else:
            cmd.warn('text="TCC motion aborted guiding"')
            return False

    return True
