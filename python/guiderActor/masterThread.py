import Queue, threading
import math, numpy, re

from guiderActor import *
import loadGprobes
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
import opscore.utility.tback as tback

import PID

from gimg.guiderImage import GuiderImageAnalysis

import pyfits
import os.path

try:
    import sm
except ImportError:
    print "Failed to import SM"
    sm = None

class GuiderState(object):
    """Save the state of the guider"""

    class Gprobe(object):
        def __init__(self, id, info, enable=True, flags=None):
            self.id = id
            self.info = info
            self.enabled = enable
            self.flags = flags

        def isEnabled(self):
            raise NotImplementedError()

        def setEnabled(self, enabled):
            if enabled:
                self.flags |= self._ENABLED
            else:
                self.flags &= ~self._ENABLED

    def __init__(self):
        self.cartridge = -1
        self.plate = -1
        self.pointing = "?"
        self.expTime = 0
        self.inMotion = False
        self.guideCmd = None
		
        self.fscanMJD = self.fscanID = -1

        self.deleteAllGprobes()

        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")

        self.pid = {}               # PIDs for various axes
        for what in ["raDec", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

        self.setDecenter("decenterRA")       #store in frameInfo
        self.setDecenter("decenterDec")      
        self.setDecenter("decenterRot")
        self.decenter = False                #store in gState

    def deleteAllGprobes(self):
        """Delete all fibers """
        self.gprobes = {}

    def setGprobeState(self, fiber, enable=True, info=None, create=False, flags=None):
        """Set a fiber's state"""

        if fiber in ("ACQUIRE", "GUIDE"):
            fiber_type = fiber
            for gp in self.gprobes.values():
                if gp.info.fiber_type == fiber_type:
                    gp.enabled = enable
        else:
            if not self.gprobes.has_key(fiber) and create:
                self.gprobes[fiber] = GuiderState.Gprobe(fiber, info, enable, flags)
            else:
                self.gprobes[fiber].enabled = enable

    def setGuideMode(self, what, enabled=True):
        if what == "axes":
            self.guideAxes = enabled
        elif what == "focus":
            self.guideFocus = enabled
        elif what == "scale":
            self.guideScale = enabled
        else:
            raise RuntimeError, ("Unknown guide mode %s" % what)

    def setCmd(self, cmd=None):
        self.guideCmd = cmd
		
    def setDecenter(self, what, value=0):
        if what == "decenterRA":
            self.decenterRA = value
        elif what == "decenterDec":
            self.decenterDec = value
        elif what == "decenterRot":
            self.decenterRot = value
        else:
            raise RuntimeError, ("Unknown decenter axis name %s" % what)

try:
    gState
except:
    gState = GuiderState()

class FrameInfo(object):
    """ Gather all info about the guiding . """

    def __init__(self):
        self.dRA = numpy.nan
        self.dDec = numpy.nan
        self.dRot = numpy.nan
        self.dFocus = numpy.nan
        self.dScale = numpy.nan

        self.filtRA = numpy.nan
        self.filtDec = numpy.nan
        self.filtRot = numpy.nan
        self.filtFocus = numpy.nan
        self.filtScale = numpy.nan

        self.offsetRA = numpy.nan
        self.offsetDec = numpy.nan
        self.offsetRot = numpy.nan
        self.offsetFocus = numpy.nan
        self.offsetScale = numpy.nan

        self.guideCameraScale = numpy.nan
        self.plugPlateScale = numpy.nan
        self.seeing = numpy.nan
        self.guideRMS= numpy.nan
        self.decenterRA = numpy.nan
        self.decenterDec = numpy.nan
        self.decenterRot = numpy.nan

def postscriptDevice(psPlotDir, frameNo, prefix=""):
    """Return the SM device to write the postscript file for guide frame frameNo"""
    return "postencap %s%d.eps" % (prefix, frameNo)

class FakeCommand(object):
    def _respond(self, tag, text):
        print "%s %s" % (tag, text)
    def warn(self, text):
        self._respond('w', text)
    def respond(self, text):
        self._respond('i', text)
    def finish(self, text):
        self._respond(':', text)
    def fail(self, text):
        self._respond('f', text)

def processOneFile(guiderFile, cartFile, plateFile):
    
    gState.setGuideMode('axes', False)
    gState.setGuideMode('focus', False)
    gState.setGuideMode('scale', False)

    cmd = FakeCommand()
    queues = dict(MASTER=Queue.Queue())

    guideStep(None, queues, cmd, cmd, guiderFile, True)

def processOneProcFile(guiderFile, cartFile, plateFile, actor=None, queues=None, cmd=None, guideCmd=None):
    gState.setGuideMode('axes', False)
    gState.setGuideMode('focus', False)
    gState.setGuideMode('scale', False)
	
    if not cmd: cmd = FakeCommand()
    if not guideCmd: guideCmd = FakeCommand()
    if not queues: queues = dict(MASTER=Queue.Queue())

    gState.setCmd(guideCmd)
    guideStep(None, queues, cmd, cmd, guiderFile, True)

def guideStep(actor, queues, cmd, inFile, oneExposure,
              plot=False, psPlot=False, sm=False):
    """ One step of the guide loop, based on the given guider file. 

    Args: (TOOOO MANY!!)
        actor      - 
"""
    sigmaToFWHM = 2.354 # conversion for a Gaussian, use this eveywhere but in ipGguide.c
    #conversion from sigma to FWHM for a JEG double Gaussian is done in ipGguide.c (sigmaToFWHMJEG = 2.468)

    psPlotDir  = "/data/gcam/scratch/"
    minStarFlux = 2000 #ADU, avoid guiding on noise spikes during acquisitions
                       #should be in photons, based on RON, Dark residual, SKY   

    actorState = guiderActor.myGlobals.actorState
    guideCmd = gState.guideCmd
    guideCmd.respond("processing=%s" % inFile)
    frameNo = int(re.search(r"([0-9]+)\.fits$", inFile).group(1))

    h = pyfits.getheader(inFile)
    flatfile = h.get('FLATFILE', None)
    flatcart = h.get('FLATCART', None)
    darkfile = h.get('DARKFILE', None)
    if not flatfile:
        guideCmd.fail("text=%s" % qstr("No flat image available"))
        gState.setCmd(None)
        return
    
    if not darkfile:
        guideCmd.fail("text=%s" % qstr("No dark image available"))
        gState.setCmd(None)
        return

    if flatcart != gState.cartridge:
        if False:
            guideCmd.fail("text=%s" % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                flatcart, gState.cartridge)))
            gState.setCmd(None)
            return
        else:
            guideCmd.warn("text=%s" % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                flatcart, gState.cartridge)))

    try:
        guideCmd.inform("text='GuiderImageAnalysis()...'")
        GI = GuiderImageAnalysis(inFile, cmd=guideCmd)
        guideCmd.inform("text='GuiderImageAnalysis.findFibers()...'")
        fibers = GI.findFibers(gState.gprobes)
        guideCmd.inform("text='GuiderImageAnalysis.findFibers() got %i fibers'" % len(fibers))
    except Exception, e:
        tback.tback("GuideTest", e)
        guideCmd.fail("text=%s" % qstr("Error in processing guide images: %s" % e))
        gState.setCmd(None)
        return

    # Object to gather all per-frame guiding info into.
    frameInfo = FrameInfo()

    # Setup to solve for the axis and maybe scale offsets.  We work consistently
    # in mm on the focal plane, only converting to angles to command the TCC
    #
    # N.B. fiber.xFocal and fiber.yFocal are the offsets of the stars
    # wrt the center of the plate in mm; fiber.xcen/star.xs are in pixels,
    # so we need a scale for the guide camera.  Nominally the guide camera
    # has the same scale as the plug plate itself, but maybe it doesn't,
    # so we'll include a possible magnification
    #
    guideCameraScale = gState.gcameraMagnification*gState.gcameraPixelSize*1e-3 # mm/pixel
    frameInfo.guideCameraScale = guideCameraScale
    frameInfo.plugPlateScale = gState.plugPlateScale

    A = numpy.matrix(numpy.zeros(3*3).reshape([3,3]))
    b = numpy.matrix(numpy.zeros(3).reshape([3,1])); b3 = 0.0

    if plot or sm:                           #setup arrays for sm
        size = len(gState.gprobes) + 1 # fibers are 1-indexed
        fiberid_np = numpy.zeros(size)
        raCenter_np = numpy.zeros(size)
        decCenter_np = numpy.zeros(size)
        dRA_np = numpy.zeros(size)
        dDec_np = numpy.zeros(size)

    guideRMS    = 0.0
    guideXRMS   = 0.0
    guideYRMS   = 0.0
    nguideRMS   = 0
    guideRaDecRMS = 0.0
    guideRaRMS = 0.0
    guideDecRMS = 0.0
    nguideRaDec = 0
    guideFitRMS = 0.0        #not calculated at present

    for fiber in fibers:
        # necessary?
        if fiber.gprobe is None:
            guideCmd.warn('text="Gprobe %d was not listed in plugmap info"' % fiber.fiberid)
            continue
        gp = fiber.gprobe
        probe = gp.info
        enabled = gp.enabled

        #
        # dx, dy are the offsets on the ALTA guider image
        #
        fiber.dx = guideCameraScale*(fiber.xs - fiber.xcen) + (probe.xFerruleOffset / 1000.)
        fiber.dy = guideCameraScale*(fiber.ys - fiber.ycen) + (probe.yFerruleOffset / 1000.)
        poserr = fiber.xyserr

        isnan = numpy.isnan
        if isnan(fiber.dx) or isnan(fiber.dy) or isnan(poserr):
            guideCmd.warn("text=%s" %
                          qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber measured=(%g, %g), nominal=(%g,%g)" % (
                              fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, probe.xCenter, probe.yCenter)))
            continue

        if fiber.flux < minStarFlux:
            guideCmd.warn("text=%s" %
                          qstr("Star in gprobe %d too faint for guiding flux %g < %g minimum flux" % (
                              fiber.fiberid, fiber.flux, minStarFlux)))
            continue
        if poserr == 0:
            guideCmd.warn("text=%s" %
                          qstr("position error is 0 for gprobe %d star=(%g, %g) fiber=(%g, %g) nominal=(%g,%g)" % (
                              fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, probe.xCenter, probe.yCenter)))
            continue


        #
        # theta is the angle to rotate (x, y) on the ALTA to (ra, alt)
        #
        # phi is the orientation of the alignment hole measured clockwise from N
        # rotation is the anticlockwise rotation from x on the ALTA to the pin
        #
        theta = 90                   # allow for 90 deg rot of camera view, should be -90 
        theta += probe.rotation # allow for intrinsic fibre rotation
        theta -= probe.phi      # allow for orientation of alignment hole

        probe.rotStar2Sky = theta # Squirrel the real angle away.
        theta = math.radians(theta)
        ct, st = math.cos(theta), math.sin(theta)
        # error in guide star position; n.b. still in mm here
        dRA   =  fiber.dx*ct + fiber.dy*st
        dDec  = -fiber.dx*st + fiber.dy*ct
        dDec *= -1

        # Apply RA & Dec user guiding offsets to mimic different xy fibers centers
        # The guiderRMS will be calculated around the new effective fiber centers
      
        if gState.decenter: 
            #convert decenter offset to mm on guider
            frameInfo.decenterRA  = gState.decenterRA*gState.plugPlateScale/3600
            frameInfo.decenterDec = gState.decenterDec*gState.plugPlateScale/3600
            frameInfo.decenterRot = gState.decenterRot            
            # apply decenter offset so that telescope (not the star) moves the decenter amount
            dRA  += frameInfo.decenterRA
            dDec += frameInfo.decenterDec
            #decenterRot applied after guide solution
        else:
            frameInfo.decenterRA = 0.0
            frameInfo.decenterRA = 0.0
            frameInfo.decenterRot = 0.0

        fiber.dRA = dRA
        fiber.dDec = dDec
        raCenter  = probe.xFocal
        decCenter = probe.yFocal

        if plot:
            try:
                fiberid_np[fiber.fiberid] = fiber.fiberid
                raCenter_np[fiber.fiberid] = raCenter
                decCenter_np[fiber.fiberid] = decCenter
                dRA_np[fiber.fiberid] = dRA
                dDec_np[fiber.fiberid] = dDec
            except IndexError, e:
                #import pdb; pdb.set_trace()
                pass

        if True:
            guideCmd.inform("probe=%d,%2d,0x%02d, %7.2f,%7.2f, %7.3f,%4.0f, %7.2f,%6.2f,%7.2f,%6.2f" % (
                frameNo, fiber.fiberid, probe.flags,
                3600.0*(fiber.dRA/gState.plugPlateScale), 3600.0*(fiber.dDec/gState.plugPlateScale),
                fiber.fwhm, probe.focusOffset,
                fiber.flux, fiber.mag, fiber.sky, fiber.skymag))

            print "%d %2d  %7.2f %7.2f  %7.2f %7.2f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %06.1f  %7.3f %7.3f %7.3f %7.3f %4.0f" % (
                frameNo,
                fiber.fiberid, dRA, dDec, fiber.dx, fiber.dy, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen,
                probe.xFocal, probe.yFocal, probe.rotStar2Sky, fiber.fwhm/sigmaToFWHM, fiber.sky, fiber.flux, fiber.mag,
                probe.focusOffset)

        if not enabled:
            continue

        #accumulate guiding errors for good stars used in fit
        guideRMS += fiber.dx**2 + fiber.dy**2
        guideXRMS += fiber.dx**2
        guideYRMS += fiber.dy**2        
        nguideRMS += 1

        b[0] += dRA
        b[1] += dDec
        b[2] += raCenter*dDec - decCenter*dRA

        A[0, 0] += 1
        A[0, 1] += 0
        A[0, 2] += -decCenter

        A[1, 0] += 0
        A[1, 1] += 1
        A[1, 2] += raCenter

        A[2, 2] += raCenter*raCenter + decCenter*decCenter
        #
        # Now scale.  We don't actually solve for scale and axis updates
        # simultanously, and we don't allow for the axis update when
        # estimating the scale. 
        #
        b3 += raCenter*dRA + decCenter*dDec

    nStar = A[0, 0]
    if nStar == 0 or gState.inMotion:
        guideCmd.warn('text="No stars are available for guiding or guiding is deferred"')
        guideRMS = 99.99
        guideXRMS = 99.99
        guideYRMS = 99.99

        GI.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)

        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.setCmd(None)
            return

        if guidingIsOK(cmd, actorState):
            queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
                                    expTime=gState.expTime))
        return
        
    A[2, 0] = A[0, 2]
    A[2, 1] = A[1, 2]
    try:
        if nStar == 1:
            guideCmd.warn('text="Only one star is usable"')
            x = b
            x[2, 0] = 0 # no rotation
            guideRMS = 99.99
            guideXRMS = 99.99
            guideYRMS = 99.99
        else:
            x = numpy.linalg.solve(A, b)

        # convert from mm to degrees
        dRA = x[0, 0]/gState.plugPlateScale
        dDec = x[1, 0]/gState.plugPlateScale
        dRot = -math.degrees(x[2, 0]) # and from radians to degrees

        #add the decenter guiding rotation offset here for now
        if gState.decenter:
            dRot += frameInfo.decenterRot

        frameInfo.dRA  = dRA
        frameInfo.dDec = dDec
        frameInfo.dRot = dRot
        print 'dRA,dDec,dRot', dRA, dDec, dRot

        offsetRa  = -gState.pid["raDec"].update(dRA)                    
        offsetDec = -gState.pid["raDec"].update(dDec)
        offsetRot = -gState.pid["rot"].update(dRot) if nStar > 1 else 0 # don't update I
        print 'offsetRA, offsetDec, offsetRot', offsetRa, offsetDec, offsetRot

        frameInfo.filtRA  = offsetRa
        frameInfo.filtDec = offsetDec
        frameInfo.filtRot = offsetRot
        frameInfo.offsetRA  = offsetRa  if gState.guideAxes else 0.0
        frameInfo.offsetDec = offsetDec if gState.guideAxes else 0.0
        frameInfo.offsetRot = offsetRot if gState.guideAxes else 0.0

        guideCmd.respond("axisError=%g, %g, %g" % (3600*dRA, 3600*dDec, 3600*dRot))
        guideCmd.respond("axisChange=%g, %g, %g, %s" % (-3600*offsetRa, -3600*offsetDec, -3600*offsetRot,
                                                        "enabled" if gState.guideAxes else "disabled"))

        #rms position error prior to this frames ccrrection                       
        try:
            guideRMS = math.sqrt(guideRMS/nguideRMS) 
            guideXRMS = math.sqrt(guideXRMS/nguideRMS)
            guideYRMS = math.sqrt(guideYRMS/nguideRMS)
        except:
            guideRMS = 99.99
            guideXRMS = 99.99
            guideYRMS = 99.99

        if gState.guideAxes:
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="offset arc %f, %f" % \
                                         (-offsetRa, -offsetDec))

            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue offset"')

            if offsetRot: 
                cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                         cmdStr="offset guide %f, %f, %g" % \
                                             (0.0, 0.0, -offsetRot))

            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue offset in rotator"')

        if sm: 
            if plot:
                try:
                    sm.device('X11')
                except:
                    guideCmd.warn('text="X display error, cannot open sm guider window"')
                    plot = False

        else:
            guideCmd.warn('text="Unable to plot as SM is not available"')
            plot = False

        if plot and psPlot:
            if not (os.path.exists(psPlotDir)):
                psPlot = False
                guideCmd.warn('text="Unable to write SM hardcopies"')

        if psPlot and not plot:
            guideCmd.warn('text="Need to enable both plot & psPlot"')           

        if plot:
            for plotdev in ("X11 -device 0", "postscript"):
                if plotdev == "postscript":
                    if not psPlot:
                        continue
                    deviceCmd = postscriptDevice(psPlotDir, frameNo)
                else:
                    deviceCmd = plotdev

                sm.device(deviceCmd)

                sm.erase(False)
                sm.limits([-400, 400], [-400, 400])
                sm.box()
                sm.frelocate(0.5, 1.03)
                sm.putlabel(5, r"\1Offsets")
                sm.xlabel(r"\2\delta Ra")
                sm.ylabel(r"\2\delta Dec")

                sm.ptype([63])
                sm.frelocate(0.85, 0.95)
                sm.putlabel(5, r"\1Frame %d" % frameNo)
                vscale = 1000 # how much to multiply position error
                sm.relocate(-350, 350)
                asec = gState.plugPlateScale/3600.0
                sm.draw(-350 + vscale*gState.plugPlateScale/3600.0, 350)
                sm.label(r"  \raise-200{\1{1 arcsec}}")

                for i in range(len(fiberid_np)):
                    if fiberid_np[i] == 0:
                        continue

                    sm.relocate(raCenter_np[i], decCenter_np[i])
                    sm.putlabel(6, r" %d" % fiberid_np[i])

                    sm.relocate(raCenter_np[i], decCenter_np[i])

                    sm.dot()

                    if not gState.gprobes[fiberid_np[i]].enabled:
                        sm.ltype(1)
                        sm.ctype(sm.CYAN)

                    sm.draw(raCenter_np[i] + vscale*dRA_np[i],
                            decCenter_np[i] + vscale*dDec_np[i])

                    sm.ltype()
                    sm.ctype()

                sm.ptype()

    except numpy.linalg.LinAlgError:
        guideCmd.warn("text=%s" % qstr("Unable to solve for axis offsets"))

    if nStar <= 1:      # don't bother with focus/scale!
        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.setCmd(None)
            return

        if guidingIsOK(cmd, actorState, force=force):
            queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
                                    expTime=gState.expTime))
        return

    # RMS guiding error 
    print "RMS guiding error= %4.3f, n stars= %d RMS_X= %4.3f, RMS_Y=%4.3f" %(guideRMS, nguideRMS, guideXRMS, guideYRMS) 
    guideCmd.inform('text=" guideRMS =%g"' % (guideRMS))

    #
    # Scale
    #
    dScale = b3/A[2, 2]
    dScaleCorrection = -dScale * 100.    #value for operators to enter manually
    offsetScale = -gState.pid["scale"].update(dScale)

    frameInfo.dScale = dScale
    frameInfo.filtScale = offsetScale
    frameInfo.offsetScale = offsetScale if gState.guideScale else 0.0

    guideCmd.respond("scaleError=%g" % (dScale))
    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                             "enabled" if gState.guideScale else "disabled"))
    guideCmd.inform('text="delta percentage scale correction =%g"' % (dScaleCorrection))
    curScale = actorState.models["tcc"].keyVarDict["scaleFac"][0]

	# There is (not terribly surprisingly) evidence of crosstalk between scale and focus adjustements.
	# So defer focus changes if we apply a scale change.
    blockFocusMove = False
		
    if gState.guideScale and abs(offsetScale) > 1e-6:
        # Clip to the motion we think is too big to apply at once.
        offsetScale = max(min(offsetScale, 5e-6), -5e-6)
        offsetScale += curScale
        cmd.warn('text="setting scale=%0.6f"' % (offsetScale))

        # Last chance to bailout.
        if offsetScale < 0.9995 or offsetScale > 1.0005:
            cmd.warn('text="NOT setting scarily large scale=%0.6f"' % (offsetScale))
        else:
            blockFocusMove = True
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="set scale=%f" % (offsetScale))
            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue scale change"')
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

    # setup arrays for sm and saving to FITS
    size = len(gState.gprobes) + 1 # fibers are 1-indexed
    x_np = numpy.zeros(size)
    xErr_np = numpy.zeros(size)
    d_np = numpy.zeros(size)

    A = numpy.matrix(numpy.zeros(2*2).reshape([2,2]))
    b = numpy.matrix(numpy.zeros(2).reshape([2,1]))

    for fiber in fibers:
        # required?
        if fiber.gprobe is None:
            continue
        gp = gState.gprobes[fiber.fiberid]
        if not gp.enabled:
            continue
        probe = gp.info

                                # FIXME -- do we want to include ACQUISITION fibers?
        rms = fiber.fwhm / sigmaToFWHM
        if isnan(rms):
            continue

        micronsPerArcsec = 1/3600.0*gState.plugPlateScale*1e3 # convert arcsec to microns
        rms *= micronsPerArcsec # in microns
        rmsErr = 1

        d = probe.focusOffset
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

        try:
            fiberid_np[fiber.fiberid] = fiber.fiberid
            x_np[fiber.fiberid] = x
            xErr_np[fiber.fiberid] = xErr
            d_np[fiber.fiberid] = d
        except IndexError, e:
            pass

    A[1, 0] = A[0, 1]
    try:
        x = numpy.linalg.solve(A, b)

        Delta = x[1, 0]/(2*C)
        try:
            rms0 = math.sqrt(x[0, 0] - C*Delta*Delta)/micronsPerArcsec
        except ValueError, e:
            rms0 = float("NaN")

        # Note sign change here.
        dFocus = -Delta*gState.dSecondary_dmm # mm to move the secondary
        offsetFocus = -gState.pid["focus"].update(dFocus)

        frameInfo.dFocus = dFocus
        frameInfo.filtFocus = offsetFocus
        frameInfo.offsetFocus = offsetFocus if gState.guideFocus else 0.0
        frameInfo.seeing = rms0*sigmaToFWHM   #in arc sec

        guideCmd.respond("seeing=%g" % (rms0*sigmaToFWHM))
        guideCmd.respond("focusError=%g" % (dFocus))
        guideCmd.respond("focusChange=%g, %s" % (offsetFocus, "enabled" if (gState.guideFocus and not blockFocusMove) else "disabled"))
        if gState.guideFocus and not blockFocusMove:
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="set focus=%f/incremental" % (offsetFocus))

            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue focus offset"')
    except numpy.linalg.LinAlgError:
        guideCmd.warn("text=%s" % qstr("Unable to solve for focus offset"))
        x = None

    if plot:
        try:
            for plotdev in ("X11 -device 1", "postscript"):
                if plotdev == "postscript":
                    if not psPlot:
                        continue
                    deviceCmd = postscriptDevice(psPlotDir, frameNo, "Focus")
                else:
                    deviceCmd = plotdev

                sm.device(deviceCmd)

                sm.erase(False)
                #
                # Bravely convert to FWHM in arcsec (brave because of the sqrt)
                #
                f = sigmaToFWHM/micronsPerArcsec
                X_np = f*numpy.sqrt(x_np)
                XErr_np = f*xErr_np/(2*numpy.sqrt(x_np))

                sm.limits(d_np, X_np)
                sm.box()
                sm.frelocate(0.5, 1.03)
                sm.putlabel(5, r"\1Focus")

                sm.xlabel(r"\2d_i \equiv{} fibre piston (\mu m)")
                sm.ylabel(r"\2\ssqrt{\sigma_i^2 - C d_i^2} (FWHM, arcsec)")

                sm.frelocate(0.85, 0.95)
                sm.putlabel(5, r"\1Frame %d" % frameNo)

                for i in range(len(fiberid_np)):
                    if fiberid_np[i] == 0:
                        continue

                    sm.relocate(d_np[i], X_np[i])
                    sm.putlabel(6, r" %d" % fiberid_np[i])

                for l in (2, 4):
                    sm.errorbar(d_np, X_np, XErr_np, l)

                #dd_np = numpy.array([-1000, 1000])
                dd_np = numpy.arange(-1000, 1000, 25)

                sm.ctype(sm.CYAN)
                if x != None: # I.e. we successfully fit for focus
                    sm.connect(dd_np, f*numpy.sqrt(x[0, 0] + x[1, 0]*dd_np))
                    sm.frelocate(0.1, 0.1); sm.label(r"\line 0 2000 \colour{default} Best fit")

                    sm.ctype(sm.MAGENTA)
                    sm.frelocate(0.1, 0.15)
                    if False:
                        sm.label(r"\line 1 2000 \colour{default} Seeing")
                        sm.ltype(1)
                        sm.connect(dd_np, 0*dd_np + rms0*sigmaToFWHM)
                    else:
                        sm.label(r"{\2\apoint 45 4 1} \colour{default} Seeing")
                        sm.relocate(0, rms0*sigmaToFWHM)
                        sm.expand(4); sm.angle(45)
                        sm.dot()
                        sm.expand(); sm.angle(0)

                sm.ltype(); sm.ctype()

                del dd_np
        except Exception, e:
            guideCmd.warn('text="plot failed: %s"' % (e))

            #
            # Write output fits file for TUI
            #
    GI.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)


def main(actor, queues):
    """Main loop for master thread"""

    threadName = "master"

    actorState = guiderActor.myGlobals.actorState
    timeout = actorState.timeout
    force = False                       # guide even if the petals are closed
    oneExposure = False                 # just take a single exposure
    plot = False                        # use SM to plot things
    psPlot = False
    fakeSpiderInstAng = None            # the value we claim for the SpiderInstAng
    
    minStarFlux = 2000 #ADU, avoid guiding on noise spikes during acquisitions
                       #FIXME defined here and in GIA

    while True:
        try:
            msg = queues[MASTER].get(timeout=timeout)

            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))

                return
            
            elif msg.type == Msg.START_GUIDING:
                if not msg.start:
                    if msg.start is None:
                        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                        continue
                        
                    if not gState.guideCmd:
                        msg.cmd.fail('text="The guider is already off"')
                        continue

                    msg.cmd.respond("guideState=stopping")
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                    gState.guideCmd.finish("guideState=off")
                    gState.setCmd(None)
                    msg.cmd.finish()
                    continue

                try:
                    expTime = msg.expTime
                    
                    if expTime >= 0 and gState.expTime != expTime:
                        gState.expTime = expTime
                        queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))
                except AttributeError:
                    pass

                force = msg.force
                oneExposure = msg.oneExposure
                plot = msg.plot
                psPlot = msg.psPlot
                fakeSpiderInstAng = msg.spiderInstAng
                
                if gState.guideCmd:
                    errMsg = "The guider appears to already be running"
                    if force:
                        msg.cmd.warn('text="%s; restarting"' % (errMsg))
                    else:
                        msg.cmd.fail('text="%s"' % (errMsg))
                        continue

                if gState.plate < 0:
                    queues[MASTER].put(Msg(Msg.FAIL, msg.cmd,
                                           text="Please tell me about your cartridge and try again"))
                    continue

                if not guidingIsOK(msg.cmd, actorState, force=force):
                    queues[MASTER].put(Msg(Msg.FAIL, msg.cmd, text=""))
                    continue

                guideCmd = msg.cmd
                gState.setCmd(guideCmd)
				
                #
                # Reset any PID I terms
                #
                for key in gState.pid.keys():
                    gState.pid[key].reset()

                guideCmd.respond("guideState=on")
                if msg.decenter: 
                    gState.decenter = True 
                else: 
                    gState.decenter = False

                queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER], expTime=gState.expTime))

            elif msg.type == Msg.REPROCESS_FILE:
                processOneProcFile(msg.filename, actor, queues, cmd=msg.cmd)
                msg.cmd.finish('text="I do hope that succeeded."')
                
            elif msg.type == Msg.READ_PLATE_FILES:
                processOneProcFile(msg.filename, actor, queues, cmd=msg.cmd)
                msg.cmd.finish('text="I do so hope that succeeded."')
                
            elif msg.type == Msg.TCC_EXPOSURE:
                queues[GCAMERA].put(Msg(Msg.EXPOSE, msg.cmd, replyQueue=queues[MASTER], 
                                        expTime=msg.expTime, forTCC=msg.forTCC, camera=msg.camera))

            elif msg.type == Msg.EXPOSURE_FINISHED:
                if msg.forTCC:
                    tccState = msg.forTCC

                    if not msg.success:
                        tccState.doreadFilename = None
                        msg.cmd.warn('text="exposure failed"')
                        msg.cmd.finish('txtForTcc=" OK"')
                        continue
                        
                    tccState.doreadFilename = msg.filename
                    ccdTemp = 0.0   # self.camera.cam.read_TempCCD()
                    msg.cmd.respond('txtForTcc=%s' % (qstr('%d %d %0.1f %0.1f %0.1f %0.1f %0.2f %d %0.2f %s' % \
                                                               (tccState.binX, tccState.binY,
                                                                tccState.ctrX, tccState.ctrY,
                                                                tccState.sizeX, tccState.sizeY,
                                                                tccState.itime, tccState.gImCamID, ccdTemp,
                                                                "image: binXY begXY sizeXY expTime camID temp"))))
                    msg.cmd.finish('txtForTcc=" OK"')
                    continue

                if not gState.guideCmd:    # exposure already finished
                    gState.inMotion = False
                    continue

                if not msg.success:
                    gState.inMotion = False
                    queues[MASTER].put(Msg(Msg.START_GUIDING, gState.guideCmd, start=False))
                    continue

                guideStep(actor, queues, msg.cmd, msg.filename, oneExposure,
                          plot=plot, psPlot=psPlot, sm=sm)
                gState.inMotion = False
                if not gState.guideCmd:    # something fatal happened in guideStep
                    continue

                #
                # Is there anything to indicate that we shouldn't be guiding?
                #
                #FIXME Is this the best place to clear decenter offsets
                #Check if offsets were changed during guiding without force
                if gState.decenter and not force:
                    cmd.fail('text="Decentred guiding must use force."')
                    gState.setDecenter("decenterRA")  #reset all to 0
                    gState.setDecenter("decenterDec")
                    gState.setDecenter("decenterRot")
                    gState.decenter = False
                    queues[MASTER].put(Msg(Msg.START_GUIDING, gState.guideCmd, start=False))
                    continue

                if not guidingIsOK(msg.cmd, actorState, force=force):
                    queues[MASTER].put(Msg(Msg.START_GUIDING, gState.guideCmd, start=False))
                    continue
                #
                # Start the next exposure
                #
                if oneExposure:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                    gState.setCmd(None)
                else:
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, gState.guideCmd, replyQueue=queues[MASTER], expTime=gState.expTime))
                
            elif msg.type == Msg.TAKE_FLAT:
                if gState.cartridge <= 0:
                    msg.cmd.fail('text="no valid cartridge is loaded"')
                    continue

                queues[GCAMERA].put(Msg(Msg.EXPOSE, msg.cmd, replyQueue=queues[MASTER], 
                                        expType="flat", expTime=msg.expTime, cartridge=gState.cartridge))

            elif msg.type == Msg.FLAT_FINISHED:
                cmd = msg.cmd
                if not msg.success:
                    cmd.fail('text="something went wrong when taking the flat"')
                    continue

                cmd.respond("processing=%s" % msg.filename)
                frameNo = int(re.search(r"([0-9]+)\.fits$", msg.filename).group(1))
                
                h = pyfits.getheader(msg.filename)
                exptype = h.get('IMAGETYP')
                if exptype != "flat":
                    cmd.fail('text="flat image processing ignoring a %s image!!"' % (exptype))
                    continue

                darkfile = h.get('DARKFILE', None)
                if not darkfile:
                    cmd.fail("text=%s" % qstr("No dark image available!!"))
                    continue

                cmd.inform("text=GuiderImageAnalysis()...")
                GI = GuiderImageAnalysis(msg.filename, cmd=cmd)
                cmd.inform("text=GuiderImageAnalysis.findFibers()...")
                try:
                    fibers = GI.findFibers(gState.gprobes)
                    flatoutname = GI.getProcessedOutputName(msg.filename) 
                    dirname, filename = os.path.split(flatoutname)
                    cmd.inform('file=%s/,%s' % (dirname, filename))
                    cmd.finish('text="flat image processing done"')
                except Exception, e:
                    tback.tback("findFibers", e)
                    cmd.fail('text="findFibers failed: %s"' % (e))

                continue
                    
            elif msg.type == Msg.FAIL:
                msg.cmd.fail('text="%s"' % msg.text);

            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllGprobes()

                gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
                gState.fscanMJD, gState.fscanID = msg.fscanMJD, msg.fscanID
                gState.boresight_ra, gState.boresight_dec = msg.boresight_ra, msg.boresight_dec
                #
                # Set the gState.gprobes array (actually a dictionary as we're not sure which fibre IDs are present)
                #
                gState.gprobes = {}
                for id, info in msg.gprobes.items():
                    # FIXABLE HACK: set broken/unplugged probes to be !exists
                    # # The fix is to unify all the probe structures
                    if info.flags & 0x3:
                        info.exists = False
                    if info.exists:
                        enabled = False if info.fiber_type == "TRITIUM" else info.enabled
                    else:
                        enabled = False

                    gState.setGprobeState(id, enable=enabled, info=info, create=True, flags=info.flags)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                
            elif msg.type == Msg.SET_PID:
                gState.pid[msg.what].setPID(Kp=msg.Kp, Ti=msg.Ti, Td=msg.Td, Imax=msg.Imax, nfilt=msg.nfilt)

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
            elif msg.type == Msg.SET_GUIDE_MODE:
                gState.setGuideMode(msg.what, msg.enable)
                #
                # Report the cartridge status
                #
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.ENABLE_FIBER:
                if gState.plate < 0:
                    msg.cmd.error("test=\"no plate is loaded\"")
                    continue
                
                gState.setGprobeState(msg.fiber, enable=msg.enable)

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
                                                    cmdStr="set scale=%.6f" % (newScale))
                if cmdVar.didFail:
                    gState.inMotion = False
                    cmd.fail('text="Failed to set scale"')
                else:
                    gState.pid['focus'].reset()
                    gState.pid['scale'].reset()
                    cmd.finish('text="scale change completed"')

            elif msg.type == Msg.SET_SCALE:
                gState.plugPlateScale = msg.plugPlateScale
                gState.gcameraMagnification = msg.gcameraMagnification
                gState.gcameraPixelSize = msg.gcameraPixelSize
                gState.dSecondary_dmm = msg.dSecondary_dmm
                
                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_TIME:
                gState.expTime = msg.expTime

                for k in gState.pid.keys():
                    gState.pid[k].setPID(dt=(gState.expTime + 5)) # "+ 5" to allow for overhead

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

                #Allow decenters to be setup prior to guiding
            elif msg.type == Msg.DECENTER:
                        gState.setDecenter("decenterRA",msg.decenterRA)
                        gState.setDecenter("decenterDec",msg.decenterDec)
                        gState.setDecenter("decenterRot",msg.decenterRot)

            elif msg.type == Msg.STATUS:
                # Try to generate status even after we have failed.
                cmd = msg.cmd if msg.cmd.alive else actor.bcast

                cmd.respond("cartridgeLoaded=%d, %d, %s, %d, %d" % (
                    gState.cartridge, gState.plate, gState.pointing, gState.fscanMJD, gState.fscanID))

                try:
                    if not msg.full:
                        if msg.finish:
                            msg.cmd.finish()
                        continue
                except AttributeError:
                    pass

                cmd.respond("guideState=%s" % ("on" if gState.guideCmd else "off"))
                cmd.inform('text="The guider is %s"' % ("running" if gState.guideCmd else "off"))
                cmd.inform('text="Decentering is %s"' % ("off" if not gState.decenter else "on"))
                
                fiberState = []
                gprobeBitsDict = {}
                for f in gState.gprobes.values():
                    if f:
                        fiberState.append("\"(%d=%s)\"" % (f.id, f.enabled))
                        gprobeBitsDict[f.id] = ("0x%x" % f.flags)

                if len(fiberState) > 0:
                    cmd.respond("gprobes=%s" % ", ".join(fiberState))

                # Some fiber IDs may be absent from gprobeBits.keys(), so make a filled list
                if gprobeBitsDict:
                    gprobeBits = [0xff,]*(1 + sorted([int(k) for k in gprobeBitsDict.keys()])[-1])
                    for k, f in gprobeBitsDict.items():
                        gprobeBits[k] = f
                    cmd.respond("gprobeBits=%s" % ", ".join(gprobeBits[1:]))
                    

                cmd.respond("guideEnable=%s, %s, %s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                cmd.respond("expTime=%g" % (gState.expTime))
                cmd.respond("scales=%g, %g, %g, %g" % (gState.plugPlateScale,
                                                           gState.gcameraMagnification, gState.gcameraPixelSize,
                                                           gState.dSecondary_dmm,))
                for w in gState.pid.keys():
                    cmd.respond("pid=%s, %g, %g, %g, %g" % (w,
                                                                gState.pid[w].Kp, gState.pid[w].Ti, gState.pid[w].Td,
                                                                gState.pid[w].Imax))
                if msg.finish:
                    cmd.finish()
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag('text="%s alive"' % threadName)
        except Exception, e:
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            actor.bcast.warn('text="%s"' % errMsg)
            gState.setCmd(False)
            # I (dstn) get infinite recursion from this...
            #tback(errMsg, e)

            #import pdb; pdb.set_trace()
            try:
                print "\n".join(tback.tback(errMsg, e)[0]) # old versions of tback return None
            except:
                pass

            try:
                msg.replyQueue.put(Msg.EXIT, cmd=msg.cmd, success=False)
            except Exception, e:
                pass

def guidingIsOK(cmd, actorState, force=False):
    """Is it OK to be guiding?"""

    if force:
        return True

    bypassed = actorState.models["sop"].keyVarDict["bypassed"]
    bypassNames = actorState.models["sop"].keyVarDict["bypassNames"]
    bypassSubsystem = dict(zip(bypassNames, bypassed))
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
        if bypassSubsystem.get("ffs", False):
            cmd.warn('text="%s; guidingIsOk failed, but ffs is bypassed in sop"' % msg)
        else:
            cmd.warn('text="%s; aborting guiding"' % msg)
            return False

#   should we allow guiding with lamps on if axes are disabled
#   check if lamps are actually ON
    ffLamp = actorState.models["mcp"].keyVarDict["ffLamp"]
    hgCdLamp = actorState.models["mcp"].keyVarDict["hgCdLamp"]
    neLamp = actorState.models["mcp"].keyVarDict["neLamp"]
    if (any(ffLamp) and not bypassSubsystem.get('ff_lamp', False)) or \
            (any(hgCdLamp) and not bypassSubstem.get('hgcd_lamp', False)) or \
            (any(neLamp) and not bypassSubsystem.get('ne_lamp', False)):
        cmd.warn('text="Calibration lamp on; aborting guiding"')
        return False

#   check if non sensed lamps are commanded ON
    uvLamp = actorState.models["mcp"].keyVarDict["uvLampCommandedOn"]
    whtLamp = actorState.models["mcp"].keyVarDict["whtLampCommandedOn"]
    if uvLamp.getValue() or whtLamp.getValue():
        cmd.warn('text="Calibration lamp commanded on; aborting guiding"')
        return False
    
    tccState = actorState.tccState
    if tccState.halted or tccState.goToNewField:
        if bypassSubsystem.get("axes", False):
            cmd.warn('text="TCC motion failed, but axis motions are bypassed in sop"')
        else:
            cmd.warn('text="TCC motion aborted guiding"')
            return False

    return True
