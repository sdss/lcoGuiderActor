import ConfigParser
import Queue, threading
import math, numpy, re
import time
import subprocess

from guiderActor import *
import loadGprobes
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
import opscore.utility.tback as tback
import opscore.utility.YPF as YPF
import RO

import PID

from gimg.guiderImage import GuiderImageAnalysis
import loadGprobes

import pyfits
import os.path
import scipy.interpolate

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
        self.centerUp = False
        self.guideCmd = None

        self.fscanMJD = self.fscanID = -1
        self.design_ha = numpy.nan
        self.deleteAllGprobes()

        self.plugPlateScale = numpy.nan
        self.dSecondary_dmm = numpy.nan
        self.gcameraPixelSize = numpy.nan
        self.gcameraMagnification = numpy.nan
        
        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")
        self.setRefractionBalance(0.0)
        
        self.pid = {}               # PIDs for various axes
        for what in ["raDec", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

        self.decenter = False                #gstate only
        self.setDecenter("decenterRA")       
        self.setDecenter("decenterDec")      
        self.setDecenter("decenterRot")
        self.decenterChanged = True
        self.decenterFocus = numpy.nan       
        self.decenterScale = numpy.nan

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

    def setRefractionBalance(self, value=0.0):
        self.refractionBalance = value
        
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

    def setScales(self, plugPlateScale=None,
                  dSecondary_dmm=None,
                  gcameraPixelSize=None,
                  gcameraMagnification=None):

        if plugPlateScale != None:
            self.plugPlateScale = plugPlateScale
        if dSecondary_dmm != None:
            self.dSecondary_dmm = dSecondary_dmm
        if gcameraPixelSize != None:
            self.gcameraPixelSize = gcameraPixelSize
        if gcameraMagnification != None:
            self.gcameraMagnification = gcameraMagnification
            
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

        self.guideRMS = numpy.nan
        self.nguideRMS = numpy.nan
        self.guideXRMS = numpy.nan
        self.guideYRMS = numpy.nan

        self.guideAzRMS = numpy.nan      #not implemented yet
        self.guideAltRMS = numpy.nan     

        self.guideFitRMS = numpy.nan     #not implemented yet
        self.nguideFitRMS = numpy.nan
        self.nrejectFitRMS = numpy.nan

        self.decenterRA = numpy.nan
        self.decenterDec = numpy.nan
        self.decenterRot = numpy.nan
        self.decenterFocus = numpy.nan
        self.decenterScale = numpy.nan

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
    def inform(self, text):
        self._respond('i', text)
    def diag(self, text):
        self._respond('d', text)
    def finish(self, text):
        self._respond(':', text)
    def fail(self, text):
        self._respond('f', text)

def setupGstate(cartFile, plateFile, cartridgeId):
    global gState

    gState = GuiderState()
    config = ConfigParser.ConfigParser()
    config.read('etc/guider.cfg')
    
    plugPlateScale = float(config.get('telescope', "scale"))
    dSecondary_dmm = float(config.get('telescope', "dSecondary_dmm"))
    gcameraPixelSize = float(config.get('gcamera', "pixelSize"))
    gcameraMagnification = float(config.get('gcamera', "magnification"))
    gState.setScales(plugPlateScale=plugPlateScale,
                     dSecondary_dmm=dSecondary_dmm,
                     gcameraMagnification=gcameraMagnification,
                     gcameraPixelSize=gcameraPixelSize)
        
    gState.setGuideMode('axes', False)
    gState.setGuideMode('focus', False)
    gState.setGuideMode('scale', False)
    gState.setRefractionBalance(0.0)
    gState.gprobes = loadGprobes.getGprobes(cartFile, plateFile, cartridgeId)
    gState.cartridge = cartridgeId
    cmd = FakeCommand()
    gState.setCmd(cmd)

def processOneFile(guiderFile):
    queues = dict(MASTER=Queue.Queue())

    guideStep(None, queues, gState.guideCmd, guiderFile, True)

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
    minStarFlux = 500 #ADU, avoid guiding on noise spikes during acquisitions
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
        guideCmd.fail('guideState="failed"; text=%s' % qstr("No flat image available"))
        gState.setCmd(None)
        return
    
    if not darkfile:
        guideCmd.fail('guideState="failed"; text=%s' % qstr("No dark image available"))
        gState.setCmd(None)
        return

    if flatcart != gState.cartridge:
        if False:
            guideCmd.fail('guideState="failed"; text=%s' % qstr("Guider flat is for cartridge %d but %d is loaded" % (
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
        guideCmd.fail('guideState="failed"; text=%s' % qstr("Error in processing guide images: %s" % e))
        gState.setCmd(None)
        tback.tback("GuideTest", e)
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
    arcsecPerMM = 3600./gState.plugPlateScale   #arcsec per mm

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
    inFocusFwhm = []

    # Grab some times for refraction correction
    longitude = -105.82045
    UTC = RO.Astro.Tm.utcFromPySec(time.time() +
                                   actorState.models["tcc"].keyVarDict["utc_TAI"][0])
    LST = RO.Astro.Tm.lastFromUT1(UTC, longitude)
    
    RAkey = actorState.models["tcc"].keyVarDict["objNetPos"][0]
    RA = RAkey.getPos()
    HA = LST-RA     # The corrections are indexed by degrees, happily.
    dHA = HA - gState.design_ha
    haLimWarn = False
    guideCmd.diag('text="LST=%0.4f RA=%0.4f HA=%0.4f desHA=%0.4f dHA=%0.4f"' %
                  (LST, RA, HA, gState.design_ha,dHA))

    # Manually set for now. CPL
    wavelength = 16600
    for fiber in fibers:
        # necessary?
        if fiber.gprobe is None:
            guideCmd.warn('text="Gprobe %d was not listed in plugmap info"' % fiber.fiberid)
            continue
        gp = fiber.gprobe
        probe = gp.info
        enabled = gp.enabled
        tooFaint = False

        # Center up on acquisition fibers only.
        if gState.centerUp and probe.fiber_type != "ACQUIRE":
            enabled = False

        if not enabled:
            guideCmd.diag('text="Gprobe %d is not enabled"' % fiber.fiberid)
            continue
            
        #
        # dx, dy are the offsets on the ALTA guider image
        #
        fiber.dx = guideCameraScale*(fiber.xs - fiber.xcen) + (probe.xFerruleOffset / 1000.)
        fiber.dy = guideCameraScale*(fiber.ys - fiber.ycen) + (probe.yFerruleOffset / 1000.)
        poserr = fiber.xyserr

        #
        # theta is the angle to rotate (x, y) on the ALTA to (ra, alt)
        #
        # phi is the orientation of the alignment hole measured clockwise from N
        # rotation is the anticlockwise rotation from x on the ALTA to the pin
        #
        theta = 90                   # allow for 90 deg rot of camera view, should be -90 
        theta += probe.rotation # allow for intrinsic fibre rotation
        try:
            theta -= probe.phi      # allow for orientation of alignment hole
        except Exception, e:
            cmd.warn('text="skipping phi-less probe %s"' % (fiber.fiberid))
            continue
        
        probe.rotStar2Sky = theta # Squirrel the real angle away.

        isnan = numpy.isnan
        #FIXME PH -- We should ignore gprobes not present on plate/pointing (MARVELS dual pointing)
        #               and ignore fibers not found in flat.
        #            However we probably want to record values of disabled fibers for diagnosis 
        if isnan(fiber.dx) or isnan(fiber.dy) or isnan(poserr):
            guideCmd.warn("text=%s" %
                          qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber measured=(%g, %g), nominal=(%g,%g)" % (
                              fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, probe.xCenter, probe.yCenter)))
            continue

        if fiber.flux < minStarFlux and enabled:
            guideCmd.warn("text=%s" %
                          qstr("Star in gprobe %d too faint for guiding flux %g < %g minimum flux" % (
                              fiber.fiberid, fiber.flux, minStarFlux)))
            tooFaint = True      #PH should we add an extra bit for this.

        if poserr == 0:
            guideCmd.warn("text=%s" %
                          qstr("position error is 0 for gprobe %d star=(%g, %g) fiber=(%g, %g) nominal=(%g,%g)" % (
                              fiber.fiberid, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen, probe.xCenter, probe.yCenter)))
            continue

        theta = math.radians(theta)
        ct, st = math.cos(theta), math.sin(theta)
        # error in guide star position; n.b. still in mm here
        dRA   =  fiber.dx*ct + fiber.dy*st
        dDec  = -fiber.dx*st + fiber.dy*ct
        dDec *= -1

        #FIXME PH -- calc dAlt and dAz for guiding diagnostics,(output as part of fiber?)
        
        # Apply refraction correction
        xOffset = 0.0
        yOffset = 0.0
        haTime = 0.0
        try:
            if wavelength in probe.haOffsetTimes:
                haTimes = probe.haOffsetTimes[wavelength]
                if dHA < haTimes[0]:
                    if not haLimWarn:
                        cmd.warn('text="dHA (%0.1f) is below interpolation table; using limit (%0.1f)"' % (dHA, haTimes[0]))
                        haLimWarn = True
                    haTime = haTimes[0]
                elif dHA > haTimes[-1]:
                    if not haLimWarn:
                        cmd.warn('text="dHA (%0.1f) is above interpolation table; using limit (%0.1f)"' % (dHA, haTimes[-1]))
                        haLimWarn = True
                    haTime = haTimes[-1]
                else:
                    haTime = dHA

                # I'm now assuming 0...offset, but it should be offset1...offset2
                xInterp = scipy.interpolate.interp1d(haTimes,
                                                     probe.haXOffsets[wavelength])
                xOffset = gState.refractionBalance * xInterp(haTime)
                yInterp = scipy.interpolate.interp1d(haTimes,
                                                     probe.haYOffsets[wavelength])
                yOffset = gState.refractionBalance * yInterp(haTime)
        except Exception, e:
            guideCmd.diag('text="failed to calc refraction offsets for %s: %s"' % (wavelength, e))
            pass

        guideCmd.inform('refractionOffset=%d,%d,%0.1f,%0.4f,%0.6f,%0.6f' % (frameNo, fiber.fiberid,
                                                                            gState.refractionBalance,
                                                                            haTime,
                                                                            xOffset*arcsecPerMM,
                                                                            yOffset*arcsecPerMM))
        dRA -= xOffset
        dDec -= yOffset
        
        # Apply RA & Dec user guiding offsets to mimic different xy fibers centers
        # The guiderRMS will be calculated around the new effective fiber centers
      
        if gState.decenter: 
            #convert decenter offset to mm on guider
            frameInfo.decenterRA  = gState.decenterRA/arcsecPerMM
            frameInfo.decenterDec = gState.decenterDec/arcsecPerMM
            frameInfo.decenterRot = gState.decenterRot            
            # apply decenter offset so that telescope moves (not the star)
            dRA  += frameInfo.decenterRA
            dDec += frameInfo.decenterDec
            #decenterRot applied after guide solution
        else:
            frameInfo.decenterRA = 0.0
            frameInfo.decenterDec = 0.0
            frameInfo.decenterRot = 0.0
            frameInfo.decenterFocus = 0.0
            frameInfo.decenterScale = 0.0

        #output the keywords only when the decenter changes
        if gState.decenterChanged: 
            guideCmd.inform("decenter=%d, %s, %7.2f, %7.2f, %7.2f, %7.2f, %7.2f" % (
            frameNo, ("enabled" if gState.decenter else "disabled"), frameInfo.decenterRA, frameInfo.decenterDec,
            frameInfo.decenterRot, frameInfo.decenterFocus, frameInfo.decenterScale))

            gState.decenterChanged = False

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

        refmag = numpy.nan
        guideCmd.inform("probe=%d,%2d,0x%02d, %7.2f,%7.2f, %7.3f,%4.0f, %7.2f,%6.2f,%6.2f, %7.2f,%6.2f" % (
            frameNo, fiber.fiberid, probe.flags,
            fiber.dRA*arcsecPerMM, fiber.dDec*arcsecPerMM,
            fiber.fwhm, probe.focusOffset,
            fiber.flux, fiber.mag, refmag, fiber.sky, fiber.skymag))
                
        print "%d %2d  %7.2f %7.2f  %7.2f %7.2f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %06.1f  %7.3f %7.3f %7.0f %7.2f %4.0f" % (
            frameNo,
            fiber.fiberid, dRA, dDec, fiber.dx, fiber.dy, fiber.xs, fiber.ys, fiber.xcen, fiber.ycen,
            probe.xFocal, probe.yFocal, probe.rotStar2Sky, fiber.fwhm/sigmaToFWHM, fiber.sky, fiber.flux, fiber.mag,
            probe.focusOffset)

        if not enabled or tooFaint:
            continue

        #Collect fwhms for good in focus stars
        #Allow for a possible small range of focus offsets
        if abs(probe.focusOffset) < 50 : inFocusFwhm.append(fiber.fwhm)

        #accumulate guiding errors for good stars used in fit
        guideRMS += fiber.dx**2 + fiber.dy**2
        guideXRMS += fiber.dx**2
        guideYRMS += fiber.dy**2        
        nguideRMS += 1
        #guideAzRMS += fiber.dAz**2
        #guideAltRMS += fiber.dAlt**2        


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

        GI.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)

        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.setCmd(None)
            return

        #if guidingIsOK(cmd, actorState):
        #    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
        #                            expTime=gState.expTime))
        return
        
    A[2, 0] = A[0, 2]
    A[2, 1] = A[1, 2]
    try:
        if nStar == 1:
            guideCmd.warn('text="Only one star is usable"')
            x = b
            x[2, 0] = 0 # no rotation
        else:
            x = numpy.linalg.solve(A, b)

        # convert from mm to degrees
        dRA = x[0, 0]/gState.plugPlateScale
        dDec = x[1, 0]/gState.plugPlateScale
        dRot = -math.degrees(x[2, 0]) # and from radians to degrees

        #add the decenter guiding rotation offset here for now (in degrees)
        if gState.decenter:
            dRot += frameInfo.decenterRot/3600.0

        frameInfo.dRA  = dRA
        frameInfo.dDec = dDec
        frameInfo.dRot = dRot
        print 'dRA,dDec,dRot', dRA, dDec, dRot

        if gState.centerUp:
            offsetRa = -dRA
            offsetDec = -dDec
            offsetRot = 0
        else:
            offsetRa  = -gState.pid["raDec"].update(dRA)                    
            offsetDec = -gState.pid["raDec"].update(dDec)
            offsetRot = -gState.pid["rot"].update(dRot) if nStar > 1 else 0 # don't update I

        print 'offsetRA, offsetDec, offsetRot', offsetRa, offsetDec, offsetRot

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
        if True:
            nFwhm = len(inFocusFwhm)
            trimLo = 1 if nFwhm > 4 else 0
            trimHi= nFwhm - trimLo
            nKept = nFwhm - 2*trimLo
            nReject = nFwhm - nKept
            meanFwhm = (sum(inFocusFwhm))/nFwhm if nFwhm>0 else numpy.nan
            tMeanFwhm = (sum(sorted(inFocusFwhm)[trimLo:trimHi]))/nKept if nKept>0 else numpy.nan
            #loKept = inFocusFwhm[trimLo]
            #hiKept = inFocusFwhm[(trimHi-1)]
            print ("FWHM: %d, %7.2f, %d, %d, %7.2f" % (
                    frameNo, tMeanFwhm, nKept, nReject, meanFwhm))

            guideCmd.inform('fwhm=%d, %7.2f, %d, %d, %7.2f' % (
                    frameNo, tMeanFwhm, nKept, nReject, meanFwhm))

            frameInfo.meanFwhm = meanFwhm
            frameInfo.tMeanFwhm = tMeanFwhm

        #rms position error prior to this frames correction
        try:
            guideRMS  = (math.sqrt(guideRMS/nguideRMS)) *arcsecPerMM
            guideXRMS = (math.sqrt(guideXRMS/nguideRMS))*arcsecPerMM
            guideYRMS = (math.sqrt(guideYRMS/nguideRMS))*arcsecPerMM
        except:
            guideRMS = numpy.nan; guideXRMS = numpy.nan; guideYRMS = numpy.nan

        frameInfo.guideRMS  = guideRMS
        frameInfo.nguideRMS = nguideRMS
        frameInfo.guideXRMS = guideXRMS
        frameInfo.guideYRMS = guideYRMS

        #FIXME PH ---Need to calculate Az and Alt RMS in arcsec
        guideAzRMS  = numpy.nan
        guideAltRMS = numpy.nan
        #frameInfo.guideAzRMS = guideAzRMS
        #frameInfo.guideAltRMS = guideAltRMS     
     
        if gState.guideAxes:
            offsetsOK = True
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="offset arc %f, %f" % \
                                         (-offsetRa, -offsetDec))
            if cmdVar.didFail:
                offsetsOK = False
                guideCmd.warn('text="Failed to issue offset"')

            if offsetRot: 
                cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                         cmdStr="offset guide %f, %f, %g" % \
                                             (0.0, 0.0, -offsetRot))
            if cmdVar.didFail:
                offsetsOK = False
                guideCmd.warn('text="Failed to issue offset in rotator"')

        elif gState.centerUp:
            # If we are in the middle of an fk5InFiber (or other TCC track/pterr),
            # adjust the calibration offsets
            doCalibOffset = actorState.models["tcc"].keyVarDict["objName"][0] == "position reference star"
            if doCalibOffset:
                guideCmd.warn('text="using arc offsets at a pointing star"')
                
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="offset arc %f, %f" % \
                                     (-offsetRa, -offsetDec))
            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue centering offset"')
            else:
                if not doCalibOffset:
                    gState.setGuideMode('axes', True)
            
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

    if nStar <= 1 or gState.centerUp:      # don't bother with focus/scale!
        GI.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)

        if oneExposure:
            queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            gState.setCmd(None)

        return

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
    # So for now defer focus changes if we apply a scale change.
    blockFocusMove = False
        
    if gState.guideScale and abs(offsetScale) > 1e-7:
        # Clip to the motion we think is too big to apply at once.
        offsetScale = 1 + max(min(offsetScale, 2e-6), -2e-6)
        offsetScale *= curScale
        # cmd.warn('text="setting scale=%0.8f"' % (offsetScale))

        # Last chance to bailout.
        if offsetScale < 0.9995 or offsetScale > 1.0005:
            cmd.warn('text="NOT setting scarily large scale=%0.8f"' % (offsetScale))
        else:
            # blockFocusMove = True
            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                     cmdStr="set scale=%.9f" % (offsetScale))
            if cmdVar.didFail:
                guideCmd.warn('text="Failed to issue scale change"')

    #Evaluate RMS on fit over fibers used in fits here
    #FIXME--PH not calculated yet
    guideFitRMS = numpy.nan
    nguideFitRMS = 0
    nguideRejectFitRMS = 0

    # RMS guiding error output has to be after scale estimation so the full fit residual can be reported

    print "RMS guiding error= %4.3f, n stars= %d RMS_Az= %4.3f, RMS_Alt=%4.3f, RMS_X= %4.3f, RMS_Y=%4.3f" %(
        guideRMS, nguideRMS, guideAzRMS, guideAltRMS, guideXRMS, guideYRMS) 
    guideCmd.inform("guideRMS=%5d,%4.3f,%4d,%4.3f,%4.3f,%4.3f,%4.3f,%4.3f,%4d,%4d" % (
        frameNo, guideRMS, nguideRMS, guideAzRMS, guideAltRMS, 
        guideXRMS, guideYRMS, guideFitRMS, nguideFitRMS, nguideRejectFitRMS))

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
                                # PH -- currently all valid enabled fibers are used so OK.
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
        guideCmd.respond("focusError=%g" % (numpy.nan))
        guideCmd.respond("focusChange=%g, %s" % (numpy.nan, "enabled" if (gState.guideFocus and not blockFocusMove) else "disabled"))
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

def loadAllProbes(cmd, gState):
    gState.allProbes = None
    try:
        path = os.path.join(os.environ['PLATEDB_DIR'],
                            'bin', 'catPlPlugMapM')
        cmd1 = "%s -c %s -m %s -p %s -f %s %s" % (path,
                                                  gState.cartridge, gState.fscanMJD,
                                                  gState.pointing, gState.fscanID,
                                                  gState.plate)
        try:
            cmd.diag('text=%s' % (qstr('running: %s' % (cmd1))))
            ret = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
            plugmapBlob, errText = ret.communicate()
        except subprocess.CalledProcessError, e:
            cmd.warn('text="failed to load plugmap file: %s"' % (e))
            return

        ypm = YPF.YPF(fromString=plugmapBlob)
        pm = ypm.structs['PLUGMAPOBJ'].asArray()


        keep = pm[numpy.where(((pm.holeType == "GUIDE") & (pm.objType == "NA"))
                              | (pm.holeType == "OBJECT"))]
        cmd.diag('text="kept %d probes"' % (len(keep)))
        gState.allProbes = keep
    except Exception, e:
        cmd.warn('text=%s' % (qstr("could not load all probe info: %s" % (e))))
    
def loadTccBlock(cmd, actorState, gState):
    try:
        cmd1 = "catPlPlugMapM -c %s -m %s -p %s -f %s %s" % (gState.cartridge, gState.fscanMJD,
                                                             gState.pointing, gState.fscanID,
                                                             gState.plate)
        scratchFile = '/tmp/v_ca1_%04d.dat' % (gState.plate)
        cmd1 += " | %s /dev/stdin > %s" % (os.path.join(os.environ['GUIDERACTOR_DIR'],
                                                        'python/guiderActor',
                                                        'convertPlPlugMap.py'),
                                           scratchFile)
        cmd.diag('text=%s' % (qstr('running: %s' % (cmd1))))
        
        cmd2 = " %s %s" % (os.path.join(os.environ['GUIDERACTOR_DIR'],
                                        'python/guiderActor',
                                        'xferBlock.py'),
                           scratchFile)
        cmd.diag('text=%s' % (qstr('running: %s' % (cmd2))))
        ret = subprocess.call(cmd1, shell=True)
        if ret < 0:
            raise RuntimeError("cat and convert job failed with %s" % (-retcode))
        ret = subprocess.call(cmd2, shell=True)
        if ret < 0:
            raise RuntimeError("xfer job failed with %s" % (-retcode))

        cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                            cmdStr="set inst=spectro/gcview=%s/keep=(scaleFac)" % (gState.plate))
        if cmdVar.didFail:
            cmd.fail('text="Failed to set inst!"')
    except Exception, e:
        cmd.warn('text=%s' % (qstr("could not load a per-cartridge instrument block: %s" % (e))))
        
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

            qlen = queues[MASTER].qsize()
            if qlen > 0 and msg.cmd:
                msg.cmd.diag("master thread has %d items after a .get()" % (qlen))
                
            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))

                return

            elif msg.type == Msg.CENTERUP:
                # Arrange for the next exposure to do a centerUp.
                if not gState.guideCmd:
                    msg.cmd.fail('text="The guider must be running in order to centerUp"')
                    continue
                else:
                    gState.centerUp = msg.cmd # Provide some way for this command to be finished.

                continue

            elif msg.type == Msg.START_GUIDING:
                if not msg.start:
                    try:
                        success = msg.success
                    except AttributeError:
                        success = True
                        
                    if msg.start is None:
                        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                        continue
                        
                    if not gState.guideCmd:
                        msg.cmd.fail('text="The guider is already off"')
                        continue

                    if success:
                        msg.cmd.respond("guideState=stopping")
                        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                        if gState.guideCmd != msg.cmd:
                            msg.cmd.finish()

                        gState.guideCmd.finish("guideState=off")
                        gState.setCmd(None)

                    else:
                        queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                        msg.cmd.respond("guideState=failed")
                        if gState.guideCmd != msg.cmd:
                            msg.cmd.fail()

                        gState.guideCmd.fail("guideState=failed")
                        gState.setCmd(None)
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
                    queues[MASTER].put(Msg(Msg.START_GUIDING, gState.guideCmd, start=False, success=False))
                    continue

                guideStep(actor, queues, msg.cmd, msg.filename, oneExposure,
                          plot=plot, psPlot=psPlot, sm=sm)
                gState.inMotion = False
                if gState.centerUp:
                    gState.centerUp.finish('')
                    gState.centerUp = False

                    # Reset any PID I terms and smoothing filters
                    for key in gState.pid.keys():
                        gState.pid[key].reset()

                    # Stuff has changed; tell STUI.
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))
                    
                if not gState.guideCmd:    # something fatal happened in guideStep
                    continue

                #Check if offsets were changed during guiding without force
                if (gState.decenter and not force):
                    cmd.fail('text="Decentred guiding must use force."')
                    gState.setDecenter("decenterRA")  #reset all to 0
                    gState.setDecenter("decenterDec")
                    gState.setDecenter("decenterRot")
                    gState.decenter = False
                    gState.decenterChanged = True
                    queues[MASTER].put(Msg(Msg.START_GUIDING, gState.guideCmd, start=False))
                    continue

                #
                # Is there anything to indicate that we shouldn't be guiding?
                #
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
                except:
                    tback.tback("findFibers", e)
                    cmd.fail('text="findFibers failed -- it probably could not find any lit fibers near their expected positions: %s"' % (e))
                    continue
                
                try:
                    flatoutname = GI.getProcessedOutputName(msg.filename) 
                    dirname, filename = os.path.split(flatoutname)
                    cmd.inform('file=%s/,%s' % (dirname, filename))
                    cmd.finish('text="flat image processing done"')
                except Exception, e:
                    tback.tback("findFibers2", e)
                    cmd.fail('text="failed to save flat: %s"' % (e))

                continue
                    
            elif msg.type == Msg.FAIL:
                msg.cmd.fail('guideState="failed"; text="%s"' % msg.text);

            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllGprobes()
                gState.setRefractionBalance(0.0)

                gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
                gState.fscanMJD, gState.fscanID = msg.fscanMJD, msg.fscanID
                gState.boresight_ra, gState.boresight_dec = msg.boresight_ra, msg.boresight_dec
                gState.design_ha = msg.design_ha
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
					
                # Build and install an instrument block for this cartridge info
                loadTccBlock(msg.cmd, actorState, gState)
                loadAllProbes(msg.cmd, gState)

                if gState.cartridge > 0 and gState.cartridge < 10:
                    gState.setRefractionBalance(1.0)

                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                
            elif msg.type == Msg.SET_PID:
                gState.pid[msg.what].setPID(Kp=msg.Kp, Ti=msg.Ti, Td=msg.Td, Imax=msg.Imax, nfilt=msg.nfilt)

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_REFRACTION:
                gState.setRefractionBalance(msg.value)

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
                for what in ('axes', 'focus', 'scale'):
                    gState.setGuideMode(what, False)
                if msg.gprobe:
                    actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=msg.cmd,
                                                                  start=True, force=True))
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
                    msg.cmd.fail("test=\"no plate is loaded\"")
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
                                                    cmdStr="set scale=%.8f" % (newScale))
                if cmdVar.didFail:
                    gState.inMotion = False
                    cmd.fail('text="Failed to set scale"')
                else:
                    gState.pid['focus'].reset()
                    gState.pid['scale'].reset()
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

                for k in gState.pid.keys():
                    gState.pid[k].setPID(dt=(gState.expTime + 5)) # "+ 5" to allow for overhead

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

                #Allow decenters to be setup prior to guiding
            elif msg.type == Msg.DECENTER:
                        gState.setDecenter("decenterRA",msg.decenterRA)
                        gState.setDecenter("decenterDec",msg.decenterDec)
                        gState.setDecenter("decenterRot",msg.decenterRot)
                        gState.decenterChanged = True
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
                    cmd.respond("pid=%s,%g,%g,%g,%g,%d" % (w, 
                                                           gState.pid[w].Kp, gState.pid[w].Ti, gState.pid[w].Td,
                                                           gState.pid[w].Imax, gState.pid[w].nfilt))
                cmd.diag('text="guideCmd=%s"' % (qstr(gState.guideCmd)))
                if gState.refractionBalance != 0.0:
                    cmd.warn('refractionBalance=%0.1f' % (gState.refractionBalance))
                else:
                    cmd.respond('refractionBalance=%0.1f' % (gState.refractionBalance))
                cmd.diag('text="design_ha=%0.1f"' % (gState.design_ha))
                
                if msg.finish:
                    cmd.finish()
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag('text="%s alive"' % threadName)
        except Exception, e:
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            if gState.guideCmd:
                gState.guideCmd.warn('text="%s"' % errMsg)
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
            (any(hgCdLamp) and not bypassSubsystem.get('hgcd_lamp', False)) or \
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
