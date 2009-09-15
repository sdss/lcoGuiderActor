import Queue, threading
import math, numpy, re

from guiderActor import *
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
from opscore.utility.tback import tback
import PID

import gcamera.pyGuide_test as pg
reload(pg)
import pyfits

try:
    import sm
except ImportError:
    print "Failed to import SM"
    sm = None

class GuiderState(object):
    """Save the state of the guider"""

    class Gprobe(object):
        def __init__(self, id, info, enable=True):
            self.id = id
            self.info = info
            self.enabled = enable

    def __init__(self):
        self.cartridge = -1
        self.plate = -1
        self.pointing = "?"
        self.expTime = 0

        self.fscanMJD = self.fscanID = -1

        self.deleteAllGprobes()

        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")

        self.pid = {}               # PIDs for various axes
        for what in ["azAlt", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

    def deleteAllGprobes(self):
        """Delete all fibers """
        self.gprobes = {}

    def setGprobeState(self, fiber, enable=True, info=None, create=False):
        """Set a fiber's state"""

        if fiber in ("ACQUIRE", "GUIDE"):
            fiber_type = fiber
            for gp in self.gprobes.values():
                if gp.info.fiber_type == fiber_type:
                    gp.enabled = enable
        else:
            if not self.gprobes.has_key(fiber) and create:
                self.gprobes[fiber] = GuiderState.Gprobe(fiber, info, enable)
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

try:
    gState
except:
    gState = GuiderState()

def main(actor, queues):
    """Main loop for master thread"""

    threadName = "master"

    actorState = guiderActor.myGlobals.actorState
    timeout = actorState.timeout
    guide_azAlt = False                 # guide in az/alt, not ra/dec
    force = False                       # guide even if the petals are closed
    oneExposure = False                 # just take a single exposure
    plot = False                        # use SM to plot things
    fakeSpiderInstAng = None            # the value we claim for the SpiderInstAng
    
    guideCmd = None                     # the Cmd that started the guide loop
    
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
                        
                    if not guideCmd:
                        msg.cmd.fail('text="The guider is already off"')
                        continue

                    msg.cmd.respond("guideState=stopping")
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=True, priority=Msg.MEDIUM))
                    guideCmd.finish("guideState=off")
                    guideCmd = None
                    continue

                try:
                    expTime = msg.expTime
                    
                    if expTime >= 0 and gState.expTime != expTime:
                        gState.expTime = expTime
                        queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))
                except AttributeError:
                    pass

                guide_azAlt = msg.azAlt
                force = msg.force
                oneExposure = msg.oneExposure
                plot = msg.plot
                fakeSpiderInstAng = msg.spiderInstAng
                
                if guideCmd:
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
                #
                # Reset any PID I terms
                #
                for key in gState.pid.keys():
                    gState.pid[key].reset()

                guideCmd.respond("guideState=starting")
                queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER], expTime=gState.expTime))

            elif msg.type == Msg.EXPOSURE_FINISHED:
                if True:
                    if not guideCmd:    # exposure already finished
                        continue

                    if not msg.success:
                        queues[MASTER].put(Msg(Msg.START_GUIDING, guideCmd, start=False))
                        continue
                    
                    guideCmd.respond("processing=%s" % msg.filename)

                    h = pyfits.getheader(msg.filename)
                    flatfile = h.get('FLATFILE', None)
                    flatcart = h.get('FLATCART', None)
                    darkfile = h.get('DARKFILE', None)
                    if not flatfile:
                        guideCmd.fail("text=%s" % qstr("No flat image available"))
                        guideCmd = None
                        continue
                    if not darkfile:
                        guideCmd.fail("text=%s" % qstr("No dark image available"))
                        guideCmd = None
                        continue
                    if flatcart != gState.cartridge:
                        if False:
                            guideCmd.fail("text=%s" % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                                flatcart, gState.cartridge)))
                            guideCmd = None
                            continue
                        else:
                            guideCmd.warn("text=%s" % qstr("Guider flat is for cartridge %d but %d is loaded" % (
                                flatcart, gState.cartridge)))
                       
                    imageObj = pg.GuideTest(dataname=msg.filename, cartridgeId=gState.cartridge, 
                                            gprobes=gState.gprobes,
                                            flatname=flatfile, darkname=darkfile, mode=1)

                    try:
                        if False:
                            imageObj.subDark()
                            correction = imageObj.addOffset()
                            imageObj.makeFinal(correction=correction)
                            fibers, stars = imageObj.finalData()
                        else:
                            fibers, stars = imageObj.runAllSteps()
                        for gcamFiber in fibers:
                            try:
                                fiber = gState.gprobes[gcamFiber.fiberid]
                            except IndexError, e:
                                guideCmd.warn("Gprobe %d was not listed in plugmap info" % star.fiberid)
                                continue

                            fiber.info.xCenter = gcamFiber.xcen
                            fiber.info.yCenter = gcamFiber.ycen
                    except Exception, e:
                        tback("GuideTest", e)
                        guideCmd.fail("text=%s" % qstr("Error in processing guide images: %s" % e))
                        guideCmd = None
                        continue

                    if guide_azAlt:
                        spiderInstAngKey = actorState.models["tcc"].keyVarDict["spiderInstAng"]
                        spiderInstAng = spiderInstAngKey[0].getPos()

                        tccPos = actorState.models["tcc"].keyVarDict["tccPos"]
                        tccAlt = tccPos[1]

                        if spiderInstAng is None or tccAlt is None:
                            txt = qstr("spiderInstAng/tccAlt is None; are we tracking?")

                            if force:
                                guideCmd.warn("text=%s" % txt)
                                tccAlt = 0.0

                                if fakeSpiderInstAng is None:
                                    guideCmd.fail('text=%s' %
                                                  qstr('You must specify a spiderInstAng to "guide" in az/alt when the telescope is not tracking"'))
                                    guideCmd = None
                                    continue
                                    
                                spiderInstAng = fakeSpiderInstAng
                            else:
                                guideCmd.fail("text=%s" % txt)
                                guideCmd = None
                                continue

                        #spiderInstAng += -90  #55077 worked
                        #spiderInstAng += 90    #ph 55081 worked
                        #spiderInstAng +=  90    #ph 55088 wrong
                        #spiderInstAng +=  0    #ph 55088
                        spiderInstAng = -spiderInstAng - 90 
                    #
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
                    
                    A = numpy.matrix(numpy.zeros(3*3).reshape([3,3]))
                    b = numpy.matrix(numpy.zeros(3).reshape([3,1])); b3 = 0.0

                    if plot:
                        if sm:
                            sm.device('x11')
                        else:
                            guideCmd.warn('text="Unable to plot as SM is not available"')

                        size = len(gState.gprobes) + 1 # fibers are 1-indexed
                        fiberid_np = numpy.zeros(size)
                        azCenter_np = numpy.zeros(size)
                        altCenter_np = numpy.zeros(size)
                        daz_np = numpy.zeros(size)
                        dalt_np = numpy.zeros(size)

                    for star in stars:
                        try:
                            fiber = gState.gprobes[star.fiberid]
                        except IndexError, e:
                            guideCmd.warn("Gprobe %d was not listed in plugmap info" % star.fiberid)
                            continue

                        try:
                            fiber.info.ra
                        except AttributeError, e:
                            guideCmd.warn("Gprobe %d was not listed in plPlugMapM file" % star.fiberid)
                            continue
                        #
                        # dx, dy are the offsets on the ALTA guider image
                        #
                        dx = guideCameraScale*(star.xs - (fiber.info.xCenter - fiber.info.xFerruleOffset))
                        dy = guideCameraScale*(star.ys - (fiber.info.yCenter - fiber.info.yFerruleOffset))

                        if dx != dx or dy != dy:
                            guideCmd.warn("text=%s" %
                                          qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber=(%g, %g)" % (
                                star.fiberid, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter)))
                            continue
                        #
                        # theta is the angle to rotate (x, y) on the ALTA to (az, alt)
                        #
                        # phi is the orientation of the alignment hole measured clockwise from N
                        # rotation is the anticlockwise rotation from x on the ALTA to the pin
                        #
                        if False:
                            theta = fiber.info.rotation - fiber.info.phi + (spiderInstAng if guide_azAlt else 0.0)
                        else:
                            print "RHL + - + -Sp   %d %.0f %0.f " % (star.fiberid,
                                                                 fiber.info.rotation, fiber.info.phi)

                            theta = 0
                            theta += fiber.info.rotation # allow for intrinsic fibre rotation
                            theta -= fiber.info.phi      # allow for orientation of alignment hole
                            if guide_azAlt:
                                theta += spiderInstAng # allow for telescope rotation

                        theta = math.radians(theta)
                        ct, st = math.cos(theta), math.sin(theta)
                        print "theta=", theta, "st=", st, "ct=", ct
                        dAz   =  dx*ct + dy*st # error in guide star position; n.b. still in mm here
                        dAlt  = -dx*st + dy*ct
                        dAz    = dAz  
                        dAlt   = -dAlt


                        if guide_azAlt:
                            ct, st = math.cos(math.radians(spiderInstAng)), math.sin(math.radians(spiderInstAng))
                            azCenter =   fiber.info.xFocal*ct + fiber.info.yFocal*st # centre of hole
                            altCenter = -fiber.info.xFocal*st + fiber.info.yFocal*ct
                        else:
                            azCenter =  fiber.info.xFocal
                            altCenter = fiber.info.yFocal

                        if plot:
                            try:
                                fiberid_np[star.fiberid] = star.fiberid
                                azCenter_np[star.fiberid] = azCenter
                                altCenter_np[star.fiberid] = altCenter
                                daz_np[star.fiberid] = dAz
                                dalt_np[star.fiberid] = dAlt
                            except IndexError, e:
                                #import pdb; pdb.set_trace()
                                pass

                        if True:
                            print "%d %2d  %7.2f %7.2f  %7.2f %7.2f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %7.3f %4.0f" % (
                                int(re.search(r"([0-9]+)\.fits$", msg.filename).group(1)),
                                star.fiberid, dAz, dAlt, dx, dy, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter,
                                fiber.info.xFocal, fiber.info.yFocal,
                                star.fwhm/2.35, fiber.info.focusOffset)


                        if not fiber.enabled:
                            continue

                        b[0] += dAz
                        b[1] += dAlt
                        b[2] += azCenter*dAlt - altCenter*dAz

                        A[0, 0] += 1
                        A[0, 1] += 0
                        A[0, 2] += -altCenter

                        A[1, 1] += 1
                        A[1, 2] += azCenter

                        A[2, 2] += azCenter*azCenter + altCenter*altCenter
                        #
                        # Now scale.  We don't actually solve for scale and axis updates
                        # simultanously, and we don't allow for the axis update when
                        # estimating the scale. 
                        #
                        b3 += azCenter*dAz + altCenter*dAlt
                        
                    nStar = A[0, 0]
                    if nStar == 0:
                        guideCmd.warn('text="No stars are available for guiding"')

                        imageObj.writeFITS(guideCmd)

                        if oneExposure:
                            queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                            guideCmd = None
                            continue

                        if guidingIsOK(msg.cmd, actorState, force=force):
                            queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
                                                    expTime=gState.expTime))
                            continue

                    A[1, 0] = A[0, 1]
                    A[2, 0] = A[0, 2]
                    A[2, 1] = A[1, 2]
                    try:
                        if nStar == 1:
                            guideCmd.warn('text="Only one star is usable; guiding on az/alt"')

                            x = b
                            x[2, 0] = 0 # no rotation
                        else:
                            x = numpy.linalg.solve(A, b)

                        # convert from mm to degrees
                        print "RHL -90 + + +"
                        dAz = x[0, 0]/gState.plugPlateScale
                        if guide_azAlt:
                            dAz /= math.cos(math.radians(tccAlt))
                        dAlt = x[1, 0]/gState.plugPlateScale
                        dRot = -math.degrees(x[2, 0]) # and from radians to degrees

                        offsetAz =  -gState.pid["azAlt"].update(dAz)                    
                        offsetAlt = -gState.pid["azAlt"].update(dAlt)
                        offsetRot = -gState.pid["rot"].update(dRot) if nStar > 1 else 0 # don't update I
#                        offsetRot = -offsetRot    #ph kluge 55081
                        offsetRot = offsetRot

                        dAzArcsec = dAz
                        offsetAzArcsec = offsetAz

                        if guide_azAlt:
                            dAzArcsec *= math.cos(math.radians(tccAlt)) # in degrees of arc
                            offsetAzArcsec *= math.cos(math.radians(tccAlt))

                        guideCmd.respond("axisError=%g, %g, %g" % (3600*dAzArcsec, 3600*dAlt, 3600*dRot))
                        guideCmd.respond("axisChange=%g, %g, %g, %s" % (3600*offsetAzArcsec,
                                                                        3600*offsetAlt, 3600*offsetRot,
                                                                        "enabled" if gState.guideAxes else "disabled"))

                        if gState.guideAxes:
                            if guide_azAlt:
                                guideCmd.inform('text="Guiding in az/alt"')
                                cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                         cmdStr="offset guide %f, %f, %f" % \
                                                         (-offsetAz, -offsetAlt, -offsetRot))
                                if cmdVar.didFail:
                                    guideCmd.warn('text="Failed to issue offset"')
                            else:
                                cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                         cmdStr="offset arc %f, %f" % \
                                                         (-offsetAz, -offsetAlt))

                                if cmdVar.didFail:
                                    guideCmd.warn('text="Failed to issue offset"')

                                if offsetRot:
                                    cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                             cmdStr="offset guide %f, %f, %g" % \
                                                             (0.0, 0.0, -offsetRot))

                                if cmdVar.didFail:
                                    guideCmd.warn('text="Failed to issue offset in rotator"')

                        if plot and sm:
                            sm.erase(False)
                            sm.limits([-400, 400], [-400, 400])
                            sm.box()
                            if guide_azAlt:
                                sm.xlabel(r"\2\delta Az")
                                sm.ylabel(r"\2\delta Alt")
                            else:
                                sm.xlabel(r"\2\delta Ra")
                                sm.ylabel(r"\2\delta Dec")
                            sm.ptype([63])
                            sm.frelocate(0.85, 0.95)
                            sm.putlabel(5, r"\1Frame " + re.search(r"([0-9]+)\.fits$", msg.filename).group(1))
#                            sm.toplabel("-SpIA-90 ")

                            vscale = 1000 # how much to multiply position error
                            sm.relocate(-350, 350)
                            asec = gState.plugPlateScale/3600.0
                            sm.draw(-350 + vscale*gState.plugPlateScale/3600.0, 350)
                            sm.label(r"  \raise-200{\1{1 arcsec}}")

                            for i in range(len(fiberid_np)):
                                if fiberid_np[i] == 0:
                                    continue

                                sm.relocate(azCenter_np[i], altCenter_np[i])
                                sm.putlabel(6, r" %d" % fiberid_np[i])

                                sm.relocate(azCenter_np[i], altCenter_np[i])

                                sm.dot()
                                
                                if not gState.gprobes[fiberid_np[i]].enabled:
                                    sm.ltype(1)
                                    sm.ctype(sm.CYAN)
                                    
                                sm.draw(azCenter_np[i] + vscale*daz_np[i],
                                        altCenter_np[i] + vscale*dalt_np[i])

                                sm.ltype()
                                sm.ctype()

                            sm.ptype()
                            
                    except numpy.linalg.LinAlgError:
                        guideCmd.warn("text=%s" % qstr("Unable to solve for axis offsets"))

                    if nStar <= 1:      # don't bother with focus/scale!
                        if oneExposure:
                            queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                            guideCmd = None
                            continue

                        if guidingIsOK(msg.cmd, actorState, force=force):
                            queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
                                                    expTime=gState.expTime))
                        continue
                    #
                    # Scale
                    #
                    dScale = b3/A[2, 2]
                    offsetScale = -gState.pid["scale"].update(dScale)

                    guideCmd.respond("scaleError=%g" % (dScale))
                    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                                             "enabled" if gState.guideScale else "disabled"))
                    #
                    # Now focus. If the ith star is d_i out of focus, and the RMS of an
                    # in-focus star would be r0, and we are Delta out of focus, we measure
                    # an RMS size R_i
                    #   R_i^2 = r0^2 + (d_i + Delta)^2
                    # i.e.
                    #   R_i^2 - d_i^2 = (r0^2 + Delta^2) + 2 Delta d_i
                    # which is a linear equation for x == R_i^2 - d_i^2
                    #
                    A = numpy.matrix(numpy.zeros(2*2).reshape([2,2]))
                    b = numpy.matrix(numpy.zeros(2).reshape([2,1]))

                    for star in stars:
                        try:
                            fiber = gState.gprobes[star.fiberid]
                        except IndexError, e:
                            guideCmd.warn("Gprobe %d was not listed in plugmap info" % star.fiberid)
                            continue
                            
                        if not fiber.enabled:
                            continue

                        try:
                            rms = star.rms
                        except AttributeError:
                            rms = star.fwhm/2.35

                        rms *= 1/3600.0*gState.plugPlateScale*1e3

                        sigma = 1

                        try:
                            ivar = 1/(2*math.pow(rms*sigma, 2)) # x's inverse variance
                        except ZeroDivisionError:
                            ivar = 1e-5

                        focalRatio = 5.0
                        d = fiber.info.focusOffset/focalRatio
                        x = rms*rms - d*d

                        b[0] += x*ivar
                        b[0] += x*d*ivar

                        A[0, 0] += ivar
                        A[0, 1] += d*ivar

                        A[1, 1] += d*d*ivar

                    A[1, 0] = A[0, 1]
                    try:
                        x = numpy.linalg.solve(A, b)
                    
                        Delta = x[1, 0]/2
                        try:
                            rms0 = math.sqrt(A[0, 0] - Delta*Delta)
                        except ValueError, e:
                            rms0 = float("NaN")

                        dFocus = Delta*gState.dSecondary_dmm # mm to move the secondary
                        offsetFocus = -gState.pid["focus"].update(dFocus)

                        guideCmd.respond("seeing=%g" % (rms0*2.35))
                        guideCmd.respond("focusError=%g" % (dFocus))
                        guideCmd.respond("focusChange=%g, %s" % (offsetFocus, "enabled" if gState.guideFocus else "disabled"))
                        if gState.guideFocus:
                            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                     cmdStr="set focus=%f/incremental" % (offsetFocus))

                            if cmdVar.didFail:
                                guideCmd.warn('text="Failed to issue focus offset"')
                    except numpy.linalg.LinAlgError:
                        guideCmd.warn("text=%s" % qstr("Unable to solve for focus offset"))
                #
                # Write output fits file for TUI
                #
                imageObj.writeFITS(guideCmd)

                #
                # Is there anything to indicate that we shouldn't be guiding?
                #
                if not guidingIsOK(msg.cmd, actorState, force=force):
                    queues[MASTER].put(Msg(Msg.START_GUIDING, guideCmd, start=False))
                    continue
                #
                # Start the next exposure
                #
                if oneExposure:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                    guideCmd = None
                else:
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER], expTime=gState.expTime))

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
                    if info.exists:
                        enabled = False if info.fiber_type == "TRITIUM" else info.enabled
                        gState.setGprobeState(id, enable=enabled, info=info, create=True)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
                
            elif msg.type == Msg.SET_PID:
                gState.pid[msg.what].setPID(Kp=msg.Kp, Ti=msg.Ti, Td=msg.Td, Imax=msg.Imax)

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
            elif msg.type == Msg.STATUS:
                msg.cmd.respond("cartridgeLoaded=%d, %d, %s, %d, %d" % (
                    gState.cartridge, gState.plate, gState.pointing, gState.fscanMJD, gState.fscanID))

                try:
                    if not msg.full:
                        if msg.finish:
                            msg.cmd.finish()
                        continue
                except AttributeError:
                    pass

                msg.cmd.inform('text="The guider is %s"' % ("running" if guideCmd else "off"))

                fiberState = []
                for f in gState.gprobes.values():
                    if f:
                        fiberState.append("\"(%d=%s)\"" % (f.id, f.enabled))

                if len(fiberState) > 0:
                    msg.cmd.respond("gprobes=%s" % ", ".join(fiberState))
                msg.cmd.respond("guideEnable=%s, %s, %s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                msg.cmd.respond("expTime=%g" % (gState.expTime))
                msg.cmd.respond("scales=%g, %g, %g, %g" % (gState.plugPlateScale,
                                                           gState.gcameraMagnification, gState.gcameraPixelSize,
                                                           gState.dSecondary_dmm,))
                for w in gState.pid.keys():
                    msg.cmd.respond("pid=%s, %g, %g, %g, %g" % (w,
                                                                gState.pid[w].Kp, gState.pid[w].Ti, gState.pid[w].Td,
                                                                gState.pid[w].Imax))
                if msg.finish:
                    msg.cmd.finish()
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag('text="%s alive"' % threadName)
        except Exception, e:
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            actor.bcast.warn('text="%s"' % errMsg)
            tback(errMsg, e)

            #import pdb; pdb.set_trace()
            try:
                print "\n".join(tback(errMsg, e)[0]) # old versions of tback return None
            except:
                pass

            guideCmd = False

            try:
                msg.replyQueue.put(Msg.EXIT, cmd=msg.cmd, success=False)
            except Exception, e:
                pass

def guidingIsOK(cmd, actorState, force=False):
    """Is it OK to be guiding?"""

    if force:
        return True

    ffsStatus = actorState.models["mcp"].keyVarDict["ffsStatus"]

    open, closed = 0, 0
    for s in ffsStatus:
        if s == None:
            cmd.warn('text="Failed to get state of flat field screen from MCP"')
            break

        open += int(s[0])
        closed += int(s[1])

    if open != 8:
        cmd.warn('text="FF petals aren\'t open; aborting guiding"')
        return False

    tccState = actorState.tccState
    if tccState.halted or tccState.goToNewField:
        cmd.warn('text="TCC motion aborted guiding"')
        print "TCC aborting", tccState.halted, tccState.slewing, tccState.goToNewField
        return False

    return True
