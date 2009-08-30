import Queue, threading
import math, numpy

from guiderActor import *
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
from opscore.utility.tback import tback
import PID

import gcamera.pyGuide_test as pg
reload(pg)
import pyfits


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

    #import pdb; pdb.set_trace()
    threadName = "master"

    actorState = guiderActor.myGlobals.actorState
    timeout = actorState.timeout
    force = False                       # guide even if the petals are closed
    oneExposure = False                 # just take a single exposure

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
                    if not guideCmd:
                        msg.cmd.fail('text="The guider is already off"')
                        continue

                    msg.cmd.respond("guideState=stopping")
                    quiet = True
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=quiet, priority=Msg.MEDIUM))
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

                force = msg.force
                oneExposure = msg.oneExposure

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
                    
                    guideCmd.respond("processing=%s" % msg.filename)

                    h = pyfits.getheader(msg.filename)
                    flatfile = h.get('FLATFILE', None)
                    flatcart = h.get('FLATCART', None)
                    darkfile = h.get('DARKFILE', None)
                    if not flatfile:
                        guideCmd.fail("text=%s" % qstr("No flat image available"))
                        continue
                    if not darkfile:
                        guideCmd.fail("text=%s" % qstr("No dark image available"))
                        continue
                    if flatcart != gState.cartridge:
                        guideCmd.fail("text=%s" % qstr("No flat file for this cartridge"))
                        continue 
                       
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
                        continue

                    spiderInstAngKey = actorState.models["tcc"].keyVarDict["spiderInstAng"]
                    spiderInstAng = spiderInstAngKey[0].getPos()

                    tccPos = actorState.models["tcc"].keyVarDict["tccPos"]
                    tccAlt = tccPos[1]
                    guideCmd.warn('text="alt=%g"' % tccAlt);

                    if spiderInstAng is None:
                        txt = qstr("spiderInstAng is None; are we tracking?")

                        if True:
                            guideCmd.warn("text=%s" % txt)
                            print "RHL Setting spiderInstAng to 0.0"
                            spiderInstAng = 0.0
                        else:
                            guideCmd.fail("text=%s" % txt)
                            continue
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

                    for star in stars:
                        try:
                            fiber = gState.gprobes[star.fiberid]
                        except IndexError, e:
                            guideCmd.warn("Gprobe %d was not listed in plugmap info" % star.fiberid)
                            continue
                            
                        if not fiber.enabled:
                            continue

                        if False:
                            print "gprobe %d star=(%g, %g) fiber=(%g, %g)" % (
                                star.fiberid, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter)
                            
                        if True:
                            theta = math.radians(fiber.info.rotation + fiber.info.phi - (spiderInstAng + 90))
                        else:
                            theta = math.radians(fiber.info.rotation + fiber.info.phi - (spiderInstAng + 90))
                            print "%2d %6.1f %6.1f %6.1f %6.1f" % \
                                  (star.fiberid, math.degrees(theta), fiber.info.rotation, fiber.info.phi,
                                   (fiber.info.rotation - fiber.info.phi - 0*(180 if star.fiberid <= 8 else 0))%360
                                   )
                        
                        dx = guideCameraScale*(star.xs - (fiber.info.xCenter - fiber.info.xFerruleOffset))
                        dy = guideCameraScale*(star.ys - (fiber.info.yCenter - fiber.info.yFerruleOffset))

                        if dx != dx or dy != dy:
                            guideCmd.warn("text=%s" %
                                          qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber=(%g, %g)" % (
                                star.fiberid, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter)))
                            continue
                        
                        ct, st = math.cos(theta), math.sin(theta)
                        dAz =   dx*ct + dy*st # error in guide star position; n.b. still in mm here
                        dAlt = -dx*st + dy*ct

                        ct, st = math.cos(math.radians(spiderInstAng*math)), math.sin(math.radians(spiderInstAng))
                        azCenter =  fiber.info.xCenter*st - fiber.info.yCenter*ct # centre of hole
                        altCenter = fiber.info.xCenter*ct + fiber.info.yCenter*st

                        if True:
                            import re
                            print "%d %2d  %7.2f %7.2f  %7.2f %7.2f  %6.1f %6.1f  %6.1f %6.1f  %7.3f %4.0f" % (
                                int(re.search(r"([0-9]+)\.fits$", msg.filename).group(1)),
                                star.fiberid, dAz, dAlt, dx, dy, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter,
                                star.fwhm/2.35, fiber.info.focusOffset)

                        b[0] += dAz
                        b[1] += dAlt
                        #b[2] += fiber.info.xCenter*dAlt - fiber.info.yCenter*dAz
                        b[2] += azCenter*dAlt - altCenter*dAz

                        A[0, 0] += 1
                        A[0, 1] += 0
                        #A[0, 2] += -fiber.info.yCenter
                        A[0, 2] += -altCenter

                        A[1, 1] += 1
                        #A[1, 2] += fiber.info.xCenter
                        A[1, 2] += azCenter

                        #A[2, 2] += fiber.info.xCenter*fiber.info.xCenter + fiber.info.yCenter*fiber.info.yCenter
                        A[2, 2] += azCenter*azCenter + altCenter*altCenter
                        #
                        # Now scale.  We don't actually solve for scale and axis updates
                        # simultanously, and we don't allow for the axis update when
                        # estimating the scale. 
                        #
                        b3 += azCenter*dAz + altCenter*dAlt
                        
                    if A[0, 0] == 0:
                        guideCmd.warn('text="No stars are available for guiding"')

                        if guidingIsOK(msg.cmd, actorState, force=force):
                            queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER],
                                                    expTime=gState.expTime))

                        continue

                    A[1, 0] = A[0, 1]
                    A[2, 0] = A[0, 2]
                    A[2, 1] = A[1, 2]
                    try:
                        x = numpy.linalg.solve(A, b)

                        #print guideCameraScale, x, A, b
                        print "x = ", x[0,0], x[1, 0], x[2, 0]
                        # convert from mm to degrees
                        dAz = -x[0, 0]/gState.plugPlateScale/math.cos(math.radians(tccAlt))
                        dAlt = x[1, 0]/gState.plugPlateScale
                        dRot = -math.degrees(x[2, 0]) # and from radians to degrees

                        offsetAz =  -gState.pid["azAlt"].update(dAz)                    
                        offsetAlt = -gState.pid["azAlt"].update(dAlt)
                        offsetRot = -gState.pid["rot"].update(dRot)

                        if False:
                            guideCmd.warn('text="Setting rot update to 0"')
                            offsetRot = 0
                        if False:
                            guideCmd.warn('text="Setting az update to 0"')
                            offsetAz = 0

                        guideCmd.respond("axisError=%g, %g, %g" % (3600*dAz, 3600*dAlt, 3600*dRot))
                        guideCmd.respond("axisChange=%g, %g, %g, %s" % (3600*offsetAz, 3600*offsetAlt, 3600*offsetRot,
                                                                        "enabled" if gState.guideAxes else "disabled"))

                        if gState.guideAxes:
                            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                     cmdStr="offset guide %f, %f, %f" % \
                                                     (-offsetAz, -offsetAlt, -offsetRot))

                            if cmdVar.didFail:
                                guideCmd.warn('text="Failed to issue offset"')
                    except numpy.linalg.LinAlgError:
                        guideCmd.warn("text=%s" % qstr("Unable to solve for axis offsets"))
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
                else:
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, replyQueue=queues[MASTER], expTime=gState.expTime))

            elif msg.type == Msg.FAIL:
                msg.cmd.fail('text="%s"' % msg.text);
            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllGprobes()

                gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
                gState.boresight_ra, gState.boresight_dec = msg.boresight_ra, msg.boresight_dec
                #
                # Set the gState.gprobes array (actually a dictionary as we're not sure which fibre IDs are present)
                #
                gState.gprobes = {}
                for id, info in msg.gprobes.items():
                    if info.exists:
                        enabled = info.enabled if info.fiber_type == "GUIDE" else False
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
                msg.cmd.inform('text="The guider is %s"' % ("running" if guideCmd else "off"))

                msg.cmd.respond("cartridgeLoaded=%d, %d, %s" % (gState.cartridge, gState.plate, gState.pointing))

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
            #import pdb; pdb.set_trace()
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            actor.bcast.warn('text="%s"' % errMsg)
            tback(errMsg, e)

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
