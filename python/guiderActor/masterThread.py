import Queue, threading
import math, numpy

from guiderActor import *
import guiderActor.myGlobals
from opscore.utility.qstr import qstr
from opscore.utility.tback import tback
import PID

class GuiderState(object):
    """Save the state of the guider"""

    class Gprobe(object):
        def __init__(self, id, info, enable=True, probeType="guide"):
            self.id = id
            self.info = info
            self.enabled = enable
            self.probeType = probeType

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
        for what in ["alt_az", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

    def deleteAllGprobes(self):
        """Delete all fibers """
        self.gprobes = {}

    def setGprobeState(self, fiber, enable=True, info=None, create=False):
        """Set a fiber's state"""

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

    timeout = guiderActor.myGlobals.actorState.timeout

    guideCmd = None                     # the Cmd that started the guide loop
    
    while True:
        try:
            msg = queues[MASTER].get(timeout=timeout)

            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform("text=\"Exiting thread %s\"" % (threading.current_thread().name))

                return
            
            elif msg.type == Msg.START_GUIDING:
                try:
                    expTime = msg.expTime
                    
                    if expTime >= 0 and gState.expTime != expTime:
                        gState.expTime = expTime
                        queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=False))
                except AttributeError:
                    pass

                if msg.start:
                    if gState.plate < 0:
                        queues[MASTER].put(Msg(Msg.FAIL, msg.cmd,
                                                text="Please tell me about your cartridge and try again"))
                        continue

                    guideCmd = msg.cmd
                    #
                    # Reset any PID I terms
                    #
                    for key in gstate.pid.keys():
                        gState.pid[k].reset()

                    guideCmd.respond("guideState=starting")
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, expTime=gState.expTime))
                else:
                    if not guideCmd:
                        msg.cmd.fail("text=\"The guider is already off\"")
                        continue

                    msg.cmd.respond("guideState=stopping")
                    quiet = True
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, msg.cmd, quiet=quiet, priority=Msg.MEDIUM))
                    guideCmd.finish("guideState=off")
                    guideCmd = None
            elif msg.type == Msg.EXPOSURE_FINISHED:
                if not msg.aborted:
                    if not guideCmd:    # exposure already finished
                        continue
                    
                    guideCmd.respond("processing=%s" % msg.filename)

                    import gcamera.pyGuide_test as pg
                    obj = pg.GuideTest(dataname=msg.filename, cartridgeId=gState.cartridge, gprobes=gState.gprobes,
                                       flatname="/data/gcam/55034/gimg-0331.fits",
                                       darkname="/data/gcam/55034/gimg-0205.fits"
                                       )

                    try:
                        #import pdb; pdb.set_trace()
                        stars = obj.runAllSteps()[1]
                    except Exception, e:
                        tback("GuideTest", e)
                        guideCmd.fail("text=%s" % qstr("Error in processing guide images: %s" % e))
                        continue

                    spiderInstAngKey = guiderActor.myGlobals.actorState.models["tcc"].keyVarDict["spiderInstAng"]
                    spiderInstAng = spiderInstAngKey[0].getPos()

                    if spiderInstAng is None:
                        msg = qstr("spiderInstAng is None; are we tracking?")

                        if True:
                            guideCmd.warn("text=%s" % msg)
                            print "RHL Setting spiderInstAng to 0.0"
                            spiderInstAng = 0.0
                        else:
                            guideCmd.fail("text=%s" % msg)
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
                            
                        if not fiber.enabled:
                            continue

                        theta = (fiber.info.rotation + fiber.info.phi + spiderInstAng)*math.pi/180
                        
                        dx = guideCameraScale*(star.xs - (fiber.info.xCenter - fiber.info.xFerruleOffset))
                        dy = guideCameraScale*(star.ys - (fiber.info.yCenter - fiber.info.yFerruleOffset))

                        if dx != dx or dy != dy:
                            guideCmd.warn("text=%s" %
                                          qstr("NaN in analysis for gprobe %d star=(%g, %g) fiber=(%g, %g)" % (
                                star.fiberid, star.xs, star.ys, fiber.info.xCenter, fiber.info.yCenter)))
                            continue
                        
                        dAz =   dx*math.cos(theta) + dy*math.sin(theta) # n.b. still in mm here
                        dAlt = -dx*math.sin(theta) + dy*math.cos(theta)

                        b[0] += dAz
                        b[1] += dAlt
                        b[2] += fiber.info.xCenter*dAlt - fiber.info.yCenter*dAz

                        A[0, 0] += 1
                        A[0, 1] += 0
                        A[0, 2] += -fiber.info.yCenter

                        A[1, 1] += 1
                        A[1, 2] += fiber.info.xCenter

                        A[2, 2] += fiber.info.xCenter*fiber.info.xCenter + fiber.info.yCenter*fiber.info.yCenter
                        #
                        # Now scale.  We don't actually solve for scale and axis updates
                        # simultanously, and we don't allow for the axis update when
                        # estimating the scale. 
                        #
                        b3 += fiber.info.xCenter*dAz + fiber.info.yCenter*dAlt
                        
                    if A[0, 0] == 0:
                        guideCmd.fail("No fibers are available for guiding")
                        continue

                    A[1, 0] = A[0, 1]
                    A[2, 0] = A[0, 2]
                    A[2, 1] = A[1, 2]
                    try:
                        x = numpy.linalg.solve(A, b)

                        # convert from mm to degrees
                        dAz =  x[0, 0]/gState.plugPlateScale*math.cos(gState.boresight_dec*math.pi/180)
                        dAlt = x[1, 0]/gState.plugPlateScale
                        dRot = x[2, 0]*180/math.pi # convert from radians to degrees

                        offsetAz = gState.pid["alt_az"].update(dAz)                    
                        offsetAlt = gState.pid["alt_az"].update(dAlt)
                        offsetRot = gState.pid["rot"].update(dRot)

                        guideCmd.respond("axisError=%g, %g, %g" % (3600*dAz, 3600*dAlt, 3600*dRot))
                        guideCmd.respond("axisChange=%g, %g, %g, %s" % (3600*offsetAz, 3600*offsetAlt, 3600*offsetRot,
                                                                        "enabled" if gState.guideAxes else "disabled"))

                        if gState.guideAxes:
                            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                     cmdStr="offset guide %f, %f, %f" % (offsetAz, offsetAlt, offsetRot))

                            if cmdVar.didFail:
                                guideCmd.warn("text=\"Failed to issue offset\"")
                    except numpy.linalg.LinAlgError:
                        guideCmd.warn("text=%s" % qstr("Unable to solve for axis offsets"))
                    #
                    # Scale
                    #
                    dScale = b3/A[2, 2]
                    offsetScale = gState.pid["scale"].update(dScale)

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
                        sigma = 1

                        ivar = 1/(2*math.pow(rms*sigma, 2)) # x's inverse variance

                        d = fiber.info.focusOffset
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
                        rms0 = math.sqrt(A[0, 0] - Delta*Delta)

                        dFocus = Delta*guideCameraScale*gState.dSecondary_dmm # mm to move the secondary
                        offsetFocus = gState.pid["focus"].update(dFocus)

                        guideCmd.respond("seeing=%g" % (rms0*2.35))
                        guideCmd.respond("focusError=%g" % (dFocus))
                        guideCmd.respond("focusChange=%g, %s" % (offsetFocus, "enabled" if gState.guideFocus else "disabled"))
                        if gState.guideFocus:
                            cmdVar = actor.cmdr.call(actor="tcc", forUserCmd=guideCmd,
                                                     cmdStr="set focus=%f/incremental" % (offsetFocus))

                            if cmdVar.didFail:
                                guideCmd.warn("text=\"Failed to issue focus offset\"")
                    except numpy.linalg.LinAlgError:
                        guideCmd.warn("text=%s" % qstr("Unable to solve for focus offset"))
                #
                # Write output fits file for TUI
                #
                

                #
                # Has the TCC done something to indicate that we should stop guiding?
                #
                tccState = guiderActor.myGlobals.actorState.tccState
                if tccState.halted or tccState.slewing:
                    guideCmd.warn("text=\"TCC motion aborted guiding\"")
                    print "TCC aborting", tccState.halted, tccState.slewing
                    queues[MASTER].put(Msg(Msg.START_GUIDING, guideCmd, start=False))
                    continue
                #
                # Start the next exposure
                #
                queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, expTime=gState.expTime))

            elif msg.type == Msg.FAIL:
                msg.cmd.fail("text=\"%s\"" % msg.text);
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
            actor.bcast.diag("text=\"master alive\"")
