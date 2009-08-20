import Queue
import numpy

from guiderActor import *
import guiderActor.myGlobals
from opscore.utility.qstr import qstr

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

def main(actor, queues):
    """Main look for master thread"""
    timeout = guiderActor.myGlobals.actorState.timeout

    gState = GuiderState()
    guideCmd = None                     # the Cmd that started the guide loop
    
    while True:
        try:
            msg = queues[MASTER].get(timeout=timeout)

            if msg.type == Msg.START_GUIDING:
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
                    guideCmd.respond("processing=%s" % msg.filename)

                    import gcamera.pyGuide_test as pg
                    # plate == 3210
                    obj = pg.GuideTest(flatname="/data/gcam/55034/gimg-0331.fits",
                                       darkname="/data/gcam/55034/gimg-0205.fits",
                                       dataname=msg.filename, cartridgeId=gState.cartridge)
                    fibers, stars = obj.runAllSteps()

                    print "RHL; not using fibers info from GuideTest"; del fibers
                    
                    spiderInstAngKey = actorState.models["tcc"].keyVarDict["spiderInstAng"]
                    spiderInstAng = spiderInstAngKey[0].getPos()
                    if True:
                        print "RHL spiderInstAng =", spiderInstAng
                        if spiderInstAng is None:
                            spiderInstAng = 0

                    if spiderInstAng is None:
                        guideCmd.warn("text=%s" % qstr("spiderInstAng is None; are we tracking?"))
                    #
                    # Setup to solve for the axis and maybe scale offsets
                    #
                    # N.B. fiber.xFocal and fiber.yFocal are the offsets of the stars
                    # wrt the center of the plate in mm; fiber.xcen/star.xs are in pixels,
                    # so we need a scale for the guide camera.  Nominally the guide camera
                    # has the same scale as the plug plate itself, but maybe it doesn't,
                    # so we'll include a possible magnification
                    #
                    guideCameraScale = gState.gcameraMagnification*gState.gcameraPixelSize*1e-3 # mm/pixel
                    
                    A = numpy.matrix(np.zeros(4*4).reshape([4,4]))
                    b = numpy.matrix(np.zeros(4).reshape([4,1]))

                    for star in stars:
                        try:
                            fiber = gState.gprobes[star.starid]
                        except IndexError, e:
                            guideCmd.warn("Gprobe %d was not listed in plugmap info" % star.starid)
                            
                        if not fiber.enabled:
                            continue
                        
                        theta = (fiber.rotation + spiderInstAngle)*math.pi/180
                        
                        dx = guideCameraScale*(fiber.xcen - fiber.xferruleOffset - star.xs)
                        dy = guideCameraScale*(fiber.ycen - fiber.yferruleOffset - star.ys)

                        dAz =   dx*math.cos(theta) + dy*math.sin(theta) # n.b. still in mm here
                        dAlt = -dx*math.sin(theta) + dy*math.cos(theta)

                        b[0] += dAz
                        b[1] += dAlt
                        b[2] += fiber.x*dAlt - fiber.y*dAz
                        b[3] += fiber.x*dAz  + fiber.y*dAlt

                        A[0, 0] += 1
                        A[0, 1] += 0
                        A[0, 2] += -dAlt
                        A[0, 3] +=  dAz

                        A[1, 0] += 0
                        A[1, 1] += 1
                        A[1, 2] += dAz
                        A[1, 3] += dAlt

                        A[2, 0] += -dAlt
                        A[2, 1] += dAz
                        A[2, 2] += dAz*dAz + dAlt*dAlt
                        A[2, 3] += 0

                        A[3, 0] += dAz
                        A[3, 1] += dAlt
                        A[3, 2] += 0
                        A[3, 3] += dAz*dAz + dAlt*dAlt
                        
                    if A[0, 0] == 0:
                        guideCmd.fail("No fibers are available for guiding")
                        continue

                    x = numpy.linalg.solve(M, b)

                    dAz =  x[0, 0]*gState.plugPlateScale # convert from mm to degrees
                    dAlt = x[1, 0]*gState.plugPlateScale
                    dRot = x[2, 0]*180/math.pi # convert from radians to degrees
                    dScale = x[3, 0]

                    posPID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetAz = posPID["P"]*dAz
                    offsetAlt = posPID["P"]*dAlt
                    offsetRot = posPID["P"]*dRot

                    guideCmd.respond("axisError=%g, %g, %g" % (dAz, dAlt, dRot))
                    guideCmd.respond("axisChange=%g, %g, %g, %s" % (offsetAz, offsetAlt, offsetRot,
                                                                    "enabled" if gState.guideAxes else "disabled"))

                    dFocus = -1e-3
                    focusPID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetFocus = focusPID["P"]*dFocus

                    guideCmd.respond("focusError=%g" % (dFocus))
                    guideCmd.respond("focusChange=%g, %s" % (offsetFocus,
                                                             "enabled" if gState.guideFocus else "disabled"))

                    scalePID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetScale = scalePID["P"]*dScale

                    guideCmd.respond("scaleError=%g" % (dScale))
                    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                                             "enabled" if gState.guideScale else "disabled"))

                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, expTime=gState.expTime))
            elif msg.type == Msg.FAIL:
                msg.cmd.fail("text=\"%s\"" % msg.text);
            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllGprobes()

                gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
                #
                # Set the gState.gprobes array (actually a dictionary as we're not sure which fibre IDs are present)
                #
                gState.gprobes = {}
                for id, info in msg.gprobes.items():
                    if info.exists:
                        import pdb; pdb.set_trace()
                        gState.setGprobeState(id, enable=info.enabled, info=info, create=True)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_GUIDE_MODE:
                gState.setGuideMode(msg.what, msg.exists)
                #
                # Report the cartridge status
                #
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

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_TIME:
                gState.expTime = msg.expTime

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
                msg.cmd.respond("plateScales=%g, %g, %g" % (gState.plugPlateScale, \
                                                            gState.gcameraMagnification, gState.gcameraPixelSize))

                if msg.finish:
                    msg.cmd.finish()
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag("text=\"master alive\"")
