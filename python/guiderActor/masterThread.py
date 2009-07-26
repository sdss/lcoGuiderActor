import Queue

from guiderActor import *
import guiderActor.myGlobals

class GuiderState(object):
    """Save the state of the guider"""

    class Fiber(object):
        def __init__(self, id, enabled=True):
            self.id = id
            self.enabled = True

    def __init__(self):
        self.cartridge = -1
        self.plate = -1
        self.pointing = "?"
        self.expTime = 0

        self.deleteAllFibers()

        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")

    def deleteAllFibers(self):
        """Delete all fibers """
        self.fibers = {}

    def setFiberState(self, fiber, enabled=True, create=False):
        """Set a fiber's state"""

        #import pdb; pdb.set_trace()
        if self.fibers[fiber] == None and create:
            self.fibers[fiber] = GuiderState.Fiber(fiber, enabled)
        else:
            self.fibers[fiber].enabled = enabled

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
                                       dataname=msg.filename, cartridgeId=7)
                    fibers, stars = obj.runAllSteps()

                    dx, dy, n = 0, 0, 0

                    for i in range(len(fibers)):
                        fiber, star = fibers[i], stars[i]
                        assert fiber.fiberid == star.fiberid
                        
                        try:
                            if gState.fibers[fiber.fiberid] and not gState.fibers[fiber.fiberid].enabled:
                                continue
                        except IndexError, e:
                            guideCmd.warn("Fiber %d was not listed in plugmap info" % fiber.fiberid)
                            
                        dx = fiber.xcen - star.xs
                        dy = fiber.ycen - star.ys
                        n += 1

                    if n == 0:
                        guideCmd.fail("No fibers are available for guiding")
                        continue

                    dAz = dx/float(n)
                    dAlt = dy/float(n)
                    dRot = 1e-5

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

                    dScale = 1e-5
                    scalePID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetScale = scalePID["P"]*dScale

                    guideCmd.respond("scaleError=%g" % (dScale))
                    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                                             "enabled" if gState.guideScale else "disabled"))

                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, expTime=gState.expTime))
            elif msg.type == Msg.FAIL:
                msg.cmd.fail("text=\"%s\"" % msg.text);
            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllFibers()

                gState.cartridge, gState.plate, gState.pointing = msg.cartridge, msg.plate, msg.pointing
                #
                # Set the gState.fibers array; note that the fiber IDs may not be 0-indexed
                #
                # find the min and max fibre listed in msg.fiberList
                #
                min, max = None, None
                for id, enabled in msg.fiberList:
                    if min is None or id < min:
                        min = id
                    if max is None or id > max:
                        max = id

                if max is not None:
                    gState.fibers = (max + 1)*[None]
                    for id, enabled in msg.fiberList:
                        gState.setFiberState(id, enabled, create=True)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.SET_GUIDE_MODE:
                gState.setGuideMode(msg.what, msg.enabled)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))

            elif msg.type == Msg.ENABLE_FIBER:
                if gState.plate < 0:
                    msg.cmd.error("test=\"no plate is loaded\"")
                    continue
                
                gState.setFiberState(msg.fiber, msg.enable)
            elif msg.type == Msg.SET_TIME:
                gState.expTime = msg.expTime

                if msg.cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, msg.cmd, finish=True))
            elif msg.type == Msg.STATUS:
                msg.cmd.respond("cartridgeLoaded=%d, %d, %s" % (gState.cartridge, gState.plate, gState.pointing))

                fiberState = []
                for f in gState.fibers:
                    if f:
                        fiberState.append("\"(%d=%s)\"" % (f.id, f.enabled))

                msg.cmd.respond("fibers=%s" % ", ".join(fiberState))
                msg.cmd.respond("guideEnable=%s, %s, %s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                msg.cmd.respond("expTime=%g" % (gState.expTime))

                if msg.finish:
                    msg.cmd.finish()
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)
        except Queue.Empty:
            actor.bcast.diag("text=\"master alive\"")
