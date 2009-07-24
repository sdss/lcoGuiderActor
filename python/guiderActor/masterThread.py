import Queue, time

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
        self.time = 0

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
                #import pdb; pdb.set_trace()
                cmd, start, expTime = msg.cmd, msg.data["start"], msg.data.get("time")
                if expTime >= 0:
                    gState.time = expTime

                if start:
                    if gState.plate < 0:
                        cmd.fail("text=\"Please tell me about your cartridge and try again\"")
                        continue

                    guideCmd = cmd

                    guideCmd.respond("guideState=starting")
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, cmd, time=gState.time))
                else:
                    if not guideCmd:
                        cmd.fail("text=\"The guider is already off\"")
                        continue

                    cmd.respond("guideState=stopping")
                    quiet = True
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, cmd, quiet=quiet, priority=msg.MEDIUM))
                    guideCmd.finish("guideState=off")
                    guideCmd = None
            elif msg.type == Msg.EXPOSURE_FINISHED:
                filename, aborted = msg.data["filename"], msg.data["aborted"]

                if not aborted:
                    guideCmd.respond("processing=%s" % filename)

                    import gcamera.pyGuide_test as pg
                    # plate == 3210
                    obj = pg.GuideTest(flatname="/data/gcam/55034/gimg-0331.fits",
                                       darkname="/data/gcam/55034/gimg-0205.fits",
                                       dataname=filename, cartridgeId=7)
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

                    guideCmd.respond("axisErr=%g, %g, %g" % (dAz, dAlt, dRot))
                    guideCmd.respond("axisChange=%g, %g, %g, %s" % (offsetAz, offsetAlt, offsetRot,
                                                                    "enabled" if gState.guideAxes else "disabled"))

                    dFocus = -1e-3
                    focusPID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetFocus = posScale["P"]*scalePID

                    guideCmd.respond("focusErr=%g" % (dFocus))
                    guideCmd.respond("focusChange=%g, %s" % (offsetFocus,
                                                             "enabled" if gState.guideFocus else "disabled"))

                    dScale = 1e-5
                    scalePID = {"P" : 0.5, "I" : 0, "D" : 0}
                    offsetScale = posScale["P"]*scalePID

                    guideCmd.respond("scaleErr=%g" % (dScale))
                    guideCmd.respond("scaleChange=%g, %s" % (offsetScale,
                                                             "enabled" if gState.guideScale else "disabled"))

                    queues[GCAMERA].put(Msg(Msg.EXPOSE, guideCmd, time=gState.time))
            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllFibers()

                cmd = msg.cmd
                gState.cartridge, gState.plate, gState.pointing = \
                                  msg.data["cartridge"], msg.data["plate"], msg.data["pointing"]
                fiberList = msg.data["fiberList"]
                #
                # Set the gState.fibers array; note that the fiber IDs may not be 0-indexed
                #
                # find the min and max fibre listed in fiberList
                #
                min, max = None, None
                for id, enabled in fiberList:
                    if min is None or id < min:
                        min = id
                    if max is None or id > max:
                        max = id

                if max is not None:
                    gState.fibers = (max + 1)*[None]
                    for id, enabled in fiberList:
                        gState.setFiberState(id, enabled, create=True)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))

            elif msg.type == Msg.SET_GUIDE_MODE:
                cmd, what, enabled = msg.cmd, msg.data["what"], msg.data["enabled"]
                
                gState.setGuideMode(what, enabled)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))

            elif msg.type == Msg.ENABLE_FIBER:
                cmd, fiber, enabled = msg.cmd, msg.data["fiber"], msg.data["enable"]

                if gState.plate < 0:
                    cmd.error("test=\"no plate is loaded\"")
                    continue
                
                gState.setFiberState(fiber, enabled)
            elif msg.type == Msg.SET_TIME:
                cmd, gState.time = msg.cmd, msg.data["time"]

                if cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
            elif msg.type == Msg.STATUS:
                cmd, finish = msg.cmd, msg.data["finish"]

                cmd.respond("cartridgeLoaded=%d, %d, %s" % (gState.cartridge, gState.plate, gState.pointing))

                fiberState = []
                for f in gState.fibers:
                    if f:
                        fiberState.append("\"(%d=%s)\"" % (f.id, f.enabled))

                cmd.respond("fibers=%s" % ", ".join(fiberState))
                cmd.respond("guideEnable=%s, %s, %s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                cmd.respond("time=%g" % (gState.time))

                if finish:
                    cmd.finish()
            else:
                raise ValueError, ("Unknown message type %d" % msg.type)
        except Queue.Empty:
            actor.bcast.diag("text=\"master alive\"")
