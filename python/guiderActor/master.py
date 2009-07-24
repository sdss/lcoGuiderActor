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
        self.plate = -1
        self.pointing = None
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

def daemon(actor, queues):
    timeout = guiderActor.myGlobals.actorState.timeout

    gState = GuiderState()
    guideCmd = None                     # the Cmd that started the guide loop
    
    while True:
        try:
            msg = queues[MASTER].get(timeout=timeout)

            if msg.type == Msg.START_GUIDING:
                #import pdb; pdb.set_trace()
                cmd, start, expTime = msg.data
                if expTime >= 0:
                    gState.time = expTime

                if start:
                    if gState.plate < 0:
                        cmd.fail("text=\"Please tell me about your cartridge and try again\"")
                        continue

                    guideCmd = cmd

                    guideCmd.respond("guideState=starting")
                    queues[GCAMERA].put(Msg(Msg.EXPOSE, (cmd, gState.time)))
                else:
                    if not guideCmd:
                        cmd.fail("text=\"The guider is already off\"")
                        continue

                    cmd.respond("guideState=stopping")
                    quiet = True
                    queues[GCAMERA].put(Msg(Msg.ABORT_EXPOSURE, (cmd, quiet), priority=msg.MEDIUM))
                    guideCmd.finish("guideState=off")
                    guideCmd = None
            elif msg.type == Msg.EXPOSURE_FINISHED:
                filename, aborted = msg.data

                if not aborted:
                    guideCmd.respond("text=\"Processing %s\"" % filename)

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
                    dFocus = -1e-3
                    dScale = 1e-5
                    guideCmd.respond("offsets=%g, %g, %g, %s" % (dAz, dAlt, dRot,
                                                                 "enabled" if gState.guideAxes else "disabled"))
                    guideCmd.respond("focusChange=%g, %s" % (dFocus,
                                                             "enabled" if gState.guideFocus else "disabled"))
                    guideCmd.respond("scaleChange=%g, %s" % (dScale,
                                                             "enabled" if gState.guideScale else "disabled"))

                    queues[GCAMERA].put(Msg(Msg.EXPOSE, (guideCmd, gState.time)))
            elif msg.type == Msg.LOAD_CARTRIDGE:
                gState.deleteAllFibers()

                cmd, gState.plate, gState.pointing, fiberList = msg.data
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
                queues[MASTER].put(Msg(Msg.STATUS, (cmd, True)))

            elif msg.type == Msg.SET_GUIDE_MODE:
                cmd, what, enabled = msg.data
                
                gState.setGuideMode(what, enabled)
                #
                # Report the cartridge status
                #
                queues[MASTER].put(Msg(Msg.STATUS, (cmd, True)))

            elif msg.type == Msg.ENABLE_FIBER:
                cmd, fiber, enabled = msg.data

                if gState.plate < 0:
                    cmd.error("test=\"no plate is loaded\"")
                    continue
                
                gState.setFiberState(fiber, enabled)
            elif msg.type == Msg.SET_TIME:
                cmd, gState.time = msg.data

                if cmd:
                    queues[MASTER].put(Msg(Msg.STATUS, (cmd, True)))
            elif msg.type == Msg.STATUS:
                cmd, finish = msg.data

                cmd.respond("plate=%d pointing=%s" % (gState.plate, gState.pointing))

                fiberState = []
                for f in gState.fibers:
                    if f:
                        fiberState.append("\"(%d=%s)\"" % (f.id, f.enabled))

                cmd.respond("fibers=%s" % ", ".join(fiberState))
                cmd.respond("axes=%s focus=%s scale=%s" % (gState.guideAxes, gState.guideFocus, gState.guideScale))
                cmd.respond("time=%g" % (gState.time))

                if finish:
                    cmd.finish()
            else:
                raise ValueError, ("Unknown message type %d" % msg.type)
        except Queue.Empty:
            actor.bcast.diag("text=\"master alive\"")
