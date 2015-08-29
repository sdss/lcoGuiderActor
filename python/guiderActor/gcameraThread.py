"""
Thread for commanding both the gcamera and ecamera.

It's ok to handle them both in one thread, because we can't use them at the same
time, so they'll never conflict with each other.
"""


import Queue, threading

from guiderActor import Msg, GCAMERA
from guiderActor import myGlobals

def expose(cmd, actorState, replyQueue, expTime, stack=1, cartridge=None, expType='expose', camera='gcamera'):
    """Take an exposure with the e/gcamera, and succeed/fail as appropriate."""
    cmd.respond('text="starting %s exposure"'%camera)
    filenameKey = actorState.models[camera].keyVarDict["filename"]

    cmdStr="{0} time={1} stack={2}".format(expType, expTime, stack)
    if expType == "flat":
        cmdStr += " cartridge={0}".format(cartridge)
        responseMsg = Msg.FLAT_FINISHED
    elif expType == "dark":
        responseMsg = Msg.DARK_FINISHED
    else:
        responseMsg = Msg.EXPOSURE_FINISHED

    # Allow for readout time for each exposure in the stack, plus a bit extra.
    timeLim = stack*expTime + stack*5 + 15

    cmd.diag('text="{0} {1} with timeLim={2}"'.format(camera, cmdStr, timeLim))
    try:
        cmdVar = actorState.actor.cmdr.call(actor=camera, cmdStr=cmdStr,
                                            keyVars=[filenameKey], timeLim=timeLim, forUserCmd=cmd)
        cmd.diag('text="{0} {1} didFail={2}"'.format(camera, cmdStr, cmdVar.didFail))
    except Exception, e:
        cmd.warn('text="{0} {1} raised {2}"'.format(camera, cmdStr, e))
        return

    if cmdVar.didFail:
        cmd.warn('text="Failed to take {0} exposure"'.format(camera))
        if cmdVar.lastReply and "Timeout" in cmdVar.lastReply.keywords:
            cmd.warn('text="{0} expose command exceeded time limit: {1}."'.format(camera, timeLim))
        replyQueue.put(Msg(responseMsg, cmd=cmd, success=False))
        return

    filename = cmdVar.getLastKeyVarData(filenameKey)[0]
    replyQueue.put(Msg(responseMsg, cmd=cmd, filename=filename, camera=camera, success=True))

def main(actor, queues):
    """Main look for thread to talk to gcamera"""

    timeout = myGlobals.actorState.timeout

    while True:
        try:
            msg = queues[GCAMERA].get(timeout=timeout)
            qlen = queues[GCAMERA].qsize()
            if qlen > 0 and msg.cmd:
                msg.cmd.diag("gcamera thread has %d items after a .get()" % (qlen))
                
            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))

                return
            
            elif msg.type == Msg.EXPOSE:
                camera = getattr(msg,'camera','gcamera')
                expType = getattr(msg,'expType','expose')
                cartridge = getattr(msg,'cartridge',None)
                stack = getattr(msg,'stack',1)
                expose(msg.cmd, myGlobals.actorState, msg.replyQueue, msg.expTime, stack=stack,
                       cartridge=cartridge, expType=expType, camera=camera)

            elif msg.type == Msg.ABORT_EXPOSURE:
                if not msg.quiet:
                    msg.cmd.respond('text="Request to abort an exposure when none are in progress"')
                with queues[GCAMERA].mutex:
                    queues[GCAMERA].queue.clear()
                
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)

        except Queue.Empty:
            actor.bcast.diag('text="gcamera alive"')
        except Exception, e:
            actor.bcast.error('text="gcamera thread got unexpected exception: %s"' % (e))
            
