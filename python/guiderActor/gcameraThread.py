import Queue, threading

from guiderActor import *
import guiderActor.myGlobals

def main(actor, queues):
    """Main look for thread to talk to gcamera"""

    timeout = guiderActor.myGlobals.actorState.timeout

    while True:
        try:
            msg = queues[GCAMERA].get(timeout=timeout)
            
            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))

                return
            
            elif msg.type == Msg.EXPOSE:
                msg.cmd.respond('text="starting exposure"')
                #
                # Take exposure
                #
                timeLim = msg.expTime + 15     # allow for readout time

                filenameKey = guiderActor.myGlobals.actorState.models["gcamera"].keyVarDict["filename"]

                cmdVar = actor.cmdr.call(actor="gcamera", cmdStr="expose time=%f" % (msg.expTime),
                                         keyVars=[filenameKey], timeLim=timeLim)
                if cmdVar.didFail:
                    msg.cmd.warn('text="Failed to take exposure"')
                    msg.replyQueue.put(Msg(Msg.EXPOSURE_FINISHED, cmd=msg.cmd, success=False))
                    continue

                filename = cmdVar.getLastKeyVarData(filenameKey)[0]

                print "Sending EXPOSURE_FINISHED to", msg.replyQueue
                msg.replyQueue.put(Msg(Msg.EXPOSURE_FINISHED, cmd=msg.cmd, filename=filename, success=True))

            elif msg.type == Msg.ABORT_EXPOSURE:
                if not msg.quiet:
                    msg.cmd.respond('text="Request to abort an exposure when none are in progress"')
                guiderActor.flushQueue(queues[GCAMERA])
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)

        except Queue.Empty:
            actor.bcast.diag('text="gcamera alive"')
