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
                    msg.cmd.inform("text=\"Exiting thread %s\"" % (threading.current_thread().name))

                return
            
            elif msg.type == Msg.EXPOSE:
                msg.cmd.respond("text=\"starting exposure\"")
                #
                # Take exposure
                #
                aborted = False

                timeLim = msg.expTime + 15     # allow for readout time

                filenameKey = guiderActor.myGlobals.actorState.models["gcamera"].keyVarDict["filename"]

                cmdVar = actor.cmdr.call(actor="gcamera", cmdStr="expose time=%f" % (msg.expTime),
                                         keyVars=[filenameKey], timeLim=timeLim)
                if cmdVar.didFail:
                    msg.cmd.warn("text=\"Failed to take exposure\"")
                    queues[MASTER].put(Msg(Msg.START_GUIDING, msg.cmd, start=False))

                    continue

                filename = cmdVar.getLastKeyVarData(filenameKey)[0]

                queues[MASTER].put(Msg(Msg.EXPOSURE_FINISHED, None, filename=filename, aborted=aborted))
                #
                # Abort logic
                #
                if False:
                    msg2 = queues[GCAMERA].get(timeout=msg.expTime) # start fake exposure

                    cmd2 = msg2.cmd
                    #import pdb; pdb.set_trace()
                    if msg2.type == Msg.ABORT_EXPOSURE:
                        aborted = True

                        quiet = msg2.data.quiet
                        guiderActor.flushQueue(queues[GCAMERA])

                        if quiet:
                            cmd2.fail("text=\"Aborting exposure\"")
                        else: 
                           cmd2.finish()
                           
                        continue
                    elif msg2.type == Msg.EXPOSE:
                        cmd2.error("text=\"Attempt to start a guide exposure while one is active\"")
                    else:
                        raise RuntimeError, ("Unexpected message type %d" % msg.type)

                    queues[MASTER].put(Msg(Msg.EXPOSURE_FINISHED, None, filename=filename, aborted=aborted))
                continue
            elif msg.type == Msg.ABORT_EXPOSURE:
                if not msg.quiet:
                    msg.cmd.respond("text=\"Request to abort an exposure when none are in progress\"")
                guiderActor.flushQueue(queues[GCAMERA])
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)

        except Queue.Empty:
            actor.bcast.diag("text=\"gcamera alive\"")
