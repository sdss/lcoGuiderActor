import Queue, time

from guiderActor import *
import guiderActor.myGlobals

def daemon(actor, queues):
    timeout = guiderActor.myGlobals.actorState.timeout

    expNum = 0
    
    while True:
        try:
            msg = queues[GCAMERA].get(timeout=timeout)
            
            if msg.type == Msg.EXPOSE:
                cmd, expTime = msg.data

                cmd.respond("text=\"starting exposure\"")
                #
                # Take exposure
                #
                aborted = False

                timeLim = expTime + 15     # allow for readout time

                filenameKey = guiderActor.myGlobals.actorState.models["gcamera"].keyVarDict["filename"]

                cmdVar = actor.cmdr.call(actor="gcamera", cmdStr="expose time=%f" % (expTime),
                                         keyVars=[filenameKey], timeLim=timeLim)
                if cmdVar.didFail:
                    cmd.fail("Failed to take exposure")
                    continue

                filename = cmdVar.getLastKeyVarData(filenameKey)

                queues[MASTER].put(Msg(Msg.EXPOSURE_FINISHED, (filename, aborted)))
                #
                # Abort logic
                #
                if False:
                    msg2 = queues[GCAMERA].get(timeout=expTime) # start fake exposure

                    cmd2 = msg2.data[0]
                    #import pdb; pdb.set_trace()
                    if msg2.type == Msg.ABORT_EXPOSURE:
                        aborted = True

                        quiet = msg2.data[1]
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

                    queues[MASTER].put(Msg(Msg.EXPOSURE_FINISHED, (filename, aborted)))
                continue
            elif msg.type == Msg.ABORT_EXPOSURE:
                cmd, quiet = msg.data
                if not quiet:
                    cmd.respond("text=\"Request to abort an exposure when none are in progress\"")
                guiderActor.flushQueue(queues[GCAMERA])
            else:
                raise ValueError, ("Unknown message type %d" % msg.type)

        except Queue.Empty:
            actor.bcast.diag("text=\"gcamera alive\"")
