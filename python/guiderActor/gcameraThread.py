import Queue, threading

from guiderActor import *
import guiderActor.myGlobals

def main(actor, queues):
    """Main look for thread to talk to gcamera"""

    timeout = guiderActor.myGlobals.actorState.timeout

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
                try:
                    camera = msg.camera
                except:
                    camera = "gcamera"

                msg.cmd.respond('text="starting exposure"')
                #
                # Take exposure
                #
                try:
                    expType = msg.expType
                except:
                    expType = "expose" 

                timeLim = msg.expTime + 45     # allow for readout time
                
                filenameKey = guiderActor.myGlobals.actorState.models[camera].keyVarDict["filename"]

                try:
                    forTCC = msg.forTCC
                except:
                    forTCC = False

                    
                cmdStr="%s time=%f" % (expType, msg.expTime)
                if expType == "flat":
                    cmdStr += " cartridge=%s" % (msg.cartridge)
                    responseMsg = Msg.FLAT_FINISHED
                else:
                    try:
                        stack = msg.stack
                        cmdStr += " stack=%d" % (stack)
                    except:
                        stack = 1

                    responseMsg = Msg.EXPOSURE_FINISHED

                actor.bcast.diag('text="%s %s with timeLim=%s"' % (camera, cmdStr, timeLim))
                try:
                    cmdVar = actor.cmdr.call(actor=camera, cmdStr=cmdStr, 
                                             keyVars=[filenameKey], timeLim=timeLim, forUserCmd=msg.cmd)
                    actor.bcast.diag('text="%s %s didFail=%s"' % (camera, cmdStr, cmdVar.didFail))
                except Exception, e:
                    actor.bcast.diag('text="%s %s raised %s"' % (camera, cmdStr, e))
                    continue

                if cmdVar.didFail:
                    msg.cmd.warn('text="Failed to take exposure"')
                    msg.replyQueue.put(Msg(responseMsg, cmd=msg.cmd, success=False, forTCC=forTCC))
                    continue

                filename = cmdVar.getLastKeyVarData(filenameKey)[0]

                print "Sending EXPOSURE_FINISHED to", msg.replyQueue
                msg.replyQueue.put(Msg(responseMsg, cmd=msg.cmd, filename=filename, success=True, forTCC=forTCC))

            elif msg.type == Msg.ABORT_EXPOSURE:
                if not msg.quiet:
                    msg.cmd.respond('text="Request to abort an exposure when none are in progress"')
                guiderActor.flushQueue(queues[GCAMERA])
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)

        except Queue.Empty:
            actor.bcast.diag('text="gcamera alive"')
        except Exception, e:
            actor.bcast.diag('text="gcamera thread got unexpected exception: %s"' % (e))
            
