"""
Thread for generating guider movies.

Calls guider_movie.py which calls fmpeg to generate the movies,
and we don't want to block on that.
"""
import Queue, threading
import subprocess
import os,os.path

import opscore.utility.tback as tback
from guiderActor import *
import guiderActor.myGlobals

class MovieMaker(object):
    """Holds the state of movie processing, including the current subprocess."""
    def __init__(self,actorState=None):
        self.popenMovie = None
        self.actorState = actorState if actorState is not None else guiderActor.myGlobals.actorState
    
    def __call__(self,msg):
        """
        Make the movie by calling guider_movie.py via subprocess.
        """
        mjd = msg.mjd
        if mjd is None:
            filedir = self.actorState.models['gcamera'].keyVarDict['dataDir'][0]
            simulating = self.actorState.models['gcamera'].keyVarDict['simulating']
            if simulating[0]:
                filedir = simulating[1]
            mjd = os.path.split(filedir)[-1]
        else:
            basedir = self.actorState.models['guider'].keyVarDict['file'][0]
            if basedir == None:
                basedir = '/data/gcam/'
            filedir = os.path.join(basedir,mjd)
        
        start = msg.start
        end = msg.end
        
        # jkp TBD: how to get the real output filename?
        #filename = # parse from stdout of subprocess?
        self.filename = os.path.join(filedir,'%s-%04d-%04d.mp4'%(mjd,start,end))
        # Takes about 1 second per frame to make the images, plus a bit more per frame to make the movie from them.
        timeLim = (end-start)*2+30
        
        cmdParts = ['guider_movie.py',filedir,str(start),str(end)]
        cmdAll = r' '.join(cmdParts)
        try:
            msg.cmd.inform('text="Making movie via: \'%s\'"'%cmdAll)
            # Using preexec_fn=os.setsid to allow:
            #    os.killpg(popen.pid,signal.SIGTERM)
            # to kill the whole process.
            self.popenMovie = subprocess.Popen(cmdParts,stdout=subprocess.PIPE,stderr=subprocess.PIPE,preexec_fn=os.setsid)
            msg.cmd.inform('moviePID=%d'%(self.popenMovie.pid,))
            stdout,stderr = self.popenMovie.communicate()
        except subprocess.CalledProcessError as e:
            msg.cmd.warn('text="%s returned %s"' % (cmdParts[0], e))
            raise e
        else:
            msg.cmd.inform('movieFile=%s'%(self.filename,))
            return self.filename
        finally:
            dbg_output = '\n'.join(('guider_movie.py STDOUT WAS:',stdout,'guider_movie.py STDERR WAS:',stderr))
            print dbg_output
    #...
    
    def check_movie(self,msg):
        """Check on the status of a currently processing movie generation."""
        pass
#...

def main(actor, queues):
    """Main for guider movie making thread."""
    
    threadName = "movie"
    
    actorState = guiderActor.myGlobals.actorState
    timeout = guiderActor.myGlobals.actorState.timeout
    movieMaker = MovieMaker(actorState)
    while True:
        try:
            msg = queues[MOVIE].get(timeout=timeout)
            qlen = queues[MOVIE].qsize()
            if qlen > 0 and msg.cmd:
                msg.cmd.diag("movie thread has %d items after a .get()" % (qlen))
                
            if msg.type == Msg.EXIT:
                if msg.cmd:
                    msg.cmd.inform('text="Exiting thread %s"' % (threading.current_thread().name))
                return
            
            elif msg.type == Msg.MAKE_MOVIE:
                try:
                    filename = movieMaker(msg)
                except subprocess.CalledProcessError as e:
                    msg.cmd.fail('text="Failed to create guider movie."')
                else:
                    # jkp: I don't think we need to send a Msg when we're done, since nothing needs to happen???
                    #msg.replyQueue.put(Msg(responseMsg, cmd=msg.cmd, filename=filename, success=True))
                    doneText = 'text="Created: %s"'%filename
                    if msg.finish:
                        msg.cmd.finish(doneText)
                    else:
                        msg.cmd.inform(doneText)

            #elif msg.type == Msg.ABORT_MOVIE:
            #    # Do something to cancel the currently running movie.
            #    guiderActor.flushQueue(queues[MOVIE])
            else:
                raise ValueError, ("Unknown message type %s" % msg.type)

        except Queue.Empty:
            actor.bcast.diag('text="guider movie alive"')
        except Exception as e:
            actor.bcast.diag('text="movie thread got unexpected exception: %s"' % (e))
            errMsg = "Unexpected exception %s in guider %s thread" % (e, threadName)
            tback.tback(errMsg, e)
    #...
#...
