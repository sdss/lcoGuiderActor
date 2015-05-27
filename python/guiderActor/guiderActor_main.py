#!/usr/bin/env python
"""An actor to run the guider"""

import re
import Queue, threading

import opscore.actor.model
import opscore.actor.keyvar

import actorcore.Actor

import gcameraThread
import masterThread
import movieThread

from guiderActor import Msg
import GuiderState
import guiderActor.myGlobals

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class State(object):
    """An object to hold globally useful state"""

    def __init__(self, actor):
        self.actor = actor
        self.dispatcher = self.actor.cmdr.dispatcher
        self.models = {}
        self.restartCmd = None

    def __str__(self):
        msg = "%s %s" % (self.actor, self.actor.cmdr.dispatcher)

        return msg

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def set_default_pids(config, gState):
    """Set the PID value defaults from the config file."""
    axes = dict(RADEC = "raDec", ROT = "rot", FOCUS = "focus", SCALE = "scale")
    for axis in config.options('PID'):
        axis = axes[axis.upper()]
        Kp, Ti_min, Ti_max, Td, Imax, nfilt = [float(v) for v in config.get('PID', axis).split()]        
        gState.set_pid_defaults(axis, Kp=Kp, Ti_min=Ti_min, Ti_max=Ti_max, Td=Td, Imax=Imax, nfilt=int(nfilt))
        gState.pid[axis].setPID(Kp=Kp, Ti=Ti_min, Td=Td, Imax=Imax, nfilt=nfilt)

def set_pid_scaling(config, gState):
    """Set the min/max altitude and the axes to scale the PID terms on."""
    gState.axes_to_scale = config.get('PID_altitude_scale', 'axes').split()
    gState.alt_min = float(config.get('PID_altitude_scale', 'min'))
    gState.alt_max = float(config.get('PID_altitude_scale', 'max'))

class Guider(actorcore.Actor.Actor):
    def __init__(self, name, configFile, debugLevel=30):
        actorcore.Actor.Actor.__init__(self, name, configFile)

        self.headURL = "$HeadURL$"

        self.logger.setLevel(debugLevel)
        #self.logger.propagate = True

        guiderActor.myGlobals.actorState = State(self)
        actorState = guiderActor.myGlobals.actorState
        actorState.gState = GuiderState.GuiderState()
        #
        # Load other actor's models so we can send it commands
        # And ours: we use the models to generate the FITS cards.
        #
        for actor in ["gcamera", "ecamera", "mcp", "platedb", "sop", "tcc", "guider", "apo"]:
            actorState.models[actor] = opscore.actor.model.Model(actor)
        #
        # spawn off the threads that sequence actions (e.g. take an exposure; move telescope)
        # and talk to the gcamera
        # 
        actorState.timeout = 60         # timeout on message queues

        Guider.startThreads(actorState, restartQueues=True)

        # Handle the hated ini file
        expTime = float(self.config.get('gcamera', "exposureTime"))
        readTime = float(self.config.get('gcamera', "binnedReadTime"))
        actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, None, expTime=expTime, readTime=readTime))
 
        plugPlateScale = float(self.config.get('telescope', "scale"))
        dSecondary_dmm = float(self.config.get('telescope', "dSecondary_dmm"))
        gcameraPixelSize = float(self.config.get('gcamera', "pixelSize"))
        gcameraMagnification = float(self.config.get('gcamera', "magnification"))
        actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_SCALE, None, 
                                                      plugPlateScale=plugPlateScale,
                                                      dSecondary_dmm=dSecondary_dmm,
                                                      gcameraMagnification=gcameraMagnification,
                                                      gcameraPixelSize=gcameraPixelSize))
        for what in self.config.options("enable"):
            enable = {"True" : True, "False" : False}[self.config.get('enable', what)]
            actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_GUIDE_MODE, None, what=what,
                                                          enable=enable))

        set_default_pids(self.config, actorState.gState)
        set_pid_scaling(self.config, actorState.gState)

        #
        # Finally start the reactor
        #
        self.run()

    @staticmethod
    def startThreads(actorState, cmd=None, restartQueues=False, restart=False):
        """Start or restart the worker threads and queues"""

        try:
            actorState.threads
        except AttributeError:
            restart = False
        
        if not restart:
            actorState.queues = {}
            actorState.threads = {}

            restartQueues = True

        newQueues = {}
        threadsToStart = []
        for tname, tid, threadModule in [("master", guiderActor.MASTER, masterThread),
                                         ("gcamera", guiderActor.GCAMERA, gcameraThread),
                                         ("movie", guiderActor.MOVIE, movieThread),
                                         ]:

            newQueues[tid] = Queue.Queue(0) if restartQueues else actorState.queues[tid]

            if restart:
                reload(threadModule)

                for t in threading.enumerate():
                    if re.search(r"^%s(-\d+)?$" % tname, t.name): # a thread of the proper type
                        guiderActor.flushQueue(actorState.queues[tid])
                        actorState.queues[tid].put(Msg(Msg.EXIT, cmd=cmd))

                        t.join(1.0)
                        if t.isAlive():
                            if cmd:
                                cmd.inform('text="Failed to kill %s"' % tname)

                def updateName(g):
                    """re.sub callback to convert master -> master-1; master-3 -> master-4"""
                    try:
                        n = int(g.group(2))
                    except TypeError:
                        n = 0

                    return "%s-%d" % (g.group(1), n + 1)

                tname = re.sub(r"^([^\d]*)(?:-(\d*))?$", updateName, actorState.threads[tid].name)

            actorState.threads[tid] = threading.Thread(target=threadModule.main, name=tname,
                                                       args=[actorState.actor, newQueues])
            actorState.threads[tid].daemon = True

            threadsToStart.append(actorState.threads[tid])
        #
        # Switch to the new queues now we've sent EXIT to the old ones
        #
        for tid, q in newQueues.items():
            actorState.queues[tid] = q

        for t in threadsToStart:
            t.start()
#
# To work
#
if __name__ == "__main__":
    guider = Guider("guider", "guiderActor", debugLevel=5)
