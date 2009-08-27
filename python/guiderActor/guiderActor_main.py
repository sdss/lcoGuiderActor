#!/usr/bin/env python
"""An actor to run the guider"""

import inspect, os, re, sys
import Queue, threading

import opscore.actor.model
import opscore.actor.keyvar

import actorcore.Actor
import actorcore.CmdrConnection

import actorkeys

import gcameraThread
import masterThread
        
from guiderActor import *
import guiderActor.myGlobals
#
# Import sdss3logging before logging if you want to use it
#
if False:
    import opscore.utility.sdss3logging as sdss3logging
import logging

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class State(object):
    """An object to hold globally useful state"""

    def __init__(self, actor):
        self.actor = actor
        self.dispatcher = self.actor.cmdr.dispatcher
        self.reactor = self.dispatcher.reactor
        self.models = {}
        self.restartCmd = None

    def __str__(self):
        msg = "%s %s %s" % (self.actor, self.actor.cmdr.dispatcher, self.dispatcher.reactor)

        return msg

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        
class TCCState(object):
    """Listen to TCC keywords to analyse when the TCC's done a move or halt that might invalidate guiding"""

    instName = None
    slewing = False
    halted = False

    def __init__(self, tccModel):
        """Register keywords that we need to pay attention to"""

        tccModel.keyVarDict["moveItems"].addCallback(self.listenToMoveItems, callNow=False)
        tccModel.keyVarDict["inst"].addCallback(self.listenToInst, callNow=False)
        tccModel.keyVarDict["slewEnd"].addCallback(self.listenToSlewEnd, callNow=False)
        tccModel.keyVarDict["tccStatus"].addCallback(self.listenToTccStatus, callNow=False)
        
    @staticmethod
    def listenToInst(keyVar):
        inst = keyVar[0]
        if inst != TCCState.instName:
            TCCState.instName = inst

    @staticmethod
    def listenToMoveItems(keyVar):
        """ Figure out if the telescope has been moved by examining the TCC's MoveItems key.

        The MoveItems keys always comes with one of:
          - Moved, indicating an immediate offset
              We guestimate the end time, and let the main loop wait on it.
          - SlewBeg, indicating the start of a real slew or a computed offset.
          - SlewEnd, indicating the end of a real slew or a computed offset
        """

        #
        # Ughh.  Look for those keys somewhere else on the same reply; Craig swears that this is safe...
        #
        key = [k for k in keyVar.reply.keywords if re.search(r"^(moved|slewBeg|slewEnd)$", k.name, re.IGNORECASE)]
        if not key:
            print "Craig lied to me"
            return
        key = key[0]
        #
        # Parse moveItems
        #
        mi = keyVar[0]

        if not mi:
            mi = "XXXXXXXXX"

        newField = False
        if mi[1] == 'Y':
            what = 'Telescope has been slewed'
            newField = True
        elif mi[3] == 'Y':
            what = 'Object offset was changed'
        elif mi[4] == 'Y':
            what = 'Arc offset was changed'
        elif mi[5] == 'Y':
            what = 'Boresight position was changed'
        elif mi[6] == 'Y':
            what = 'Rotation was changed'
        elif mi[8] == 'Y':
            what = 'Calibration offset was changed'
        else:
            # Ignoring:
            #    mi[0] - Object name
            #    mi[2] - Object magnitude
            #    mi[7] - Guide offset
            return
        #
        # respond to what's been ordered
        #
        if key.name == "Moved":         # Has an uncomputed offset just been issued?
            # telescope has started moving, and may still be moving
            TCCState.goToNewField = False
            TCCState.halted = False
            TCCState.moved = True
            TCCState.slewing = False

        elif key.name == "SlewBeg":     # Has a slew/computed offset just been issued? 
            TCCState.goToNewField = newField
            TCCState.halted = False
            TCCState.moved = False
            TCCState.slewing = True

        elif key.name == "SlewEnd":     # A computed offset finished
            listenToSlewEnd(None)
        else:
            print "Impossible key:", key.name

        print "TCCState.slewing", TCCState.slewing

    @staticmethod
    def listenToSlewEnd(keyVar):
        """Listen for SlewEnd commands"""

        TCCState.goToNewField = False
        TCCState.halted = False
        TCCState.moved = False
        TCCState.slewing = False
        

    @staticmethod
    def listenToTccStatus(keyVar):
        """ Figure out if the telescope has been halted by examining the TCC's TCCStatus key.
        """

        axisStat = keyVar.valueList[0]

        if axisStat is not None and 'H' in axisStat:
            TCCState.goToNewField = False
            TCCState.halted = True

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getActorState():
    print "RHL actorState", actorState
    return guiderActor.myGlobals.actorState

class Guider(actorcore.Actor.Actor):
    def __init__(self, name, configFile, debugLevel=30):
        actorcore.Actor.Actor.__init__(self, name, configFile)

        self.logger.setLevel(debugLevel)
        #self.logger.propagate = True

        self.cmdr = actorcore.CmdrConnection.Cmdr(name, self)
        self.cmdr.connectionMade = self.connectionMade
        self.cmdr.connect()

        guiderActor.myGlobals.actorState = State(self)
        actorState = guiderActor.myGlobals.actorState
        #
        # Load other actor's models so we can send it commands
        #
        for actor in ["gcamera", "mcp", "platedb", "tcc"]:
            actorState.models[actor] = opscore.actor.model.Model(actor)
        #
        # spawn off the threads that sequence actions (e.g. take an exposure; move telescope)
        # and talk to the gcamera
        # 
        actorState.timeout = 60         # timeout on message queues

        Guider.startThreads(actorState, startQueues=True)
        #
        # Start listening to the TCC's keywords that announce that it's done a move or halt
        # that might invalidate guiding
        #
        actorState.tccState = TCCState(actorState.models["tcc"])
        #
        # Handle the hated ini file
        #
        import ConfigParser
        
        try:
            expTime = float(self.config.get('gcamera', "exposureTime"))
        except ConfigParser.NoOptionError:
            expTime = 10.0

        actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, None, expTime=expTime))
 
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

        for what in self.config.options("PID"):
            Kp, Ti, Td, Imax = [float(v) for v in self.config.get('PID', what).split()]
            actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_PID, None, what=what,
                                                          Kp=Kp, Ti=Ti, Td=Td, Imax=Imax))

        #
        # Finally start the reactor
        #
        self.run()

    def connectionMade(self):
        self.bcast.warn("Guider is connected.")
        
    @staticmethod
    def startThreads(actorState, cmd=None, startQueues=False, restart=False):
        """Start or restart the worker threads and queues"""

        try:
            actorState.threads
        except AttributeError:
            restart = False
        
        if not restart:
            actorState.queues = {}
            actorState.threads = {}

            startQueues = True

        for tname, tid, threadModule in [("master", guiderActor.MASTER, masterThread),
                                         ("gcamera", guiderActor.GCAMERA, gcameraThread),
                                         ]:

            if startQueues:
                actorState.queues[tid] = Queue.PriorityQueue(0)

            if restart:
                reload(threadModule)

                for t in threading.enumerate():
                    if re.search(r"^%s(-\d+)?$" % tname, t.name): # a thread of the proper type
                        guiderActor.flushQueue(actorState.queues[tid])
                        actorState.queues[tid].put(Msg(Msg.EXIT, cmd=cmd))

                        t.join(1.0)

                def updateName(g):
                    """re.sub callback to convert master -> master-1; master-3 -> master-4"""
                    try:
                        n = int(g.group(2))
                    except TypeError:
                        n = 0

                    return "%s-%d" % (g.group(1), n + 1)

                tname = re.sub(r"^([^\d]*)(?:-(\d*))?$", updateName, actorState.threads[tid].name)

            actorState.threads[tid] = threading.Thread(target=threadModule.main, name=tname,
                                                       args=[actorState.actor, actorState.queues])
            actorState.threads[tid].start()
#
# To work
#
if __name__ == "__main__":
    guider = Guider("guider", "guiderActor", debugLevel=5)
