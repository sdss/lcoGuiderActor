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

    def __str__(self):
        msg = "%s %s %s" % (self.actor, self.actor.cmdr.dispatcher, self.dispatcher.reactor)

        return msg

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getActorState():
    print "RHL actorState", actorState
    return guiderActor.myGlobals.actorState

class Guider(actorcore.Actor.Actor):
    def __init__(self, name, configFile, debugLevel=30):
        actorcore.Actor.Actor.__init__(self, name, configFile)

        self.logger.setLevel(debugLevel)
        #self.logger.propagate = True

        self.cmdr = actorcore.CmdrConnection.Cmdr('%s.%s' % (name, name), self)
        self.cmdr.connectionMade = self.connectionMade
        self.cmdr.connect()

        guiderActor.myGlobals.actorState = State(self)
        actorState = guiderActor.myGlobals.actorState
        #
        # Load other actor's models so we can send it commands
        #
        for actor in ["gcamera", "platedb", "tcc"]:
            actorState.models[actor] = opscore.actor.model.Model(actor)
        #
        # spawn off the threads that sequence actions (e.g. take an exposure; move telescope)
        # and talk to the gcamera
        # 
        actorState.timeout = 60         # timeout on message queues

        actorState.queues = {}
        actorState.threads = {}

        actorState.queues[MASTER] = Queue.PriorityQueue(0)
        actorState.queues[GCAMERA] = Queue.PriorityQueue(0)
        
        actorState.threads[GCAMERA] = threading.Thread(target=gcameraThread.main, args=[self, actorState.queues])
        actorState.threads[GCAMERA].start()
        
        actorState.threads[MASTER] = threading.Thread(target=masterThread.main, args=[self, actorState.queues])
        actorState.threads[MASTER].start()
        #
        # Handle the hated ini file
        #
        import ConfigParser
        
        try:
            expTime = float(self.config.get('gcamera', "exposureTime"))
        except ConfigParser.NoOptionError:
            expTime = 10.0

        actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, None, expTime=expTime))
 
        plugPlateScale = float(self.config.get('plugPlate', "scale"))
        gcameraPixelSize = float(self.config.get('gcamera', "pixelSize"))
        gcameraMagnification = float(self.config.get('gcamera', "magnification"))
        actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_SCALE, None, 
                                                      plugPlateScale=plugPlateScale,
                                                      gcameraMagnification=gcameraMagnification,
                                                      gcameraPixelSize=gcameraPixelSize))
        #
        # Finally start the reactor
        #
        self.run()

    def connectionMade(self):
        self.bcast.warn("Guider is connected.")
        
#
# To work
#
if __name__ == "__main__":
    guider = Guider("guider", "guiderActor", debugLevel=5)
