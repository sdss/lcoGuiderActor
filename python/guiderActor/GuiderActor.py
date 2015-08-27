#!/usr/bin/env python
"""An actor to run the guider"""

import abc

import opscore.actor.model
import opscore.actor.keyvar

import actorcore.Actor

import gcameraThread
import masterThread
import movieThread

import GuiderState
import guiderActor.myGlobals

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

class GuiderActor(actorcore.Actor.SDSSActor):
    """Manage the threads that calculate guiding corrections and gcamera commands."""
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def newActor(location=None,**kwargs):
        """Return the version of the actor based on our location."""

        location = GuiderActor._determine_location(location)
        if location == 'apo':
            return GuiderActorAPO('guider',productName='guiderActor',**kwargs)
        elif location == 'lco':
            return GuiderActorLCO('guider',productName='guiderActor',**kwargs)
        else:
            raise KeyError("Don't know my location: cannot return a working Actor!")

    def __init__(self, name, configFile, debugLevel=30):
        actorcore.Actor.Actor.__init__(self, name, configFile)

        self.headURL = "$HeadURL: svn+ssh://sdss3svn@sdss3.org/repo/ops/actors/guiderActor/branches/lco/python/guiderActor/guiderActor_main.py $"

        self.logger.setLevel(debugLevel)

        guiderActor.myGlobals.actorState = actorcore.Actor.State(self)
        actorState = guiderActor.myGlobals.actorState
        actorState.gState = GuiderState.GuiderState()
        gState = actorState.gState

        # Load other actor's models so we can send it commands
        # And ours: we use the models to generate the FITS cards.
        for actor in ["gcamera", "ecamera", "mcp", "platedb", "sop", "tcc", "guider", "apo"]:
            actorState.models[actor] = opscore.actor.model.Model(actor)

        actorState.timeout = 60         # timeout on message queues

        self.threadList = [("master", guiderActor.MASTER, masterThread),
                           ("gcamera", guiderActor.GCAMERA, gcameraThread),
                           ("movie", guiderActor.MOVIE, movieThread),]

        expTime = float(self.config.get('gcamera', "exposureTime"))
        readTime = float(self.config.get('gcamera', "binnedReadTime"))
        masterThread.set_time(gState,expTime,readTime=readTime)

        gState.plugPlateScale = float(self.config.get('telescope', "scale"))
        gState.dSecondary_dmm = float(self.config.get('telescope', "dSecondary_dmm"))
        gState.longitude = float(self.config.get('telescope', "longitude"))
        gState.focalRatio = float(self.config.get('telescope', "focalRatio"))
        gState.gcameraPixelSize = float(self.config.get('gcamera', "pixelSize"))
        gState.gcameraMagnification = float(self.config.get('gcamera', "magnification"))

        for what in self.config.options("enable"):
            enable = {"True" : True, "False" : False}[self.config.get('enable', what)]
            gState.setGuideMode(what, enable)

        set_default_pids(self.config, gState)
        set_pid_scaling(self.config, gState)


class GuiderActorAPO(GuiderActor):
    """APO version of this actor."""
    location='APO'

class GuiderActorLCO(GuiderActor):
    """LCO version of this actor."""
    location='LCO'
