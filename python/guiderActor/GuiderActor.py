#!/usr/bin/env python
"""An actor to run the guider"""

import abc

import opscore.actor.model
import opscore.actor.keyvar

import actorcore.Actor

import gcameraThread
import masterThread
import movieThread

import guiderActor
import GuiderState
from guiderActor import myGlobals

import warnings


def set_default_pids(config, gState):
    """Set the PID value defaults from the config file."""
    axes = dict(RADEC="raDec", ROT="rot", FOCUS="focus", SCALE="scale")
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


def set_telescope(config, gState):
    """Set values related to the telescope from the config file."""
    gState.plugPlateScale = float(config.get('telescope', "scale"))
    gState.dSecondary_dmm = float(config.get('telescope', "dSecondary_dmm"))
    gState.longitude = float(config.get('telescope', "longitude"))
    gState.focalRatio = float(config.get('telescope', "focalRatio"))


def set_gcamera(config, gState):
    """Set values related to the guide camera from the config file."""
    expTime = float(config.get('gcamera', "exposureTime"))
    readTime = float(config.get('gcamera', "binnedReadTime"))
    masterThread.set_time(gState, expTime, 1, readTime)
    gState.gcameraPixelSize = float(config.get('gcamera', "pixelSize"))
    gState.gcameraMagnification = float(config.get('gcamera', "magnification"))


class GuiderActor(actorcore.Actor.SDSSActor):
    """Manage the threads that calculate guiding corrections and gcamera commands."""
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def newActor(location=None, **kwargs):
        """Return the version of the actor based on our location."""
        location = GuiderActor._determine_location(location)
        if location == 'APO':
            return GuiderActorAPO('guider', productName='guiderActor', **kwargs)
        elif location == 'LCO':
            return GuiderActorLCO('guider', productName='guiderActor', **kwargs)
        elif location == 'LOCAL':
            return GuiderActorLocal('guider', productName='guiderActor', **kwargs)
        else:
            raise KeyError("Don't know my location: cannot return a working Actor!")

    def __init__(self, name, debugLevel=30, productName=None, makeCmdrConnection=True):
        actorcore.Actor.Actor.__init__(self, name, productName=productName, makeCmdrConnection=makeCmdrConnection)

        self.headURL = ('$HeadURL: svn+ssh://sdss3svn@sdss3.org/repo/ops/actors/guiderActor/branches/lco/'
                        'python/guiderActor/guiderActor_main.py $')

        self.logger.setLevel(debugLevel)
        self.logger.propagate = True

        #guiderActor.myGlobals.actorState = actorcore.Actor.ActorState(self)
        #actorState = guiderActor.myGlobals.actorState
        #self.actorState = actorState
        #actorState.gState = GuiderState.GuiderState()
        #actorState.actorConfig = self.config
        #gState = actorState.gState

        # Define thread list
        self.threadList = [("master", guiderActor.MASTER, masterThread),
                           ("gcamera", guiderActor.GCAMERA, gcameraThread),
                           ("movie", guiderActor.MOVIE, movieThread), ]

        # Load other actor's models so we can send it commands
        # And ours: we use the models to generate the FITS cards.
        self.models = {}
        for actor in ["gcamera", "ecamera", "mcp", "platedb", "sop", "tcc", "guider", "apo"]:
            self.models[actor] = opscore.actor.model.Model(actor)

        self.actorState = actorcore.Actor.ActorState(self, self.models)
        self.actorState.gState = GuiderState.GuiderState()
        self.actorState.actorConfig = self.config
        myGlobals.actorState = self.actorState
        self.actorState.timeout = 60  # timeout on message queues
        gState = self.actorState.gState

        for what in self.config.options("enable"):
            enable = {"True": True, "False": False}[self.config.get('enable', what)]
            gState.setGuideMode(what, enable)

        set_default_pids(self.config, gState)
        set_pid_scaling(self.config, gState)
        set_telescope(self.config, gState)
        set_gcamera(self.config, gState)

        # Sets the value of the bigFiberRadius
        gState.bigFiberRadius = float(self.config.get('gprobes', 'bigFiberRadius'))
        gState.zeropoint = float(self.config.get('telescope', 'zeroPoint'))

    def getLoadedCartridge(self, cmd, actor, actorState=None):
        """Returns the value of instrumentNum from actor."""

        instrumentNumKey = actorState.models[actor].keyVarDict["instrumentNum"]
        cmdVar = actorState.actor.cmdr.call(actor=actor,
                                            forUserCmd=cmd,
                                            cmdStr="info",
                                            keyVars=[instrumentNumKey])

        if cmdVar.didFail:
            cmd.fail('text=\"Failed to ask {0} for info on cartridges\"'.format(actor))
            return

        loadedCartridge = cmdVar.getLastKeyVarData(instrumentNumKey)[0]

        return loadedCartridge


class GuiderActorAPO(GuiderActor):
    """APO version of this actor."""
    location = 'APO'

    def guidingIsOK(self, cmd, actorState, force=False):
        """Is it OK to be guiding?"""

        if force:
            return True

        bypassedNames = actorState.models["sop"].keyVarDict["bypassedNames"]

        ffsStatus = actorState.models["mcp"].keyVarDict["ffsStatus"]
        open, closed = 0, 0
        for s in ffsStatus:
            if s is None:
                cmd.warn('text="Failed to get state of flat field screen from MCP"')
                break

            open += int(s[0])
            closed += int(s[1])

        if open != 8:
            msg = "FF petals aren\'t all open"
            if 'ffs' in bypassedNames:
                cmd.warn('text="%s; guidingIsOk failed, but ffs is bypassed in sop"' % msg)
            else:
                cmd.warn('text="%s; aborting guiding"' % msg)
                return False

        # This lets guiderImageAnalysis know to ignore dark frames.
        actorState.bypassDark = 'guider_dark' in bypassedNames

    #   should we allow guiding with lamps on if axes are disabled
    #   check if lamps are actually ON
        ffLamp = actorState.models["mcp"].keyVarDict["ffLamp"]
        hgCdLamp = actorState.models["mcp"].keyVarDict["hgCdLamp"]
        neLamp = actorState.models["mcp"].keyVarDict["neLamp"]
        if (any(ffLamp) and 'lamp_ff' not in bypassedNames) or \
           (any(hgCdLamp) and 'lamp_hgcd' not in bypassedNames) or \
           (any(neLamp) and 'lamp_ne' not in bypassedNames):
            cmd.warn('text="Calibration lamp on; aborting guiding"')
            return False

    #   check if non sensed lamps are commanded ON
        uvLamp = actorState.models["mcp"].keyVarDict["uvLampCommandedOn"]
        whtLamp = actorState.models["mcp"].keyVarDict["whtLampCommandedOn"]
        if uvLamp.getValue() or whtLamp.getValue():
            cmd.warn('text="Calibration lamp commanded on; aborting guiding"')
            return False

        tccModel = actorState.models['tcc']
        axisCmdState = tccModel.keyVarDict['axisCmdState']
        if any(x.lower() != 'tracking' for x in axisCmdState):
            if 'axes' in bypassedNames:
                cmd.warn('text="TCC motion failed, but axis motions are bypassed in sop"')
            else:
                cmd.warn('text="TCC motion aborted guiding"')
                return False

        # Checks the state of the primary and secondary mirrors
        mirror2CmdState = tccModel.keyVarDict['secState'][0]
        if mirror2CmdState.lower() != 'done':
            cmd.warn('text="Secondary mirror state: {0}. Aborted guiding."'
                     .format(mirror2CmdState))
            return False

        primaryCmdState = tccModel.keyVarDict['primState'][0]
        if primaryCmdState.lower() != 'done':
            cmd.warn('text="Primary mirror state: {0}. Aborted guiding."'.format(primaryCmdState))
            return False

        return True

    def getLoadedCartridge(self, cmd, actorState):
        """Returns the number of the cart loaded in the telescope.

        At APO, this information is provided by the MCP.

        """

        return super(GuiderActorAPO, self).getLoadedCartridge(
            cmd, 'mcp', actorState=actorState)


class GuiderActorLCO(GuiderActor):
    """LCO version of this actor."""
    location = 'LCO'

    def guidingIsOK(self, cmd, actorState, force=False):
        """Is it OK to be guiding?"""

        if force:
            return True

        bypassedNames = actorState.models["sop"].keyVarDict["bypassedNames"]
        # This lets guiderImageAnalysis know to ignore dark frames.

        actorState.bypassDark = 'guider_dark' in bypassedNames

        tccModel = actorState.models['tcc']
        axisCmdState = tccModel.keyVarDict['axisCmdState']

        if None in axisCmdState:
            cmd.warn('text="Could not get axisCmdState. Is the TCC connected?"')
            return False

        raState, decState, rotState = map(lambda xx: xx.lower(), axisCmdState)

        # At LCO rotator does not guide, so we check that it is halted.
        if raState != 'tracking' or decState != 'tracking' or rotState != 'halted':
            if 'axes' in bypassedNames:
                cmd.warn('text="TCC motion failed, but axis motions are bypassed in sop"')
            else:
                cmd.warn('text="TCC motion aborted guiding"')
                return False

        # Checks the state of the secondary mirror. Fails if it's moving.
        mirror2CmdState = tccModel.keyVarDict['secState'][0]
        if mirror2CmdState.lower() != 'done':
            cmd.warn('text="Secondary mirror state: {0}. Aborted guiding."'
                     .format(mirror2CmdState))
            return False
        
        # LCOHACK: disabling this part untils the scaling ring is connected
        # Checks the state of the scaling ring. Fails if it's moving.
        # scaleCmdState = tccModel.keyVarDict['threadringState'][0]
        # if scaleCmdState.lower() != 'done':
        #     cmd.warn('text="Scaling ring state: {0}. Aborted guiding."'.format(scaleCmdState))
        #     return False

        return True

    def getLoadedCartridge(self, cmd, actorState):
        """Returns the number of the cart loaded in the telescope.

        At LCO, this information is provided by the TCC.

        """

        return super(GuiderActorLCO, self).getLoadedCartridge(
            cmd, 'tcc', actorState=actorState)


class GuiderActorLocal(GuiderActor):
    """Test version of this actor. In prnciple, inherits from GuiderActorAPO."""
    location = 'LOCAL'

    def guidingIsOk(self, cmd, actorState, force=False):
        """Is it OK to be guiding?"""

        return True

    def getLoadedCartridge(self, cmd, actorState):
        """Returns the number of the cart loaded in the telescope."""

        pass
