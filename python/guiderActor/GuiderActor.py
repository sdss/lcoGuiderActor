#!/usr/bin/env python
"""An actor to run the guider"""

import abc
import os
from distutils import version

import actorcore.Actor
import gcameraThread
import guiderActor
import GuiderState
import masterThread
import movieThread
import opscore.actor.keyvar
import opscore.actor.model
from guiderActor import myGlobals


def set_default_pids(config, gState):
    """Set the PID value defaults from the config file."""
    axes = dict(RADEC='raDec', ROT='rot', FOCUS='focus', SCALE='scale')
    for axis in config.options('PID'):
        axis = axes[axis.upper()]
        Kp, Ti_min, Ti_max, Td, Imax, nfilt, ncorr = [float(v)
                                                      for v in config.get('PID', axis).split()]
        gState.set_pid_defaults(axis, Kp=Kp, Ti_min=Ti_min, Ti_max=Ti_max,
                                Td=Td, Imax=Imax, nfilt=int(nfilt), ncorr=ncorr)
        gState.pid[axis].setPID(Kp=Kp, Ti=Ti_min, Td=Td, Imax=Imax, nfilt=int(nfilt),
                                ncorr=ncorr)


def set_pid_scaling(config, gState):
    """Set the min/max altitude and the axes to scale the PID terms on."""
    gState.axes_to_scale = config.get('PID_altitude_scale', 'axes').split()
    gState.alt_min = float(config.get('PID_altitude_scale', 'min'))
    gState.alt_max = float(config.get('PID_altitude_scale', 'max'))


def set_telescope(config, gState):
    """Set values related to the telescope from the config file."""
    gState.plugPlateScale = float(config.get('telescope', 'scale'))
    gState.dSecondary_dmm = float(config.get('telescope', 'dSecondary_dmm'))
    gState.longitude = float(config.get('telescope', 'longitude'))
    gState.focalRatio = float(config.get('telescope', 'focalRatio'))


def set_gcamera(config, gState):
    """Set values related to the guide camera from the config file."""
    expTime = float(config.get('gcamera', 'exposureTime'))
    readTime = float(config.get('gcamera', 'binnedReadTime'))
    masterThread.set_time(gState, expTime, 1, readTime)
    gState.gcameraPixelSize = float(config.get('gcamera', 'pixelSize'))
    gState.gcameraMagnification = float(config.get('gcamera', 'magnification'))


class GuiderActor(actorcore.Actor.SDSSActor):
    """Manage the threads that calculate guiding corrections and gcamera commands."""
    __metaclass__ = abc.ABCMeta

    @staticmethod
    def newActor(location=None, **kwargs):
        """Return the version of the actor based on our location."""
        location = GuiderActor._determine_location(location)
        if location == 'APO':
            raise ValueError('this actor cannot run at APO.')
        elif location == 'LCO':
            return GuiderActorLCO('guider', productName='guiderActor', **kwargs)
        elif location == 'LOCAL':
            return GuiderActorLocal('guider', productName='guiderActor', **kwargs)
        else:
            raise KeyError("Don't know my location: cannot return a working Actor!")

    def __init__(self, name, debugLevel=30, productName=None, makeCmdrConnection=True):

        actorcore.Actor.Actor.__init__(self, name, productName=productName,
                                       productDir=(os.path.dirname(__file__) + '/../../'),
                                       makeCmdrConnection=makeCmdrConnection)

        self.version = guiderActor.__version__

        self.logger.setLevel(debugLevel)
        self.logger.propagate = True

        self.check_versions()

        # guiderActor.myGlobals.actorState = actorcore.Actor.ActorState(self)
        # actorState = guiderActor.myGlobals.actorState
        # self.actorState = actorState
        # actorState.gState = GuiderState.GuiderState()
        # actorState.actorConfig = self.config
        # gState = actorState.gState

        # Define thread list
        self.threadList = [('master', guiderActor.MASTER, masterThread),
                           ('gcamera', guiderActor.GCAMERA, gcameraThread),
                           ('movie', guiderActor.MOVIE, movieThread), ]

        # Load other actor's models so we can send it commands
        # And ours: we use the models to generate the FITS cards.
        self.models = {}
        for actor in ['gcamera', 'ecamera', 'platedb', 'sop', 'tcc', 'guider']:
            self.models[actor] = opscore.actor.model.Model(actor)

        self.actorState = actorcore.Actor.ActorState(self, self.models)
        self.actorState.gState = GuiderState.GuiderState()
        self.actorState.actorConfig = self.config
        myGlobals.actorState = self.actorState
        self.actorState.timeout = 60  # timeout on message queues
        gState = self.actorState.gState

        self.actorState.bypassDark = False

        for what in self.config.options('enable'):
            enable = {'True': True, 'False': False}[self.config.get('enable', what)]
            gState.setGuideMode(what, enable)

        set_default_pids(self.config, gState)
        set_pid_scaling(self.config, gState)
        set_telescope(self.config, gState)
        set_gcamera(self.config, gState)

        # Sets the value of the bigFiberRadius
        gState.bigFiberRadius = float(self.config.get('gprobes', 'bigFiberRadius'))
        gState.zeropoint = float(self.config.get('telescope', 'zeroPoint'))

    def getLoadedCartridge(self, cmd, actor, command='info', actorState=None):
        """Returns the value of instrumentNum from actor."""

        instrumentNumKey = actorState.models[actor].keyVarDict['instrumentNum']

        cmdVar = actorState.actor.cmdr.call(actor=actor,
                                            forUserCmd=cmd,
                                            cmdStr=command,
                                            keyVars=[instrumentNumKey])

        if cmdVar.didFail:
            cmd.warn('text=\"Failed to ask {0} for info on cartridges\"'
                     .format(actor))
            return

        loadedCartridge = actorState.models[actor].keyVarDict['instrumentNum'][0]

        return loadedCartridge

    def check_versions(self):
        """Checks whether library versions are ok.

        This is not a mandatory method, but specific locations can override it if needed.

        """

        pass


class GuiderActorLCO(GuiderActor):
    """LCO version of this actor."""
    location = 'LCO'
    astropy_max_version = '1.1.2'

    def guidingIsOK(self, cmd, actorState, force=False):
        """Is it OK to be guiding?"""

        if force:
            return True

        # Forces the TCC to output the status of the axes
        cmdVar = actorState.actor.cmdr.call(actor='tcc',
                                            forUserCmd=cmd,
                                            cmdStr='device status')

        if cmdVar.didFail:
            cmd.warn('text=\"Failed to get tcc status\"')
            return False

        bypassedNames = actorState.models['sop'].keyVarDict['bypassedNames']
        # This lets guiderImageAnalysis know to ignore dark frames.

        actorState.bypassDark = 'guider_dark' in bypassedNames

        tccModel = actorState.models['tcc']
        axisCmdState = tccModel.keyVarDict['axisCmdState']

        if None in axisCmdState:
            cmd.warn('text="Could not get axisCmdState. Is the TCC connected?"')
            return False

        raState, decState, rotState = map(lambda xx: xx.lower(), axisCmdState)

        # At LCO rotator does not guide, so we check that it is halted.
        if raState != 'tracking' or decState != 'tracking' or rotState != 'tracking':
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

        # Checks the state of the scaling ring. Fails if it's moving.
        scaleCmdState = tccModel.keyVarDict['threadringState'][0]
        if scaleCmdState.lower() != 'done':
            cmd.warn('text="Scaling ring state: {0}. Aborted guiding."'.format(scaleCmdState))
            return False

        return True

    def getLoadedCartridge(self, cmd, actorState):
        """Returns the number of the cart loaded in the telescope.

        At LCO, this information is provided by the TCC.

        """

        return super(GuiderActorLCO, self).getLoadedCartridge(
            cmd, 'tcc', command='device status scale', actorState=actorState)

    def check_versions(self):
        """Check that astropy is installed and it's the right version."""

        try:
            import astropy
        except ImportError:
            raise ImportError('cannot import astropy. '
                              'Please, make sure astropy <= {0} '
                              'is installed.'.format(self.astropy_max_version))

        if (version.StrictVersion(astropy.__version__) >
                version.StrictVersion(self.astropy_max_version)):
            raise ValueError('astropy needs to be <= {0}'.format(self.astropy_max_version))


class GuiderActorLocal(GuiderActor):
    """Test version of this actor. In prnciple, inherits from GuiderActorAPO."""
    location = 'LOCAL'

    def guidingIsOk(self, cmd, actorState, force=False):
        """Is it OK to be guiding?"""

        return True

    def getLoadedCartridge(self, cmd, actorState):
        """Returns the number of the cart loaded in the telescope."""

        pass
