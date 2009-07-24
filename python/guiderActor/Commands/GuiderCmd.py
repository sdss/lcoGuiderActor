#!/usr/bin/env python

""" Wrap top-level ICC functions. """

import pdb
import logging
import pprint
import sys
import ConfigParser

import opscore.protocols.validation as validation
import opscore.protocols.keysformat as keysformat
import opscore.protocols.keys as keys
import opscore.protocols.types as types

import Commands.CmdSet
from opscore.utility.qstr import qstr

from guiderActor import *
import guiderActor
import guiderActor.myGlobals as myGlobals

class GuiderCmd(Commands.CmdSet.CmdSet):
    """ Wrap commands to the guider actor"""

    def __init__(self, actor):
        Commands.CmdSet.CmdSet.__init__(self, actor)
        
        #
        # Set the keyword dictionary
        #
        self.keys = keys.KeysDictionary("guidercommands", (1, 1),
                                        keys.Key("cartridge", types.Int(), help="A cartridge ID"),
                                        keys.Key("fibers", types.Int()*(1,None), help="A list of fibers"),
                                        keys.Key("plate", types.Int(), help="A plugplate ID"),
                                        keys.Key("pointing", types.String(),
                                                 help="A pointing for the given plugplate"),
                                        keys.Key("time", types.Float(), help="Exposure time for guider"),
                                        )

        keys.CmdKey.setKeys(self.keys)
        #
        # Declare commands
        #
        self.vocab = [
            #("guide", "(on) <time>", self.guideOn),
            #("guide", "(off)", self.guideOff),
            ("guide", "(on|off) [<time>]", self.guide),
            ("setTime", "<time>", self.setTime),
            ("disableFibers", "<fibers>", self.disableFibers),
            ("enableFibers", "<fibers>", self.enableFibers),
            ("loadCartridge", "<cartridge> <plate> [<pointing>]", self.loadCartridge),
            ('ping', '', self.ping),
            ('axes', '(on|off)', self.axes),
            ('focus', '(on|off)', self.focus),
            ('scale', '(on|off)', self.scale),
            ('status', '', self.status),
            ]
    #
    # Define commands' callbacks
    #
    def disableFibersImpl(self, cmd, enable=True):
        """Disable a set of fibers"""

        for f in cmd.cmd.keywords["fibers"].values:
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.ENABLE_FIBER, cmd, fiber=f, enable=enable))

        self.status(cmd)                # finishes this command

    def disableFibers(self, cmd):
        """Disable a set of fibers"""

        self.disableFibersImpl(cmd, enable=False)

    def enableFibers(self, cmd):
        """Enable a set of fibers"""

        self.disableFibersImpl(cmd, enable=True)

    def setTime(self, cmd):
        """Set the exposure time"""

        time = cmd.cmd.keywords["time"].values[0]
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, cmd, time=time))

    if False:
        def guideOn(self, cmd):
            """Turn guiding on"""

            expTime = cmd.cmd.keywords["time"].values[0]
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd, start=True, time=expTime))

        def guideOff(self, cmd):
            """Turn guiding off"""

            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd, start=False))
    else:
        def guide(self, cmd):
            """Turn guiding on or off"""

            on = "on" in cmd.cmd.keywords
            expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else None
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd, start=on, time=expTime))

    def loadCartridge(self, cmd):
        """Load a cartridge"""

        cartridge = cmd.cmd.keywords["cartridge"].values[0]
        plate = cmd.cmd.keywords["plate"].values[0]
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        #
        # Fake reading the plugmap file
        #
        fiberList = []
        for f in range(1, 12):
            fiberList.append((f, True))
        
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.LOAD_CARTRIDGE, cmd,
                                                                cartridge=cartridge, plate=plate, pointing=pointing,
                                                                fiberList=fiberList))

    def ping(self, cmd):
        """ Top-level "ping" command handler. Query the actor for liveness/happiness. """

        cmd.finish('text="pong"')

    def correctionImpl(self, cmd, what):
        """Turn guiding something on or off"""

        on = "on" in cmd.cmd.keywords
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_GUIDE_MODE, cmd, what=what, enabled=on))

    def axes(self, cmd):
        """Turn guiding the plate axes on or off"""

        self.correctionImpl(cmd, "axes")

    def focus(self, cmd):
        """Turn guiding the plate focus on or off"""

        self.correctionImpl(cmd, "focus")

    def scale(self, cmd):
        """Turn guiding the plate scale on or off"""

        self.correctionImpl(cmd, "scale")

    def status(self, cmd):
        """Return guide status status"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STATUS, cmd, finish=True))
