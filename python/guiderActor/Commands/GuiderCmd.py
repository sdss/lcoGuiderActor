#!/usr/bin/env python

""" Wrap top-level ICC functions. """

import pdb
import logging
import re, sys
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
        self.keys = keys.KeysDictionary("guider_guider", (1, 1),
                                        keys.Key("cartridge", types.Int(), help="A cartridge ID"),
                                        keys.Key("fibers", types.Int()*(1,None), help="A list of fibers"),
                                        keys.Key("pointing", types.String(),
                                                 help="A pointing for the given plugplate"),
                                        keys.Key("expTime", types.Float(), help="Exposure time for guider"),
                                        )
        #
        # Declare commands
        #
        self.vocab = [
            ("guide", "(on|off) [<expTime>]", self.guide),
            ("setExpTime", "<expTime>", self.setExpTime),
            ("disableFibers", "<fibers>", self.disableFibers),
            ("enableFibers", "<fibers>", self.enableFibers),
            ("loadCartridge", "<cartridge> [<pointing>]", self.loadCartridge),
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
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.ENABLE_FIBER, cmd=cmd, fiber=f, enable=enable))

        self.status(cmd)                # finishes this command

    def disableFibers(self, cmd):
        """Disable a set of fibers"""

        self.disableFibersImpl(cmd, enable=False)

    def enableFibers(self, cmd):
        """Enable a set of fibers"""

        self.disableFibersImpl(cmd, enable=True)

    def setExpTime(self, cmd):
        """Set the exposure time"""

        expTime = cmd.cmd.keywords["expTime"].values[0]
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, cmd=cmd, expTime=expTime))

    def guide(self, cmd):
        """Turn guiding on or off"""

        on = "on" in cmd.cmd.keywords
        expTime = cmd.cmd.keywords["expTime"].values[0] if "expTime" in cmd.cmd.keywords else None
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd,
                                                                start=on, expTime=expTime))

    def loadCartridge(self, cmd):
        """Load a cartridge"""

        cartridge = cmd.cmd.keywords["cartridge"].values[0]
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        #
        # Get the plate from the plateDB
        #
        actorState = guiderActor.myGlobals.actorState

        plateKey = actorState.models["platedb"].keyVarDict["pointingInfo"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb",
                                            cmdStr="loadCartridge cartridge=%d pointing=%s" % (cartridge, pointing),
                                            keyVars=[plateKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup plate corresponding to %d/%s\"" % (cartridge, pointing))
            return

        plate = cmdVar.getLastKeyVarData(plateKey)[0]
        #
        # Lookup the valid gprobes
        #
        gprobeKey = actorState.models["platedb"].keyVarDict["gprobe"]
        gprobesInUseKey = actorState.models["platedb"].keyVarDict["gprobesInUse"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb",
                                            cmdStr="getGprobes cartridge=%d" % (cartridge),
                                            keyVars=[gprobeKey, gprobesInUseKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup gprobes for cartridge %d\"" % (cartridge))
            return

        enabled = {}
        for el in cmdVar.getLastKeyVarData(gprobesInUseKey):
            mat = re.search(r"^\((\d+)\s*=\s*(True|False)\s*\)$", el)
            id, isEnabled = int(mat.group(1)), (mat.group(2) == 'True')

            enabled[id] = isEnabled

        gprobes = []
        for cartridgeID, gpID, exists, rotation, focusOffset in cmdVar.getKeyVarData(gprobeKey):
            gprobes.append((gpID, exists, enabled.get(gpID, False), rotation, focusOffset))

        spiderInstAngKey = actorState.models["tcc"].keyVarDict["spiderInstAng"]
        print "RHL", spiderInstAngKey[0].getPos()

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                                                                cartridge=cartridge, plate=plate, pointing=pointing,
                                                                gprobes=gprobes))

    def ping(self, cmd):
        """ Top-level "ping" command handler. Query the actor for liveness/happiness. """

        cmd.finish('text="pong"')

    def correctionImpl(self, cmd, what):
        """Turn guiding something on or off"""

        on = "on" in cmd.cmd.keywords
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_GUIDE_MODE, cmd=cmd, what=what, enabled=on))

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

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STATUS, cmd=cmd, finish=True))
