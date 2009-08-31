#!/usr/bin/env python

""" Wrap top-level ICC functions. """

import pdb
import logging
import re, sys
import threading

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

from guiderActor import *
import guiderActor
import guiderActor.myGlobals as myGlobals
import masterThread
import gcameraThread

class GuiderCmd(object):
    """ Wrap commands to the guider actor"""

    class GprobeInfo(object):
        """Capture information about a guider probe"""
        def __init__(self, exists, enabled, xCenter, yCenter, radius, rotation,
                     xFerruleOffset, yFerruleOffset, focusOffset, fiber_type):
            self.exists = exists
            self.enabled = enabled
            self.xCenter = xCenter
            self.yCenter = yCenter
            self.radius = radius
            self.xFerruleOffset = xFerruleOffset
            self.yFerruleOffset = yFerruleOffset
            self.rotation = rotation
            self.focusOffset = focusOffset
            self.fiber_type = fiber_type

    def __init__(self, actor):
        self.actor = actor
        #
        # Declare keys that we're going to use
        #
        self.keys = keys.KeysDictionary("guider_guider", (1, 1),
                                        keys.Key("cartridge", types.Int(), help="A cartridge ID"),
                                        keys.Key("fibers", types.Int()*(1,None), help="A list of fibers"),
                                        keys.Key("pointing", types.String(),
                                                 help="A pointing for the given plugplate"),
                                        keys.Key("time", types.Float(), help="Exposure time for guider"),
                                        keys.Key("force", help="Force requested action to happen"),
                                        keys.Key("gprobes", types.Enum("acquire", "guide"), help="Type of gprobe"),
                                        keys.Key("oneExposure", help="Take just one exposure"),
                                        keys.Key("Kp", types.Float(), help="Proportional gain"),
                                        keys.Key("Ti", types.Float(), help="Integral time"),
                                        keys.Key("Td", types.Float(), help="Derivative time"),
                                        keys.Key("Imax", types.Float(), help="|maximum value of I| (-ve to disable)"),
                                        keys.Key("geek", help="Show things that only some of us love"),
                                        keys.Key("plot", help="Plot things"),
                                        )
        #
        # Declare commands
        #
        self.vocab = [
            ("on", "[<time>] [force] [oneExposure] [plot]", self.guideOn),
            ("off", "", self.guideOff),
            ("setExpTime", "<time>", self.setExpTime),
            ("setPID", "(azAlt|rot|focus|scale) <Kp> [<Ti>] [<Td>] [<Imax>]", self.setPID),
            ("disable", "<fibers>|<gprobes>", self.disableFibers),
            ("enable", "<fibers>|<gprobes>", self.enableFibers),
            ("loadCartridge", "<cartridge> [<pointing>]", self.loadCartridge),
            ("showCartridge", "", self.showCartridge),
            ('ping', '', self.ping),
            ('restart', '', self.restart),
            ('axes', '(on|off)', self.axes),
            ('focus', '(on|off)', self.focus),
            ('scale', '(on|off)', self.scale),
            ('status', "[geek]", self.status),
            ]
    #
    # Define commands' callbacks
    #
    def disableFibersImpl(self, cmd, enable=True):
        """Disable a set of fibers"""

        if "fibers" in cmd.cmd.keywords:
            for f in cmd.cmd.keywords["fibers"].values:
                myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.ENABLE_FIBER,
                                                                        cmd=cmd, fiber=f, enable=enable))
        elif "gprobes" in cmd.cmd.keywords:
            gprobeType = cmd.cmd.keywords["gprobes"].values[0].upper()
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.ENABLE_FIBER, cmd=cmd,
                                                                    fiber=gprobeType, enable=enable))

        self.status(cmd)                # finishes this command

    def disableFibers(self, cmd):
        """Disable a set of fibers"""

        self.disableFibersImpl(cmd, enable=False)

    def enableFibers(self, cmd):
        """Enable a set of fibers"""

        self.disableFibersImpl(cmd, enable=True)

    def setExpTime(self, cmd):
        """Set the exposure time"""

        expTime = cmd.cmd.keywords["time"].values[0]
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, cmd=cmd, expTime=expTime))

    def setPID(self, cmd):
        """Set something's PID coefficients"""

        what = None
        for k in ["azAlt", "rot", "focus", "scale"]:
            if k in cmd.cmd.keywords:
                what = k
                break

        if not what:
            cmd.fail("text=\"Impossible condition in setPID\"")

        Kp = cmd.cmd.keywords["Kp"].values[0]
        Ti = cmd.cmd.keywords["Ti"].values[0] if "Ti" in cmd.cmd.keywords else 0
        Td = cmd.cmd.keywords["Td"].values[0] if "Td" in cmd.cmd.keywords else 0
        Imax = cmd.cmd.keywords["Imax"].values[0] if "Imax" in cmd.cmd.keywords else 0

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_PID, cmd=cmd, what=what,
                                                                Kp=Kp, Ti=Ti, Td=Td, Imax=Imax))

    def guideOn(self, cmd):
        """Turn guiding on"""

        force = "force" in cmd.cmd.keywords
        oneExposure = "oneExposure" in cmd.cmd.keywords
        plot = "plot" in cmd.cmd.keywords
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else None
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd,
                                                                start=True, expTime=expTime,
                                                                force=force, oneExposure=oneExposure,
                                                                plot=plot))
    def guideOff(self, cmd):
        """Turn guiding off"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd, start=False))

    def loadCartridge(self, cmd):
        """Load a cartridge"""

        cartridge = cmd.cmd.keywords["cartridge"].values[0]
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        #
        # Cartridge ID of 0 means that no cartridge is loaded
        #
        if cartridge == 0:
            gprobes = {}
            plate = 0
            boresight_ra = float("NaN")
            boresight_dec = float("NaN")
            #
            # Send that information off to the master thread
            #
            myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                                                                cartridge=cartridge, plate=plate, pointing=pointing,
                                                                boresight_ra=boresight_ra, boresight_dec=boresight_dec,
                                                                gprobes=gprobes))
            return
        #
        # Get the plate from the plateDB
        #
        actorState = guiderActor.myGlobals.actorState

        pointingInfoKey = actorState.models["platedb"].keyVarDict["pointingInfo"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="loadCartridge cartridge=%d pointing=%s" % (cartridge, pointing),
                                            keyVars=[pointingInfoKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup plate corresponding to %d/%s\"" % (cartridge, pointing))
            return

        plate = cmdVar.getLastKeyVarData(pointingInfoKey)[0]
        boresight_ra = cmdVar.getLastKeyVarData(pointingInfoKey)[3]
        boresight_dec = cmdVar.getLastKeyVarData(pointingInfoKey)[4]
        #
        # Lookup the valid gprobes
        #
        gprobeKey = actorState.models["platedb"].keyVarDict["gprobe"]
        gprobesInUseKey = actorState.models["platedb"].keyVarDict["gprobesInUse"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
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

        gprobes = {}
        for cartridgeID, gpID, exists, xCenter, yCenter, radius, rotation, \
                xFerruleOffset, yFerruleOffset, focusOffset, fiber_type in cmdVar.getKeyVarData(gprobeKey):
            gprobes[gpID] = GuiderCmd.GprobeInfo(exists, enabled.get(gpID, False), xCenter, yCenter, radius, rotation,
                                                 xFerruleOffset, yFerruleOffset, focusOffset, fiber_type)

        #
        # Add in the plate/fibre geometry from plPlugMapM
        #
        guideInfoKey = actorState.models["platedb"].keyVarDict["guideInfo"]
        plPlugMapMKey = actorState.models["platedb"].keyVarDict["plPlugMapM"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="getGprobesPlateGeom cartridge=%d" % (cartridge),
                                            keyVars=[guideInfoKey, plPlugMapMKey])
        if cmdVar.didFail:
            cmd.fail("text=%s" % qstr("Failed to lookup gprobes's geometry for cartridge %d" % (cartridge)))
            return

        assert int(cmdVar.getLastKeyVarData(plPlugMapMKey)[0]) == plate
        fscanMJD = cmdVar.getLastKeyVarData(plPlugMapMKey)[1]
        fscanID = cmdVar.getLastKeyVarData(plPlugMapMKey)[2]

        for el in cmdVar.getKeyVarData(guideInfoKey):
            i = 0
            id = int(el[0]); i += 1
            
            try:
                gprobes[id].ra = float(el[i]); i += 1
                gprobes[id].dec = float(el[i]); i += 1
                gprobes[id].xFocal = float(el[i]); i += 1
                gprobes[id].yFocal = float(el[i]); i += 1
                gprobes[id].phi = float(el[i]); i += 1
                gprobes[id].throughput = float(el[i]); i += 1

            except KeyError:
                cmd.warn("text=\"Unknown fiberId %d from plugmap file (%s)\"" % (id, ", ".join(el[1:])))
                continue
        #
        # Send that information off to the master thread
        #
        myGlobals.actorState.queues[guiderActor.MASTER].put(
            Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                cartridge=cartridge, plate=plate, pointing=pointing,
                fscanMJD=fscanMJD, fscanID=fscanID,
                boresight_ra=boresight_ra, boresight_dec=boresight_dec,
                gprobes=gprobes))

    def ping(self, cmd):
        """ Top-level 'ping' command handler. Query the actor for liveness/happiness. """

        cmd.finish('text="pong"')

    def restart(self, cmd):
        """Restart the worker threads"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd, start=None))

        actorState = myGlobals.actorState

        if actorState.restartCmd:
            actorState.restartCmd.finish("text=\"Nunc dimittis servum tuum Domine\"")
            actorState.restartCmd = None

        actorState.actor.startThreads(actorState, cmd, restart=True)
        #
        # We can't finish this command now as the threads may not have died yet,
        # but we can remember to clean up _next_ time we restart
        #
        cmd.inform("text=\"Restarting threads\"")
        actorState.restartCmd = cmd

    def correctionImpl(self, cmd, what):
        """Turn guiding something on or off"""

        on = "on" in cmd.cmd.keywords
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_GUIDE_MODE, cmd=cmd, what=what, enable=on))

    def axes(self, cmd):
        """Turn guiding the plate axes on or off"""

        self.correctionImpl(cmd, "axes")

    def focus(self, cmd):
        """Turn guiding the plate focus on or off"""

        self.correctionImpl(cmd, "focus")

    def scale(self, cmd):
        """Turn guiding the plate scale on or off"""

        self.correctionImpl(cmd, "scale")

    def showCartridge(self, cmd, full=True):
        """Reveal the identity of the current cartridge"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STATUS, cmd=cmd, full=False, finish=True))

    def status(self, cmd, full=True):
        """Return guide status status"""

        if "geek" in cmd.cmd.keywords:
            for t in threading.enumerate():
                cmd.inform('text="%s"' % t)

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STATUS, cmd=cmd, finish=True))
