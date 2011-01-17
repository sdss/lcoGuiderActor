#!/usr/bin/env python

""" Wrap top-level ICC functions. """

import pdb
import logging
import os, re, sys
import numpy
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

    GOOD   =  0x0                       # N.b. these are repeated in PlatedbCmd.py. Caveat editor
    BROKEN =  0x1                       # Worse, they are repeated in 
    NOSTAR =  0x2
    DISABLE = 0x4
    UNKNOWN = 0xff                      # shouldn't ever happen.

    class GprobeInfo(object):
        """Capture information about a guider probe"""
        def __init__(self, exists, enabled, xCenter, yCenter, radius, rotation,
                     xFerruleOffset, yFerruleOffset, focusOffset, fiber_type, flags):
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
            self.rotStar2Sky = numpy.nan
            self.flags = flags

    def __init__(self, actor):
        self.actor = actor
        #
        # Declare keys that we're going to use
        #
        self.keys = keys.KeysDictionary("guider_guider", (2, 1),
                                        keys.Key("cartridge", types.Int(), help="A cartridge ID"),
                                        keys.Key("fscanId", types.Int(), help="The fscanId identifying a plate scanning"),
                                        keys.Key("mjd", types.Int(), help="The MJD when a plate was scanned"),
                                        keys.Key("plate", types.Int(), help="A plugplate ID"),
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
                                        keys.Key("nfilt", types.Int(), help="number of input readings to filter with."),
                                        keys.Key("geek", help="Show things that only some of us love"),
                                        keys.Key("plot", help="Plot things"),
                                        keys.Key("psPlot", help="Save copy of plots"),
                                        keys.Key("display", types.String(), help="DISPLAY variable to use"),
                                        keys.Key("spiderInstAng", types.Float(), help="Value to use for spiderInstAng"),
                                        keys.Key("cartfile", types.String(), help="cartridge file"),
                                        keys.Key("plugfile", types.String(), help="plugmap file"),
                                        keys.Key("file", types.String(), help="guider file"),
                                        keys.Key("decenter", types.String(), help="apply decenter offsets to guiding"),
                                        keys.Key("decenterRA", types.Float(), help="Telescope absolute offset for guiding in RA arcsec"),
                                        keys.Key("decenterDec", types.Float(), help="Telescope absolute offset for guiding in Dec arcsec"),
                                        keys.Key("decenterRot", types.Float(), help="Telescope absolute offset for guiding in Rot"),
                                        keys.Key("scale", types.Float(), help="Current scale from \"tcc show scale\""),
                                        keys.Key("delta", types.Float(), help="Delta scale (percent)"),
                                       )
        #
        # Declare commands
        #
        self.vocab = [
            ("on", "[<time>] [force] [oneExposure] [decenter] [plot] [psPlot] [<display>] [<spiderInstAng>]", self.guideOn),
            ("off", "", self.guideOff),
            ("setExpTime", "<time>", self.setExpTime),
            ("setPID", "(raDec|rot|focus|scale) <Kp> [<Ti>] [<Td>] [<Imax>] [nfilt]", self.setPID),
            ("disable", "<fibers>|<gprobes>", self.disableFibers),
            ("enable", "<fibers>|<gprobes>", self.enableFibers),
            ("loadCartridge", "[<cartridge>] [<pointing>] [<plate>] [<mjd>] [<fscanId>] [force]", self.loadCartridge),
            ("showCartridge", "", self.showCartridge),
            ("loadPlateFiles", "<cartfile> <plugfile>", self.loadPlateFiles),
            ("reprocessFile", "<file>", self.reprocessFile),
            ("flat", "[<time>]", self.flat),
            ('ping', '', self.ping),
            ('restart', '', self.restart),
            ('axes', '(on|off)', self.axes),
            ('focus', '(on|off)', self.focus),
            ('scale', '(on|off)', self.scale),
            ('status', "[geek]", self.status),
            ('centerUp', "", self.centerUp),
            ("setScale", "<delta>|<scale>", self.setScale),
            ("scaleChange", "<delta>|<scale>", self.scaleChange),
            ('setDecenter', "[<decenterRA>] [<decenterDec>] [<decenterRot>]", self.setDecenter),
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

    def flat(self, cmd):
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else 0.5
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.TAKE_FLAT, cmd=cmd,
                                                                start=True, expTime=expTime))

    def setPID(self, cmd):
        """Set something's PID coefficients"""

        what = None
        for k in ["raDec", "rot", "focus", "scale"]:
            if k in cmd.cmd.keywords:
                what = k
                break

        if not what:
            cmd.fail("text=\"Impossible condition in setPID\"")

        Kp = cmd.cmd.keywords["Kp"].values[0]
        Ti = cmd.cmd.keywords["Ti"].values[0] if "Ti" in cmd.cmd.keywords else 0
        Td = cmd.cmd.keywords["Td"].values[0] if "Td" in cmd.cmd.keywords else 0
        Imax = cmd.cmd.keywords["Imax"].values[0] if "Imax" in cmd.cmd.keywords else 0
        nfilt = cmd.cmd.keywords["nfilt"].values[0] if "nfilt" in cmd.cmd.keywords else 0

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_PID, cmd=cmd, what=what,
                                                                Kp=Kp, Ti=Ti, Td=Td, Imax=Imax,
                                                                nfilt=nfilt))

    def guideOn(self, cmd):
        """Turn guiding on"""

        force = "force" in cmd.cmd.keywords
        oneExposure = "oneExposure" in cmd.cmd.keywords
        plot = "plot" in cmd.cmd.keywords
        psPlot = "psPlot" in cmd.cmd.keywords
        display = cmd.cmd.keywords["display"].values[0] if "display" in cmd.cmd.keywords else None
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else None
        spiderInstAng = cmd.cmd.keywords["spiderInstAng"].values[0] if "spiderInstAng" in cmd.cmd.keywords else None
        decenter = True if "decenter" in cmd.cmd.keywords else None

        if spiderInstAng and not force:
            cmd.fail('text="You may only specify a spiderInstAng with force, for debugging during the day"')
            return

        if decenter and not force:
            cmd.fail('text="You must specify decenter with force, its for guider engineering test only "')
            return

        if display:
            os.environ["DISPLAY"] = display

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd,
                                                                start=True, expTime=expTime,
                                                                spiderInstAng=spiderInstAng,
                                                                force=force, oneExposure=oneExposure,
                                                                plot=plot, psPlot=psPlot, decenter=decenter))
    def guideOff(self, cmd):
        """Turn guiding off"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd, start=False))

    def centerUp(self, cmd):
        """Force a single XY offset"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.CENTERUP, cmd=cmd))

    def reprocessFile(self, cmd):
        """Reprocess a single file."""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.REPROCESS_FILE, cmd=cmd, 
                                                                filename=cmd.cmd.keywords["filename"].values[0]))

    def loadPlateFiles(self, cmd):
        """Read in cartridge and plugmap files. """

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.READ_PLATE_FILES, cmd=cmd,
                                                                plugfile=cmd.cmd.keywords["plugfile"].values[0],
                                                                cartfile=cmd.cmd.keywords["cartfile"].values[0]))

    def loadCartridge(self, cmd):
        """Load a cartridge

If the cartridge ID is omitted the currently-mounted cartridge is used.  It's an error to specify a cartridge
that isn't actually mounted (unless you specify force)
"""

        force = "force" in cmd.cmd.keywords
        cartridge = cmd.cmd.keywords["cartridge"].values[0] if "cartridge" in cmd.cmd.keywords else -1
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        #
        # If they specify a plate explicitly, we'll bypass the active table and give them what they want
        #
        plate = str(cmd.cmd.keywords["plate"].values[0]) if "plate" in cmd.cmd.keywords else None
        mjd = cmd.cmd.keywords["mjd"].values[0] if "mjd" in cmd.cmd.keywords else None
        fscanId = cmd.cmd.keywords["fscanId"].values[0] if "fscanId" in cmd.cmd.keywords else None
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
        # Check that the claimed cartridge is actually on the telescope
        #
        actorState = guiderActor.myGlobals.actorState

        instrumentNumKey = actorState.models["mcp"].keyVarDict["instrumentNum"]
        cmdVar = actorState.actor.cmdr.call(actor="mcp", forUserCmd=cmd,
                                            cmdStr="info", keyVars=[instrumentNumKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to ask mcp for info on cartridges\"")
            return

        loadedCartridge = cmdVar.getLastKeyVarData(instrumentNumKey)[0]
        cmd.inform("text=\"Cartridge %s is on the telescope\"" % loadedCartridge)

        if cartridge < 0:
            cartridge = loadedCartridge
            
        if loadedCartridge != cartridge:
            msg = "Expected cartridge %s, but %s is loaded" % (cartridge, loadedCartridge)
            if force:
                cmd.warn("text=\"%s\"" % (msg + "; proceeding"))
            else:
                cmd.fail("text=\"%s\"" % msg)
                return
        #
        # Get the plate from the plateDB
        #
        pointingInfoKey = actorState.models["platedb"].keyVarDict["pointingInfo"]
        extraArgs = ""
        if plate: extraArgs += " plate=%s" % (plate)
        if mjd: extraArgs += " mjd=%s" % (mjd)
        if fscanId: extraArgs += " fscanId=%s" % (fscanId)
        
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="loadCartridge cartridge=%d pointing=%s %s" % \
                                                (cartridge, pointing, extraArgs),
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
        extraArgs = ""
        if plate: extraArgs += "plate=%s" % (plate)
        gprobeKey = actorState.models["platedb"].keyVarDict["gprobe"]
        gprobesInUseKey = actorState.models["platedb"].keyVarDict["gprobesInUse"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="getGprobes cartridge=%d pointing=%s %s" % \
                                                (cartridge, pointing, extraArgs),
                                            keyVars=[gprobeKey, gprobesInUseKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup gprobes for cartridge %d\"" % (cartridge))
            return

        enabled = {}; flags = {}
        for el in cmdVar.getLastKeyVarData(gprobesInUseKey):
            mat = re.search(r"^\((\d+)\s*=\s*(\S+)\s*\)$", el)
            id, flags[id] = int(mat.group(1)), int(mat.group(2), 16)

            enabled[id] = True if flags[id] == GuiderCmd.GOOD else False

        gprobes = {}
        for cartridgeID, gpID, exists, xCenter, yCenter, radius, rotation, \
                xFerruleOffset, yFerruleOffset, focusOffset, fiber_type in cmdVar.getKeyVarData(gprobeKey):
            gprobes[gpID] = GuiderCmd.GprobeInfo(exists, enabled.get(gpID, False), xCenter, yCenter, radius,
                                                 rotation, xFerruleOffset, yFerruleOffset, focusOffset,
                                                 fiber_type, 
                                                 flags.get(gpID, GuiderCmd.UNKNOWN) | GuiderCmd.NOSTAR)

        #
        # Add in the plate/fibre geometry from plPlugMapM
        #
        guideInfoKey = actorState.models["platedb"].keyVarDict["guideInfo"]
        plPlugMapMKey = actorState.models["platedb"].keyVarDict["plPlugMapM"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="getGprobesPlateGeom cartridge=%d %s" % (cartridge, extraArgs),
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

            if id < 0:                  # invalid; typically -9999
                continue
            
            try:
                gprobes[id].ra = float(el[i]); i += 1
                gprobes[id].dec = float(el[i]); i += 1
                gprobes[id].xFocal = float(el[i]); i += 1
                gprobes[id].yFocal = float(el[i]); i += 1
                gprobes[id].phi = float(el[i]); i += 1
                gprobes[id].throughput = float(el[i]); i += 1
                gprobes[id].flags &= ~GuiderCmd.NOSTAR

            except KeyError:
                cmd.warn("text=\"Unknown fiberId %d from plugmap file (%s)\"" % (id, ", ".join([str(e) for e in el[1:]])))
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

        self.actor.sendVersionKey(cmd)
        cmd.finish('text="pong"')

    def restart(self, cmd):
        """Restart the worker threads"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd, start=None))

        actorState = myGlobals.actorState

        if actorState.restartCmd:
            actorState.restartCmd.finish("text=\"Nunc dimittis servum tuum Domine\"")
            actorState.restartCmd = None

        actorState.actor.startThreads(actorState, cmd, restart=True, restartQueues=True)
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

    def scaleChange(self, cmd):
        """Alias for setScale """
        self.setScale(cmd)

    def setScale(self, cmd):
        """Change telescope scale by a factor of (1 + 0.01*delta), or to scale """

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.CHANGE_SCALE, cmd=cmd, finish=True))

    def status(self, cmd, full=True):
        """Return guide status status"""

        self.actor.sendVersionKey(cmd)
        if "geek" in cmd.cmd.keywords:
            for t in threading.enumerate():
                cmd.inform('text="%s"' % t)

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STATUS, cmd=cmd, finish=True))

    def setDecenter(self, cmd):
        """guide at an absolute offset location"""
        #require elsewhere that guiding be run with force to allow decentered guiding
        #for now Decenter rot is around (RA+decenterRA, Dec+decenterDec)
        decenterRA  = cmd.cmd.keywords["decenterRA"].values[0] if "decenterRA" in cmd.cmd.keywords else 0
        decenterDec = cmd.cmd.keywords["decenterDec"].values[0] if "decenterDEC" in cmd.cmd.keywords else 0
        decenterRot = cmd.cmd.keywords["decenterRot"].values[0] if "decenterRot" in cmd.cmd.keywords else 0
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.DECENTER, cmd=cmd, decenterRA=decenterRA,
                                                                decenterDec=decenterDec, decenterRot=decenterRot))
