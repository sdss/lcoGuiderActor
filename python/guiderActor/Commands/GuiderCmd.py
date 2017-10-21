#!/usr/bin/env python
"""Define the available commands and how to handle them."""

import os
import threading

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr
import opscore.utility.YPF as YPF

from guiderActor import Msg, GuiderState
import guiderActor
import guiderActor.myGlobals as myGlobals

from fitPlateProfile import DuPontMeasurement, DuPontProfile, plt#, MeasRadii, DirThetaMapDuPont


import numpy as np


class GuiderCmd(object):
    """ Wrap commands to the guider actor"""
    def __init__(self, actor):
        """
        Declares keys that this actor uses, and available commands that can be sent to it.

        actor is the actor that this is part of (guiderActor, in this case).
        """
        self.actor = actor
        #
        # Declare keys that we're going to use
        #
        self.keys = keys.KeysDictionary(
            "guider_guider", (3, 1),
            keys.Key("cartridge", types.Int(), help="A cartridge ID"),
            keys.Key("fscanId", types.Int(), help="The fscanId identifying a plate scanning"),
            keys.Key("mjd", types.Int(), help="The MJD when a plate was scanned"),
            keys.Key("plate", types.Int(), help="A plugplate ID"),
            keys.Key("fibers", types.Int()*(1, None), help="A list of fibers"),
            keys.Key("probe", types.Int(), help="A probe ID, 1-indexed"),
            keys.Key("gprobe", types.Int(), help="A probe ID, 1-indexed"),
            keys.Key("fromProbe", types.Int(), help="A probe ID, 1-indexed"),
            keys.Key("fromGprobe", types.Int(), help="A probe ID, 1-indexed"),
            keys.Key("pointing", types.String(), help="A pointing for the given plugplate"),
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
            keys.Key("cartfile", types.String(), help="cartridge file"),
            keys.Key("plugfile", types.String(), help="plugmap file"),
            keys.Key("file", types.String(), help="guider file"),
            keys.Key("decenterRA", types.Float(), help="Telescope absolute offset for guiding in RA arcsec"),
            keys.Key("decenterDec", types.Float(), help="Telescope absolute offset for guiding in Dec arcsec"),
            keys.Key("decenterRot", types.Float(), help="Telescope absolute offset for guiding in Rot"),
            keys.Key("ditherPos", types.String(),  help="Named MaNGA guider dither position"),
            keys.Key("scale", types.Float(), help="Current scale from \"tcc show scale\""),
            keys.Key("delta", types.Float(), help="Delta scale (percent)"),
            keys.Key("stack", types.Int(), help="number of itime gcamera integrations to request per exposure."),
            keys.Key("corrRatio", types.Float(), help="How much refraction correction to apply (0..)"),
            keys.Key("plateType", types.String(), help="Name of the current plateType (survey concatenation)"),
            keys.Key("surveyMode", types.String(), help="Name of the current surveyMode"),
            keys.Key("movieMJD", types.String(), help="The MJD that we want to generate the movie for."),
            keys.Key("start", types.Int(), help="Guider frame number to start the movie at."),
            keys.Key("end", types.Int(), help="Guider frame number to end the movie at."),
            keys.Key("bin", types.Int(), help="bin factor for exposure"),
            keys.Key('gprobeX', types.Int(), help='fake gprobe x position (binned)'),
            keys.Key('gprobeY', types.Int(), help='fake gprobe y position (binned)'),
            keys.Key('gprobeRadius', types.Int(), help='fake gprobe radius (binned)'),
            keys.Key('gprobeRA', types.Float(), help='fake gprobe RA'),
            keys.Key('gprobeDec', types.Float(), help='fake gprobe declination'))

        #
        # Declare commands
        #
        self.vocab = [
            ("on", "[<time>] [force] [oneExposure] [<stack>]", self.guideOn),
            ("off", "", self.guideOff),
            ("setExpTime", "<time> [<stack>]", self.setExpTime),
            ("setPID", "(raDec|rot|focus|scale) <Kp> [<Ti>] [<Td>] [<Imax>] [nfilt]", self.setPID),
            ("disable", "<fibers>|<gprobes>", self.disableFibers),
            ("enable", "<fibers>|<gprobes>", self.enableFibers),
            # fake this, get info directly from plPlugMap file
            ("loadCartridge", "[<cartridge>] [<pointing>] [<plate>] [<mjd>] [<fscanId>] [force]", self.loadCartridge),
            ("showCartridge", "", self.showCartridge),
            ("loadPlateFiles", "<cartfile> <plugfile>", self.loadPlateFiles),
            ("reprocessFile", "<file>", self.reprocessFile),
            ('bias', '[<stack>]', self.bias),
            ("flat", "[<time>]", self.flat),
            ("dark", "[<time>] [<stack>]", self.dark),
            ('ping', '', self.ping),
            ('restart', '', self.restart),
            ('axes', '(on|off)', self.axes),
            ('focus', '(on|off)', self.focus),
            ('scale', '(on|off)', self.scale),
            ('status', "[geek]", self.status),
            ('centerUp', "", self.centerUp),
            ('fk5InFiber', "[<probe>] [<time>]", self.fk5InFiber),
            ('starInFiber', "[<probe>] [<gprobe>] [<fromProbe>] [<fromGprobe>]", self.starInFiber),
            ("setScale", "<delta>|<scale>", self.setScale),
            ("scaleChange", "<delta>|<scale>", self.scaleChange),
            ('decenter', '(on|off)', self.decenter),
            ('setDecenter', "[<decenterRA>] [<decenterDec>] [<decenterRot>]", self.setDecenter),
            ('mangaDither', "<ditherPos>", self.mangaDither),
            ('setRefractionBalance', "[<corrRatio>] [<plateType>] [<surveyMode>]", self.setRefractionBalance),
            ('makeMovie', '[<movieMJD>] <start> <end>', self.makeMovie),
            ('findstar', '[<time>] [<bin>]', self.ecam_findstar),
            ('sendfield', '', self.sendfield),
            ('loadFakeGProbe', '[<gprobeX>] [<gprobeY>] [<gprobeRadius>] [<gprobeRA>] [<gprobeDec>]', self.loadFakeGProbe),
            ('oneStep', '[<time>]', self.oneStep)
            ]

        self.fakeCartNumber = 99

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
        """Disable a set of fibers, by probe number."""

        self.disableFibersImpl(cmd, enable=False)

    def enableFibers(self, cmd):
        """Enable a set of fibers, by probe number."""

        self.disableFibersImpl(cmd, enable=True)

    def setExpTime(self, cmd):
        """Set the exposure time"""

        expTime = cmd.cmd.keywords["time"].values[0]
        stack = cmd.cmd.keywords["stack"].values[0] if "stack" in cmd.cmd.keywords else 1
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_TIME, cmd=cmd, expTime=expTime, stack=stack))

    def flat(self, cmd):
        """Take, and process, a guider flat."""
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else 0.5
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.TAKE_FLAT, cmd=cmd,
                                                                expTime=expTime))

    def bias(self, cmd):
        """Takes a stack of bias images."""

        stack = (cmd.cmd.keywords['stack'].values[0]
                 if 'stack' in cmd.cmd.keywords else 7)

        if stack < 5:
            cmd.fail('text="bias requires a minimum stack of 5 images."')
            return

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.TAKE_BIAS,
                                                                cmd=cmd,
                                                                stack=stack))

    def dark(self, cmd):
        """
        Take, and process, a guider dark.
        Recommended: guider dark time=15 stack=9
        Minimum for a processed dark: guider dark time=10 stack=5
        """
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else 15
        stack = cmd.cmd.keywords["stack"].values[0] if "stack" in cmd.cmd.keywords else 3
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.TAKE_DARK, cmd=cmd,
                                                                expTime=expTime, stack=stack))

    def setPID(self, cmd):
        """Set something's PID coefficients"""

        axis = None
        for k in ["raDec", "rot", "focus", "scale"]:
            if k in cmd.cmd.keywords:
                axis = k
                break

        if not axis:
            cmd.fail("text=\"Impossible condition in setPID\"")

        # TODO: should Imax and nfilt default to zero? (JSG)
        Kp = cmd.cmd.keywords["Kp"].values[0]
        Ti = cmd.cmd.keywords["Ti"].values[0] if "Ti" in cmd.cmd.keywords else 0
        Td = cmd.cmd.keywords["Td"].values[0] if "Td" in cmd.cmd.keywords else 0
        Imax = cmd.cmd.keywords["Imax"].values[0] if "Imax" in cmd.cmd.keywords else 0
        nfilt = cmd.cmd.keywords["nfilt"].values[0] if "nfilt" in cmd.cmd.keywords else 0
        ncorr = cmd.cmd.keywords['ncorr'].values[0] if 'ncorr' in cmd.cmd.keywords else 1

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_PID, cmd=cmd, axis=axis,
                                                                Kp=Kp, Ti=Ti, Td=Td, Imax=Imax,
                                                                nfilt=nfilt, ncorr=ncorr))

    def guideOn(self, cmd):
        """Turn guiding on"""

        force = "force" in cmd.cmd.keywords
        oneExposure = "oneExposure" in cmd.cmd.keywords
        expTime = cmd.cmd.keywords["time"].values[0] if "time" in cmd.cmd.keywords else None
        stack = cmd.cmd.keywords["stack"].values[0] if "stack" in cmd.cmd.keywords else 1
        camera = 'ecamera' if myGlobals.actorState.gState.plateType == 'ecamera' else 'gcamera'

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd,
                                                                expTime=expTime, stack=stack, camera=camera,
                                                                force=force, oneExposure=oneExposure))

    def guideOff(self, cmd):
        """Turn guiding off"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STOP_GUIDING, cmd=cmd))

    def centerUp(self, cmd):
        """Force a single XY offset"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.CENTERUP, cmd=cmd))

    def fk5InFiber(self, cmd):
        """Have the TCC put a bright star down a given probe"""

        actorState = guiderActor.myGlobals.actorState
        probe = cmd.cmd.keywords['probe'].values[0] if 'probe' in cmd.cmd.keywords else None
        expTime = cmd.cmd.keywords["time"].values[0] if 'time' in cmd.cmd.keywords else 0.1
        stack = cmd.cmd.keywords["stack"].values[0] if 'stack' in cmd.cmd.keywords else 1

        # Force up an image-only guide loop
        for what in ["scale", "focus", "axes"]:
            actorState.queues[guiderActor.MASTER].put(Msg(Msg.SET_GUIDE_MODE, cmd=cmd, what=what, enable=False))
        actorState.queues[guiderActor.MASTER].put(Msg(Msg.START_GUIDING, cmd=cmd, oneExposure=False,
                                                      expTime=expTime, stack=stack, force=True))
        if probe:
            cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                                cmdStr="set ptErrProbe=%d" % (probe))
            if cmdVar.didFail:
                cmd.fail("text=\"Failed to set the pointing error probe to %s\"" % (probe))
                return

        cmdVar = actorState.actor.cmdr.call(actor="tcc", forUserCmd=cmd,
                                            cmdStr="track/pterr")
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to move to a bright star\"")
            return

        cmd.finish("text='There should be a bright star in probe'")

    def loadAllProbes(self, cmd):
        pass

    def starInFiber(self, cmd):
        """ Put a star down a given probe """

        probe = cmd.cmd.keywords['probe'].values[0] if 'probe' in cmd.cmd.keywords else None
        gprobe = cmd.cmd.keywords['gprobe'].values[0] if 'gprobe' in cmd.cmd.keywords else None
        if (probe is None and gprobe is None) or (probe is not None and gprobe is not None):
            cmd.fail('text="exactly one destination probe must specified"')
            return

        fromProbe = cmd.cmd.keywords["fromProbe"].values[0] if 'fromProbe' in cmd.cmd.keywords else None
        fromGprobe = cmd.cmd.keywords["fromGprobe"].values[0] if 'fromGprobe' in cmd.cmd.keywords else None
        if (fromProbe is not None and fromGprobe is not None):
            cmd.fail('text="no more than one source probe can be specified"')
            return

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STAR_IN_FIBER, cmd=cmd,
                                                                probe=probe, gprobe=gprobe,
                                                                fromProbe=fromProbe, fromGprobe=fromGprobe))

    def reprocessFile(self, cmd):
        """Reprocess a single file."""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.REPROCESS_FILE, cmd=cmd,
                                                                filename=cmd.cmd.keywords["file"].values[0]))

    def loadPlateFiles(self, cmd):
        """Read in cartridge and plugmap files. """

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.READ_PLATE_FILES, cmd=cmd,
                                                                plugfile=cmd.cmd.keywords["plugfile"].values[0],
                                                                cartfile=cmd.cmd.keywords["cartfile"].values[0]))

    def loadCartridge(self, cmd):
        """
        Load a cartridge.
        If the cartridge ID is omitted the currently-mounted cartridge is used.
        Error if cartridge that isn't actually mounted is specified (unless force is also given).
        """

        queue = myGlobals.actorState.queues[guiderActor.MASTER]

        force = "force" in cmd.cmd.keywords
        cartridge = cmd.cmd.keywords["cartridge"].values[0] if "cartridge" in cmd.cmd.keywords else -1
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        #
        # If they specify a plate explicitly, we'll bypass the active table and give them what they want
        #
        plate = str(cmd.cmd.keywords["plate"].values[0]) if "plate" in cmd.cmd.keywords else None
        mjd = cmd.cmd.keywords["mjd"].values[0] if "mjd" in cmd.cmd.keywords else None
        fscanId = cmd.cmd.keywords["fscanId"].values[0] if "fscanId" in cmd.cmd.keywords else None

        # Cartridge ID of 0 means that no cartridge is loaded
        if cartridge == 0:
            gprobes = {}
            plate = 0
            boresight_ra = float("NaN")
            boresight_dec = float("NaN")
            design_ha = ("NaN")
            # Send that information off to the master thread
            queue.put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                      cartridge=cartridge, plate=plate, pointing=pointing,
                      boresight_ra=boresight_ra, boresight_dec=boresight_dec,
                      design_ha=design_ha,
                      gprobes=gprobes))
            return
        #
        # Check that the claimed cartridge is actually on the telescope
        #
        actorState = guiderActor.myGlobals.actorState

        loadedCartridge = self.actor.getLoadedCartridge(cmd, actorState)
        if loadedCartridge is None:
            if not force:
                cmd.fail('text="failed retrieving the number of the cart on the telescope."')
                return
            else:
                cmd.warn('text="failed retrieving the number of the cart on the telescope. Proceeding with force."')
        else:
            cmd.inform("text=\"Cartridge %s is on the telescope\"" % loadedCartridge)

        # Only auto-select the cart if a plate was not specified.
        if cartridge is not None and cartridge < 0 and plate is None:
            cartridge = loadedCartridge

        if loadedCartridge != cartridge:
            msg = "Expected cartridge %s, but %s is loaded" % (cartridge, loadedCartridge)
            if force:
                cmd.warn("text=\"%s\"" % (msg + "; proceeding"))
            else:
                cmd.fail("text=\"%s\"" % msg)
                return

        # cart 19 is the engineering camera, and has no info in platedb.
        if cartridge == 19:
            # don't do anything but clear the gprobes and output status.
            gState = actorState.gState
            gState.deleteAllGprobes()
            gState.cartridge = cartridge
            gState.plate = 0
            gState.pointing = pointing
            gState.plateType = 'ecamera'
            gState.surveyMode = None
            queue.put(Msg(Msg.STATUS, cmd, finish=True))
            return

        # Get the plate from the plateDB
        pointingInfoKey = actorState.models["platedb"].keyVarDict["pointingInfo"]
        extraArgs = ""
        if plate:
            extraArgs += " plate=%s" % (plate)
        if mjd:
            extraArgs += " mjd=%s" % (mjd)
        if fscanId:
            extraArgs += " fscanId=%s" % (fscanId)

        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="loadCartridge cartridge=%d pointing=%s %s" %
                                            (cartridge, pointing, extraArgs),
                                            keyVars=[pointingInfoKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup plate corresponding to %d/%s\"" % (cartridge, pointing))
            return

        plate = cmdVar.getLastKeyVarData(pointingInfoKey)[0]
        boresight_ra = cmdVar.getLastKeyVarData(pointingInfoKey)[3]
        boresight_dec = cmdVar.getLastKeyVarData(pointingInfoKey)[4]
        design_ha = cmdVar.getLastKeyVarData(pointingInfoKey)[5]
        survey = cmdVar.getLastKeyVarData(pointingInfoKey)[8]
        surveyMode = cmdVar.getLastKeyVarData(pointingInfoKey)[9]
        if design_ha < 0:
            design_ha += 360

        # Lookup the valid gprobes
        gprobeKey = actorState.models["platedb"].keyVarDict["gprobe"]
        gprobesInUseKey = actorState.models["platedb"].keyVarDict["gprobesInUse"]
        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="getGprobes cartridge=%d pointing=%s %s" %
                                            (cartridge, pointing, extraArgs),
                                            keyVars=[gprobeKey, gprobesInUseKey])
        if cmdVar.didFail:
            cmd.fail("text=\"Failed to lookup gprobes for cartridge %d\"" % (cartridge))
            return

        # Unpack the various platedb guider keys into a Probe instance for each probe
        # NOTE: ordered so that we first set the gprobebits, then fill in the rest of the values.
        # as otherwise the gprobebits would overwrite some of the state we set.
        gprobes = {}
        for key in cmdVar.getLastKeyVarData(gprobesInUseKey):
            probeId, flags = key.strip('()').split('=')
            gprobes[int(probeId)] = GuiderState.GProbe(int(probeId))
            gprobes[int(probeId)].gprobebits = int(flags, 16)

        for key in cmdVar.getKeyVarData(gprobeKey):
            try:
                gprobes[key[1]].from_platedb_gprobe(key)
            except (KeyError, ValueError), e:
                cmd.warn('text=%s' % e)
                cmd.warn('text="Unknown probeId %s from platedb.gprobe. %s"' % (probeId, str(key)))
                continue

        # Add in the plate/fibre geometry from plPlugMapM
        plPlugMapMKey = actorState.models["platedb"].keyVarDict["plPlugMapM"]
        guideInfoKey = actorState.models["platedb"].keyVarDict["guideInfo"]

        cmdVar = actorState.actor.cmdr.call(actor="platedb", forUserCmd=cmd,
                                            cmdStr="getGprobesPlateGeom %s" % (extraArgs),
                                            keyVars=[guideInfoKey, plPlugMapMKey])
        if cmdVar.didFail:
            cmd.fail("text=%s" % qstr("Failed to lookup gprobes's geometry for cartridge %d" % (cartridge)))
            return
        assert int(cmdVar.getLastKeyVarData(plPlugMapMKey)[0]) == plate
        fscanMJD = cmdVar.getLastKeyVarData(plPlugMapMKey)[1]
        fscanID = cmdVar.getLastKeyVarData(plPlugMapMKey)[2]

        # unpack the platedb guideInfo keys into the probe
        for key in cmdVar.getKeyVarData(guideInfoKey):
            try:
                gprobes[key[0]].from_platedb_guideInfo(key)
            except (KeyError, ValueError), e:
                cmd.warn('text=%s' % e)
                cmd.warn('text="Unknown probeId %d from plugmap file. %s"' % (key[0], str(key)))
                continue

        # Add in the refraction functions from plateGeomCoeffs
        #
        # I'm not sure how to get numeric pointing IDs, but it turns out that
        # shared plates will only ever have one pointing.
        pointingID = 1
        if pointing != 'A':
            cmd.warn('text="pointing name is %s, but we are using pointing #1. This is probably OK."' % (pointing))

        self.addGuideOffsets(cmd, plate, pointingID, gprobes)

        # LCOHACK: test adding profile (focus) errors based on guide x/yFocal
        self.add_prof_offsets(cmd, plate, gprobes)

        # Send that information off to the master thread
        #
        queue.put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                  cartridge=cartridge, plate=plate, pointing=pointing,
                  fscanMJD=fscanMJD, fscanID=fscanID,
                  boresight_ra=boresight_ra, boresight_dec=boresight_dec,
                  design_ha=design_ha, survey=survey, surveyMode=surveyMode,
                  gprobes=gprobes))

    def addGuideOffsets(self, cmd, plate, pointingID, gprobes):
        """
        Read in the new (needed for APOGEE/MARVELS) plateGuideOffsets interpolation arrays.
        """

        # Get .par file name in the platelist product.
        # plates/0046XX/004671/plateGuideOffsets-004671-p1-l16600.par
        for wavelength in (16600,):
            try:
                path = os.path.join(os.environ['PLATELIST_DIR'],
                                    'plates',
                                    '%04dXX' % (int(plate/100)),
                                    '%06d' % (plate),
                                    'plateGuideOffsets-%06d-p%d-l%05d.par' % (plate, pointingID, wavelength))
            except:
                cmd.warn('text="no refraction corrections for plate %d at %dA, could not locate PLATELIST_DIR"' % (plate, wavelength))
                continue
            if not os.path.exists(path):
                cmd.warn('text="no refraction corrections for plate %d at %dA"' % (plate, wavelength))
                continue

            try:
                ygo = YPF.YPF(path)
                guideOffsets = ygo.structs['HAOFFSETS'].asObjlist()
                cmd.inform('text="loaded guider coeffs for %dA from %s"' % (wavelength, path))
            except Exception, e:
                cmd.warn('text="failed to read plateGuideOffsets file %s: %s"' % (path, e))
                continue

            for gpID, gProbe in gprobes.items():
                if gProbe.fiber_type == 'TRITIUM':
                    continue

                offset = [o for o in guideOffsets if o.holetype == "GUIDE" and o.iguide == gpID]
                if len(offset) != 1:
                    cmd.warn('text="no or too many (%d) guideOffsets for probe %s"' % (len(offset), gpID))
                    continue

                gProbe.haOffsetTimes[wavelength] = offset[0].delha
                gProbe.haXOffsets[wavelength] = offset[0].xfoff
                gProbe.haYOffsets[wavelength] = offset[0].yfoff


    def add_prof_offsets(self, cmd, plate, gprobes):
        dpf = DuPontProfile()
        dpf.getProfileFromDB(int(plate))
        gprobe_ids = sorted(gprobes)
        for gprobe_id in gprobe_ids:
            if not gprobes[gprobe_id].exists:
                #probably gprobe 17 or tritium?
                continue
            xFocal = gprobes[gprobe_id].xFocal
            yFocal = gprobes[gprobe_id].yFocal
            profOffset = dpf.getErr(xFocal,yFocal)
            gprobes[gprobe_id].focusOffset += profOffset*1000.0 # mm to microns

    def add_cmm_offsets(self, cmd, plate, fscan_id, pointing, gprobes):
        """Gets the CMM offsets from a file in etc and adds them as gprobe x/yFocal."""

        cmm_file = os.path.realpath(
            os.path.join(os.path.dirname(__file__),
                         '../../../etc/plate{0:04d}_errs.txt'.format(plate)))

        if not os.path.exists(cmm_file):
            cmd.warn('text="no CMM measurements found."')
            return

        cmm_data = np.loadtxt(cmm_file)

        n_guides = 16
        pointing_num = ['a', 'b', 'c', 'd'].index(pointing.lower())
        pre_index = n_guides * 3 * pointing_num
        fiberid_range = pre_index + np.arange(1 + (fscan_id - 1) * n_guides,
                                              1 + fscan_id * n_guides)

        cmd.inform('text="loading CMM offsets for range {0}."'.format(fiberid_range))

        cmm_errors = cmm_data[:, [2, 3]][fiberid_range]

        gprobe_ids = sorted(gprobes)
        for gprobe_id in gprobe_ids:
            if gprobe_id > len(cmm_errors):
                cmd.warn('text="gprobe={0}: no CMM measurements found."'.format(gprobe_id))
                continue
            cmm_x_offset, cmm_y_offset = cmm_errors[gprobe_id - 1]
            cmd.inform('text="gprobe={0}: CMM errors x={1:+.3e}, y={2:+.3e}"'
                       .format(gprobe_id, cmm_x_offset, cmm_y_offset))
            gprobes[gprobe_id].xFocal += cmm_x_offset
            gprobes[gprobe_id].yFocal += cmm_y_offset

        return

    def setRefractionBalance(self, cmd):
        """Set refraction balance to a specific correction ratio, or based on plateType/surveyMode."""
        keywords = cmd.cmd.keywords
        corrRatio = keywords["corrRatio"].values[0] if 'corrRatio' in keywords else None
        plateType = keywords["plateType"].values[0] if 'plateType' in keywords else None
        surveyMode = keywords["surveyMode"].values[0] if 'surveyMode' in keywords else None

        myGlobals.actorState.queues[guiderActor.MASTER].put(
            Msg(Msg.SET_REFRACTION, corrRatio=corrRatio, plateType=plateType, surveyMode=surveyMode, cmd=cmd))

    def ping(self, cmd):
        """ Top-level 'ping' command handler. Query the actor for liveness/happiness. """

        self.actor.sendVersionKey(cmd)
        cmd.finish('text="pong"')

    def restart(self, cmd):
        """Restart the worker threads"""

        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.STOP_GUIDING, cmd=cmd))

        actorState = myGlobals.actorState

        cmd.inform('text="Restarting threads or at least _trying to_)"')

        # We can't finish this command after all the threads have died, 'cuz we might not get there.
        actorState.actor.startThreads(actorState, actorState.actor.bcast, restart=True, restartQueues=True)

        cmd.finish("text=\"Nunc dimittis servum tuum Domine\"")

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

    def decenter(self, cmd):
        """Enable/disable decentered guiding."""
        on = "on" in cmd.cmd.keywords
        masterQueue = myGlobals.actorState.queues[guiderActor.MASTER]
        masterQueue.put(Msg(Msg.DECENTER, cmd=cmd, enable=on))

    def setDecenter(self, cmd):
        """Specify absolute offset location for decentered guiding."""
        keywords = cmd.cmd.keywords
        # for now Decenter rot is around (RA+decenterRA, Dec+decenterDec)
        decenters = {}
        decenters['decenterRA'] = keywords["decenterRA"].values[0] if "decenterRA" in keywords else 0
        decenters['decenterDec'] = keywords["decenterDec"].values[0] if "decenterDEC" in keywords else 0

        # Though these are currently available, we don't want to use them.
        if "decenterRot" in keywords:
            cmd.fail('Guider cannot apply a decenter in Rotation (yet).')
            return

        masterQueue = myGlobals.actorState.queues[guiderActor.MASTER]
        masterQueue.put(Msg(Msg.DECENTER, cmd=cmd, decenters=decenters))

    def mangaDither(self, cmd):
        """Specify a particular manga dither position for decentered guiding."""
        # ra, dec, rot
        dithers = {'N': {'decenterRA': -0.417, 'decenterDec': +0.721, 'decenterRot': 0.0},
                   'S': {'decenterRA': -0.417, 'decenterDec': -0.721, 'decenterRot': 0.0},
                   'E': {'decenterRA': +0.833, 'decenterDec': 0., 'decenterRot': 0.0},
                   'C': {'decenterRA': 0., 'decenterDec': 0., 'decenterRot': 0.0}}
        ditherPos = cmd.cmd.keywords['ditherPos'].values[0]
        try:
            decenters = dithers[ditherPos]
            decenters['mangaDither'] = ditherPos
        except KeyError:
            cmd.fail("text=%s" % qstr("Failed to parse manga dither position: %s" % ditherPos))
        else:
            masterQueue = myGlobals.actorState.queues[guiderActor.MASTER]
            masterQueue.put(Msg(Msg.DECENTER, cmd=cmd, decenters=decenters))

    def makeMovie(self, cmd):
        """Create a movie of guider images in /data/gcam/movieMJD from a range of exposures from start to end."""
        mjd = cmd.cmd.keywords['movieMJD'].values[0] if 'movieMJD' in cmd.cmd.keywords else None
        start = cmd.cmd.keywords['start'].values[0]
        end = cmd.cmd.keywords['end'].values[0]
        movieQueue = myGlobals.actorState.queues[guiderActor.MOVIE]
        movieQueue.put(Msg(Msg.MAKE_MOVIE, cmd=cmd,
                           mjd=mjd, start=start, end=end, finish=True))

    def ecam_findstar(self, cmd):
        """
        Take one ecam exposure, reduce it, and output the stars found therein.
        """
        time = cmd.cmd.keywords['time'].values[0] if 'time' in cmd.cmd.keywords else 5

        # TODO: Can't change ecam binning yet!
        bin = cmd.cmd.keywords['bin'].values[0] if 'bin' in cmd.cmd.keywords else 1

        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        queue.put(Msg(Msg.START_GUIDING, cmd=cmd, expTime=time, oneExposure=True,
                  bin=bin, camera='ecamera'))

    def sendfield(self, cmd):
        """Commands the TCC to go to the RA/Dec of the loaded cart."""

        actorState = guiderActor.myGlobals.actorState
        gState = actorState.gState
        ra, dec = gState.boresight_ra, gState.boresight_dec

        cmdVar = actorState.actor.cmdr.call(
            actor='tcc', forUserCmd=cmd,
            cmdStr='target {0:.5f}, {1:.5f} icrs'.format(ra, dec))

        if cmdVar.didFail:
            cmd.warn('text="Failed to send field."')

        return

    def loadFakeGProbe(self, cmd):
        """Loads a fake cart (99) with one synthetic gprobe."""

        queue = myGlobals.actorState.queues[guiderActor.MASTER]

        gprobeX = (cmd.cmd.keywords['gprobeX'].values[0]
                   if 'gprobeX' in cmd.cmd.keywords else None)
        gprobeY = (cmd.cmd.keywords['gprobeY'].values[0]
                   if 'gprobeY' in cmd.cmd.keywords else None)
        gprobeRadius = (cmd.cmd.keywords['gprobeRadius'].values[0]
                        if 'gprobeRadius' in cmd.cmd.keywords else None)

        gprobeRA = (cmd.cmd.keywords['gprobeRA'].values[0]
                    if 'gprobeRA' in cmd.cmd.keywords else None)
        gprobeDec = (cmd.cmd.keywords['gprobeDec'].values[0]
                     if 'gprobeDec' in cmd.cmd.keywords else None)

        if gprobeX is None or gprobeY is None or gprobeRadius is None:
            cmd.fail('text=%s' % qstr('Missing gprobe information.'))
            return

        # Fakes cart, plate, and plugging information.
        cartridge = self.fakeCartNumber
        plate = -99
        pointing = 'A'
        fscanMJD = -999
        fscanID = 1
        design_ha = 0.
        survey = 'NONE'
        surveyMode = 'NONE'

        GOOD_PROBE = '0x00'

        gprobes = {}
        gprobes[1] = GuiderState.GProbe(int(1))
        gprobes[1].gprobebits = int(GOOD_PROBE, 16)

        # Fakes the gprobeKey information
        gprobeKey = [
            self.fakeCartNumber,
            1,                    # gprobe_id
            True,                 # exists
            gprobeX,
            gprobeY,
            gprobeRadius,
            0,                    # Rotation
            0,                    # x_ferrule_offset
            0,                    # y_ferrule_offset
            0,                    # focus_offset
            'ACQUIRE'             # fiber_type
            ]

        gprobes[1].from_platedb_gprobe(gprobeKey)

        # Fakes gprobeInfoKey
        gprobeInfoKey = [
            1,
            gprobeRA,
            gprobeDec,
            0,            # xFocal
            0,            # yFocal
            90,           # phi
            0             # throughput
            ]

        gprobes[1].from_platedb_guideInfo(gprobeInfoKey)

        queue.put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                  cartridge=cartridge, plate=plate, pointing=pointing,
                  fscanMJD=fscanMJD, fscanID=fscanID,
                  boresight_ra=gprobeRA, boresight_dec=gprobeDec,
                  design_ha=design_ha, survey=survey, surveyMode=surveyMode,
                  gprobes=gprobes))

    def oneStep(self, cmd):
        """Centers stars in acquisition gprobes."""

        expTime = (cmd.cmd.keywords['time'].values[0]
                   if 'time' in cmd.cmd.keywords else 0.5)
        myGlobals.actorState.queues[guiderActor.MASTER].put(
            Msg(Msg.ONESTEP, cmd=cmd, expTime=expTime))
