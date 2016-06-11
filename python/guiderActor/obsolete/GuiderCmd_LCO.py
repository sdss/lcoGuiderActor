#!/usr/bin/env python
"""Define the available commands for APO and how to handle them."""

import os
import math
import itertools

import numpy as np

from sdss.utilities import yanny

import opscore.protocols.keys as keys
import opscore.protocols.types as types
from opscore.utility.qstr import qstr
import opscore.utility.YPF as YPF

from guiderActor.Commands import GuiderCmd

from guiderActor import Msg, GuiderState
import guiderActor
import guiderActor.myGlobals as myGlobals


def hardingRotation(apoMeasuredRot):
    """Convert Dan Long's measured rotation angle
    To the angle that Paul Harding determined
    """
    return float(apoMeasuredRot) + 18.


def getGprobeKeys():
    """Output a list of gprobeKey where gprobeKey itself is also a list
    gprobeKey needs specific ordering for use by GuiderState, here's code
    that uses gprobeKey:

        self._check_id(gprobeKey[1],'platedb.gprobe')
        self.broken = False if gprobeKey[2] else True
        if self.broken:
            self.disabled = True
        self.xCenter = gprobeKey[3]
        self.yCenter = gprobeKey[4]
        self.radius = gprobeKey[5]
        self.rotation = gprobeKey[6]
        self.xFerruleOffset = gprobeKey[7]
        self.yFerruleOffset = gprobeKey[8]
        self.focusOffset = gprobeKey[9]
        self.fiber_type = gprobeKey[10]
    """
    prodDir = os.environ['GUIDERACTOR_DIR']
    gprobes = yanny.yanny(os.path.join(prodDir, 'etc/gcamFiberInfo_LCO.par'), np=True)['GPROBE']
    gprobeFields = [
        "cartridgeId",
        "gProbeId",
        "exists",
        "xcen",
        "ycen",
        "radius",
        "rot",
        "xferruleOffset",
        "yferruleOffset",
        "focusOffset",
        "fiberType",
    ]  # order matters
    # throw out last value (we don't want tritium source),
    # castMap = (int,)*3 + (float,)*3 +(hardingRotation,) +(float,)*3 + (str,)
    # note removing harding rotation for Nov Run
    castMap = (int,)*3 + (float,)*7 + (str,)
    # LCOHACK: warning!!! remove the [:-1] if the tritium source gets removed from the gcamFiberInfo file!!!
    gprobeKeysNumpy = np.asarray([gprobes[key][:-1] for key in gprobeFields]).T
    assert gprobeKeysNumpy.shape == (16, len(gprobeFields))
    # convert to list of mixed types (numpy got messy here)
    gprobeKeys = []
    for gProbeKey in gprobeKeysNumpy:
        gprobeKeys.append([castFunc(val) for val, castFunc in itertools.izip(gProbeKey, castMap)])
    return gprobeKeys


def getGuideInfoKey(gProbeId, guideNumber, plYanny):
    """
    gProbeId: int 1-16 (17 not uesd)
    guideNumber: fiberNumber corresponding to this gProbeId in the plYanny
    plYanny: a yanny-parsed plPlugMap file

    Output a list of guideInfoKey (list of list):
    guideInfoKey looks like this:
        self._check_id(guideInfoKey[0],'platedb.guideInfo')
        self.ra = guideInfoKey[1]
        self.dec = guideInfoKey[2]
        self.xFocal = guideInfoKey[3]
        self.yFocal = guideInfoKey[4]
        self.phi = guideInfoKey[5]
        self.throughput = guideInfoKey[6]

    note throughput value is divided by 65535.0 and phi always = 0
    in PlatedbCmd.py so I will do the same here

    list output = [gprobId, ra, dec, xFocal, yFocal, phi, throughput]
    """
    assert gProbeId in range(1, 17)
    guideInfo = [gProbeId]
    alignXYpos = []
    # find the values corresponding to guideNumber in the yanny file
    # search the file for object holeType=GUIDE and fiberID = guideNumber
    objs = plYanny["PLUGMAPOBJ"]
    indexGuide = np.where((objs["holeType"] == "GUIDE") & (objs["fiberId"] == guideNumber))
    indexAlign = np.where((objs["holeType"] == "ALIGNMENT") & (objs["fiberId"] == guideNumber))
    attrList = ["ra", "dec", "xFocal", "yFocal", "phi", "throughput"]  # order matters
    for attr in attrList:
        if attr == "phi":
            assert len(alignXYpos) == 2  # paranoia
            alignX, alignY = alignXYpos
            guideX = guideInfo[3]
            guideY = guideInfo[4]
            phi = 90. - math.atan2(float(alignY) - float(guideY), float(alignX) - float(guideX))*180/math.pi
            guideInfo.append(phi)
        elif attr == "throughput":
            guideInfo.append(float(objs[attr][indexGuide]/65535.0))
        else:
            guideInfo.append(float(objs[attr][indexGuide]))
        # tag on alignment x/yFocal
        if attr in ["xFocal", "yFocal"]:
            # xfocal comes first always
            alignXYpos.append(float(objs[attr][indexAlign]))
    return guideInfo


class GuiderCmd_LCO(GuiderCmd.GuiderCmd):

    # globals here, only need definition once
    gprobekeys = getGprobeKeys()
    validpointings = ["A", "B", "C", "D"]

    def __init__(self, actor):
        # initialize from the superclass
        super(GuiderCmd_LCO, self).__init__(actor)

        # Define some new command keywords
        self.keys = keys.KeysDictionary("guider_guiderlco", (3, 0),
                                        keys.Key("plate", types.Int(), help="A plugplate ID"),
                                        keys.Key("pointing", types.String(), help="A pointing for the given plugplate"),
                                        keys.Key("fiberPos", types.Int(), help="A fiber position, 1-indexed"),
                                        keys.Key("pmDir", types.String(), help="Directory where the plPlugMapP files are.")
                                        )
        # Define new commands for APO
        self.vocab = [
            ("fakeCartridge", "<plate> <fiberPos> <pointing> [<pmDir>]", self.fakeCartridge),
            ]

    def fakeCartridge(self, cmd):
        """
        Load a cartridge but grab all info from plPlugMap and gcamFiberInfo_LCO files on disk (rather than from the db)
        """
        plate = int(cmd.cmd.keywords["plate"].values[0])
        fiberPos = int(cmd.cmd.keywords["fiberPos"].values[0]) if "fiberPos" in cmd.cmd.keywords else 1
        if fiberPos not in [1, 2, 3]:
            cmd.fail("text=\"fiberPos parameter must be 1, 2 or 3 in fakeCartridge\"")
            return
        pointing = cmd.cmd.keywords["pointing"].values[0] if "pointing" in cmd.cmd.keywords else "A"
        if pointing not in self.validpointings:
            cmd.fail("text=\"pointing parameter must be one of %s in fakeCartridge\"" % ", ".join(self.validpointings))
            return
        plPlugMapDir = cmd.cmd.keywords["pmDir"].values[0] if "pmDir" in cmd.cmd.keywords else "/data/plPlugMap/"

        queue = myGlobals.actorState.queues[guiderActor.MASTER]

        # get the path to the plPlugMap file
        pmFilePath = os.path.join(plPlugMapDir, "plPlugMapP-%i%s.par" % (plate, "" if pointing == "A" else pointing))
        # parse the plPlugMapFile
        plYanny = yanny.yanny(pmFilePath, np=True)
        # may use plPlugMap.getPointingInfo() to return updated keyword
        pointingIndex = self.validpointings.index(pointing)
        boresight_ra = float(plYanny["raCen"])
        boresight_dec = float(plYanny["decCen"])
        # design_ha = float(plYanny["ha"][pointingIndex])
        design_ha = 0.0
        survey = "APOGEE-2"
        surveyMode = "APOGEE lead"
        if design_ha < 0:
            design_ha += 360

        # Lookup the valid gprobes
        # the correct gprobe number corresponds to one of 3 possible holes
        # we have chosen to plug for each fiber
        # 1-16 correspond to pointing 1 fiberPos 1
        # 13-32 pointing 1 fiberPos 2
        # ...
        # ...
        # 49 - 64 pointing 2 fiberPos 1...
        # 65 - 80 pointing 2 fiberPos 2 ... you get the picture
        # determine the starting guide number
        #
        guideStartNum = 48*(pointingIndex) + 16*(fiberPos-1) + 1
        guideNums = np.arange(guideStartNum, guideStartNum+16)  # 16 guide fibers
        gprobes = {}
        for guideNum, gProbeKey in itertools.izip(guideNums, self.gprobekeys):
            gprobeId = int(gProbeKey[1])  # index 0 is cart 1 is gProbeId
            # alignmentPos is [xFocal, yFocal]
            cmd.diag('text="probeNum,guideNum: {},{}"'.format(gprobeId, guideNum))
            guideInfoKey = getGuideInfoKey(gprobeId, guideNum, plYanny)
            gProbe = GuiderState.GProbe(id=gprobeId, gprobeKey=gProbeKey, guideInfoKey=guideInfoKey)
            # need to explicitly set the gprobebits
            if not bool(gProbeKey[2]):
                # set flag to bad (defaults to good)
                # gProbeKey[2] is exists
                gProbe.broken = True
            gprobes[gprobeId] = gProbe

        #
        # I'm not sure how to get numeric pointing IDs, but it turns out that
        # shared plates will only ever have one pointing.
        pointingID = 1

        self.addGuideOffsets(cmd, plate, pointingID, gprobes)

        # Send that information off to the master thread
        #
        queue.put(Msg(Msg.LOAD_CARTRIDGE, cmd=cmd,
                  cartridge=1, plate=plate, pointing=pointing,
                  fscanMJD=-1, fscanID=-1,
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

