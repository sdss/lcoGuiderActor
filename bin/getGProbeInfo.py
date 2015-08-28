#!/usr/bin/env python
import os
import itertools
import numpy

from guiderActor.Commands.GuiderCmd import getGprobeKeys, getGuideInfoKey, GuiderCmd, yanny, GuiderState

gprobekeys = getGprobeKeys()

def getGprobes(plate, pointing, fiberPos):
    plPlugMapDir = "/data/plPlugMap/"
    pmFilePath = os.path.join(plPlugMapDir, "plPlugMapP-%i%s.par"%(plate, "" if pointing == "A" else pointing))
    plYanny = yanny.yanny(pmFilePath, np=True)
    pointingIndex = GuiderCmd.validpointings.index(pointing)
    guideStartNum = 48*(pointingIndex) + 16*(fiberPos-1) + 1
    guideNums = numpy.arange(guideStartNum, guideStartNum+16) # 16 guide fibers
    gprobes = {}
    for guideNum, gProbeKey in itertools.izip(guideNums, gprobekeys):
        gprobeId = int(gProbeKey[1]) # index 0 is cart 1 is gProbeId
        guideInfoKey = getGuideInfoKey(gprobeId, guideNum, plYanny)
        gProbe = GuiderState.GProbe(id=gprobeId, gprobeKey=gProbeKey, guideInfoKey=guideInfoKey)
        # need to explicitly set the gprobebits
        if not bool(gProbeKey[2]):
            # set flag to bad (defaults to good)
            # gProbeKey[2] is exists
            gProbe.broken = True
        gprobes[gprobeId] = gProbe
    return gprobes, guideNums

if __name__ == "__main__":
    plate = 8645
    pointing = "B"
    fiberPos = 3
    guideNumber = 11
    gprobes, guideNums = getGprobes(plate, pointing, fiberPos)
    print "plate=%i, pointing=%s, fiberPos=%i"%(plate, pointing, fiberPos)
    print "guideNumber: %i, expandedNumber: %i"%(guideNumber, guideNums[guideNumber-1])
    print "xFocal: %.6f  yFocal: %.6f"%(gprobes[guideNumber].xFocal, gprobes[guideNumber].yFocal)