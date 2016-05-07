#!/usr/bin/env python
# encoding: utf-8
#
# create_plPlugMapM_LCO.py
#
# Created by José Sánchez-Gallego on 6 May 2016.
# Licensed under a 3-clause BSD license.
#
# Revision history:
#    6 May 2016 J. Sánchez-Gallego
#       Initial version


from __future__ import division
from __future__ import print_function
import sys
import os
from sdss.utilities import yanny
import numpy as np
import string


cart = 20
mjd = 57514

template = """
EVILSCAN
fscanVersion $HeadURL: https://svn.sdss.org/repo/operations/general/idlmapper/v6_0_7/src/evilscan.c $
pluggers     Jose
plateId      {plateID}
fscanMJD     {mjd}
fscanId      {fscanID}
fscanDate    Fri May  6 12:00:00 2016
fscanFile    fiberScan-{plateID}-57514-{fscanID:02d}.par
fscanMode    interpolated
fscanSpeed   400
fscanRows    960
fscanCols    960
fscanBias    45.000000
motorId1     35
motorId2     31
motorId3     3
cartridgeId {cart}
fmapVersion NOCVS:v6_0_7
idlutilsVersion v5_5_17
idlVersion 7.1
"""


def addHeader(header, file, plateID, fscanID):
    """Adds a header to a yanny file after the commented section."""

    fscanData = template.format(plateID=plateID, mjd=mjd, fscanID=fscanID,
                                cart=cart)

    header = [''] + header + fscanData.splitlines() + ['']

    fileLines = open(file, 'r').read().splitlines()

    for ii, line in enumerate(fileLines):
        if not line.startswith('#'):
            break

    for line in header[::-1]:
        fileLines.insert(ii, line)

    unit = open(file, 'w')
    for line in fileLines:
        unit.write(line + '\n')
    unit.close()

    return


def doPointing(plateID, pointing, fscanIDs):
    """Creates the plPlugMapM files for a specific pointing."""

    pointingName = '' if pointing == 'A' else pointing
    filename = os.path.join(
        os.environ['PLATELIST_DIR'], 'plates',
        '{0:06d}'.format(plateID)[:-2] + 'XX', '{0:06d}'.format(plateID),
        'plPlugMapP-{0:d}{1}.par'.format(plateID, pointingName))

    assert os.path.exists(filename)

    yannyFile = yanny.yanny(filename, np=True)
    rawFile = open(filename, 'r').read().splitlines()

    plPlugMapObj = yannyFile['PLUGMAPOBJ']

    # Manually retrieves the header from the raw lines.
    header = []
    for line in rawFile:
        if line.strip().startswith('typedef enum {'):
            break
        header.append(line)

    # Loops over fscanIDs
    for fscanID in fscanIDs:

        # Calculates the range of fiberIds that correspond to this pointing
        # and fscanID.
        pointingNum = string.uppercase.index(pointingName)
        preIndex = 16 * 3 * pointingNum
        fiberID_range = preIndex + np.arange(1 + (fscanID - 1) * 16,
                                             1 + fscanID * 16)

        outFileName = 'plPlugMapM-{0}{3}-{1}-{2:02d}.par'.format(
            plateID, mjd, fscanID, pointingName)

        # Gets a list of the holes that we should keep for this fscanID
        validHoles = plPlugMapObj[(plPlugMapObj['holeType'] == 'LIGHT_TRAP') |
                                  (np.in1d(plPlugMapObj['fiberId'],
                                           fiberID_range))]

        enums = {'holeType': ['HOLETYPE', yannyFile._enum_cache['HOLETYPE']],
                 'objType': ['OBJTYPE', yannyFile._enum_cache['OBJTYPE']]}
        yanny.write_ndarray_to_yanny(outFileName, validHoles, enums=enums,
                                     structname='PLUGMAPOBJ')

        # Replaces the guidenums with the slice used for this file
        for ii, line in enumerate(header):
            if line.startswith('guidenums' + str(pointingNum + 1)):
                guides = str('guidenums' + str(pointingNum + 1)) + ' ' + \
                    ' '.join(map(str, fiberID_range))
                header[ii] = guides
                break

        addHeader(header, outFileName, plateID, fscanID)


def create_plPlugMapM_LCO(plateIDs):
    """Converts the plPlugMapP files for a `plateID` into a plPlugMapM.

    This scripts converts the plPlugMapP files for a `plateID` into a
    series of plPlugMapM files for LCO guider testing.

    Each guider commissioning plate contains four pointings, each with 3
    sets of 16 guiding stars. We distribute those 48 holes in different
    fScanIDs.

    """

    pointings = ['A', 'B', 'C', 'D']
    fscanIDs = [1, 2, 3]

    for plateID in plateIDs:
        for pointing in pointings:
            doPointing(plateID, pointing, fscanIDs)

    return


if __name__ == '__main__':

    if len(sys.argv) == 2:
        plateIDs = [int(sys.argv[1])]
    else:
        plateIDs = range(int(sys.argv[1]), int(sys.argv[2]) + 1)

    create_plPlugMapM_LCO(plateIDs)
