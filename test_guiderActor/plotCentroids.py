#!/usr/bin/env python
# encoding: utf-8
#
# plotCentroids.py
#
# Created by José Sánchez-Gallego on 29 Apr 2016.
# Licensed under a 3-clause BSD license.
#
# Revision history:
#    29 Apr 2016 J. Sánchez-Gallego
#       Initial version


from __future__ import division
from __future__ import print_function
from matplotlib import pyplot as plt
from astropy.io import fits
from guiderTester import getTestFile
from PyGuide import findStars, CCDInfo
import numpy as np
import json

# This script plots the contours for a processed guider image and overlays
# the star centroids as measured using a number of methods.

# Loads a processed gimg generated using guiderActor and PyGuide.
ff = fits.open(getTestFile('gcam', 57357, 'processed/proc-gimg-0040.fits.gz'))
data = np.clip(ff[0].data, 400, 10500)

plt.contour(data, levels=np.arange(500, 10000, 1000),
            origin='lower', colors='k', linewidths=0.1)

# Plots the centres of the centroids as calculated for the previous image.
xstar = ff[-1].data['xstar']
ystar = ff[-1].data['ystar']
guider_PyGuide = plt.scatter(xstar, ystar, marker='x', color='r', s=2, lw=0.2,
                             label='guiderActor with PyGuide')

# Idem but adding +1 to x/ystar
guider_PyGuider_1 = plt.scatter(
    xstar + 1., ystar + 1., marker='x', color='g', s=2, lw=0.2,
    label='guiderActor with PyGuide x/y + 1')

# Loads measurements from IRAF.
dataIRAF_file = getTestFile('gcam', 57357, 'expected_7660-57356-1.json')
dataIRAF = json.load(open(dataIRAF_file))

# Removes 0.5 pixels to make measurements 0-indexed.
xystar_imexam = np.array(dataIRAF['57357']['imexam_xstar_ystar']) - 0.5
iraf = plt.scatter(xystar_imexam[:, 0], xystar_imexam[:, 1], marker='x',
                   color='b', s=2, lw=0.2, label='IRAF')

# Finally, we call PyGuide without any sophisticated process, just passing the
# mask and the CCD parameters.

bias = np.median(fits.getdata(getTestFile('gcam', 57357,
                                          'processed/proc-gimg-0001.fits.gz')))
readNoise = 10.4
ccdGain = 1.4
ccdInfo = CCDInfo(bias, readNoise, ccdGain)

stars = findStars(ff[0].data, ff[1].data, None, ccdInfo)
centroids = np.array([star.xyCtr for star in stars[0]])

PyGuide_raw = plt.scatter(centroids[:, 0], centroids[:, 1], marker='x',
                          color='cyan', s=2, lw=0.2,
                          label='PyGuide unadulterated')

plt.xlabel('x [pixels]')
plt.ylabel('y [pixels]')

plt.legend(frameon=False, loc='upper left', markerscale=5, fontsize=8)

plt.savefig('centroids.pdf')

# Outputs some stistics comparing IRAF and PyGuide

iraf_centroids = xystar_imexam[:-1, :]
iraf_centroids_sorted = iraf_centroids[np.argsort(iraf_centroids[:, 0])]

pyguide_centroids_sorted = centroids[np.argsort(centroids[:, 0])]
diff_PyGuide_IRAF = np.abs(1. - pyguide_centroids_sorted /
                           iraf_centroids_sorted)
print('Relative difference between IRAF and PyGuide:')
print('    Max: {0:.3g}'.format(np.max(diff_PyGuide_IRAF)))
print('    Mean: {0:.3g}'.format(np.mean(diff_PyGuide_IRAF)))

guiderActor_1 = np.array([xstar, ystar]).T
guiderActor_1_sorted = guiderActor_1[np.argsort(guiderActor_1[:, 0])]
diff_PyGuide_guiderActor_1 = (np.delete(pyguide_centroids_sorted, -2, axis=0) -
                              guiderActor_1_sorted[:-2])

print('\nAbsolute difference between PyGuide and guiderActor+1:')
print('    Min: {0:.3g}'.format(np.min(diff_PyGuide_guiderActor_1)))
print('    Max: {0:.3g}'.format(np.max(diff_PyGuide_guiderActor_1)))
print('    Mean: {0:.3g}'.format(np.mean(diff_PyGuide_guiderActor_1)))
