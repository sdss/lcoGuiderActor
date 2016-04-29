#!/usr/bin/env python
# encoding: utf-8
"""

test_guiderImage.py

Created by José Sánchez-Gallego on 4 Apr 2016.
Licensed under a 3-clause BSD license.

Revision history:
    4 Apr 2016 J. Sánchez-Gallego
      Initial version

"""

from __future__ import division
from __future__ import print_function
import os
import glob
import unittest
import pyfits
import numpy as np
import guiderTester
from guiderActor import GuiderState


def getTestImagePaths(dir, mjd, file):
    """Return the infile and outfile associated with this filename."""

    infile = guiderTester.getTestImage(dir, mjd, file)
    outfile = guiderTester.getTestImage(dir, mjd, 'proc-' + file,
                                        raiseError=False)

    return infile, outfile


@guiderTester.skipIfNoGuiderImages
class TestGuiderImage(guiderTester.GuiderTester, unittest.TestCase):
    def setUp(self):
        self.verbose = True
        super(TestGuiderImage, self).setUp()

    def tearDown(self):
        if hasattr(self, 'outFile'):
            self._remove_file(self.outFile)

    @staticmethod
    def tearDownClass():
        testDataPath = guiderTester.guiderImagesPath
        dataRegEx = os.path.join(testDataPath, 'data', '?cam/*/proc-*.fits.gz')
        for filename in glob.glob(dataRegEx):
            os.remove(filename)

    def _check_overwriting(self, inFile, outFile, analyze, args=[]):
        """If the proc- file exists, we should only read it, not rewrite it.

        This fails if we try to write to it.

        """

        os.chmod(outFile, 0444)
        analyze(inFile, *args, cmd=self.cmd)
        os.chmod(outFile, 0644)
        self._remove_file(outFile)

    def test_analyzeDark(self):
        """Test GuiderImageAnalysis.analyzeDark()"""

        inFile, outFile = getTestImagePaths('gcam', 57357, 'gimg-0001.fits.gz')
        self.outFile = outFile  # for tearDown
        self.gi.analyzeDark(inFile, cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile))
        self._check_overwriting(inFile, outFile, self.gi.analyzeDark)

    def test_analyzeFlat(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""

        self.init_probes(mjd=57357, plateid=7660, fscan_mjd=57356, fscan_id=1)

        inFile, self.outFile = getTestImagePaths('gcam', 57357,
                                                 'gimg-0003.fits.gz')

        flatExpect = guiderTester.getTestImage('gcam', 57357,
                                               'expect-gimg-0003.fits.gz')

        self.gi.analyzeFlat(inFile, self.gState.gprobes, cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outFile))

        # First, we compare with the proc-gimg file generated using the C
        # code. We will tell us if our figure are in the ballpark.
        # TODO: the current expected results file is missing fiber 11 because
        # it's damaged. If we fix the code to deal with that per #2357 we'll
        # need a test that knows about that fiber!
        result = pyfits.open(self.outFile)
        expect = pyfits.open(flatExpect)

        for column in result[6].data.names:
            if column == 'fiber_type':
                # Trouble testing numpy string arrays
                self.assertTrue(
                    (result[6].data[column] == expect[6].data[column]).all(),
                    '{0} has a mismatch'.format(column))

            elif column not in ['rotStar2Sky', 'ugriz', 'ref_mag']:
                # TODO: remove the if rotStar2Sky once we have a better test!
                # The problem is rotStar2Sky was never saved for old flats.
                # TODO: don't have ugriz data for these to test,
                # but could add some. Would have to pull it from the database.
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column],
                    err_msg='{0} has a mismatch'.format(column))

        self._check_overwriting(inFile, self.outFile,
                                self.gi.analyzeFlat, [self.gState.gprobes])

    def test_call(self):
        """Test GuiderImageAnalysis.__call__()"""

        self.init_probes(mjd=57357, plateid=7660, fscan_mjd=57356, fscan_id=1)

        inFile, self.outFile = getTestImagePaths('gcam', 57357,
                                                 'gimg-0040.fits.gz')
        dataExpect = guiderTester.getTestImage('gcam', 57357,
                                               'expect-gimg-0040.fits.gz')

        # disable the acquisition probes, since the observers did so.
        self.gState.gprobes[3].disabled = True
        self.gState.gprobes[11].disabled = True

        frameInfo = GuiderState.FrameInfo(40, 1, 2, 3)
        self._call_gi(inFile)
        self.gi.writeFITS(self.actorState.models, self.cmd, frameInfo,
                          self.gState.gprobes)

        # First we check with the values generated with the previous
        # implementation of the code.
        result = pyfits.open(self.outFile)
        expect = pyfits.open(dataExpect)

        for column in result[6].data.names:
            err_msg = '{0} has a mismatch'.format(column)
            if column == 'fiber_type':
                # Trouble testing numpy string arrays
                self.assertTrue(
                    (result[6].data[column] == expect[6].data[column]).all(),
                    err_msg)
            elif column in ['xstar', 'ystar']:
                # Checks that star positions are ok within 1%.
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column], rtol=1e-2,
                    err_msg=err_msg)
            elif column == 'fwhm':
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column], rtol=1e-1,
                    err_msg=err_msg)
            elif column in ['rotStar2Sky', 'ugriz', 'ref_mag', 'dx', 'dy',
                            'dRA', 'dDec']:
                # TODO: the problem is rotStar2Sky was not saved for old flats.
                # TODO: don't have ugriz data for these to test,
                # but could add some. Would have to pull it from the database.
                # TODO: dx/dy are currently calculated in masterThreaad, but
                # they could be calculated as part of guiderImage(),
                # if we passed in the frameInfo data.
                # NOTE: is the last point still true?

                continue

            else:
                # FIXME: This test fails for flux, flux, mag, sky, and poserr.
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column],
                    err_msg=err_msg)


if __name__ == '__main__':
    verbosity = 2

    suite = None
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)
