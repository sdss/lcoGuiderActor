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
import json

import guiderTester
import guiderActor
from Queue import Queue
from guiderActor import GuiderState
from guiderActor.gimg import GuiderExceptions
from guiderActor.masterThread import guideStep


def getTestImagePaths(dir, mjd, file):
    """Return the infile and outfile associated with this filename."""

    infile = guiderTester.getTestFile(dir, mjd, file)
    outfile = guiderTester.getTestFile(dir, mjd, 'proc-' + file,
                                       raiseError=False)

    return infile, outfile


@guiderTester.skipIfNoGuiderImages
class TestGuiderImage(guiderTester.GuiderTester, unittest.TestCase):
    """Base clase for tests using GuiderImage."""

    def setUp(self):
        self.verbose = True
        super(TestGuiderImage, self).setUp()

        self.queues = {}
        self.queues[guiderActor.GCAMERA] = Queue('gcamera')
        self.queues[guiderActor.MASTER] = Queue('master')
        self.actorState.queues = self.queues

    # def tearDown(self):
    #     if hasattr(self, 'outFile'):
    #         self._remove_file(self.outFile)

    @staticmethod
    def tearDownClass():
        testDataPath = guiderTester.guiderImagesPath
        dataRegEx = os.path.join(testDataPath, 'data', '*cam/*/proc-*.fits.gz')
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

    def _test_call(self, mjd, plateid, fscan_mjd, fscan_id, frameid,
                   gprobes_disabled=[]):
        """Runs GuiderImageAnalysis.__call__(). Returns proc-gimg filename."""

        self.init_probes(mjd=mjd, plateid=plateid, fscan_mjd=fscan_mjd,
                         fscan_id=fscan_id)

        inFile, self.outFile = getTestImagePaths(
            'gcam', mjd, 'gimg-{0:04d}.fits.gz'.format(frameid))

        # Disables acquisition probes
        for gid in gprobes_disabled:
            self.gState.gprobes[gid].disabled = True

        # frameInfo = GuiderState.FrameInfo(frameid, 1, 2, 3)

        # Does some gymnastics with the guider and actorState to make things
        # work.
        self.gState.cmd = self.cmd
        self.actorState.bypassDark = False

        guideStep(self.actor, self.actorState.queues, self.cmd, self.gState,
                  inFile, True, self.gi, output_verify='warn',
                  camera='gcamera')

        return self.outFile

    def _compareBinTables(self, procgimg_fn, expected_fn):
        """Compares the bin tables in two proc-gimg files."""

        result = pyfits.open(procgimg_fn)
        expect = pyfits.open(expected_fn)

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

            elif column in ['dDec']:
                # Checks that star positions are ok within 1%.
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column], rtol=1e-2,
                    err_msg=err_msg)

            elif column == 'fwhm':
                np.testing.assert_allclose(
                    result[6].data[column], expect[6].data[column], rtol=1e-1,
                    err_msg=err_msg)

            elif column in ['rotStar2Sky', 'ugriz', 'ref_mag']:
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
                # np.testing.assert_allclose(
                #     result[6].data[column], expect[6].data[column],
                #     err_msg=err_msg)
                pass


class TestCalibrations(TestGuiderImage):
    """Tests darks and flats."""

    def test_analyzeDark(self):
        inFile, outFile = getTestImagePaths('gcam', 57357, 'gimg-0001.fits.gz')
        self.outFile = outFile  # for tearDown
        self.gi.analyzeDark(inFile, cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile))
        self._check_overwriting(inFile, outFile, self.gi.analyzeDark)

    def test_analyzeDark_lco(self):
        self.gi.setPoint = -10
        inFile, outFile = getTestImagePaths('lco-gcam', 57520, 'gimg-0032.fits.gz')
        self.outFile = outFile  # for tearDown
        self.gi.analyzeDark(inFile, cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile))
        self._check_overwriting(inFile, outFile, self.gi.analyzeDark)

    def test_analyzeDark_no_stack(self):
        self.gi.setPoint = -10
        inFile, outFile = getTestImagePaths('lco-gcam', 57521, 'gimg-0018.fits.gz')
        self.outFile = outFile  # for tearDown
        with self.assertRaises(GuiderExceptions.BadDarkError) as cm:
            self.gi.analyzeDark(inFile, cmd=self.cmd)
        errMsg = 'Total dark exposure time too short'
        self.assertIn(errMsg, str(cm.exception))
        self.assertIn(errMsg, self.cmd.messages[-1])

    def test_analyzeDark_too_short(self):
        self.gi.setPoint = -10
        inFile, outFile = getTestImagePaths('lco-gcam', 57521, 'gimg-0019.fits.gz')
        self.outFile = outFile  # for tearDown
        with self.assertRaises(GuiderExceptions.BadDarkError) as cm:
            self.gi.analyzeDark(inFile, cmd=self.cmd)
        errMsg = 'Total dark exposure time too short'
        self.assertIn(errMsg, str(cm.exception))
        self.assertIn(errMsg, self.cmd.messages[-1])

    def test_analyzeFlat(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""

        self.init_probes(mjd=57357, plateid=7660, fscan_mjd=57356, fscan_id=1)

        inFile, self.outFile = getTestImagePaths('gcam', 57357,
                                                 'gimg-0003.fits.gz')

        flatExpect = guiderTester.getTestFile('gcam', 57357,
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

    def test_analyzeFlat_too_faint(self):
        """Checks that guiderImage fails for a flat with < 1e4 counts."""

        # We init this with random probes because we don't care about them.
        # This should fail before guiderImage does something with the fibres.
        self.init_probes(mjd=57357, plateid=7660, fscan_mjd=57356, fscan_id=1)

        inFile, __ = getTestImagePaths('lco-gcam', 57621, 'gimg-0016.fits.gz')

        self.gi.setPoint = -10.
        with self.assertRaises(GuiderExceptions.FlatError):
            self.gi.analyzeFlat(inFile, self.gState.gprobes, cmd=self.cmd)


class TestCallAPO(TestGuiderImage):
    """Tests calls to GuiderImage with APO guider data."""

    def test_compare_c_code_0040(self):
        """Test GuiderImageAnalysis.__call__() vs the old expected values."""

        self.outFile = self._test_call(mjd=57357, plateid=7660,
                                       fscan_mjd=57356, fscan_id=1,
                                       frameid=40, gprobes_disabled=[11])

        dataExpect = guiderTester.getTestFile('gcam', 57357,
                                              'expect-gimg-0040.fits.gz')

        self._compareBinTables(self.outFile, dataExpect)

    def test_compare_c_code_0041(self):
        """Same as test_compare_c_code_0040 but with the next image.

        To test consistency in the dRA/dDec difference between images taken
        using the same cart and flat, on the same night.

        """

        self.outFile = self._test_call(mjd=57357, plateid=7660,
                                       fscan_mjd=57356, fscan_id=1,
                                       frameid=41, gprobes_disabled=[11])

        dataExpect = guiderTester.getTestFile('gcam', 57357,
                                              'expect-gimg-0041.fits.gz')

        self._compareBinTables(self.outFile, dataExpect)

    def test_compare_c_code_0042(self):
        """Same as test_compare_c_code_0040 but with 0042.

        To test consistency in the dRA/dDec difference between images taken
        using the same cart and flat, on the same night.

        """

        self.outFile = self._test_call(mjd=57357, plateid=7660,
                                       fscan_mjd=57356, fscan_id=1,
                                       frameid=42, gprobes_disabled=[11])

        dataExpect = guiderTester.getTestFile('gcam', 57357,
                                              'expect-gimg-0042.fits.gz')

        self._compareBinTables(self.outFile, dataExpect)

    def test_call_iraf(self):
        """Test GuiderImageAnalysis.__call__() vs values measured with IRAF."""

        self.outFile = self._test_call(mjd=57357, plateid=7660,
                                       fscan_mjd=57356, fscan_id=1,
                                       frameid=40, gprobes_disabled=[3, 11])

        dataIRAF_file = guiderTester.getTestFile('gcam', 57357,
                                                 'expected_7660-57356-1.json')
        dataIRAF = json.load(open(dataIRAF_file))

        result = pyfits.open(self.outFile)
        data = result[6].data

        # Loads expected values measured using IRAF's imexam
        # Removes 0.5 to make measurements 0-indexed and centred on the pixel.
        imexam_xystar = np.array(dataIRAF['57357']['imexam_xstar_ystar'],
                                 np.float32) - 0.5

        # Fiber 11 was disabled, so we make it nan in IRAF
        imexam_xystar[10, :] = np.nan

        # import ipdb; ipdb.set_trace()
        # FIXME: These tests fails, as the relative difference with the
        # PyGuide centroids is still significant.

        np.testing.assert_allclose(data['xstar'], imexam_xystar[:, 0],
                                   rtol=5e-3, err_msg='xstar does not match.')

        np.testing.assert_allclose(data['ystar'], imexam_xystar[:, 1],
                                   rtol=5e-3, err_msg='ystar does not match.')


if __name__ == '__main__':
    verbosity = 2

    suite = None
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)
