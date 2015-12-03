#!/usr/bin/env python
"""
Test the behavior of guider flats, including finding fibers..
"""
import os
import glob
import unittest
import pyfits
import numpy as np

from actorcore import TestHelper

import guiderTester

from guiderActor import GuiderState
from guiderActor.gimg import GuiderExceptions

def _inout_names(dir,file):
    """Return the infile and outfile associated with this filename."""
    infile = os.path.join(dir,file)
    outfile = os.path.join(dir,'proc-'+file)
    return infile,outfile

class TestGuiderImage(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        self.verbose = True
        # calibration files
        self.inDarkFile,self.outDarkFile = _inout_names('gcam','gimg-0001.fits.gz')
        self.inFlatFile,self.outFlatFile = _inout_names('gcam','gimg-0003.fits.gz')
        self.inFlatEcamFile,self.outFlatEcamFile = _inout_names('ecam','gimg-0003.fits.gz')

        # bad files
        self.saturatedFile = 'gcam/gimg-0004.fits.gz'
        self.badReadFile1 = 'gcam/gimg-0859.fits.gz'
        self.badReadFile2 = 'gcam/gimg-0519.fits.gz'

        # good data files
        self.inDataFile,self.outDataFile = _inout_names('gcam','gimg-0040.fits.gz')

        # values for comparison after processing good data files.
        self.flatExpect = 'gcam/expect-gimg-0003.fits.gz'
        self.dataExpect = 'gcam/expect-gimg-0040.fits.gz'

        super(TestGuiderImage,self).setUp()
        self._init_probes_gimg_0040()

    def tearDown(self):
        if hasattr(self,'outFile'):
            self._remove_file(self.outFile)

    @staticmethod
    def tearDownClass():
        for filename in glob.glob('?cam/proc-*.fits.gz'):
            os.remove(filename)

    def _check_overwriting(self, inFile, outFile, analyze, args=[]):
        """
        If the proc- file exists, we should only read it, not rewrite it.
        This fails if we try to write to it.
        """
        os.chmod(outFile,0444)
        analyze(inFile,*args,cmd=self.cmd)
        os.chmod(outFile,0644)
        self._remove_file(outFile)
        
    def test_analyzeDark(self):
        """Test GuiderImageAnalysis.analyzeDark()"""
        inFile = self.inDarkFile
        outFile = self.outDarkFile
        self.outFile = outFile # for tearDown
        self.gi.analyzeDark(inFile, cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile),'analyzeDark file write: %s'%outFile)
        self._check_overwriting(inFile,outFile,self.gi.analyzeDark)
    
    def test_analyzeFlat(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""
        inFile = self.inFlatFile
        self.outFile = self.outFlatFile
        self.gi.analyzeFlat(inFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outFile),'analyzeFlat file write: %s'%self.outFile)

        # TBD: I know this is a dumb test: we want to do a better analysis
        # (maybe by hand in IRAF or something?) of the images and use that for
        # our comparison here.
        # NOTE TBD: In particular, the current expected results file is missing
        # fiber 11 because it's damaged. If we fix the code to deal with that
        # per #2357 we'll need a test that knows about that fiber!
        result = pyfits.open(self.outFile)
        expect = pyfits.open(self.flatExpect)
        for x in result[6].data.names:
            if x == 'fiber_type':
                # Trouble testing numpy string arrays
                self.assertTrue((result[6].data[x] == expect[6].data[x]).all(),'%s has a mismatch'%x)
            elif x != 'rotStar2Sky' and x != 'ugriz' and x != 'ref_mag':
                # !!!!!!!!!!!!!!!
                # TBD: remove the if rotStar2Sky once we have a better test!
                # The problem is rotStar2Sky was never saved for old flats.
                # !!!!!!!!!!!!!!!
                # TBD: don't have ugriz data for these to test, but could add some.
                # Would have to pull it from the database.
                # !!!!!!!!!!!!!!!
                np.testing.assert_allclose(result[6].data[x],expect[6].data[x],err_msg='%s has a mismatch'%x)
        self._check_overwriting(inFile,self.outFile,self.gi.analyzeFlat,[self.gState.gprobes])

    def test_analyzeFlat_ecam(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""
        inFile = self.inFlatEcamFile
        outFile = self.outFlatEcamFile
        self.outFile = outFile # for tearDown
        self.gi.camera = 'ecamera'
        self.gi.analyzeFlat(inFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile),'analyzeFlat_ecam file write: %s'%outFile)
        # TBD: best test is probably to check that the flat is median ~1.
        self.assertFalse('Need to create a test for this!')
        self._check_overwriting(inFile,outFile,self.gi.analyzeFlat,[self.gState.gprobes])

    def test_call(self):
        """Test GuiderImageAnalysis.__call__()"""

        # disable the acquisition probes, since the observers did so.
        self.gState.gprobes[3].disabled = True
        self.gState.gprobes[11].disabled = True

        self.outFile = self.outDataFile # for tearDown
        frameInfo = GuiderState.FrameInfo(40,1,2,3)
        fibers = self._call_gi(self.inDataFile)
        self.gi.writeFITS(self.actorState.models, self.cmd, frameInfo, self.gState.gprobes)

        # TBD: I know this is a dumb test: we want to do a better analysis
        # (maybe by hand in IRAF or something?) of the images and use that for
        # our comparison here.
        result = pyfits.open(self.outFile)
        expect = pyfits.open(self.dataExpect)

        for x in result[6].data.names:
            if x == 'fiber_type':
                # Trouble testing numpy string arrays
                self.assertTrue((result[6].data[x] == expect[6].data[x]).all(),'%s has a mismatch'%x)
            elif x == 'xstar' or x == 'ystar':
                # !!!!!!!!!!!!!!!
                # TBD: remove once I know I have optimal fit values to compare to.
                # TBD: while I sort out the exact position differences...
                np.testing.assert_allclose(result[6].data[x],expect[6].data[x],rtol=1e-2,err_msg='%s has a mismatch'%x)
            elif x == 'fwhm':
                # !!!!!!!!!!!!!!!
                # TBD: remove once I know I have optimal fit values to compare to.
                # TBD: while I sort out the exact fwhm differences
                np.testing.assert_allclose(result[6].data[x],expect[6].data[x],rtol=1e-1,err_msg='%s has a mismatch'%x)
            elif x != 'rotStar2Sky' and x != 'ugriz' and x != 'ref_mag' and \
                 x != 'dx' and x != 'dy' and x != 'dRA' and x != 'dDec':
                # !!!!!!!!!!!!!!!
                # TBD: remove the if rotStar2Sky once we have a better test!
                # The problem is rotStar2Sky was never saved for old flats.
                # !!!!!!!!!!!!!!!
                # TBD: don't have ugriz data for these to test, but could add some.
                # Would have to pull it from the database.
                # !!!!!!!!!!!!!!!
                # TBD: dx/dy are currently calculated in masterThreaad, but
                # they could be calculated as part of guiderImage(), if we passed
                # in the frameInfo data.
                np.testing.assert_allclose(result[6].data[x],expect[6].data[x],err_msg='%s has a mismatch'%x)

        self.assertFalse(True,'make a test!')

    def test_findStars_ecam(self):
        self.assertFalse(True,"make some tests for this!")
    
    def _temp_run(self,filename,errorText):
        """Help with temperature tests."""
        header = pyfits.getheader(filename)
        self.gi.cmd = self.cmd
        result = self.gi._check_ccd_temp(header)
        self.assertFalse(result)
        self.assertEqual(self.cmd.levels,'w')
        self.assertEqual('text="%s"'%errorText, self.cmd.messages[-1])

    def test_badSetPoint_dark(self):
        """Test what happens when the ccdtemp is outside the setPoint spec for a dark."""
        self.gi.setPoint = self.setPoint_bad
        self._temp_run(self.inDarkFile,'CCD temp signifcantly different (>3.0) from setPoint: -40.1, expected -35.0')

    def test_badSetPoint_image(self):
        """Test what happens when the ccdtemp is outside the setPoint spec for an image."""
        self.gi.setPoint = self.setPoint_bad
        self.gi.darkTemperature = self.setPoint_good
        self._temp_run(self.inDataFile,'CCD temp signifcantly different (>3.0) from setPoint: -40.1, expected -35.0')

    def test_badSetPoint_dark_temp(self):
        """Test what happens when the ccdtemp is outside the dark temp spec for an image."""
        self.gi.darkTemperature = self.setPoint_bad
        self._temp_run(self.inDataFile,'CCD temp signifcantly different (>3.0) from dark temp: -40.1, expected -35.0')

    def test_no_setPoint(self):
        """Test what happens when the setPoint is None."""
        #self.gi.darkTemperature = self.setPoint_bad
        self.gi.setPoint = None
        header = pyfits.getheader(self.inDataFile)
        self.gi.cmd = self.cmd
        with self.assertRaises(GuiderExceptions.GuiderError):
            self.gi._check_ccd_temp(header)

    def test_saturatedImage(self):
        """Test GuiderImageAnalysis.__call__() on a completely saturated image."""
        self.assertRaises(GuiderExceptions.GuiderError,self._call_gi,self.saturatedFile)
        self.assertTrue(self.cmd.levels[-2:] == 'we')
        self.assertTrue('Fully saturated' in self.cmd.messages[-1])
        
    def test_badRead(self):
        """Confirm that bad gcamera reads are rejected immediately."""
        self.gi.cmd = self.cmd
        self.assertRaises(GuiderExceptions.BadReadError,self.gi._pre_process,self.badReadFile1,binning=2)
        self.assertRaises(GuiderExceptions.BadReadError,self.gi._pre_process,self.badReadFile2,binning=2)
    
    def test_dither_headers(self):
        """check that the dither keywords get into the header."""
        self.actorState.models['guider'] = TestHelper.Model('guider',TestHelper.guiderState['guiderOnDecenter'])
        self.gi.cmd = self.cmd
        objectname = 'somename'
        hdu = pyfits.open(self.inDataFile)[0]
        frameInfo = GuiderState.FrameInfo(-1,1,2,3)
        self.gi.fillPrimaryHDU(self.cmd,self.actorState.models,hdu,frameInfo,objectname)
        self.assertEqual(hdu.header['MGDPOS'],'N')
#...

if __name__ == '__main__':
    verbosity = 2

    suite = None
    # to test just one piece
    # suite = unittest.TestLoader().loadTestsFromName('test_guiderImage.TestGuiderImage.test_call')
    # suite = unittest.TestLoader().loadTestsFromName('test_guiderImage.TestGuiderImage.test_analyzeFlat')
    # suite = unittest.TestLoader().loadTestsFromName('test_guiderImage.TestGuiderImage.test_analyzeFlat_ecam')
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestGuiderImage)
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)

