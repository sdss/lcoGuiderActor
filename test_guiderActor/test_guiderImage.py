#!/usr/bin/env python
"""
Test the behavior of guider flats, including finding fibers..
"""
import os
import unittest
import pyfits

from actorcore import TestHelper

import guiderTester

from guiderActor import GuiderState
from guiderActor.gimg import GuiderExceptions

class TestGuiderImage(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        self.verbose = True
        # calibration files
        self.inDarkFile = 'gcam/gimg-0001.fits.gz'
        self.outDarkFile = 'gcam/proc-'+self.inDarkFile
        self.inFlatFile = 'gcam/gimg-0003.fits.gz'
        self.outFlatFile = 'gcam/proc-'+self.inFlatFile
        self.inFlatEcamFile = 'ecam/gimg-0003.fits.gz'
        self.outFlatEcamFile = 'ecam/proc-'+self.inFlatFile

        # bad files
        self.saturatedFile = 'gcam/gimg-0004.fits.gz'
        self.badReadFile1 = 'gcam/gimg-0859.fits.gz'
        self.badReadFile2 = 'gcam/gimg-0519.fits.gz'

        # good data files
        self.inDataFile = 'gcam/gimg-0040.fits.gz'
        self.outDataFile = 'gcam/proc-'+self.inFlatFile

        # values for comparison after processing good data files.
        self.inDataResults = 'something.fits'
        super(TestGuiderImage,self).setUp()

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
        self.gi.analyzeDark(inFile, cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile),'analyzeDark file write')
        self._check_overwriting(inFile,outFile,self.gi.analyzeDark)
    
    def test_analyzeFlat(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""
        inFile = self.inFlatFile
        outFile = self.outFlatFile
        self.gi.analyzeFlat(inFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile),'analyzeFlat file write')
        # TBD: some way to test that the correct fibers are identified?
        for name,i in self.probeNames.items():
        #    #self.assertEqual(i,self.gi.flatFibers[i].fiberid)
            print i,self.gi.flatFibers[i].fiberid
        self._check_overwriting(inFile,outFile,self.gi.analyzeFlat,[self.gState.gprobes])

    def test_analyzeFlat_ecam(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""
        inFile = self.inFlatEcamFile
        outFile = self.outFlatEcamFile
        self.gi.camera = 'ecamera'
        self.gi.analyzeFlat(inFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(outFile),'analyzeFlat file write')
        # TBD: best test is probably to check that the flat is median ~1.
        self.assertFail('Need to create a test for this!')
        self._check_overwriting(inFile,outFile,self.gi.analyzeFlat,[self.gState.gprobes])

    def test_call(self):
        """Test GuiderImageAnalysis.__call__()"""
        fibers = self._call_gi(self.inDataFile)


        for name,i in self.probeNames.items():
            self.assertEqual(i,fibers[i].fiberid-1)
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
        """Test what happens when the setPoint `is None."""
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
    unittest.main(verbosity=2)

