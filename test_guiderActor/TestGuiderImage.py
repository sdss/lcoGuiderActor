#!/usr/bin/env python
"""
Test the behavior of guider flats, including finding fibers..
"""
import os
import unittest
import re
import pyfits

import guiderTester

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState
from guiderActor.gimg import GuiderExceptions

class TestGuiderImage(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        
        self.inDarkFile = 'gimg-0001.fits.gz'
        self.outDarkFile = 'proc-'+self.inDarkFile
        self.inFlatFile = 'gimg-0003.fits.gz'
        self.outFlatFile = 'proc-'+self.inFlatFile
        self.inDataFile = 'gimg-0040.fits.gz'
        self.outDataFile = 'proc-'+self.inFlatFile
        self.saturatedFile = 'gimg-0004.fits.gz'
        self.badReadFile1 = 'gimg-0859.fits.gz'
        self.badReadFile2 = 'gimg-0519.fits.gz'
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
        self.gi.analyzeDark(self.inDarkFile,cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outDarkFile),'analyzeDark file write')
        self._check_overwriting(self.inDarkFile,self.outDarkFile,self.gi.analyzeDark)
    
    def test_analyzeFlat(self):
        """Test GuiderImageAnalysis.analyzeFlat()"""
        self.gi.analyzeFlat(self.inFlatFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outFlatFile),'analyzeFlat file write')
        # TBD: some way to test that the correct fibers are identified?
        for name,i in self.probeNames.items():
        #    #self.assertEqual(i,self.gi.flatFibers[i].fiberid)
            print i,self.gi.flatFibers[i].fiberid
        self._check_overwriting(self.inFlatFile,self.outFlatFile,self.gi.analyzeFlat,[self.gState.gprobes])

    def test_call(self):
        """Test GuiderImageAnalysis.__call__()"""
        fibers = self._call_gi(self.inDataFile)
        for name,i in self.probeNames.items():
            self.assertEqual(i,fibers[i].fiberid-1)
        self.assertFalse(True,'make a test!')
    
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
    
    def test_badSetPoint_image(self):
        """Test what happens when the ccdtemp is outside the dark temp spec for an image."""
        self.gi.darkTemperature = self.setPoint_bad
        self._temp_run(self.inDataFile,'CCD temp signifcantly different (>3.0) from dark temp: -40.1, expected -35.0')

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
#...

if __name__ == '__main__':
    unittest.main()

