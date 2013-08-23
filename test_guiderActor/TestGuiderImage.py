#!/usr/bin/env python
"""
Test the behavior of guider flats, including finding fibers..
"""
import os
import unittest

import guiderTester

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState

class TestGuiderImage(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        self.inDarkFile = 'gimg-0001.fits.gz'
        self.outDarkFile = 'proc-'+self.inDarkFile
        self.inFlatFile = 'gimg-0002.fits.gz'
        self.outFlatFile = 'proc-'+self.inFlatFile
        self.inDataFile = 'gimg-0040.fits.gz'
        self.outDataFile = 'proc-'+self.inFlatFile
        super(TestGuiderImage,self).setUp()
    
    def test_analyzeDark(self):
        self.gi.analyzeDark(self.inDarkFile,cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outDarkFile),'analyzeDark file write')
        # If the proc- file exists, we should only read it, not rewrite it.
        # This fails if it tries to write to it.
        os.chmod(self.outDarkFile,0444)
        self.gi.analyzeDark(self.inDarkFile,cmd=self.cmd)
        os.chmod(self.outDarkFile,0644)
        os.remove(self.outDarkFile)
    
    def test_analyzeFlat(self):
        self.gi.analyzeFlat(self.inFlatFile,self.gState.gprobes,cmd=self.cmd)
        self.assertTrue(os.path.exists(self.outFlatFile),'analyzeFlat file write')
        
        # TBD: some way to test that the correct fibers are identified?
        for name,i in self.probeNames.items():
            self.assertEqual(i,self.gi.flatFibers[i].fiberid)
        
        # If the proc- file exists, we should only read it, not rewrite it.
        # This fails if it tries to write to it.
        os.chmod(self.outFlatFile,0444)
        self.gi.analyzeFlat(self.inFlatFile,self.gState.gprobes,cmd=self.cmd)
        os.chmod(self.outFlatFile,0644)
        os.remove(self.outFlatFile)

    def test_call(self):
        fibers = self.gi(self.cmd,self.inDataFile,self.gState.gprobes,self.setPoint_good)
        for name,i in self.probeNames.items():
            self.assertEqual(i,fibers[i].fiberid-1)
        
        self.assertFalse(True,'make a test!')

    def test_badSetPoint(self):
        self.gi(self.cmd,self.inDataFile,self.gState.gprobes,self.setPoint_bad)
        self.assertFalse(True,'make a test!')
#...

if __name__ == '__main__':
    unittest.main()

