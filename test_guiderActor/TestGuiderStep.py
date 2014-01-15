#!/usr/bin/env python
"""
Test the phases of guiderStep in masterThread.py
"""
import unittest
import os

import guiderTester

from guiderActor.gimg import guiderImage
from guiderActor import masterThread, GuiderState

class TestGuiderStep(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        self.centerUpIn = 'gimg-0024.fits.gz'
        self.centerUpOut = 'proc-'+self.centerUpIn
        self.guidingIn = 'gimg-0040.fits.gz'
        self.guidingOut = 'proc-'+self.guidingIn
        super(TestGuiderStep,self).setUp()
    
    def tearDown(self):
        self._remove_file(self.centerUpOut)
        self._remove_file(self.guidingOut)
        super(TestGuiderStep,self).tearDown()
    
    def test_check_fiber_guiding(self):
        self.gState.centerUp = False
        self.fibers = self.gi(self.cmd,self.guidingIn,self.gState.gprobes,setPoint=self.setPoint_good)
        for name,i in self.probeNames.items():
            probe = self.gState.gprobes[i]
            fiber = self.fibers[i-1]
            result = masterThread._check_fiber(fiber,self.gState,self.cmd)
            #print name,i,result
            if 'disabled' in name or 'broken' in name:
                self.assertFalse(result,name)
            else:
                self.assertTrue(result,name)
    
    def test_check_fiber_centerUp(self):
        self.gState.centerUp = True
        self.fibers = self.gi(self.cmd,self.centerUpIn,self.gState.gprobes,setPoint=self.setPoint_good)
        for name,i in self.probeNames.items():
            probe = self.gState.gprobes[i]
            fiber = self.fibers[i-1]
            result = masterThread._check_fiber(fiber,self.gState,self.cmd)
            #print name,i,result
            if 'disabled' in name or 'broken' in name or 'acquire' not in name:
                self.assertFalse(result,name)
            else:
                self.assertTrue(result,name)
#...

if __name__ == '__main__':
    unittest.main()
