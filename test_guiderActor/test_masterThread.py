#!/usr/bin/env python
"""
Test the behavior of the various guider master thread commands.
"""
import os
import unittest
import pyfits

import guiderTester

from guiderActor.gimg import guiderImage
from guiderActor import masterThread, GuiderState

from guiderActor.gimg import GuiderExceptions


class TestMasterThread(guiderTester.GuiderTester,unittest.TestCase):
    """Test specific masterThread commands."""
    def setUp(self):
        super(TestMasterThread,self).setUp()
    
    def test_load_cartridge(self):
        pass

class TestGuiderStep(guiderTester.GuiderTester,unittest.TestCase):
    """Test the phases of guiderStep in masterThread"""
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


class TestDecenter(guiderTester.GuiderTester,unittest.TestCase):
    def test_set_decenter_enable(self):
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        self.assertTrue(self.gState.decenter)
    def test_set_decenter_disable(self):
        decenters = {'decenterRA':1,'decenterDec':2,'mangaDither':'N'}
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        masterThread.set_decenter(self.cmd, {}, self.gState, False)
        self.assertFalse(self.gState.decenter)
        self.assertEqual(self.gState.decenterRA,0)
        self.assertEqual(self.gState.decenterDec,0)
        self.assertEqual(self.gState.mangaDither,'C')
    
    def _set_decenter_ok(self,decenters):
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        self.assertEqual(self.gState.decenterRA,decenters.get('decenterRA',0))
        self.assertEqual(self.gState.decenterDec,decenters.get('decenterDec',0))
        self.assertEqual(self.gState.mangaDither,decenters.get('mangaDither','?'))
        self.assertIn(self.cmd,self.gState.decenterCmd)
    def test_set_decenter_new(self):
        decenters = {'decenterRA':1,'decenterDec':2}
        self._set_decenter_ok(decenters)
    def test_set_decenter_mangaDither(self):
        decenters = {'decenterRA':1,'decenterDec':2,'mangaDither':'N'}
        self._set_decenter_ok(decenters)
    
    def test_set_decenter_cannot_change(self):
        decenters = {'decenterRA':1,'decenterDec':2}
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        self.assertTrue(self.cmd.didFail)
        self.assertFalse(self.gState.decenter)
        self.assertEqual(self.gState.decenterRA,0)
        self.assertEqual(self.gState.decenterDec,0)
        self.assertEqual(self.gState.mangaDither,'C')

if __name__ == '__main__':
    unittest.main(verbosity=2)
