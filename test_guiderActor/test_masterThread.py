#!/usr/bin/env python
"""
Test the behavior of the various guider master thread commands.
"""
import unittest
from actorcore import TestHelper

import guiderTester

from guiderActor import masterThread

# from guiderActor import GuiderState
# from guiderActor.gimg import GuiderExceptions


class TestMasterThread(guiderTester.GuiderTester,unittest.TestCase):
    """Test specific masterThread commands."""
    def setUp(self):
        super(TestMasterThread,self).setUp()
    
    def test_load_cartridge(self):
        pass

class TestGuiderStep(guiderTester.GuiderTester,unittest.TestCase):
    """Test the phases of guiderStep in masterThread"""
    def setUp(self):
        self.centerUpIn = 'data/gimg-0024.fits.gz'
        self.centerUpOut = 'data/proc-'+self.centerUpIn
        self.guidingIn = 'data/gimg-0040.fits.gz'
        self.guidingOut = 'data/proc-'+self.guidingIn
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

class TestSetRefraction(guiderTester.GuiderTester,unittest.TestCase):
    def _set_refraction(self, corrRatio, plateType, surveyMode, expect):
        self.gState.refractionBalance = -100 # ensure it's always different to start.
        masterThread.set_refraction(self.cmd, self.gState, corrRatio, plateType, surveyMode)
        self.assertEqual(self.gState.refractionBalance, expect)
    def test_corrRatio_1(self):
        self._set_refraction(1, None, None, 1)
    def test_apogee(self):
        self._set_refraction(None, 'APOGEE', None, 1)
    def test_boss(self):
        self._set_refraction(None, 'eBOSS', None, 0)


class TestGuidingIsOK(guiderTester.GuiderTester,unittest.TestCase):
    def _guidingIsOK(self, expect, nWarn=0, force=False):
        result = masterThread.guidingIsOK(self.cmd, self.actorState, force=force)
        self.assertEqual(result, expect)
        self._check_levels(0, 0, nWarn, 0)

    def test_force(self):
        self._guidingIsOK(True, force=True)

    def test_boss_science(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(True)

    def test_ffs_closed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['all_off'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_ffs_closed_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['all_off'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs'])
        self._guidingIsOK(True, 1)

    def test_fflamp_on(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['flats'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_fflamp_on_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['flats'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs','lamp_ff'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(True, 1)
    def test_arclamps_on(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['arcs'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_arclamps_on_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['arcs'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs','lamp_hgcd','lamp_ne'])
        self._guidingIsOK(True, 1)

    def test_tcc_halted(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['halted'])
        self._guidingIsOK(False, 1)
    def test_tcc_motion_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['halted'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['axes'])
        self._guidingIsOK(True, 1)

class TestEcam(guiderTester.GuiderTester,unittest.TestCase):
    """Tests for turning ecam processing on and off."""
    def test_ecam_on(self):
        self.fail()

if __name__ == '__main__':
    unittest.main(verbosity=2)
