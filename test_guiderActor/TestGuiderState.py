#!/usr/bin/env python
"""
Test the GuiderState module using unittest.
"""
import unittest
import numpy as np

from actorcore import TestHelper

import guiderTester

from guiderActor import GuiderState

gprobeKey = {}
gprobeKey['good'] = [10,1,True,100,100,10,0,-.1,-.1,0,'GUIDE']
gprobeKey['acquire'] = [10,2,True,200,200,10,0,-.1,-.1,0,'ACQUIRE']
gprobeKey['tritium'] = [10,3,True,300,300,10,0,-.1,-.1,0,'TRITIUM']
gprobeKey['aboveFocus'] = [10,4,True,400,400,10,0,-.1,-.1,-400,'GUIDE']
gprobeKey['belowFocus'] = [10,5,True,500,500,10,0,-.1,-.1,+400,'GUIDE']
gprobeKey['broken'] = [10,6,False,600,600,10,0,-.1,-.1,0,'GUIDE']

ugriz = np.array([13,13.5,14,14.5,15])

class TestGuiderState(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        super(TestGuiderState,self).setUp()
        self.gState = GuiderState.GuiderState()
        #myGlobals.actorState.models = {'tcc':Model('axePos',[20,30,40])}
        for k,v in gprobeKey.items():
            self.gState.gprobes[v[1]] = GuiderState.GProbe(gprobeKey=v)
            self.gState.gprobes[v[1]].ugriz = ugriz
    
    def test_aboveFocus(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.aboveFocus
            if 'aboveFocus' in name:
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_belowFocus(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.belowFocus
            if 'belowFocus' in name:
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_atFocus(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.atFocus
            if 'aboveFocus' in name or 'belowFocus' in name:
                self.assertFalse(value,name)
            else:
                self.assertTrue(value,name)
    def test_exists(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.exists
            if 'broken' not in name:
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_enabled(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.enabled
            if 'tritium' in name or 'broken' in name:
                self.assertFalse(value,name)
            else:
                self.assertTrue(value,name)
    def test_disabled(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.disabled
            if 'tritium' in name or 'broken' in name:
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    
    def _setDecenter(self,decenters,enable=None,new=True):
        decenters = {'decenterRA':1,'decenterDec':2}
        self.gState.setDecenter(decenters,self.cmd,enable)
        self.assertEqual(self.gState.decenterRA,decenters.get('decenterRA',0))
        self.assertEqual(self.gState.decenterDec,decenters.get('decenterDec',0))
        self.assertEqual(self.gState.decenterRot,decenters.get('decenterRot',0))
        self.assertEqual(self.gState.mangaDither,decenters.get('mangaDither','?'))
        if new:
            self.assertIn(self.cmd,self.gState.decenterCmd)
        else:
            self.assertEqual(self.gState.decenterCmd,[])
    def test_setDecenter_on(self):
        self._setDecenter({},enable=True)
    def test_setDecenter_off(self):
        self.gState.setDecenter({},self.cmd,True)
        self._setDecenter({},enable=False)
    def test_setDecenter_new(self):
        decenters = {'decenterRA':1,'decenterDec':2}
        self._setDecenter(decenters)
    def test_setDecenter_mangaDither(self):
        decenters = {'decenterRA':1,'decenterDec':2,'mangaDither':'C'}
        self._setDecenter(decenters)
    def test_setDecenter_same(self):
        decenters = {'decenterRA':1,'decenterDec':2}
        self._setDecenter(decenters)
        self.gState.finish_decenter()
        # this prevents "this command has already finished"
        self.cmd = TestHelper.Cmd()
        self._setDecenter(decenters,new=False)
        self.assertTrue(self.cmd.finished)

    def test_FrameInfo(self):
        plugPlateScale = 10.
        arcsecPerMM = 3600./plugPlateScale
        gcameraPixelSize = 15e-6
        gcameraMagnification = 10.
        guideCameraScale = gcameraMagnification * gcameraPixelSize * 1e-3
        frameInfo = GuiderState.FrameInfo(1,arcsecPerMM,guideCameraScale,plugPlateScale)
        self.assertAlmostEqual(frameInfo.micronsPerArcsec,plugPlateScale/3.6)

    def _setRefractionBalance(self, survey, expect):
        self.gState.setRefractionBalance(survey)
        self.assertEqual(self.gState.refractionBalance, expect)

    def test_setRefractionBalance_eBOSS(self):
        self._setRefractionBalance('eBOSS', 0)
    def test_setRefractionBalance_APOGEE(self):
        self._setRefractionBalance('APOGEE', 1)
    def test_setRefractionBalance_APOGEE2(self):
        self._setRefractionBalance('APOGEE-2', 1)
    def test_setRefractionBalance_MaNGA(self):
        self._setRefractionBalance('MaNGA', 0)
    def test_setRefractionBalance_ApogeeManga(self):
        self._setRefractionBalance('APOGEE&MaNGA', 0)
#...

if __name__ == '__main__':
    unittest.main()
