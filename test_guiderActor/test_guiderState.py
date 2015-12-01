#!/usr/bin/env python
"""
Test the GuiderState module using unittest.
"""
import os
import unittest
import numpy as np

from actorcore import TestHelper

import guiderTester

from guiderActor import GuiderState
from guiderActor.GuiderActor import set_default_pids, set_pid_scaling

gprobeKey = {}
gprobeKey['good'] = [10,1,True,100,100,10,0,-.1,-.1,0,'GUIDE']
gprobeKey['acquire'] = [10,2,True,200,200,10,0,-.1,-.1,0,'ACQUIRE']
gprobeKey['tritium'] = [10,3,True,300,300,10,0,-.1,-.1,0,'TRITIUM']
gprobeKey['aboveFocus'] = [10,4,True,400,400,10,0,-.1,-.1,-400,'GUIDE']
gprobeKey['belowFocus'] = [10,5,True,500,500,10,0,-.1,-.1,+400,'GUIDE']
gprobeKey['broken'] = [10,6,False,600,600,10,0,-.1,-.1,0,'GUIDE']
# gprobesInUse = ["(1=0x0)","(2=0x0)","(3=0x2)","(4=0x0)","(5=0x0)","(6=0x1)"]
# NOTE: the above gprobesInUse would translate to these integer bits:
gprobebits = [0,0,2,0,0,1]
guideInfoKey = [1,10.,20.,30.0,40.0,50,0.00]

ugriz = np.array([13,13.5,14,14.5,15])

class TestGuiderState(guiderTester.GuiderTester,unittest.TestCase):
    def setUp(self):
        
        super(TestGuiderState,self).setUp()

        # have to set up the PID scaling values for APO
        self.config.read(os.path.expandvars('$GUIDERACTOR_DIR/etc/guider_APO.cfg'))
        set_default_pids(self.config, self.gState)
        set_pid_scaling(self.config, self.gState)

        #myGlobals.actorState.models = {'tcc':Model('axePos',[20,30,40])}
        for k,v in gprobeKey.items():
            guideInfoKey[0] = v[1]
            self.gState.gprobes[v[1]] = GuiderState.GProbe(gprobeKey=v,
                                                           guideInfoKey=guideInfoKey,
                                                           gprobeBits=gprobebits[v[1]-1])
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
            if 'broken' not in name and 'tritium' not in name:
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

    def _setRefractionBalance(self, plateType, surveyMode, expect):
        self.gState.setRefractionBalance(plateType, surveyMode)
        self.assertEqual(self.gState.refractionBalance, expect)

    def test_setRefractionBalance_eBOSS(self):
        self._setRefractionBalance('eBOSS', 'None', 0)
    def test_setRefractionBalance_APOGEE(self):
        self._setRefractionBalance('APOGEE', 'None', 1)
    def test_setRefractionBalance_APOGEE2(self):
        self._setRefractionBalance('APOGEE-2', 'None', 1)
    def test_setRefractionBalance_MaNGADither(self):
        self._setRefractionBalance('MaNGA', 'MaNGA Dither', 0)
    def test_setRefractionBalance_MaNGAStare(self):
        self._setRefractionBalance('MaNGA', 'MaNGA Stare', 0)
    def test_setRefractionBalance_ApogeeLead(self):
        self._setRefractionBalance('APOGEE-2&MaNGA', 'APOGEE lead', 1)
    def test_setRefractionBalance_ApogeeMangaDither(self):
        self._setRefractionBalance('APOGEE-2&MaNGA', 'MaNGA Dither', 0)
    def test_setRefractionBalance_ApogeeMangaStare(self):
        self._setRefractionBalance('APOGEE-2&MaNGA', 'MaNGA Stare', 0)

    def test_set_pid_defaults(self):
        values = dict(Kp=1, Ti=2, Td=3, Imax=4, nfilt=5)
        self.gState.set_pid_defaults('raDec',**values)
        for x in values:
            self.assertEqual(self.gState.pid_defaults['raDec'][x], values[x])

    def test_reset_pid_terms_all(self):
        for axis in ['raDec','rot','focus','scale']:
            self.gState.pid[axis]._x = 10000000
        self.gState.reset_pid_terms()
        for axis in ['raDec','rot','focus','scale']:
            self.assertIsNone(self.gState.pid[axis]._x, "%s was not reset"%axis)

    def test_reset_pid_terms_rot(self):
        for axis in ['raDec','rot','focus','scale']:
            self.gState.pid[axis]._x = 10000000
        axis = 'rot'
        self.gState.reset_pid_terms([axis])
        self.assertIsNone(self.gState.pid[axis]._x, "%s was not reset"%axis)
        for axis in ['raDec','focus','scale']:
            self.assertIsNotNone(self.gState.pid[axis]._x, "%s should not have reset"%axis)

    def _scale_pid_with_alt(self, alt, expect_Ti):
        result = self.gState.scale_pid_with_alt(alt)
        for axis in expect_Ti:
            self.assertEqual(self.gState.pid[axis].Ti, expect_Ti[axis], "%s does not match"%axis)
        return result
    def test_scale_pid_with_alt_low(self):
        expect_Ti = {'raDec':250, 'rot':250, 'scale':0, 'focus':0}
        result = self._scale_pid_with_alt(50, expect_Ti)
        self.assertFalse(result)
    def test_scale_pid_with_alt_mid(self):
        expect_Ti = {'raDec':175, 'rot':250, 'scale':0, 'focus':0}
        result = self._scale_pid_with_alt(70, expect_Ti)
        self.assertTrue(result)
    def test_scale_pid_with_alt_high(self):
        expect_Ti = {'raDec':100, 'rot':250, 'scale':0, 'focus':0}
        result = self._scale_pid_with_alt(85, expect_Ti)
        self.assertTrue(result)

    def test_update_pid_time_first_try(self):
        result = self.gState.update_pid_time('raDec',2)
        self.assertIsNone(result)
        self.assertEqual(self.gState.pid_time['raDec'], 2)

    def test_update_pid_time(self):
        self.gState.update_pid_time('raDec',2)
        result = self.gState.update_pid_time('raDec',3)
        self.assertEqual(result, 1)
        self.assertEqual(self.gState.pid_time['raDec'], 3)


    def test_output_pid(self):
        self.gState.output_pid(self.cmd)
        self._check_cmd(0,4,0,0,False)
#...

if __name__ == '__main__':
    unittest.main(verbosity=2)
