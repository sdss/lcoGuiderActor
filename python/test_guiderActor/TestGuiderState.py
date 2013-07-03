#!/usr/bin/env python
"""
Test the GuiderState module using unittest.
"""
import unittest
import numpy as np

from guiderActor import GuiderState

gprobeKey = {}
gprobeKey['good'] = [10,1,True,100,100,10,0,-.1,-.1,0,'GUIDE']
gprobeKey['acquire'] = [10,2,True,200,200,10,0,-.1,-.1,0,'ACQUIRE']
gprobeKey['tritium'] = [10,3,True,300,300,10,0,-.1,-.1,0,'TRITIUM']
gprobeKey['aboveFocus'] = [10,4,True,400,400,10,0,-.1,-.1,-400,'GUIDE']
gprobeKey['belowFocus'] = [10,5,True,500,500,10,0,-.1,-.1,+400,'GUIDE']
gprobeKey['broken'] = [10,6,False,600,600,10,0,-.1,-.1,0,'GUIDE']

ugriz = np.array([13,13.5,14,14.5,15])

class Model(object):
    """quick replacement for Model in opscore/actorcore."""
    def __init__(self,key=None,value=None):
        self.keyVarDict = {}
        if key:
            self.keyVarDict[key] = value
#...
class State(object):
    """quick replacement for guiderActor_main.State"""
    def __init__(self,actor):
        self.actor = actor
        self.models = {}
#...

class TestGuiderState(unittest.TestCase):
    def setUp(self):
        self.gState = GuiderState.GuiderState()
        # need a tcc.axePos value
        #GuiderState.myGlobals.actorState = guiderActor_main.State(None)
        GuiderState.myGlobals.actorState = State(None)
        GuiderState.myGlobals.actorState.models = {'tcc':Model('axePos',[20,30,40])}
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
    def test_aboveFocus(self):
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
    
    def test_FrameInfo(self):
        plugPlateScale = 10.
        arcsecPerMM = 3600./plugPlateScale
        gcameraPixelSize = 15e-6
        gcameraMagnification = 10.
        guideCameraScale = gcameraMagnification * gcameraPixelSize * 1e-3
        frameInfo = GuiderState.FrameInfo(1,arcsecPerMM,guideCameraScale,plugPlateScale)
        self.assertAlmostEqual(frameInfo.micronsPerArcsec,plugPlateScale/3.6)
#...

if __name__ == '__main__':
    unittest.main()
