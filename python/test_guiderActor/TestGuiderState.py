#!/usr/bin/env python
"""
Test the GuiderState module using unittest.
"""
import unittest2 as unittest
from guiderActor import GuiderState
#from guiderActor import guiderActor_main
import numpy as np

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
            if name == 'aboveFocus':
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_aboveFocus(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.belowFocus
            if name == 'belowFocus':
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_atFocus(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.atFocus
            if name == 'aboveFocus' or name == 'belowFocus':
                self.assertFalse(value,name)
            else:
                self.assertTrue(value,name)
    def test_exists(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.exists
            if name != 'broken':
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
    def test_enabled(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.enabled
            if name == 'tritium' or name == 'broken':
                self.assertFalse(value,name)
            else:
                self.assertTrue(value,name)
    def test_disabled(self):
        for name in gprobeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            value = probe.disabled
            if name == 'tritium' or name == 'broken':
                self.assertTrue(value,name)
            else:
                self.assertFalse(value,name)
#...

if __name__ == '__main__':
    unittest.main()
