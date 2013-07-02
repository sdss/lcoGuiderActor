#!/usr/bin/env python
"""
Test the phases of guiderStep in masterThread.py
"""
from guiderActor import masterThread

class Cmd(object):
    def inform(self,txt):
        print 'i',txt
    def diag(self,txt):
        print 'd',txt
    def warn(self,txt):
        print 'w',txt

gProbeKey = {}
guideInfoKey = {}
gprobeKey['guide'] = [11,1,True,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey['guide'] = [1,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]
gprobeKey['guide_disabled'] = [11,2,True,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey['guide_disabled'] = [2,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]
gprobeKey['acquire'] = [11,11,True,391.00,119.50,28.50,329.00,5.80,0.80,0.00,'ACQUIRE']
guideInfoKey['acquire'] = [11,216.3289,53.1114,-22.3338,40.5621,43.757,0.00]
gprobeKey['acquire_disabled'] = [11,7,True,391.00,119.50,28.50,329.00,5.80,0.80,0.00,'ACQUIRE']
guideInfoKey['acquire_disabled'] = [7,216.3289,53.1114,-22.3338,40.5621,43.757,0.00]

class TestGuiderStep(unittest.TestCase):
    def setUp(self):
        self.cmd = Cmd()
        self.gi = guiderImage.GuiderImageAnalysis(-40)
        gState = GuiderState.GuiderState()
        for k in gprobeKey:
            gk = gprobeKey[k]
            gik = guideInfoKey[k])
            gState.gprobes[gk[1]] = GuiderState.GProbe(gk[1],gprobeKey=gk,guideInfoKey=gik]
        self.gState = gState
        
    def test_check_fiber_normal(self):
        self.gState.centerUp = False
        for name in gProbeKey:
            probe = self.gState.gprobes[gprobeKey[name][1]]
            masterThread._check_fiber(fiber,)
    
    def test_check_fiber_centerUp(self):
        self.gState.centerUp = True
#...
