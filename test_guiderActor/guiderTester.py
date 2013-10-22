"""
Helps with setting up and tearing down guideractor tests.

Example:
    class TestGuiderThing(guiderTester.GuiderTester,unittest.TestCase):
        def setUp(self):
            #do some stuff
            super(TestGuiderImage,self).setUp())
    
        def test_guiderImageThingy(self):
            self.call_gi(filename)
    
    if __name__ == '__main__':
        unittest.main()
"""

import unittest
import os

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState

class Cmd(object):
    def __init__(self):
        """Save the level of any messages that pass through."""
        self.levels = ''
        self.messages = []
    def _msg(self,txt,level):
        print level,txt
        self.levels += level
        self.messages.append(txt)
    def inform(self,txt):
        self._msg(txt,'i')
    def diag(self,txt):
        self._msg(txt,'d')
    def warn(self,txt):
        self._msg(txt,'w')
    def fail(self,txt):
        self._msg(txt,'f')
    def error(self,txt):
        self._msg(txt,'e')

gprobeKey = {}
guideInfoKey = {}

# TBD: need to create a fake tritium star too...

gprobeKey['guide'] = [11,1,True,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey['guide'] = [1,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]
gprobeKey['guide_disabled'] = [11,4,True,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey['guide_disabled'] = [4,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]

gprobeKey['guide_broken'] = [11,3,False,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey['guide_broken'] = [3,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]

gprobeKey['acquire'] = [11,11,True,391.00,119.50,28.50,329.00,5.80,0.80,0.00,'ACQUIRE']
guideInfoKey['acquire'] = [11,216.3289,53.1114,-22.3338,40.5621,43.757,0.00]
gprobeKey['acquire_disabled'] = [11,2,True,391.00,119.50,28.50,329.00,5.80,0.80,0.00,'ACQUIRE']
guideInfoKey['acquire_disabled'] = [2,216.3289,53.1114,-22.3338,40.5621,43.757,0.00]

class GuiderTester(object):
    """
    guiderActor test suites should subclass this and unittest, in that order.
    """
    def setUp(self):
        """Populate fake guide probes, etc."""
        self.setPoint_good = -40
        self.setPoint_bad = -35
        self.cmd = Cmd()
        self.gi = guiderImage.GuiderImageAnalysis(self.setPoint_good)
        gState = GuiderState.GuiderState()
        self.probeNames = {}
        for name in gprobeKey:
            gk = gprobeKey[name]
            gik = guideInfoKey[name]
            self.probeNames[name] = gk[1]
            gState.gprobes[gk[1]] = GuiderState.GProbe(gk[1],gprobeKey=gk,guideInfoKey=gik)
            if 'disabled' in name:
                gState.gprobes[gk[1]].disabled = True
        self.gState = gState
    
    def tearDown(self):
        pass
    
    def _call_gi(self,filename,setpoint=None,args=[]):
        """Use this to simplify calls to guiderImageAnalysis."""
        if setpoint is None: setpoint = self.setPoint_good
        return self.gi(self.cmd,filename,self.gState.gprobes,setpoint,*args)

    def _remove_file(self,filename):
        if os.path.exists(filename):
            os.remove(filename)
    
