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

import os
import ConfigParser

from actorcore import TestHelper

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState
import guiderActor.myGlobals as myGlobals

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


def updateModel(name,model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name,model)


def initialize_pid(config, gState):
    """Initialize the PID classes and set default PID values."""
    for axis in config.options('PID'):
        axis = dict(RADEC = "raDec", ROT = "rot", FOCUS = "focus", SCALE = "scale")[axis.upper()]
        Kp, Ti, Td, Imax, nfilt = [float(v) for v in config.get('PID', axis).split()]
        gState.set_pid_defaults(axis, Kp=Kp, Ti=Ti, Td=Td, Imax=Imax, nfilt=nfilt)
        gState.pid[axis].setPID(Kp=Kp, Ti=Ti, Td=Td, Imax=Imax, nfilt=nfilt)

class GuiderTester(TestHelper.ActorTester):
    """
    guiderActor test suites should subclass this and unittest, in that order.
    """
    def setUp(self):
        """Populate fake guide probes, etc."""
        self.verbose = True
        self.name = 'guider'
        super(GuiderTester,self).setUp()
        myGlobals.actorState = self.actorState
        self.setPoint_good = -40
        self.setPoint_bad = -35
        self.gi = guiderImage.GuiderImageAnalysis(self.setPoint_good)
        gState = GuiderState.GuiderState()

        self.config = ConfigParser.ConfigParser()
        self.config.read('../etc/guider.cfg')
        initialize_pid(self.config, gState)

        self.probeNames = {}
        for name in gprobeKey:
            gk = gprobeKey[name]
            gik = guideInfoKey[name]
            self.probeNames[name] = gk[1]
            gState.gprobes[gk[1]] = GuiderState.GProbe(gk[1],gprobeKey=gk,guideInfoKey=gik)
            if 'disabled' in name:
                gState.gprobes[gk[1]].disabled = True
        self.gState = gState
    
    def _call_gi(self,filename,setpoint=None,args=[]):
        """Use this to simplify calls to guiderImageAnalysis."""
        if setpoint is None: setpoint = self.setPoint_good
        return self.gi(self.cmd,filename,self.gState.gprobes,setpoint,*args)

    def _remove_file(self,filename):
        if os.path.exists(filename):
            os.remove(filename)
    
