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
import unittest

from actorcore import TestHelper

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState
import guiderActor.myGlobals as myGlobals
from guiderActor.GuiderActor import set_default_pids, set_pid_scaling

gprobeKey = {}
guideInfoKey = {}

# These keys are for the guider files for cart 13, plate 7660 on MJD 57356.
platedb_gprobe = [[13,1,True,226.00,448.00,8.50,58.00,-8.40,-9.80,400.00,'GUIDE'],
    [13,2,True,159.50,309.50,8.50,44.00,-10.60,13.10,0.00,'GUIDE'],
    [13,3,True,129.50,139.50,28.50,41.00,-16.30,6.90,0.00,'ACQUIRE'],
    [13,4,True,226.00,102.50,8.50,67.00,2.20,-2.30,-400.00,'GUIDE'],
    [13,5,True,226.00,275.50,8.50,329.00,13.30,7.00,0.00,'GUIDE'],
    [13,6,True,226.00,188.50,8.50,76.00,8.30,-6.40,-400.00,'GUIDE'],
    [13,7,True,92.50,243.00,8.50,343.00,-2.10,-11.40,0.00,'GUIDE'],
    [13,8,True,226.00,362.00,8.50,3.00,6.90,-2.70,400.00,'GUIDE'],
    [13,9,True,303.50,362.00,8.50,328.00,-5.50,2.70,400.00,'GUIDE'],
    [13,10,True,303.50,275.50,8.50,80.00,11.10,-1.00,0.00,'GUIDE'],
    [13,11,True,401.00,139.50,28.50,236.00,-10.90,-15.70,0.00,'ACQUIRE'],
    [13,12,True,303.50,188.50,8.50,18.00,-4.90,7.00,-400.00,'GUIDE'],
    [13,13,True,371.00,400.50,8.50,357.00,8.80,-2.00,0.00,'GUIDE'],
    [13,14,True,303.50,102.50,8.50,65.00,13.30,3.00,-400.00,'GUIDE'],
    [13,15,True,438.00,333.00,8.50,5.00,-3.80,-4.80,0.00,'GUIDE'],
    [13,16,True,303.50,448.00,8.50,61.00,8.40,2.10,400.00,'GUIDE'],
    [13,17,True,129.50,411.00,1.00,0.00,0.00,0.00,0.00,'TRITIUM']]
platedb_gprobesInUse = ["(1=0x0)","(2=0x0)","(3=0x0)","(4=0x0)","(5=0x0)","(6=0x0)","(7=0x0)","(8=0x0)","(9=0x0)","(10=0x0)","(11=0x0)","(12=0x0)","(13=0x0)","(14=0x0)","(15=0x0)","(16=0x0)","(17=0x2)"]
platedb_plPlugMapM = [7660,57356,1]
platedb_guideInfo = [[1,6.7933,25.9177,-269.7172,11.4905,247.181,0.00],
    [2,7.4312,25.1582,-145.5555,-154.9257,172.294,0.00],
    [3,8.3402,25.9051,33.3775,7.3405,245.313,0.00],
    [4,9.0094,24.9979,165.8246,-189.7583,156.619,0.00],
    [5,9.3816,25.3645,238.6046,-109.3392,192.807,0.00],
    [6,7.7632,24.6309,-80.4990,-270.2283,183.272,0.00],
    [7,8.2122,24.6483,8.4350,-266.4980,184.951,0.00],
    [8,8.8256,24.6050,130.0581,-275.7856,180.771,0.00],
    [9,9.6424,25.9622,288.6084,21.4059,52.3773,0.00],
    [10,8.7879,26.5346,120.4569,144.7009,-3.10519,0.00],
    [11,8.2969,26.3027,24.8438,93.8529,19.7762,0.00],
    [12,7.5231,26.7023,-125.7832,181.2737,-19.563,0.00],
    [13,7.0452,26.2957,-219.5712,93.3688,19.994,0.00],
    [14,8.9030,27.0475,142.3818,256.7546,9.33557,0.00],
    [15,7.9536,26.8680,-41.9554,217.0998,27.1801,0.00],
    [16,7.6758,27.0000,-95.8673,246.0751,14.1412,0.00]]


gprobeKey['guide'] = [5,1,True,216.00,428.00,8.50,51.00,-6.30,-2.80,400.00,'GUIDE']
guideInfoKey['guide'] = [1,329.3098,39.483133,-223.103,-86.4817,203.093,0.00]
gprobeKey['guide_disabled'] = [5,4,True,216.00,82.50,8.50,169.00,-11.40,-12.30,-400.00,'GUIDE']
guideInfoKey['guide_disabled'] = [4,329.3098,39.483133,162.202,-143.16,177.588,0.00]

gprobeKey['guide_broken'] = [11,3,False,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'ACQUIRE']
guideInfoKey['guide_broken'] = [3,329.3098,39.483133,-32.3079,-59.4912,215.239,0.00]

gprobeKey['acquire'] = [11,11,True,391.00,119.50,28.50,257.00,15.80,12.30,0.00,'ACQUIRE']
guideInfoKey['acquire'] = [11,329.3098,39.483133,5.5082,31.1517,47.9917,0.00]
gprobeKey['acquire_disabled'] = [11,3,True,119.50,119.50,28.50,287.00,15.90,-8.70,0.00,'ACQUIRE']
guideInfoKey['acquire_disabled'] = [3,329.3098,39.483133,-32.3079,-59.4912,215.239,0.00]


def updateModel(name,model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name,model)

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
        self.gi = guiderImage.GuiderImageAnalysis(self.setPoint_good,self.actorState.actor.location)
        gState = GuiderState.GuiderState()
        myGlobals.actorState.gState = gState

        self.config = ConfigParser.ConfigParser()
        self.config.read(os.path.expandvars('$GUIDERACTOR_DIR/etc/guider.cfg'))
        set_default_pids(self.config, gState)
        set_pid_scaling(self.config, gState)

        self.probeNames = {}
        for name in gprobeKey:
            gk = gprobeKey[name]
            gik = guideInfoKey[name]
            self.probeNames[name] = gk[1]
            gState.gprobes[gk[1]] = GuiderState.GProbe(gk[1],gprobeKey=gk,guideInfoKey=gik)
            if 'disabled' in name:
                gState.gprobes[gk[1]].disabled = True
        self.gState = gState
    
    def _init_probes(self):
        """
        Initialize the "real" guide probes.
        Use with tests that e.g. call guiderImageAnalysis on the MJD 57356 files.
        """

        for probe,info,bits in zip(platedb_gprobe,platedb_guideInfo,platedb_gprobesInUse):
            junk,bits = bits.strip('()').split('=')
            self.gState.gprobes[probe[1]] = GuiderState.GProbe(probe[1],gprobeKey=probe,guideInfoKey=info,gprobeBits=int(bits,16))
        tritium = platedb_gprobe[16]
        self.gState.gprobes[tritium[1]] = GuiderState.GProbe(tritium[1],gprobeKey=tritium,gprobeBits=2)

    def _call_gi(self,filename,setpoint=None,args=[]):
        """Use this to simplify calls to guiderImageAnalysis."""
        if setpoint is None: setpoint = self.setPoint_good
        return self.gi(self.cmd,filename,self.gState.gprobes,setpoint,*args)

    def _remove_file(self,filename):
        if os.path.exists(filename):
            os.remove(filename)


class GuiderThreadTester(GuiderTester,unittest.TestCase):
    """
    guiderActor Thread test suites should subclass this and unittest, in that order.
    """
    def __init__(self, *args, **kwargs):
        """Load up the cmd calls for this test class, so _check_cmd will use them."""
        unittest.TestCase.__init__(self, *args, **kwargs)
        # -1 is the test function, -2 is test class, -3 (or 0) should be main
        class_name = self.id().split('.')[-2]
        self._load_cmd_calls(class_name)
        # lets us see really long list/list diffs
        self.maxDiff = None

