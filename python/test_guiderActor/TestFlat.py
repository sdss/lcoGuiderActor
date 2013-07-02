#!/usr/bin/env python
"""
Test the behavior of guider flats, including finding fibers..
"""
from guiderActor import gimg.guiderImage
from guiderActor import GuiderState

class Cmd(object):
    def inform(self,txt):
        print 'i',txt
    def diag(self,txt):
        print 'd',txt
    def warn(self,txt):
        print 'w',txt

gprobeKey1 = [11,1,True,216.00,428.00,8.50,237.00,-14.10,-3.30,400.00,'GUIDE']
guideInfoKey1 = [1,214.6303,52.5876,-247.4540,-70.2827,210.383,0.00]
gprobeKey2 = [11,11,True,391.00,119.50,28.50,329.00,5.80,0.80,0.00,'ACQUIRE']
guideInfoKey2 = [11,216.3289,53.1114,-22.3338,40.5621,43.757,0.00]

class TestGuiderState(unittest.TestCase):
    def setUp(self):
        self.cmd = Cmd()
        self.gi = guiderImage.GuiderImageAnalysis(-40)
        gState = GuiderState.GuiderState()
        gState.gprobes[1] = GuiderState.GProbe(1,gprobeKey=gprobeKey1,guideInfoKey=guideInfoKey1)
        gState.gprobes[2] = GuiderState.GProbe(2,gprobeKey=gprobeKey2,guideInfoKey=guideInfoKey2)
        self.gState = gState
    
    def test_analyzeFlat(self):
        gi.analyzeFlat('gimg-0002.fits.gz',self.gState.gprobes,cmd=self.cmd)
    
#...

if __name__ == '__main__':
    unittest.main()

