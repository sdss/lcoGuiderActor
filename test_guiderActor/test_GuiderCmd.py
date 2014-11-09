"""
Tests of the GuiderCmd functions, without executing any threads.

Each of these tests should confirm that a Guidermd command call calls the
correct queue with the correct parameters.

If these tests work correctly, each masterThread function should work
correctly when called via a GuiderCmd (assuming test_masterThread clears).
"""
import unittest
from Queue import Queue

import guiderActor
import guiderActor.myGlobals as myGlobals
from guiderActor.Commands import GuiderCmd

from actorcore import TestHelper
import guiderTester

class GuiderCmdTester(guiderTester.GuiderTester):
    def setUp(self):
        self.verbose = True
        self.actor = TestHelper.FakeActor('guider','guiderActor')
        super(GuiderCmdTester,self).setUp()
        self.timeout = 1
        # Do this after super setUp, as that's what creates actorState.
        myGlobals.actorState.queues = {}
        myGlobals.actorState.queues[guiderActor.MASTER] = Queue('master')
        self.guiderCmd = GuiderCmd.GuiderCmd(self.actor)

class TestDecentering(GuiderCmdTester,unittest.TestCase):
    """mangaDither, decenter on/off, setDecenter."""
    def _mangaDither(self,expect,args):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        msg = self._run_cmd('mangaDither %s'%(args),queue)
        self.assertEqual(msg.type,guiderActor.Msg.DECENTER)
        self.assertEqual(msg.decenters['decenterRA'],expect['decenterRA'])
        self.assertEqual(msg.decenters['decenterDec'],expect['decenterDec'])
        self.assertEqual(msg.decenters['decenterRot'],expect['decenterRot'])
        self.assertEqual(msg.decenters['mangaDither'],expect['mangaDither'])
    def test_mangaDither_C(self):
        expect = {'decenterRA':0,
                  'decenterDec':0,
                  'decenterRot':0,
                  'mangaDither':'C'}
        self._mangaDither(expect,'ditherPos=C')
    def test_mangaDither_N(self):
        expect = {'decenterRA':-0.417,
                  'decenterDec':0.721,
                  'decenterRot':0,
                  'mangaDither':'N'}
        self._mangaDither(expect,'ditherPos=N')
    def test_mangaDither_S(self):
        expect = {'decenterRA':-0.417,
                  'decenterDec':-0.721,
                  'decenterRot':0,
                  'mangaDither':'S'}
        self._mangaDither(expect,'ditherPos=S')
    def test_mangaDither_E(self):
        expect = {'decenterRA':0.833,
                  'decenterDec':0,
                  'decenterRot':0,
                  'mangaDither':'E'}
        self._mangaDither(expect,'ditherPos=E')
    
    def _decenter(self,expect,args):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        msg = self._run_cmd('decenter %s'%(args),queue)
        self.assertEqual(msg.type,guiderActor.Msg.DECENTER)
        self.assertEqual(msg.enable,expect['on'])
        self.assertIsNone(getattr(msg,'finish',None))
    def test_decenter_on(self):
        expect = {'on':True}
        self._decenter(expect,'on')
    def test_decenter_off(self):
        expect = {'on':False}
        self._decenter(expect,'off')

    def _setDecenter(self,expect,args,didFail=False):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        if didFail:
            with self.assertRaises(AttributeError):
                msg = self._run_cmd('setDecenter %s'%(args), None, empty=True)
                self.assertEquals(msg,None)
                self._check_cmd(0,0,0,0, True, True)
                self.assertTrue(msg.finish)
        else:
            msg = self._run_cmd('setDecenter %s'%(args), queue)
            self.assertEqual(msg.type,guiderActor.Msg.DECENTER)
            self.assertIsNone(getattr(msg,'finish',None))
            self.assertEqual(msg.decenters['decenterRA'],expect.get('decenterRA',0))
            self.assertEqual(msg.decenters['decenterDec'],expect.get('decenterDec',0))
    def test_setDecenter_ra(self):
        expect = {'decenterRA':10}
        self._setDecenter(expect,'decenterRA=10')
    def test_setDecenter_dec(self):
        expect = {'decenterDec':10}
        self._setDecenter(expect,'decenterDec=10')
    def test_setDecenter_rot(self):
        self._setDecenter({},'decenterRot=10',didFail=True)

class TestSetRefractionBalance(GuiderCmdTester,unittest.TestCase):
    def _setRefractionBalance(self, args, expect):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        msg = self._run_cmd('setRefractionBalance %s'%(args),queue)
        self.assertEqual(msg.type, guiderActor.Msg.SET_REFRACTION)
        self.assertEqual(msg.corrRatio, expect.get('corrRatio',None))
        self.assertEqual(msg.plateType, expect.get('plateType',None))
        self.assertEqual(msg.surveyMode, expect.get('surveyMode',None))

    def test_corrRatio_1(self):
        args = 'corrRatio=1'
        expect = {'corrRatio':1}
        self._setRefractionBalance(args,expect)
    def test_survey(self):
        args = 'plateType="APOGEE-2&MaNGA" surveyMode="APOGEE lead" '
        expect = {'plateType':'APOGEE-2&MaNGA',
                  'surveyMode':'APOGEE lead'}
        self._setRefractionBalance(args,expect)

class TestEcam(GuiderCmdTester,unittest.TestCase):
    def _ecamOn(self, args, expect={}):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        msg = self._run_cmd('ecamOn %s'%(args),queue)
        self.assertEqual(msg.type, guiderActor.Msg.ECAM_ON)
        self.assertEqual(msg.time, expect.get('time',5))
    def test_ecamOn(self):
        self._ecamOn('')
    def test_ecamOn_time10(self):
        self._ecamOn('time=11', {'time':11})

    def test_ecamOff(self):
        queue = myGlobals.actorState.queues[guiderActor.MASTER]
        msg = self._run_cmd('ecamOff',queue)
        self.assertEqual(msg.type, guiderActor.Msg.ECAM_OFF)

if __name__ == '__main__':
    verbosity = 2
    
    suite = None
    # to test just one piece
    #suite = unittest.TestLoader().loadTestsFromTestCase(TestClassifyCartridge)
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)
