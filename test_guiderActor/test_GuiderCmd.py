"""
Tests of the GuiderCmd functions, without executing any threads.

Each of these tests should confirm that a Guidermd command call calls the
correct queue with the correct parameters.

If these tests work correctly, each masterThread function should work
correctly when called via a GuiderCmd (assuming test_masterThread clears).
"""
import os
import unittest
from Queue import Queue

import guiderActor
import guiderActor.myGlobals as myGlobals
import guiderTester
from actorcore import TestHelper
from guiderActor.Commands import GuiderCmd


class GuiderCmdTester(guiderTester.GuiderTester, unittest.TestCase):
    def setUp(self):
        self.verbose = True
        super(GuiderCmdTester, self).setUp()
        self.timeout = 1
        # Do this after super setUp, as that's what creates actorState.
        self.queues = {}
        self.queues[guiderActor.MASTER] = Queue('master')
        myGlobals.actorState.queues = self.queues
        self.guiderCmd = GuiderCmd.GuiderCmd(self.actor)


class GuiderAPOCmdTester(GuiderCmdTester):
    def setUp(self):
        self.name = 'guider'
        self.actor = TestHelper.FakeActor(self.name, self.name + 'Actor', location='APO')
        super(GuiderAPOCmdTester, self).setUp()
        self.guidercmd = self.actor.commandSets['GuiderCmd_APO']

    def test_ping(self):
        """Ping just finishes."""
        self._run_cmd('ping', None)
        self._check_cmd(0, 1, 0, 0, True)


class GuiderLCOCmdTester(GuiderCmdTester):
    def setUp(self):
        self.name = 'guider'
        self.actor = TestHelper.FakeActor(self.name, self.name + 'Actor', location='LCO')
        super(GuiderLCOCmdTester, self).setUp()
        self.guidercmd = self.actor.commandSets['GuiderCmd_LCO']

    def test_ping(self):
        """Ping just finishes."""
        self._run_cmd('ping', None)
        self._check_cmd(0, 1, 0, 0, True)


class GuiderLocalCmdTester(GuiderCmdTester):
    def setUp(self):
        self.name = 'guider'
        self.actor = TestHelper.FakeActor(self.name, self.name + 'Actor', location='LOCAL')
        super(GuiderLocalCmdTester, self).setUp()
        self.guidercmd = self.actor.commandSets['GuiderCmd_LOCAL']

    def test_ping(self):
        """Ping just finishes."""
        self._run_cmd('ping', None)
        self._check_cmd(0, 1, 0, 0, True)


class TestDecentering(GuiderCmdTester, unittest.TestCase):
    """mangaDither, decenter on/off, setDecenter."""

    def _mangaDither(self, expect, args):
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('mangaDither %s' % (args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.DECENTER)
        self.assertEqual(msg.decenters['decenterRA'], expect['decenterRA'])
        self.assertEqual(msg.decenters['decenterDec'], expect['decenterDec'])
        self.assertEqual(msg.decenters['decenterRot'], expect['decenterRot'])
        self.assertEqual(msg.decenters['mangaDither'], expect['mangaDither'])

    def test_mangaDither_C(self):
        expect = {'decenterRA': 0,
                  'decenterDec': 0,
                  'decenterRot': 0,
                  'mangaDither': 'C'}
        self._mangaDither(expect, 'ditherPos=C')

    def test_mangaDither_N(self):
        expect = {'decenterRA': -0.417,
                  'decenterDec': 0.721,
                  'decenterRot': 0,
                  'mangaDither': 'N'}
        self._mangaDither(expect, 'ditherPos=N')

    def test_mangaDither_S(self):
        expect = {'decenterRA': -0.417,
                  'decenterDec': -0.721,
                  'decenterRot': 0,
                  'mangaDither': 'S'}
        self._mangaDither(expect, 'ditherPos=S')

    def test_mangaDither_E(self):
        expect = {'decenterRA': 0.833,
                  'decenterDec': 0,
                  'decenterRot': 0,
                  'mangaDither': 'E'}
        self._mangaDither(expect, 'ditherPos=E')

    def _decenter(self, expect, args):
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('decenter %s' % (args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.DECENTER)
        self.assertEqual(msg.enable, expect['on'])
        self.assertIsNone(getattr(msg, 'finish', None))

    def test_decenter_on(self):
        expect = {'on': True}
        self._decenter(expect, 'on')

    def test_decenter_off(self):
        expect = {'on': False}
        self._decenter(expect, 'off')

    def _setDecenter(self, expect, args, didFail=False):
        queue = self.queues[guiderActor.MASTER]
        if didFail:
            with self.assertRaises(AttributeError):
                msg = self._run_cmd('setDecenter %s' % (args), None, empty=True)
                self.assertEquals(msg, None)
                self._check_cmd(0, 0, 0, 0, True, True)
                self.assertTrue(msg.finish)
        else:
            msg = self._run_cmd('setDecenter %s' % (args), queue)
            self.assertEqual(msg.type, guiderActor.Msg.DECENTER)
            self.assertIsNone(getattr(msg, 'finish', None))
            self.assertEqual(msg.decenters['decenterRA'], expect.get('decenterRA', 0))
            self.assertEqual(msg.decenters['decenterDec'], expect.get('decenterDec', 0))

    def test_setDecenter_ra(self):
        expect = {'decenterRA': 10}
        self._setDecenter(expect, 'decenterRA=10')

    def test_setDecenter_dec(self):
        expect = {'decenterDec': 10}
        self._setDecenter(expect, 'decenterDec=10')

    def test_setDecenter_rot(self):
        self._setDecenter({}, 'decenterRot=10', didFail=True)


class TestSetRefractionBalance(GuiderCmdTester, unittest.TestCase):
    def _setRefractionBalance(self, args, expect):
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('setRefractionBalance %s' % (args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.SET_REFRACTION)
        self.assertEqual(msg.corrRatio, expect.get('corrRatio', None))
        self.assertEqual(msg.plateType, expect.get('plateType', None))
        self.assertEqual(msg.surveyMode, expect.get('surveyMode', None))

    def test_corrRatio_1(self):
        args = 'corrRatio=1'
        expect = {'corrRatio': 1}
        self._setRefractionBalance(args, expect)

    def test_survey(self):
        args = 'plateType="APOGEE-2&MaNGA" surveyMode="APOGEE lead" '
        expect = {'plateType': 'APOGEE-2&MaNGA',
                  'surveyMode': 'APOGEE lead'}
        self._setRefractionBalance(args, expect)


class TestFakeCartridge(GuiderCmdTester, unittest.TestCase):
    testAttrs = [
        'ra',
        'dec',
        'xFocal',
        'yFocal',
        'xCenter',
        'yCenter',
        'radius',
        'rotation',
        'xFerruleOffset',
        'yFerruleOffset',
        'focusOffset',
        'fiber_type']

    def setUp(self):
        self.name = 'guider'
        self.actor = TestHelper.FakeActor(self.name, self.name + 'Actor', location='LCO')
        super(TestFakeCartridge, self).setUp()

    def _fakeCartridge(self, args, expect=None):
        if expect is None:
            expect = {}
        queue = self.queues[guiderActor.MASTER]
        pmDir = os.path.expandvars('$GUIDERACTOR_DIR/test_guiderActor/gcam')
        msg = self._run_cmd('fakeCartridge pmDir=%s %s' % (pmDir, args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.LOAD_CARTRIDGE)
        self.assertEqual(msg.plate, expect.get('plate', None))
        self.assertEqual(msg.pointing, expect.get('pointing', None))
        for gprobeGot in msg.gprobes.itervalues():
            for gprobeId, gprobeExpect in expect['gprobes'].iteritems():
                if gprobeGot.id == gprobeId:
                    for attr in self.testAttrs:
                        if attr in gprobeExpect:
                            if isinstance(getattr(gprobeGot, attr), str):
                                self.assertEqual(
                                    getattr(
                                        gprobeGot, attr), gprobeExpect.get(
                                        attr, None), '{} for gprobeId: {}'.format(
                                        attr, gprobeId))
                            else:
                                # NOTE: grumble, float comparisons, grumble...
                                self.assertAlmostEqual(
                                    getattr(
                                        gprobeGot, attr), gprobeExpect.get(
                                        attr, None), 6, msg='{} for gprobeId: {}'.format(
                                        attr, gprobeId))

    def test_fakeCartridge_8641(self):
        # gProbe1 = guiderActor.GuiderState.GProbe(
        #     gprobeKey = [1, 1, 1, 216.0, 428.0, 8.5, 272.0, -4.5, 3.20000004768, 400.0, 'GUIDE'],
        #     guideInfoKey = [1, 282.49004, -41.867939, -120.56083, 287.07792, 0.0, 0.0]
        #     )
        gprobes = {
            1: {'xCenter': 170.5, 'yCenter': 412.5, 'ra': 282.49004, 'dec': -41.867939},
        }
        expect = {
            'plate': 8641,
            'pointing': 'A',
            'fiberPos': 1,
            'gprobes': gprobes}
        self._fakeCartridge(
            'plate={plate} pointing={pointing} fiberPos={fiberPos}'.format(
                **expect), expect)

    def test_fakeCartridge_8646(self):
        gprobes = {
            1: {'ra': 7.713199, 'dec': -29.765444, 'xFocal': -61.316109, 'yFocal': 253.6546, 'xCenter': 170.5, 'yCenter': 412.5,
                'radius': 8.5, 'rotation': 272.0, 'xFerruleOffset': -4.5, 'yFerruleOffset': 3.2, 'focusOffset': 400.0, 'fiber_type': 'GUIDE'},
            2: {'ra': 7.865471, 'dec': -29.11857, 'xFocal': -105.42248, 'yFocal': 39.318875, 'xCenter': 140.5, 'yCenter': 264.5,
                'radius': 8.5, 'rotation': 268.0, 'xFerruleOffset': -6.9, 'yFerruleOffset': 0.3, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            3: {'ra': 7.7258693, 'dec': -28.697817, 'xFocal': -65.42454, 'yFocal': -99.718591, 'xCenter': 154.5, 'yCenter': 90.5,
                'radius': 27, 'rotation': 132.0, 'xFerruleOffset': 9.2, 'yFerruleOffset': 1.7, 'focusOffset': 0.0, 'fiber_type': 'ACQUIRE'},
            4: {'ra': 7.9546143, 'dec': -28.722066, 'xFocal': -131.74361, 'yFocal': -91.585178, 'xCenter': 254.5, 'yCenter': 76.5,
                'radius': 8.5, 'rotation': 278.0, 'xFerruleOffset': -6.6, 'yFerruleOffset': 4.2, 'focusOffset': -400.0, 'fiber_type': 'GUIDE'},
            5: {'ra': 7.563513, 'dec': -28.37965, 'xFocal': -18.484631, 'yFocal': -205.19272, 'xCenter': 213.5, 'yCenter': 245.5,
                'radius': 8.5, 'rotation': 271.0, 'xFerruleOffset': -4.9, 'yFerruleOffset': 8.8, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            6: {'ra': 8.091459, 'dec': -29.331564, 'xFocal': -170.55195, 'yFocal': 110.10227, 'xCenter': 236.5, 'yCenter': 162.5,
                'radius': 8.5, 'rotation': 251.0, 'xFerruleOffset': 14.9, 'yFerruleOffset': -1.4, 'focusOffset': -400.0, 'fiber_type': 'GUIDE'},
            7: {'ra': 8.109592, 'dec': -29.023161, 'xFocal': -176.20022, 'yFocal': 8.1160002, 'xCenter': 90.5, 'yCenter': 182.5,
                'radius': 8.5, 'rotation': 267.0, 'xFerruleOffset': 3.1, 'yFerruleOffset': 6.4, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            8: {'ra': 8.0936182, 'dec': -28.763194, 'xFocal': -172.06145, 'yFocal': -77.863252, 'xCenter': 192.5, 'yCenter': 328.5,
                'radius': 8.5, 'rotation': 344.0, 'xFerruleOffset': -0.8, 'yFerruleOffset': -1.1, 'focusOffset': 400.0, 'fiber_type': 'GUIDE'},
            9: {'ra': 7.0428383, 'dec': -28.166214, 'xFocal': 133.73854, 'yFocal': -276.40956, 'xCenter': 268.5, 'yCenter': 347.5,
                'radius': 8.5, 'rotation': 264.0, 'xFerruleOffset': -9.3, 'yFerruleOffset': -0.8, 'focusOffset': 400.0, 'fiber_type': 'GUIDE'},
            10: {'ra': 6.955436, 'dec': -28.951372, 'xFocal': 157.45553, 'yFocal': -15.700456, 'xCenter': 289.5, 'yCenter': 264.5,
                 'radius': 8.5, 'rotation': 273.0, 'xFerruleOffset': -7.1, 'yFerruleOffset': 8.2, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            11: {'ra': 7.1099969, 'dec': -29.118803, 'xFocal': 112.50913, 'yFocal': 39.422145, 'xCenter': 418.5, 'yCenter': 156.5,
                 'radius': 27, 'rotation': 90.0, 'xFerruleOffset': -6.6, 'yFerruleOffset': -2.6, 'focusOffset': 0.0, 'fiber_type': 'ACQUIRE'},
            12: {'ra': 7.1206053, 'dec': -29.200405, 'xFocal': 109.37635, 'yFocal': 66.367385, 'xCenter': 308.5, 'yCenter': 182.5,
                 'radius': 8.5, 'rotation': 245.0, 'xFerruleOffset': -10.6, 'yFerruleOffset': -3.5, 'focusOffset': -400.0, 'fiber_type': 'GUIDE'},
            13: {'ra': 7.1268997, 'dec': -29.703818, 'xFocal': 107.35518, 'yFocal': 233.32247, 'xCenter': 324.5, 'yCenter': 402.5,
                 'radius': 8.5, 'rotation': 224.0, 'xFerruleOffset': -7.5, 'yFerruleOffset': 6.6, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            14: {'ra': 6.524485, 'dec': -28.842894, 'xFocal': 283.34771, 'yFocal': -50.92314, 'xCenter': 329.5, 'yCenter': 97.5,
                 'radius': 8.5, 'rotation': 253.0, 'xFerruleOffset': -8.3, 'yFerruleOffset': 5.0, 'focusOffset': -400.0, 'fiber_type': 'GUIDE'},
            15: {'ra': 6.444669, 'dec': -28.884399, 'xFocal': 306.63774, 'yFocal': -36.988933, 'xCenter': 407.5, 'yCenter': 354.5,
                 'radius': 8.5, 'rotation': 198.0, 'xFerruleOffset': -1.1, 'yFerruleOffset': -0.9, 'focusOffset': 0.0, 'fiber_type': 'GUIDE'},
            16: {'ra': 6.4775823, 'dec': -29.128048, 'xFocal': 296.27414, 'yFocal': 43.766755, 'xCenter': 293.5, 'yCenter': 428.0,
                 'radius': 8.5, 'rotation': 0.0, 'xFerruleOffset': -6.5, 'yFerruleOffset': -12.7, 'focusOffset': 400.0, 'fiber_type': 'GUIDE'},
        }
        # apply the "Harding Rotation Factor"
        for g in gprobes:
            gprobes[g]['rotation'] += 18.
        expect = {
            'plate': 8646,
            'pointing': 'D',
            'fiberPos': 3,
            'gprobes': gprobes}
        self._fakeCartridge(
            'plate={plate} pointing={pointing} fiberPos={fiberPos}'.format(
                **expect), expect)


class TestGuideOnOff(GuiderCmdTester, unittest.TestCase):
    def _guideOn(self, args, expect=None):
        if expect is None:
            expect = {}
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('on %s' % (args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.START_GUIDING)
        self.assertEqual(msg.expTime, expect.get('time', None))
        self.assertEqual(msg.stack, expect.get('stack', 1))
        self.assertEqual(msg.oneExposure, expect.get('oneExposure', False))
        self.assertEqual(msg.force, expect.get('force', False))
        self.assertEqual(msg.camera, expect.get('camera', 'gcamera'))

    def test_guideOn(self):
        self._guideOn('')

    def test_guideOn_force(self):
        self._guideOn('force', {'force': True})

    def test_guideOn_time(self):
        self._guideOn('time=15', {'time': 15})

    def test_guideOn_oneExposure(self):
        self._guideOn('oneExposure', {'oneExposure': True})

    def test_guideOn_stack(self):
        self._guideOn('stack=3', {'stack': 3})

    def test_guideOn_ecam(self):
        self.gState.plateType = 'ecamera'
        self._guideOn('', {'camera': 'ecamera'})

    def test_guideOff(self):
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('off', queue)
        self.assertEqual(msg.type, guiderActor.Msg.STOP_GUIDING)


class TestEcam(GuiderCmdTester, unittest.TestCase):
    def _findstar(self, args, expect={}):
        queue = self.queues[guiderActor.MASTER]
        msg = self._run_cmd('findstar %s' % (args), queue)
        self.assertEqual(msg.type, guiderActor.Msg.START_GUIDING)
        self.assertEqual(msg.expTime, expect.get('time', 5))
        self.assertEqual(msg.bin, expect.get('bin', 1))
        self.assertEqual(msg.oneExposure, True)
        self.assertEqual(msg.camera, 'ecamera')

    def test_findstar(self):
        self._findstar('')

    def test_findstar_nondefault(self):
        self._findstar('time=11 bin=2', {'time': 11, 'bin': 2})


if __name__ == '__main__':
    verbosity = 2

    suite = None
    # to test just one piece
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestClassifyCartridge)
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)
