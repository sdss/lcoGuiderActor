#!/usr/bin/env python
"""
Test the behavior of the various guider master thread commands.
"""
import unittest
from Queue import Queue

import guiderActor
import guiderActor.myGlobals as myGlobals
import guiderTester
from actorcore import Actor, TestHelper
from guiderActor import masterThread
from opscore.actor import KeyVarDispatcher, Model


class TestMasterThread(guiderTester.GuiderTester, unittest.TestCase):
    """Test specific masterThread commands."""

    def setUp(self):
        super(TestMasterThread, self).setUp()

    def test_load_cartridge(self):
        pass


@guiderTester.skipIfNoGuiderImages
class TestGuiderStep(guiderTester.GuiderTester, unittest.TestCase):
    """Test the phases of guiderStep in masterThread"""

    def setUp(self):
        self.centerUpIn = guiderTester.getTestImage('data', 57357,
                                                    'gimg-0024.fits.gz')
        self.centerUpOut = guiderTester.getTestImage('data', 57357,
                                                     'proc-' + self.centerUpIn)
        self.guidingIn = guiderTester.getTestImage('data', 57357,
                                                   'gimg-0040.fits.gz')
        self.guidingOut = guiderTester.getTestImage('data', 57357,
                                                    'proc-' + self.guidingIn)
        super(TestGuiderStep, self).setUp()

    def tearDown(self):
        self._remove_file(self.centerUpOut)
        self._remove_file(self.guidingOut)
        super(TestGuiderStep, self).tearDown()

    def test_check_fiber_guiding(self):
        self.gState.centerUp = False
        self.fibers = self.gi(
            self.cmd,
            self.guidingIn,
            self.gState.gprobes,
            setPoint=self.setPoint_good)
        for name, i in self.probeNames.items():
            # probe = self.gState.gprobes[i]
            fiber = self.fibers[i - 1]
            result = masterThread._check_fiber(fiber, self.gState, self.cmd)
            # print name,i,result
            if 'disabled' in name or 'broken' in name:
                self.assertFalse(result, name)
            else:
                self.assertTrue(result, name)

    def test_check_fiber_centerUp(self):
        self.gState.centerUp = True
        self.fibers = self.gi(
            self.cmd,
            self.centerUpIn,
            self.gState.gprobes,
            setPoint=self.setPoint_good)
        for name, i in self.probeNames.items():
            # probe = self.gState.gprobes[i]
            fiber = self.fibers[i - 1]
            result = masterThread._check_fiber(fiber, self.gState, self.cmd)
            # print name,i,result
            if 'disabled' in name or 'broken' in name or 'acquire' not in name:
                self.assertFalse(result, name)
            else:
                self.assertTrue(result, name)


class TestDecenter(guiderTester.GuiderTester, unittest.TestCase):
    def test_set_decenter_enable(self):
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        self.assertTrue(self.gState.decenter)

    def test_set_decenter_disable(self):
        decenters = {'decenterRA': 1, 'decenterDec': 2, 'mangaDither': 'N'}
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        masterThread.set_decenter(self.cmd, {}, self.gState, False)
        self.assertFalse(self.gState.decenter)
        self.assertEqual(self.gState.decenterRA, 0)
        self.assertEqual(self.gState.decenterDec, 0)
        self.assertEqual(self.gState.mangaDither, 'C')

    def _set_decenter_ok(self, decenters):
        masterThread.set_decenter(self.cmd, {}, self.gState, True)
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        self.assertEqual(self.gState.decenterRA, decenters.get('decenterRA', 0))
        self.assertEqual(self.gState.decenterDec, decenters.get('decenterDec', 0))
        self.assertEqual(self.gState.mangaDither, decenters.get('mangaDither', '?'))
        self.assertIn(self.cmd, self.gState.decenterCmd)

    def test_set_decenter_new(self):
        decenters = {'decenterRA': 1, 'decenterDec': 2}
        self._set_decenter_ok(decenters)

    def test_set_decenter_mangaDither(self):
        decenters = {'decenterRA': 1, 'decenterDec': 2, 'mangaDither': 'N'}
        self._set_decenter_ok(decenters)

    def test_set_decenter_cannot_change(self):
        decenters = {'decenterRA': 1, 'decenterDec': 2}
        masterThread.set_decenter(self.cmd, decenters, self.gState, None)
        self.assertTrue(self.cmd.didFail)
        self.assertFalse(self.gState.decenter)
        self.assertEqual(self.gState.decenterRA, 0)
        self.assertEqual(self.gState.decenterDec, 0)
        self.assertEqual(self.gState.mangaDither, 'C')


class TestSetRefraction(guiderTester.GuiderTester, unittest.TestCase):
    def _set_refraction(self, corrRatio, plateType, surveyMode, expect):
        self.gState.refractionBalance = -100  # ensure it's always different to start.
        masterThread.set_refraction(self.cmd, self.gState, corrRatio, plateType, surveyMode)
        self.assertEqual(self.gState.refractionBalance, expect)

    def test_corrRatio_1(self):
        self._set_refraction(1, None, None, 1)

    def test_apogee(self):
        self._set_refraction(None, 'APOGEE', None, 1)

    def test_boss(self):
        self._set_refraction(None, 'eBOSS', None, 0)


class TestStartStopGuider(guiderTester.GuiderTester, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # can only configure the dispatcher once.
        if Model.dispatcher is None:
            Model.setDispatcher(KeyVarDispatcher())
        Actor.setupRootLogger = TestHelper.setupRootLogger

    def _close_port(self):
        """
        Close the connection: requires handling the deferred. Details here:
            https://jml.io/pages/how-to-disconnect-in-twisted-really.html
        """
        deferred = self.actor.commandSources.port.stopListening()
        deferred.callback(None)

    def tearDown(self):
        # have to clear any actors that were registered previously.
        Model._registeredActors = set()
        super(TestStartStopGuider, self).tearDown()

    def setUp(self):
        self.actor = guiderActor.GuiderActor.GuiderActor.newActor(
            location='apo', makeCmdrConnection=False)
        self.actor.runInReactorThread = False
        super(TestStartStopGuider, self).setUp()
        self.addCleanup(self._close_port)
        # Do this after super setUp, as that's what creates actorState.
        self.queues = {}
        self.queues[guiderActor.GCAMERA] = Queue('gcamera')
        self.queues[guiderActor.MASTER] = Queue('master')
        myGlobals.actorState.queues = self.queues

    def _start_guider(self, nWarn=0, **kwargs):
        guiderTester.updateModel('mcp', TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc', TestHelper.tccState['tracking'])
        self.gState.cartridge = self.actorState.models['guider'].keyVarDict['cartridgeLoaded'][0]

        masterThread.start_guider(
            self.cmd,
            self.gState,
            self.actorState,
            self.actorState.queues,
            **kwargs)

        msg = self.queues[guiderActor.MASTER].get(False)
        self.assertEqual(msg.type, guiderActor.Msg.STATUS)
        self.assertFalse(msg.finish)

        self.assertEqual(self.gState.cmd, self.cmd)

        startFrame = self.actorState.models['gcamera'].keyVarDict['nextSeqno'][0] + 1
        self.assertEqual(self.gState.startFrame, startFrame)

        msg = self.queues[guiderActor.GCAMERA].get(False)
        self.assertEqual(msg.type, guiderActor.Msg.EXPOSE)
        self.assertEqual(msg.camera, kwargs.get('camera', 'gcamera'))
        self.assertEqual(msg.expTime, kwargs.get('expTime', 5))
        self.assertEqual(msg.stack, kwargs.get('stack', 1))
        self._check_cmd(0, 1, nWarn, 0, False)

    def test_start_guider(self):
        self._start_guider()

    def test_start_guider_stack_2(self):
        self._start_guider(stack=2)

    def test_start_guider_expTime_10(self):
        self._start_guider(expTime=10)

    def test_start_guider_already_running_force(self):
        self.gState.cmd = TestHelper.Cmd()
        self._start_guider(force=True, nWarn=1)

    def test_start_guider_ecamera(self):
        self._start_guider(camera='ecamera')

    def test_start_guider_no_plate(self):
        self.gState.cartridge = -1
        masterThread.start_guider(self.cmd, self.gState, self.actorState, self.actorState.queues)
        self.assertTrue(self.queues[guiderActor.GCAMERA].empty())
        self.assertTrue(self.queues[guiderActor.MASTER].empty())
        self.assertIsNone(self.gState.cmd)
        self._check_cmd(0, 0, 0, 0, True, True)

    def test_start_guider_not_ok_to_guide(self):
        guiderTester.updateModel('mcp', TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc', TestHelper.tccState['stopped'])
        self.gState.cartridge = self.actorState.models['guider'].keyVarDict['cartridgeLoaded'][0]
        masterThread.start_guider(self.cmd, self.gState, self.actorState, self.actorState.queues)
        self.assertTrue(self.queues[guiderActor.GCAMERA].empty())
        self.assertTrue(self.queues[guiderActor.MASTER].empty())
        self.assertIsNone(self.gState.cmd)
        self._check_cmd(0, 0, 1, 0, True, True)

    def test_start_guider_already_running(self):
        oldCmd = TestHelper.Cmd()
        self.gState.cmd = oldCmd
        masterThread.start_guider(self.cmd, self.gState, self.actorState, self.actorState.queues)
        self.assertTrue(self.queues[guiderActor.GCAMERA].empty())
        self.assertTrue(self.queues[guiderActor.MASTER].empty())
        self.assertEqual(self.gState.cmd, oldCmd)
        self._check_cmd(0, 0, 0, 0, True, True)

    def _stop_guider(self, success=True):
        self.gState.cmd = self.cmd
        masterThread.stop_guider(
            self.cmd,
            self.gState,
            self.actorState,
            self.actorState.queues,
            1234,
            success)
        msg = self.queues[guiderActor.GCAMERA].get(False)
        self.assertEqual(msg.type, guiderActor.Msg.ABORT_EXPOSURE)
        self.assertEqual(msg.quiet, True)
        self.assertIsNone(self.gState.cmd)

    def test_stop_guider(self):
        self._stop_guider()
        self._check_cmd(0, 3, 0, 0, True, didFail=False)

    def test_stop_guider_with_movie(self):
        self.gState.startFrame = 10
        self._stop_guider()
        self._check_cmd(0, 3, 0, 0, True, didFail=False)
        self.assertIsNone(self.gState.startFrame)

    def test_stop_guider_failed(self):
        self._stop_guider(success=False)
        self._check_cmd(0, 0, 0, 0, True, didFail=True)

    def test_stop_guider_already_off(self):
        masterThread.stop_guider(
            self.cmd,
            self.gState,
            self.actorState,
            self.actorState.queues,
            1234,
            True)
        self._check_cmd(0, 0, 0, 0, True, didFail=True)


class TestSetTime(guiderTester.GuiderTester, unittest.TestCase):
    def test_set_time(self):
        expTime = 5
        stack = 1
        readTime = 1
        masterThread.set_time(self.gState, expTime, stack=stack, readTime=readTime)
        self.assertEqual(self.gState.expTime, expTime)
        self.assertEqual(self.gState.stack, stack)
        self.assertEqual(self.gState.readTime, readTime)
        for pid in self.gState.pid.values():
            self.assertEqual(pid.dt, (expTime + readTime) * stack + 5)


if __name__ == '__main__':
    verbosity = 2

    suite = None
    # to test just one piece
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestStartStopGuider)
    # suite = unittest.TestLoader().loadTestsFromTestCase(TestSetTime)
    if suite:
        unittest.TextTestRunner(verbosity=verbosity).run(suite)
    else:
        unittest.main(verbosity=verbosity)
