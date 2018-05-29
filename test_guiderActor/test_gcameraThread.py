#!/usr/bin/env python
"""
Test the behavior of the various guider gcamera thread commands.
"""
import unittest
from Queue import Queue

import guiderActor
import guiderTester
from guiderActor import gcameraThread


class TestExpose(guiderTester.GuiderThreadTester, unittest.TestCase):
    def setUp(self):
        super(TestExpose, self).setUp()
        self.replyQueue = Queue('fakeReply')
        self.fail_on_no_cmd_calls = True  # we need cmd_calls for all of these.

    def _expose(self, expTime, didFail=False, expect={}, **kwargs):
        gcameraThread.expose(self.cmd, self.actorState, self.replyQueue, expTime, **kwargs)
        msg = self._queue_get(self.replyQueue)
        camera = kwargs.get('camera', 'gcamera')
        self.assertEqual(msg.type, expect.get('expType', guiderActor.Msg.EXPOSURE_FINISHED))
        self.assertEqual(msg.filename, self.actorState.models[camera].keyVarDict['filename'][0])
        self.assertEqual(msg.camera, camera)
        self._check_cmd(1, 1, 0, 0, False, didFail)

    def test_expose_10(self):
        self._expose(10)

    def test_expose_10_stack_2(self):
        self._expose(10, stack=2)

    def test_expose_ecamera_10(self):
        self._expose(10, camera='ecamera')

    def test_flat_5(self):
        self._expose(5, expType='flat', cartridge=10, expect={
                     'expType': guiderActor.Msg.FLAT_FINISHED})

    def test_dark_10_stack_9(self):
        self._expose(10, stack=9, expType='dark', expect={
                     'expType': guiderActor.Msg.DARK_FINISHED})

    def test_expose_10_fails(self):
        self.cmd.failOn = 'gcamera expose time=10 stack=1'
        expTime = 10
        gcameraThread.expose(self.cmd, self.actorState, self.replyQueue, expTime)
        msg = self._queue_get(self.replyQueue)
        self.assertEqual(msg.type, guiderActor.Msg.EXPOSURE_FINISHED)
        self._check_cmd(1, 1, 1, 0, False, True)


if __name__ == '__main__':
    verbosity = 2

    unittest.main(verbosity=verbosity)
