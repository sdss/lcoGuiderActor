#!/usr/bin/env python
"""unittests for guiderActor itself"""

import unittest

from actorcore import Actor, ICC
from opscore.actor import Model, KeyVarDispatcher

from actorcore import TestHelper

from guiderActor import GuiderActor

import guiderTester

logDirBase = 'temp/'

class TestGuiderActor(guiderTester.GuiderTester,unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # can only configure the dispatcher once.
        Model.setDispatcher(KeyVarDispatcher())
        Actor.setupRootLogger = TestHelper.setupRootLogger
        ICC.makeOpsFileLogger = TestHelper.fakeOpsFileLogger

    def setUp(self):
        # have to clear any actors that were registered previously.
        Model._registeredActors = set()

    def tearDown(self):
        # close the connection: requires handling the deferred.
        # see here for more details:
        #     https://jml.io/pages/how-to-disconnect-in-twisted-really.html
        if getattr(self,'guider',None) is not None:
            deferred = self.guider.commandSources.port.stopListening()
            deferred.callback(None)

    def test_init_apo(self):
        self.guider = GuiderActor.GuiderActor.newActor(location='apo',makeCmdrConnection=False)
        self.assertIsInstance(self.guider,GuiderActor.GuiderAPO)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        # self.assertIsInstance(self.guider.actorState,GuiderActor.State)
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_APO',logged)

    def test_init_lco(self):
        self.guider = GuiderActor.GuiderActor.newActor(location='lco',makeCmdrConnection=False)
        self.assertIsInstance(self.guider,GuiderActor.GuiderLCO)
        # self.assertIsInstance(self.guider.actorState,GuiderActor.State)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_LCO',logged)

if __name__ == '__main__':
    verbosity = 2
    
    unittest.main(verbosity=verbosity)

