#!/usr/bin/env python
"""unittests for guiderActor itself"""

import unittest

from actorcore import Actor
from opscore.actor import Model, KeyVarDispatcher
from actorcore import TestHelper

from guiderActor import GuiderActor

logDirBase = 'temp/'

# NOTE: for this particular set of unittests, you need unittest.TestCase as the
# first parent, in order for setUp to work correctly and allow TestHelper.logBuffer
# to do what it's supposed to. I haven't tracked down exactly why, but it's a
# simple change.
class TestGuiderActor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # can only configure the dispatcher once.
        Model.setDispatcher(KeyVarDispatcher())
        Actor.setupRootLogger = TestHelper.setupRootLogger

    def setUp(self):
        super(TestGuiderActor,self).setUp()
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
        self.assertIsInstance(self.guider,GuiderActor.GuiderActorAPO)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        print logged
        self.assertIsInstance(self.guider.actorState,Actor.ActorState)
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_APO',logged)

    def test_init_lco(self):
        self.guider = GuiderActor.GuiderActor.newActor(location='lco',makeCmdrConnection=False)
        self.assertIsInstance(self.guider,GuiderActor.GuiderActorLCO)
        self.assertIsInstance(self.guider.actorState,Actor.ActorState)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_LCO',logged)

if __name__ == '__main__':
    verbosity = 2
    
    unittest.main(verbosity=verbosity)

