#!/usr/bin/env python
"""unittests for guiderActor itself"""

import unittest

from actorcore import Actor
from opscore.actor import Model, KeyVarDispatcher
from actorcore import TestHelper

import guiderTester

from guiderActor import GuiderActor

logDirBase = 'temp/'

# NOTE: for this particular set of unittests, you need unittest.TestCase as the
# first parent, in order for setUp to work correctly and allow TestHelper.logBuffer
# to do what it's supposed to. I haven't tracked down exactly why, but it's a
# simple change.
class GuiderActorTester(unittest.TestCase):
    """Parent class for tests that actually need a proper guiderActor instance."""
    @classmethod
    def setUpClass(cls):
        # can only configure the dispatcher once.
        if Model.dispatcher is None:
            Model.setDispatcher(KeyVarDispatcher())
        Actor.setupRootLogger = TestHelper.setupRootLogger

    def setUp(self):
        super(GuiderActorTester,self).setUp()
        self.addCleanup(self._close_port)

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
        super(GuiderActorTester,self).tearDown()


class TestGuiderActor(GuiderActorTester):
    def test_init_apo(self):
        self.actor = GuiderActor.GuiderActor.newActor(location='apo',makeCmdrConnection=False)
        self.assertIsInstance(self.actor,GuiderActor.GuiderActorAPO)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        print logged
        self.assertIsInstance(self.actor.actorState,Actor.ActorState)
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_APO',logged)

    def test_init_lco(self):
        self.actor = GuiderActor.GuiderActor.newActor(location='lco',makeCmdrConnection=False)
        self.assertIsInstance(self.actor,GuiderActor.GuiderActorLCO)
        self.assertIsInstance(self.actor.actorState,Actor.ActorState)
        self.assertEqual(TestHelper.logBuffer.basedir,'/data/logs/actors/guider')
        logged = TestHelper.logBuffer.getvalue()
        self.assertIn('attaching command set GuiderCmd',logged)
        # self.assertIn('attaching command set CameraCmd_LCO',logged)


class TestGuidingIsOK_APO(GuiderActorTester,guiderTester.GuiderTester):
    def setUp(self):
        self.actor = GuiderActor.GuiderActor.newActor(location='apo',makeCmdrConnection=False)
        super(TestGuidingIsOK_APO,self).setUp()
        # need to explicitly call this, to get cmd,actorState,etc. set up.
        guiderTester.GuiderTester.setUp(self)

    def _guidingIsOK(self, expect, nWarn=0, force=False):
        result = self.actor.guidingIsOK(self.cmd, self.actorState, force=force)
        self.assertEqual(result, expect)
        self._check_levels(0, 0, nWarn, 0)

    def test_force(self):
        self._guidingIsOK(True, force=True)

    def test_boss_science(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(True)

    def test_ffs_closed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['all_off'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_ffs_closed_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['all_off'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs'])
        self._guidingIsOK(True, 1)

    def test_fflamp_on(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['flats'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_fflamp_on_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['flats'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs','lamp_ff'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(True, 1)
    def test_arclamps_on(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['arcs'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self._guidingIsOK(False, 1)
    def test_arclamps_on_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['arcs'])
        guiderTester.updateModel('tcc',TestHelper.tccState['tracking'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['ffs','lamp_hgcd','lamp_ne'])
        self._guidingIsOK(True, 1)

    def test_tcc_halted(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['halted'])
        self._guidingIsOK(False, 1)
    def test_tcc_motion_bypassed(self):
        guiderTester.updateModel('mcp',TestHelper.mcpState['boss_science'])
        guiderTester.updateModel('tcc',TestHelper.tccState['halted'])
        self.actorState.models['sop'].keyVarDict['bypassedNames'].set(['axes'])
        self._guidingIsOK(True, 1)



if __name__ == '__main__':
    verbosity = 2
    
    unittest.main(verbosity=verbosity)

