#!/usr/bin/env python
"""Test the speed of guider image processing."""
import cProfile
import os

from actorcore import TestHelper

from guiderActor import masterThread
from guiderActor.gimg import guiderImage
from guiderActor import GuiderState
import guiderActor.myGlobals as myGlobals

import guiderTester

def updateModel(name,model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name,model)

darkFile = 'data/gimg-0001.fits.gz'
flatFile = 'data/gimg-0003.fits.gz'
dataFile = 'data/gimg-0040.fits.gz'

helper = guiderTester.GuiderTester()
helper.setUp()
actor = TestHelper.FakeActor('guider','guiderActor')
updateModel('mcp',TestHelper.mcpState['science'])
helper.actorState.bypassDark = False

# pre-process the dark and flat.
helper.gi.analyzeDark(darkFile,cmd=helper.cmd)
helper.gi.analyzeFlat(flatFile,helper.gState.gprobes,cmd=helper.cmd)
helper.gState.guideCmd = helper.cmd

def timeit_helper():
    masterThread.guideStep(actor,None,helper.cmd,helper.gState,dataFile,False,helper.gi)
    os.remove('data/proc-gimg-0040.fits.gz')

# Paste this line for a multi-pass timing run in ipython, after pasting the above.
# %timeit timeit_helper()


prof = cProfile.Profile()
#result = prof.runcall(helper.gi,helper.cmd,dataFile,helper.gState.gprobes,-40)
result = prof.runcall(masterThread.guideStep,actor,None,helper.cmd,helper.gState,dataFile,False,helper.gi)
prof.dump_stats('guideStep.profile')
