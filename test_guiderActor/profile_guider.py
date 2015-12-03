#!/usr/bin/env python
"""Test the speed of guider image processing."""
import cProfile
import pstats
import os,glob

from actorcore import TestHelper

from guiderActor import masterThread
from guiderActor.gimg import guiderImage
import guiderActor.myGlobals as myGlobals

import guiderTester

def updateModel(name,model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name,model)

darkFile = 'gcam/gimg-0001.fits.gz'
flatFile = 'gcam/gimg-0003.fits.gz'
dataFile = 'gcam/gimg-0040.fits.gz'

helper = guiderTester.GuiderTester()
helper.setUp()
# force the gi instance to be at location=APO.
helper.gi = guiderImage.GuiderImageAnalysis(helper.setPoint_good,'APO')
# setup to use all 16 "real" probes
helper._init_probes_gimg_0040()
updateModel('mcp',TestHelper.mcpState['boss_science'])
helper.actorState.bypassDark = False

# pre-process the dark and flat.
helper.gi.analyzeDark(darkFile,cmd=helper.cmd)
helper.gi.analyzeFlat(flatFile,helper.gState.gprobes,cmd=helper.cmd)
helper.gState.cmd = helper.cmd

def timeit_helper():
    masterThread.guideStep(helper.actor,None,helper.cmd,helper.gState,dataFile,False,helper.gi)
    os.remove('gcam/proc-gimg-0040.fits.gz')

# Paste this line for a multi-pass timing run in ipython, after pasting the above.
# %timeit timeit_helper()

def cleanup():
    for filename in glob.glob('?cam/proc-*.fits.gz'):
        os.remove(filename)

prof = cProfile.Profile()
#result = prof.runcall(helper.gi,helper.cmd,dataFile,helper.gState.gprobes,-40)
result = prof.runcall(masterThread.guideStep,helper.actor,None,helper.cmd,helper.gState,dataFile,False,helper.gi)
cleanup()
prof.dump_stats('guideStep.profile')
stats = pstats.Stats('guideStep.profile')
stats.strip_dirs()
stats.sort_stats('cumtime').print_stats(25)
