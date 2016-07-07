#!/usr/bin/env python
"""Test the speed of guider image processing."""
import cProfile
import pstats
import os
import glob

from actorcore import TestHelper

from guiderActor import masterThread
from guiderActor.gimg import guiderImage
import guiderActor.myGlobals as myGlobals

import guiderTester


def updateModel(name, model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name, model)

mjd = 57357

darkFile = guiderTester.getTestFile('gcam/', mjd, 'gimg-0001.fits.gz')
flatFile = guiderTester.getTestFile('gcam/', mjd, 'gimg-0003.fits.gz')
dataFile = guiderTester.getTestFile('gcam/', mjd, 'gimg-0040.fits.gz')

gTester = guiderTester.GuiderTester()
gTester.attachCmdSets = False
gTester.setUp()
# force the gi instance to be at location=APO.
gTester.gi = guiderImage.GuiderImageAnalysis(gTester.setPoint_good, 'APO')
# setup to use all 16 "real" probes
gTester.init_probes(mjd, 7660, 57356, 1, camera='gcam')
updateModel('mcp', TestHelper.mcpState['boss_science'])
gTester.actorState.bypassDark = False

# pre-process the dark and flat.
gTester.gi.analyzeDark(darkFile, cmd=gTester.cmd)
gTester.gi.analyzeFlat(flatFile, gTester.gState.gprobes, cmd=gTester.cmd)
gTester.gState.cmd = gTester.cmd


def timeit_gTester():
    masterThread.guideStep(gTester.actor, None, gTester.cmd,
                           gTester.gState, dataFile, False, gTester.gi)
    os.remove('gcam/proc-gimg-0040.fits.gz')

# Paste this line for a multi-pass timing run in ipython, after pasting the above.
# %timeit timeit_gTester()


def cleanup():
    for filename in glob.glob('?cam/proc-*.fits.gz'):
        os.remove(filename)

prof = cProfile.Profile()
# result = prof.runcall(gTester.gi,gTester.cmd,dataFile,gTester.gState.gprobes,-40)
result = prof.runcall(masterThread.guideStep, gTester.actor, None, gTester.cmd,
                      gTester.gState, dataFile, False, gTester.gi)
cleanup()
prof.dump_stats('guideStep.profile')
stats = pstats.Stats('guideStep.profile')
stats.strip_dirs()
stats.sort_stats('cumtime').print_stats(25)
