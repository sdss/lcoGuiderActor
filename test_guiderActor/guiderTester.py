"""
Helps with setting up and tearing down guideractor tests.

Example:
    class TestGuiderThing(guiderTester.GuiderTester,unittest.TestCase):
        def setUp(self):
            #do some stuff
            super(TestGuiderImage,self).setUp())

        def test_guiderImageThingy(self):
            self.call_gi(filename)

    if __name__ == '__main__':
        unittest.main()
"""

import os
import unittest

from actorcore import TestHelper

from guiderActor.gimg import guiderImage
from guiderActor import GuiderState
import guiderActor.myGlobals as myGlobals
from guiderActor import GuiderActor
import json


guiderImagesPath = os.path.realpath(os.environ['GUIDERIMAGES_DIR'])


def updateModel(name, model):
    """Update the named actorState model with new parameters."""
    myGlobals.actorState.models[name] = TestHelper.Model(name, model)


def getGProbeTesterData(filename=None, mjd=None, plateid=None,
                        fscan_mjd=None, fscan_id=None, camera='gcam'):
    """Returns the GuiderTester gprobe configuration as a dictionary."""

    assert (filename or (plateid and fscan_mjd and fscan_id) or
            (not plateid and not fscan_mjd and not fscan_id)), \
        ('cannot find gprobe data for plateid={0}, fscan_mjd={1}, '
         'fscan_id={2}'.format(plateid, fscan_mjd, fscan_id))

    assert camera in ['gcam', 'ecam']

    if not filename:
        if not plateid:
            filename = os.path.join(guiderImagesPath,
                                    'gprobes_default.json')
        else:
            filename = os.path.join(
                guiderImagesPath,
                '{0}/{1}/gprobes_{2}-{3}-{4}.json'.format(camera, mjd, plateid,
                                                          fscan_mjd, fscan_id))

    return json.load(open(filename))


# Decorator to skip a test if the guider images product is not set up.
skipIfNoGuiderImages = unittest.skipIf(
    not os.path.exists(guiderImagesPath), 'guiderImages product not set up.')


def getTestFile(directory, mjd, filename, raiseError=True):
    """Returns the full path in the guiderImages product for a test file."""

    filePath = os.path.join(guiderImagesPath,
                            directory, str(mjd), filename)
    if raiseError and not os.path.exists(filePath):
        raise ValueError('filename {0} does not exist in guiderImages'
                         .format(filePath))

    return filePath


class GuiderTester(TestHelper.ActorTester):
    """Test suites should subclass this and unittest, in that order."""

    def setUp(self):
        """Populate fake guide probes, etc."""

        self.verbose = True
        self.name = 'guider'
        super(GuiderTester, self).setUp()

        myGlobals.actorState = self.actorState

        self.setPoint_good = -40
        self.setPoint_bad = -35

        if self.actorState.actor.location is None:
            self.actorState.actor.location = 'APO'

        self.gi = guiderImage.GuiderImageAnalysis(
            self.setPoint_good, self.actorState.actor.location)

        gState = GuiderState.GuiderState()
        self.gState = gState
        myGlobals.actorState.gState = self.gState

        GuiderActor.set_default_pids(self.actor.config, self.gState)
        GuiderActor.set_pid_scaling(self.actor.config, self.gState)
        GuiderActor.set_telescope(self.actor.config, self.gState)
        GuiderActor.set_gcamera(self.actor.config, self.gState)

        # Initialises default gprobe and guideInfo keys
        guiderTesterData = getGProbeTesterData()
        gprobeKey = guiderTesterData['gprobeKey']
        guideInfoKey = guiderTesterData['guideInfoKey']

        for name in gprobeKey:
            gk = gprobeKey[name]
            gik = guideInfoKey[name]
            gState.gprobes[gk[1]] = GuiderState.GProbe(gk[1], gprobeKey=gk,
                                                       guideInfoKey=gik)
            if 'disabled' in name:
                gState.gprobes[gk[1]].disabled = True

    def init_probes(self, mjd, plateid, fscan_mjd, fscan_id, camera='gcam'):
        """Initialize the "real" guide probes.

        Use with tests that e.g. call guiderImageAnalysis
        on the MJD 57356 files.

        """

        guiderTesterData = getGProbeTesterData(
            mjd=mjd, plateid=plateid, fscan_mjd=fscan_mjd,
            fscan_id=fscan_id, camera=camera)

        platedb_gprobe = guiderTesterData['platedb_gprobe']
        platedb_guideInfo = guiderTesterData['platedb_guideInfo']
        platedb_gprobesInUse = guiderTesterData['platedb_gprobesInUse']

        for probe, info, bits in zip(platedb_gprobe, platedb_guideInfo,
                                     platedb_gprobesInUse):

            __, bits = bits.strip('()').split('=')

            self.gState.gprobes[probe[1]] = GuiderState.GProbe(
                probe[1], gprobeKey=probe, guideInfoKey=info,
                gprobeBits=int(bits, 16))

        tritium = platedb_gprobe[16]
        self.gState.gprobes[tritium[1]] = GuiderState.GProbe(
            tritium[1], gprobeKey=tritium, gprobeBits=2)

    def _call_gi(self, filename, setpoint=None, args=[]):
        """Use this to simplify calls to guiderImageAnalysis."""

        if setpoint is None:
            setpoint = self.setPoint_good

        return self.gi(self.cmd, filename, self.gState.gprobes,
                       setpoint, *args)

    def _remove_file(self, filename):
        if os.path.exists(filename):
            os.remove(filename)


class GuiderThreadTester(GuiderTester, unittest.TestCase):
    """Test suites should subclass this and unittest, in that order."""

    def __init__(self, *args, **kwargs):
        """Load up the cmd calls for this test class, so _check_cmd will use them."""
        unittest.TestCase.__init__(self, *args, **kwargs)
        # -1 is the test function, -2 is test class, -3 (or 0) should be main
        class_name = self.id().split('.')[-2]
        self._load_cmd_calls(class_name)
        # lets us see really long list/list diffs
        self.maxDiff = None
