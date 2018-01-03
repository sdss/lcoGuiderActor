"""
Classes related to the current state of the guider.
"""

import numpy as np
import math

import PID
from guiderActor import myGlobals

# gprobebits
# To help manage the guide probe status bits.
# WARNING: ensure the usage in guiderActor, platedbActor, STUI are consistent with this!
GOOD   =  0x00     # A happy, working, in-use probe has all bits set to 0.
BROKEN =  0x01     # a known broken probe, labeled as such in plPlugMap
NOSTAR =  0x02     # probe with no star in plPlugMap (e.g. tritium)
DISABLED = 0x04     # probe that has been disabled by the observers (e.g. no star present, double star, wrong position)
ABOVEFOCUS = 0x08  # fiber is above the focal plane
BELOWFOCUS = 0x10  # fiber is below the focal plane
TOOFAINT = 0x20    # observed star is too faint to be reliably used for guiding
UNKNOWN = 0xff     # shouldn't ever happen

class GProbe(object):
    """
    Contains information about a single guide probe.

    ugriz is an array of fiber magnitudes  (through 2arcsec fibers) from in plPlugMapP.par.
    When set, it computes self.ref_mag, which is the synthetic predicted magnitude for this fiber.

    GProbe flag bits are set via the corresponding property:
        broken, disabled (enabled), noStar, notExist, atFocus (aboveFocus,belowFocus), toofaint
    and gProbe.good will tell you if all bits are in the OK state.
    """
    def __init__(self, id=-9999, gprobeKey=None, guideInfoKey=None, gprobeBits=GOOD):
        """
        Init this probe using the appropriate keywords.
        kwargs:
            gprobeKey: platedb.gprobe keyword
            guideInfoKey: platedb.guideInfo keyword
            gprobebits (int): platedb.gprobesInUse value for this probe
        """
        self.id = id
        self._bits = gprobeBits
        self._ugriz = [np.nan,]*5
        self.ref_mag = np.nan
        if gprobeKey is not None:
            self.from_platedb_gprobe(gprobeKey)
        if guideInfoKey is not None:
            self.from_platedb_guideInfo(guideInfoKey)

    def checkFocus(self):
        """Set the above/below focus bits based on the focusOffset value."""
        # allow a small range of allowed focus offsets.
        if self.focusOffset < -150:
            self.aboveFocus = True
        elif self.focusOffset > 150:
            self.belowFocus = True
        else:
            self.aboveFocus = False
            self.belowFocus = False

    def from_platedb_gprobe(self,gprobeKey):
        """
        Fill in data from the platedb.gprobe key.
        Expects a list, as output by CmdVar.getKeyVarData()
        """

        # Gets the additional grpobe rotation due to the CCD chip rotation.
        actorConfig = myGlobals.actorState.actorConfig
        rotationCCD = actorConfig.getfloat('gcamera', 'ccdRotation')

        self._check_id(gprobeKey[1],'platedb.gprobe')
        self.broken = False if gprobeKey[2] else True
        if self.broken:
            self.disabled = True
        self.xCenter = gprobeKey[3]
        self.yCenter = gprobeKey[4]
        self.radius = gprobeKey[5]
        self.rotation = gprobeKey[6] + rotationCCD
        self.xFerruleOffset = gprobeKey[7]
        self.yFerruleOffset = gprobeKey[8]
        self.focusOffset = gprobeKey[9]
        self.fiber_type = gprobeKey[10]
        self.haOffsetTimes = {}
        self.haXOffsets = {}
        self.haYOffsets = {}
        self.checkFocus()

    def from_platedb_guideInfo(self,guideInfoKey):
        """
        Fill in data from the platedb.guideInfo keyword.
        Expects a list, as output by CmdVar.getKeyVarData()
        """
        self._check_id(guideInfoKey[0],'platedb.guideInfo')
        self.ra = guideInfoKey[1]
        self.dec = guideInfoKey[2]
        self.xFocal = guideInfoKey[3]
        self.yFocal = guideInfoKey[4]
        self.phi = guideInfoKey[5]
        self.throughput = guideInfoKey[6]
        self.set_rotStar2Sky()

    def set_rotStar2Sky(self):
        """Compute rotStar2Sky from phi and rotation, if they're available."""
        # rotStar2Sky is the angle to rotate (x, y) on the camera to (ra, alt)
        # phi is the orientation of the alignment hole measured clockwise from N
        # rotation is the anticlockwise rotation from x on the camera to the pin
        try:
            self.rotStar2Sky = 90 + self.rotation - self.phi
        except:
            self.rotStar2Sky = np.nan

    def _check_id(self,id,fromName):
        """
        Verify that the id is correct for this probe,
        or set it if it hasn't been set yet.

        fromName is the name of the actorkey the new id came from.
        """
        if self.id == -9999:
            self.id = id
        elif id != self.id:
            raise ValueError("%s id does not match current probe id!"%fromName)
        else:
            # otherwise, everything's fine.
            pass
    #...

    def _unset(self,bit):
        """Set bit to 0."""
        self._bits = self._bits & (~bit)

    def _set(self,bit):
        """Set bit to 1."""
        self._bits = self._bits | bit

    @property
    def good(self):
        """True if all probe bits except [above|below]Focus are 0."""
        return True if not (self._bits & ~ABOVEFOCUS & ~BELOWFOCUS) else False

    @property
    def exists(self):
        """True if this probe is not broken and has a star defined."""
        return False if self.broken else True

    @property
    def disabled(self):
        """True if this probe is disabled (negation of enabled)."""
        return True if (self._bits & DISABLED) else False
    @disabled.setter
    def disabled(self,value):
        self._set(DISABLED) if value else self._unset(DISABLED)

    @property
    def enabled(self):
        """True if this probe is enabled (negation of disabled)."""
        return False if (self._bits & DISABLED) else True
    @enabled.setter
    def enabled(self,value):
        self._unset(DISABLED) if value else self._set(DISABLED)

    @property
    def broken(self):
        """True if this probe is broken."""
        return True if (self._bits & BROKEN) else False
    @broken.setter
    def broken(self,value):
        self._set(BROKEN) if value else self._unset(BROKEN)

    @property
    def noStar(self):
        """True if this probe has no star defined in the plugmap."""
        return True if (self._bits & NOSTAR) else False
    @noStar.setter
    def noStar(self,value):
        self._set(NOSTAR) if value else self._unset(NOSTAR)

    @property
    def atFocus(self):
        """True if this fiber is at the focal plane."""
        return not ((self._bits & ABOVEFOCUS) | (self._bits & BELOWFOCUS))

    @property
    def aboveFocus(self):
        """True if this fiber is above the focal plane."""
        return True if (self._bits & ABOVEFOCUS) else False
    @aboveFocus.setter
    def aboveFocus(self,value):
        if value:
            self._unset(BELOWFOCUS) # can't be both above and below focus!
        self._set(ABOVEFOCUS) if value else self._unset(ABOVEFOCUS)

    @property
    def belowFocus(self):
        """True if this fiber is below the focal plane."""
        return True if (self._bits & BELOWFOCUS) else False
    @belowFocus.setter
    def belowFocus(self,value):
        if value:
            self._unset(ABOVEFOCUS) # can't be both above and below focus!
        self._set(BELOWFOCUS) if value else self._unset(BELOWFOCUS)

    @property
    def tooFaint(self):
        """True if this star in this fiber is too faint to use for guiding."""
        return True if (self._bits & TOOFAINT) else False
    @tooFaint.setter
    def tooFaint(self,value):
        self._set(TOOFAINT) if value else self._unset(TOOFAINT)

    @property
    def gprobebits(self):
        """The bits for this probe's status (int)."""
        return self._bits
    @gprobebits.setter
    def gprobebits(self,value):
        if not isinstance(value,int):
            raise ValueError("gprobebits must be set as an integer!")
        else:
            self._bits = value

    @property
    def tritium(self):
        return self.fiber_type == 'TRITIUM'

    @property
    def ugriz(self):
        '''
        The 2" fiber magnitudes of the object in this fiber.
        Computes the synthetic predicted reference magnitude (self.ref_mag) when set.
        '''
        return self._ugriz
    @ugriz.setter
    def ugriz(self,value):
        """
        Compute the reference magnitude for this probe's target.

        The magnitude the guider should measure for this star/fiber at
        the current telescope position. Guider effective wavelength is
        5400A, so calculate guidermag from g and r. Then correct for
        atmospheric extinction.

        APO atmospheric extinction coeff at airmass=1 taken from table 3 of the
        ubercal paper, Padmanabhan et al. 2008 ApJ.
        with the k0 value for the filters, g:0.17, r:0.10

        Color terms for transformation from r + a1 + a2*(g-r) + a3*(g-r)**2 = ref_mag
        a1=0.535, a2=0.506, a3=-0.0312
        Coefficient are from Masayuki, Gunn&Strkyer stds observed through modeled guider & g,r passbands
        """
        self._ugriz = value
        actorState = myGlobals.actorState
        k0_g = 0.17
        k0_r = 0.10
        a1 = 0.535
        a2 = 0.506
        a3 = -0.0312

        #get airmass from tcc only gives alt = tcc.axePos[2]
        try:
            zd = 90. - actorState.models["tcc"].keyVarDict["axePos"][1]
        except TypeError:
            zd = 90. # if the tcc doesn't return a proper position.
        # TODO: zd=0 never occurs for tracking, but need to test for zd=0 for simulate
        airmass = 1./math.cos(math.radians(zd))
        gobs = value[1] + airmass*k0_g
        robs = value[2] + airmass*k0_r
        self.ref_mag = robs + a1 + a2*(gobs-robs) + a3*(gobs-robs)**2
        #self.ref_mag = (gobs+robs)/2. #jkp TODO: placeholder
#...

class GuiderState(object):
    """
    The current state of the guider.

    Contains information about the currently loaded cartridge, the target of
    each guide probe, custom parameters (decentering, stacking, exposure time, etc.).

    Does not know about any gcamera exposures: that data is stored in FrameInfo.
    """

    def __init__(self):
        self.cartridge = -1
        self.plate = -1
        self.pointing = "?"
        self.plateType = "?"
        self.surveyMode = "?"
        self.expTime = 0
        self.readTime = 0
        self.stack = 1
        self.inMotion = False
        self.centerUp = False
        self.cmd = None
        self.startFrame = None # guider frame number to start movie at.

        self.fscanMJD = self.fscanID = -1
        self.design_ha = np.nan

        self.plugPlateScale = np.nan
        self.dSecondary_dmm = np.nan
        self.gcameraPixelSize = np.nan
        self.gcameraMagnification = np.nan
        self.longitude = np.nan
        self.focalRatio = np.nan

        # Start with all fibers
        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")
        self.refractionBalance = 0

        # Will contain [id]:gProbe pairs
        self.gprobes = {}

        self.bigFiberRadius = np.nan

        # PIDs for various axes, and their default, on-initialization values.
        self.pid = {}
        self.pid_defaults = {}
        self.pid_time = {}  # times that each pid is updated, to track delta-t.
        for axis in ["raDec", "rot", "scale", "focus"]:
            self.pid[axis] = PID.PID(self.expTime, 0, 0, 0)
            self.pid_defaults[axis] = {}
            self.pid_time[axis] = None

        # The min and max altitudes to scale PID terms between for the given axes.
        self.axes_to_scale = []
        self.alt_min = 90.
        self.alt_max = 90.

        # reset the decenter positions.
        self.clearDecenter()
    #...

    def deleteAllGprobes(self):
        """Delete all fibers """
        self.gprobes = {}

    def setGprobeState(self, fiber, enable=True):
        """
        Set the enable/disable state of either one fiber, or all fibers of type fiber.

        fiber types: ACQUIRE GUIDE TRITIUM
        If an integer, must refer to a currently loaded probe.
        """
        if fiber in ("ACQUIRE", "GUIDE", "TRITIUM"):
            fiber_type = fiber
            for gp in self.gprobes.values():
                if gp.info.fiber_type == fiber_type:
                    gp.enabled = enable
        else:
            self.gprobes[fiber].enabled = enable

    def setGuideMode(self, mode, enabled=True):
        """Enable a guide mode, from "axes", "focus", or "scale"."""
        if mode == "axes":
            self.guideAxes = enabled
        elif mode == "focus":
            self.guideFocus = enabled
        elif mode == "scale":
            self.guideScale = enabled
        else:
            raise RuntimeError, ("Unknown guide mode %s" % mode)

    def setRefractionBalance(self, plateType, surveyMode):
        """Set the refraction balance based on the survey name."""
        if plateType == 'APOGEE' or plateType == 'APOGEE-2' or plateType == 'APOGEE-2S':
            self.refractionBalance = 1
        elif (plateType == 'APOGEE&MaNGA' or plateType == 'APOGEE-2&MaNGA') and \
             (surveyMode == 'APOGEE lead'):
            self.refractionBalance = 1
        else:
            self.refractionBalance = 0

    def setDecenter(self, decenters, cmd, enable):
        """
        Set decenter[RA,Dec,Rot] and mangaDither to the values in decenters.
        Save the cmd that sent it, if it is a different position from the current one.
        """
        if enable is not None:
            self.decenter = enable

        newRA = decenters.get('decenterRA',0)
        newDec = decenters.get('decenterDec',0)
        newRot = decenters.get('decenterRot',0)
        # Store the cmd, if we'll have to actually apply a new location,
        # or just finish if the decenter is identical to the current one.
        if self.decenterRA != newRA or \
           self.decenterDec != newDec or \
           self.decenterRot != newRot or \
           enable is not None:
            self.decenterCmd.append(cmd)
        else:
            cmd.finish()
            return
        self.decenterRA = newRA
        self.decenterDec = newDec
        self.decenterRot = newRot
        # if this isn't the 'C' position, unless mangaDither is specified
        # we don't know what the actual dither location is!
        if newRA == 0 and newDec == 0 and newRot == 0:
            default = 'C'
        else:
            default = '?'
        self.mangaDither = decenters.get('mangaDither',default)

    def finish_decenter(self):
        """Finish any pending decenter cmds."""
        for cmd in self.decenterCmd:
            cmd.finish('')
        self.decenterCmd = []

    def clearDecenter(self):
        """Clear all decenter information."""

        self.decenterCmd = []
        self.decenter = False
        self.mangaDither = 'C'
        self.decenterRA = 0
        self.decenterDec = 0
        self.decenterRot = 0
        self.decenterFocus = np.nan
        self.decenterScale = np.nan

    def set_pid_defaults(self, axis, **kwargs):
        """Set the default pid values for axis."""
        for key in kwargs:
            self.pid_defaults[axis][key] = kwargs[key]

    def reset_pid_terms(self, terms=None):
        """Reset all PID terms, or just those listed."""
        if terms == None:
            terms = self.pid.keys()
        for key in terms:
            self.pid[key].reset()

    def scale_pid_with_alt(self, alt):
        """
        Change the PID coefficients with altitude, so we are smoothly increasing
        Ti as altitude increases.
        Return True if we changed the PID values.
        """
        def twoPointForm(x0,y0,x1,y1,x):
            """Equation of a line from two points"""
            return (y1-y0)/(x1-x0) * (x-x0) + y0

        old_Ti = self.pid['raDec'].Ti

        # use .get() with a reasonable value incase set_pid_defaults wasn't called.
        Ti_low = self.pid_defaults['raDec'].get('Ti_min',200)
        Ti_high = self.pid_defaults['raDec'].get('Ti_max',200)
        if alt < self.alt_min:
            Ti = Ti_low
        elif alt > self.alt_max:
            Ti = Ti_high
        else:  # between low and high limits...
            Ti = twoPointForm(self.alt_min, Ti_low, self.alt_max, Ti_high, alt)
        if old_Ti != Ti:
            for axis in self.axes_to_scale:
                self.pid[axis].setPID(Ti=Ti)
            return True
        else:
            return False

    def update_pid_time(self, axis, newtime):
        """Update the time associated with axis, and output dt (new-old)."""
        try:
            dt = newtime - self.pid_time[axis]
        except TypeError:
            dt = None
        self.pid_time[axis] = newtime
        return dt

    def output_pid(self, cmd=None):
        if cmd is None:
            cmd = self.cmd
        for axis in self.pid.keys():
            cmd.respond("pid=%s,%g,%g,%g,%g,%d" % (axis,
                        self.pid[axis].Kp, self.pid[axis].Ti, self.pid[axis].Td,
                        self.pid[axis].Imax, self.pid[axis].nfilt))


class FrameInfo(object):
    """
    Holds data about the most recently read image frame.
    """
    def __init__(self,frameNo,arcsecPerMM,guideCameraScale,plugPlateScale):
        """Sets all parameters to NaN, so that they at least exist."""
        self.frameNo = frameNo

        self.dRA = np.nan
        self.dDec = np.nan
        self.dRot = np.nan
        self.dFocus = np.nan
        self.dScale = np.nan

        self.filtRA = np.nan
        self.filtDec = np.nan
        self.filtRot = np.nan
        self.filtFocus = np.nan
        self.filtScale = np.nan

        self.offsetRA = np.nan
        self.offsetDec = np.nan
        self.offsetRot = np.nan
        self.offsetFocus = np.nan
        self.offsetScale = np.nan

        self.guideCameraScale = guideCameraScale
        self.plugPlateScale = plugPlateScale
        self.arcsecPerMM = arcsecPerMM
        self.micronsPerArcsec = 1/3600.0*plugPlateScale*1e3 # convert arcsec to microns
        self.seeing = np.nan

        # conversion for a Gaussian, use this eveywhere but in ipGguide.c
        # conversion from sigma to FWHM for a JEG double Gaussian is done in ipGguide.c (sigmaToFWHMJEG = 2.468)
        self.sigmaToFWHM = 2.354
        #ADU, avoid guiding on noise spikes during acquisitions
        #should be in photons, based on RON, Dark residual, SKY
        self.minStarFlux = 500

        self.guideRMS = 0.
        self.nguideRMS = 0
        self.guideXRMS = 0.
        self.guideYRMS = 0.
        self.guideRaRMS = 0.
        self.guideDecRMS = 0.
        self.inFocusFwhm = []

        self.guideAzRMS = np.nan      #not implemented yet
        self.guideAltRMS = np.nan

        self.guideFitRMS = np.nan     #not implemented yet
        self.nguideFitRMS = np.nan
        self.nrejectFitRMS = np.nan

        self.decenterRA = np.nan
        self.decenterDec = np.nan
        self.decenterRot = np.nan
        self.decenterFocus = np.nan
        self.decenterScale = np.nan

        self.refractionBalance = np.nan
        self.wavelength = np.nan
        self.dHA = np.nan

        self.A = np.matrix(np.zeros(3*3).reshape([3,3]))
        self.b = np.matrix(np.zeros(3).reshape([3,1]))
        self.b3 = 0
        self.guideAxes = False
        self.guideFocus = False
        self.guideScale = False

    def setGuideMode(self,gState=None):
        """Save the guide* values from gState, or reset them to False."""
        self.guideAxes = gState.guideAxes
        self.guideFocus = gState.guideFocus
        self.guideScale = gState.guideScale

    def setDecenter(self,gState=None):
        """Fill the decenter values from gState, or reset them to 0."""
        if gState:
            self.decenterRA  = gState.decenterRA
            self.decenterDec = gState.decenterDec
            self.decenterRot = gState.decenterRot
            self.decenterFocus = gState.decenterFocus
            self.decenterScale = gState.decenterScale
        else:
            self.decenterRA = 0.0
            self.decenterDec = 0.0
            self.decenterRot = 0.0
            self.decenterFocus = 0.0
            self.decenterScale = 0.0
#...
