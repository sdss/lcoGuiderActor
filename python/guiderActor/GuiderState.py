"""
Classes related to the current state of the guider.
"""

import numpy
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
    def __init__(self,id=-9999,gprobeKey=None,guideInfo=None):
        """Pass the contents of the platedb.gprobe and/or platedb.guideInfo keyword to initialize"""
        self.id = id
        self._bits = GOOD
        self._ugriz = [numpy.nan,]*5
        self.ref_mag = numpy.nan
        if gprobeKey is not None:
            self.from_platedb_gprobe(gprobeKey)
        if guideInfo is not None:
            self.from_platedb_guideInfo(guideInfoKey)
    
    def checkFocus(self):
        """Set the above/below focus bits based on the focusOffset value."""
        # allow a small range of allowed focus offsets.
        if self.focusOffset < -50:
            self.aboveFocus = True
        elif self.focusOffset > 50:
            self.belowFocus = True
        else:
            self.aboveFocus = False
            self.belowFocus = False
        
    def checkTritium(self):
        """If this probe is labeled a tritium star, disable it."""
        if self.fiber_type == 'TRITIUM':
            self.disabled =True

    def from_platedb_gprobe(self,gprobeKey):
        """
        Fill in data from the platedb.gprobe key.
        Expects a list, as output by CmdVar.getKeyVarData()
        """
        self._check_id(gprobeKey[1],'platedb.gprobe')
        self.broken = False if gprobeKey[2] else True
        if self.broken:
            self.disabled = True
        self.xCenter = gprobeKey[3]
        self.yCenter = gprobeKey[4]
        self.radius = gprobeKey[5]
        self.rotation = gprobeKey[6]
        self.xFerruleOffset = gprobeKey[7]
        self.yFerruleOffset = gprobeKey[8]
        self.focusOffset = gprobeKey[9]
        self.fiber_type = gprobeKey[10]
        self.checkTritium()
        self.rotStar2Sky = numpy.nan
        self.haOffsetTimes = {}
        self.haXOffsets = {}
        self.haYOffsets = {}
        self.checkFocus()
    #...
    
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
    #...

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
        return False if ((self.broken) or (self.noStar)) else True

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
        zd = 90. - actorState.models["tcc"].keyVarDict["axePos"][1]
        # TBD: zd=0 never occurs for tracking, but need to test for zd=0 for simulate
        airmass = 1./math.cos(math.radians(zd))
        gobs = value[1] + airmass*k0_g
        robs = value[2] + airmass*k0_r
        self.ref_mag = robs + a1 + a2*(gobs-robs) + a3*(gobs-robs)**2
        #self.ref_mag = (gobs+robs)/2. #jkp TBD: placeholder
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
        self.expTime = 0
        self.stack = 1
        self.inMotion = False
        self.centerUp = False
        self.guideCmd = None

        self.fscanMJD = self.fscanID = -1
        self.design_ha = numpy.nan

        self.plugPlateScale = numpy.nan
        self.dSecondary_dmm = numpy.nan
        self.gcameraPixelSize = numpy.nan
        self.gcameraMagnification = numpy.nan
        
        # Start with all fibers 
        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")
        self.setRefractionBalance(0.0)
        
        # Will contain [id]:gProbe pairs
        self.gprobes = {}
        
        self.pid = {}               # PIDs for various axes
        for what in ["raDec", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

        # reset the decenter positions.
        self.decenter = False                #gstate only
        self.setDecenter("decenterRA",0.)
        self.setDecenter("decenterDec",0.)
        self.setDecenter("decenterRot",0.)
        self.decenterChanged = True
        self.decenterFocus = numpy.nan
        self.decenterScale = numpy.nan
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

    def setRefractionBalance(self, value=0.0):
        self.refractionBalance = value
        
    def setCmd(self, cmd=None):
        self.guideCmd = cmd

    def setDecenter(self, axis, value=0):
        """Set axis="decenter[RA,Dec,Rot]" to value."""
        if axis == "decenterRA":
            self.decenterRA = value
        elif axis == "decenterDec":
            self.decenterDec = value
        elif axis == "decenterRot":
            self.decenterRot = value
        else:
            raise RuntimeError, ("Unknown decenter axis name %s" % axis)

    def setScales(self, plugPlateScale=None,
                  dSecondary_dmm=None,
                  gcameraPixelSize=None,
                  gcameraMagnification=None):

        if plugPlateScale != None:
            self.plugPlateScale = plugPlateScale
        if dSecondary_dmm != None:
            self.dSecondary_dmm = dSecondary_dmm
        if gcameraPixelSize != None:
            self.gcameraPixelSize = gcameraPixelSize
        if gcameraMagnification != None:
            self.gcameraMagnification = gcameraMagnification

