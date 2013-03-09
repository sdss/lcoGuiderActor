"""
Classes related to the current state of the guider.
"""

import numpy as np
import math 

from guiderActor import *

# gprobebits
# To help manage the guide probe status bits.
# WARNING: ensure the usage in platedb and STUI are consistent with this!
GOOD   =  0x00     # A happy, working, in-use probe has all bits set to 0.
BROKEN =  0x01     # a known broken probe, labeled as such in plPlugMap
NOSTAR =  0x02     # probe with no star in plPlugMap (e.g. tritium)
DISABLE = 0x04     # probe that has been disabled by the observers (e.g. no star present, double star, wrong position)
ABOVEFOCUS = 0x08  # star in probe is out of focus, above focal plane
BELOWFOCUS = 0x10  # star in probe is out of focus, below focal plane
UNKNOWN = 0xff     # shouldn't ever happen

class ProbeInfo(object):
    """
    Contains information about a single guide probe.
    
    ugriz is an array of fiber magnitudes  (through 2arcsec fibers) from in plPlugMapP.par.
    
    Probe flag bits are set via the corresponding property:
        broken, disabled, noStar, notExist, outOfFocus
    and probeInfo.good will tell you if all bits are in the OK state.
    """
    def __init__(self,gprobeKey=None,guideInfo=None):
        """Pass the contents of the platedb.probe keyword to initialize"""
        self.id = -9999
        self._bits = GOOD
        if gprobeKey is not None:
            self.from_gprobeKey(gprobeKey)
        if guideInfo is not None:
            self.from_guideInfoKey(guideInfoKey)
    
    def checkTritium(self):
        """If this probe is labeled a tritium star, disable it."""
        if self.fiberType == 'TRITIUM':
            self.disable =True
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
    
    def from_platedb_gprobe(self,gprobeKey):
        """
        Fill in data from the platedb.gprobe key.
        Expects a list, as output by CmdVar.getKeyVarData()
        """
        self._check_id(gprobeKey[1],'platedb.gprobe')        
        self.broken = False if gprobeKey[2] else self.broken = True
        self.xCenter = gprobeKey[3]
        self.yCenter = gprobeKey[4]
        self.radius = gprobeKey[5]
        self.rotation = gprobeKey[6]
        self.xFerruleOffset = gprobeKey[7]
        self.yFerruleOffset = gprobeKey[8]
        self.focusOffset = gprobeKey[9]
        self.fiberType = gprobeKey[10]
        self.checkTritium()
        self.rotStar2Sky = numpy.nan
        self.haOffsetTimes = {}
        self.haXOffsets = {}
        self.haYOffsets = {}
        self.ugriz = np.nan
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

    def _unset(self,bit):
        """Set bit to 0."""
        self._bits = self._bits & (~bit)
    
    def _set(self,bit):
        """Set bit to 1."""
        self._bits = self._bits | bit
    
    @property
    def good(self):
        """True if all probe bits are OK."""
        return True if not self._bits else False

    @property
    def disabled(self):
        """True if this probe is disabled (negation of enabled)."""
        return (self._bits & DISABLE)
    @disabled.setter
    def disabled(self,value):
        self._set(DISABLE) if value else self._unset(DISABLE)

    @property
    def enabled(self):
        """True if this probe is enabled (negation of disabled)."""
        return not (self._bits & DISABLE)
    @enabled.setter
    def enabled(self,value):
        self._unset(DISABLE) if value else self._set(DISABLE)

    @property
    def broken(self):
        """True if this probe is broken."""
        return (self._bits & BROKEN)
    @broken.setter
    def broken(self,value):
        self._set(BROKEN) if value else self._unset(BROKEN)
    
    @property
    def noStar(self):
        """True if this probe has no star visible."""
        return (self._bits & NOSTAR)
    @noStar.setter
    def noStar(self,value):
        self._set(NOSTAR) if value else self._unset(NOSTAR)
    
    @property
    def outOfFocus(self):
        """True if the star in this probe is out of focus."""
        return (self._bits & OUTOFFOCUS)
    @outOfFocus.setter
    def outOfFocus(self,value):
        self._set(OUTOFFOCUS) if value else self._unset(OUTOFFOCUS)

    @property
    def gprobebits(self):
        """The hex bitstring for this probe's status."""
        return "0x%x"self._bits
    @gprobebits.setter
    def gprobebits(self,value):
        if not isinstance(value,int):
            raise ValueError("gprobebits must be set as an integer!")
        else:
            self._bits = value

    def get_ref_mag(self):
        """
        Return the reference magnitude for this probe's target.

        The magnitude the guider should measure for this star/fiber at
        the current telescope position. Guider effective wavelength is
        5400A, so calculate guidermag from g and r. Then correct for
        atmospheric extinction.
        
        APO atmospheric extinction coeff at airmass=1 taken from table 3 of the
        ubercal paper, Padmanabhan et al. 2008 ApJ.
        with the k0 value for the filters, g:0.17, r:0.10

        Color terms for transformation from a*g + b*r = guidermag 
        a=xx, b=yy
        """
        actorState = myGlobals.actorState
        k0_g = 0.17
        k0_r = 0.10
        #get airmass from tcc only gives alt = tcc.axePos[2] 
        zd = 90. - actorState.models["tcc"].keyVarDict["axcPos"][1]
        # TBD: zd=0 never occurs for tracking, but need to test for zd=0 for simulate
        airmass = 1./math.cos(math.radians(zd))
        gobs = ugriz[1] + airmass*k0_g
        robs = ugriz[2] + airmass*k0_r
        #guidermag = xx*gobs + yy*robs
        return guidermag
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
        
        # Will contain [id]:probeInfo pairs
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

    def setGprobeState(self, fiber, enable=True, info=None, create=False):
        """Set a fiber's state"""

        if fiber in ("ACQUIRE", "GUIDE"):
            fiber_type = fiber
            for gp in self.gprobes.values():
                if gp.info.fiber_type == fiber_type:
                    gp.enabled = enable
        else:
            if not self.gprobes.has_key(fiber) and create:
                self.gprobes[fiber] = ProbeInfo(info)
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
            self.gcameraMagnification = gcameraMagnificatio

