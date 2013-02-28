"""
Classes related to the current state of the guider.
"""

import numpy as np
import math 

class ProbeInfo(object):
    """
    Contains information about a single guide probe.
    """
    def __init__(self, exists, enabled, xCenter, yCenter, radius, rotation,
                 xFerruleOffset, yFerruleOffset, focusOffset, fiber_type, flags):
        self.exists = exists
        self.enabled = enabled
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.radius = radius
        self.xFerruleOffset = xFerruleOffset
        self.yFerruleOffset = yFerruleOffset
        self.rotation = rotation
        self.focusOffset = focusOffset
        self.fiber_type = fiber_type
        self.rotStar2Sky = numpy.nan
        self.flags = flags
        self.haOffsetTimes = {}
        self.haXOffsets = {}
        self.haYOffsets = {}
        self.ugriz = np.nan

    def set_magnitude(self,ugriz):
        """
        Set the fibermag (ugriz array) for this probe.
        Fibermag originates in plPlugMapP.par, mag through 2arcsec fiber. 
        """
        self.ugriz = ugriz
    
    def get_ref_mag(self):
        """
        Reference magnitude for this probe's target is mag the guider 
        should measure for this star/fiber at the current telescope position
        Guider effective wavelength is 5400A, so calc guidermag from g and r.
        Then correct for atmospheric extinction
        
        APO atmospheric extinction coeff at airmass=1 taken from ubercal paper
        Padmanabhan etal. 2008 ApJ.  g:k0=0.17,  r:0.10

        Color terms for transformation from a*g + b*r = guidermag 
        a=xx, b=yy
        """
        roughmag = (ugriz[2] + ugriz[1])/2
        #get airmass form tcc  only gives alt = tcc.axePos[2] 
        zd = 90. - actorState.models["tcc"].keyVarDict["axcPos"][1]
        #zd=0 never occurs for tracking, but need to test for zd=0 for simulate
        airmass = 1./math.cos(math.radians(zd))
        gobs = ugriz[1] + airmass*0.17
        robs = ugriz[1] + airmass*0.10
        #guidermag = xx*gobs + yy*robs
        return roughmag
#...

class GuiderState(object):
    """
    The current state of the guider.
    
    Contains information about the currently loaded cartridge, the target of
    each guide probe, custom parameters (decentering, stacking, exposure time, etc.).
    
    Does not know about any gcamera exposures: that's in FrameInfo.
    """

    class Gprobe(object):
        """Contains information about a single guide probe."""
        def __init__(self, id, probeInfo, enable=True, flags=None):
            self.id = id
            self.probeInfo = probeInfo
            self.enabled = enable
            self.flags = flags

        def isEnabled(self):
            raise NotImplementedError()

        def setEnabled(self, enabled):
            if enabled:
                self.flags |= self._ENABLED
            else:
                self.flags &= ~self._ENABLED

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
        self.deleteAllGprobes()

        self.plugPlateScale = numpy.nan
        self.dSecondary_dmm = numpy.nan
        self.gcameraPixelSize = numpy.nan
        self.gcameraMagnification = numpy.nan
        
        self.setGuideMode("axes")
        self.setGuideMode("focus")
        self.setGuideMode("scale")
        self.setRefractionBalance(0.0)
        
        self.pid = {}               # PIDs for various axes
        for what in ["raDec", "rot", "scale", "focus"]:
            self.pid[what] = PID.PID(self.expTime, 0, 0, 0)

        self.decenter = False                #gstate only
        self.setDecenter("decenterRA")       
        self.setDecenter("decenterDec")      
        self.setDecenter("decenterRot")
        self.decenterChanged = True
        self.decenterFocus = numpy.nan
        self.decenterScale = numpy.nan

    def deleteAllGprobes(self):
        """Delete all fibers """
        self.gprobes = {}

    def setGprobeState(self, fiber, enable=True, info=None, create=False, flags=None):
        """Set a fiber's state"""

        if fiber in ("ACQUIRE", "GUIDE"):
            fiber_type = fiber
            for gp in self.gprobes.values():
                if gp.info.fiber_type == fiber_type:
                    gp.enabled = enable
        else:
            if not self.gprobes.has_key(fiber) and create:
                self.gprobes[fiber] = GuiderState.Gprobe(fiber, info, enable, flags)
            else:
                self.gprobes[fiber].enabled = enable

    def setGuideMode(self, what, enabled=True):
        if what == "axes":
            self.guideAxes = enabled
        elif what == "focus":
            self.guideFocus = enabled
        elif what == "scale":
            self.guideScale = enabled
        else:
            raise RuntimeError, ("Unknown guide mode %s" % what)

    def setRefractionBalance(self, value=0.0):
        self.refractionBalance = value
        
    def setCmd(self, cmd=None):
        self.guideCmd = cmd

    def setDecenter(self, what, value=0):
        if what == "decenterRA":
            self.decenterRA = value
        elif what == "decenterDec":
            self.decenterDec = value
        elif what == "decenterRot":
            self.decenterRot = value
        else:
            raise RuntimeError, ("Unknown decenter axis name %s" % what)

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

