"""
Old code for generating some guider plots with SM.

Torn from the guts of masterThread.guideStep().

!!!!!
jkp 2012.01.23: May not be fully tested
!!!!!

"""
import os.path
import numpy

try:
	import sm
except ImportError:
	print "Failed to import SM"
	sm = None

class PlotGuider(object):
    """Store parameters for plotting, and create requested guider plots."""
    def __init__(self,size,psPlot=False,psPlotDir="/data/gcam/scratch/"):
        """Pass the number of valid guide probes to set up the arrays."""
        self.psPlotDir = psPlotDir
        self.psPlot = psPlot
        self.size = size
        self.fiberid_np = numpy.zeros(size)
        self.raCenter_np = numpy.zeros(size)
        self.decCenter_np = numpy.zeros(size)
        self.dRA_np = numpy.zeros(size)
        self.dDec_np = numpy.zeros(size)
        self.x_np = numpy.zeros(size)
        self.xErr_np = numpy.zeros(size)
        self.d_np = numpy.zeros(size)
    #...
    
    def checkPlot(self):
        """
        Check that everything is ok for us to plot.
        Returns True if everything's ok, otherwise a message to
        send at warning level.
        """
        if sm:
            try:
                sm.device('X11')
            except:
                return "X display error, cannot open SM guider window"
        else:
            return "Unable to plot as SM is not available"
        
        if plot and psPlot:
            if not (os.path.exists(self.psPlotDir)):
                self.psPlot = False
                return "Unable to write SM hardcopies"
    
        return True
    #...

    def postscriptDevice(self, psPlotDir, frameNo, prefix=""):
        """Return the SM device to write the postscript file for guide frame frameNo"""
        return "postencap %s%d.eps" % (prefix, frameNo)    
    
    def plotOffsets(self,device,gState,frameInfo):
        """Generate a plot of the fiber offsets"""
        if device == "postscript":
            if not self.psPlot:
                return
            deviceCmd = self.postscriptDevice(self.psPlotDir, frameInfo.frameNo)
        else:
            deviceCmd = device

        sm.device(deviceCmd)

        sm.erase(False)
        sm.limits([-400, 400], [-400, 400])
        sm.box()
        sm.frelocate(0.5, 1.03)
        sm.putlabel(5, r"\1Offsets")
        sm.xlabel(r"\2\delta Ra")
        sm.ylabel(r"\2\delta Dec")

        sm.ptype([63])
        sm.frelocate(0.85, 0.95)
        sm.putlabel(5, r"\1Frame %d" % frameInfo.frameNo)
        vscale = 1000 # how much to multiply position error
        sm.relocate(-350, 350)
        asec = gState.plugPlateScale/3600.0
        sm.draw(-350 + vscale*gState.plugPlateScale/3600.0, 350)
        sm.label(r"  \raise-200{\1{1 arcsec}}")

        for i in range(len(self.fiberid_np)):
            if self.fiberid_np[i] == 0:
                continue

            sm.relocate(self.raCenter_np[i], self.decCenter_np[i])
            sm.putlabel(6, r" %d" % self.fiberid_np[i])

            sm.relocate(plotaArams.raCenter_np[i], self.decCenter_np[i])

            sm.dot()

            if not gState.gprobes[self.fiberid_np[i]].enabled:
                sm.ltype(1)
                sm.ctype(sm.CYAN)

            sm.draw(self.raCenter_np[i] + vscale*self.dRA_np[i],
                    self.decCenter_np[i] + vscale*self.dDec_np[i])

            sm.ltype()
            sm.ctype()

        sm.ptype()
    #...

    def plotFWHM(self,device,frameInfo):
        """Plot the FWHM of each fiber."""
        if device == "postscript":
            if not self.psPlot:
				return
            deviceCmd = postscriptDevice(self.psPlotDir, frameInfo.frameNo, "Focus")
        else:
            deviceCmd = device

        sm.device(deviceCmd)

        sm.erase(False)
        #
        # Bravely convert to FWHM in arcsec (brave because of the sqrt)
        #
        f = frameInfo.sigmaToFWHM/frameInfo.micronsPerArcsec
        X_np = f*numpy.sqrt(self.x_np)
        XErr_np = f*xErr_np/(2*numpy.sqrt(self.x_np))

        sm.limits(self.d_np, X_np)
        sm.box()
        sm.frelocate(0.5, 1.03)
        sm.putlabel(5, r"\1Focus")

        sm.xlabel(r"\2d_i \equiv{} fibre piston (\mu m)")
        sm.ylabel(r"\2\ssqrt{\sigma_i^2 - C d_i^2} (FWHM, arcsec)")

        sm.frelocate(0.85, 0.95)
        sm.putlabel(5, r"\1Frame %d" % frameInfo.frameNo)

        for i in range(len(self.fiberid_np)):
            if fiberid_np[i] == 0:
                continue

            sm.relocate(self.d_np[i], X_np[i])
            sm.putlabel(6, r" %d" % self.fiberid_np[i])

        for l in (2, 4):
            sm.errorbar(self.d_np, X_np, XErr_np, l)

        #dd_np = numpy.array([-1000, 1000])
        dd_np = numpy.arange(-1000, 1000, 25)

        sm.ctype(sm.CYAN)
        if x != None: # I.e. we successfully fit for focus
            sm.connect(dd_np, f*numpy.sqrt(x[0, 0] + x[1, 0]*dd_np))
            sm.frelocate(0.1, 0.1); sm.label(r"\line 0 2000 \colour{default} Best fit")

            sm.ctype(sm.MAGENTA)
            sm.frelocate(0.1, 0.15)
            if False:
                sm.label(r"\line 1 2000 \colour{default} Seeing")
                sm.ltype(1)
                sm.connect(dd_np, 0*dd_np + rms0*frameInfo.sigmaToFWHM)
            else:
                sm.label(r"{\2\apoint 45 4 1} \colour{default} Seeing")
                sm.relocate(0, rms0*frameInfo.sigmaToFWHM)
                sm.expand(4); sm.angle(45)
                sm.dot()
                sm.expand(); sm.angle(0)

        sm.ltype(); sm.ctype()

        del dd_np
    #...
#...
