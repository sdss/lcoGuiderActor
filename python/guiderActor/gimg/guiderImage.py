"""
For processing guider images, indepdendent of commanding the telescope.

Can be used to batch-process old guider images to extract more information
from them in bulk.
"""

import os.path
from operator import attrgetter
import ctypes
import datetime
import re

import pyfits
import numpy as np

from scipy.ndimage.morphology import binary_closing, binary_dilation, binary_erosion
from scipy.ndimage.measurements import label, center_of_mass

import PyGuide

import GuiderExceptions
import actorcore.utility.fits as actorFits
from opscore.utility.tback import tback
from opscore.utility.qstr import qstr

class Fiber(object):
    """A guider fiber and the star image seen through it."""
    def __init__(self, fiberid, xc=np.nan, yc=np.nan, r=np.nan, illr=np.nan, label=-1):
        self.fiberid = fiberid
        self.xcen = xc
        self.ycen = yc
        self.radius = r
        self.illrad = illr

        self.xs = np.nan
        self.ys = np.nan
        self.xyserr = np.nan
        self.fwhm = np.nan
        self.sky = np.nan
        self.skymag = np.nan
        self.flux = np.nan
        self.mag = np.nan

        self.fwhmErr = np.nan
        self.dx = np.nan
        self.dy = np.nan
        self.dRA = np.nan
        self.dDec = np.nan

        self.gProbe = None
        
        # The label of this fiber in the labeled, masked image (from scipy.ndimage.label)
        self.label = label

    def __str__(self):
        return ('Fiber id %i: center %g,%g; star %g,%g; radius %g' %
                (self.fiberid, self.xcen, self.ycen, self.xs, self.ys, self.radius))

    def set_fake(self):
        self.xcen = np.nan

    def is_fake(self):
        return np.isnan(self.xcen)
#...

# The following are ctypes classes for interaction with the ipGguide.c code.
class REGION(ctypes.Structure):
    _fields_ = [("nrow", ctypes.c_int),
                ("ncol", ctypes.c_int),
                ("rows_s16", ctypes.POINTER(ctypes.POINTER(ctypes.c_int16)))]

class MASK(ctypes.Structure):
    _fields_ = [("nrow", ctypes.c_int),
                ("ncol", ctypes.c_int),
                ("rows", ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)))]

class FIBERDATA(ctypes.Structure):
    _fields_ = [("g_nfibers", ctypes.c_int),
                ("g_fid", ctypes.POINTER(ctypes.c_int32)),
                ("g_xcen", ctypes.POINTER(ctypes.c_double)),
                ("g_ycen", ctypes.POINTER(ctypes.c_double)),
                ("g_fibrad", ctypes.POINTER(ctypes.c_double)),
                ("g_illrad", ctypes.POINTER(ctypes.c_double)),
                ("g_xs", ctypes.POINTER(ctypes.c_double)),
                ("g_ys", ctypes.POINTER(ctypes.c_double)),
                ("flux", ctypes.POINTER(ctypes.c_double)),
                ("sky", ctypes.POINTER(ctypes.c_double)),
                ("fwhm", ctypes.POINTER(ctypes.c_double)),
                ("poserr", ctypes.POINTER(ctypes.c_double)),
                ("g_readnoise", ctypes.c_double),
                ("g_npixmask", ctypes.c_int)]

# Must match ipGguide.h
FWHM_BAD = 99.99

def np_array_to_REGION(A):
    H, W = A.shape
    ptrtype = ctypes.POINTER(ctypes.c_int16)
    rows = (ptrtype*H)(*[row.ctypes.data_as(ptrtype) for row in A])
    return REGION(H, W, rows)

def np_array_to_MASK(A):
    H, W = A.shape
    ptrtype = ctypes.POINTER(ctypes.c_uint8)
    rows = (ptrtype*H)(*[row.ctypes.data_as(ptrtype) for row in A])
    return MASK(H, W, rows)

def bin_image(img, BIN):
    """Return an image rebinned by BINxBIN pixels."""
    binned = np.zeros((img.shape[0]/BIN, img.shape[1]/BIN), np.float32)
    for i in range(BIN):
        for j in range(BIN):
            binned += img[i::BIN,j::BIN]
    binned /= (BIN*BIN)
    return binned


class GuiderImageAnalysis(object):
    '''
    A class for analyzing the images taken by the guiding camera.
    
    guiderImageAnalysis = GuiderImageAnalysis()
    # for each new exposure:
        gimg_filename = 'gimg-0400.fits'
        # cmd must have methods "inform(string)", "warn(string)".
        fibers = guiderImageAnalysis(gimg_filename, gState.gprobes, cmd)
        # run guider control loop...
        guiderImageAnalyze.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)
    '''

    # According to
    #  http://sdss3.apo.nmsu.edu/opssoft/guider/ProcessedGuiderImages.html
    #   0 = good
    mask_saturated = 1
    mask_badpixels = 2
    mask_masked    = 4 # ie, outside the guide fiber

    def __init__(self,setPoint):
        """
        New GuiderImageAnalysis instances are ready to accept files for processing.
        setPoint is the current gcamera temperature set point.
        """
        self.setPoint = setPoint
        self.outputDir = ''
        self.camera = 'gcamera'
        # set during findStars():
        self.fibers = None
        self.guiderImage = None
        self.guiderHeader = None
        self.maskImage = None
        self.bypassDark = False

        # So that we don't have to re-open dark and flat files every time.
        self.currentDarkName = ''
        self.currentFlatName = ''
        self.processedDark = None
        self.processedFlat = None
        self.darkTemperature = None

        # Print debugging?
        self.printDebug = False
        
        # amount we let the gcamera temperature vary from the setPoint
        self.deltaTemp = 3.0
        
        # "Big" (acquisition) fibers are bigger than this pixel
        # radius.  The older SDSS cartridges don't declare (in the
        # gcamFiberInfo.par file) the larger fibers to be ACQUIRE,
        # though they are declared to have radii of 14.1 pixels (vs
        # 8.5 for the GUIDE fibers).  We therefore cut on this radius.
        self.bigFiberRadius = 12.

        # Saturation level.
        #need to make use of full 64k image for bright marvels guide stars.
        #A solution was to scale by 2 data into and out of the C code
        #Problems also with the rotation of postage stamps
        # jkp: note that the saturation level is a bit below the
        # full well level of the chip, to prevent column bleedthrough.
        #self.saturationLevel = 0xA000     #41k
        self.saturationLevel = 0xF000    #62,000
                
        # The value to replace saturated pixels by.
        #self.saturationReplacement = 0xA000
        self.saturationReplacement = 0xF000 #62,000

        # The factor by which guider images are binned down.
        # That is, unbinned (flat) images are this factor bigger in
        # each dimension.
        self.binning = 2

        # The pixel scale of the guider camera, when binned by "binning"
        # In arcsec/pixel
        self.pixelscale = 0.428

        # The photometric zero-point for (g + r)/2 band was 25.34
        # This was calibrated vs MJD 55246, with AZ~=73 deg, airmass~=1.046
        # Masayuki zero average point 25.70 for use for his color transform
        self.zeropoint = 25.70
    #...
    
    def __call__(self, cmd, gimgfn, gprobes, setPoint, bypassDark=False, camera='gcamera'):
        """
        Calls findStars to process gimgfn/gprobes and return found fibers.
        
        gimgfn is the unprocessed gcamera file to process.
        gprobes is from GuiderState: it's a dict of probeId to GProbe object.
        cmd is a Commander object, to allow messaging (diag/inform/warn).
        setPoint is the current gcamera temperature set point.
        set bypassDark to ignore guider dark frames, and not do dark subtraction.
        camera can be either 'gcamera' to find fibers for guiding, or
        'ecamera' to find the brightest star for a pointing model.
        """
        self.gimgfn = gimgfn
        # don't like the re, but it works (stolen from masterThread.guideStep)
        self.frameNo = int(re.search(r"([0-9]+)\.fits.*$", gimgfn).group(1))
        self.setPoint = setPoint
        self.bypassDark = bypassDark
        self.cmd = cmd
        self.camera = camera
        
        return self.findStars(gprobes)
    #...

    def pixels2arcsec(self, pix):
        """Convert pix to arcseconds, using the pixelscale."""
        return pix * self.pixelscale

    def flux2mag(self, flux, exptime):
        """Convert flux to magnitude, using exptime."""
        if exptime == 0:
            return -99
        return -2.5 * np.log10(flux / exptime) + self.zeropoint

    def find_bias_level(self,image,binning=1):
        """
        Find the bias level of the image.
        Set binning to the number of binned pixels in x and y.
        
        If there is no overscan region (expected to be 24 extra columns unbinned)
        then we use a kludged overscan.
        """
        # The overscan region is an extra 24 columns (12 after binning)
        if image.shape[1] == 1048/binning:
            # subtracting the overall median value should be good enough
            # NOTE: The the "inner" portion of the bias region is most
            # representative of the rest of the chip.
            # TBD: Which is correct here?
            #bias = np.median(image[:,(1024/binning):(1038/binning)])
            # TBD: The upper is what Renbin wanted,
            # TBD: the lower is what makes sense from over-exposed flats and the ecam.
            bias = np.median(image[:,(1039/binning):])
        else:
            # find bias = BIAS_PERCENTILE (ipGguide.h) = (100 - 70%)
            self.cmd.warn('text=%s'%qstr("Cheating with bais level! No overscan was found!"))
            ir = image.ravel()
            I = np.argsort(ir)
            bias = ir[I[int(0.3 * len(ir))]]
        self.imageBias = bias
    #...
    
    def ensureLibraryLoaded(self):
        """
        Load C library that does the fluxing, etc.: lib/libguide.so -> self.libguide
        See src/ipGguide.c for the actual calculations.
        """
        path = os.path.expandvars("$GUIDERACTOR_DIR/lib/libguide.so")
        libguide = ctypes.CDLL(path)
        if not libguide:
            self.cmd.error('text=%s'%qstr('Failed to load "libguide.so" from %s ($GUIDERACTOR_DIR/lib/libguide.so)' % path))
        libguide.gfindstars.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(FIBERDATA), ctypes.c_int]
        libguide.gfindstars.restype = ctypes.c_int
        libguide.fiberdata_new.argtypes = [ctypes.c_int]
        libguide.fiberdata_new.restype = ctypes.POINTER(FIBERDATA)
        libguide.fiberdata_free.argtypes = [ctypes.POINTER(FIBERDATA)]
        libguide.fiberdata_free.restype = None
        libguide.rotate_region.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(REGION), ctypes.c_float]
        libguide.rotate_region.restype = None
        libguide.rotate_mask.argtypes = [ctypes.POINTER(MASK), ctypes.POINTER(MASK), ctypes.c_float]
        libguide.rotate_mask.restype = None
        self.libguide = libguide

    def findDarkAndFlat(self, gimgfn, fitsheader):
        """ findDarkAndFlat(...)

        Returns the filenames containing the dark and flat images for
        the given guider-camera image.

        These are listed in the gimg-####.fits header; this method exists
        to make testing easier, and to make path name magic more explicit.

        DARKFILE= '/data/gcam/55205/gimg-0003.fits.gz'
        FLATFILE= '/data/gcam/55205/gimg-0224.fits.gz'
        """
        # files prior to MJD 56465 have the dark/flat without .gz in the header.
        darkfile = fitsheader['DARKFILE']
        if not os.path.exists(darkfile):
            darkfile = darkfile+'.gz'
        flatfile = fitsheader.get('FLATFILE', None)
        if flatfile is not None:
            if not os.path.exists(flatfile):
                flatfile = flatfile+'.gz'
        return darkfile,flatfile

    def getProcessedOutputName(self, imgfn):
        """Return the name of the file that we will save the processed results to."""
        
        (dirname, filename) = os.path.split(imgfn)
        outname = 'proc-%s' % (filename)
        if self.outputDir:
            thedir = self.outputDir
        else:
            thedir = dirname
        procpath = os.path.join(thedir, outname)
        return procpath

    def addPixelWcs(self, header, wcsName=""):
        """Add a WCS that sets the bottom left pixel's centre to be (0.5, 0.5)"""
        header.update("CRVAL1%s" % wcsName, 0, "(output) Column pixel of Reference Pixel")
        header.update("CRVAL2%s" % wcsName, 0, "(output) Row pixel of Reference Pixel")
        header.update("CRPIX1%s" % wcsName, 0.5, "Column Pixel Coordinate of Reference")
        header.update("CRPIX2%s" % wcsName, 0.5, "Row Pixel Coordinate of Reference")
        header.update("CTYPE1%s" % wcsName, "LINEAR", "Type of projection")
        header.update("CTYPE1%s" % wcsName, "LINEAR", "Type of projection")
        header.update("CUNIT1%s" % wcsName, "PIXEL", "Column unit")
        header.update("CUNIT2%s" % wcsName, "PIXEL", "Row unit")

    def fillPrimaryHDU(self, cmd, models, imageHDU, frameInfo, objectname):
        """ Add in all the site and environment keywords. """
        try:
            imageHDU.header.update('SEEING', frameInfo.seeing if np.isfinite(frameInfo.seeing) else 0.0,
                                   'Estimate of current seeing, arcsec fwhm')
            imageHDU.header.update('OBJECT', objectname, '')
            imageHDU.header.update('GCAMSCAL', frameInfo.guideCameraScale, 'guide camera plate scale (mm/pixel)')
            imageHDU.header.update('PLATSCAL', frameInfo.plugPlateScale, 'plug plate scale (mm/degree)')
                        # Do this first, before we need the models.
            guiderCards = self.getGuideloopCards(cmd, frameInfo)
            actorFits.extendHeader(cmd, imageHDU.header, guiderCards)
            self.addPixelWcs(imageHDU.header)

            plateCards = actorFits.plateCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, plateCards)
            guiderCards = actorFits.guiderCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, guiderCards)
            apoCards = actorFits.apoCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, apoCards)

        except Exception as e:
            self.cmd.error('text=%s'%qstr('!!!!! failed to fill out primary HDU  !!!!! (%s)' % (e)))

            raise e
    def getGuideloopCards(self, cmd, frameInfo):
        defs = (('dRA', 'DRA', 'measured offset in RA, deg'),
                ('dDec', 'DDec', 'measured offset in Dec, deg'),
                ('dRot', 'DRot', 'measured rotator offset, deg'),
                ('dFocus', 'DFocus', 'measured focus offset, um '),
                ('dScale', 'DScale', 'measured scale offset, %'),
                ('filtRA', 'FILTRA', 'filtered offset in RA, deg'),
                ('filtDec', 'FILTDec', 'filtered offset in Dec, deg'),
                ('filtRot', 'FILTRot', 'filtered rotator offset, deg'),
                ('filtFocus', 'FILTFcus', 'filtered focus offset, um '),
                ('filtScale', 'FILTScle', 'filtered scale offset, %'),
                ('offsetRA', 'OFFRA', 'applied offset in RA, deg'),
                ('offsetDec', 'OFFDec', 'applied offset in Dec, deg'),
                ('offsetRot', 'OFFRot', 'applied rotator offset, deg'),
                ('offsetFocus', 'OFFFocus', 'applied focus offset, um'),
                ('offsetScale', 'OFFScale', 'applied scale offset, %'),
                ('guideRMS',    'gdRMS',    'RMS guiding error total, arcsec'),
                ('nguideRMS',    'ngdRMS',   'N stars used for RMS'),
                ('guideXRMS',    'gdXRMS',   'CCD X component of guiding RMS, arcsec'),
                ('guideYRMS',    'gdYRMS',   'CCD Y component of guiding RMS, arcsec'),
                ('guideAzRMS',    'gdAzRMS',  'Az component of guiding RMS error, arcsec'),
                ('guideAltRMS', 'gdAltRMS', 'Alt component of guiding RMS error, arcsec'),
                ('guideRaRMS',    'gdRaRMS',  'RA component of guiding RMS error, arcsec'),
                ('guideDecRMS', 'gdDecRMS', 'Dec component of guiding RMS error, arcsec'),
                ('guideFitRMS',  'gdFRMS',  'RMS of fit to guide star posn, arcsec'),
                ('nguideFitRMS', 'ngdFRMS', 'N stars used for fit RMS'),
                ('decenterRA',    'dcnRA',   'applied user supplied offset in RA, arcsec'),
                ('decenterDec',   'dcnDec',  'applied user supplied offset in Dec, arcsec'),
                ('decenterRot',   'dcnRot',  'applied user supplied rotator offset, arcsec'),
                ('decenterFocus', 'dcnFcus', 'applied user supplied focus offset, um'),
                ('decenterScale', 'dcnScle', 'applied user supplied scale offset, %' ),
                ('refractionBalance','refrBal','specified refraction balance between (0,1)'),
                )
                #TBD: FIXME PH --- do we change to 1e6 units for scale
        cards = []
        for name, fitsName, comment in defs:
            try:
                val = None
                val = getattr(frameInfo, name)
                #print name, val
                if np.isnan(val):
                    val = -99999.9 # F.ing F.TS
                c = actorFits.makeCard(cmd, fitsName, val, comment)
                cards.append(c)
            except Exception as e:
                self.cmd.warn('text=%s'%qstr('failed to make guider card %s=%s (%s)' % (name, val, e)))
        return cards

    def getStampHDUs(self, fibers, bg, image, mask):
        """Get the FITS HDUs for the image stamps."""
        if len(fibers) == 0:
            return [pyfits.ImageHDU(np.array([[]]).astype(np.int16)),
                    pyfits.ImageHDU(np.array([[]]).astype(np.uint8))]
        self.ensureLibraryLoaded()
        r = int(np.ceil(max([f.radius for f in fibers])))
        stamps = []
        maskstamps = []
        for f in fibers:
            xc = int(f.xcen + 0.5)
            yc = int(f.ycen + 0.5)
            rot = -f.gProbe.rotStar2Sky
            self.cmd.diag('text=%s'%qstr("rotating fiber %d at (%d,%d) by %0.1f degrees" % (f.fiberid, xc, yc, rot)))
            # Rotate the fiber image...
            stamp = image[yc-r:yc+r+1, xc-r:xc+r+1].astype(np.int16)
            rstamp = np.zeros_like(stamp)
            self.libguide.rotate_region(np_array_to_REGION(stamp),
                                        np_array_to_REGION(rstamp),
                                        rot)
            stamps.append(np.flipud(rstamp))
            # Rotate the mask image...
            stamp = mask[yc-r:yc+r+1, xc-r:xc+r+1].astype(np.uint8)
            #print 'stamp values:', unique(stamp.ravel())
            # Replace zeros by 255 so that after rotate_region (which puts
            # zeroes in the "blank" regions) we can replace them.
            stamp[stamp == 0] = 255
            rstamp = np.zeros_like(stamp)
            self.libguide.rotate_mask(np_array_to_MASK(stamp),
                                      np_array_to_MASK(rstamp),
                                      rot)
            # "blank" regions become masked-out pixels.
            rstamp[rstamp == 0] = GuiderImageAnalysis.mask_masked
            # Reinstate the zeroes.
            rstamp[rstamp == 255] = 0
            #print 'rotated stamp values:', unique(rstamp.ravel())
            maskstamps.append(np.flipud(rstamp))

        # Stack the stamps into one image.
        stamps = np.vstack(stamps)
        maskstamps = np.vstack(maskstamps)
        # Replace zeroes by bg
        stamps[stamps==0] = bg
        # Also replace masked regions by bg.
        # HACKery....
        #stamps[maskstamps > 0] = maximum(bg - 500, 0)
        stamps[maskstamps > 0] = bg
        #print 'bg=', bg
        return [pyfits.ImageHDU(stamps), pyfits.ImageHDU(maskstamps)]

    def _get_basic_hdulist(self, image, primhdr, bg):
        """Return an hdulist with the most basic header keywords filled in."""
        imageHDU = pyfits.PrimaryHDU(image, header=primhdr)
        imageHDU.header.update('IMGBACK', bg, 'crude background for entire image. For displays.')
        imageHDU.header.update('OVERSCAN', self.imageBias, 'Bias level of this image.')
        # TBD: sdssfmt should go away, but having it missing breaks the current STUI's 
        # ability to show plateview. putting it back for now, just to solve
        # that without releasing a new STUI.
        imageHDU.header.update('SDSSFMT', 'GPROC 1 4', 'type major minor version for this file')
        hdulist = pyfits.HDUList()
        hdulist.append(imageHDU)
        return hdulist

    def _getProcGimgHDUList(self, primhdr, gprobes, fibers, image, mask, stampImage=None):
        """Generate an HDU list to be inserted into the proc-file header."""
        if stampImage is None:
            stampImage = image
        bg = np.median(image[mask == 0])
        hdulist = self._get_basic_hdulist(image,primhdr,bg)

        try:
            # List the fibers by fiber id.
            fiberids = sorted(gprobes.keys())
            sfibers = []
            myfibs = dict([(f.fiberid,f) for f in fibers])
            for fid in fiberids:
                if fid in myfibs:
                    sfibers.append(myfibs[fid])
                else:
                    # Put in NaN fake fibers for ones that weren't found.
                    fake = Fiber(fid)
                    fake.gProbe = gprobes[fid]
                    fake.gProbe.disabled = True
                    sfibers.append(fake)
            fibers = sfibers

            # Split into small and big fibers.
            # Also create the columns "stampInds" and "stampSizes"
            # for the output table.
            smalls = []
            bigs = []
            stampInds = []
            stampSizes = []
            for f in fibers:
                if f.is_fake():
                    stampSizes.append(0)
                    stampInds.append(-1)
                elif f.radius < self.bigFiberRadius:
                    # must do this before appending to smalls!
                    stampInds.append(len(smalls))
                    stampSizes.append(1)
                    smalls.append(f)
                else:
                    stampInds.append(len(bigs))
                    stampSizes.append(2)
                    bigs.append(f)

            hdulist.append(pyfits.ImageHDU(mask))
            hdulist += self.getStampHDUs(smalls, bg, stampImage, mask)
            hdulist += self.getStampHDUs(bigs, bg, stampImage, mask)

            # jkp TBD: rework this to make it more legible/easily extensible.
            pixunit = 'guidercam pixels (binned)'
            gpinfofields = [
                           # python attr, FITScol (if diff), FITStype, NIL val, unit
                           ('exists',         None,    'L',   np.nan,     None),
                           ('enabled',        None,    'L',   np.nan,     None),
                           ('gprobebits',     None,    'B',   np.nan,     None),
                           ('xFocal',         None,    'E',   np.nan,     'plate mm'),
                           ('yFocal',         None,    'E',   np.nan,     'plate mm'),
                           ('radius',         None,    'E',   np.nan,     pixunit),
                           ('xFerruleOffset', None,    'E',   np.nan,     None),
                           ('yFerruleOffset', None,    'E',   np.nan,     None),
                           ('rotation',       None,    'E',   np.nan,     None),
                           ('rotStar2Sky',    None,    'E',   np.nan,     None),
                           ('focusOffset',    None,    'E',   np.nan,     'micrometers'),
                           ('fiber_type',      None,    'A20', np.nan,     None),
                           ('ugriz',          None,    '5E',  [np.nan,]*5, 'mag'),
                           ('ref_mag',        None,    'E',   np.nan,     'synthetic predicted fiber mag'),
                           ]

            ffields = [
                      # FITScol,  variable, FITStype, unit
                      ('fiberid', None,     'I', None),
                      ('xCenter', 'xcen',   'E', pixunit),
                      ('yCenter', 'ycen',   'E', pixunit),
                      ('xstar',   'xs',     'E', pixunit),
                      ('ystar',   'ys',     'E', pixunit),
                      ('dx',      None,     'E', 'residual in mm, guidercam frame'),
                      ('dy',      None,     'E', 'residual in mm, guidercam frame'),
                      ('dRA',     None,     'E', 'residual in mm, plate frame'),
                      ('dDec',    None,     'E', 'residual in mm, plate frame'),
                      ('fwhm',    None,     'E', 'arcsec'),
                      ('flux',    None,     'E', 'total DN'),
                      ('mag',     None,     'E', 'mag'),
                      ('sky',     None,     'E', 'DN/pixel'),
                      ('skymag',  None,     'E', 'mag/(arcsec^2)'),
                      ('poserr',  'xyserr', 'E', pixunit),
                      ]

            # TBD: FIXME -- rotStar2Sky -- should check with "hasattr"...
            cols = []
            for name,fitsname,fitstype,nilval,units in gpinfofields:
                if fitsname is None:
                    fitsname = name
                cols.append(pyfits.Column(name=fitsname, format=fitstype, unit=units,
                                          array=np.array([getattr(f.gProbe, name, nilval) for f in fibers])))
            for name,atname,fitstype,units in ffields:
                cols.append(pyfits.Column(name=name, format=fitstype, unit=units,
                                          array=np.array([getattr(f, atname or name) for f in fibers])))

            cols.append(pyfits.Column(name='stampSize', format='I', array=np.array(stampSizes)))
            cols.append(pyfits.Column(name='stampIdx', format='I', array=np.array(stampInds)))

            hdulist.append(pyfits.new_table(cols))
        except Exception as e:
            self.cmd.warn('text=%s'%qstr('could not create header for proc- guider file: %s' % (e,)))
            tback('guiderImage write', e)
            raise e
        return hdulist


    def writeFITS(self, models, cmd, frameInfo, gprobes, output_verify='warn'):
        """
        Write a fits file containing the processed results for this exposure.
        """
        if not self.fibers and self.camera != 'ecamera':
            raise Exception('must call findStars() before writeFITS()')
        image = self.guiderImage
        hdr = self.guiderHeader

        procpath = self.getProcessedOutputName(self.gimgfn)
        objectname = os.path.splitext(self.gimgfn)[0]

        try:
            if self.camera == 'gcamera':
                hdulist = self._getProcGimgHDUList(hdr, gprobes, self.fibers, image, self.maskImage)
            elif self.camera == 'ecamera':
                bg = np.median(image)#TBD: this is a poor choice for star-filled ecam images!
                hdulist = self._get_basic_hdulist(image, hdr, bg)
                hdulist.append(pyfits.ImageHDU(self.maskImage))
            imageHDU = hdulist[0]
            self.fillPrimaryHDU(cmd, models, imageHDU, frameInfo, objectname)
            directory,filename = os.path.split(procpath)
            actorFits.writeFits(cmd, hdulist, directory, filename, doCompress=True, chmod=0644, checksum=True, output_verify=output_verify)
            self.cmd.inform('file=%s/,%s' % (directory, filename))
        except Exception as e:
            cmd.error('text="failed to write FITS file %s: %r"' % (procpath, e))
            raise e
    
    def _check_ccd_temp(self,header):
        """Return True if the gcamera CCDTEMP is within deltaTemp of setPoint and darkTemperature."""
        def tempCheck(temp1,temp2,delta):
            return ((temp1 - delta) < temp2 < (temp1 + delta))

        ccdtemp = header['CCDTEMP']
        imageType = header['IMAGETYP']
        result = True
        warnText = 'CCD temp signifcantly different (>%.1f) from %s: %.1f, expected %.1f'
        try:
            if not tempCheck(self.setPoint,ccdtemp,self.deltaTemp):
                self.cmd.warn('text=%s'%qstr(warnText%(self.deltaTemp, 'setPoint', ccdtemp, self.setPoint)))
                result = False
        except TypeError:
            msg = 'unknown error when checking setPoint (%s) against exposure ccdtemp (%s).'%(self.setPoint,ccdtemp)
            # this problem sometimes occurs after the guider is restarted.
            if self.setPoint is None:
                msg = 'setPoint is None: issue e/gcamera status and try again.'
            raise GuiderExceptions.GuiderError(msg)

        # redundant for darks, irrelevant for flats (we don't dark subtract them)
        if imageType != 'dark' and imageType != 'flat' and not tempCheck(self.darkTemperature,ccdtemp,self.deltaTemp):
            self.cmd.warn('text=%s'%qstr(warnText%(self.deltaTemp, 'dark temp', ccdtemp, self.darkTemperature)))
            result = False
        return result
    
    def _pre_process(self,filename,binning=1):
        """
        Initial checks and processing on any kind of exposure.
        
        Returns image,hdr,sat if everything goes well, raises exceptions if not.
        """
        assert(filename)
        self.ensureLibraryLoaded()
        
        # files prior to MJD 56465 have the dark/flat without .gz in the header.
        if not os.path.exists(filename):
            filename += '.gz'
        # Load guider-cam image.
        self.cmd.diag('text=%s'%qstr('Reading guider-cam image %s' % filename))
        try:
            image,hdr = pyfits.getdata(filename,0,header=True)
        except IOError:
            raise GuiderExceptions.GuiderError('File not found: %s'%filename)
        
        # Find saturation level pre-bias.
        sat = (image.astype(int) >= self.saturationLevel)
        nSat = sat.sum()
        
        image = self.applyBias(image,binning)
        # Occasionally there is a bad read from the camera.
        # In this case, the bias level is ~35,000, and the stddev is low.
        # We can just reject such frames, as they are useless.
        if self.imageBias > 5000:
            self.cmd.warn('text=%s'%qstr('Bad guider read! This exposure is mangled and will not be used.'))
            raise GuiderExceptions.BadReadError

        if nSat > 0:
            self.cmd.warn('text=%s'%qstr('Guider raw exposure has %i saturated pixels.' % nSat))
        if nSat > 4000:
            self.cmd.error('text=%s'%qstr('Fully saturated! Please reduce exposure time or wait for the excess light to go away.'))
            raise GuiderExceptions.GuiderError
        image[sat] = self.saturationReplacement
        
        return image,hdr,sat
    #...

    def _find_stars_ecam(self, image, mask):
        """Find the stars in a processed ecamera image."""
        # TBD: readnoise and ccd gain should be set in the config file,
        # or even better, written by gcameraICC and read from the header.
        readNoise = 10.4 # electrons RMS, from: http://www.ccd.com/alta_f47.html
        ccdGain = 1
        ccdInfo = PyGuide.CCDInfo(self.imageBias,readNoise,ccdGain,)
        # set a high threshold, since we only care about obviously bright stars.
        # TBD: this threshold could be a configurable parameter, so the observer
        # can adjust it in the Pointing Data STUI script, and relax it when focusing.
        stars = PyGuide.findStars(image,mask,None,ccdInfo,thresh=5)
        try:
            star = stars[0][0]
        except:
            return None,None
        try:
            shape = PyGuide.StarShape.starShape(image,mask,star.xyCtr,100)
        except:
            shape = PyGuide.StarShape.StarShapeData(False,'failed to compute star shape')
        return star,shape

    def findStars(self, gprobes):
        """
        Identify the centers of the stars in the fibers.
        Assumes self.gimgfn is set to the correct file name.
        
        gState is a GuiderState instance, containing a gprobes dict.
        
        Analyze a single image from the guider camera, finding fiber
        centers and star centers.

        The fiber centers are found by referring to the associated guider
        flat.  Since this flat is shared by many guider-cam images, the
        results of analyzing it are cached in a "proc-gimg-" file for the flat.

        Returns a list of fibers; also sets several fields in this object.

        The list of fibers contains an entry for each fiber found.
        """
        image,hdr,sat = self._pre_process(self.gimgfn,binning=self.binning)

        exptime = hdr.get('EXPTIME', 0)

        (darkFileName, flatFileName) = self.findDarkAndFlat(self.gimgfn, hdr)
        image = self.applyDark(image,darkFileName,exptime)

        # Check after we've loaded the dark, to ensure the dark temperature was set.
        self._check_ccd_temp(hdr)
        self.cmd.diag('text=%s'%qstr('Using flat image: %s' % flatFileName))
        if flatFileName != self.currentFlatName:
            try:
                self.analyzeFlat(flatFileName, gprobes)
            except GuiderExceptions.FlatError as e:
                # e.g.: no fibers could be found in the flat
                self.cmd.warn('text=%s'%qstr('Error processsing flat!'))
                raise e
        fibers = [f for f in self.flatFibers if not f.is_fake()]
        # mask the saturated pixels with the appropriate value.
        mask = self.flatMask.copy()
        mask[sat] |= GuiderImageAnalysis.mask_saturated

        # Divide by the flat (avoiding NaN where the flat is zero)
        image /= (self.flatImage + (self.flatImage == 0)*1)
        self.cmd.diag('text=%s'%qstr('After flattening: image range: %g to %g' % (image.min(), image.max())))

        # NOTE: jkp: post-flat fielding, we need to re-check for saturated pixels and remask them
        sat_flat = (image.astype(int) >= self.saturationLevel)
        image[sat_flat] = self.saturationReplacement
        mask[sat_flat] |= GuiderImageAnalysis.mask_saturated

        # Save the processed image
        self.guiderImage = image/2.0  #PH***quick Kluge needs to be corrected
        self.guiderHeader = hdr
        self.maskImage = mask

        # TBD: once we've converted to use PyGuide, we can get rid of these kludges!

        # Prepare image for Gunn C code PSF fitting in "gfindstars"...
        # PH..Kluge to conserve full dynamic range, scale image to fit into signed int in C code.
        # Have to scale up the outputs.
        img = image[:]/2.0     #

        # Zero out parts of the image that are masked out.
        # In this mask convention, 0 = good, >0 is bad.
        # Mark negative pixels
        badpixels = (img < 0) & (mask == 0)
        mask[badpixels] |= GuiderImageAnalysis.mask_badpixels
        # Blank out masked pixels.
        img[mask > 0] = 0

        if self.camera == 'ecamera':
            # mask the overscan too, since we're keeping it around for monitoring.
            mask[:,(1039/self.binning):] |= GuiderImageAnalysis.mask_saturated
            ecam_mask = mask == GuiderImageAnalysis.mask_saturated
            star,shape = self._find_stars_ecam(image, ecam_mask)
            if star is None:
                self.cmd.warn('text="ecam_star=None in exposure %d!"'%(self.frameNo))
                return []
            if not shape.isOK:
                self.cmd.warn('text="Could not determine star shape: %s"'%shape.msgStr)
            self.cmd.inform('ecam_star={},{:.2f},{:.2f},{:.3f},{:.2f},{:.1f}'.format(self.frameNo,star.xyCtr[0],star.xyCtr[1],shape.fwhm,shape.bkgnd,shape.ampl))
            return [] # no fibers to return
        else:
            # The "img16" object must live until after gfindstars() !
            img16 = img.astype(np.int16)
            c_image = np_array_to_REGION(img16)

            goodfibers = [f for f in fibers if not f.is_fake()]
            c_fibers = self.libguide.fiberdata_new(len(goodfibers))

            for i,f in enumerate(goodfibers):
                c_fibers[0].g_fid[i] = f.fiberid
                c_fibers[0].g_xcen[i] = f.xcen
                c_fibers[0].g_ycen[i] = f.ycen
                c_fibers[0].g_fibrad[i] = f.radius
                # FIXME ??
                c_fibers[0].g_illrad[i] = f.radius
            # TBD: FIXME --
            #c_fibers.readnoise = ...

            # mode=1: data frame; 0=spot frame
            mode = 1
            res = self.libguide.gfindstars(ctypes.byref(c_image), c_fibers, mode)
            # SH_SUCCESS is this following nutty number...
            if np.uint32(res) == np.uint32(0x8001c009):
                self.cmd.diag('text=%s'%qstr('gfindstars returned successfully.'))
            else:
                self.cmd.warn('text=%s'%qstr('gfindstars() returned an error code: %08x (%08x; success=%08x)' % (res, np.uint32(res), np.uint32(0x8001c009))))

            # pull star positions out of c_fibers, stuff outputs...
            for i,f in enumerate(goodfibers):
                f.xs     = c_fibers[0].g_xs[i]
                f.ys     = c_fibers[0].g_ys[i]
                f.xyserr = c_fibers[0].poserr[i]
                fwhm     = c_fibers[0].fwhm[i]
                if fwhm != FWHM_BAD:
                    f.fwhm = self.pixels2arcsec(fwhm)
                # else leave fwhm = nan.
                # TBD: FIXME -- figure out good units -- mag/(pix^2)?
                f.sky = (c_fibers[0].sky[i])*2.0    #correct for image div by 2 for Ggcode
                f.flux = (c_fibers[0].flux[i])*2.0
                if f.flux > 0:
                    f.mag = self.flux2mag(f.flux, exptime)
                # else leave f.mag = nan.
            self.libguide.fiberdata_free(c_fibers)

            self.fibers = fibers
            return fibers
    
    def applyBias(self,image,binning):
        """Apply a bias correction to the image, and return the image and bias level."""
        self.find_bias_level(image,binning=binning)
        self.cmd.diag('text=%s'%qstr('subtracting bias level: %g' % self.imageBias))
        image -= self.imageBias
        return image
        
    def applyDark(self,image,darkFileName,exptime):
        """Dark subtract the current image, and return the result"""
        # Create and process the dark image if this is the first time through, or a new dark exposure
        self.cmd.diag('text=%s'%qstr('Using dark image: %s' % darkFileName))
        if darkFileName != self.currentDarkName and not self.bypassDark:
            self.analyzeDark(darkFileName)
        
        # scale and subtract the dark
        if not self.bypassDark:
            image -= self.processedDark * exptime
        else:
            self.cmd.inform('text=%s'%qstr('Bypassing dark subtraction.'))
        return image
    
    def readProcessedFlat(self, flatFileName, gprobes):
        """
        Read processed flatFileName and return (flat, mask, fibers).
        NOTE, returns a list of fibers the same length as 'gprobes';
        some will have xcen=ycen=NaN; test with fiber.is_fake()
        """
        flatfits = pyfits.open(flatFileName)
        flat = flatfits[0].data
        if self.camera == 'ecamera':
            # TBD: eventually we'll actually compute a mask for the ecam flats
            # for now, just assume all pixels are good.
            mask = np.zeros_like(flat,np.uint8)
            return flat, mask, []
        
        mask = flatfits[1].data
        table = flatfits[6].data
        x = table.field('xcenter')
        y = table.field('ycenter')
        radius = table.field('radius')
        fiberid = table.field('fiberid')
        fibers = []
        if gprobes is not None:
            for xi,yi,ri,fi in zip(x,y,radius,fiberid):
                f = Fiber(fi, xi, yi, ri, 0)
                f.gProbe = gprobes.get(fi)
                fibers.append(f)

        return (flat, mask, fibers)

    def readProcessedDark(self, darkFileName):
        """Read darkFileName and return the data."""
        data = pyfits.open(darkFileName)
        self.darkTemperature = data[0].header['CCDTEMP']
        return data[0].data
    
    def analyzeDark(self, darkFileName, cmd=None, setPoint=None):
        """
        Open a dark file, process it, and save the processd dark as
        self.processedDark, and its temperature as self.darkTemperature.
        """
        if cmd is not None:
            self.cmd = cmd
        if setPoint is not None:
            self.setPoint = setPoint

        darkout = self.getProcessedOutputName(darkFileName)
        if os.path.exists(darkout):
            self.cmd.inform('text=%s'%qstr('Reading processed dark-field from %s' % darkout))
            try:
                self.processedDark = self.readProcessedDark(darkout)
                return
            except:
                self.cmd.warn('text=%s'%qstr('Failed to read processed dark-field from %s; regenerating it.' % darkout))
        
        image,hdr,sat = self._pre_process(darkFileName,binning=self.binning)
        # Fail on bad CCD temp here, since we really need the dark to be at the setPoint.
        if not self._check_ccd_temp(hdr):
            raise GuiderExceptions.BadDarkError
        else:
            self.darkTemperature = hdr['CCDTEMP']
        
        # NOTE: darks are binned.
        # there are very few hot pixels and bulk dark is < 0.02 e/sec at -40C
        
        # Check if it's a good dark
        darkMean = image.mean()
        # TBD: if we use the "outer" bias region (per the ecam data), we need
        # ~40 for a bad dark; for Renbin's preferred valeus ~20 is fine.
        if darkMean > 40:
            self.cmd.warn('text=%s'%qstr('Much more signal in the dark than expected: %6.2f')%darkMean)
            raise GuiderExceptions.BadDarkError
        exptime = hdr['EXPTIME']
        stack = hdr.get('STACK',1)
        exptimen = hdr.get('EXPTIMEN',exptime)
        # We didn't institute the long dark stacking until 2013.06.23,
        # so don't check darks before that date.
        expdate = hdr.get('DATE-OBS').split()[0] # only need date, not time.
        expdate = datetime.datetime.strptime(expdate,'%Y-%m-%d')
        pre_stacking_date = datetime.datetime(2013,6,23)
        if expdate > pre_stacking_date and (stack < 3 or exptimen < 30):
            self.cmd.warn('text=%s'%qstr('Total dark exposure time too short: minimum stack of 3 with total time > 30s .'))
            raise GuiderExceptions.BadDarkError
        if (exptime < 0.5):    #proxy for zero second exposure
           self.cmd.warn('text=%s' % qstr("Dark image less than 0.5 sec"))
           raise GuiderExceptions.BadDarkError
        
        # Convert the dark into a 1-second equivalent exposure.
        # NOTE: we really do want to divide by exptime, not exptimen,
        # because the median-stack behaves like it was exptime long, just with
        # better noise properties.
        image /= exptime
        hdr.update('ORIGEXPT',exptime,'Original, unscaled exposure time.')
        hdr.update('EXPTIME',1.,'dark scaled to 1 second equivalent exposure.')
        
        #Write the dark image
        directory,filename = os.path.split(darkout)
        hdu = pyfits.PrimaryHDU(image,hdr)
        actorFits.writeFits(cmd,hdu,directory,filename,doCompress=True,chmod=0644)
        
        self.processedDark = image
        self.currentDarkName = darkFileName
    #...
    
    def _find_fibers_in_flat(self, image, flatFileName, gprobes, hdr):
        """Identify the fibers in a flat image."""
        # TBD: FIXME -- we KNOW the area (ie, the number of pixels) of the
        # fibers -- we could choose the threshold appropriately.
        # NOTE, that's not always true, we sometimes lose fibers, even
        # acquisition fibers which are pretty big.

        # HACK -- find threshold level inefficiently.
        # This is the threshold level that was used in "gfindfibers".
        i2 = image.copy().ravel()
        I = np.argsort(i2)
        # median
        med = i2[I[len(I)/2]]
        # pk == 99th percentile
        #pk = i2[I[int(len(I)*0.99)]]
        pk = i2[I[int(len(I)*0.998)]]
        #pk = i2[I[int(len(I)*0.998)]]  # change percentile based on N big fibers, N small fibers
        thresh = (med + pk)/2.

        # Threshold image
        T = (image > thresh)
        # Fix nicks in fibers, and fill holes (necessary for the acquisition fibers)
        T = binary_closing(T, iterations=10)

        # Find the background in regions below the threshold.
        background = np.median(image[np.logical_not(T)])
        self.cmd.diag('text=%s'%qstr('Background in flat %s: %g' % (flatFileName, background)))
    
        # TBD: PH, should we make the mask only after the labeled regions that
        # match fibers have been identified?
        # Make the mask a bit smaller than the thresholded fibers.
        mask = binary_erosion(T, iterations=3)

        # Make an annular/ring mask around thresholded fibers
        np.ringmask = binary_dilation(T, iterations=5)

        # Label connected components.
        (fiber_labels,nlabels) = label(T)

        # TBD: FIXME -- the following could be made a bit quicker by using:
        #    objs = find_objects(fiber_albels, nlabels)
        # ==> "objs" is a list of (slice-rows, slice-cols)

        BIN = self.binning

        fibers = []
        self.cmd.diag('text=%s'%qstr('%d components' % (nlabels)))
        for i in range(1, nlabels+1):
            # find pixels labelled as belonging to object i.
            obji = (fiber_labels == i)
            # ri,ci = nonzero(obji)
            npix = obji.sum()
            # center_of_mass returns row,column... swap to x,y
            (yc,xc) = center_of_mass(obji)
            # x,y,radius / BIN because the flat is unbinned pixels, but we want
            # to report in binned pixels.
            # The 0.25 pixel offset makes these centroids agree with gfindstar's
            # pixel coordinate convention.
            self.cmd.diag('text=%s'%qstr('fiber %d (%d,%g,%g,%g)' % (i,npix,xc,yc,np.sqrt(npix/np.pi))))
            fibers.append(Fiber(-1, xc/BIN - 0.25, yc/BIN - 0.25, np.sqrt(npix/np.pi)/BIN, -1, label=i))

        # Match up the fibers with the known probes.

        # Filter fibers that are not within X% of the size of a probe.
        # This should drop small objects (eg, cosmic rays)
        proberads = np.unique([p.radius for p in gprobes.values()])
        keepfibers = []
        for f in fibers:
            dr = abs((f.radius - proberads)/proberads).min()
            # Typically < 7%
            if dr < 0.25:
                keepfibers.append(f)
            else:
                self.cmd.inform('text=%s'%qstr('Rejecting fiber at (%i,%i) in unbinned pixel coords, with radius %g: too far from expected radii %s' %
                          (f.xcen, f.ycen, f.radius, '{' + ', '.join(['%g'%r for r in proberads]) + '}')))
        fibers = keepfibers

        if len(fibers) == 0:
            self.cmd.error('text=%s'%qstr('Failed to find any fibers in guider flat!'))
            raise GuiderExceptions.NoFibersFoundError

        # Find a single x,y offset by testing possibly corresponding
        # pairs of fibers/probes and checking how many fibers are
        # where you expect them to be.
        maxhits = 0
        best = None
        FX = np.array([f.xcen for f in fibers])
        FY = np.array([f.ycen for f in fibers])
        PX = np.array([p.xCenter for p in gprobes.values()])
        PY = np.array([p.yCenter for p in gprobes.values()])
        (H,W) = image.shape
        for i,f in enumerate(fibers):
            for j,(px,py) in enumerate(zip(PX,PY)):
                dx,dy = f.xcen - px, f.ycen - py
                # If this (dx,dy) would put some of the expected probes outside the image, skip it.
                if (min(PX + dx) < 0 or min(PY + dy) < 0 or
                    max(PX + dx) > W/BIN or max(PY + dy) > H/BIN):
                    continue
                # Look for matches...
                fmatch = []
                for k,pp in gprobes.items():
                    # Find the distance from this probe to each fiber...
                    D = np.sqrt(((pp.yCenter + dy) - FY)**2 + ((pp.xCenter + dx) - FX)**2)
                    a = np.argmin(D)
                    # If the nearest one is within range...
                    if D[a] < fibers[a].radius:
                        # Fiber "a" is probe id "k".
                        # Compute the dx,dy implied by this match...
                        mydx = FX[a] - (pp.xCenter + dx)
                        mydy = FY[a] - (pp.yCenter + dy)
                        fmatch.append((a,k, mydx, mydy))

                nhits = len(fmatch)
                if nhits > maxhits:
                    # FIXME -- check for conflicts in 'fmatch' here?
                    maxhits = nhits
                    best = fmatch
                if nhits == len(gprobes):
                    break

        if best is None:
            # How can this happen?  No fibers or no probes...
            self.cmd.warn('text=%s'%qstr("This can't happen?  No matched fibers/probes."))
            raise GuiderExceptions.FlatError

        fmatch = best
        # Check for duplicates...
        fmap = {}
        finvmap = {}
        for fi,probei,dx,dy in fmatch:
            if fi in fmap:
                self.cmd.warn('text=%s'%qstr('Fiber %i wants to match to probe %i and %i.' % (fi, fmap[fi], probei)))
                continue
            if probei in fmap.values():
                self.cmd.warn('text=%s'%qstr('Fiber %i wants to be matched to already-matched probe %i.' % (fi, probei)))
                continue
            fmap[fi] = probei
            finvmap[int(probei)] = fi
        dx = np.mean([x for fi,probei,x,y in fmatch])
        dy = np.mean([y for fi,probei,x,y in fmatch])
        self.cmd.inform('text=%s'%qstr('Matched %i fibers, with dx,dy = (%g,%g)' % (len(fmap), dx, dy)))

        # Record the fiber id, and remember the gprobe...
        for i,f in enumerate(fibers):
            if i not in fmap:
                self.cmd.diag('text=%s'%qstr('probe %d (%d) not matched; skipping...' % (i, f.fiberid)))
                continue
            f.fiberid = fmap[i]
            f.gProbe = gprobes[f.fiberid]
        # Filter out those with no match...
        fibers = [f for i,f in enumerate(fibers) if i in fmap and f.fiberid >= 0]

        # TBD: PH should we create the masks here rather than above?
        #   nonreal fibers have been removed
        #   set non real labels to 0? to modify the mask

        # Reorder fibers by fiberid.
        fibers.sort(key=attrgetter('fiberid'))
        for f in fibers:
            self.cmd.diag('text=%s'%qstr('Fiber id %i at (%.1f, %.1f)' % (f.fiberid, f.xcen-dx, f.ycen-dy)))

        # Create the processed flat image.
        # NOTE: jkp: using float32 to keep the fits header happier.
        flat = np.zeros_like(image).astype(np.float32)
        all_median = np.empty(len(fibers),dtype=np.float32)

        # TBD: Paul Harding
        # For consistent mags of the guide stars independent of cartridges or fibers
        # we need to normalize the flats in a robust way.
        # Assume that acquisition fibers all have a similar throughput
        # The normalization needs to be insensitive to unplugged fibers, dirt on fibers etc. 
        
        # For now try a median over the 14 guide fibers of
        # the median of the unmasked pixel values in each fiber

        #PH to be strictly correct we should use a local bkg around the fiber
 
        #find the median of each fiber
        # have to enumerate, because we could be missing fibers,
        for i,fiber in enumerate(fibers):
            obji = (fiber_labels == fiber.label) # find pixels belonging to object i.
            objflat = (image[obji] - background)   #should calc a local bkg here using ring mask
            flat[obji] += objflat
            # Do not use acquisition fibers, which have higher throughput than guide fibers.
            if fiber.gProbe.fiber_type != 'GUIDE':
                continue
            # so can't use fiberid as an index.
            all_median[i] = np.median(objflat)

        flatscale = np.median(all_median)
        flat /= flatscale

        # Bin down the mask and flat.
        bmask = np.zeros((mask.shape[0]/BIN, mask.shape[1]/BIN), bool)
        for i in range(BIN):
            for j in range(BIN):
                bmask = np.logical_or(bmask, mask[i::BIN,j::BIN] > 0)

        # Use the mask conventions of ipGguide.c...
        mask = np.where(bmask, 0, GuiderImageAnalysis.mask_masked)
        mask = mask.astype(np.uint8)

        flat = bin_image(flat, BIN)
        # Now clip the flat to be within a reasonable range.
        # Have to do it post-binning, otherwise the fiber edges get wonky.
        # This prevents dividing by really big/small numbers when using the flat.
        np.clip(flat, 0.5, 1.5, out=flat)
        
        #int16 for the rotation code in C
        binnedImage = flat.astype(np.int16)

        # TBD: DX?  DY?  Any other stats in here?
        # we don't want the bzero header keyword set.
        # it somehow is ending up in the raw flat files, but not the object files.
        del hdr['BZERO']

        # For writing fiber postage stamps, fill in rotation fields.
        for p in gprobes.values():
            if not hasattr(p, 'rotStar2Sky'):
                p.rotStar2Sky = 90 + p.rotation - p.phi

        hdulist = self._getProcGimgHDUList(hdr, gprobes, fibers, flat, mask, stampImage=binnedImage)
        if hdulist is None:
            self.cmd.warn('text=%s'%qstr('Failed to create processed flat file'))
        
        return hdulist,gprobes

    def _process_ecam_flat(self, image, flatFileName, hdr):
        """
        Process an engineering camera flat.
        Determine the mask, check flux levels, normalize.
        """
        image = bin_image(image,self.binning)

        # not a terrible choice to show the nominal flat level.
        bg = np.median(image)
        # normalize the flat to ~1.
        image /= bg
        hdulist = self._get_basic_hdulist(image, hdr, bg)
        return hdulist

    def analyzeFlat(self, flatFileName, gprobes, cmd=None, setPoint=None):
        """
        Return (flat,mask,fibers): with the processed flat, mask to apply
        to the image and flat to mask everything but the fibers, and a list
        of the fiber number to x,y position.
        
        NOTE: returns a list of fibers the same length as 'gprobes';
        some will have xcen=ycen=NaN; test with fiber.is_fake()
        """
        if cmd is not None:
            self.cmd = cmd
        if setPoint is not None:
            self.setPoint = setPoint

        flatout = self.getProcessedOutputName(flatFileName)
        if self.camera == 'ecamera':
            #TBD: kludge, since we don't gzip unprocessed ecam files, because IRAF.
            flatout = flatout + '.gz'
        directory,filename = os.path.split(flatout)
        
        if os.path.exists(flatout):
            self.cmd.inform('text=%s'%qstr('Reading processed flat-field from %s' % flatout))
            try:
                self.flatImage,self.flatMask,self.flatFibers = self.readProcessedFlat(flatout, gprobes)
                return
            except:
                self.cmd.warn('text=%s'%qstr('Failed to read processed flat-field from %s; regenerating it.' % flatout))

        image,hdr,sat = self._pre_process(flatFileName,binning=1)
        self._check_ccd_temp(hdr)
       
        # NOTE: we are not dark-subtracting the flats,
        # because they are so short and have ~20k counts and are unbinned.

        if self.camera == 'gcamera':
            hdulist,gprobes = self._find_fibers_in_flat(image, flatFileName, gprobes, hdr)
        elif self.camera == 'ecamera':
            hdulist = self._process_ecam_flat(image, flatFileName, hdr)
            gprobes = None
        
        actorFits.writeFits(cmd,hdulist,directory,filename,doCompress=True,chmod=0644)
        # Now read that file we just wrote...
        self.flatImage,self.flatMask,self.flatFibers = self.readProcessedFlat(flatout, gprobes)
        self.currentFlatName = flatFileName
