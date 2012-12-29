import os.path
import numpy
from numpy import *
import ctypes
import pyfits
from scipy.ndimage.morphology import binary_closing, binary_dilation, binary_fill_holes, binary_erosion
from scipy.ndimage.measurements import label, center_of_mass, find_objects

import actorcore.utility.fits as actorFits
from opscore.utility.tback import tback

# A guider fiber and the star image seen through it.
class fiber(object):
    def __init__(self, fibid, xc=numpy.nan, yc=numpy.nan, r=numpy.nan, illr=numpy.nan):
        self.fiberid = fibid
        self.xcen = xc
        self.ycen = yc
        self.radius = r
        self.illrad = illr

        self.xs = numpy.nan
        self.ys = numpy.nan
        self.xyserr = numpy.nan
        self.fwhm = numpy.nan
        self.sky = numpy.nan
        self.skymag = numpy.nan
        self.flux = numpy.nan
        self.mag = numpy.nan

        self.fwhmErr = numpy.nan
        self.dx = numpy.nan
        self.dy = numpy.nan
        self.dRA = numpy.nan
        self.dDec = numpy.nan

        self.gprobe = None

    def __str__(self):
        return ('Fiber id %i: center %g,%g; star %g,%g; radius %g' %
                (self.fiberid, self.xcen, self.ycen, self.xs, self.ys, self.radius))

    def set_fake(self):
        self.xcen = numpy.nan

    def is_fake(self):
        return isnan(self.xcen)


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

def numpy_array_to_REGION(A):
    H, W = A.shape
    ptrtype = ctypes.POINTER(ctypes.c_int16)
    rows = (ptrtype*H)(*[row.ctypes.data_as(ptrtype) for row in A])
    return REGION(H, W, rows) 

def numpy_array_to_MASK(A):
    H, W = A.shape
    ptrtype = ctypes.POINTER(ctypes.c_uint8)
    rows = (ptrtype*H)(*[row.ctypes.data_as(ptrtype) for row in A])
    return MASK(H, W, rows) 

class GuiderImageAnalysis(object):
    '''
    A class for analyzing the images taken by the guiding camera.

    gimg_filename = 'gimg-0400.fits'
    GI = GuiderImageAnalysis(gimg_filename, cmd)
    # cmd must have methods "inform(string)", "warn(string)".
    fibers = GI.findFibers(gState.gprobes)

    # run guider control loop...

    GI.writeFITS(actorState.models, guideCmd, frameInfo, gState.gprobes)

    '''

    # According to
    #  http://sdss3.apo.nmsu.edu/opssoft/guider/ProcessedGuiderImages.html
    #   0 = good
    mask_saturated = 1
    mask_badpixels = 2
    mask_masked    = 4 # ie, outside the guide fiber

    def __init__(self, gimgfn, cmd=None):
        self.gimgfn = gimgfn
        self.outputDir = None
        # set during findFibers():
        self.fibers = None
        self.guiderImage = None
        self.guiderHeader = None
        self.maskImage = None

        if cmd:
            self.informFunc = cmd.inform
            self.warnFunc = cmd.warn
            self.debugFunc = cmd.diag
        else:
            self.informFunc = None
            self.warnFunc = None
            self.debugFunc = None

        # Print debugging?
        self.printDebug = False

        # "Big" (acquisition) fibers are bigger than this pixel
        # radius.  The older SDSS cartridges don't declare (in the
        # gcamFiberInfo.par file) the larger fibers to be ACQUIRE,
        # though they are declared to have radii of 14.1 pixels (vs
        # 8.5 for the GUIDE fibers).  We therefore cut on this radius.
        #
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

        # The photometric zero-point for (g + r)/2 band
        # This was calibrated vs MJD 55246, with AZ~=73 deg, airmass~=1.046
        self.zeropoint = 25.34
        
    def pixels2arcsec(self, pix):
        return pix * self.pixelscale

    def flux2mag(self, flux, exptime, fiber):
        if exptime == 0:
            return -99
        return -2.5 * log10(flux / exptime) + self.zeropoint

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
            # NOTE: can only use the last few columns, as the first few
            # include some  counts that have bled through.
            bias = numpy.median(image[:,(1040/binning):])
        else:
            # find bias = BIAS_PERCENTILE (ipGguide.h) = (100 - 70%)
            ir = image.ravel()
            I = argsort(ir)
            bias = ir[I[int(0.3 * len(ir))]]
        return bias
    #...
    
    def ensureLibraryLoaded(self):
        # Load C library "ipGguide.c"
        path = os.path.expandvars("$GUIDERACTOR_DIR/lib/libguide.so")
        libguide = ctypes.CDLL(path)
        if not libguide:
            self.warn('Failed to load "libguide.so" from %s ($GUIDERACTOR_DIR/lib/libguide.so)' % path)
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

    # For ease of testing...
    def warn(self, s, hasKey=False):
        if self.warnFunc is None:
            print 'warning:', s
        else:
            if not hasKey:
                s = 'text="%s"' % (s)
            self.warnFunc(s)
    def inform(self, s, hasKey=False):
        if self.informFunc is None:
            print 'info:', s
        else:
            if not hasKey:
                s = 'text="%s"' % (s)
            self.informFunc(s)

    def debug(self, s, hasKey=False):
        if self.debugFunc is None:
            print 'debug:', s
        else:
            if not hasKey:
                s = 'text="%s"' % (s)
            self.debugFunc(s)

    def findDarkAndFlat(self, gimgfn, fitsheader):
        ''' findDarkAndFlat(...)

        Returns the filenames containing the dark and flat images for
        the given guider-camera image.

        These are listed in the gimg-####.fits header; this method exists
        to make testing easier, and to make path name magic more explicit.

        DARKFILE= '/data/gcam/55205/gimg-0003.fits'
        FLATFILE= '/data/gcam/55205/gimg-0224.fits'
        '''
        return (fitsheader['DARKFILE'], fitsheader.get('FLATFILE', None))

    def setOutputDir(self, dirnm):
        self.outputDir = dirnm

    def getProcessedOutputName(self, imgfn):
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
            imageHDU.header.update('OBJECT', objectname, '')
            imageHDU.header.update('GCAMSCAL', frameInfo.guideCameraScale, 'guide camera plate scale (mm/pixel)')
            imageHDU.header.update('PLATSCAL', frameInfo.plugPlateScale, 'plug plate scale (mm/degree)')
                        # Do this first, before we need the models.
            guiderCards = self.getGuideloopCards(cmd, frameInfo)
            actorFits.extendHeader(cmd, imageHDU.header, guiderCards)
            self.addPixelWcs(imageHDU.header)

            tccCards = actorFits.tccCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, tccCards)
            mcpCards = actorFits.mcpCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, mcpCards)
            plateCards = actorFits.plateCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, plateCards)
        except Exception, e:
            self.warn('!!!!! failed to fill out primary HDU  !!!!! (%s)' % (e))

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
                ('guideFitRMS',  'gdFRMS',  'RMS of fit to guide star posn, arcsec'),
                ('nguideFitRMS', 'ngdFRMS', 'N stars used for fit RMS'),
                ('decenterRA',    'dcnRA',   'applied user supplied offset in RA, arcsec'),
                ('decenterDec',   'dcnDec',  'applied user supplied offset in Dec, arcsec'),
                ('decenterRot',   'dcnRot',  'applied user supplied rotator offset, arcsec'),
                ('decenterFocus', 'dcnFcus', 'applied user supplied focus offset, um'),
                ('decenterScale', 'dcnScle', 'applied user supplied scale offset, %' ))
                #FIXME PH --- do we change to 1e6 units for scale
        cards = []
        for name, fitsName, comment in defs:
            try:
                val = None
                val = getattr(frameInfo, name)
                print name, val
                if isnan(val):
                    val = -99999.9 # F.ing F.TS
                c = actorFits.makeCard(cmd, fitsName, val, comment)
                cards.append(c)
            except Exception, e:
                self.warn('failed to make guider card %s=%s (%s)' % (name, val, e))
        return cards

    def getStampHDUs(self, fibers, bg, image, mask):
        if len(fibers) == 0:
            return [pyfits.ImageHDU(array([[]]).astype(int16)),
                    pyfits.ImageHDU(array([[]]).astype(uint8))]
        self.ensureLibraryLoaded()
        r = int(ceil(max([f.radius for f in fibers])))
        stamps = []
        maskstamps = []
        for f in fibers:
            xc = int(f.xcen + 0.5)
            yc = int(f.ycen + 0.5)
            rot = -f.gprobe.info.rotStar2Sky
            self.debug("rotating fiber %d at (%d,%d) by %0.1f degrees" % (f.fiberid, xc, yc, rot))
            # Rotate the fiber image...
            stamp = image[yc-r:yc+r+1, xc-r:xc+r+1].astype(int16)
            rstamp = zeros_like(stamp)
            self.libguide.rotate_region(numpy_array_to_REGION(stamp),
                                        numpy_array_to_REGION(rstamp),
                                        rot)
            stamps.append(numpy.flipud(rstamp))
            # Rotate the mask image...
            stamp = mask[yc-r:yc+r+1, xc-r:xc+r+1].astype(uint8)
            #print 'stamp values:', unique(stamp.ravel())
            # Replace zeros by 255 so that after rotate_region (which puts
            # zeroes in the "blank" regions) we can replace them.
            stamp[stamp == 0] = 255
            rstamp = zeros_like(stamp)
            self.libguide.rotate_mask(numpy_array_to_MASK(stamp),
                                      numpy_array_to_MASK(rstamp),
                                      rot)
            # "blank" regions become masked-out pixels.
            rstamp[rstamp == 0] = GuiderImageAnalysis.mask_masked
            # Reinstate the zeroes.
            rstamp[rstamp == 255] = 0
            #print 'rotated stamp values:', unique(rstamp.ravel())
            maskstamps.append(numpy.flipud(rstamp))

        # Stack the stamps into one image.
        stamps = vstack(stamps)
        maskstamps = vstack(maskstamps)
        # Replace zeroes by bg
        stamps[stamps==0] = bg
        # Also replace masked regions by bg.
        # HACKery....
        #stamps[maskstamps > 0] = maximum(bg - 500, 0)
        stamps[maskstamps > 0] = bg
        print 'bg=', bg
        return [pyfits.ImageHDU(stamps), pyfits.ImageHDU(maskstamps)]

    def _getProcGimgHDUList(self, primhdr, gprobes, fibers, image, mask, stampImage=None):
        if stampImage is None:
            stampImage = image

        #bg = median(image[numpy.isfinite(image)])
        bg = median(image[mask == 0])

        imageHDU = pyfits.PrimaryHDU(image, header=primhdr)
        imageHDU.header.update('SDSSFMT', 'GPROC 1 4', 'type major minor version for this file')
        imageHDU.header.update('IMGBACK', bg, 'crude background for entire image. For displays.')
        imageHDU.header.update('OVERSCAN', self.imageBias, 'Bias level of this image.')

        try:
            # List the fibers by fiber id.
            fiberids = gprobes.keys()
            fiberids.sort()
            sfibers = []
            myfibs = dict([(f.fiberid,f) for f in fibers])
            for fid in fiberids:
                if fid in myfibs:
                    sfibers.append(myfibs[fid])
                else:
                    # Put in NaN fake fibers for ones that weren't found.
                    fake = fiber(fid)
                    fake.gprobe = gprobes[fid]
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

            hdulist = pyfits.HDUList()
            hdulist.append(imageHDU)
            hdulist.append(pyfits.ImageHDU(mask))
            hdulist += self.getStampHDUs(smalls, bg, stampImage, mask)
            hdulist += self.getStampHDUs(bigs, bg, stampImage, mask)

            pixunit = 'guidercam pixels (binned)'
            gpfields = [
                        # FITScol,      FITStype, unit
                        ('enabled',        'L',   None),
                        ('flags',          'L',   None),
                        ]
            gpinfofields = [
                           # python attr, FITScol (if diff), FITStype, NIL val, unit
                           ('exists',         None,    'L',   numpy.nan,     None),
                           ('xFocal',         None,    'E',   numpy.nan,     'plate mm'),
                           ('yFocal',         None,    'E',   numpy.nan,     'plate mm'),
                           ('radius',         None,    'E',   numpy.nan,     pixunit),
                           ('xFerruleOffset', None,    'E',   numpy.nan,     None),
                           ('yFerruleOffset', None,    'E',   numpy.nan,     None),
                           ('rotation',       None,    'E',   numpy.nan,     None),
                           ('rotStar2Sky',    None,    'E',   numpy.nan,     None),
                           ('focusOffset',    None,    'E',   numpy.nan,     'micrometers'),
                           ('fiber_type',     None,    'A20', numpy.nan,     None),
                           ('mag',            'ugriz', '5E',  [numpy.nan]*5, 'mag'),
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
                      ('mag',     None,     'E', '(g+r)/2 mag'),
                      ('sky',     None,     'E', 'DN/pixel'),
                      ('skymag',  None,     'E', 'mag/(arcsec^2)'),
                      ('poserr',  'xyserr', 'E', pixunit),
                      ]

            # FIXME -- rotStar2Sky -- should check with "hasattr"...
            cols = []
            for name,fitstype,units in gpfields:
                cols.append(pyfits.Column(name=name, format=fitstype, unit=units,
                                          array=numpy.array([getattr(f.gprobe, name) for f in fibers])))
            for name,fitsname,fitstype,nilval,units in gpinfofields:
                if fitsname is None:
                    fitsname = name
                cols.append(pyfits.Column(name=fitsname, format=fitstype, unit=units,
                                          array=numpy.array([getattr(f.gprobe.info, name, nilval) for f in fibers])))
            for name,atname,fitstype,units in ffields:
                cols.append(pyfits.Column(name=name, format=fitstype, unit=units,
                                          array=numpy.array([getattr(f, atname or name) for f in fibers])))

            cols.append(pyfits.Column(name='stampSize', format='I', array=numpy.array(stampSizes)))
            cols.append(pyfits.Column(name='stampIdx', format='I', array=numpy.array(stampInds)))

            hdulist.append(pyfits.new_table(cols))
        except Exception, e:
            self.warn('could not write proc- guider file: %s' % (e,))
            tback('Narf', e)
            return None
        return hdulist


    def writeFITS(self, models, cmd, frameInfo, gprobes):
        """
        Write a fits file containing the processed results for this exposure.
        """
        if not self.fibers:
            raise Exception('must call findFibers() before writeFITS()')
        fibers = self.fibers
        gimg = self.guiderImage
        hdr = self.guiderHeader

        procpath = self.getProcessedOutputName(self.gimgfn)
        objectname = os.path.splitext(self.gimgfn)[0]

        try:
            hdulist = self._getProcGimgHDUList(hdr, gprobes, fibers, gimg, self.maskImage)
            imageHDU = hdulist[0]
            imageHDU.header.update('SEEING', frameInfo.seeing if numpy.isfinite(frameInfo.seeing) else 0.0,
                           'Estimate of current seeing, arcsec fwhm')
            self.fillPrimaryHDU(cmd, models, imageHDU, frameInfo, objectname)
            hdulist.writeto(procpath, clobber=True)
            dirname, filename = os.path.split(procpath)
            self.inform('file=%s/,%s' % (dirname, filename), hasKey=True)
        except Exception, e:
            cmd.warn('text="failed to write FITS file %s: %r"' % (procpath, e))

    def findFibers(self, gprobes):
        ''' findFibers(gprobes)

        * gprobes is from masterThread.py -- it must be a dict from
          fiber id (int) to Gprobe objects.  The "enabled" and "info"
          fields are used; "info" must have "radius", "xCenter", "yCenter",
          and "rotStar2Sky" fields.

        Analyze a single image from the guider camera, finding fiber
        centers and star centers.

        The fiber centers are found by referring to the associated guider
        flat.  Since this flat is shared by many guider-cam images, the
        results of analyzing it are cached in a "proc-gimg-" file for the flat.

        Returns a list of fibers; also sets several fields in this object.

        The list of fibers contains an entry for each fiber found.
        '''
        assert(self.gimgfn)
        self.ensureLibraryLoaded()

        # Load guider-cam image.
        self.debug('Reading guider-cam image %s' % self.gimgfn)
        p = pyfits.open(self.gimgfn)
        image = p[0].data
        hdr = p[0].header

        #print 'image', image.min(), image.max()
        #print 'sat', self.saturationLevel
        #print image >= self.saturationLevel

        sat = (image.astype(int) >= self.saturationLevel)
        if any(sat):
            self.warn('the exposure has %i saturated pixels' % sum(sat))
        image[sat] = self.saturationReplacement

        # Get dark and flat
        (darkfn, flatfn) = self.findDarkAndFlat(self.gimgfn, hdr)
        cartridgeId = int(hdr['FLATCART'])

        # FIXME -- we currently don't do anything with the dark frame.
        #   we'll need hdr['EXPTIME'] if we do.
        # FIXME -- filter probes here for !exists, !tritium ?
        exptime = hdr.get('EXPTIME', 0)
        print 'Exposure time', exptime

        if hdr['IMAGETYP'] == 'flat':
            flatfn = self.gimgfn
            self.debug('Analysing flat image %s' % (flatfn))
            (flat, mask, fibers) = self.analyzeFlat(flatfn, cartridgeId, gprobes)
            flatoutname = self.getProcessedOutputName(flatfn)
            return fibers

        self.debug('Using flat image %s' % flatfn)
        X = self.analyzeFlat(flatfn, cartridgeId, gprobes)
        if X is None:
            # TBD: jkp: should we ever get here? What kind of crazy result is this?
            return None
        (flat, mask, fibers) = X
        fibers = [f for f in fibers if not f.is_fake()]
        # mask the saturated pixels with the appropriate value.
        mask[sat] |= GuiderImageAnalysis.mask_saturated

        bias = self.find_bias_level(image,binning=self.binning)
        self.inform('subtracting bias level: %g' % bias)
        image -= bias
        self.imageBias = bias

        # Divide by the flat (avoiding NaN where the flat is zero)
        image /= (flat + (flat == 0)*1)
        self.debug('After flattening: image range: %g to %g' % (image.min(), image.max()))

        # NOTE: jkp: post-flat fielding, we need to re-check for saturated pixels and remask them
        sat_flat = (image.astype(int) >= self.saturationLevel)
        image[sat_flat] = self.saturationReplacement
        mask[sat_flat] |= GuiderImageAnalysis.mask_saturated

        # Save the processed image
        self.guiderImage = image/2.0  #PH***quick Kluge needs to be corrected
        self.guiderHeader = hdr
        self.maskImage = mask

        # Prepare image for Gunn C code PSF fitting in "gfindstars"...
        # PH..Kluge to conserve full dynamic range, scale image to fit into signed int in C code.
        # Have to scale up the outputs.  
        img = image[:]/2.0     #
        # The "img16" object must live until after gfindstars() !

        # Zero out parts of the image that are masked out.
        # In this mask convention, 0 = good, >0 is bad.
        # Mark negative pixels
        badpixels = (img < 0) & (mask == 0)
        mask[badpixels] |= GuiderImageAnalysis.mask_badpixels
        # Blank out masked pixels.
        img[mask > 0] = 0
        img16 = img.astype(int16)
        #self.inform('Image (i16) range: %i to %i' % (img16.min(), img16.max()))
        c_image = numpy_array_to_REGION(img16)

        goodfibers = [f for f in fibers if not f.is_fake()]
        #print '%i fibers; %i good.' % (len(fibers), len(goodfibers))
        c_fibers = self.libguide.fiberdata_new(len(goodfibers))

        for i,f in enumerate(goodfibers):
            c_fibers[0].g_fid[i] = f.fiberid
            c_fibers[0].g_xcen[i] = f.xcen
            c_fibers[0].g_ycen[i] = f.ycen
            c_fibers[0].g_fibrad[i] = f.radius
            # FIXME ??
            c_fibers[0].g_illrad[i] = f.radius
        # FIXME --
        #c_fibers.readnoise = ...

        # mode=1: data frame; 0=spot frame
        mode = 1
        res = self.libguide.gfindstars(ctypes.byref(c_image), c_fibers, mode)
        # SH_SUCCESS is this following nutty number...
        if numpy.uint32(res) == numpy.uint32(0x8001c009):
            self.debug('gfindstars returned successfully.')
        else:
            self.warn('gfindstars() returned an error code: %08x (%08x; success=%08x)' % (res, numpy.uint32(res), numpy.uint32(0x8001c009)))

        # pull star positions out of c_fibers, stuff outputs...
        for i,f in enumerate(goodfibers):
            f.xs     = c_fibers[0].g_xs[i]
            f.ys     = c_fibers[0].g_ys[i]
            f.xyserr = c_fibers[0].poserr[i]
            fwhm     = c_fibers[0].fwhm[i]
            if fwhm != FWHM_BAD:
                f.fwhm = self.pixels2arcsec(fwhm)
            # else leave fwhm = nan.
            # FIXME -- figure out good units -- mag/(pix^2)?
            f.sky    = (c_fibers[0].sky[i])*2.0    #correct for image div by 2 for Ggcode 
            f.flux = (c_fibers[0].flux[i])*2.0     
            if f.flux > 0:
                f.mag    = self.flux2mag(f.flux, exptime, f)
            # else leave f.mag = nan.
        self.libguide.fiberdata_free(c_fibers)

        self.fibers = fibers
        return fibers

    @staticmethod
    def binImage(img, BIN):
        binned = zeros((img.shape[0]/BIN, img.shape[1]/BIN), numpy.float32)
        for i in range(BIN):
            for j in range(BIN):
                binned += img[i::BIN,j::BIN]
        binned /= (BIN*BIN)
        return binned
        

    # Returns (flat, mask, fibers)
    # NOTE, returns a list of fibers the same length as 'gprobes';
    # some will have xcen=ycen=NaN; test with fiber.is_fake()
    @staticmethod
    def readProcessedFlat(flatfn, gprobes, stamps=False):
        p = pyfits.open(flatfn)
        flat = p[0].data
        mask = p[1].data
        # small fiber stamps & masks
        # big fiber stamps & masks
        table = p[6].data
        x = table.field('xcenter')
        y = table.field('ycenter')
        radius = table.field('radius')
        fiberid = table.field('fiberid')
        fibers = []
        for xi,yi,ri,fi in zip(x,y,radius,fiberid):
            f = fiber(fi, xi, yi, ri, 0)
            f.gprobe = gprobes.get(fi)
            fibers.append(f)

        if stamps:
            # add stamp,mask fields to fibers.
            smallstamps = p[2].data
            smallmasks  = p[3].data
            bigstamps = p[4].data
            bigmasks  = p[5].data

            # 1=small, 2=big
            stampsize = table.field('stampSize')
            # 0,1,... for each of small and big stamps.
            stampindex = table.field('stampIdx')

            # Assert consistency...
            assert(smallstamps.shape == smallmasks.shape)
            assert(bigstamps.shape == bigmasks.shape)
            SH,SW = smallstamps.shape
            BH,BW = bigstamps.shape
            assert(SW * sum(stampsize == 1) == SH)
            assert(BW * sum(stampsize == 2) == BH)
            assert(all(stampindex[stampsize == 1] >= 0))
            assert(all(stampindex[stampsize == 1] < SH/SW))
            assert(all(stampindex[stampsize == 2] >= 0))
            assert(all(stampindex[stampsize == 2] < BH/BW))

            for i,f in enumerate(fibers):
                # Small
                if stampsize[i] == 1:
                    si = stampindex[i]
                    f.stamp     = smallstamps[si*SW:(si+1)*SW, :]
                    f.stampmask = smallmasks [si*SW:(si+1)*SW, :]
                elif stampsize[i] == 2:
                    si = stampindex[i]
                    f.stamp     = bigstamps[si*BW:(si+1)*BW, :]
                    f.stampmask = bigmasks [si*BW:(si+1)*BW, :]

        return (flat, mask, fibers)

    def analyzeDark(self, darkfn, cartridgeId):
        """Open a dark file, process it, and return the dark data. """
        darkout = self.getProcessedOutputName(darkfn)

        if os.path.exists(darkout):
            self.inform('Reading processed flat-field from %s' % flatout)
            try:
                return GuiderImageAnalysis.readProcessedDark(darkout, gprobes, stamps)
            except:
                self.warn('Failed to read processed dark-field from %s; regenerating it.' % darkout)

        # darks are binned.
        bias = self.find_bias_level(img,binning=self.binning)
        self.inform('subtracting bias level: %g' % bias)
        img -= bias
        self.imageBias = bias
    #...
    
    def analyzeFlat(self, flatfn, cartridgeId, gprobes, stamps=False):
        """
        Return (flat,mask,fibers): with the processed flat, mask to apply
        to the image and flat to mask everything but the fibers, and a list
        of the fiber number to x,y position.
        
        NOTE: returns a list of fibers the same length as 'gprobes';
        some will have xcen=ycen=NaN; test with fiber.is_fake()

        """
        flatout = self.getProcessedOutputName(flatfn)

        if os.path.exists(flatout):
            self.inform('Reading processed flat-field from %s' % flatout)
            try:
                return GuiderImageAnalysis.readProcessedFlat(flatout, gprobes, stamps)
            except:
                self.warn('Failed to read processed flat-field from %s; regenerating it.' % flatout)

        img = pyfits.getdata(flatfn)

        bias = self.find_bias_level(img,binning=1)
        self.inform('subtracting bias level: %g' % bias)
        img -= bias
        self.imageBias = bias

        # FIXME -- we KNOW the area (ie, the number of pixels) of the
        # fibers -- we could choose the threshold appropriately.
        # NOTE, that's not always true, we sometimes lose fibers, even
        # acquisition fibers which are pretty big.

        # HACK -- find threshold level inefficiently.
        # This is the threshold level that was used in "gfindfibers".
        i2 = img.copy().ravel()
        I = argsort(i2)
        # median
        med = i2[I[len(I)/2]]
        # pk == 99th percentile
        #pk = i2[I[int(len(I)*0.99)]]
        pk = i2[I[int(len(I)*0.998)]]
        #pk = i2[I[int(len(I)*0.998)]]  # change percentile based on N big fibers, N small fibers
        thresh = (med + pk)/2.

        # Threshold image
        T = (img > thresh)
        # Fix nicks in fibers, and fill holes (necessary for the acquisition fibers)
        T = binary_closing(T, iterations=10)

        # Find the background in regions below the threshold.
        background = median(img[logical_not(T)].ravel())
        self.debug('Background in flat %s: %g' % (flatfn, background))
    
        # Make the mask a bit smaller than the thresholded fibers.
        mask = binary_erosion(T, iterations=3)

        # Label connected components.
        (fiber_labels,nlabels) = label(T)

        # FIXME -- the following could be made a bit quicker by using:
        #    objs = find_objects(fiber_albels, nlabels)
        # ==> "objs" is a list of (slice-rows, slice-cols)

        BIN = self.binning

        fibers = []
        self.debug('%d components' % (nlabels))
        for i in range(1, nlabels+1):
            # find pixels labelled as belonging to object i.
            obji = (fiber_labels == i)
            # ri,ci = nonzero(obji)
            # print 'fiber', i, 'extent', ci.min(), ci.max(), ri.min(), ri.max()
            npix = sum(obji)
            # center_of_mass returns row,column... swap to x,y
            (yc,xc) = center_of_mass(obji)
            # x,y,radius / BIN because the flat is unbinned pixels, but we want
            # to report in binned pixels.
            # The 0.25 pixel offset makes these centroids agree with gfindstar's
            # pixel coordinate convention.
            self.debug('fiber %d (%d,%g,%g,%g)' % (i,npix,xc,yc,sqrt(npix/pi)))
            fibers.append(fiber(-1, xc/BIN - 0.25, yc/BIN - 0.25, sqrt(npix/pi)/BIN, -1))

        # Match up the fibers with the known probes.

        # Filter fibers that are not within X% of the size of a probe.
        # This should drop small objects (eg, cosmic rays)
        proberads = unique([p.info.radius for p in gprobes.values()])
        keepfibers = []
        for f in fibers:
            dr = abs((f.radius - proberads)/proberads).min()
            # Typically < 7%
            if dr < 0.25:
                keepfibers.append(f)
            else:
                self.inform('Rejecting fiber at (%i,%i) in unbinned pixel coords, with radius %g: too far from expected radii %s' %
                          (f.xcen, f.ycen, f.radius, '{' + ', '.join(['%g'%r for r in proberads]) + '}'))
        fibers = keepfibers

        if len(fibers) == 0:
            self.warn('Failed to find any fibers in guider flat!')
            return None

        # Find a single x,y offset by testing possibly corresponding
        # pairs of fibers/probes and checking how many fibers are
        # where you expect them to be.
        maxhits = 0
        best = None
        FX = array([f.xcen for f in fibers])
        FY = array([f.ycen for f in fibers])
        PX = array([p.info.xCenter for p in gprobes.values()])
        PY = array([p.info.yCenter for p in gprobes.values()])
        (H,W) = img.shape
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
                    D = sqrt(((pp.info.yCenter + dy) - FY)**2 + ((pp.info.xCenter + dx) - FX)**2)
                    a = argmin(D)
                    # If the nearest one is within range...
                    if D[a] < fibers[a].radius:
                        # Fiber "a" is probe id "k".
                        # Compute the dx,dy implied by this match...
                        mydx = FX[a] - (pp.info.xCenter + dx)
                        mydy = FY[a] - (pp.info.yCenter + dy)
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
            self.warn("This can't happen?  No match fibers/probes")
            return None

        fmatch = best
        # Check for duplicates...
        fmap = {}
        finvmap = {}
        for fi,probei,dx,dy in fmatch:
            if fi in fmap:
                self.warn('Fiber %i wants to match to probe %i and %i.' % (fi, fmap[fi], probei))
                continue
            if probei in fmap.values():
                self.warn('Fiber %i wants to be matched to already-matched probe %i.' % (fi, probei))
                continue
            fmap[fi] = probei
            finvmap[int(probei)] = fi
        dx = mean([dx for fi,probei,dx,dy in fmatch])
        dy = mean([dy for fi,probei,dx,dy in fmatch])
        self.inform('Matched %i fibers, with dx,dy = (%g,%g)' % (len(fmap), dx, dy))
        # import pdb; pdb.set_trace()

        # Record the fiber id, and remember the gprobe...
        for i,f in enumerate(fibers):
            if i not in fmap:
                self.debug('probe %d (%d) not matched; skipping...' % (i, f.fiberid))
                continue
            f.fiberid = fmap[i]
            f.gprobe = gprobes[f.fiberid]
        # Filter out those with no match...
        fibers = [f for i,f in enumerate(fibers) if i in fmap and f.fiberid >= 0]

        # Reorder fibers by fiberid.
        I = argsort(array([f.fiberid for f in fibers]))
        fibers = [fibers[i] for i in I]
        for f in fibers:
            self.debug('Fiber id %i at (%.1f, %.1f)' % (f.fiberid, f.xcen-dx, f.ycen-dy))

        # Create the processed flat image.
        # NOTE: jkp: using float32 to keep the fits header happier.
        flat = zeros_like(img).astype(numpy.float32)
        all_mean = numpy.empty(nlabels,dtype=numpy.float32)
        for i in range(1, nlabels+1):
            # find pixels belonging to object i.
            obji = (fiber_labels == i)
            # FIXME -- flat -- normalize by max? median? sum? of this fiber
            # FIXME -- clipping, etc.
            objflat = (img[obji]-background) / median(img[obji]-background)
            flat[obji] += objflat
            all_mean[i-1] = objflat.mean()
            print i, 'objflat (min,max,med):', objflat.min(), objflat.max(), all_mean[i-1]

        # Bin down the mask and flat.
        bmask = zeros((mask.shape[0]/BIN, mask.shape[1]/BIN), bool)
        for i in range(BIN):
            for j in range(BIN):
                bmask = logical_or(bmask, mask[i::BIN,j::BIN] > 0)

        # Use the mask conventions of ipGguide.c...
        mask = where(bmask, 0, GuiderImageAnalysis.mask_masked)
        mask = mask.astype(uint8)

        flat = GuiderImageAnalysis.binImage(flat, BIN)
        # Now clip the flat to be within a reasonable range.
        # Have to do it post-binning, otherwise the fiber edges get wonky.
        # This prevents dividing by really big/small numbers when using the flat.
        numpy.clip(flat, 0.5, 1.5, out=flat)

        #int16 for the rotation code in C
        binimg = GuiderImageAnalysis.binImage(img, BIN).astype(int16)

        # Write output file... copy header cards from input image.
        # DX?  DY?  Any other stats in here?
        hdr = pyfits.getheader(flatfn)
        # we don't want the bzero header keyword set.
        # it somehow is ending up in the raw flat files, but not the object files.
        del hdr['BZERO']
        pixunit = 'binned guider-camera pixels'

        # For writing fiber postage stamps, fill in rotation fields.
        for p in gprobes.values():
            if not hasattr(p.info, 'rotStar2Sky'):
                p.info.rotStar2Sky = 90 + p.info.rotation - p.info.phi

        hdulist = self._getProcGimgHDUList(hdr, gprobes, fibers, flat, mask, stampImage=binimg)
        if hdulist is None:
            self.warn('Failed to create processed flat file')
            return (flat, mask, fibers)

        hdulist.writeto(flatout, clobber=True)
        # Now read that file we just wrote...

        return GuiderImageAnalysis.readProcessedFlat(flatout, gprobes, stamps)
