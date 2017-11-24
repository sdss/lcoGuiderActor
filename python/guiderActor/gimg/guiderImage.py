"""
For processing guider images, indepdendent of commanding the telescope.

Can be used to batch-process old guider images to extract more information
from them in bulk.
"""

import os.path
from operator import attrgetter
import datetime
import re

from astropy.io import fits
import numpy as np

from scipy.ndimage.morphology import binary_closing, binary_dilation, binary_erosion
from scipy.ndimage.measurements import label, center_of_mass
from scipy.ndimage import interpolation

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

        self.reset_star_values()

        self.gProbe = None

        # The label of this fiber in the labeled, masked image (from scipy.ndimage.label)
        self.label = label

    def reset_star_values(self):
        """Reset all values that are measured from a star image."""
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

    def __str__(self):
        return ('Fiber id %i: center %g,%g; star %g,%g; radius %g' %
                (self.fiberid, self.xcen, self.ycen, self.xs, self.ys, self.radius))

    def set_fake(self):
        self.xcen = np.nan

    def is_fake(self):
        return np.isnan(self.xcen)

# Must match ipGguide.h
FWHM_BAD = 99.99


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

    def __init__(self, setPoint, location, bigFiberRadius, zeropoint):
        """
        New GuiderImageAnalysis instances are ready to accept files for processing.

        Args:
            setPoint (float): the current gcamera temperature set point.
            location (str): LCO or APO, so we know what header keywords to use.
        """
        self.setPoint = setPoint
        self.location = location

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

        self.biasFile = None
        self.biasLevel = None

        # amount we let the gcamera temperature vary from the setPoint
        self.deltaTemp = 3.0

        # Clip the flat to lie between these values.
        # Prevents dividing by really big/small numbers when using the flat.
        self.flat_clip = (0.5,1.5)

        # Saturation level.
        #need to make use of full 64k image for bright marvels guide stars.
        # jkp: note that the saturation level is a bit below the
        # full well level of the chip, to prevent column bleedthrough.
        #self.saturationLevel = 0xA000     #41k
        self.saturationLevel = 0xF000    #62,000

        # The value to replace saturated pixels by.
        # self.saturationReplacement = 0xA000
        self.saturationReplacement = 0xF000 #62,000

        # TODO: everything below here should come either from the gimg files,
        # or from the config files. See #2487
        # ===================================================================
        # "Big" (acquisition) fibers are bigger than this pixel
        # radius.  The older SDSS cartridges don't declare (in the
        # gcamFiberInfo.par file) the larger fibers to be ACQUIRE,
        # though they are declared to have radii of 14.1 pixels (vs
        # 8.5 for the GUIDE fibers).  We therefore cut on this radius.
        self.bigFiberRadius = bigFiberRadius

        # The factor by which guider images are binned down.
        # That is, unbinned (flat) images are this factor bigger in
        # each dimension.
        # TODO: this is confusing. Either we should change the name of this
        # attribute or use the binning from the gimgs (JSG).
        self.binning = 2

        # The pixel scale of the guider camera, when binned by "binning"
        # In arcsec/pixel. This is set from the gimg.
        self.pixelscale = np.nan

        self.zeropoint = zeropoint

        # These values are set from the gimg when GuiderImageAnalysis is called
        self.readNoise = np.nan
        self.ccdGain = np.nan

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

    def flux2mag(self, flux):
        """Convert flux to magnitude, using self.exptime."""
        if self.exptime == 0:
            return -99
        return -2.5 * np.log10(flux / self.exptime) + self.zeropoint

    def find_bias_level(self, image, hdr, filename, binning=1):
        """
        Find the bias level of the image.
        Set binning to the number of binned pixels in x and y.

        If there is no overscan region (expected to be 24 extra columns
        unbinned) uses a bias image.
        """
        # The overscan region is an extra 24 columns (12 after binning)
        if image.shape[1] == 1048/binning:
            # subtracting the overall median value should be good enough
            # NOTE: The the "inner" portion of the bias region is most
            # representative of the rest of the chip.
            # TODO: Which is correct here?
            #bias = np.median(image[:,(1024/binning):(1038/binning)])
            # TODO: The upper is what Renbin wanted,
            # TODO: the lower is what makes sense from over-exposed flats and the ecam.
            bias = np.median(image[:,(1039/binning):])
        else:
            # At LCO the gimgs don't have overscan but we have bias images.
            hdr_bias_fn = hdr['BIASFILE']
            if hdr_bias_fn == self.biasFile and self.biasLevel is not None:
                bias = self.biasLevel
                self.cmd.diag('text="Using stored bias level={0:.1f}"'
                              .format(bias))
            else:
                bias_fn = os.path.join(os.path.dirname(filename),
                                       os.path.basename(hdr_bias_fn))
                if not os.path.exists(bias_fn):
                    raise GuiderExceptions.BadBiasError(
                        'Could not find a bias image to use.')

                # TODO: for now I'm using the entire image, but probably it
                # could be trimmed a bit to avoid edge problems. Alternatively,
                # and better, we could use the whole image. (JSG)
                self.cmd.diag('text="reading bias level from {0}"'
                              .format(bias_fn))
                bias = np.median(fits.getdata(bias_fn))
                self.biasLevel = bias
                self.biasFile = bias_fn
                self.cmd.diag('text="bias level={0:.1f}"'.format(bias))

        self.imageBias = bias

        return bias
    #...

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
            darkfile = os.path.join(os.path.dirname(gimgfn),
                                    os.path.basename(darkfile))
        if not os.path.exists(darkfile):
            darkfile = darkfile + '.gz'

        flatfile = fitsheader.get('FLATFILE', None)

        if flatfile is not None:
            if not os.path.exists(flatfile):
                flatfile = os.path.join(os.path.dirname(gimgfn),
                                        os.path.basename(flatfile))
            if not os.path.exists(flatfile):
                flatfile = flatfile + '.gz'

        return darkfile, flatfile

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
        header["CRVAL1%s" % wcsName] = (0, "(output) Column pixel of Reference Pixel")
        header["CRVAL2%s" % wcsName] = (0, "(output) Row pixel of Reference Pixel")
        header["CRPIX1%s" % wcsName] = (0.5, "Column Pixel Coordinate of Reference")
        header["CRPIX2%s" % wcsName] = (0.5, "Row Pixel Coordinate of Reference")
        header["CTYPE1%s" % wcsName] = ("LINEAR", "Type of projection")
        header["CTYPE1%s" % wcsName] = ("LINEAR", "Type of projection")
        header["CUNIT1%s" % wcsName] = ("PIXEL", "Column unit")
        header["CUNIT2%s" % wcsName] = ("PIXEL", "Row unit")

    def fillPrimaryHDU(self, cmd, models, imageHDU, gState, frameInfo, objectname):
        """ Add in all the site and environment keywords. """
        try:
            imageHDU.header['SEEING'] = (
                frameInfo.seeing if np.isfinite(frameInfo.seeing) else 0.0,
                'Estimate of current seeing, arcsec fwhm')
            imageHDU.header['OBJECT'] = (objectname, '')
            imageHDU.header['GCAMSCAL'] = (frameInfo.guideCameraScale, 'guide camera plate scale (mm/pixel)')
            imageHDU.header['PLATSCAL'] = (frameInfo.plugPlateScale, 'plug plate scale (mm/degree)')
                        # Do this first, before we need the models.
            guiderCards = self.getGuideloopCards(cmd, frameInfo, gState)
            actorFits.extendHeader(cmd, imageHDU.header, guiderCards)
            self.addPixelWcs(imageHDU.header)

            plateCards = actorFits.plateCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, plateCards)
            guiderCards = actorFits.guiderCards(models, cmd=cmd)
            actorFits.extendHeader(cmd, imageHDU.header, guiderCards)

            if self.location == 'APO':
                apoCards = actorFits.apoCards(models, cmd=cmd)
                actorFits.extendHeader(cmd, imageHDU.header, apoCards)

        except Exception as e:
            self.cmd.error('text=%s'%qstr('!!!!! failed to fill out primary HDU  !!!!! (%s)' % (e)))
            raise e

    def getGuideloopCards(self, cmd, frameInfo, gState):
        defs = (('guideAxes','guideAxe','guiding on Ra/Dec/Rot?'),
                ('guideFocus','guideFoc','guiding on Focus?'),
                ('guideScale','guideScl','guiding on Scale?'),
                ('dRA', 'DRA', 'measured offset in RA, deg'),
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
                #TODO: FIXME PH --- do we change to 1e6 units for scale
        cards = []
        for name, fitsName, comment in defs:
            try:
                val = None
                val = getattr(frameInfo, name)
                if np.isnan(val):
                    val = -99999.9 # F.ing F.TS
                c = actorFits.makeCard(cmd, fitsName, val, comment)
                cards.append(c)
            except Exception as e:
                self.cmd.warn('text=%s'%qstr('failed to make guider card %s=%s (%s)' % (name, val, e)))

        # LCOHACK: adds PID terms to the header
        for axis in ['raDec', 'rot', 'scale', 'focus']:
            for pid_coeff in ['Kp', 'Ti', 'Td', 'Imax', 'nfilt', 'ncorr']:
                pid_value = getattr(gState.pid[axis], pid_coeff)
                fitsName = '{0}_{1}'.format(axis, pid_coeff)
                if len(fitsName) > 8:
                    fitsName = fitsName[0:8]
                comment = 'PID {0} {1} coefficient'.format(axis, pid_coeff)
                c = actorFits.makeCard(cmd, fitsName, pid_value, comment)
                cards.append(c)

        return cards

    def rotate_one_fiber(self, image, mask, fiber, r):
        """
        Rotate one fiber image stamp (radius r) by the specified angle and
        return the image stamp and mask stamp.
        """
        xc = int(fiber.xcen + 0.5)
        yc = int(fiber.ycen + 0.5)
        rot = -fiber.gProbe.rotStar2Sky
        self.cmd.diag('text=%s'%qstr("rotating fiber %d at (%d,%d) by %0.1f degrees" % (fiber.fiberid, xc, yc, rot)))
        # Rotate the fiber image
        stamp_slice = np.s_[yc-r:yc+r+1, xc-r:xc+r+1]
        # NOTE: we need to rotate by -rot, since interpolation.rotate goes in
        # the opposite direction of the old rotator.
        rstamp = interpolation.rotate(image[stamp_slice],-rot,cval=self.rot_cval,reshape=False)
        # Now rotate the mask, filling with the masked value.
        # Use order=0 for this one, so we don't create "fake" mask values along the edges.
        # This should be maximally pessimistic: only points with no mask in their
        # interpolated value will be unmasked.
        rmaskstamp = interpolation.rotate(mask[stamp_slice],-rot,cval=self.mask_masked,reshape=False,order=0)

        return np.flipud(rstamp), np.flipud(rmaskstamp)

    def getStampHDUs(self, fibers, bg, image, mask):
        """Get the FITS HDUs for the image stamps."""
        if len(fibers) == 0:
            return [fits.ImageHDU(np.array([[]]).astype(np.int16)),
                    fits.ImageHDU(np.array([[]]).astype(np.uint8))]
        r = int(np.ceil(max([f.radius for f in fibers])))
        stamps = []
        maskstamps = []
        for fiber in fibers:
            stamp, maskstamp = self.rotate_one_fiber(image,mask,fiber,r)
            stamps.append(stamp)
            maskstamps.append(maskstamp)

        # Stack the stamps into one image.
        stamps = np.vstack(stamps)
        maskstamps = np.vstack(maskstamps)
        # Replace zeroes by bg
        stamps[stamps==0] = bg
        # Also replace masked regions by bg.
        stamps[maskstamps > 0] = bg
        return [fits.ImageHDU(stamps), fits.ImageHDU(maskstamps)]

    def _get_basic_hdulist(self, image, primhdr, bg):
        """Return an hdulist with the most basic header keywords filled in."""
        imageHDU = fits.PrimaryHDU(image, header=primhdr)
        imageHDU.header['IMGBACK'] = (bg, 'crude background for entire image. For STUI.')
        imageHDU.header['OVERSCAN'] = (self.imageBias, 'Bias level of this image.')
        imageHDU.header['SDSSFMT'] = ('{}PROC 1 4'.format(self.camera[0].upper()), 'guider file version for STUI: v1 has 7 HDUs.')
        hdulist = fits.HDUList()
        hdulist.append(imageHDU)
        return hdulist

    def _get_HDU7_fields(self, fibers, stampSizes, stampInds):
        """Return a list containing the data for HDU7 (the fiber data table)."""

        # jkp TODO: rework this to make it more legible/easily extensible.
        pixunit = 'guidercam pixels (binned)'
        gpinfofields = [
                       # python attr, FITScol (if diff), FITStype, NIL val, unit
                       ('exists',         None,    'L',   np.nan,     None),
                       ('enabled',        None,    'L',   np.nan,     None),
                       ('gprobebits',     None,    'B',   np.nan,     None),
                       ('ra',             None,    'E',   np.nan,     None),
                       ('dec',            None,    'E',   np.nan,     None),
                       ('phi',            None,    'E',   np.nan,     None),
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

        cols = []
        for name,fitsname,fitstype,nilval,units in gpinfofields:
            if fitsname is None:
                fitsname = name
            cols.append(fits.Column(name=fitsname, format=fitstype, unit=units,
                                      array=np.array([getattr(f.gProbe, name, nilval) for f in fibers])))
        for name,atname,fitstype,units in ffields:
            cols.append(fits.Column(name=name, format=fitstype, unit=units,
                                      array=np.array([getattr(f, atname or name) for f in fibers])))

        cols.append(fits.Column(name='stampSize', format='I', array=np.array(stampSizes)))
        cols.append(fits.Column(name='stampIdx', format='I', array=np.array(stampInds)))

        return cols

    def _get_stamp_data(self, gprobes, fibers):
        """Return the data about the different postage stamps, so we can make images."""

        # TODO: this could almost certainly be made much cleaner, possibly
        # TODO: by resetting each gProbe and/or fiber when findStars is called.

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

        # Split into small and big fibers.
        # Also create the columns "stampInds" and "stampSizes"
        # for the output table.
        smalls = []
        bigs = []
        stampInds = []
        stampSizes = []
        for f in sfibers:
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
        return sfibers,bigs,smalls,stampSizes,stampInds

    def _getProcGimgHDUList(self, primhdr, gprobes, fibers, image, mask):
        """Generate an HDU list to be inserted into the proc-file header."""

        bg = np.median(image[mask == 0])
        hdulist = self._get_basic_hdulist(image,primhdr,bg)

        try:
            fibers,bigs,smalls,stampSizes,stampInds = self._get_stamp_data(gprobes,fibers)

            hdulist.append(fits.ImageHDU(mask))
            hdulist += self.getStampHDUs(smalls, bg, image, mask)
            hdulist += self.getStampHDUs(bigs, bg, image, mask)

            hdulist.append(fits.new_table(self._get_HDU7_fields(fibers, stampSizes, stampInds)))

        except Exception as e:
            self.cmd.warn('text=%s'%qstr('could not create header for proc- guider file: %s' % (e,)))
            tback('guiderImage write', e)
            raise e
        return hdulist

    def writeFITS(self, models, cmd, frameInfo, gState, output_verify='warn'):
        """
        Write a fits file containing the processed results for this exposure.
        """

        gprobes = gState.gprobes

        if not self.fibers and self.camera != 'ecamera':
            raise Exception('must call findStars() before writeFITS()')
        image = self.guiderImage
        hdr = self.guiderHeader

        procpath = self.getProcessedOutputName(self.gimgfn)
        objectname = os.path.splitext(self.gimgfn)[0]

        try:
            if self.camera == 'gcamera':
                hdulist = self._getProcGimgHDUList(hdr, gprobes, self.fibers, image, self.maskImage)
                doCompress = True
            elif self.camera == 'ecamera':
                bg = np.median(image)#TODO: this is a poor choice for star-filled ecam images!
                hdulist = self._get_basic_hdulist(image, hdr, bg)
                hdulist.append(fits.ImageHDU(self.maskImage))
                doCompress = False
            imageHDU = hdulist[0]
            self.fillPrimaryHDU(cmd, models, imageHDU, gState, frameInfo, objectname)
            directory,filename = os.path.split(procpath)
            actorFits.writeFits(cmd, hdulist, directory, filename,
                                doCompress=doCompress, chmod=0644,
                                checksum=True, output_verify=output_verify)
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

        # files prior to MJD 56465 have the dark/flat without .gz in the header.
        if not os.path.exists(filename):
            filename += '.gz'
        # Load guider-cam image.
        self.cmd.diag('text=%s'%qstr('Reading guider-cam image %s' % filename))
        try:
            # Read the image as uint (if available), with the BZERO/BSCALE scaling as uint,
            # but then convert to float32 for processing.
            # This solves a problem with over-subtraction of bias level
            # in the LCO files causing uint overflow into large values.
            image,hdr = fits.getdata(filename,0,header=True,uint=True)
            image = np.array(image, dtype=np.float32)
        except IOError:
            raise GuiderExceptions.GuiderError('File not found: %s'%filename)

        # Find saturation level pre-bias.
        sat = (image.astype(int) >= self.saturationLevel)
        nSat = sat.sum()

        image = self.applyBias(image, hdr, filename, binning)
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

        # Sets the pixelScale, gain, and read noise from the header. In case
        # the values are not in the header (for re-reduction of old images),
        # we use old values.

        if 'PIXELSC' in hdr:
            self.pixelscale = hdr['PIXELSC']
        else:
            if self.location == 'APO':
                self.pixelscale = 0.428
            elif self.location == 'LCO':
                self.pixelscale = 0.2834
            self.cmd.warn('text={0}'.format(
                qstr('PIXELSC not found in image header. Using {0:.3f} arcsec.'
                     .format(self.pixelscale))))

        if 'READNOIS' in hdr:
            self.readNoise = hdr['READNOIS']
        else:
            self.readNoise = 10.4
            self.cmd.warn('text={0}'.format(
                qstr('READNOIS not found in image header. Using {0:.1f}.'
                     .format(self.pixelscale))))

        if 'GAIN' in hdr:
            self.ccdGain = hdr['GAIN']
        else:
            self.ccdGain = 1.4
            self.cmd.warn('text={0}'.format(
                qstr('GAIN not found in image header. Using {0:.1f}.'
                     .format(self.pixelscale))))

        return image,hdr,sat
    #...

    def _find_stars_ecam(self, image, mask):
        """Find the stars in a processed ecamera image."""
        ccdInfo = PyGuide.CCDInfo(self.imageBias,self.readNoise,self.ccdGain)
        # set a high threshold, since we only care about obviously bright stars.
        # TODO: this threshold could be a configurable parameter, so the observer
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

    def _set_fiber_star(self, fiber, stars, image, mask, stampFrameCoords):
        """Fill in the fiber values for this star, safely."""

        try:
            # findStars returns in order of brightness, so take the first one.
            star = stars[0]
        except IndexError:
            self.cmd.warn('text=%s'%qstr('No star found via PyGuide.findStars for fiber %d.'%(fiber.fiberid)))
            return

        if not star.isOK:
            self.cmd.warn('text=%s'%qstr('Danger: this should never happen! PyGuide.findStars failed on fiber %d with: %s.'%(fiber.fiberid,star.msgStr)))
            return

        # shift the centers back into the data frame.
        fiber.xs = stampFrameCoords[0] + star.xyCtr[0]
        fiber.ys = stampFrameCoords[1] + star.xyCtr[1]
        fiber.xyserr = np.hypot(*star.xyErr)

        try:
            shape = PyGuide.StarShape.starShape(image, mask, star.xyCtr, 100)
            if not shape.isOK:
                self.cmd.warn('text=%s'%qstr('PyGuide.starShape failed on fiber %d with: %s.'%(fiber.fiberid,star.msgStr)))
                return

            fiber.fwhm = self.pixels2arcsec(shape.fwhm)
            fiber.sky = shape.bkgnd
            # counts is the total, so we need to remove the integrated background.
            fiber.flux = star.counts - shape.bkgnd*star.pix
            if fiber.flux > 0:
                fiber.mag = self.flux2mag(fiber.flux)
        except Exception as e:
            self.cmd.warn('text=%s'%qstr('PyGuide.starShape failed on fiber %d with Exception: %s.'%(fiber.fiberid,e)))

    def _find_stars_gcam(self, image, mask, fibers):
        """Find the stars in a processed gcamera image."""
        ccdInfo = PyGuide.CCDInfo(self.imageBias,self.readNoise,self.ccdGain)
        # findStars wants False for regions of good data
        good_mask = (mask & self.mask_masked) != 0
        saturated = (mask & self.mask_saturated) != 0

        for fiber in fibers:
            fiber.reset_star_values()
            xc = fiber.xcen
            yc = fiber.ycen
            r = fiber.radius

            # Calculates the stamp slice.
            y0, y1 = np.rint(yc - r).astype(np.int), np.rint(yc + r + 1).astype(np.int)
            x0, x1 = np.rint(xc - r).astype(np.int), np.rint(xc + r + 1).astype(np.int)

            # This defines the origin of the coordinates of the stamp with
            # respect to the general data frame. We'll use it to convert
            # back from stamp centroid coordinates to data frame ones.
            stampFrameCoords = (x0, y0)

            stamp = np.s_[y0:y1, x0:x1]

            tritium_offset = []  # Stores flexure measurements using the tritium source(s)

            # use a medium threshold, since the stars might not be that bright when acquiring
            try:
                stars = PyGuide.findStars(image[stamp], good_mask[stamp], saturated[stamp], ccdInfo, thresh=2)[0]
            except Exception as e:
                self.cmd.warn('text=%s'%qstr('PyGuide.findStars failed on fiber %d with: %s.'%(fiber.fiberid,e)))
            else:

                # For not trititum sources we call _set_fiber_star to set xstar/ystar in the
                # fibre to the values measured by PyGuide. For trititum sources, we compare it
                # with the value in the fibre and determine an offset
                if not fiber.gprobe.tritium:
                    self._set_fiber_star(fiber, stars, image[stamp], good_mask[stamp],
                                         stampFrameCoords)
                else:
                    if len(stars) == 0:
                        self.cmd.warn('text="PyGuide.findStars did not '
                                      'find a source for fiber {:d}."'.format(fiber.fiberid))
                        continue

                    xs = stampFrameCoords[0] + stars[0].xyCtr[0]
                    ys = stampFrameCoords[1] + stars[0].xyCtr[1]

                    tritium_offset.append([xs - fiber.xCenter, ys - fiber.yCenter])

        if len(tritium_offset) > 0:
            tritium_offset = tritium_offset.mean(axis=0)
            if np.any(np.abs(tritium_offset) > 1):
                self.cmd.warn('text="tritium sources detected an offset of '
                              '({:.1f}, {:.1f}) pix."'.format(tritium_offset[0],
                                                              tritium_offset[1]))

                for fiber in fibers:
                    fiber.xCenter += tritium_offset[0]
                    fiber.yCenter += tritium_offset[1]

        self.fibers = fibers

    def findStars(self, gprobes):
        """
        Identify the centers of the stars in the fibers.
        Assumes self.gimgfn is set to the correct file name.

        Analyze a single image from the guider camera, finding fiber
        centers and star centers.

        The fiber centers are found by referring to the associated guider
        flat.  Since this flat is shared by many guider-cam images, the
        results of analyzing it are cached in a "proc-gimg-" file for the flat.

        Return a list of fibers; also sets several fields in this object.
        The list of fibers contains an entry for each fiber found.

        Args:
            gprobes is the dictionary of GuideProbes from GuiderState.

        """
        image,hdr,sat = self._pre_process(self.gimgfn,binning=self.binning)

        self.exptime = hdr.get('EXPTIME', 0)

        (darkFileName, flatFileName) = self.findDarkAndFlat(self.gimgfn, hdr)
        image = self.applyDark(image, darkFileName)

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
        # mask the saturated pixels with the appropriate value.
        # NOTE: Why to do this now if we are going to repeat it later?
        mask = self.flatMask.copy()
        mask[sat] |= self.mask_saturated

        # Divide by the flat (avoiding NaN where the flat is zero)
        image /= (self.flatImage + (self.flatImage == 0)*1)
        self.cmd.diag('text=%s'%qstr('After flattening: image range: %g to %g' % (image.min(), image.max())))

        # NOTE: jkp: post-flat fielding, we need to re-check for saturated pixels and remask them
        sat_flat = (image.astype(int) >= self.saturationLevel)
        image[sat_flat] = self.saturationReplacement
        mask[sat_flat] |= self.mask_saturated

        # Save the processed image
        self.guiderImage = image
        self.guiderHeader = hdr
        self.maskImage = mask

        # Zero out parts of the image that are masked out.
        # In this mask convention, 0 = good, >0 is bad.
        # Mark negative pixels
        badpixels = (image < 0) & (mask == 0)
        mask[badpixels] |= self.mask_badpixels

        # If the percentage of bad pixels is > 5%, there may have been a
        # flip in the lenses. We issue a warning.
        if np.sum(badpixels) / float(np.size(self.guiderImage)) > 0.05:
            self.cmd.warn('text="more than 5% of the pixels are marked bad. '
                          'Maybe take another guider flat?"')

        # Blank out masked pixels.
        image[mask > 0] = 0

        # images use a different rotation fill constant than flats.
        self.rot_cval = 0

        if self.camera == 'ecamera':
            # mask the overscan too, since we're keeping it around for monitoring.
            mask[:,(1039/self.binning):] |= self.mask_saturated
            ecam_mask = mask == self.mask_saturated
            star,shape = self._find_stars_ecam(image, ecam_mask)
            if star is None:
                self.cmd.warn('text="ecam_star=None in exposure %d!"'%(self.frameNo))
                return []
            if not shape.isOK:
                self.cmd.warn('text="Could not determine star shape: %s"'%shape.msgStr)
            self.cmd.inform('ecam_star={},{:.2f},{:.2f},{:.3f},{:.2f},{:.1f}'.format(self.frameNo,star.xyCtr[0],star.xyCtr[1],shape.fwhm,shape.bkgnd,shape.ampl))
            return [] # no fibers to return
        else:
            fibers = [f for f in self.flatFibers if not f.is_fake()]
            self._find_stars_gcam(image, mask, fibers)
            return self.fibers

    def applyBias(self, image, hdr, filename, binning):
        """Apply a bias correction to the image, and return the image and bias level."""
        self.find_bias_level(image, hdr, filename, binning=binning)
        self.cmd.diag('text=%s'%qstr('subtracting bias level: %g' % self.imageBias))
        image = image - self.imageBias
        return image

    def applyDark(self,image,darkFileName):
        """Dark subtract the current image, and return the result"""
        # Create and process the dark image if this is the first time through, or a new dark exposure
        self.cmd.diag('text=%s'%qstr('Using dark image: %s' % darkFileName))
        if darkFileName != self.currentDarkName and not self.bypassDark:
            self.analyzeDark(darkFileName)

        # scale and subtract the dark
        if not self.bypassDark:
            image -= self.processedDark * self.exptime
        else:
            self.cmd.inform('text=%s'%qstr('Bypassing dark subtraction.'))
        return image

    def readProcessedFlat(self, flatFileName, gprobes):
        """
        Read processed flatFileName and return (flat, mask, fibers).
        NOTE, returns a list of fibers the same length as 'gprobes';
        some will have xcen=ycen=NaN; test with fiber.is_fake()
        """
        flatfits = fits.open(flatFileName)
        flat = flatfits[0].data
        if self.camera == 'ecamera':
            # TODO: eventually we'll actually compute a mask for the ecam flats
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
        data = fits.open(darkFileName)
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
            raise GuiderExceptions.BadDarkError('CCD temp significantly different from setPoint.')
        else:
            self.darkTemperature = hdr['CCDTEMP']

        # NOTE: darks are binned.
        # there are very few hot pixels and bulk dark is < 0.02 e/sec at -40C

        # First, check that the dark is stacked and long enough to have some signal.
        exptime = hdr['EXPTIME']
        stack = hdr.get('STACK', 1)
        exptimen = hdr.get('EXPTIMEN', exptime)
        # We didn't institute the long dark stacking until 2013.06.23,
        # so don't check darks before that date.
        expdate = hdr.get('DATE-OBS').split()[0]  # only need date, not time.
        expdate = datetime.datetime.strptime(expdate, '%Y-%m-%d')
        pre_stacking_date = datetime.datetime(2013, 6, 23)
        if expdate > pre_stacking_date and (stack < 3 or exptimen < 30):
            errMsg = 'Total dark exposure time too short: minimum stack of 3 with total time > 30s .'
            self.cmd.warn('text=%s' % qstr(errMsg))
            raise GuiderExceptions.BadDarkError(errMsg)

        # Check if it's a good dark
        darkMean = image.mean()
        # TODO: if we use the "outer" bias region (per the ecam data), we need
        # ~40 for a bad dark; for Renbin's preferred values ~20 is fine.
        if darkMean > 40:
            errMsg = 'Much more signal in the dark than expected: %6.2f'%darkMean
            self.cmd.warn('text=%s'%qstr(errMsg))
            raise GuiderExceptions.BadDarkError(errMsg)

        # Convert the dark into a 1-second equivalent exposure.
        # NOTE: we really do want to divide by exptime, not exptimen,
        # because the median-stack behaves like it was exptime long, just with
        # better noise properties.
        image /= exptime
        hdr['ORIGEXPT'] = (exptime, 'Original, unscaled exposure time.')
        hdr['EXPTIME'] = (1., 'dark scaled to 1 second equivalent exposure.')

        #Write the dark image
        directory,filename = os.path.split(darkout)
        hdu = fits.PrimaryHDU(image,hdr)
        actorFits.writeFits(cmd,hdu,directory,filename,doCompress=True,chmod=0644)

        self.processedDark = image
        self.currentDarkName = darkFileName
    #...

    def _find_fibers_in_flat(self, image, flatFileName, gprobes, hdr):
        """Identify the fibers in a flat image."""
        # TODO: FIXME -- we KNOW the area (ie, the number of pixels) of the
        # fibers -- we could choose the threshold appropriately.
        # NOTE, that's not always true, we sometimes lose fibers, even
        # acquisition fibers which are pretty big.

        # Finds threshold level.
        median = np.median(image)
        pk = np.percentile(image, 99.8, interpolation='higher')  # We could use 'linear'
        thresh = (median + pk) / 2.

        # Threshold image
        T = (image > thresh)
        # Fix nicks in fibers, and fill holes (necessary for the acquisition fibers)
        T = binary_closing(T, iterations=10)

        # Find the background in regions below the threshold.
        background = np.median(image[np.logical_not(T)])
        self.cmd.diag('text=%s'%qstr('Background in flat %s: %g' % (flatFileName, background)))

        # TODO: PH, should we make the mask only after the labeled regions that
        # match fibers have been identified?
        # Make the mask a bit smaller than the thresholded fibers.
        mask = binary_erosion(T, iterations=3)

        # Checks that the fibres have a good average number of counts
        image_masked = image[mask]
        fibers_mean = np.median(image_masked)
        if fibers_mean < 1e4:
            self.cmd.error('text="average number of counts in fibres < 1e4."')
            raise GuiderExceptions.FlatError

        self.cmd.warn('text="average number of counts in fibres = {0:.1f}"'.format(fibers_mean))

        # Make an annular/ring mask around thresholded fibers
        # NOTE: where is this used? Why is it assigned to a Numpy attribute.
        np.ringmask = binary_dilation(T, iterations=5)

        # Label connected components.
        (fiber_labels,nlabels) = label(T)

        # TODO: FIXME -- the following could be made a bit quicker by using:
        #    objs = find_objects(fiber_albels, nlabels)
        # ==> "objs" is a list of (slice-rows, slice-cols)

        BIN = self.binning

        fibers = []
        self.cmd.diag('text=%s'%qstr('%d components' % (nlabels)))
        for i in range(1, nlabels + 1):
            # find pixels labelled as belonging to object i.
            obji = (fiber_labels == i)
            # ri,ci = nonzero(obji)
            npix = obji.sum()

            # center_of_mass returns row,column... swap to x,y
            (yc, xc) = center_of_mass(obji)
            radius = np.sqrt(npix / np.pi)

            # x,y,radius / BIN because the flat is unbinned pixels, but we want
            # to report in binned pixels.
            # The 0.25 pixel offset makes these centroids agree with gfindstar's
            # pixel coordinate convention.

            self.cmd.diag('text={0}'.format(qstr('fiber %d (%d,%g,%g,%g)' %
                                                 (i, npix, xc, yc, radius))))

            # Creates the Fiber objects
            xx = xc / BIN - 0.25
            yy = yc / BIN - 0.25
            rr = np.sqrt(npix / np.pi) / BIN  # Radius of the fiber
            fibers.append(Fiber(-1, xx, yy, rr, -1, label=i))

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
        PX = np.array([p.xCenter for p in gprobes.values() if not gprobe.tritium])
        PY = np.array([p.yCenter for p in gprobes.values() if not gprobe.tritium])
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

        # TODO: PH should we create the masks here rather than above?
        #   nonreal fibers have been removed
        #   set non real labels to 0? to modify the mask

        # Reorder fibers by fiberid.
        fibers.sort(key=attrgetter('fiberid'))
        for f in fibers:
            self.cmd.diag('text=%s'%qstr('Fiber id %i at (%.1f, %.1f)' % (f.fiberid, f.xcen-dx, f.ycen-dy)))

        good_mask = (mask & self.mask_masked) != 0
        saturated = (mask & self.mask_saturated) != 0

        # Now we deal with tritium/LED sources. We use PyGuide since these
        # sources are point-like.
        for gprobe in gprobes:

            if not gprobe.tritium:
                continue

            ccdInfo = PyGuide.CCDInfo(self.imageBias, self.readNoise, self.ccdGain)

            gprobe_xc = gprobe.xCenter
            gprobe_yc = gprobe.xCenter
            gprobe_r = gprobe.radius

            # Calculates a stamp slice around the gprobe centre
            y0 np.rint(gprobe_yc - gprobe_r).astype(np.int)
            y1 = np.rint(gprobe_yc + gprobe_r + 1).astype(np.int)

            x0 = np.rint(gprobe_xc - gprobe_r).astype(np.int)
            x1 = np.rint(gprobe_xc + gprobe_r + 1).astype(np.int)

            stampFrameCoords = (x0, y0)
            stamp_slice = np.s_[y0:y1, x0:x1]

            stamp = image[stamp_slice]

            # For now we use empty masks. Maybe this will need to change.
            mask = saturated = np.zeros(image.shape)

            try:
                stars = PyGuide.findStars(stamp, mask, saturated, ccdInfo, thresh=2)[0]
                star = stars[0]
                assert star.isOK
            except Exception as ee:
                self.cmd.warn('text="failed to find centroid for tritium gprobe {}: {}"'
                              .format(gprobe.id, ee))
            else:
                fiber_x = stampFrameCoords[0] + star.xyCtr[0]
                fiber_y = stampFrameCoords[1] + star.xyCtr[1]
                fiber = Fiber(gprobe.id, fiber_x, fiber_y, 0, -1, label=-1)
                fiber.gProbe = gprobe

                fibers.append(fiber)

                self.cmd.inform('text="found tritium source for '
                                'gprobe {} at ({:.1f}, {:.1f})"'.format(gprobe.id,
                                                                        fiber_x,
                                                                        fiber_y))

        # Create the processed flat image.
        # NOTE: jkp: using float32 to keep the fits header happier.
        flat = np.zeros_like(image).astype(np.float32)
        all_median = np.empty(len(fibers),dtype=np.float32)

        # TODO: Paul Harding
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
        mask = np.where(bmask, 0, self.mask_masked)
        mask = mask.astype(np.uint8)

        flat = bin_image(flat, BIN)
        # Now clip the flat to be within a reasonable range.
        # Have to do it post-binning, otherwise the fiber edges get wonky.
        # This prevents dividing by really big/small numbers when using the flat.
        np.clip(flat, self.flat_clip[0], self.flat_clip[1], out=flat)

        # images use a different rotation fill constant than flats.
        self.rot_cval = self.flat_clip[0]

        # TODO: DX?  DY?  Any other stats in here?
        # we don't want the bzero header keyword set.
        # it somehow is ending up in the raw flat files, but not the object files.
        if 'BZERO' in hdr:
            del hdr['BZERO']

        hdulist = self._getProcGimgHDUList(hdr, gprobes, fibers, flat, mask)
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

    def create_fibers_fake_gprobe(self, image, flatFileName, gprobes, hdr):
        """Processes a flat that has a single synthetic gprobe."""

        self.rot_cval = self.flat_clip[0]

        binning = self.binning
        mask = np.zeros((image.shape[0] / self.binning,
                         image.shape[1] / self.binning))

        fibers = []
        # For each fake gprobe, we create a fiber matching position and radius.
        for gprobeID in gprobes:

            gprobe = gprobes[gprobeID]

            xCenter = gprobe.xCenter
            yCenter = gprobe.yCenter
            radius = gprobe.radius

            fiber = Fiber(gprobe.id, xCenter, yCenter, radius, -1,
                          label=gprobe.id)
            fiber.gProbe = gprobe
            fibers.append(fiber)

            # Adds the fiber to the mask
            yGrid, xGrid = np.ogrid[-yCenter:mask.shape[0] - yCenter,
                                    -xCenter:mask.shape[1] - xCenter]
            mask_fiber = xGrid ** 2 + yGrid ** 2 <= radius ** 2
            mask[mask_fiber] = True

            # We print the unbinned parameters of the fake fibre
            xx = (xCenter + 0.25) * binning
            yy = (yCenter + 0.25) * binning
            rr = radius * binning

            self.cmd.diag(
                'text={0}'.format(qstr('fiber {0:d} ({1:g},{2:g},{3:g})'
                                       .format(gprobe.id, xx, yy, rr))))

        # For now, let's just make the flat the original image, binned to
        # the desired binning.
        # LCOHACK: apply background and flatscale here
        flat = bin_image(image, self.binning)

        mask = np.where(mask, 0, self.mask_masked)
        mask = mask.astype(np.uint8)

        hdulist = self._getProcGimgHDUList(hdr, gprobes, fibers, flat, mask)

        if hdulist is None:
            self.cmd.warn('text=%s' %
                          qstr('Failed to create processed flat file'))

        return hdulist, gprobes

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
        if self.camera == 'ecamera' and '.gz' not in flatout:
            #TODO: kludge, since we currently don't gzip unprocessed ecam files, because IRAF.
            flatout = flatout + '.gz'
        directory,filename = os.path.split(flatout)

        if os.path.exists(flatout):
            self.cmd.inform('text=%s'%qstr('Reading processed flat-field from %s' % flatout))
            try:
                self.flatImage,self.flatMask,self.flatFibers = self.readProcessedFlat(flatout, gprobes)
                self.currentFlatName = flatFileName
                return
            except:
                self.cmd.warn('text=%s'%qstr('Failed to read processed flat-field from %s; regenerating it.' % flatout))

        image,hdr,sat = self._pre_process(flatFileName,binning=1)
        self._check_ccd_temp(hdr)

        # NOTE: we are not dark-subtracting the flats,
        # because they are so short and have ~20k counts and are unbinned.

        if self.camera == 'gcamera':
            # If cart is 99, we don't want to detect fibres from the image. We
            # create a fibre that matches the gprobe.
            if hdr['FLATCART'] == 99:
                hdulist, gprobes = self.create_fibers_fake_gprobe(image, flatFileName,
                                                                  gprobes, hdr)
            else:
                hdulist, gprobes = self._find_fibers_in_flat(image, flatFileName,
                                                             gprobes, hdr)

        elif self.camera == 'ecamera':
            hdulist = self._process_ecam_flat(image, flatFileName, hdr)
            gprobes = None

        actorFits.writeFits(cmd,hdulist,directory,filename,doCompress=True,chmod=0644)
        # Now read that file we just wrote...
        self.flatImage,self.flatMask,self.flatFibers = self.readProcessedFlat(flatout, gprobes)
        self.currentFlatName = flatFileName
