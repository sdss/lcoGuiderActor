import ctypes
import numpy as np
import os
import pyfits

import YPF

import pdb

ACQFIBER = 1
GUIDEFIBER = 2
POSFIBER = 3

class Cartridge(object):
    # Per some SDSS-1 image, extended with some dummies.
    fiberInfoType = np.dtype([('gp_xcen','float32'),
                              ('gp_ycen','float32'),
                              ('gp_rot','float32'),
                              ('gp_radius','float32'),
                              ('gp_xFerruleOffset','float32'),
                              ('gp_yFerruleOffset','float32'),
                              ('gp_focusOffset','float32'),
                              ('gp_exists', 'int8'),
                              ('gp_fiberType', 'int8'),
                              
                              ('fiberCenterX', 'float32'),
                              ('fiberCenterY', 'float32'),
                              ('starCenterX', 'float32'),
                              ('starCenterY', 'float32'),

                              ('starCenterErrMajor', 'float32'),
                              ('starCenterErrMinor', 'float32'),
                              ('starCenterErrCross', 'float32'),
                              
                              ])


    def __init__(self, cartID):
        self.cartID = cartID
        self._loadInfo()

    def _loadInfo(self):
        prodDir = os.environ['GCAMERA_DIR']
        gprobes = YPF.YPFreadone(os.path.join(prodDir, 'etc/gcamFiberInfo-55027.par'),
                                 'GPROBE')
        w = np.where(gprobes.cartridgeId == self.cartID)
        gprobes = gprobes[w]
        nprobes = len(gprobes)

        fiberInfo = np.zeros(nprobes, dtype=self.fiberInfoType)

        # Add a fiber type field -- this is bad.
        bigfibers = np.where(gprobes.radius > 20)
        fiberInfo['gp_fiberType'][bigfibers] = ACQFIBER

        smallfibers = np.where((gprobes.radius < 20) & (gprobes.radius > 2))
        fiberInfo['gp_fiberType'][smallfibers] = GUIDEFIBER

        posfiber = np.where(gprobes.radius < 2)
        fiberInfo['gp_fiberType'][posfiber] = POSFIBER
        
        realfibers = np.where(gprobes.radius > 2)
        nfibers = len(realfibers[0])
        
        for n in ('exists', 'xcen', 'ycen', 'radius', 'rot',
                  'xFerruleOffset', 'yFerruleOffset', 'focusOffset'):
            fiberInfo['gp_'+n] = gprobes[n]

        # Some number are unmeasured. Provide signal...
        fiberInfo['gp_xFerruleOffset'] = 0
        fiberInfo['gp_yFerruleOffset'] = 0
        fiberInfo['gp_rot'] = 0
        fiberInfo['gp_focusOffset'] = 0

        fiberInfo['gp_xcen'] /= 2
        fiberInfo['gp_ycen'] /= 2

        # pdb.set_trace()
        
        fiberInfo['gp_xFerruleOffset'][realfibers] = np.random.uniform(low=-2.0, high=2.0,
                                                                       size=(nfibers)).astype('f4')
        fiberInfo['gp_yFerruleOffset'][realfibers] = np.random.uniform(low=-3.0, high=1.0,
                                                                       size=(nfibers)).astype('f4')
        fiberInfo['gp_rot'][realfibers] = np.random.uniform(low=-90.0, high=90.0,
                                                            size=(nfibers)).astype('f4')
        fiberInfo['gp_focusOffset'][0:4] = 100
        fiberInfo['gp_focusOffset'][4:8] = -100

        self.fiberInfo = fiberInfo
        
    def getInfo(self):
        return self.fiberInfo
        
class PullTest(object):
    def __init__(self, filename, cartID, plugmap):
        """ """

        # pdb.set_trace()
        self.rawImage = pyfits.getdata(filename).astype('u2')

        self.makeFakeRawMetadata(filename)
        cart = Cartridge(cartID)
        self.fiberInfo = cart.getInfo()
        self.getPlugmap(plugmap)
        # Lie, cheat, steal
        w = np.where(self.fiberInfo['gp_fiberType'] == ACQFIBER)
        self.bigFiberSize = self.fiberInfo[w[0][0]]['gp_radius'] * 2
        w = np.where(self.fiberInfo['gp_fiberType'] == GUIDEFIBER)
        self.smallFiberSize = self.fiberInfo[w[0][0]]['gp_radius'] * 2
        
        self.gutils = None
        try:
            self._setupForCcalls()
        except Exception, e:
            print "failed to load external C library: ", e

    def _setupForCcalls(self):
        """ Load the external C library, and define the prototypes. """
        self.gutils = ctypes.CDLL(os.path.expandvars('$GCAMERA_DIR/lib/gimg.so'))

        # void grot(float, short **, short **, int, int);
        self.gutils.grot.restype = None
        self.gutils.grot.argtypes = [ctypes.c_float,
                                     ctypes.POINTER(ctypes.POINTER(ctypes.c_short)),
                                     ctypes.POINTER(ctypes.POINTER(ctypes.c_short)),
                                     ctypes.c_int, ctypes.c_int]
        
        # void maskrot(float, char **, char **, int, int);
        self.gutils.maskrot.restype = None
        self.gutils.maskrot.argtypes = [ctypes.c_float,
                                        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
                                        ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
                                        ctypes.c_int, ctypes.c_int]

    def makeFakeRawMetadata(self, filename=None):
        """ Create the basic exposure metadata. """

        fd = pyfits.getheader(filename)
        md = {}
        # Copy most of the input FITS cards
        for c in ('exptime', 'date-obs', 'filename', 'ccdtemp',
                  'binx', 'biny', 'begx', 'begy', 'fullx', 'fully'):
            if c in fd:
                md[c] = fd[c]
            else:
                print "WARNING: card %s is missing" % (c)

        # Overwrite the timesys card with what we will actually use.
        md['timesys'] = 'TAI'

        self.metadata = md
    
    def makeRawMetadata(self):
        """ Create the basic exposure metadata. """
        raise UnimplementedError()
    
    def findCenters(self):
        """ Dummy centroider: copies through the seed positions. """

        self.fiberInfo['fiberCenterX'] = self.fiberInfo['gp_xcen'] + 0.1
        self.fiberInfo['fiberCenterY'] = self.fiberInfo['gp_ycen'] - 0.2

    def calcOffsets(self):
        realfibers = np.where(self.fiberInfo['gp_fiberType'] != POSFIBER)
        nfibers = len(realfibers[0])
        self.fiberInfo['starCenterX'][realfibers] = (self.fiberInfo['fiberCenterX'][realfibers] + 
                                                     np.random.uniform(nfibers) * -0.2)
        self.fiberInfo['starCenterY'][realfibers] = (self.fiberInfo['fiberCenterY'][realfibers] +
                                                     np.random.uniform(nfibers) * -0.2)
        self.fiberInfo['starCenterErrMajor'][realfibers] = 0.1
        self.fiberInfo['starCenterErrMinor'][realfibers] = 0.1
        self.fiberInfo['starCenterErrCross'][realfibers] = 0.1
        
    def _rowPtrs(self, img, c_type):
        """ Uh-oh. The Mirella routines take arrays of pointers to the
            rows. Make such a thing from a 2-d image.

        Taken from the ctypes tutorial. Is that a strange __call__ in a ctypes type object?
        """
        ptrType = ctypes.POINTER(c_type)
        rowPtrs = (ptrType*len(img))(*[row.ctypes.data_as(ptrType) for row in img])

        return rowPtrs
    
    def rotateStamp(self, img, theta, isMask=False):
        """ rotate an image by theta degrees. Calls the external grot() C function. """

        if not self.gutils:
            return False

        if isMask:
            rotFunc = self.gutils.maskrot
            pixType = ctypes.c_char
        else:
            rotFunc = self.gutils.grot
            pixType = ctypes.c_short
            
        rotImg = np.zeros(img.shape, dtype=img.dtype)
        rotFunc(theta, 
                self._rowPtrs(img, pixType),
                self._rowPtrs(rotImg, pixType),
                img.shape[1], img.shape[0])
        return rotImg
        
    def getOneStamp(self, fullImage, xc, yc, rad, rot, isMask=False):
            img = fullImage[yc-rad:yc+rad+1,
                            xc-rad:xc+rad+1]
            rotImg = self.rotateStamp(img, rot, isMask=isMask)

            return rotImg
        
    def getStamps(self, fullImage, flist, rad, isMask=False):
        """ Generates a list of postage stamps. Note that they are 2*r+1 pixels square. """
        
        n = len(flist)
        stamps = [None] * n

        for i in range(n):
            xc = int(flist['fiberCenterX'][i])
            yc = int(flist['fiberCenterY'][i])
            rot = flist['gp_rot'][i]

            stamps[i] = self.getOneStamp(fullImage, xc, yc, rad, rot, isMask=isMask)
        return stamps

    def getMaskStamps(self, fullMaskImage, flist, r, nplanes):
        n = len(flist)
        stamps = [None] * n

        for p in range(nplanes):
            # Use the regular image rotation function.
            maskPlane = ((fullMaskImage & (1 << p)) >> p).astype('int16')
            planeStamps = self.getStamps(maskPlane, flist, r)

            for i in range(n):
                stampPlane = (planeStamps[i] >= 0.5).astype('uint8') << p

                if p == 0:
                    stamps[i] = stampPlane
                else:
                    stamps[i] |= stampPlane

        return stamps
            
    def getBigStamps(self):
        bigFiberInfo = self.fiberInfo[np.where(self.fiberInfo['gp_fiberType'] == ACQFIBER)]
        self.bigStamps = self.getStamps(self.cleanImage, bigFiberInfo, self.bigFiberSize/2)
        self.bigStampMasks = self.getMaskStamps(self.maskImage, bigFiberInfo, self.bigFiberSize/2, 3)
    def getSmallStamps(self):
        smallFiberInfo = self.fiberInfo[np.where(self.fiberInfo['gp_fiberType'] == GUIDEFIBER)]
        self.smallStamps = self.getStamps(self.cleanImage, smallFiberInfo, self.smallFiberSize/2)
        self.smallStampMasks = self.getMaskStamps(self.maskImage, smallFiberInfo, self.smallFiberSize/2, 3)
            
    def getMaskImage(self):
        """ Stuff the 1-bit mask planes into an 8-bit image. """

        maskImg = (self.satMask != 0).astype('uint8')
        maskImg |= (self.badMask != 0) << 1
        maskImg |= (self.maskedMask != 0) << 2

        return maskImg
        
    def genMasks(self):
        """ Generate dummy pixel masks. """
        self.satMask = (self.cleanImage > 500)
        self.badMask = (self.cleanImage < 0)
        self.maskedMask = (self.cleanImage < 50)

        self.maskImage = self.getMaskImage()
        
    def calcAll(self):
        self.findCenters()
        self.genMasks()
        self.getBigStamps()
        self.getSmallStamps()
        self.calcOffsets()

    def getSmallStampImage(self):
        return np.vstack(self.smallStamps)
    def getBigStampImage(self):
        return np.vstack(self.bigStamps)

    def getSmallMaskStampImage(self):
        return np.vstack(self.smallStampMasks)
    def getBigMaskStampImage(self):
        return np.vstack(self.bigStampMasks)

    def addHeaderCards(self, hdr, metadata):
        names = list(metadata.keys())
        names.sort()

        for n in names:
            hdr.update(n, metadata[n], 'no comment, sorry')
        
    def getFiberHDU(self):
        fibers = self.fiberInfo[self.fiberInfo['fiberType'] != POSFIBER]
        return pyfits.BinTableHDU(fibers)
                                      
    def getObjectHDU(self):
        cols = []
        cols.append(pyfits.Column(name='ra', format='E', array=self.plugmap.ra))
        cols.append(pyfits.Column(name='dec', format='E', array=self.plugmap.dec))
        cols.append(pyfits.Column(name='xFocal', format='E', array=self.plugmap.xFocal))
        cols.append(pyfits.Column(name='yFocal', format='E', array=self.plugmap.yFocal))
        cols.append(pyfits.Column(name='spectrographId', format='E', array=self.plugmap.spectrographId))
        cols.append(pyfits.Column(name='fiberId', format='E', array=self.plugmap.fiberId))

        return pyfits.new_table(cols)

    def attachISRImage(self, cleanImage, maskImage):
        self.cleanimage = cleanImage
        self.makeImage = maskImage

    def attachProbeInfo(self, probeInfo):
        self.probeInfo = probeInfo

    def writeFITS(self, filename):
        
        # The cleaned image
        imageHDU = pyfits.PrimaryHDU(self.cleanImage)
        self.addHeaderCards(imageHDU.header, self.metadata)
        
        # The mask planes.
        maskHDU = pyfits.ImageHDU(self.maskImage)

        # The small fiber postage stamps
        smallStampHDU = pyfits.ImageHDU(self.getSmallStampImage())
        smallMaskStampHDU = pyfits.ImageHDU(self.getSmallMaskStampImage())

        # The big fiber postage stamps
        bigStampHDU = pyfits.ImageHDU(self.getBigStampImage())
        bigMaskStampHDU = pyfits.ImageHDU(self.getBigMaskStampImage())

        # Measured quantities.
        fiberHDU = self.getFiberHDU()

        # Object quantities.
        objectHDU = self.getObjectHDU()

        # Pile 'em all together.
        hlist = pyfits.HDUList()
        hlist.append(imageHDU)
        hlist.append(maskHDU)
        hlist.append(smallStampHDU)
        hlist.append(smallMaskStampHDU)
        hlist.append(bigStampHDU)
        hlist.append(bigMaskStampHDU)
        hlist.append(fiberHDU)
        hlist.append(objectHDU)
        
        try:
            hlist.writeto(filename)
        except Exception, e:
            raise

def main(infile, cartID, plugmap):
    infile = os.path.expandvars(infile)
    infile = os.path.expanduser(infile)
    dirname, filename = os.path.split(infile)
    outfile = os.path.join(dirname, "proc-" + filename)
    
    t = PullTest(infile, cartID, plugmap)
    t.calcAll()
    t.writeFITS(outfile)

if __name__ == "__main__":
    import sys
    
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    
