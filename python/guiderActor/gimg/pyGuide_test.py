#!/usr/bin/python -i

import ctypes
import numpy as np
import os
import pyfits
import pdb
import time

import actorcore.utility.fits as actorFits
reload(actorFits)

class REGION(ctypes.Structure):
	_fields_ = [("nrow", ctypes.c_int),
				("ncol", ctypes.c_int),
				("rows_s16", ctypes.POINTER(ctypes.POINTER(ctypes.c_int16)))]
				
class GHIST(ctypes.Structure):
	_fields_ = [("ghist_medn", ctypes.c_int),
				("ghist_peak", ctypes.c_int),
				("ghist_ref", ctypes.c_int),
				("ghist_refl", ctypes.c_int),
				("ghistarr", ctypes.POINTER(ctypes.c_int))]
				
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
				("g_mag", ctypes.POINTER(ctypes.c_double)),
				("g_fwhm", ctypes.POINTER(ctypes.c_double)),
				("g_poserr", ctypes.POINTER(ctypes.c_double)),
				("g_fitbkgrd", ctypes.POINTER(ctypes.c_double)),
				("g_readnoise", ctypes.c_double),
				("g_npixmask", ctypes.c_int)]



#only needs to be pased to gfindfibers?
class PLATEDATA(ctypes.Structure):
	_fields_ = [("gProbeId", ctypes.POINTER(ctypes.c_int)),
				("exists", ctypes.POINTER(ctypes.c_int)),
				("xcen", ctypes.POINTER(ctypes.c_double)),
				("ycen", ctypes.POINTER(ctypes.c_double)),
				("radius", ctypes.POINTER(ctypes.c_double))]
				
				
#return these classes:
class fiber:
	def __init__(self, fibid, xc, yc, r, illr):
		self.fiberid = fibid
		self.xcen = xc
		self.ycen = yc
		self.radius = r
		self.illrad = illr

class star:
	def __init__(self, fibid, xstar, ystar, poserr, fwhm, strsky, mag):
		self.fiberid = fibid
		self.xs = xstar
		self.ys = ystar
		self.xyserr = poserr
		self.fwhm = fwhm
		self.starsky = strsky
		self.mag = mag

class GuideTest(object):
	def __init__(self, dataname, cartridgeId, gprobes, mode=0, flatname=None, darkname=None, cmd=None):
		self.NFIBERS = 17 + 1	# allow for 1-indexed arrays
		self.cmd = cmd

		self.dataname = dataname
		self.rawdata = pyfits.getdata(dataname)

                # We can't cope with bright pixels. Hack a fix.
                w = np.where(self.rawdata >= 0x8000)
                if len(w[0] > 0):
                        if cmd:
                                cmd.warn('text="the exposure has saturated pixels"')
                        self.rawdata[w[0],w[1]] = 0x7fff
                        
		self.cleandata = self.rawdata.astype('Int16')

		h = pyfits.getheader(dataname) 
		self.dataExptime = h.get('EXPTIME', 1.0)

		if flatname:
			self.rawflat = pyfits.getdata(flatname)
			self.cleanflat = self.rawflat.astype('Int16')
		else:
			self.cleanflat = self.cleandata * 0 + 1

		if darkname:
			self.rawdark = pyfits.getdata(darkname)
			self.cleandark = self.rawdark.astype('Int16')
			h = pyfits.getheader(darkname) 
			self.darkExptime = h.get('EXPTIME', 1.0)
		else:
			self.cleandark = self.cleandata * 0
		self.cleandark = self.cleandark * 0

		self.mode = mode # ???? 
		self.ipGguide = None
		self.ipGguide = ctypes.CDLL(os.path.expandvars("$GUIDERACTOR_DIR/lib/libguide.so"))
		
		self.CID = cartridgeId
		self.gprobes = gprobes
		self.platedata = self.makePlatedata(gprobes, self.CID)
		
		self.ipGguide.gtranspose.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(REGION)]
		self.ipGguide.gmakehist.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(GHIST)]
		self.ipGguide.gmakeflat.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(REGION), ctypes.POINTER(GHIST)]
		self.ipGguide.gextendmask.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(MASK), ctypes.c_int]
		
		self.ipGguide.gfindfibers.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(MASK), 
						      ctypes.POINTER(FIBERDATA), ctypes.POINTER(GHIST), ctypes.POINTER(PLATEDATA)]
		
		self.ipGguide.gfindstars.argtypes = [ctypes.POINTER(REGION), ctypes.POINTER(FIBERDATA), ctypes.c_int]
		# void grot(float, short **, short **, int, int);
		self.ipGguide.grot.restype = None
		self.ipGguide.grot.argtypes = [ctypes.c_float,
					       ctypes.POINTER(ctypes.POINTER(ctypes.c_short)),
					       ctypes.POINTER(ctypes.POINTER(ctypes.c_short)),
					       ctypes.c_int, ctypes.c_int]
        
		# void maskrot(float, char **, char **, int, int);
		self.ipGguide.maskrot.restype = None
		self.ipGguide.maskrot.argtypes = [ctypes.c_float,
						  ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
						  ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
						  ctypes.c_int, ctypes.c_int]

		self.flatarr = np.zeros(self.cleanflat.shape, 'Int16')
		self.darkarr = np.zeros(self.cleandark.shape, 'Int16')
		self.dark = self._getRegion(self.darkarr)
		self.flattener = self._getRegion(self.flatarr)
		self.flatreg = self._getRegion(self.cleanflat)
		self.maskarr = np.zeros(self.cleandata.shape, 'uint8')
		self.mask = MASK(self.cleandata.shape[0], self.cleandata.shape[1], 
				 self._rowPtrs(self.maskarr, ctypes.c_ubyte))
		
		self.readnoise = 5.0
		self.pixmask = 0
		self.fid = np.zeros([self.NFIBERS], 'int32')
		self.xcen = np.zeros([self.NFIBERS], 'double')
		self.ycen = np.zeros([self.NFIBERS], 'double')
		self.fibrad = np.zeros([self.NFIBERS], 'double')
		self.illrad = np.zeros([self.NFIBERS], 'double')
		self.xs = np.zeros([self.NFIBERS], 'double')
		self.ys = np.zeros([self.NFIBERS], 'double')
		self.mag = np.zeros([self.NFIBERS], 'double')
		self.fwhm = np.zeros([self.NFIBERS], 'double')
		self.poserr = np.zeros([self.NFIBERS], 'double')
		self.fitbkgrd = np.zeros([self.NFIBERS], 'double')
		self.ghistarr = np.zeros([32767], 'Int32')
		
		
		self.phist = GHIST(0, 0, 0, 0, self._makePtr(self.ghistarr, ctypes.c_int32))
		self.pfib = FIBERDATA(0,
				self._makePtr(self.fid, ctypes.c_int32),
				self._makePtr(self.xcen, ctypes.c_double), 
				self._makePtr(self.ycen, ctypes.c_double), 
				self._makePtr(self.fibrad, ctypes.c_double), 
				self._makePtr(self.illrad, ctypes.c_double), 
				self._makePtr(self.xs, ctypes.c_double),
				self._makePtr(self.ys, ctypes.c_double), 
				self._makePtr(self.mag, ctypes.c_double), 
				self._makePtr(self.fwhm, ctypes.c_double), 
				self._makePtr(self.poserr, ctypes.c_double), 
				self._makePtr(self.fitbkgrd, ctypes.c_double),
				0.0, 
				0)
				
		self.doubleflat = np.zeros(self.cleanflat.shape, 'double')
		
		self.offset = None
		
	def makePlatedata(self, gprobes, cartridgeID):
		
		self.pprobeId =  np.zeros([self.NFIBERS], int)
		self.pexists =  np.zeros([self.NFIBERS], int)
		self.pxcen =  np.zeros([self.NFIBERS], 'double')
		self.pycen =  np.zeros([self.NFIBERS], 'double')
		self.pradius = np.zeros([self.NFIBERS], 'double')

		for i in range(self.NFIBERS): # the max fibre ID is NFIBERS + 1; the tritium star
			probe = gprobes.get(i, None)
			if not probe:
				self.pexists[i] = False
				continue
			else:
				info = probe.info
				self.pexists[i] = True
				self.pprobeId[i] =  i
				self.pxcen[i] = info.xCenter
				self.pycen[i] = info.yCenter
				self.pradius[i] = info.radius

		self.pprobeId = self.pprobeId.astype(ctypes.c_int)
		self.pexists = self.pexists.astype(ctypes.c_int)
		self.pxcen = self.pxcen.astype(ctypes.c_double)
		self.pycen = self.pycen.astype(ctypes.c_double)
		self.pradius = self.pradius.astype(ctypes.c_double)
		platedata = PLATEDATA(self._makePtr(self.pprobeId, ctypes.c_int), 
								self._makePtr(self.pexists, ctypes.c_int), 
								self._makePtr(self.pxcen, ctypes.c_double),
								self._makePtr(self.pycen, ctypes.c_double),
								self._makePtr(self.pradius, ctypes.c_double))
		return platedata
			
			
	def _makePtr(self, arr, ctype):
		length = len(arr)
		#pdb.set_trace()
		arr = arr.ctypes.data_as(ctypes.POINTER(ctype))
		#ctypes.resize(arr, ctypes.sizeof(ctype*length))
		return arr
	
	def _rowPtrs(self, img, c_type):
		ptrType = ctypes.POINTER(c_type)
		rowPtrs = (ptrType*len(img))(*[row.ctypes.data_as(ptrType) for row in img])

		return rowPtrs
	
	def _getRegion(self, arr):
		x,y = np.shape(arr)
		reg = REGION(x, y, self._rowPtrs(arr, ctypes.c_int16))
		
		return reg
		
	def binImage(self, arr, binsize):
		binshape = tuple([i/binsize for i in np.shape(arr)])
		binpic = np.zeros(binshape)
		for x in range(len(binpic)):
			for y in range(len(binpic[x])):
				mean = (np.sum(arr[binsize*x:binsize*x+binsize, binsize*y:binsize*y+binsize]))/(binsize**2)
				binpic[x][y] = mean
				
		return binpic

	def checkBin(self, im):
		x,y = np.shape(im)
		binx = x/512
		biny = y/512
		if binx>1 and biny>1:
			if binx != biny:
				print "unsymmetric binning, abort"
				return -1
			else:
				binim = self.binImage(im, binx)
		else:
			binim = im
		return binim
	
	def subOverscan(self):
		""" subtract a dummy overscan of 40 for now, will have to depend on overscan later """
		magicBias = 0

		self.cleanflat -= magicBias
		self.cleandata -= magicBias
		# self.cleandark = self.cleandark-magicBias
	    	
	def makeHist(self):
		# This scales the short int flat into regular int values
		
		self.cleanflat = self.cleanflat
		self.flatreg = self._getRegion(self.cleanflat)
		self.ipGguide.gmakehist(ctypes.byref(self.flatreg), ctypes.byref(self.phist))
		
	def makeFlat(self):
		t1 = time.time()
		self.doubleflat = self.cleanflat.astype('double')
		medn = self.phist.ghist_medn
		peak = self.phist.ghist_peak
		ref = self.phist.ghist_ref
		norm = peak - medn
		thresh = ref - medn
		spotx = self.pxcen[16]*2
		spoty = self.pycen[16]*2
		#pdb.set_trace()
		spotr = self.NFIBERS
		self.doubleflat = np.where(self.doubleflat>thresh, norm/(self.doubleflat - medn), 0)
		self.cleanflat = np.where(self.cleanflat>thresh, self.cleanflat, 0)
		t2 = time.time()
		print "time to make flat:", t2-t1
		
		if(self.pexists[16]):
			spotbox = self.doubleflat[spotx-spotr:spotx+spotr, spoty-spotr:spoty+spotr]
			for x in range(len(spotbox)):
				for y in range(len(spotbox)):
					#print self.cleanflat[x][y]
					if ((spotr-x)**2+(spotr-y)**2) < spotr**2:
						self.doubleflat[spoty+y-spotr][spotx+x-spotr] = 1
			t3 = time.time()
			print "time to add spot:", t3-t2
	
	def makeMask(self):
		self.cleanflat = (self.checkBin(self.cleanflat)).astype("int16")
		self.flattener = self._getRegion(self.cleanflat)
		self.ipGguide.gextendmask(ctypes.byref(self.flattener), ctypes.byref(self.mask), 5)
		
	def makeFiberdata(self):
		#need to calculate average offset....
		self.cleanflat = (self.checkBin(self.cleanflat)).astype("int16")
		self.flattener = self._getRegion(self.cleanflat)
		self.ipGguide.gfindfibers(ctypes.byref(self.flattener), ctypes.byref(self.mask), \
			ctypes.byref(self.pfib), ctypes.byref(self.phist), ctypes.pointer(self.platedata))
		
	def addOffset(self):
		""" Add x,y correction derived from the tritium star. """
		correction = (0,0)
		if self.pexists[16]:
			realx = self.xcen[:-1]
			realy = self.ycen[:-1]
			print realx.sum(), realy.sum()
			idealx = self.pxcen[:-1]
			idealy = self.pycen[:-1]
			print idealx.sum(), idealy.sum()
			xoffset = (realx.sum()-idealx.sum())/len(realx)
			yoffset = (realy.sum()-idealy.sum())/len(realy)
			
			correction = (xoffset, yoffset)
			rad = 6
			spot = self.doubleflat[self.pxcen[16]*2+rad:self.pxcen[16]*2-rad, self.pycen[16]*2+rad:self.pycen[16]*2-rad]
			self.doubleflat[self.pxcen[16]*2+rad:self.pxcen[16]*2-rad, self.pycen[16]*2+rad:self.pycen[16]*2-rad] = 0
			newx = int(self.pxcen[16]*2 + xoffset)
			newy = int(self.pycen[16]*2 + yoffset)
			self.doubleflat[newx+rad:newx-rad, newy+rad:newy-rad] = spot
		return correction
			
	def subDark(self):
		pass
		scaledark = (self.dataExptime/self.darkExptime)*self.cleandark
		self.cleandata = self.cleandata - scaledark.astype('i2')

	def makeFinal(self, correction = (0,0)):
		#pdb.set_trace()
		offset = None
		tolerance = 1
		self.doubleflat = self.checkBin(self.doubleflat).astype('double')
		### DISABLE flattening. 
		# self.cleandata = (self.cleandata*self.doubleflat).astype('Int16')
		self.datareg = self._getRegion(self.cleandata)
		self.ipGguide.gfindstars(ctypes.byref(self.datareg), ctypes.byref(self.pfib), self.mode)
		if self.mode == 0:
			xcorrect, ycorrect = correction
			if ((self.xs[16]-xcorrect-self.xcen[16])**2+(self.ys[16]-ycorrect-self.xcen[16])**2) > tolerance**2:
				offset = (self.xs[16]-xcorrect-self.xcen[16], self.ys[16]-ycorrect-self.xcen[16])
		return offset
    	
	def UNUSED_transpose(self, reginarr):
		xreginarr, yreginarr = np.shape(reginarr)
		regin = REGION(xreginarr , yreginarr, self._rowPtrs(new_numpy_array, ctypes.c_ushort))
		new_numpy_array = np.zeros([xreginarr,yreginarr], 'Int16')
		regout = REGION(xreginarr, yreginarr, self._rowPtrs(new_numpy_array, ctypes.c_ushort))
		self.ipGguide.gtranspose(ctypes.pointer(regin), ctypes.pointer(regout))
		return new_numpy_array  
		
	def finalData(self):
		self.fibers = []
		self.stars = []
		for i in range(self.pfib.g_nfibers):
			newfib = fiber(self.fid[i], self.xcen[i], self.ycen[i], self.fibrad[i], self.illrad[i])
			newstar = star(self.fid[i], self.xs[i], self.ys[i], self.poserr[i], self.fwhm[i], self.fitbkgrd[i], self.mag[i])
			self.fibers.append(newfib)
			self.stars.append(newstar)
		return self.fibers, self.stars

	def rotateStamp(self, img, theta, isMask=False):
		""" rotate an image by theta degrees. Calls the external grot() C function. """

		if isMask:
			rotFunc = self.ipGguide.maskrot
			pixType = ctypes.c_ubyte
		else:
			rotFunc = self.ipGguide.grot
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
		# print "img, rotimg: %s, %s" % (str(img.shape), str(rotImg.shape))
		return np.flipud(rotImg)
        
	def getStamps(self, fullImage, probeList, rad, isMask=False):
		""" Generates a list of postage stamps. Note that they are 2*r+1 pixels square. """
		
		stamps = []
		for i in probeList:
			info = self.gprobes[i].info
			xc = int(info.xCenter + 0.5)
			yc = int(info.yCenter + 0.5)
			rot = -info.rotStar2Sky
                        print "rotating probe image %d at (%d,%d) by %0.1f degrees" % (i,xc,yc,rot)
			stamps.append(self.getOneStamp(fullImage, xc, yc, rad, rot, isMask=isMask))
		return stamps

	def getStampImages(self, probeTypes, maskImage, byRadHack=None):
		fiberList = []
		for i in range(self.pfib.g_nfibers):
			probe = self.gprobes.get(i, None)
                        if not probe:
                                continue
                        if byRadHack:
                                if int(probe.info.radius) in byRadHack:
                                        fiberList.append(i)
                        else:
                                if probe.info.fiber_type in probeTypes:
                                        fiberList.append(i)

		if len(fiberList) == 0:
			return (np.zeros([0,0],dtype=self.cleandata.dtype),
                                np.zeros([0,0],dtype=maskImage.dtype),
                                [])
		
		rad = reduce(max, [self.gprobes[i].info.radius for i in fiberList])
		rad = int(rad+0.5)
		imageStamps = self.getStamps(self.cleandata, fiberList, rad)
		maskStamps = self.getStamps(maskImage, fiberList, rad, isMask=True)

                imageStamps = np.vstack(imageStamps)
                maskStamps = np.vstack(maskStamps)
                
                # Fill the stamp background
                med = np.median(self.cleandata)
                w = np.where(imageStamps == 0)
                imageStamps[w[0],w[1]] = med
                
		return imageStamps, maskStamps, fiberList

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
		
	def getMaskImage(self):
		""" Stuff the 1-bit mask planes into an 8-bit image. """

		satMask = (self.rawdata >= 0x7FFF)
		badMask = (self.rawdata == 0) # No idea what is bad. Is there a map?
		maskedMask = (self.rawdata * 0) # Must regenerate image from mask object

		maskImg = (satMask != 0).astype('uint8')
		maskImg |= (badMask != 0) << 1
		maskImg |= (maskedMask != 0) << 2

		return maskImg

	def getAwfulFITSColumn(self, name, npType, fitsType, objs):
		""" Generate a pyfits Column filled with attributes objs.info

                Args:
                     name     - the name of the attribute to look for in objs.info
                     npType   - the numpy type of the given attribute
                     fitsType - the FITS type we want.
                     objs     - the objects whose .info we extract from.

                Returns:
                     - a pyfits Column
                """
                
		col = np.zeros(len(objs), dtype=npType)
		i = 0
		for o in objs.values():
			probeInfo = o.info
                        try:
                                col[i] = getattr(probeInfo,name)
                        except:
                                pass
			i += 1

		return pyfits.Column(name=name, format=fitsType, array=col)
		
	def fillAwfulFITSColumn(self, name, npType, fitsType, nobj, objs):
		""" Fill a pyfits Column with attributes from objs.info

                Args:
                     name     - the name of the attribute to look for in objs.info
                     npType   - the numpy type of the given attribute
                     fitsType - the FITS type we want.
                     nobj     - the number of rows we generate
                     objs     - the objects whose .info we extract from.

                Returns:
                     - a pyfits Column
                """
                
		col = np.zeros(nobj, dtype=npType)
		for o in objs.values():
			fid = o.fiberid-1
                        try:
                                col[fid] = getattr(o,name)
                        except:
                                pass

		return pyfits.Column(name=name, format=fitsType, array=col)
		
	def getProbeHDU(self, cmd, starInfo, smallIds, bigIds):
		""" Return an HDU containing all the probe+star info. 

		Jayzus. Where is all this stuff squirrelled away?
		"""

		probeFields = (('exists','u1','L'),
			       ('enabled','u1','L'),
                               ('xFocal','f4','E'),
                               ('yFocal','f4','E'),
			       ('xCenter','f4','E'),
			       ('yCenter','f4','E'),
			       ('radius','f4','E'),
			       ('xFerruleOffset','f4','E'),
			       ('yFerruleOffset','f4','E'),
			       ('rotation','f4','E'),
			       ('rotStar2Sky','f4','E'),
			       ('focusOffset','f4','E'),
			       ('fiber_type','S20','A20'))

		starFields = (('xstar','f4','E'),
			      ('ystar','f4','E'),
			      ('dx','f4','E'),
			      ('dy','f4','E'),
			      ('dRA','f4','E'),
			      ('dDec','f4','E'),
			      ('fwhm','f4','E'),
			      ('poserr','f4','E'))

		cols = []
		for f in probeFields:
			name, npType, fitsType = f
			col = self.getAwfulFITSColumn(name, npType, fitsType, self.gprobes)
			cols.append(col)

		nobj = len(self.gprobes)
		for f in starFields:
			name, npType, fitsType = f
			col = self.fillAwfulFITSColumn(name, npType, fitsType, nobj, starInfo)
			cols.append(col)

                # Index columns indicating which stampHDU and which image within that
                stampSize = np.zeros(nobj, dtype='i2')-1
                stampIdx = np.zeros(nobj, dtype='i2')-1
                for i in range(len(smallIds)):
                        idx = smallIds[i]
                        stampSize[idx-1] = 1
                        stampIdx[idx-1] = i
                for i in range(len(bigIds)):
                        idx = bigIds[i]
                        stampSize[idx-1] = 2
                        stampIdx[idx-1] = i
                cols.append(pyfits.Column(name='stampSize', format='I', array=stampSize))
                cols.append(pyfits.Column(name='stampIdx', format='I', array=stampIdx))

		hdu = pyfits.new_table(pyfits.ColDefs(cols))
		return hdu

	def getGuideloopCards(self, cmd, frameInfo):
		cards = []
		
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
                        ('offsetScale', 'OFFScale', 'applied scale offset, %'))
		for name, fitsName, comment in defs:
			try:
				val = getattr(frameInfo,name)
				if val != val:
					val = -99999.9 # F.ing F.TS
				c = actorFits.makeCard(cmd, fitsName, val, comment)
				cards.append(c)
			except Exception, e:
				cmd.warn('text="failed to make guider card %s=%s (%s)"' % (
						name, val, e))

		return cards

	def fillPrimaryHDU(self, cmd, models, imageHDU, frameInfo, filename):
		""" Add in all the site and environment keywords. """

		try:
			imageHDU.header.update("OBJECT", os.path.splitext(filename)[0], "")
			imageHDU.header.update("GCAMSCAL", frameInfo.guideCameraScale, "guide camera plate scale (mm/pixel)")
			imageHDU.header.update("PLATSCAL", frameInfo.plugPlateScale, "plug plate scale (mm/degree)")
			tccCards = actorFits.tccCards(models, cmd=cmd)
			actorFits.extendHeader(cmd, imageHDU.header, tccCards)

			mcpCards = actorFits.mcpCards(models, cmd=cmd)
			actorFits.extendHeader(cmd, imageHDU.header, mcpCards)

			plateCards = actorFits.plateCards(models, cmd=cmd)
			actorFits.extendHeader(cmd, imageHDU.header, plateCards)

			guiderCards = self.getGuideloopCards(cmd, frameInfo)
			actorFits.extendHeader(cmd, imageHDU.header, guiderCards)

			self.addPixelWcs(imageHDU.header)
		except Exception, e:
			cmd.warn('text="!!!!! failed to fill out primary HDU  !!!!! (%s)"' % (e))
		
	def writeFITS(self, models, cmd, frameInfo, starInfo):
		try:
			rawHeader = pyfits.getheader(self.dataname)
		except Exception, e:
			cmd.warn("text='could not read raw guider file %s: %s'" % (self.dataname, e,))
			rawHeader = None
		dirname, filename = os.path.split(self.dataname)
		outname = "proc-%s" % (filename)
		procpath = os.path.join(dirname, outname)

		# Start with the raw guider header.
		imageHDU = pyfits.PrimaryHDU(self.cleandata, header=rawHeader)
                imageHDU.header.update("SDSSFMT", "GPROC 1 0", "type major minor version for this file")
		self.fillPrimaryHDU(cmd, models, imageHDU, frameInfo, filename)
                
		try:
			# The mask planes.
			maskImage = self.getMaskImage()
			maskHDU = pyfits.ImageHDU(maskImage)

			# The small fiber postage stamps
			stampImage, stampMaskImage, smallIds = self.getStampImages(probeTypes=('GUIDE', 'TRITIUM'),
                                                                                   maskImage=maskImage,
                                                                                   byRadHack=(8,1))
			smallStampHDU = pyfits.ImageHDU(stampImage)
			smallMaskStampHDU = pyfits.ImageHDU(stampMaskImage)

			# The big fiber postage stamps
			# The old cartridge big fibers are not declared as ACQUIRE, so hack in
			# a radius test. Find a better way, Loomis.
			stampImage, stampMaskImage, bigIds = self.getStampImages(probeTypes=('ACQUIRE',),
                                                                                 maskImage=maskImage,
                                                                                 byRadHack=(14,28))
			bigStampHDU = pyfits.ImageHDU(stampImage)
			bigMaskStampHDU = pyfits.ImageHDU(stampMaskImage)

			# probe input&Measured quantities.
                        #import pdb; pdb.set_trace()
			probeHDU = self.getProbeHDU(cmd, starInfo, smallIds, bigIds)

			# Pile 'em all together.
			hlist = pyfits.HDUList()
			hlist.append(imageHDU)
			hlist.append(maskHDU)
			hlist.append(smallStampHDU)
			hlist.append(smallMaskStampHDU)
			hlist.append(bigStampHDU)
			hlist.append(bigMaskStampHDU)
			hlist.append(probeHDU)
        
			hlist.writeto(procpath, clobber=True)
		except Exception, e:
			cmd.warn("text='could not write proc- guider file: %s'" % (e,))
			# import pdb; pdb.set_trace()
			return

		cmd.inform('file=%s,%s' % (dirname, outname))

	def runAllSteps(self):
		#import pdb; pdb.set_trace()
		self.subOverscan()
		self.makeHist()
		self.makeFlat()
		self.makeMask()
		self.makeFiberdata()
		c = self.addOffset()
		#self.subDark()
		offset = self.makeFinal(correction = c)

		return self.finalData()


def main(argv):
	gt = GuideTest("dummy_flat.fits", "darksum.fits", "dummy_star.fits", 10, mode = 1)
	fibers, stars = gt.runAllSteps()
	return fibers, stars
	
if __name__ == "__main__":
    import sys
    
    data = main("")		
