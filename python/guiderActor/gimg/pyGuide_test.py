import ctypes
import numpy as np
import os
import pyfits
import pdb
import time
import re

import actorcore.utility.fits as actorFits
reload(actorFits)


class REGION(ctypes.Structure):
	_fields_ = [("nrow", ctypes.c_int),
				("ncol", ctypes.c_int),
				("rows_s16", ctypes.POINTER(ctypes.POINTER(ctypes.c_int16)))]
				
class GHIST(ctypes.Structure):
	_fields_ = [("ghist_medn", ctypes.c_int),
				("ghist_peak", ctypes.c_int),
				("ghist_ref",  ctypes.c_int),
				("ghist_refl", ctypes.c_int),
				("ghist_refgf",ctypes.c_int),
                                ("ghist_bias", ctypes.c_int),
				("ghistarr",   ctypes.POINTER(ctypes.c_int))]
				
class MASK(ctypes.Structure):
	_fields_ = [("nrow", ctypes.c_int),
				("ncol", ctypes.c_int),
				("rows", ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)))]
				
                    #output values form gunn code
class FIBERDATA(ctypes.Structure):
	_fields_ = [("g_nfibers",  ctypes.c_int),
                    ("g_fid",      ctypes.POINTER(ctypes.c_int32)),
                    ("g_xcen",     ctypes.POINTER(ctypes.c_double)),
                    ("g_ycen",     ctypes.POINTER(ctypes.c_double)),
                    ("g_fibrad",   ctypes.POINTER(ctypes.c_double)),
                    ("g_illrad",   ctypes.POINTER(ctypes.c_double)),
                    ("g_xs",       ctypes.POINTER(ctypes.c_double)),
                    ("g_ys",       ctypes.POINTER(ctypes.c_double)),
                    ("g_mag",      ctypes.POINTER(ctypes.c_double)),
                    ("g_fwhm",     ctypes.POINTER(ctypes.c_double)),
                    ("g_poserr",   ctypes.POINTER(ctypes.c_double)),
                    ("g_fitbkgrd", ctypes.POINTER(ctypes.c_double)),
                    ("g_fibercts", ctypes.POINTER(ctypes.c_double)),
                    ("g_skycts",   ctypes.POINTER(ctypes.c_double)),
                    ("g_rmswidth", ctypes.POINTER(ctypes.c_double)),
                    ("g_readnoise",ctypes.c_double),
                    ("g_npixmask", ctypes.c_int)]

                    #input from gcamFiberInfo
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
	def __init__(self, fibid, xstar, ystar, poserr, fwhm, strsky, mag, fibercts, skycts, rmswidth):
		self.fiberid = fibid
		self.xs = xstar
		self.ys = ystar
		self.xyserr = poserr
		self.fwhm = fwhm
		self.starsky = strsky
		self.mag = mag
                self.fibercts = fibercts
                self.skycts = skycts
                self.rmswidth = rmswidth

class GuideTest(object):
	def __init__(self, dataname, cartridgeId, gprobes, mode=0, flatname=None, proc = False, darkname=None, cmd=None):
		self.NFIBERS = 17 + 1	# allow for 1-indexed arrays
		self.tritiumTest = False
		self.spotx = 0
		self.spoty = 0
		self.cmd = cmd
                if flatname:
                    self.flatname = flatname
                    self.flatid =  flatNo = re.search(r"([0-9]+)\.fits$", self.flatname).group(1)
                    self.rawflat = pyfits.getdata(flatname)   
                    self.cleanflat = self.rawflat.astype('Int16')
                    self.procdir = os.path.dirname(self.flatname)

                if dataname:
                    self.dataname = dataname
                    self.rawdata = (pyfits.getdata(dataname)).astype('double')

                # Need to do this test after flattening as well
                # We can't cope with bright pixels. Hack a fix.  
                    w = np.where(self.rawdata >= 0x8000)
                    if len(w[0] > 0):
                        if cmd:
                                cmd.warn('text="the raw exposure has saturated pixels"')
                    self.rawdata[w[0],w[1]] = 0x7fff

                    self.cleandata = self.rawdata.astype('Int16')
                    dh = pyfits.getheader(dataname) 
                    self.dataExptime = dh.get('EXPTIME', 1.0)


#		if darkname:
                if False:
			self.rawdark = pyfits.getdata(darkname)
			self.cleandark = self.rawdark.astype('Int16')
			dkh = pyfits.getheader(darkname) 
			self.darkExptime = dkh.get('EXPTIME', 1.0)
		else:
			self.cleandark = self.cleandata * 0


		self.mode = mode # ???? 
		self.ipGguide = None
		self.ipGguide = ctypes.CDLL(os.path.expandvars("$GUIDERACTOR_DIR/lib/libguide.so"))
		
		self.CID = cartridgeId
		self.gprobes = gprobes
		self.platedata = self.makePlatedata(gprobes, self.CID)
		
                #argument type definitions for CTYPES call to Gunn C code subroutines
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
					       ctypes.c_int16, ctypes.c_int16]
        
		# void maskrot(float, char **, char **, int, int);
		self.ipGguide.maskrot.restype = None
		self.ipGguide.maskrot.argtypes = [ctypes.c_float,
						  ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
						  ctypes.POINTER(ctypes.POINTER(ctypes.c_ubyte)),
						  ctypes.c_int, ctypes.c_int]

		self.darkarr = np.zeros(self.cleandark.shape, 'Int16')
		self.dark = self._getRegion(self.darkarr)
		
		self.readnoise = 5.0      #should come form gcameraICC info file
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
		self.fibercts = np.zeros([self.NFIBERS], 'double')
		self.skycts = np.zeros([self.NFIBERS], 'double')
		self.rmswidth = np.zeros([self.NFIBERS], 'double')
		self.ghistarr = np.zeros([32767], 'Int32')
		
		
		self.phist = GHIST(0, 0, 0, 0, 0,0, self._makePtr(self.ghistarr, ctypes.c_int32))

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
				self._makePtr(self.fibercts, ctypes.c_double),
				self._makePtr(self.skycts, ctypes.c_double),
				self._makePtr(self.rmswidth, ctypes.c_double),
				0.0, 
				0)
				

                self.offset = None
		
	def makePlatedata(self, gprobes, cartridgeID):
		self.pprobeId =  np.zeros([self.NFIBERS], int)
		self.pexists =  np.zeros([self.NFIBERS], int)
		self.pxcen =  np.zeros([self.NFIBERS], 'double')
		self.pycen =  np.zeros([self.NFIBERS], 'double')
		self.pradius = np.zeros([self.NFIBERS], 'double')

                #need to decide if should pass tritium star fiber to C code

                #need to handle non existing fibers better, dont want to match to fibers
                #delcared not to exist in the info file
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
		
	def _getMask(self, arr):
		x,y = np.shape(arr)
		reg = MASK(x, y, self._rowPtrs(arr, ctypes.c_ubyte))
		
		return reg

	def slowbinImage(self, arr, binsize):
		binshape = tuple([i/binsize for i in np.shape(arr)])
		binpic = np.zeros(binshape)
		for x in range(len(binpic)):
			for y in range(len(binpic[x])):
				mean = (np.sum(arr[binsize*x:binsize*x+binsize, binsize*y:binsize*y+binsize]))/(binsize**2)
				binpic[x][y] = mean
				
		return binpic

	def binImage(self, arr, binsize):
                '''faster binning 2d image only '''
                t1 = time.time()
                x,y = np.shape(arr)
                x /= binsize ; y /= binsize
                binpic = arr * 1  #create a copy
		binpic.shape = int(x),binsize,int(y),binsize
		binpic = (binpic.sum(axis=3)).sum(axis=1)/(binsize*binsize) 
                t2 = time.time()
                print "time to binImage:", t2-t1                				
#                check image binned correctly
#                procname = "fb-gimg-%(flat)s.fits" %{'flat':self.flatid}
#                hdu = pyfits.PrimaryHDU(binpic)
#                hdu.writeto(procname, clobber = True)

		return binpic


	def checkBin(self, im):
#               fixme : need to read image and flat size from a cgf file                
		x,y = np.shape(im)
		binx = x/512
		biny = y/512
                if binx != biny:
                    print "asymmetric binning, abort"
                    return -1

		if binx==2:
                    binim = self.binImage(im, binx)
		elif binx==1:
                    binim = im
                else:
                    print "incorrect flat size",x,y, "abort"
                    return -1

		return binim	
	def subOverscan(self):
		""" subtract a dummy overscan of 40 for now, will have to depend on overscan later """
		magicBias = 0

		self.cleanflat -= magicBias
		self.cleandata -= magicBias
		# self.cleandark = self.cleandark-magicBias
	    	
	def makeHist(self, data):
		# This scales the short int flat into regular int values
		
		
		self.flatreg = self._getRegion(data)
		self.ipGguide.gmakehist(ctypes.byref(self.flatreg), ctypes.byref(self.phist))
		
	def makeFlat(self):
            #pyflattener is for flattening in python
            #maskflat integer flat is for calculating masks in c
            #note cant add tritium spot until have fiberdata

            #default to processing flat each time for now
                if True:

                    self.makeHist(self.cleanflat) # unbinned
                    #use pseudo bias from hist (60th percentile down, median is 55th % down )
                    #to bias subtract python flattener
                    bias = self.phist.ghist_bias
                    medn = self.phist.ghist_medn - bias
                    peak = self.phist.ghist_peak - bias
                    ref  = self.phist.ghist_ref  - bias 

                    norm = peak - medn
                    thresh = ref - medn
                    
                    self.temp = self.cleanflat.astype('double') - bias
                    self.pyflattener = (self.checkBin(self.temp))                    
                    self.pyflattener = np.where(self.pyflattener>thresh, (peak)/(self.pyflattener), 0)

                    # Create a non-normalized thresholded flat for makemask.c
                    # Mask is expanded around the thresholded regions, this hopefully
                    # circularizes wierd masked regions on bad fibers 
                    # prior to finding the fiber centroids.
                    # Bias level is not subtracted 

                    self.temp = (self.checkBin(self.cleanflat)).astype('int16')
                    self.temp = np.where(self.temp>thresh, self.temp, 0)
                    self.flattener = self._getRegion(self.temp)

                    #setup an empty mask array here to ensure its the same shape as the flattener
                    #a mask is needed for centroiding fibers in the flat

                    self.maskarr = np.zeros(self.cleandata.shape, 'uint8')
                    self.mask = MASK(self.cleandata.shape[0], self.cleandata.shape[1], 
				 self._rowPtrs(self.maskarr, ctypes.c_ubyte))

                    print 'making mask'
                    self.makeMask()
                    print 'making fiber data'
                    self.makeFiberdata()

                    self.writeProcFlat()
                else:
                    #this is where we would reuse the existing flat
                    #need to create a class in masterThread to keep track of flat info
                    print FLATID
                    print self.cleanflat.shape
                    

        def tempWriteFits(self,prefix,imagenum,array):
                ''' write fits images for debugging tests '''
                procdir = os.path.dirname(self.flatname)
                procname = procdir + prefix + "-" + "%(num)s.fits" %{'num':imagenum}
                hdu = pyfits.PrimaryHDU(array)
                hdu.writeto(procname, clobber = True)
                
               
	
	def makeMask(self):
                  
		self.ipGguide.gextendmask(ctypes.byref(self.flattener), ctypes.byref(self.mask), 5)
		
	def makeFiberdata(self):

		self.ipGguide.gfindfibers(ctypes.byref(self.flattener), ctypes.byref(self.mask), \
			ctypes.byref(self.pfib), ctypes.byref(self.phist), ctypes.pointer(self.platedata))

                #save the fiber data 
		self.fibers = []
		for i in range(self.pfib.g_nfibers):
			newfib = fiber(self.fid[i], self.xcen[i], self.ycen[i], self.fibrad[i], self.illrad[i])
			self.fibers.append(newfib)

		#get average offset from pfib.xcen/ycen and move spot by that amount...
		#only done if there is no processed flat
	
		
        def measureFlexure(self):
                """Keep track of the tritum stars location on images """
                #setup a reference position from the first 5 guider images
                #then track the shift there after
                #output tritium star location keywords, Xtrtm, Ytrtm,
                #output tritium star correction keywords XtrtmCorr, YtrtmCorr
                #output information about when to shift the flat

		#if self.flexcorrect            #this used to be the old mode flag
		#	xflex, yflex = flexure
                #
                #check differential value of flexure above some noise threshold
                #       dxflex = self.xs[spotID] - ( xflex + xspotref)
                #       dxfley = self.ys[spotID] - ( yflex + yspotref)
		#	if ((dxflex**2+dyflex**2) > tolerance**2:
		#		xflex += dxflex ; yflex += dyflex
                #
		#	if ((dxflex**2+dyflex**2) > flattolerance**2:
                #          SPOTMOVED = True 
                flexure = 0.0, 0,0                
                return flexure

	def addOffset(self, flexure):
		""" Correct fiber centroids for tritium star flexure. """
                #Correct fiber centroid coords for the flexure
                #This should be a slow change 
                
		pass


        def maskTritium(self):
            """ Set the flat to 1 around the actual location of the tritium star"""
            #compare the fiber.info and measured fiber locations on the flat 
            #to define offset for tritium mask location                

            #The tritium star region of the flat needs to be set to 1.0
            #so we don't flatfield away the tritium spot in the images
            #The tritium star is too faint on the flat to locate precisely 


            #calculate the mean output fiber block offset
            #only want to use coord diffs of fibers actually found
            #worth keeping track of, if someting changes there may be a problem with the cartridge

           #pdb.set_trace()

            realx = self.xcen[:-1] # add where exists
            realy = self.ycen[:-1]
            print realx.sum(), realy.sum()
            idealx = self.pxcen[:-1]
            idealy = self.pycen[:-1]
            print idealx.sum(), idealy.sum()
            xoffset = (realx.sum()-idealx.sum())/len(realx)
            yoffset = (realy.sum()-idealy.sum())/len(realy)

            buttonOffset = (xoffset, yoffset)  # bin image coords


            self.tritiumTest = False
            for id, probe in self.gprobes.items():
                #print "id, prob", id, probe.info.fiber_type
                if probe.info.fiber_type == 'TRITIUM':          #there is by defn only 1 tritium fiber
                    spotx = (probe.info.xCenter*2 + xoffset)    #convert to unbinned coords
                    spoty = (probe.info.yCenter*2 + yoffset) 
                    spotr = probe.info.radius * 2
                    self.tritiumTest = True
                    break



            if self.tritiumTest:
                    rad = 15      #fixme
                    self.pyflattener[spotx-spotr:spotx+spotr, spoty-spotr:spoty+spotr] = 1.


 
        def writeProcFlat(self):
                #write processed flat... need to add check for processed flat
                procdir = os.path.dirname(self.flatname)
                flatNo = re.search(r"([0-9]+)\.fits$", self.flatname).group(1)
                procname = "/proc-gimg-%(flat)s.fits" %{'flat':flatNo}
                rawHeader = pyfits.getheader(self.flatname)
                #add extra keywords:normalization, etc
                hdu = pyfits.PrimaryHDU(self.pyflattener, header=rawHeader)
                hdu.writeto(procdir+procname, clobber = True)


        def subDark(self):
		pass
		scaledark = (self.dataExptime/self.darkExptime)*self.cleandark
		self.cleandata = self.cleandata - scaledark.astype('i2')

	def makeFinal(self):
		offset = None
		tolerance = 1
#		self.pyflattener = self.checkBin(self.pyflattener).astype('double')
		self.cleandata -= self.phist.ghist_bias
		self.cleandata = (self.cleandata*self.pyflattener).astype('Int16')
		self.datareg = self._getRegion(self.cleandata)
		self.ipGguide.gfindstars(ctypes.byref(self.datareg), ctypes.byref(self.pfib), self.mode)

#               write fits file to check image flattened  correctly
                procdir = os.path.dirname(self.dataname)
                imageNo = re.search(r"([0-9]+)\.fits$", self.dataname).group(1)
                procname = procdir+"/proc-gimg-%s.fits" %imageNo

                hdu = pyfits.PrimaryHDU(self.cleandata)
                hdu.writeto(procname, clobber = True)

    	
	def UNUSED_transpose(self, reginarr):
		xreginarr, yreginarr = np.shape(reginarr)
		regin = REGION(xreginarr , yreginarr, self._rowPtrs(new_numpy_array, ctypes.c_ushort))
		new_numpy_array = np.zeros([xreginarr,yreginarr], 'Int16')
		regout = REGION(xreginarr, yreginarr, self._rowPtrs(new_numpy_array, ctypes.c_ushort))
		self.ipGguide.gtranspose(ctypes.pointer(regin), ctypes.pointer(regout))
		return new_numpy_array  
		
	def finalData(self):

#		self.fibers = []
		self.stars = []
		for i in range(self.pfib.g_nfibers):
			newfib = fiber(self.fid[i], self.xcen[i], self.ycen[i], self.fibrad[i], self.illrad[i])
			newstar = star(self.fid[i],   self.xs[i],       self.ys[i],  self.poserr[i], 
                                       self.fwhm[i],  self.fitbkgrd[i], self.mag[i], self.fibercts[i],
                                       self.skycts[i],self.rmswidth)




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
			rot = -info.rotStar2Sky        #NOTE: -ve of rot in gcamFiberInfo.par
                        print "rotating probe image %d at (%d,%d) by %0.1f degrees" % (i,xc,yc,rot)
			stamps.append(self.getOneStamp(fullImage, xc, yc, rad, rot, isMask=isMask))
		return stamps

	def getStampImages(self, probeTypes, maskImage, byRadHack=None, fillBackground=None):
		fiberList = []
		for i in self.gprobes.keys():
			probe = self.gprobes[i]
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
		if fillBackground != None:
			w = np.where(imageStamps == 0)
			imageStamps[w[0],w[1]] = fillBackground
                
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

	def getAwfulFITSColumn(self, name, npType, fitsType, objs, inInfo=True):
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
		for i, o in objs.items():
			fid = i-1
			if inInfo:
				o = o.info
                        try:
                                col[fid] = getattr(o,name)
                        except:
                                pass

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
		for i, o in objs.items():
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

		probeFields = (('exists','u1','L',True),
			       ('enabled','u1','L',False),
			       ('flags','u1','L',False),
                               ('xFocal','f4','E',True),
                               ('yFocal','f4','E',True),
			       ('xCenter','f4','E',True),
			       ('yCenter','f4','E',True),
			       ('radius','f4','E',True),
			       ('xFerruleOffset','f4','E',True),
			       ('yFerruleOffset','f4','E',True),
			       ('rotation','f4','E',True),
			       ('rotStar2Sky','f4','E',True),
			       ('focusOffset','f4','E',True),
			       ('fiber_type','S20','A20',True))

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
			name, npType, fitsType, inInfo = f
			col = self.getAwfulFITSColumn(name, npType, fitsType, self.gprobes, inInfo=inInfo)
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

                imageBackground = np.median(self.cleandata)

		# Start with the raw guider header.
		imageHDU = pyfits.PrimaryHDU(self.cleandata, header=rawHeader)
                imageHDU.header.update("SDSSFMT", "GPROC 1 2", "type major minor version for this file")
                imageHDU.header.update("IMGBACK", imageBackground, "crude background for entire image. For displays.")
		imageHDU.header.update("SEEING", frameInfo.seeing if frameInfo.seeing == frameInfo.seeing else 0.0, 
				       "Estimate of current seeing, arcsec fwhm")

		self.fillPrimaryHDU(cmd, models, imageHDU, frameInfo, filename)
                
		try:
			# The mask planes.
			maskImage = self.getMaskImage()
			maskHDU = pyfits.ImageHDU(maskImage)

			# The small fiber postage stamps
			stampImage, stampMaskImage, smallIds = self.getStampImages(probeTypes=('GUIDE', 'TRITIUM'),
                                                                                   maskImage=maskImage,
                                                                                   byRadHack=(8,1),
										   fillBackground=imageBackground)
			smallStampHDU = pyfits.ImageHDU(stampImage)
			smallMaskStampHDU = pyfits.ImageHDU(stampMaskImage)

			# The big fiber postage stamps
			# The old cartridge big fibers are not declared as ACQUIRE, so hack in
			# a radius test. Find a better way, Loomis.
			stampImage, stampMaskImage, bigIds = self.getStampImages(probeTypes=('ACQUIRE',),
                                                                                 maskImage=maskImage,
                                                                                 byRadHack=(14,28),
										 fillBackground=imageBackground)
			bigStampHDU = pyfits.ImageHDU(stampImage)
			bigMaskStampHDU = pyfits.ImageHDU(stampMaskImage)

			# probe input&Measured quantities.
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

		cmd.inform('file=%s/,%s' % (dirname, outname))

	def runAllSteps(self):
		#import pdb; pdb.set_trace()
#                print 'subtracting overscan'
#		self.subOverscan()
#		print 'making histogram'
#		self.makeHist(self.cleanflat)
#		print 'makeing flat'
		self.makeFlat()
#		print 'making mask'
#		self.makeMask()
#		print 'making fiber data'
#		self.makeFiberdata()
                print 'setting tritium spot area to 1' 
                self.maskTritium()
#		print 'adding offset'
#		c = self.addOffset()
		#self.subDark()
		print 'get final data (find stars)'
		self.makeHist(self.cleandata)
		self.makeFinal()

		return self.finalData()


def main(argv):
	gt = GuideTest("dummy_flat.fits", "darksum.fits", "dummy_star.fits", 10, mode = 1)
	fibers, stars = gt.runAllSteps()
	return fibers, stars
	
if __name__ == "__main__":
    import sys
    
    data = main("")		
