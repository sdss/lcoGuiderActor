import os
import os.path
import re
import sys

import matplotlib
matplotlib.use('Agg')

from numpy import *
from numpy.random import *

from scipy.ndimage.measurements import *

import pyfits

from gimg.guiderImage import *
from gimg.guiderImagePlots import *
from YPF import *

from tests import getGprobes, ducky

#class TestGuiderImageAnalysis(GuiderImageAnalysis):
#pass

def testPixelConventions():
	bg = 1500
	bgsig = sqrt(bg)

	W,H = 100,100
	# fiber:
	fx,fy = 49.5,52
	fr = 8.5 * 2
	fflux = 500
	flatflux = 10000
	# star:
	sx,sy = 49.5,fy
	ssig = 2.0 * 2
	sflux = 100000

	flatfn = 'test1-flat.fits'
	imgfn  = 'test1-img.fits'

	X,Y = meshgrid(arange(W),arange(H))
	flat = zeros((H,W), int16)
	flat += bg
	flat += flatflux * (((X-fx)**2 + (Y-fy)**2) < fr**2)

	pyfits.PrimaryHDU(flat).writeto(flatfn, clobber=True)

	gprobes = {}
	gp = ducky()
	gp.enabled = True
	gp.flags = 0
	gp.info = ducky()
	gp.info.exists = True
	gp.info.fiber_type = 'GUIDE'
	gp.info.xCenter = 25
	gp.info.yCenter = 25
	gp.info.radius = fr/2.
	gp.info.xFerruleOffset = 0
	gp.info.yFerruleOffset = 0
	gp.info.rotation = 0
	gp.info.phi = 0
	gp.info.rotStar2Sky = 0
	gp.info.focusOffset = 0
	gp.info.ra = 0
	gp.info.dec = 0
	gp.info.xFocal = 0
	gp.info.yFocal = 0
	gprobes[1] = gp

	clf()

	#sxstep = 0.05
	#sxstep = 0.25
	sxstep = 1.
	#SXs = arange(34.0, 66.05, sxstep)
	#SXs = arange(20.0, 70.05, sxstep)
	#SXs = arange(45.0, 53.05, sxstep)
	#SXs = arange(47.0, 51.05, 0.05)
	#SXs = arange(47.0, 51.05, 0.25)
	SXs = array([fx])
	bgsig = 0.

	allsx = []
	allsy = []

	allims = []
	allimxs = []

	for i,sx in enumerate(SXs):
		imgfn = 'test1-img-%.2f.fits' % sx

		img = zeros((H,W))
		img += bg + bgsig * standard_normal(img.shape)
		img += fflux * (((X-fx)**2 + (Y-fy)**2) < fr**2)
		img += sflux * 1./(2.*pi*ssig**2) * exp(-((X-sx)**2 + (Y-sy)**2)/(2.*ssig**2))

		S = sflux * 1./(2.*pi*ssig**2) * exp(-((X-sx)**2 + (Y-sy)**2)/(2.*ssig**2))
		#print 'star flux cm', center_of_mass(S)
		#print 'star flux', S.min(), S.max(), S.sum()
		binimg = GuiderImageAnalysis.binImage(img, 2)
		#print 'binimg:', binimg.min(), binimg.max()
		binimg = binimg.astype(int16)
		#print 'binned cm', center_of_mass(binimg)
		#print 'binimg:', binimg.min(), binimg.max()
		imhdu = pyfits.PrimaryHDU(binimg)
		h = imhdu.header
		h.update('FLATFILE', flatfn)
		h.update('DARKFILE', '')
		h.update('FLATCART', 1000)
		imhdu.writeto(imgfn, clobber=True)
		#os.system('an-fitstopnm -i %s -r -v | pnmtopng > %s' % (imgfn, imgfn.replace('.fits','.png')))

		if i % 3 == 0:
			#imshow(binimg, origin='lower', interpolation='nearest',
			#	   extent=[SXs[i]
			allims.append(binimg)
			allimxs.append(SXs[i])

		GI = GuiderImageAnalysis(imgfn)
		GI.setOutputDir('test-outputs')
		fibers = GI.findFibers(gprobes)
		assert(len(fibers) == 1)
		print 'fiber:', fibers[0]
		allsx.append(fibers[0].xs)
		allsy.append(fibers[0].ys)

	allsx = array(allsx)
	allsy = array(allsy)

	SXs /= 2.
	allimxs = array(allimxs)
	allimxs /= 2.

	#plot(SXs, SXs/2-0.25, '-', color='0.5')
	plot(SXs, SXs-0.25, '-', color='0.5')
	plot(SXs, allsx, 'r.')
	#plot(SXs, allsy, 'b.')
	#for sx in [48,49,49.5,50]:
	#	axvline(sx, color='0.5')
	#	I = argmin(abs(SXs - sx))
	#	axhline(allsx[I], color='0.5')

	a = axis()
	imw = (a[1]-a[0]) * 0.05 * 1.3
	imh = (a[3]-a[2]) * 0.05 * 1.3 # adjust for plot aspect ratio
	#imw = 1./(len(allims))
	#imh = 1./(len(allims))
	for i,(x,im) in enumerate(zip(allimxs, allims)):
		imshow(im[15:35,15:35], origin='lower', interpolation='nearest',
			   #transform=gca().transAxes,
			   transform=gcf().transFigure,
			   #extent=[i*imw, (i+1)*imw, 0, imh])
			   extent=[x, x+imw, a[2], a[2]+imh])
		gray()
	axis(a)

	xlabel('True star x')
	ylabel('Found star x')

	title('slope %.3f' % ((allsx[-1]-allsx[0])/(SXs[-1]-SXs[0])))
	savefig('sx.png')


if __name__ == '__main__':
	os.environ['GUIDERACTOR_DIR'] = '..'
	fiberinfofn = '../etc/gcamFiberInfo.par'

	#testPixelConventions()
	#sys.exit(0)
	
	GI = GuiderImageAnalysis(None)
	GI.setOutputDir('test-outputs')


	testInputs = [
		(55243,    1,    7, 10, 'plPlugMapM-3650-55242-01.par'),
		(55243,    2,    7, 10, 'plPlugMapM-3650-55242-01.par'),
		(55243,    5,    7, 10, 'plPlugMapM-3650-55242-01.par'),
		(55243,    6,    7, 10, 'plPlugMapM-3650-55242-01.par'),
		(55243,   85,    7, 10, 'plPlugMapM-3650-55242-01.par'),
		(55243,   86,   87, 12, 'plPlugMapM-3681-55242-01.par'),
		(55243,  288,   87, 12, 'plPlugMapM-3681-55242-01.par'),
		(55243,  289,  290, 13, 'plPlugMapM-3781-55241-02.par'),
		(55243,  469,  290, 13, 'plPlugMapM-3781-55241-02.par'),
		(55243,  470,  471, 17, 'plPlugMapM-3782-55242-01.par'),
		(55243,  784,  471, 17, 'plPlugMapM-3782-55242-01.par'),
		(55243,  785,  786, 16, 'plPlugMapM-3774-55241-01.par'),
		(55243,  926,  786, 16, 'plPlugMapM-3774-55241-01.par'),
		(55243,  927,  928, 11, 'plPlugMapM-3852-55242-01.par'),
		(55243, 1193,  928, 11, 'plPlugMapM-3852-55242-01.par'),

		(55244,    2,    4, 10, 'plPlugMapM-3650-55242-01.par'),
		(55244,    3,    4, 10, 'plPlugMapM-3650-55242-01.par'),
		(55244,  146,    4, 10, 'plPlugMapM-3650-55242-01.par'),
		(55244,  147,  148, 13, 'plPlugMapM-3657-55244-01.par'),
		(55244,  424,  148, 13, 'plPlugMapM-3657-55244-01.par'),
		(55244,  425,  426, 12, 'plPlugMapM-3682-55244-01.par'),
		(55244,  645,  426, 12, 'plPlugMapM-3682-55244-01.par'),
		(55244,  646,  647, 17, 'plPlugMapM-3782-55242-01.par'),
		(55244,  776,  647, 17, 'plPlugMapM-3782-55242-01.par'),
		(55244,  777,  778, 16, 'plPlugMapM-3774-55241-01.par'),
		(55244, 1018,  778, 16, 'plPlugMapM-3774-55241-01.par'),
		(55244, 1019, 1020, 11, 'plPlugMapM-3879-55244-02.par'),
		(55244, 1548, 1020, 11, 'plPlugMapM-3879-55244-02.par'),
		]

	alldx = []
	alldy = []

	figure(figsize=(6,6))
	for (mjd, gimgnum, compnum, cart, plugmapbasefn) in testInputs:

		#if mjd != 55244:
		#	continue

		fn = '../testfiles/%i/gimg-%04i.fits' % (mjd, gimgnum)
		plugmapfn = '../testfiles/%i/' % mjd + plugmapbasefn

		plate = int(re.match(r'plPlugMapM-(\d*)-', plugmapbasefn).group(1))

		comparefn = '../testfiles/%i/proc-gimg-%04i.fits' % (mjd, compnum)

		print
		print 'MJD %i  --  Gimg %04i  --  cart %i  --  plugmap %s  --  plate %i' % (mjd, gimgnum, cart, plugmapbasefn, plate)
		print

		gprobes = getGprobes(fiberinfofn, plugmapfn, cart)
		for k,v in gprobes.items():
			d = ducky()
			d.info = v
			d.enabled = True
			d.flags = 0
			gprobes[k] = d
	
		(flat, mask, fibers) = GI.analyzeFlat(fn, cart, gprobes, stamps=True)
		fibers = [f for f in fibers if not f.is_fake()]

		print 'Got %i fibers' % len(fibers)
		for f in fibers:
			print '  ', f

		outbase = '%i-%04i' % (mjd, gimgnum)

		compare = False
		if comparefn and os.path.exists(comparefn):
			compare = True
			p = pyfits.open(comparefn)[6].data
			xc = p.field('xCenter')
			yc = p.field('yCenter')
			for f in fibers:
				f.xs = xc[f.fiberid - 1]
				f.ys = yc[f.fiberid - 1]

			alldx += [f.xs - f.xcen for f in fibers]
			alldy += [f.ys - f.ycen for f in fibers]

		clf()
		plot_fibers(flat, fibers, centerxy=True, showaxes=True,
					starxy=compare)
		subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
						wspace=0.25, hspace=0.25)

		text(0.5, 0.98,
			 'MJD %i  -  Gimg %04i  -  Plate %i' % (mjd, gimgnum, plate),
			 transform=gcf().transFigure, fontsize=12,
			 horizontalalignment='center', verticalalignment='top')

		savefig(outbase + '.png')

		#clf()
		#plot_fiber_stamps(fibers)
		#savefig('55243-0001-stamps.png')

		#clf()
		#subplot(111)
		#subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
		#				wspace=0.05, hspace=0.05)
		#imshow(flat, origin='lower', interpolation='nearest')
		#gray()
		#savefig('55243-0001-flat.png')
	
	clf()
	subplot(2,1,1)
	title('(Dustin - As-Run) fiber centroids')
	hist(array(alldx).ravel(), 20)
	xlabel('dx (pixels)')
	subplot(2,1,2)
	hist(array(alldy).ravel(), 20)
	xlabel('dy (pixels)')
	subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.9,
					wspace=0.2, hspace=0.2)
	savefig('dxdy-hist.png')

	clf()
	subplot(111)
	title('(Dustin - As-Run) fiber centroids')
	plot(array(alldx).ravel(), array(alldy).ravel(), 'b.')
	xlabel('dx (pixels)')
	ylabel('dy (pixels)')
	subplots_adjust(left=0.2, right=0.95, bottom=0.1, top=0.9,
					wspace=0.2, hspace=0.2)
	savefig('dxdy.png')
