# wget "http://antwrp.gsfc.nasa.gov/apod/image/1007/FourPlanetSunset_hao.jpg"
# jpegtopnm FourPlanetSunset_hao.jpg | pnmcut -height 500 | ppmtopgm | pnmtofits > cut.fits

import sys
import pyfits
from scipy.ndimage.filters import median_filter, gaussian_filter
from scipy.ndimage.measurements import label, find_objects
from scipy.ndimage.morphology import binary_dilation
from numpy import median, array, argsort, flatnonzero, meshgrid, arange, ones, sum

# given points at (0, f0), (1, f1), (2, f2), assume there is a
# quadratic passing through the three points; return the peak position
# and value of the quadratic:
#
# f = a x^2 + b x + c
#
# df/dx = 2ax + b = 0  =>  x* = -b/2a
#
# f* = a x*^2 + b x* + c
#
def dcen3(f0, f1, f2):
	a = 0.5 * (f2 - 2*f1 + f0)
	if a == 0.:
		print 'a==0 in dcen3'
		return 1
	b = f1 - a - f0
	xc = -0.5 * b / a
	if (xc < 0.0) or (xc > 2.0):
		print 'xc=%g in dcen3' % xc
		return 1
	return xc

def dcen3x3(img):
	# rows
	mx0 = dcen3(img[0,0], img[0,1], img[0,2])
	mx1 = dcen3(img[1,0], img[1,1], img[1,2])
	mx2 = dcen3(img[2,0], img[2,1], img[2,2])

	# cols
	my0 = dcen3(img[0,0], img[1,0], img[2,0])
	my1 = dcen3(img[0,1], img[1,1], img[2,1])
	my2 = dcen3(img[0,2], img[1,2], img[2,2])

	# Fit straight line to peak positions along the rows...
	# x = (y-1) mx + bx */
	bx = (mx0 + mx1 + mx2) / 3.
	mx = (mx2 - mx0) / 2.
	# ... and along the columns...
	# y = (x-1) my + by
	by = (my0 + my1 + my2) / 3.
	my = (my2 - my0) / 2.

	# find intersection
	xc = (mx * (by - my - 1.) + bx) / (1. + mx * my)
	yc = (xc - 1.) * my + by

	# check that we are in the box
	if xc < 0 or xc > 2 or yc < 0 or yc > 2:
		print 'xc,yc = %g,%g in dcen3x3' % (xc,yc)
		return (1,1)
	return (xc,yc)


def find_single_star(img, medfiltsize = 10, psfsigma = 1.0,
					 min_npix = 10, noise_est_margin=10, noise_est_dpix=5,
					 grow = 5):
	# median filter to remove background
	bg = median_filter(img, size=medfiltsize)
	img = img - bg

	# convolve by Gaussian to smooth noise
	img = gaussian_filter(img, psfsigma)

	# estimate noise in the smoothed image using the Blanton trick
	# of looking at pixel pairs separated by 5 pixels
	(H,W) = img.shape
	# Avoid the edges
	margin,dp = noise_est_margin,noise_est_dpix
	noise = median(abs((img[margin:H-(margin+dp), margin:W-(margin+dp)] - img[(margin+dp):H-margin, (margin+dp):W-margin]).ravel()))
	print 'Estimated noise:', noise

	# binary search for a threshold that produces only one peak.
	siglo,sighi = 10,10000
	iters = 0
	while True:
		iters += 1
		if iters > 100:
			print 'Failed to find a threshold that yielded exactly one object with %i pixels above threshold.' % (min_npix)
			return None
		nsigma = (siglo + sighi) / 2.
		thresh = nsigma * noise
		# Find pixels above threshold
		detected = (img > thresh)
		# Label connected components
		(labels, nobjs) = label(detected)
		slices = find_objects(labels, nobjs)
		# Count number of pixels above threshold in each connected component
		npix = array([sum((labels[s] == i+1).ravel()) for i,s in enumerate(slices)])
		# Find components with more than min_npix pixels above threshold
		I = (npix >= min_npix)
		print ('Trying nsigma=%g, threshold=%g ==> detected %i objects, %i with more than %i pixels above threshold' %
			   (nsigma, thresh, nobjs, sum(I), min_npix))

		if sum(I) == 0:
			# No objects detected -- lower threshold.
			sighi = nsigma
		elif sum(I) == 1:
			# just right.
			break
		else:
			# raise threshold.
			siglo = nsigma

	print 'N sigma:', nsigma

	# Which component is it?
	i = flatnonzero(I)[0]
	# Pull out the region of interest.
	slc = slices[i]
	(sy,sx) = slc
	# Grow region of interest
	# The +1 here isn't actually necessary.
	sx = slice(sx.start - (grow+1), sx.stop + (grow+1))
	sy = slice(sy.start - (grow+1), sy.stop + (grow+1))
	slc = (sy,sx)
	x0,y0 = sx.start, sy.start
	subimg = img[slc]
	subobj = (labels[slc] == i+1)
	subobj = binary_dilation(subobj, iterations=grow, structure=ones((3,3)))
	(SH,SW) = subobj.shape

	# Compute the centroids.
	# "subobj" is a binary mask image.
	norm = sum(subobj * subimg)
	mx = sum(sum(subobj * subimg, axis=0) * arange(SW)) / norm
	my = sum(sum(subobj * subimg, axis=1) * arange(SH)) / norm
	ox = x0 + mx
	oy = y0 + my

	# Starting at the centroid position, find the peak of the paraboloid
	# using Blanton's dcen3x3.
	ix,iy = int(round(mx)), int(round(my))
	xc,yc = dcen3x3(subimg[iy-1:iy+2, ix-1:ix+2])
	xc += ix - 1 + x0
	yc += iy - 1 + y0

	ix,iy = int(round(xc)), int(round(yc))
	
	return dict(centroid = (xc,yc),
				background = bg[iy, ix],
				flux = norm,
				subimg = subimg,
				submask = subobj,
				offset = (x0,y0),
				)



if __name__ == '__main__':
	img = pyfits.open('cut.fits')[0].data.astype(float)

	D = find_single_star(img)

	subimg = D['subimg']
	subobj = D['submask']
	x0,y0 = D['offset']
	xc,yc = D['centroid']
	flux = D['flux']
	bg = D['background']

	print 'Flux', flux
	print 'Background', bg
	
	from pylab import *

	clf()
	subplot(1,2,1)
	imshow(subimg, interpolation='nearest', origin='lower')
	a=axis()
	plot([xc-x0], [yc-y0], 'r.')
	axis(a)
	subplot(1,2,2)
	imshow(subobj, interpolation='nearest', origin='lower')
	savefig('sub.png')

	clf()
	imshow(img, interpolation='nearest', origin='lower')
	a=axis()
	plot([xc], [yc], 'r.')
	axis(a)
	colorbar()
	savefig('obj.png')

	axis([xc-10, xc+10, yc-10, yc+10])
	savefig('obj2.png')

	
