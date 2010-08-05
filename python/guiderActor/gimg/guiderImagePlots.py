from numpy import isfinite
from pylab import *
from matplotlib import rc
rc('xtick', labelsize=8)
rc('ytick', labelsize=8)

def plot_fiber_stamps(fibers, showaxes=False):
	clf()
	fibers = [f for f in fibers if hasattr(f, 'stamp')]
	for i,f in enumerate(fibers):
		subplot(4,4,i+1)
		#print 'f.stamp:', f.stamp.shape, f.stamp.min(), f.stamp.max()
		imshow(f.stamp, origin='lower', interpolation='nearest')
		gray()
		if not showaxes:
			axis('off')
	if not showaxes:
		subplots_adjust(left=0, right=1, bottom=0, top=1,
						wspace=0.05, hspace=0.05)
	else:
		subplots_adjust(left=0.1, right=1, bottom=0.05, top=1,
						wspace=0.25, hspace=0.25)
	

def plot_fibers(image, fibers, showaxes=False, centerxy=False, starxy=False, R=0):
	clf()
	for i,f in enumerate(fibers):
		subplot(4,4,i+1)
		if R == 0:
			if f.radius > 10:
				N = 60
			else:
				N = 21
		else:
			N = R
		y = int(f.ycen)
		x = int(f.xcen)
		stamp = image[y-N/2:y+N/2+1, x-N/2:x+N/2+1]
		#print 'stamp range:', stamp.min(), stamp.max()
		a = [x-N/2-0.5, x+N/2+0.5, y-N/2-0.5, y+N/2+0.5]
		imshow(stamp, origin='lower', interpolation='nearest',
			   extent=a)
		gray()
		if centerxy:
			axhline(f.ycen, color='c')
			axvline(f.xcen, color='c')
		if starxy:
			if isfinite(f.ys):
				axhline(f.ys, color='r')
			if isfinite(f.xs):
				axvline(f.xs, color='r')
		axis(a)
		if not showaxes:
			axis('off')
		axis(a)
		#print axis()
	if not showaxes:
		subplots_adjust(left=0, right=1, bottom=0, top=1,
						wspace=0.05, hspace=0.05)
	else:
		subplots_adjust(left=0.1, right=1, bottom=0.05, top=1,
						wspace=0.25, hspace=0.25)

