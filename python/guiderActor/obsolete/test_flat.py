import os
import os.path
import re
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib import rc

from numpy import *
from numpy.random import *

from scipy.ndimage.measurements import *

import pyfits

from gimg.guiderImage import *
from gimg.guiderImagePlots import *
from YPF import *

from tests import getGprobes, getModels, getCommand, ducky

from astrometry.util.pyfits_utils import *
from astrometry.util.file import *
import pylab

rc('xtick', labelsize=10)
rc('ytick', labelsize=10)

def savefig(fn):
    print 'Saving', fn
    pylab.savefig(fn)

def makeTestImage(bg, bgsig, fx, fy, fr, fflux, sx, sy, ssig, sflux,
                  W, H):
    X,Y = meshgrid(arange(W),arange(H))
    img = zeros((H,W))
    img += bg + bgsig * standard_normal(img.shape)
    infiber = (((X-fx)**2 + (Y-fy)**2) < fr**2)
    img += fflux * infiber
    img += sflux * 1./(2.*pi*ssig**2) * exp(-((X-sx)**2 + (Y-sy)**2)/(2.*ssig**2)) * infiber
    return img

def measureFiber(bg, bgsig, fx, fy, fr, fflux,
                 sx, sy, ssig, sflux, W, H,
                 imgfn, flatfn, gprobes):
    img = makeTestImage(bg, bgsig, fx, fy, fr, fflux, sx, sy, ssig, sflux, W, H)
    binimg = GuiderImageAnalysis.binImage(img, 2)
    binimg = binimg.astype(int16)
    imhdu = pyfits.PrimaryHDU(binimg)
    h = imhdu.header
    h.update('FLATFILE', flatfn)
    h.update('DARKFILE', '')
    h.update('FLATCART', 1000)
    imhdu.writeto(imgfn, clobber=True)
    #os.system('an-fitstopnm -i %s -r -v | pnmtopng > %s' % (imgfn, imgfn.replace('.fits','.png')))

    #if i % 3 == 0:
    #    #imshow(binimg, origin='lower', interpolation='nearest',
    #    #       extent=[SXs[i]
    #    allims.append(binimg)
    #    allimxs.append(SXs[i])

    GI = GuiderImageAnalysis(imgfn)
    GI.setOutputDir('test-outputs')
    fibers = GI.findFibers(gprobes)
    assert(len(fibers) == 1)
    return fibers[0]

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
    sflux = 100000.

    flatfn = 'test1-flat.fits'
    imgfn  = 'test1-img.fits'

    flat = makeTestImage(bg, bgsig, fx, fy, fr, flatflux, 0, 0, 1, 0, W, H)

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
    for gp in gprobes.values():
        gp.info.rotStar2Sky = 90 + gp.info.rotation - gp.info.phi
    gprobes[1] = gp

    clf()

    sxstep = 1.
    #SXs = arange(34.0, 66.05, sxstep)
    SXs = array([fx])

    #allims = []
    #allimxs = []

    # arcsec/(binned)pix
    pixelscale = 0.428

    if True:
        truefwhm = 1.5
        # magic 2 because we bin by 2 below.
        ssig = 2. * truefwhm / pixelscale / 2.35
        sxstep = 1.
        #SXs = arange(34.0, 66.05, sxstep)
        SXs = arange(30.0, 69.05, sxstep)

        fibers = []
        for i,sx in enumerate(SXs):
            imgfn = 'test1-img-%i.fits' % i
            fiber = measureFiber(bg, bgsig, fx, fy, fr, fflux,
                                 sx, sy, ssig, sflux, W, H,
                                 imgfn, flatfn, gprobes)
            fibers.append(fiber)
            print 'sx:', sx, 'measured:', fiber.xs
        allsx = array([f.xs for f in fibers])

        clf()
        lo,hi = 14*pixelscale, 35*pixelscale
        plot([lo,hi],[lo,hi], 'b-')
        plot(SXs/2. * pixelscale, (allsx+0.25) * pixelscale, 'r.')
        #axvline((fx - fr)/2.*pixelscale, color='0.5')
        #axvline((fx + fr)/2.*pixelscale, color='0.5')
        #axvline(2.2, color='0.5')
        title('Synthetic image test: position')
        xlabel('True star position (arcsec)')
        ylabel('Star position measured by guider (arcsec)')
        axis([lo,hi,lo,hi])
        a=axis()
        imw = (a[1]-a[0]) * 0.1 * 1.3
        imh = (a[3]-a[2]) * 0.1 * 1.3
        for i in [3, 13, 27, 35]:
            ix = SXs[i]/2*pixelscale
            im = pyfits.open('test1-img-%i.fits' % i)[0].data
            dirn = 1 if (ix > 11) else -1
            plot([ix,ix],[11+dirn*imh/2,ix],'b--')
            imshow(im[16:36,15:35], origin='lower', interpolation='nearest',
                   extent=[ix-imw/2, ix+imw/2, 11-imh/2, 11+imh/2])
                #transform=gcf().transFigure,
                #extent=[i*imw, (i+1)*imw, 0, imh])
                #extent=[x, x+imw, a[2], a[2]+imh])
        axis(a)
        savefig('sx.png')
        savefig('sx.pdf')


    if True:
        sx = fx
        truefwhm = arange(0.03, 3.0, 0.03)
        # magic 2 because we bin by 2 below.
        starsigmas = 2. * truefwhm / pixelscale / 2.35

        fibers = []
        for i,ssig in enumerate(starsigmas):
            imgfn = 'test1-img-%i.fits' % i
            fiber = measureFiber(bg, bgsig, fx, fy, fr, fflux,
                                 sx, sy, ssig, sflux, W, H,
                                 imgfn, flatfn, gprobes)
            fibers.append(fiber)
            print 'sigma:', ssig, ', measured fwhm', fiber.fwhm

        allfwhm = array([f.fwhm for f in fibers])
        allmag = array([f.mag for f in fibers])

        clf()
        plot(truefwhm, allfwhm, 'r.')
        #plot([0.8,2.2],[0.8,2.2], 'b.-')
        plot([0,3],[0,3], 'b.-')
        axvline(0.8, color='0.5')
        axvline(2.2, color='0.5')
        print 'slope:', allfwhm[-1] / truefwhm[-1]
        axis([0, 3, 0, 3])
        title('Synthetic image test: FWHM')
        xlabel('True FWHM (arcsec)')
        ylabel('FWHM measured by guider (arcsec)')

        a=axis()
        imw = (a[1]-a[0]) * 0.1 * 1.3
        imh = (a[3]-a[2]) * 0.1 * 1.3
        vmin,vmax = None,None
        for i in [25, 50, 75]:
            ix = truefwhm[i]
            im = pyfits.open('test1-img-%i.fits' % i)[0].data
            if vmin is None:
                vmin,vmax = im.min(), im.max()
                vmax = vmin + (vmax - vmin)*0.5
            y = 0.5
            dirn = 1 if (ix > y) else -1
            plot([ix,ix],[y+dirn*imh/2,ix],'b--')
            imshow(im[16:36,15:35], origin='lower', interpolation='nearest',
                   extent=[ix-imw/2, ix+imw/2, y-imh/2, y+imh/2],
                   vmin=vmin, vmax=vmax)
        axis(a)
        savefig('fwhm.png')
        savefig('fwhm.pdf')

    return


    # arcsec
    fwhm = 1.0
    ssig = 2. * fwhm / pixelscale / 2.35

    truepeakflux = exp(linspace(log(1000), log(32000), 32))
    # bin?
    trueflux = truepeakflux * 2.*pi*ssig**2

    fibers = []
    for i,sflux in enumerate(trueflux):
        imgfn = 'test1-img-%i.fits' % i
        fiber = measureFiber(bg, bgsig, fx, fy, fr, fflux,
                             sx, sy, ssig, sflux, W, H,
                             imgfn, flatfn, gprobes)
        fibers.append(fiber)
        print 'sigma:', ssig, ', measured mag', fiber.mag

    allmag = array([f.mag for f in fibers])
    clf()
    #plot(truepeakflux, allmag, 'r.')
    semilogx(truepeakflux, allmag, 'r.')
    xlabel('True peak flux')
    ylabel('mag')
    savefig('mag2.png')


    return

    SXs /= 2.
    allimxs = array(allimxs)
    allimxs /= 2.

    #plot(SXs, SXs/2-0.25, '-', color='0.5')
    plot(SXs, SXs-0.25, '-', color='0.5')
    plot(SXs, allsx, 'r.')
    #plot(SXs, allsy, 'b.')
    #for sx in [48,49,49.5,50]:
    #    axvline(sx, color='0.5')
    #    I = argmin(abs(SXs - sx))
    #    axhline(allsx[I], color='0.5')

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
    savefig('sx.pdf')



def test_fwhm_flux():
    bg = 1500
    bgsig = sqrt(bg)
    W,H = 100,100
    # fiber:
    fx,fy = 49.5,52
    fr = 8.5 * 2
    fflux = 500
    # star:
    sx,sy = 49.5,fy
    ssig = 2.0 * 2
    sflux = 100000
    img = makeTestImage(bg, bgsig, fx, fy, fr, fflux, sx, sy, ssig, sflux, W, H)



class TestGI(GuiderImageAnalysis):
    def findDarkAndFlat(self, gimgfn, fitsheader):
        flat = fitsheader['FLATFILE']
        flat = os.path.basename(flat)
        return (None, os.path.join(self.basedir, flat))

def flux_calibration():
    #GI = GuiderImageAnalysis(None)
    GI = TestGI(None)
    GI.setOutputDir('test-outputs')
    #GI.printDebug = True

    # get plPlugMap files from speclog.

    testInputs = [
        (55246,  None, None, None, 'plPlugMapM-3684-55246-01.par'),
        #(55246,  1370, None, None, 'plPlugMapM-3684-55246-01.par'),
        #(55247,    2, None, None, 'plPlugMapM-3651-55245-01.par'),
        #(55247,    7, None, None, 'plPlugMapM-3651-55245-01.par'),
        #(55247,   69, None, None, 'plPlugMapM-3651-55245-01.par'),
        #(55247,  175, None, None, 'plPlugMapM-3665-55246-01.par'),
        #(55247,  176, None, None, 'plPlugMapM-3665-55246-01.par'),
        #(55247,  423, None, None, 'plPlugMapM-3685-55247-01.par'),
        #(55247,  835, None, None, 'plPlugMapM-3768-55247-01.par'),
        #(55247, 1028, None, None, 'plPlugMapM-3854-55247-01.par'),
        ]
    for i,t in enumerate(testInputs):
        (mjd, gimgnum, compnum, cart, fn) = t
        pathnm = '../testfiles/%i/' % mjd + fn
        txt = read_file(pathnm)
        cart = int(re.search(r'cartridgeId (\d*)', txt).group(1))
        testInputs[i] = (mjd, gimgnum, compnum, cart, fn)
        print fn, cart

    TI = testInputs[0]
    testInputs = [(TI[0],gimgnum,TI[2],TI[3],TI[4]) for gimgnum in range(533,600,2)]
    #testInputs = [(TI[0],gimgnum,TI[2],TI[3],TI[4]) for gimgnum in range(533,600,10)]
    #testInputs = [(TI[0],gimgnum,TI[2],TI[3],TI[4]) for gimgnum in range(533,534,10)]

    figure(figsize=(6,6))

    allmags = []
    allg = []
    allr = []

    for (mjd, gimgnum, compnum, cart, plugmapbasefn) in testInputs:
        dirnm = '../testfiles/%i/' % mjd
        fn = dirnm + 'gimg-%04i.fits' % gimgnum
        GI.basedir = dirnm
        #print 'Reading file', fn
        plugmapfn = dirnm + plugmapbasefn
        plate = int(re.match(r'plPlugMapM-(\d*)-', plugmapbasefn).group(1))

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
        for gp in gprobes.values():
            gp.info.rotStar2Sky = 90 + gp.info.rotation - gp.info.phi

        GI.gimgfn = fn
        fibers = GI.findFibers(gprobes)
        #print 'Found fibers:', fibers

        frameInfo = ducky()
        frameInfo.seeing = 1.0
        frameInfo.guideCameraScale = 1.0 #GI.pixelscale
        frameInfo.plugPlateScale = 1.0

        GI.writeFITS(getModels(), getCommand(), frameInfo, gprobes)

        #for f in fibers:
        #    print 'Fiber ID', f.fiberid, 'mag', f.mag, 'fwhm', f.fwhm, 'center (%i,%i)' % (int(f.xcen), int(f.ycen)), 'type', f.gprobe.info.fiber_type

        if False:
            clf()
            plot([f.xcen for f in fibers], [f.ycen for f in fibers], 'ro')
            for f in fibers:
                text(f.xcen, f.ycen, 'mag %.1f, fwhm %.2f' %(f.mag, f.fwhm), ha='center', fontsize=8)
            savefig('where.png')
            clf()

        mag = array([f.mag for f in fibers])
        gmag = array([f.gprobe.info.mag[1] for f in fibers])
        rmag = array([f.gprobe.info.mag[2] for f in fibers])

        allmags += list(mag)
        allg += list(gmag)
        allr += list(rmag)
        
    mag = array(allmags)
    gmag = array(allg)
    rmag = array(allr)

    clf()
    #plot(gmag, mag, 'g.')
    #plot(rmag, mag, 'r.')
    grmag = (gmag+rmag)/2.
    #plot(grmag, mag, 'k.')
    #plot(grmag + 0.02*numpy.random.normal(size=grmag.shape), mag, 'k.')
    plot(grmag + 0.01*numpy.random.normal(size=grmag.shape), mag, 'b.', alpha=0.5)
    u = unique(grmag)
    meanmag,stdmag = [],[]
    for m in u:
        I = (grmag == m)
        meanmag.append(mean(mag[I]))
        stdmag.append(std(mag[I]))
    errorbar(u, meanmag, yerr=stdmag, fmt=None)

    print 'Mean offset (g band):', mean(mag - gmag)
    print 'Mean offset (r band):', mean(mag - rmag)
    print 'Mean offset ((g+r)/2 band):', mean(mag - grmag)
    print 'Median offset (g band):', median(mag - gmag)
    print 'Median offset (r band):', median(mag - rmag)
    print 'Median offset ((g+r)/2 band):', median(mag - grmag)

    #xlabel('g / r mag')
    #ylabel('measured mag')
    title('Guide camera zero-point calibration')
    xlabel('Guide star magnitude, (g+r)/2 (mag)')
    ylabel('Calibrated mag (mag)')
    maglo,maghi = min(grmag)-0.25, max(grmag)+0.25
    plot([maglo,maghi], [maglo,maghi],'-', color='0.5')
    axis([maglo,maghi,maglo,maghi])
    savefig('zeropoint.png')
    savefig('zeropoint.pdf')

    #clf()
    #plot(grmag, mag, 'k.')
    #xlabel('(g+r)/2 mag')
    #ylabel('measured mag')
    #savefig('gr.png')
    clf()
    xlabel('g-mag')
    ylabel('r-mag')
    plot(gmag-mag, rmag-mag, 'k.')
    savefig('resids.png')


def acq_fiber_centers():
    ### Test...
    acqs = [f for f in fibers if f.radius > 10]
    f = acqs[0]
    r = int(f.radius * 2) + 5
    y = int(f.ycen * 2)
    x = int(f.xcen * 2)
    img = origflat[y-r:y+r+1, x-r:x+r+1]
    a = [x-r-0.5, x+r+0.5, y-r-0.5, y+r+0.5]
    clf()
    imshow(img,  origin='lower', interpolation='nearest', extent=a)
    axhline(f.ycen*2+0.5, color='r')
    axvline(f.xcen*2+0.5, color='r')
    gca().add_artist(Circle([f.xcen*2+0.5,f.ycen*2+0.5], radius=f.radius*2, fc='none', ec='r'))
    savefig('tst.png')

    i2 = origflat.copy().ravel()
    I = argsort(i2)
    med = i2[I[len(I)/2]]
    pk = i2[I[int(len(I)*0.99)]]
    thresh = (med + pk)/2.

    clf()
    imshow(img>thresh,  origin='lower', interpolation='nearest', extent=a)
    axhline(f.ycen*2+0.5, color='r')
    axvline(f.xcen*2+0.5, color='r')
    gca().add_artist(Circle([f.xcen*2+0.5,f.ycen*2+0.5], radius=f.radius*2, fc='none', ec='r'))
    savefig('tst2.png')

    from scipy.ndimage.filters import gaussian_laplace, gaussian_gradient_magnitude
    #L = gaussian_laplace(img, 2.0)
    L = gaussian_gradient_magnitude(img, 1.0)
    clf()
    imshow(L,  origin='lower', interpolation='nearest', extent=a)
    axhline(f.ycen*2+0.5, color='r')
    axvline(f.xcen*2+0.5, color='r')
    gca().add_artist(Circle([f.xcen*2+0.5,f.ycen*2+0.5], radius=f.radius*2, fc='none', ec='r'))
    savefig('tst3.png')

    #C = zeros_like(img)
    X,Y = meshgrid(arange(x-r, x+r+1), arange(y-r, y+r+1))
    xc,yc = f.xcen*2+0.5, f.ycen*2+0.5
    rad = f.radius * 2
    sig = 1.5

    C = exp(-0.5 * (sqrt((X-xc)**2 + (Y-yc)**2) - rad)**2 / sig**2)
    clf()
    imshow(C,  origin='lower', interpolation='nearest', extent=a)
    axhline(f.ycen*2+0.5, color='r')
    axvline(f.xcen*2+0.5, color='r')
    gca().add_artist(Circle([f.xcen*2+0.5,f.ycen*2+0.5], radius=f.radius*2, fc='none', ec='r'))
    savefig('tst4.png')

    allxyr = []
    for i in range(10):
        rstep = 1./(2.**i)
        xstep = ystep = rstep

        for k in range(10):
            for dr,dx,dy in [ (rstep,0,0), (0,xstep,0), (0,0,ystep) ]:
                for j in range(10):
                    C0 = exp(-0.5 * (sqrt((X-(xc   ))**2 + (Y-(yc   ))**2) - (rad      ))**2 / sig**2)
                    Cp = exp(-0.5 * (sqrt((X-(xc+dx))**2 + (Y-(yc+dy))**2) - (rad+dr))**2 / sig**2)
                    Cm = exp(-0.5 * (sqrt((X-(xc-dx))**2 + (Y-(yc-dy))**2) - (rad-dr))**2 / sig**2)
                    s0 = (C0 * L).sum()
                    sp = (Cp * L).sum()
                    sm = (Cm * L).sum()
                    print
                    print 'x,y', xc,yc, 'radius:', rad
                    print 'dx,dy', dx,dy, 'dradius:', dr
                    print '  -:', sm
                    print '  0:', s0
                    print '  p:', sp

                    if s0 > sp and s0 > sm:
                        break
                    if sp > sm:
                        rad += dr
                        xc += dx
                        yc += dy
                    else:
                        rad -= dr
                        xc -= dx
                        yc -= dy
                    allxyr.append((xc,yc,rad))

    C = exp(-0.5 * (sqrt((X-xc)**2 + (Y-yc)**2) - rad)**2 / sig**2)
    clf()
    imshow(C,  origin='lower', interpolation='nearest', extent=a)
    axhline(f.ycen*2+0.5, color='r')
    axvline(f.xcen*2+0.5, color='r')
    gca().add_artist(Circle([f.xcen*2+0.5,f.ycen*2+0.5], radius=f.radius*2, fc='none', ec='r'))
    #for x,y,r in allxyr:
    #    gca().add_artist(Circle([x,y], radius=r, fc='none', ec='r'))
    savefig('tst5.png')

    clf()
    imshow(img,  origin='lower', interpolation='nearest', extent=a)
    axhline(yc, color='r')
    axvline(xc, color='r')
    gca().add_artist(Circle([xc,yc], radius=rad, fc='none', ec='r'))
    savefig('tst6.png')

def resids(pfn, tt):
    X = unpickle_from_file(pfn)

    # for (xf,yf,dra,ddec,en,gdra,gddec,gdrot,gdscale,gdfocus) in X:

    gras    = array([x[5] for x in X])
    gdecs   = array([x[6] for x in X])
    grots   = array([x[7] for x in X])
    gscales = array([x[8] for x in X])
    gfoci   = array([x[9] for x in X])
    
    gresids = []
    resids2 = []
    residsra = []
    residsdec = []
    residsra2 = []
    residsdec2 = []

    arcsec_per_mm = 3600./217.7

    for (xf,yf,dra,ddec,en,gdra,gddec,gdrot,gdscale,gdfocus) in X:
        # corrections at guide fibers, in arcsec
        rf = sqrt(xf**2 + yf**2)
        cx = gdra *3600. + -yf*deg2rad(gdrot)*arcsec_per_mm + xf*gdscale*.01*arcsec_per_mm
        cy = gddec*3600. +  xf*deg2rad(gdrot)*arcsec_per_mm + yf*gdscale*.01*arcsec_per_mm

        I = logical_and(isfinite(dra), isfinite(ddec), en)

        dra_as = arcsec_per_mm * dra[I]
        ddec_as = arcsec_per_mm * ddec[I]
        cx = cx[I]
        cy = cy[I]

        print
        print 'before correction:', sqrt(mean(dra_as**2)), sqrt(mean(ddec_as**2))
        #print 'after correction:', sqrt(mean((dra_as - cx)**2)), sqrt(mean((ddec_as - cy)**2))

        print 'after correction pos: (%g,%g)' % (gdra,gddec), sqrt(mean((dra_as - gdra*3600.)**2)), sqrt(mean((ddec_as - gddec*3600.)**2))
        print 'after correction rot: (%g)' % gdrot, sqrt(mean((dra_as - -yf[I]*deg2rad(gdrot)*arcsec_per_mm)**2)), sqrt(mean((ddec_as - xf[I]*deg2rad(gdrot)*arcsec_per_mm)**2))
        print 'after correction scale (%g):' % gdscale, sqrt(mean((dra_as - xf[I]*gdscale*0.01*arcsec_per_mm)**2)), sqrt(mean((ddec_as - yf[I]*gdscale*0.01*arcsec_per_mm)**2))

        print 'before correction:', sqrt(mean(dra_as**2 + ddec_as**2))
        print 'after correction pos: (%g,%g)' % (gdra,gddec), sqrt(mean((dra_as - gdra*3600.)**2 + (ddec_as - gddec*3600.)**2))
        print 'after correction rot: (%g)' % gdrot, sqrt(mean((dra_as - -yf[I]*deg2rad(gdrot)*arcsec_per_mm)**2 + (ddec_as - xf[I]*deg2rad(gdrot)*arcsec_per_mm)**2))
        print 'after correction scale (%g):' % gdscale, sqrt(mean((dra_as - xf[I]*gdscale*0.01*arcsec_per_mm)**2 + (ddec_as - yf[I]*gdscale*0.01*arcsec_per_mm)**2))
        print 'after correction all:', sqrt(mean((dra_as - cx)**2 + (ddec_as - cy)**2))

        if False:
            clf()
            plot(xf[I], yf[I], 'r.')
            axhline(0, color='0.5')
            axvline(0, color='0.5')
            sc = 100
            for x,y,dx,dy in zip(xf[I], yf[I], dra_as, ddec_as):
                arrow(x,y,sc*dx,sc*dy)
            arrow(0,0,sc*gdra*3600.,sc*gddec*3600.)
            for x,y,dx,dy in zip(xf[I], yf[I], dra_as, ddec_as):
                arrow(x,y,sc*(dx-gdra*3600),sc*(dy-gddec*3600.), ec='b')
            axis('scaled')
            savefig('err-radec.png')

            clf()
            plot(xf[I], yf[I], 'r.')
            axhline(0, color='0.5')
            axvline(0, color='0.5')
            for x,y,dx,dy in zip(xf[I], yf[I], dra_as, ddec_as):
                arrow(x,y,sc*dx,sc*dy)
            for x,y,dx,dy in zip(xf[I], yf[I], dra_as - yf[I]*deg2rad(gdrot)*arcsec_per_mm, ddec_as - -xf[I]*deg2rad(gdrot)*arcsec_per_mm):
                arrow(x,y,sc*dx, sc*dy, ec='b')
            axis('scaled')
            savefig('err-rot.png')


        gresids.append(sqrt(dra_as**2 + ddec_as**2))

        print 'dra[I]*arcsec_per_mm:', dra[I]*arcsec_per_mm
        print 'mean(dra[I]*arcsec_per_mm):', mean(dra[I]*arcsec_per_mm)
        print 'mean(gras)*3600:', mean(gras)*3600.
        print 'ddec[I]*arcsec_per_mm:', ddec[I]*arcsec_per_mm
        print 'mean(ddec[I]*arcsec_per_mm):', mean(ddec[I]*arcsec_per_mm)
        print 'mean(gdecs)*3600.:', mean(gdecs)*3600.
        
        resids2.append(sqrt((dra[I]*arcsec_per_mm - mean(gras)*3600.)**2 + (ddec[I]*arcsec_per_mm - mean(gdecs)*3600.)**2))
        residsra2.append(dra[I]*arcsec_per_mm - mean(gras)*3600.)
        residsdec2.append(ddec[I]*arcsec_per_mm - mean(gdecs)*3600.)
        residsra.append(dra[I]*arcsec_per_mm)
        residsdec.append(ddec[I]*arcsec_per_mm)


    clf()
    for x,y in zip(residsra,residsra2):
        X = arange(len(x))
        plot(X, x, 'r.')
        plot(X+0.5, y, 'b.')
    savefig('resids-ra.png')
    clf()
    for x,y in zip(residsdec,residsdec2):
        plot(x, 'r.')
        plot(y, 'b.')
    savefig('resids-dec.png')

    clf()
    allresids = hstack(gresids)
    rms = sqrt(mean(allresids[5:] ** 2))
    
    (meanr,stdr) = ([mean(x) for x in gresids],
                    [std(x) for x in gresids])
    #for i,x in enumerate(gresids):
    #    #plot(i*ones_like(x), x*arcsec_per_mm, 'r.')
    errorbar(range(len(meanr)), meanr, yerr=stdr)
    ylabel('Residuals at guide stars (arcsec)')
    xlabel('Guider frame #')
    title(tt)
    ylim(ymax=0.8)
    xlim(xmax=84)
    print 'RMS', rms
    axhline(rms, color=(0.3,0.3,1))
    text(40, 0.7, 'RMS = %.2f arcsec' % rms, horizontalalignment='center', fontsize=12)
    savefig('resids.png')
    savefig('resids.pdf')

    clf()
    for i,x in enumerate(gresids):
        plot(i*ones_like(x), x, 'r.', alpha=0.5)
    ylabel('Residuals at guide stars (arcsec)')
    xlabel('Guider frame #')
    ylim(ymax=1)
    xlim(xmax=84)
    title(tt)
    axhline(rms, color=(0.3,0.3,1))
    text(40, 0.7, 'RMS = %.2f arcsec' % rms, horizontalalignment='center', fontsize=12)
    savefig('all-resids.png')
    savefig('all-resids.pdf')

    errorbar(range(len(meanr)), meanr, yerr=stdr)
    ylim(ymax=1)
    savefig('all-resids2.png')

    

    clf()
    (meanr,stdr) = ([mean(x) for x in resids2],
                    [std(x) for x in resids2])
    errorbar(range(len(meanr)), meanr, yerr=stdr)
    ylabel('Detrended residuals at guide stars (arcsec)')
    ylim(ymax=1)
    savefig('detrend-resids.png')

    clf()
    for i,x in enumerate(resids2):
        plot(i*ones_like(x), x, 'r.', alpha=0.5)
    ylabel('Detrended residuals at guide stars (arcsec)')
    savefig('detrend-all-resids.png')

    
    clf()
    #plot(gras*3600., gdecs*3600., 'r.')
    ymn,ymx,xmx = -0.3,0.3,0.4
    plot(gras*3600., gdecs*3600., 'r.')
    axhline(0, color='0.5')
    axvline(0, color='0.5')
    xlabel('RA guiding correction (arcsec)')
    ylabel('Dec guiding correction (arcsec)')
    title(tt)
    xlim(-0.2, 0.4)
    ylim(-0.3, 0.3)
    savefig('radec.png')
    savefig('radec.pdf')

    clf()
    plot(grots*3600., 'r.')
    xlabel('Guider frame #')
    ylabel('Guider rotation correction (arcsec)')
    axhline(0, color='0.5')
    ylim(-2.5,2.5)
    title(tt)
    savefig('rot.png')
    savefig('rot.pdf')

    clf()
    plot(gscales * 1e4, 'r.')
    ylabel('Guider scale correction (ppm)')
    axhline(0, color='0.5')
    savefig('scale.png')

    #clf()
    #plot(gras, gdecs, 'r.')
    #savefig('g.png')


def summarize(nums, pfn):
    X = []
    for n in nums:
        fn = 'proc-gimg-%04i.fits' % n
        if not os.path.exists(fn):
            continue
        hdr = pyfits.getheader(fn)
        tab = table_fields(pyfits.getdata(fn, 6))
        xf,yf = tab.xfocal, tab.yfocal
        dra,ddec = tab.dra, tab.ddec
        # DRA, DDEC, DROT, DFOCUS, DSCALE
        X.append((xf,yf,dra,ddec,hdr['DRA'], hdr['DDEC'], hdr['DROT'], hdr['DSCALE'], hdr['DFOCUS']))
    pickle_to_file(X, pfn)

if __name__ == '__main__':
    os.environ['GUIDERACTOR_DIR'] = '..'
    fiberinfofn = '../etc/gcamFiberInfo.par'

    flux_calibration()
    sys.exit(0)
    testPixelConventions()
    resids('3666-2.pickle', 'MJD 55259 -- Plate 3666')
    #resids('55247-3854-1.pickle')
    #summarize(range(1029,1463), '55247-3854-1.pickle')
    #summarize(range(63, 123), '3666-1.pickle')

    GI = TestGI(None)
    GI.setOutputDir('test-outputs')


    
    #sys.exit(0)

    testInputs = [
        (55205,     4,   50, 13, 'plPlugMapM-3615-55201-09.par'),
        #(55205,     4,    5, 13, 'plPlugMapM-3615-55201-09.par'),
        #(55243,    1,    7, 10, 'plPlugMapM-3650-55242-01.par'),
        #(55244,  424,  148, 13, 'plPlugMapM-3657-55244-01.par'),
        ]
    nothing = [
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
    #    ]
    #testInputs = [
        (55244,    2,    4, 10, 'plPlugMapM-3650-55242-01.par'),
        (55244,    3,    4, 10, 'plPlugMapM-3650-55242-01.par'),
        (55244,  146,    4, 10, 'plPlugMapM-3650-55242-01.par'),
        (55244,  147,  148, 13, 'plPlugMapM-3657-55244-01.par'),
        #(55244,  424,  148, 13, 'plPlugMapM-3657-55244-01.par'),
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

        outdir = 'test-outputs/%i' % mjd
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        GI.setOutputDir(outdir)

        dirnm = '../testfiles/%i/' % mjd
        fn = dirnm + 'gimg-%04i.fits' % gimgnum
        GI.basedir = dirnm
        print 'Reading file', fn
        plugmapfn = dirnm + plugmapbasefn
        plate = int(re.match(r'plPlugMapM-(\d*)-', plugmapbasefn).group(1))

        comparefn = None
        if compnum is not None:
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
    
        X = GI.analyzeFlat(fn, cart, gprobes, stamps=True)
        if X is None:
            print 'FAILED!'
            continue
        (flat, mask, fibers) = X
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

        #newcenters = [(f.xcen,f.ycen) for f in fibers]
        #oldcenters = [(f.xs,f.ys) for f in fibers]

        origflat = pyfits.getdata(fn)
        halfflat = GuiderImageAnalysis.binImage(origflat, 2)

        if False:
            clf()
            plot_fibers(halfflat, fibers)
            subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
                            wspace=0.05, hspace=0.05)
            text(0.5, 0.98,
                 'MJD %i  -  Gimg %04i  -  Plate %i' % (mjd, gimgnum, plate),
                 transform=gcf().transFigure, fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
            savefig(outbase + '-flat.png')

            clf()
            ax = plot_fibers(halfflat, fibers, dim=0.2)
            for f,a in zip(fibers,ax):
                a.add_artist(Circle([f.xs, f.ys], radius=f.radius, fc='none', ec='r', lw=1.5, alpha=0.5))
            subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
                            wspace=0.05, hspace=0.05)
            text(0.5, 0.98,
                 'MJD %i  -  Gimg %04i  -  Plate %i - Old Centroiding' % (mjd, gimgnum, plate),
                 transform=gcf().transFigure, fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
            savefig(outbase + '-old.png')

            clf()
            ax = plot_fibers(halfflat, fibers, dim=0.2)
            for f,a in zip(fibers,ax):
                a.add_artist(Circle([f.xcen, f.ycen], radius=f.radius, fc='none',
                                    #ec=(0.1,0.1,1),
                                    ec='r', lw=1.5, alpha=0.5))
            subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
                            wspace=0.05, hspace=0.05)
            text(0.5, 0.98,
                 'MJD %i  -  Gimg %04i  -  Plate %i - New Centroiding' % (mjd, gimgnum, plate),
                 transform=gcf().transFigure, fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
            savefig(outbase + '-new.png')

            clf()
            plot_fibers(halfflat, fibers, centerxy=True, showaxes=True,
                        starxy=compare, circles=True)
            subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
                            wspace=0.25, hspace=0.25)
            text(0.5, 0.98,
                 'MJD %i  -  Gimg %04i  -  Plate %i' % (mjd, gimgnum, plate),
                 transform=gcf().transFigure, fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
            savefig(outbase + '.png')

        if True and comparefn and os.path.exists(comparefn):
            img = pyfits.open(comparefn)[0].data
            p = pyfits.open(comparefn)[6].data
            xs = p.field('xStar')
            ys = p.field('yStar')
            for f in fibers:
                f.xs = xs[f.fiberid - 1]
                f.ys = ys[f.fiberid - 1]
            clf()
            plot_fibers(img, fibers, circles=True)#, starxy=True)
            subplots_adjust(left=0.05, right=0.99, bottom=0.03, top=0.95,
                            wspace=0.05, hspace=0.05)
            text(0.5, 0.98,
                 'MJD %i  -  Gimg %04i  -  Plate %i' % (mjd, compnum, plate),
                 transform=gcf().transFigure, fontsize=12,
                 horizontalalignment='center', verticalalignment='top')
            savefig(outbase + '-stars.png')
            sys.exit(0)

        #acq_fiber_centers()


        #clf()
        #plot_fiber_stamps(fibers)
        #savefig('55243-0001-stamps.png')

        #clf()
        #subplot(111)
        #subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95,
        #                wspace=0.05, hspace=0.05)
        #imshow(flat, origin='lower', interpolation='nearest')
        #gray()
        #savefig('55243-0001-flat.png')
    
    if len(alldx):
        clf()
        subplot(2,1,1)
        title('Difference between old and new fiber centroids')
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
        title('Difference between old and new fiber centroids')
        axhline(0, color='0.5', lw=1)
        axvline(0, color='0.5', lw=1)
        plot(array(alldx).ravel(), array(alldy).ravel(), 'b.')
        xlabel('dx (pixels)')
        ylabel('dy (pixels)')
        subplots_adjust(left=0.2, right=0.95, bottom=0.1, top=0.9,
                        wspace=0.2, hspace=0.2)
        axis([-1.5, 1.5, -1.5, 1.5])
        savefig('dxdy.png')
