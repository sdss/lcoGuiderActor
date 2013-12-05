#!/usr/bin/env python
"""
Generate a movie from processed guider files to show the guider behavior.
"""
import subprocess
import tempfile
import shutil
import os.path
import sys

import pyfits
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import numpy as np

from opscore.utility import assembleImage

gimgbase = 'proc-gimg-%04d.fits.gz'
tempbase = 'temp-gimg-%04d.png'

def asinh(inputArray, scale_min=None, scale_max=None, non_linear=2.0):
    """
    Performs asinh scaling of the input numpy array.
    Taken from:
        http://dept.astro.lsa.umich.edu/~msshin/science/code/Python_fits_image/img_scale.py
    """
    
    imageData=np.array(inputArray, copy=True)
    if scale_min == None:
        scale_min = imageData.min()
    if scale_max == None:
        scale_max = imageData.max()
    factor = np.arcsinh((scale_max - scale_min)/non_linear)
    indices0 = np.where(imageData < scale_min)
    indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
    indices2 = np.where(imageData > scale_max)
    imageData[indices0] = 0.0
    imageData[indices2] = 1.0
    imageData[indices1] = np.arcsinh((imageData[indices1] - scale_min)/non_linear)/factor
    
    return imageData
#...

def makeGProbeName(gprobeNum, gprobeBits):
    """Construct a guide probe name from its number and gProbeBits
    
    Inputs:
    gprobeNum: guide probe number (an integer, though a string will do)
    gprobeBits: guide probe bits; if None then the above/below focus suffix is not added
    """
    aboveFocus = '+'
    belowFocus = '-'
    if gprobeBits is None:
        suffixStr = ""
    elif gprobeBits & 1<<3 != 0:
        suffixStr = aboveFocus
    elif gprobeBits & 1<<4 != 0:
        suffixStr = belowFocus
    else:
        suffixStr = ""
    return str(gprobeNum) + suffixStr

class ImageMaker(object):
    """Initialize with the file to read. Call to get an image saved to that filename."""
    def __init__(self,infile):
        np.seterr(all='ignore')
        self.cmap1 = 'gray'
        self.cmap2 = 'gist_ncar'
        self.aspect = 'auto' # normal vs. equal?
        
        ErrPixPerArcSec = 40 # pixels per arcsec of error on the plug plate 
        scale_min = 5.
        scale_max = 30000.
        
        self.width = 512*2. # 10px buffer
        self.height = 512.
        self.dpi = 100.
        self.qscale = 10.
        self.assembler = assembleImage.AssembleImage(1,0)
        
        self.index = os.path.splitext(infile)[0].split('-')[-1]
        data = pyfits.open(infile)
        self.plate = data[0].header['PLATEID']
        self.cart =  data[0].header['CARTID']
            
        
        self.guiderView = asinh(data[0].data,scale_min=5.,scale_max=scale_max,non_linear=10.)
        self.guiderSat = np.nonzero(data[1].data == 1)
        self.guiderBad = np.nonzero(data[1].data == 2)
        
        self.plateInfo = self.assembler(data)
        plateView = self.plateInfo.plateImageArr
        zeros = plateView == 0
        self.plateView = asinh(plateView,scale_min=-5.,scale_max=scale_max,non_linear=10.)
        self.plateView[zeros] = -999
        mask = self.plateInfo.plateMaskArr
        self.plateView = np.ma.masked_array(self.plateView, mask == 4)
        
        hdr = data[0].header
        self.seeing = hdr['seeing']
        self.offset = tuple(np.array((hdr['dra'],hdr['ddec'],hdr['drot']))*3600.)
        self.focus = hdr['dfocus']
        self.scale = hdr['filtscle']
        self.rms = hdr['gdrms']
        self.time = hdr['date-obs']
        self.ccdtemp = hdr['ccdtemp']
        
        self._guide_locations()
        # dictionaries to hold information about the probes in each view
        self._labels(self.plateInfo.stampList)
    #...
    
    def _labels(self,stampList):
        """Make the labels for each guide probe for both plots."""
        self.guideLabels = {}
        for i in range(1,17):
            self.guideLabels[i] = [self.guiderLocations[i],str(i),self.guiderLocations[i],'disabled',None]
        self.plateLabels = dict([[i,[]] for i in range(1,17)])
        vectors = []
        for stamp in self.plateInfo.stampList:
            self.guideLabels[stamp.gpNumber] = self._guide_labels(stamp)
            self.plateLabels[stamp.gpNumber] = self._plate_labels(stamp)
            center = stamp.decImCtrPos
            error = stamp.starRADecErrArcSec
            vectors.append(np.array([center[0],center[1],error[0],error[1]]))
        self.vectors = np.array(vectors).T

    def _plate_labels(self,stamp):
        """Return the positions, probe names, and fwhm for the plate view."""
        boxWidth = stamp.image.shape[0] / 2.0
        probeName = makeGProbeName(stamp.gpNumber, stamp.gpBits)
        center = stamp.decImCtrPos
        loc1 = center + (boxWidth + 2,5)
        loc2 = center - (boxWidth + 5,30)
        fwhm = stamp.fwhmArcSec
        if not stamp.gpEnabled:
            disabled = [np.array([center+(boxWidth,boxWidth),center-(boxWidth,boxWidth)]).T,
                        np.array([center+(boxWidth,-boxWidth),center-(boxWidth,-boxWidth)]).T]
        else:
            disabled = None
        return loc1,probeName,loc2,fwhm,disabled
    
    def _guide_labels(self,stamp):
        """Return the location positions1/2, probe names, and fwhm for the guider view."""
        center = stamp.gpCtr
        boxWidth = stamp.image.shape[0] / 2.0
        probeName = makeGProbeName(stamp.gpNumber, stamp.gpBits)
        fwhm = stamp.fwhmArcSec
        loc1 = center-(boxWidth+4,-5)
        loc2 = center-(0,boxWidth+17)
        if not stamp.gpEnabled:
            disabled = [np.array([center+(boxWidth,boxWidth),center-(boxWidth,boxWidth)]).T,
                        np.array([center+(boxWidth,-boxWidth),center-(boxWidth,-boxWidth)]).T]
        else:
            disabled = None
        return loc1,probeName,loc2,'%4.2f"'%fwhm,disabled
    
    def _guide_locations(self):
        """Makes a list of approximate fiber locations in the guider view."""
        self.guiderLocations = {} # fiber number: location pairs
        # two (near) vertical rows of 5 each
        self.guiderLocations[1] = [228,449]
        self.guiderLocations[8] = [229,362]
        self.guiderLocations[5] = [228,275]
        self.guiderLocations[6] = [228,189]
        self.guiderLocations[4] = [227,103]
        self.guiderLocations[16] = [306,449]
        self.guiderLocations[9]  = [306,362]
        self.guiderLocations[10] = [306,274]
        self.guiderLocations[12] = [306,189]
        self.guiderLocations[14] = [306,101]
        # acquisition
        self.guiderLocations[3] = [129,140]
        self.guiderLocations[11] = [401,137]
        # four others at other locations
        self.guiderLocations[2] = [159,310]
        self.guiderLocations[7] = [92,243]
        self.guiderLocations[13] = [374,400]
        self.guiderLocations[15] = [441,331]
        # make numpy arrays out of all of them
        for x in self.guiderLocations:
            self.guiderLocations[x] = np.array(self.guiderLocations[x])
    
    def __call__(self,outfile,cmap_name=None):
        """Save an image to outfile."""        
        if not cmap_name:
            cmap_name = self.cmap2
        cmap = cm.cmap_d[self.cmap1]
        cmap2 = cm.cmap_d[cmap_name]
        cmap2._init()
        cmap2._lut[0,:] = (0,0,0,1)
        cmap2.set_bad((0,0.15,0)) # mark the masked pixels (not fibers) with dark green
        cmap2.set_under((0.05,0.05,0.05)) # mark everything outside the fibers as dull-gray
        
        fig = plt.figure(figsize=(self.width/self.dpi,self.height/self.dpi),dpi=self.dpi,frameon=False)

        ax1 = plt.Axes(fig,(0,0,.5,1),frame_on=False)
        ax1.set_axis_off()
        fig.add_axes(ax1)
        ax1.imshow(self.guiderView,aspect=self.aspect,origin='lower',cmap=cmap,interpolation='nearest',vmin=0,vmax=1)
        # jkp TBD: scatter doesn't work right for single pixels in matplotlib < 1.1.1
        # see: https://github.com/matplotlib/matplotlib/pull/695
        if self.guiderSat[0].shape > 0:
            ax1.scatter(self.guiderSat[1]+0.5,self.guiderSat[0]+0.5,s=1,marker=',',color='magenta',edgecolor='none',antialiased=False)
        if self.guiderBad[0].shape > 0:
            ax1.scatter(self.guiderBad[1],self.guiderBad[0],s=1,marker=',',color='red')
        ax1.axis((0,512,0,512))
        ax1.text(10,490,'TAI=%s'%self.time,color='white')
        ax1.text(10,470,'frameNo=%s'%self.index,color='white')
        ax1.text(10,450,'seeing=%4.2f'%self.seeing,color='white')
        ax1.text(10,10,'cart=%2d, plate=%d'%(self.cart,self.plate),color='white')
        ax1.autoscale(False)
        for i,label in self.guideLabels.items():
            loc = label[0]
            ax1.text(loc[0],loc[1],label[1],color='white',horizontalalignment='right',fontsize=10)
            loc = label[2]
            ax1.text(loc[0],loc[1],label[3],color='green',horizontalalignment='center',fontsize=10)
            if label[4] != None:
                ax1.plot(label[4][0][0],label[4][0][1],color='red')
                ax1.plot(label[4][1][0],label[4][1][1],color='red')
                
        # Be careful: Axes takes: [left,bottom,width,height]!
        # Note those last two are *not* "right","top"!
        ax2 = plt.Axes(fig,(.5,0,.5,1),frame_on=False)
        ax2.set_axis_off()
        fig.add_axes(ax2)
        ax2.imshow(self.plateView,aspect=self.aspect,origin='lower',cmap=cmap2,interpolation='nearest',vmin=0,vmax=1)
        ax2.quiver(self.vectors[0],self.vectors[1],self.vectors[2],self.vectors[3],color=(0,1,0),width=0.003,edgecolor='none',headwidth=0,units='width',scale=self.qscale)
        ax2.text(10,490,'offset=%4.3f", %4.3f", %4.3f"'%self.offset,color='white')
        ax2.text(10,470,'focus=%4.2f$\mu$m'%self.focus,color='white')
        # jkp TBD: Do I need to multiply this by 100? Schlegel says yes,
        # but the fits header says its in %...
        ax2.text(10,450,'scale=%4.2e%%'%(self.scale),color='white')
        ax2.text(10,430,'RMSerr=%4.2f'%(self.rms),color='white')
        ax2.text(25,410,'1"',color=(0,1,0),)
        ax2.quiver(10,390,1,0,color=(0,1,0),width=0.003,edgecolor='none',headwidth=0,units='width',scale=self.qscale)
        # orientation
        ax2.arrow(470,450,0,30,facecolor='gray',edgecolor='gray',head_width=5, head_length=5)
        ax2.arrow(470,450,30,0,facecolor='gray',edgecolor='gray',head_width=5, head_length=5)
        ax2.text(470,490,'N',color='gray')
        ax2.text(508,450,'E',color='gray')
        ax2.autoscale(False)
        for i,label in self.plateLabels.items():
            try:
                loc = label[0]
                ax2.text(loc[0],loc[1],label[1],color='white',fontsize=10)
                # jkp: skip the vertical lines for fwhm, as they are distracting.
                #loc = label[2]
                #ax2.quiver(loc[0],loc[1],0,label[3],color='green',width=0.003,edgecolor='none',headwidth=0,units='width',scale=self.qscale)
                if label[4] != None:
                    ax2.plot(label[4][0][0],label[4][0][1],color='red')
                    ax2.plot(label[4][1][0],label[4][1][1],color='red')
            except IndexError:
                continue

        fig.set_size_inches(self.width/self.dpi,self.height/self.dpi)
        plt.savefig(outfile,pad_inches=0,dpi=self.dpi)
        plt.close() # close it, so we don't eat too much RAM.
    #...
#...

def make_images(gimgdir,start,end,tempdir,cmap=None,verbose=False):
    """Generate all jpegs from the processed guider files in gimgdir from start to end."""
    # not guaranteed to have every number from start to end, so we have to count.
    count = 0
    files = []
    for i in range(start,end):
        infile = os.path.join(gimgdir,gimgbase%i)
        if os.path.exists(infile):
            outfile = os.path.join(tempdir,tempbase%count)
            imageMaker = ImageMaker(infile)
            imageMaker(outfile,cmap_name=cmap)
            if verbose:
                print '%s -> %s'%(infile,outfile)
            files.append(outfile)
            count += 1
    return count,files
#...

def make_movie(indir,outfile,framerate,verbose=False):
    """Create the movie with ffmpeg, from files in tempdir."""
    inpath = os.path.join(indir,tempbase)
    # NOTE: the order of ffmpeg arguments *REALLY MATTERS*.
    # Reorder them at your own peril!
    # more notes:
    # profile main and pix_fmt yuv420p are for quicktime compatibility.
    # threads 1 so we don't eat up all the processors on hub25m.
    # b:v 10000k is 10MB/second video rate.
    # ffmpeg verbose levels are listed here:
    # http://superuser.com/questions/326629/how-can-i-make-ffmpeg-be-quieter-less-verbose
    if verbose:
        ffmpeg_verbose = 'info'
    else:
        ffmpeg_verbose = 'error'
    cmd = ['ffmpeg',
           '-v',ffmpeg_verbose,
           '-f','image2',
           '-y',
           '-framerate',framerate,
           '-i',inpath,
           '-vcodec','libx264',
           '-b:v','10000k',
           '-profile:v','main',
           '-pix_fmt','yuv420p',
           '-threads','1',
           '-r',framerate,
           os.path.join(outfile)]
    print 'Running:',r' '.join(cmd)
    subprocess.check_call(r' '.join(cmd),shell=True)
#...

def do_work(opts,gimgdir,start,end):
    """Create temp directory, create images, make movie, and cleanup."""
    retval = 0
    tempdir = tempfile.mkdtemp()
    if opts.verbose:
        print 'Writing image files to:',tempdir
    mjd = gimgdir.split('/')
    # /foo/bar/ splits to have [-1] == '', so check for that
    mjd = mjd[-1] if mjd[-1] != '' else mjd[-2]
    if opts.outdir is not None:
        outdir = opts.outdir
    else:
        outdir = gimgdir
    outfile = os.path.join(outdir,'%s-%04d-%04d.mp4'%(mjd,start,end))
    try:
        import time
        time0 = time.time()
        # uncomment these lines and comment out the normal call to get a profile dump
        #import cProfile
        #prof = cProfile.Profile()
        #count,files = prof.runcall(make_images,gimgdir,start,end,tempdir,cmap=opts.cmap)
        #prof.dump_stats('images.profile')
        count,files = make_images(gimgdir,start,end,tempdir,cmap=opts.cmap,verbose=opts.verbose)
        time1 = time.time()
        if opts.verbose:
            print 'Seconds to make %d pngs: %5.1f'%(count,time1-time0)
            print 'Writing movie to:',outfile
        make_movie(tempdir,outfile,opts.framerate,opts.verbose,)
        time2 = time.time()
        if opts.verbose:
            print 'Seconds to make movie: %5.1f'%(time2-time1)
    except Exception as e:
        retval = -2
        print 'Error producing movie:',e
        import traceback
        traceback.print_exc()
    else:
        print 'Wrote:',outfile
        retval = 0
    finally:
        if opts.verbose:
            print "Cleaning up",tempdir
        shutil.rmtree(tempdir)
        return retval
#...

def main(argv=None):
    from optparse import OptionParser
    if argv is None: argv = sys.argv[1:]
	
    usage = '%prog [OPTIONS] DIR STARTNUM ENDNUM'
    usage += '\n\nGenerate a movie of the guider images in DIR from STARTNUM to ENDNUM.'
    usage += '\nResulting movie is written to DIR as an x264 compressed .mp4.'
    usage += '\n\nDIR example: /data/gcam/56233'
    usage += '\nSTARNUM ENDNUM example: 19 303'
    parser = OptionParser(usage)
    parser.add_option('-r','--framerate',dest='framerate',default='10',
                      help='Frame rate of output video (%default).')
    #parser.add_option('--raw',dest='raw',default='""',
    #                  help='Raw commands to pass on to ffmpeg directly (%default)')
    parser.add_option('--cmap',dest='cmap',default='hsv',
                      help='Colormap used to make the images (log10 scaled) from the fits data (%default).')
    parser.add_option('--outdir',dest='outdir',default=None,
                      help='Directory to write resulting movie to (default to DIR).')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',
                      help='Be verbose with progress (%default).')

    # need options?
    (opts,args) = parser.parse_args(args=argv)

    try:
        gimgdir = args[0]
        start = int(args[1])
        end = int(args[2])
    except IndexError,ValueError:
        parser.error('Need DIR STARTNUM ENDNUM. Pass -h or --help for more information.')
        return -1
	
    return do_work(opts,gimgdir,start,end)
#...

if __name__ == "__main__":
    sys.exit(main())
