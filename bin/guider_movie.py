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

import AssembleImage

gimgbase = 'proc-gimg-%04d.fits'
tempbase = 'temp-gimg-%04d.png'
width = 512*2. # 10px buffer
height = 512.
dpi = 100.
assembler = AssembleImage.AssembleImage(1,0)

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

def one_image(infile,outfile,count,cmap=None):
    """Generate a single jpeg from a processed guider fits file."""
    aspect = 'normal' # normal vs. equal?
    index = os.path.splitext(infile)[0].split('-')[-1]
    
    data = pyfits.open(infile)
    plate = data[0].header['PLATEID']
    cart =  data[0].header['CARTID']
    guiderView = asinh(data[0].data,non_linear=5.)
    plateInfo = assembler(data)
    plateView = plateInfo.plateImageArr
    zeros = plateView == 0
    plateView = asinh(plateView,non_linear=5.0)
    plateView[zeros] = 0
    # TBD: should do something to mark the saturated (1) and bad (2) pixels?
    # STUI labels them with colors, but that only works with a grayscale cmap.
    #plateView[plateInfo.plateMaskArr == 1] = ??
    #plateView[plateInfo.plateMaskArr == 2] = ??
    
    #import pdb
    #pdb.set_trace()
    if not cmap:
        cmap = 'hsv_r'
    cmap = cm.cmap_d[cmap]
    cmap2 = cm.gist_ncar
    cmap2._init()
    cmap2._lut[0,:] = (0,0,0,1)
    
    fig = plt.figure(figsize=(width/dpi,height/dpi),dpi=dpi,frameon=False)

    ax1 = plt.Axes(fig,(0,0,.5,1),frame_on=False)
    ax1.set_axis_off()
    fig.add_axes(ax1)
    ax1.imshow(guiderView,aspect=aspect,origin='lower',cmap=cmap)
    ax1.text(20,470,'frameNo=%s'%index)
    ax1.text(20,450,'plate=%d'%plate)
    ax1.text(20,430,'cart=%2d'%cart)
    
    # Careful: Axes takes: [left,bottom,width,height]!
    # Note those last two are *not* "right","top"!
    ax2 = plt.Axes(fig,(.5,0,.5,1),frame_on=False)
    ax2.set_axis_off()
    fig.add_axes(ax2)
    ax2.imshow(plateView,aspect=aspect,origin='lower',cmap=cmap2)
    # TBD: some useful text on this frame?
    #ax2.text(40,400,index)
    
    fig.set_size_inches(width/dpi,height/dpi)
    plt.savefig(outfile,pad_inches=0,dpi=dpi)
    print 'Wrote:',outfile
#...

def make_images(gimgdir,start,end,tempdir,cmap=None):
    """Generate all jpegs from the processed guider files in gimgdir from start to end."""
    # not guaranteed to have every number from start to end, so we have to count.
    count = 0
    files = []
    for i in range(start,end):
        infile = os.path.join(gimgdir,gimgbase%i)
        if os.path.exists(infile):
            outfile = os.path.join(tempdir,tempbase%count)
            one_image(infile,outfile,count)
            files.append(outfile)
            count += 1
    return count,files
#...

def make_movie(opts,indir,outfile):
    """Create the movie with ffmpeg, from files in tempdir."""
    inpath = os.path.join(indir,tempbase)
    # profile main and yuv420p are for quicktime compatibility.
    # remove the '-v panic' to get it to print out all of its encoding steps.
    cmd = ['ffmpeg',
           '-v','fatal',
           '-f','image2',
           '-y',
           '-r',opts.framerate,
           '-i',inpath,
           '-vcodec','libx264',
           '-b:v','10000k',
           '-profile:v','main',
           '-pix_fmt','yuv420p',
           '-threads','1',
           os.path.join(outfile)]
    print 'Running:',r' '.join(cmd)
    subprocess.check_call(r' '.join(cmd),shell=True)
#...

def do_work(opts,gimgdir,start,end):
    """Create temp directory, create images, make movie, and cleanup."""
    tempdir = tempfile.mkdtemp()
    print 'Writing image files to:',tempdir
    mjd = gimgdir.split('/')
    # /foo/bar/ splits to have [-1] == '', so check for that
    mjd = mjd[-1] if mjd[-1] != '' else mjd[-2]
    outfile = '%s-%04d-%04d.mp4'%(mjd,start,end)
    try:
        import time
        time0 = time.time()
        # uncomment these lines and comment out the normal call to get a profile dump
        #import cProfile
        #prof = cProfile.Profile()
        #count,files = prof.runcall(make_images,gimgdir,start,end,tempdir,cmap=opts.cmap)
        #prof.dump_stats('images.profile')
        count,files = make_images(gimgdir,start,end,tempdir,cmap=opts.cmap)
        time1 = time.time()
        print 'Seconds to make %d pngs: %5.1f'%(count,time1-time0)
        print 'Writing movie to:',outfile
        make_movie(opts,tempdir,outfile)
        time2 = time.time()
        print 'Seconds to make movie: %5.1f'%(time2-time1)
    except Exception as e:
        print 'Error producing movie:',e
        import traceback
        traceback.print_exc()
    finally:
        print "Cleaning up",tempdir
        shutil.rmtree(tempdir)
        pass
#...

def main(argv=None):
    from optparse import OptionParser
    if argv is None: argv = sys.argv[1:]
	
    usage = '%prog [OPTIONS] DIR STARTNUM ENDNUM'
    usage += '\n\nGenerate a movie of the guider images in DIR from STARTNUM to ENDNUM.'
    usage += '\nDefaults to generating the images in the guider  frame.'
    usage += '\n\nDIR example: /data/gcam/56233'
    usage += '\nSTARNUM ENDNUM example: 19 303'
    parser = OptionParser(usage)
    parser.add_option('-r','--framerate',dest='framerate',default='10',
                      help='Frame rate of output video (%default).')
    parser.add_option('--raw',dest='raw',default='""',
                      help='Raw commands to pass on to ffmpeg directly (%default)')
    parser.add_option('--cmap',dest='cmap',default='hsv',
                      help='Colormap used to make the images (log10 scaled) from the fits data (%default).')

    # need options?
    (opts,args) = parser.parse_args(args=argv)

    try:
        gimgdir = args[0]
        start = int(args[1])
        end = int(args[2])
    except IndexError,ValueError:
        parser.error('Need DIR STARTNUM ENDNUM. Pass -h or --help for more information.')
	
    do_work(opts,gimgdir,start,end)
#...

if __name__ == "__main__":
    sys.exit(main())
