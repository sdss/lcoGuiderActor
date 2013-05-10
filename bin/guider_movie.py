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
import numpy as np

gimgbase = 'proc-gimg-%04d.fits'
tempbase = 'temp-gimg-%04d.png'
width = 512*2+10 # 10px buffer
height = 512
dpi = 72

def asinh(inputArray, scale_min=None, scale_max=None, non_linear=2.0):
	"""
    Performs asinh scaling of the input numpy array.
    taken from:
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
    data = pyfits.open(infile)
    if not cmap:
        cmap = 'hsv'
    cmap = cm.cmap_d[cmap]
    fig = plt.figure(figsize=(width/dpi,height/dpi))
    axg = fig.add_subplot(121)
    axg.set_axis_off()    
    axg.imshow(asinh(data[0].data),aspect='equal',origin='lower',cmap=cmap)
    axp = fig.add_subplot(212)
    plt.savefig(outfile,bbox_inches='tight',dpi=dpi)
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
    #'-v','panic',
    # profile main and yuv420p are for quicktime compatibility.
    cmd = ['ffmpeg',
           '-f','image2',
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
    mjd = os.path.split(gimgdir)
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
        #shutil.rmtree(tempdir)
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
