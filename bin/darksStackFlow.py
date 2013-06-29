#!/usr/bin/env python

"""
2013-06-27 by EM
darksStackFlow.py 56400 56470 -n 10

Stack dark files for n-days, save them and their images in some directory 

Hystory: 

2013-06-27 EM:  initial version created. 
2013-06-28 EM:  added cycle range(mjd1,mjd2, nfiles); 
       added png images per request. 

"""
import argparse
import sys, os
import glob 
import pyfits, numpy, Image

def  getDatList(mjds):   
   # list of lists of files for each mjd
   datLists = [glob.glob("/data/gcam/%s/dark-????.dat*" % (m)) for m in mjds]
   # join list of mjd-lists to one single list with all *.dat files
   datList=[]; 
   for ll in datLists: datList.extend(ll) 
   return datList

def getDarkList(datList):
   darkList=[]
   for i,f in enumerate(datList):   # read files *.dat  to get darks filenames  
      dd = open(f, "r")
      darkName=dd.readline()
      darkName=darkName[:-1] # remove last \n
      darkName=darkName[9:]   # remove first "dark" word
      dd.close() 
      if not os.path.exists(darkName):  # check if such file exist on disk
         print "no file found", darkName
         continue
      darkList.append(darkName)
   return darkList

def do_one_stack(mjds):
    # create the list of dat files with dark pattern 
    datList=getDatList(mjds)
    if len(datList) == 0: 
        print "     no darks"     
        return False
    else: 
        print "     %s darks"% len(datList) 
        
    # create the list of dark files -darkList 
    darkList=getDarkList(datList)
    if len(darkList) == 0: 
        print " something wrong, no dark files but they shoudl be here " 
        return False

    # stack dark images to one cube
    hdu1= pyfits.PrimaryHDU(numpy.zeros([512,524]))
    for i,dd in enumerate(darkList):   # read fits darks
        hdulist=pyfits.open(dd,'readonly')  # check is this file dark? 
        hdr = hdulist[0].header
        dat=hdulist[0].data  # got numpy array
        hdulist.close()
        if hdr.get('IMAGETYP') != "dark":
            print "File is not dark, something wrong:  %s" % dd
            sys.exit()
        if i == 0:  stack=dat
        else: stack=numpy.dstack([stack,dat])
        hdu1.header.add_comment("added %s" % dd)  # add dark name to header comment 
        if args.list:
            print "      %4i  %s  median=%s " % (i,dd, numpy.median(dat))

# calcumate median image using stack images
    median=numpy.median(stack, axis=2) 

# save median of stack to fits file    
    if args.fits: 
        result="%s+%s.fits"%(mjds[0],args.nfiles)  #  fits file name
        hdu1.data=median
        hdu1.writeto(result,clobber= True,)

#  save mean image to png as JP asked
# Rescale to 0-255 and convert to uint8
    if args.image: 
        image="%s+%s.png"%(mjds[0],args.nfiles)  #  fits file name
        min=1780;  max=1900
        rescaled=((median-min)*255.0/(max-min)).astype(numpy.uint8)
        im = Image.fromarray(rescaled)
        im.save(image)

if __name__ == "__main__":
    desc = 'Stack guider darks /data/gcam/[mjd1]/gimg-*.fits for n nights'
    usage=' %(prog)s [OPTIONS] mjd1 mjd2'
    parser = argparse.ArgumentParser(description=desc,usage=usage)
    parser.add_argument("mjd1", help="mjd1 to start", type=int)
    parser.add_argument("mjd2", help="mjd2 to end", type=int)
    parser.add_argument('-n', '--nfiles', help="number of nights to stack (def=7)",
        default=7, type=int)
    parser.add_argument('-f', '--fits', help="save as fits ", action='store_true')
    parser.add_argument('-i', '--image', help="save as png image", action='store_true')
    parser.add_argument('-l', '--list', help="list files for stack", action='store_true')
    args = parser.parse_args()   
    
    line="-"*60
    print line
    for i,m in enumerate(range(args.mjd1, args.mjd2, args.nfiles)): 
        mjds=range(m, m+args.nfiles)
        print "%3i  mjds=%s" % (i,mjds)
        q=do_one_stack(mjds)
        print line 
    print ""
        
        