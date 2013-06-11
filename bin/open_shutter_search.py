#!/usr/bin/env python

## refinement with two functions
##  open_shutter_search.py 56376 56379 -o "dd.txt"

"""
Search a list of directories on hub25m or sos3  for examples of bad gimg frames.
Assumes any gimg frame with a mean value > 20,000 counts and a 
standard deviation < 2000 counts is potentially bad.
Dumps the list of potential bad frames to a specified file.
"""
import argparse
import sys, os
import glob 
import pyfits, numpy

def prnLine(m, mjd, ndirs='-',  nfiles='-',  nflats='-',  nbad='-', proc='  -'): 
      ss="%4s  %5s  %1s  %4s  %4s   %4s   %3s" % (m, mjd, ndirs, nfiles, nflats, nbad, proc)
      print ss

def ds9topy(p): return [p[1],p[0]]      
def getStat(npData, p,dp): 
    mean=npData[(p[0]-dp):(p[0]+dp),(p[1]-dp):(p[1]+dp)].mean()
    std =npData[(p[0]-dp):(p[0]+dp),(p[1]-dp):(p[1]+dp)].std()
    max=npData[(p[0]-dp):(p[0]+dp),(p[1]-dp):(p[1]+dp)].max()
    return mean, std,max
  
def do_one_dir(m, mjd, outfile, listfiles=False, cutoff=10000):
    """Scan the files in one directory."""
    mdir="/data/gcam/%s/" % (mjd)
    nflats=0;   nbad=0;  nbadList=[]

# check if girectory exist    
    if  os.path.exists(mdir)  != True:
        prnLine(m, mjd )
        return nflats, nbad
    ndirs="y" 
     
# check files
    files = glob.glob(mdir+'/gimg-*.fits')
    nfiles= len(files)
    if nfiles==0:
         prnLine(m, mjd, ndirs=ndirs, nfiles=nfiles) 
         return nflats, nbad
         
    for f  in sorted(files) :
         hdulist=pyfits.open(f,'readonly')
         hdr = hdulist[0].header
         data=hdulist[0].data
         hdulist.close()
         
         imType=hdr['IMAGETYP']
         if imType.strip() != args.imtype:
              continue
         nflats=nflats+1
         npData=numpy.array(data)
#       numrows, numcols = npData.shape
#    flat  npData.shape  == (1024, 1048), other imtypes resize to flat size. 
         if npData.shape != (1024, 1048):
              npData= numpy.resize(npData, (1024, 1048))
         
# area selection for data analysis
         Qp=[84,970]    
         Fp=[459,891]   # for A
         Ap=[459,970]  # top left
         Bp=[615,970]  # top right
         Cp=[459,99]  #  bottom left as A 
         Dp=[615,99]  #  bottom right as B 

# statistics
         dp=10  # square size
         Qmean, Qstd, Qmax =  getStat(npData, ds9topy(Qp), dp)
         Fmean, Fstd, Fmax =  getStat(npData, ds9topy(Fp), dp)
         Amean, Astd, Amax =  getStat(npData, ds9topy(Ap), dp)
         Bmean, Bstd, Bmax =  getStat(npData, ds9topy(Bp), dp)   
         Cmean, Cstd, Cmax =  getStat(npData, ds9topy(Cp), dp)   
         Dmean, Dstd, Dmax =  getStat(npData, ds9topy(Dp), dp)   

     # if bad data 
         def norm(Xmax):  return (Xmax-Qmean)/Qstd 
         Qt=max(norm(Amax), norm(Bmax), norm(Cmax), norm(Dmax))
         if Qt > args.cutoff: 
             qq="-bad-"
             nbad= nbad+1
             nbadList.append([f, Qmean, Qstd, Qt, qq])   
         else: qq=""
          
     # write to file                
         ss1="%s %8.1f %8.1f " % (f, Qmean, Qstd)
         ss2=ss1+"%6.2f  %6.2f  %6.2f  %6.2f " % (norm(Amax), norm(Bmax), norm(Cmax), norm(Dmax))
         ss3=ss2+"%6.2f %s " % (Qt, qq)
         outfile.write("%s \n" % (ss3))
         outfile.flush()    
         
         if f=="/data/gcam/56378/gimg-0665.fits" :
            print f
            print Qmean, Qstd, Qmax
            print Amean, Astd, Amax
            print Bmean, Bstd, Bmax
            print Cmean, Cstd, Cmax
            print Dmean, Dstd, Dmax
            
#   ------- end for 

# print summary for this mjd  
    if nflats==0: 
       prnLine(m, mjd, ndirs=ndirs, nfiles=nfiles, nflats=0)    
    else:  
       proc=nbad*100.0/nflats        
       prnLine(m, mjd, ndirs=ndirs, nfiles=nfiles, nflats=nflats, nbad=nbad, proc="%3i"%proc) 
          
# and print bad files if option -l selected  
    if listfiles: 
       for fbad in nbadList:
            ss="%s   %8.2f   %8.2f  %8.2f  %s " %(fbad[0], fbad[1],fbad[2], fbad[3],fbad[4])
            print "      ",ss

    return nflats, nbad  
   
if __name__ == "__main__":
    desc = 'Check guider flats in /data/gcam/[mjd1-mjd2]'
    usage=' %(prog)s [OPTIONS] mjd1  mjd2'
    parser = argparse.ArgumentParser(description=desc,usage=usage)
    parser.add_argument("mjd1", help="mjd to start", type=int)
    parser.add_argument("mjd2", help="mjd to end", type=int)
    parser.add_argument('-l', '--list', help='print the list of bad files', action='store_true')
    parser.add_argument('-c', '--cutoff', help='''if median of selected portion of frame 
        more then this level, we classify frame as bad''', default=1, type=float)
    parser.add_argument('-o', '--outfile', default="bad_gimg.txt", help='file to save bad frame names')
    parser.add_argument('-i', '--imtype', default="flat", help='flat or object, default flat ')

    args = parser.parse_args()    
    
    nflatsTot=0;   nbadTot=0; 
    mjdsList=sorted([args.mjd1, args.mjd2])
    mjds=range(mjdsList[0], mjdsList[1]+1)
    line = 50*"-"  
    listfiles=args.list
    outfile = open(args.outfile,'w')
    outfile.write("File %s  Dark   Fiber   A    B      Lev    Qty\n" %('-'*27))
    outfile.flush() 

    print line
    print "   i   mjd  dir total  imtype  bad     % "
    print line    
    for m, mjd in enumerate(mjds):
        nflats, nbad= do_one_dir(m, mjd, outfile, listfiles=listfiles, cutoff=args.cutoff)  
        nflatsTot=nflatsTot+nflats
        nbadTot=nbadTot+nbad
    print line
 
    if  nflatsTot == 0:  pers="n/a"
    else:  pers="%s" % int(nbadTot*100.0/nflatsTot)+'%'
    print "Summary for guider data:  %s" % args.imtype
    print "MJDs %s - %s,  cutoff=%s" % (args.mjd1, args.mjd2, args.cutoff)
    print "Sum imtype = %s, Sum bad = %s,  percent = %s"  % (nflatsTot, nbadTot, pers)
    print line
    print ""
#...
