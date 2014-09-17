#!/usr/bin/env python
"""Extract some header values from the past for analysis."""

import glob
import argparse
import os
import sys

import fitsio
#from astropy.io import fits
import pyfits
import numpy as np

class Headers(object):
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.keys = ['AZ',
                     'ALT',
                     'IPA',
                     'SPA',
                     'ARCOFFX',
                     'ARCOFFY',
                     'GUIDOFFR',
                     'DRA',
                     'DDEC',
                     'DROT',
                     'DFOCUS',
                     'DSCALE',
                     'OFFRA',
                     'OFFDEC',
                     'OFFROT',
                     'OFFFOCUS',
                     'OFFSCALE',
                     'GDRMS',
                     'NGDRMS',
                     'GDXRMS',
                     'GDYRMS',
                     'FILENAME',
                     'DATE-OBS']

    def __call__(self, files, outfilename):
        """Extract the headers defined in init, and write the result as FITS."""
        self.get_all(files)
        self.write(outfilename)

    def get_one(self, header):
        """Return the requested values from one file header."""
        result = []
        for key in self.keys:
            result.append(header[key])
        return result

    def get_all(self, files):
        """Get the headers from all requested files."""
        result = []
        for f in files:
            header = fitsio.read_header(f, 0)
            #header = pyfits.getheader(f)
            # don't want darks or flats
            if 'object' in header['IMAGETYP']:
                if self.verbose:
                    print "Extracting:", f
                line = self.get_one(header)
                result.append(line)
            else:
                if self.verbose:
                    print "Ignoring:", f
        self.data = np.rec.fromrecords(result,names=self.keys)

    def write(self, filename):
        """Write the result to a fits file."""
        #hdu = pyfits.BinTableHDU(self.data)
        #hdu.writeto(filename)
        hdu = fitsio.write(filename, self.data)
        if self.verbose:
            print "Wrote to:", filename

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog=os.path.basename(sys.argv[0]))
    parser.add_argument('FILEGLOB', metavar='FILEGLOB', type=str,
                        help='A quoted glob regex for the files to process.')
    parser.add_argument('OUTFILE', metavar='OUTFILE', type=str,
                        help='The output fits filename to write to.')
    parser.add_argument('-v','--verbose',action='store_true',dest='verbose',
                        help='Print lots of extra output.')

    args = parser.parse_args()

    files = glob.glob(args.FILEGLOB)
    if files == []:
        print "Error: no files found in glob %s"%args.FILEGLOB
        return
    headers = Headers(verbose=args.verbose)
    headers(sorted(files), args.OUTFILE)

if __name__ == "__main__":
    main()
