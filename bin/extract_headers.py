#!/usr/bin/env python
"""Extract some header values from the past for analysis."""

import argparse
import glob
import os
import sys

import fitsio
import numpy as np
from astropy.table import Table


class Headers(object):
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.keys = [['AZ', float],
                     ['ALT', float],
                     ['IPA', float],
                     ['SPA', float],
                     ['RA', float],
                     ['DEC', float],
                     ['FOCUS', float],
                     ['ARCOFFX', float],
                     ['ARCOFFY', float],
                     ['GUIDOFFR', float],
                     ['DRA', float],
                     ['DDEC', float],
                     ['DROT', float],
                     ['DFOCUS', float],
                     ['DSCALE', float],
                     ['OFFRA', float],
                     ['OFFDEC', float],
                     ['OFFROT', float],
                     ['OFFFOCUS', float],
                     ['OFFSCALE', float],
                     ['GDRMS', float],
                     ['NGDRMS', float],
                     ['GDXRMS', float],
                     ['GDYRMS', float],
                     ['FILENAME', str],
                     ['DATE-OBS', str],
                     ['PLATEID', int], ]

    def __call__(self, files, outfilename):
        """Extract the headers defined in init, and write the result as FITS."""
        self.get_all(files)
        self.write(outfilename)

    def get_one(self, header):
        """Return the requested values from one file header."""
        return [type(header.get(key, np.nan)) for key, type in self.keys]

    def get_all(self, files):
        """Get the headers from all requested files."""
        result = []
        for f in files:
            frameno = int(f.split('-')[-1].split('.')[0])
            header = fitsio.read_header(f, 0)
            #header = pyfits.getheader(f)
            # don't want darks or flats
            if 'object' in header['IMAGETYP']:
                if self.verbose:
                    print 'Extracting:', f
                line = self.get_one(header)
                line.insert(0, frameno)
                result.append(line)
            else:
                if self.verbose:
                    print 'Ignoring:', f
            del header
        names = [k[0] for k in self.keys]
        types = [k[1] for k in self.keys]
        for i, t in enumerate(types):
            if t == str:
                types[i] = '|S40'
        names.insert(0, 'FRAMENO')
        types.insert(0, int)
        # self.data = np.rec.fromrecords(result,dtype=np.dtype(zip(names,types)))
        # self.data = np.array(result,dtype=np.dtype(zip(names,types)))
        self.data = Table(rows=result, names=names, dtype=types)

    def write(self, filename):
        """Write the result to a fits file."""
        # hdu = fits.BinTableHDU(self.data)
        # hdu.writeto(filename, clobber=True)
        # hdu = fitsio.write(filename, self.data, clobber=True)
        self.data.write(filename, format='fits')
        if self.verbose:
            print 'Wrote to:', filename


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog=os.path.basename(sys.argv[0]))
    parser.add_argument('FILEGLOB', metavar='FILEGLOB', type=str,
                        help='A quoted glob regex for the files to process.')
    parser.add_argument('OUTFILE', metavar='OUTFILE', type=str,
                        help='The output fits filename to write to.')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose',
                        help='Print lots of extra output.')

    args = parser.parse_args()

    files = glob.glob(args.FILEGLOB)
    if files == []:
        print 'Error: no files found in glob %s' % args.FILEGLOB
        return
    headers = Headers(verbose=args.verbose)
    headers(sorted(files), args.OUTFILE)


if __name__ == '__main__':
    main()
