#!/usr/bin/env python
"""
Search a list of directories on hub25m for examples of bad gimg frames.
Assumes any gimg frame with a mean value > 20,000 counts and a 
standard deviation < 2000 counts is potentially bad.
Dumps the list of potential bad frames to a specified file.
"""
import glob
import pyfits
from optparse import OptionParser, OptionGroup
import sys

def do_one_dir(directory,outfile):
	"""Scan the files in one directory."""
	files = glob.glob(directory+'/gimg-*.fits')
	for f in sorted(files):
		data = pyfits.open(f)[0].data
		if data.mean() > 20000 and data.std() < 2000:
			print 'Potentially bad:',f
			outfile.write(f+'\n')
			outfile.flush()
#...

def main(argv=None):
	if argv is None: argv = sys.argv[1:]
	usage = '%prog [OPTIONS] DIR1 [DIR2 [DIR3 ...]]'
	usage += '\n\nDIR example: /data/gcam/5623*'
	parser = OptionParser(usage)
	parser.add_option('--outfile',dest='outfile',default='bad_gimg.txt',help='write bad file list to this (%default)')
	# need options?
	(opts,args) = parser.parse_args(args=argv)

	if len(args) == 0:
		parser.error('need to pass some directories')
	else:
		directories = args
	
	outfile = open(opts.outfile,'w')
	for directory in directories:
		print 'processing:',directory
		do_one_dir(directory,outfile)
#...

if __name__ == "__main__":
    sys.exit(main())
