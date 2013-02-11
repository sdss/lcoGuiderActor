#!/usr/bin/env python

import os
import sys
from ftplib import FTP

tcc = 'tcc25m-p.apo.nmsu.edu'

def xferFile(pathname):
    dir, fname = os.path.split(pathname)
    os.chdir(dir)
    try:
        ftp = FTP(tcc)
        ftp.login(user='tcc', passwd='no$$change')
        ftp.cwd('tinst')
        ftp.storlines('STOR %s' % (fname), open(fname))
        ftp.close()
        print("updated %s on %s" % (fname, tcc))
    except:
        print("failed to copy %s to %s" % (fname, tcc))

def main():
    fname = sys.argv[-1]
    xferFile(fname)
#...

if __name__ == "__main__":
    main()

