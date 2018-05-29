#!/usr/bin/env python
# encoding: utf-8
#
# guider_movie_allnight.py
#
# Created by José Sánchez-Gallego on 28 Mar 2017.


from __future__ import absolute_import, division, print_function

import argparse
import glob
import os
import sys
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter

import astropy.io.fits as fits

from guider_movie import do_work


class FakeOps(object):
    """Creates a fake opts object to be sent to guider_movie.do_work."""

    def __init__(self):
        self.verbose = False
        self.outdir = None
        self.cmap = 'hsv'
        self.skipcals = False
        self.framerate = '10'


def get_intervals(frame_list):
    """Returns a list of intervals in a list of numbers.

    See http://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list

    """

    intervals = []
    for kk, gg in groupby(enumerate(frame_list), lambda ii_xx: ii_xx[0] - ii_xx[1]):
        intervals.append(map(itemgetter(1), gg))

    return intervals


def guider_movie_allnight(mjd, ecam=False):
    """Creates guider movies for a whole MJD, by cart. Skips intervals of fewer than 5 frames."""

    opts = FakeOps()

    carts = OrderedDict(((ii, []) for ii in range(21, 26)))

    path = os.path.join('/data', 'gcam' if ecam is False else 'ecam', str(mjd))
    files = sorted(glob.glob(os.path.join(path, 'proc-*.fits.gz')))

    print('Collecting data ...')

    for fn in files:

        frame_no = int(os.path.basename(fn).split('-')[-1].split('.')[0])

        header = fits.getheader(fn)

        if 'CARTID' not in header or header['IMAGETYP'].strip().lower() != 'object':
            continue

        cart_id = header['CARTID']

        if cart_id < 21:
            continue
        carts[cart_id].append(frame_no)

    for cart in carts:
        print('Cart', cart, '...')
        intervals = get_intervals(carts[cart])
        for interval in intervals:
            if len(interval) > 5:
                print('Creating movie for frames', interval[0], 'to', interval[-1])
                do_work(opts, path, interval[0], interval[-1])
            else:
                print('Skipping frames', interval)

    return


def main():

    parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description='Processes a full night of guider '
                                                 'images and creates guider movies for each cart.')

    parser.add_argument('mjd', metavar='MJD', type=int, help='The MJD to process.')
    parser.add_argument('-e', '--ecam', dest='ecam', default=False, action='store_true',
                        help='If set, creates movies from the ecam directory for the MJD.')

    args = parser.parse_args()

    guider_movie_allnight(**vars(args))


if __name__ == '__main__':
    main()
