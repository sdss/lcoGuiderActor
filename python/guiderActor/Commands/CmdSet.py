#!/usr/bin/env python

""" CmdSet.py -- a vocabulary of MC commands. """

class CmdSet(object):
    def __init__(self, icc):

        self.icc = icc

        self.help = ()
        self.vocab = {}
        self.keys = ()
