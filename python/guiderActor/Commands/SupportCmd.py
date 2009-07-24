#!/usr/bin/env python

""" Wrap top-level ICC functions. """

import pdb
import logging
import pprint
import sys
import ConfigParser

import opscore.protocols.validation as validation
import opscore.protocols.keys as keys
import opscore.protocols.types as types

import Commands.CmdSet
from opscore.utility.qstr import qstr

class SupportCmd(Commands.CmdSet.CmdSet):
    """ Wrap commands to the alerts actor"""
    
    def __init__(self, icc):
        Commands.CmdSet.CmdSet.__init__(self, icc)

        #
        # Set the keyword dictionary
        #
        self.keys = keys.KeysDictionary("alertscommands", (1, 1),
                                        keys.Key("cmds", types.String()*(1,None),
                                                 help="A list of command modules to reload"),
                                        )

        keys.CmdKey.setKeys(self.keys)

        self.vocab = [
            ('reload', '[<cmds>]', self.reloadCommands),
            ('reloadConfiguration', '', self.reloadConfiguration),
            ('exitexit', '', self.exitCmd),
            ]
        
    def reloadCommands(self,cmd):
        """ If cmds defined, define the listed commands, other wise reload all command sets. """

        if 'cmds' in cmd.cmd.keywords:
            # Load the specified module
            commands = cmd.cmd.keywords['cmds'].values
            for command in commands:
                cmd.respond('text="Attaching %s."' % (command))
                self.icc.attachCmdSet(command)
        else:
            # Load all
            cmd.respond('text="Attaching all command sets."')
            self.icc.attachAllCmdSets()

        cmd.finish('')
    
    def reloadConfiguration(self,cmd):
        """ Reload the configuration. """
        cmd.respond('text="Reparsing the configuration file."')
        logging.warn("reading config file %s", self.icc.configFile)
        self.icc.config = ConfigParser.ConfigParser()
        self.icc.config.read(self.icc.configFile)
        cmd.finish('')
    
    def exitCmd(self, cmd):
        """ Brutal exit when all else has failed. """
        from twisted.internet import reactor

        reactor.stop()
        sys.exit(0)
