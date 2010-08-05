""" Wrap TCC gcam functions. """

import pdb
import logging
import os, re, sys
import time
import numpy
import threading

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

class TCCCmd(object):
    '''
 For commands from the tcc, we need to support:
    init
    setcam N
    doread S Xbin Ybin Xctr Yctr Xsize Ysize
    findstars N Xctr Yctr Xsize Ysize XpredFWHM YpredFWHM
    
findstars            1      171.0      171.0     1024.0     1024.0        3.5        3.5
  yields:
i RawText="findstars            1      171.0      171.0     1024.0     1024.0        3.5        3.5"
i RawText="3 3   213.712 144.051   5.73 4.90 77.5   192.1 5569.1 328.0   0.008 0.008   0"
: RawText=" OK"

  but must be preceded by a doread:

doread       8.00     3     3      171.0      171.0     1024.0     1024.0
  
    '''

    def __init__(self, actor):
        self.actor = actor

        self.keys = keys.KeysDictionary("guider_tccgcam", (1, 1))
        self.vocab = [
            ("tccCmd", "@raw", self.tccCmd),
            ("OK", "", self.OK)
            ]

        # These are the TCC command we implement, and are the 1st word in
        # the @raw tccCmd argument
        self.tccCmds = dict(init=self.doTccInit,
                            setcam=self.doTccSetcam,
                            doread=self.doTccDoread,
                            showstatus=self.doTccShowstatus,
                            findstars=self.doTccFindstars)

        # The ID the TCC expects us to be. We pick up the real value
        # from the setcam N command
        self.GCamId = -1

        self.GImName = "AltaE6ISwear"
        self.size = [1024, 1024]
        
    def echoToTcc(self, cmd, ret):
        """ Pass all the ret lines back to the tcc. """

        for i in range(len(ret)-1):
            cmd.respond('txtForTcc=%s' % (qstr(ret[i])))
        cmd.finish('txtForTcc=%s' % (qstr(ret[-1])))

    def echoCmdToTcc(self, cmd):
        """ echo our command back to the tcc. """

        cmd.respond('txtForTcc=%s' % (qstr(cmd.cmd.keywords['raw'].values[0])))

    def splitTccCmd(self, cmd):
        rawCmd = cmd.cmd.keywords['raw'].values[0]
        return rawCmd.split()
    
    def tccCmd(self, cmd):
        """ Accept a raw tcc command, saved in the @raw argument. """

        tccCmdStr = self.splitTccCmd(cmd)
        if len(tccCmdStr) == 0:         # ignore blank lines.
            return;
        tccCmd = tccCmdStr[0]
        
        if not tccCmd in self.tccCmds:
            cmd.inform('txtForTcc="UNIMPLEMENTED gcam command: %s"' % (tccCmd))
            cmd.warn('text="UNIMPLEMENTED gcam command: %s"' % (tccCmd))
            self.OK(cmd)

        self.echoCmdToTcc(cmd)
        self.tccCmds[tccCmd](cmd)
        
    def OK(self, cmd):
        """ A TOTAL hack, for recovering after we fail in the middle of a TCC command. This
        merely generates a fake completion for the TCC. Even though we say OK, the command 
        has certainly failed, and the TCC will recognize this -- because it has not 
        seen the expected command response. 
        """
        
        cmd.finish('txtForTcc=" OK"')
    
    def doTccInit(self, cmd):
        """ Clean up/stop/initialize ourselves. """

        cmd.finish('txtForTcc="OK"')
        
    def doTccSetcam(self, cmd):
        ''' Respond to a tcc 'setcam N' command.

        This is intended to follow an 'init' command, and truly configures a GImCtrl camera. We,
        however, do not need it, so we gin up a fake response.

        A sample response for the NA2 camera is:
        setcam 1
        1 \"PXL1024\" 1024 1024 16 -11.11 244 \"camera: ID# name sizeXY bits/pixel temp lastFileNum\"
         OK
        '''

        # Parse off the camera number:
        tccCmd = self.splitTccCmd(cmd)
        gid = int(tccCmd[-1])
        self.GImCamID = gid
        
        lastImageNum = 'nan' # self.camera.lastImageNum() # 'nan'
        
        cmd.respond('txtForTcc=%s' % (qstr('%d "%s" %d %d %d nan %s "%s"' % \
                                           (gid, self.GImName,
                                            self.size[0], self.size[1], 16,
                                            lastImageNum,
                                            "camera: ID# name sizeXY bits/pixel temp lastFileNum"))))
        cmd.finish('txtForTcc=" OK"')

    def doTccDoread(self, cmd):
        """ doread S Xbin Ybin Xctr Yctr Xsize Ysize

        e.g.:
        doread       8.00     3     3      171.0      171.0     1024.0     1024.0
        """

        tccCmd = self.splitTccCmd(cmd)
        try:
            itime = float(tccCmd[1])
            xbin, ybin, xCtr, yCtr, xSize, ySize = tccCmd[2:]
        except Exception, e:
            cmd.warn('text="doTccRead barfed on %s: %s"' % (tccCmd, e))
            self.OK(cmd)
            return

        cmd.warn('text="doTccDoRead itime=%01.f"' % (itime))
        self.OK(cmd)

    def doTccFindstars(self, cmd):
        """ findstars N Xctr Yctr Xsize Ysize XpredFWHM YpredFWHM
    
        e.g.:
        findstars            1      171.0      171.0     1024.0     1024.0        3.5        3.5
        
        yields:

i RawText="findstars            1      171.0      171.0     1024.0     1024.0        3.5        3.5"
i RawText="3 3   213.712 144.051   5.73 4.90 77.5   192.1 5569.1 328.0   0.008 0.008   0"
: RawText=" OK"
        """

        tccCmd = self.splitTccCmd(cmd)
        try:
            itime = float(tccCmd[1])
            xCtr, yCtr, xSize, ySize, xPred, yPred = map(float(tccCmd[2:]))
        except Exception, e:
            cmd.warn('text="doTccRead barfed on %s: %s"' % (tccCmd, e))
            self.OK(cmd)
            return

        cmd.warn('text="doTccFindstars not really implemented"')
        cmd.warn('txtForTcc="%d %d %0.3f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.3f %0.3f %d' % \
                     (3, 3,   213.712, 144.051,   5.73, 4.90, 77.5,   192.1, 5569.1, 328.0,   0.008, 0.008,   0))
        self.OK(cmd)
        
    def doTccShowstatus(self, cmd):
        ''' Respond to a tcc 'showstatus' command.
        
showstatus
1 "PXL1024" 1024 1024 16 -26.02 2 "camera: ID# name sizeXY bits/pixel temp lastFileNum"
1 1 0 0 0 0 nan 0 nan "image: binXY begXY sizeXY expTime camID temp"
8.00 1000 params: boxSize (FWHM units) maxFileNum
 OK

        '''
        temp = 0.0    # self.camera.cam.read_TempCCD(),
        cmd.respond("txtForTcc=%s" % (qstr('%d "%s" %d %d %d %0.2f nan "%s"' % \
                                           (self.GImCamID, self.GImName,
                                            self.size[0], self.size[1], 16,
                                            temp,
                                            "camera: ID# name sizeXY bits/pixel temp lastFileNum"))))
        time.sleep(0.1)
        cmd.respond('txtForTcc=%s' % (qstr('%d %d %d %d %d %d nan 0 nan "%s"' % \
                                           (1, 1, 0, 0, 0, 0,
                                            "image: binXY begXY sizeXY expTime camID temp"))))
        time.sleep(0.1)
        cmd.respond('txtForTcc=%s' % (qstr('8.00 1000 "%s"' % \
                                           ("params: boxSize (FWHM units) maxFileNum"))))
        self.OK(cmd)
        
