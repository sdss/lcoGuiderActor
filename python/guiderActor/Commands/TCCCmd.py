""" Wrap TCC gcam functions. """

import pdb
import logging
import os, re, sys
import time
import numpy
import threading
import pyfits

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

# I hate this global stuff...
from guiderActor import Msg
import guiderActor
import guiderActor.myGlobals as myGlobals

import gimg.findsinglestar as findstar
reload(findstar)

cameraNames = {3:'ecamera',
               6:'gcamera'}

class TCCShimState(object):
    def __init__(self):
        self.reset()


        ###### HACK CAMID=3 UNTIL TCC FIX
    def reset(self, camID=3, lastFileNum='nan'):
        # The ID the TCC expects us to be. We pick up the real value
        # from the setcam N command
        self.gImCamID = camID
        self.camera = cameraNames[camID]
        self.gImName = "AltaE6ISwear+%s" % (self.camera)
        self.ccdSize = [1024, 1024]

        self.doreadFilename = None
        self.lastFileNum = lastFileNum

        self.itime = 1.0
        self.binX = self.binY = 2
        self.ctrX = (self.ccdSize[0]+0.5)/self.binX / 2
        self.ctrY = (self.ccdSize[1]+0.5)/self.binY / 2
        self.sizeX = self.ccdSize[0] / self.binX
        self.sizeY = self.ccdSize[1] / self.binY

    def startExposure(self, itime, binX, binY, ctrX, ctrY, sizeX, sizeY):
        self.doreadFilename = None
        self.itime = itime
        self.binX = binX
        self.binY = binY
        self.ctrX = ctrX
        self.ctrY = ctrY
        self.sizeX = sizeX if (sizeX > 0) else self.ccdSize[0]/binX
        self.sizeY = sizeY if (sizeY > 0) else self.ccdSize[1]/binY
        
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
        
        self.tccState = TCCShimState()

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
        
    def OK(self, cmd, fail=False):
        """ A TOTAL hack, for recovering after we fail in the middle of a TCC command. This
        merely generates a fake completion for the TCC. Even though we say OK, the command 
        has certainly failed, and the TCC will recognize this -- because it has not 
        seen the expected command response. 
        """
    
        if fail:
            cmd.warn('text=%s' % (qstr(fail)))
        cmd.finish('txtForTcc=" OK"')
    
    def notOK(self, cmd, why):
        self.OK(cmd, fail=why)
    
    def doTccInit(self, cmd):
        """ Clean up/stop/initialize ourselves. """

        self.tccState.reset()
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

        tccState = self.tccState
        tccState.reset(camID=gid)
        
        cmd.respond('txtForTcc=%s' % (qstr('%d "%s" %d %d %d nan %s "%s"' % \
                                           (gid, self.tccState.gImName,
                                            tccState.ccdSize[0], tccState.ccdSize[1], 16,
                                            tccState.lastFileNum,
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
            binX, binY = map(int, tccCmd[2:4])
            ctrX, ctrY, sizeX, sizeY = map(float, tccCmd[4:])
        except Exception, e:
            cmd.warn('text="doTccRead barfed on %s: %s"' % (tccCmd, e))
            self.OK(cmd)
            return

        self.tccState.startExposure(itime, binX, binY, ctrX, ctrY, sizeX, sizeY)
        myGlobals.actorState.queues[guiderActor.MASTER].put(Msg(Msg.TCC_EXPOSURE, cmd=cmd,
                                                                expTime=itime,
                                                                forTCC=self.tccState, 
                                                                camera=self.tccState.camera))

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
            ctrX, ctrY, sizeX, sizeY, predX, predY = map(float, tccCmd[2:])
        except Exception, e:
            self.notOK(cmd, "doTccRead barfed on %s: %s" % (tccCmd, e))
            return
        
        tccState = self.tccState
        if not tccState.doreadFilename:
            self.notOK(cmd, "no doread file available to process!")
            return

        cmd.diag('text="doTccFindstars processing %s"' % (tccState.doreadFilename))

        try:
            img = pyfits.getdata(tccState.doreadFilename)
        except Exception, e:
            self.notOK(cmd, "could not read image file %s: %s" % (tccState.doreadFilename, e))
            return

        try:
            cmd.inform('file="%s","%s"' % os.path.split(tccState.doreadFilename))
        except Exception, e:
            cmd.warn('text="failed to send file key: %s"' % (e))

        try:
            # Subframe before sending -- CPL
            star = findstar.find_single_star(img)
        except Exception, e:
            self.notOK(cmd, "could not find a star in %s: %s" % (tccState.doreadFilename, e))
            return

        if not star:
            self.notOK(cmd, "could not find a star in %s" % (tccState.doreadFilename))
            return

        starX, starY = star['centroid']
        cmd.warn('txtForTcc="%d %d %0.3f %0.3f %0.2f %0.2f %0.1f %0.1f %0.1f %0.1f %0.3f %0.3f %d' % \
                     (tccState.binX, tccState.binY,
                      starX, starY,  
                      2.5, 2.5, 77.5,   192.1, 5569.1, 328.0,   0.008, 0.008,   0))
        
        self.OK(cmd)
        
    def doTccShowstatus(self, cmd):
        ''' Respond to a tcc 'showstatus' command.
        
showstatus
1 "PXL1024" 1024 1024 16 -26.02 2 "camera: ID# name sizeXY bits/pixel temp lastFileNum"
1 1 0 0 0 0 nan 0 nan "image: binXY begXY sizeXY expTime camID temp"
8.00 1000 params: boxSize (FWHM units) maxFileNum
 OK

        '''

        tccState = self.tccState
        temp = 0.0    # self.camera.cam.read_TempCCD(),
        cmd.respond("txtForTcc=%s" % (qstr('%d "%s" %d %d %d %0.2f %s "%s"' % \
                                           (tccState.gImCamID, tccState.gImName,
                                            tccState.ccdSize[0], tccState.ccdSize[1], 16,
                                            temp, tccState.lastFileNum,
                                            "camera: ID# name sizeXY bits/pixel temp lastFileNum"))))
#        time.sleep(0.4)
        cmd.respond('txtForTcc=%s' % (qstr('%d %d %0.1f %0.1f %0.1f %0.1f %d %0.1f "%s"' % \
                                               (tccState.binX, tccState.binY,
                                                tccState.ctrX, tccState.ctrY,
                                                tccState.sizeX, tccState.sizeY,
                                                tccState.itime, tccState.gImCamID, temp,
                                                "image: binXY begXY sizeXY expTime camID temp"))))
#        time.sleep(0.4)
        cmd.respond('txtForTcc=%s' % (qstr('8.00 9999 "%s"' % \
                                           ("params: boxSize (FWHM units) maxFileNum"))))
        self.OK(cmd)
        
