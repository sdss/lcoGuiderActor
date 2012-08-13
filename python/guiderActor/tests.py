import Queue
import math
import numpy
import os

from guiderActor import *
import guiderActor.myGlobals as myGlobals
#from guiderActor.Commands.GuiderCmd import GprobeInfo
import masterThread

from YPF import YPF

# Generic duck-type object
class ducky(object):
    pass

def mycall(actor=None, forUserCmd=None, cmdStr=None):
    print 'call: actor=', actor, 'forUserCmd=', forUserCmd, 'cmdStr=', cmdStr
    r = ducky()
    r.didFail = False
    return r

def mydiag(s):
    print 'diag:', s

def mywarn(s):
    print 'warn:', s

def myinform(s):
    print 'inform:', s

def myfail(s):
    print 'fail:', s

def myrespond(s):
    print 'respond:', s

def getGprobes(fiberinfofn, plugmapfn, cartridge):
    par = YPF(fiberinfofn)
    probes = par.structs['GPROBE'].asObjlist()
    probes = [p for p in probes if p.cartridgeId == cartridge]
    guideprobes = [p for p in probes if p.fiberType == 'GUIDE' or p.fiberType == 'ACQUIRE']
    tritiumprobes = [p for p in probes if p.fiberType == 'TRITIUM']

    plugmapfile = YPF(plugmapfn)
    plugmap = plugmapfile.structs['PLUGMAPOBJ'].asObjlist()
    aholes = [p for p in plugmap if p.holeType == 'ALIGNMENT']
    gholes = [p for p in plugmap if p.holeType == 'GUIDE']

    print '%i probes' % len(probes)
    print '  %i guide/capture probes' % len(guideprobes)
    print '  %i tritium probes' % len(tritiumprobes)
    print '%i alignment holes' % len(aholes)
    print '%i guide holes' % len(gholes)

    gprobes = {}
    for (p, ghole, ahole) in zip(guideprobes + tritiumprobes, gholes + [None], aholes + [None]):
        info = ducky()
        info.enabled = True
        info.exists = p.exists
        info.fiber_type = p.fiberType
        info.flags = 0
        info.xCenter = p.xcen
        info.yCenter = p.ycen
        info.radius = p.radius
        info.xFerruleOffset = p.xferruleOffset
        info.yFerruleOffset = p.yferruleOffset
        info.rotation = p.rot
        info.focusOffset = p.focusOffset
        if p.fiberType in ['GUIDE', 'ACQUIRE']:
            info.ra  = ghole.ra
            info.dec = ghole.dec
            info.xFocal = ghole.xFocal
            info.yFocal = ghole.yFocal
            info.phi = 90 - math.atan2(ahole.yFocal - ghole.yFocal,
                                       ahole.xFocal - ghole.xFocal) * 180/math.pi
            info.mag = ghole.mag

        else:
            info.ra = 0
            info.dec = 0
            info.xFocal = 0
            info.yFocal = 0
            info.phi = 0
            info.rotStar2Sky = numpy.nan


        gprobes[p.gProbeId] = info
    return gprobes


def getModels():
    models = {}
    tcc = ducky()
    tcc.keyVarDict = {}
    models['tcc'] = tcc
    mcp = ducky()
    mcp.keyVarDict = {'ffsStatus':[[8,0]]}
    models['mcp'] = mcp
    return models


def getCommand():
    cmd = ducky()
    cmd.inform = myinform
    cmd.warn = mywarn
    cmd.fail = myfail
    cmd.respond = myrespond
    return cmd

if __name__ == '__main__':
    os.environ['GUIDERACTOR_DIR'] = '..'

    mq = Queue.Queue(0)
    queues = {MASTER: mq,
              GCAMERA: Queue.Queue(0),
              }

    actor = ducky()
    actor.cmdr = ducky()
    actor.cmdr.call = mycall
    actor.bcast = ducky()
    actor.bcast.warn = mywarn
    actor.bcast.diag = mydiag

    myGlobals.actorState = ducky()
    myGlobals.actorState.timeout = 0
    myGlobals.actorState.models = getModels()
    myGlobals.actorState.tccState = ducky()
    myGlobals.actorState.tccState.halted = False
    myGlobals.actorState.tccState.goToNewField = False

    cmd = getCommand()

    # Simple test
    #m = Msg(Msg.EXIT, cmd)
    #mq.put(m)
    #masterThread.main(actor, queues)

    cartridge = 13

    gprobes = getGprobes('../etc/gcamFiberInfo.par',
                         '../testfiles/plPlugMapM-3615-55201-09.par',
                         cartridge)
    plate=1

    m = Msg(Msg.LOAD_CARTRIDGE, cmd, cartridge=cartridge, plate=plate, pointing=0,
            fscanMJD=0, fscanID=0, boresight_ra=0, boresight_dec=0,
            gprobes=gprobes)
    mq.put(m)

    m = Msg(Msg.SET_SCALE, cmd,
            plugPlateScale=217.7358,
            gcameraMagnification=1.00015,
            gcameraPixelSize=0.026 * 1e3,
            # FIXME
            dSecondary_dmm=1.0)
    mq.put(m)

    m = Msg(Msg.SET_TIME, cmd, expTime=30)
    mq.put(m)

    m = Msg(Msg.SET_PID, cmd, what='raDec', Kp=0.6, Ti=0, Imax=-1, Td=0.)
    mq.put(m)
    m = Msg(Msg.SET_PID, cmd, what='rot', Kp=0.5, Ti=0, Imax=-1, Td=0.)
    mq.put(m)
    m = Msg(Msg.SET_PID, cmd, what='focus', Kp=0.1, Ti=0, Imax=-1, Td=0.)
    mq.put(m)
    m = Msg(Msg.SET_PID, cmd, what='scale', Kp=0.5, Ti=0, Imax=-1, Td=0.)
    mq.put(m)

    m = Msg(Msg.START_GUIDING, cmd, start=True, force=False, oneExposure=False,
            plot=True, psPlot=True, spiderInstAng=0, expTime=0)
    mq.put(m)

    m = Msg(Msg.EXPOSURE_FINISHED, cmd,
            filename='../testfiles/modified-55205-gimg-0005.fits',
            success=True)
    mq.put(m)

    m = Msg(Msg.EXIT, cmd)
    mq.put(m)

    masterThread.main(actor, queues)

    gq = queues[GCAMERA]
    print 'Gcam queue:', gq
    print 'contents:'
    while not gq.empty():
        m = gq.get()
        print '  ', m
    print 'end of contents'

    print 'Master queue:', mq
    print 'contents:'
    while not mq.empty():
        m = mq.get()
        print '  ', m
    print 'end of contents'

