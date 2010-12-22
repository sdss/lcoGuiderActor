#!/usr/bin/env python

import os, sys
import math

import numpy
import opscore.utility.YPF as YPF

"""
typedef struct {
   int objId[5];
   HOLETYPE holeType;
   double ra;
   double dec;
   float mag[5];
   float starL;
   float expL;
   float deVaucL;
   OBJTYPE objType;
   double xFocal;
   double yFocal;
   int spectrographId;
   int fiberId;
   int throughput;
   int primTarget;
   int secTarget;
} PLUGMAPOBJ;
                
"""

"""
 ! Guide Probe(s)
 ! probe number (in range [1,20]), does probe exist? (T=yes,F=no)
 ! if and only if the probe exists then this data must follow in order:
 ! center of guide probe on guide image (x,y unbinned pix)
 ! mininum guide image coordinates (x,y unbinned pix)
 ! maximum guide image coordinates (x,y unbinned pix)
 ! position of guide camera probe w.r.t. rotator (x,y deg on sky)
 ! angle from the guide image x axis to rotator x axis (deg)
 GProbe         1    T  ! number, exists?
                512.000000000000        512.000000000000             ! center (x,y unbinned pix)
                0.000000000000000E+000  0.000000000000000E+000        ! mininum (x,y unbinned pix)
                1024.00000000000        1024.00000000000             ! maximum (x,y unbinned pix)
                -0.223999000000000       1.607000000000000E-003        ! pos. w.r.t. rotator (x,y deg on sky)
                0.630000000000000             ! angle of rot w.r.t. image (deg)
                                                                                              
"""

def writeGProbe(cartInfo, probeInfo):
    """ Write out a single GProbe entry for a GCView block. """

    mmToDeg = 1/217.7358
    pcen = numpy.array([cartInfo['xcen'], cartInfo['ycen']])
    prad = cartInfo['radius'] * math.sqrt(2.0)
    print "! radius = %0.1f" % (cartInfo['radius'])
    print "GProbe  %d T" % (probeInfo['fiberId'])
    print "    %0.1f %0.1f" % (pcen[0], pcen[1])
    print "    %0.1f %0.1f" % (pcen[0]-prad, pcen[1]-prad)
    print "    %0.1f %0.1f" % (pcen[0]+prad, pcen[1]+prad)
    print "    %0.9f %0.9f" % (-probeInfo['yFocal'] * mmToDeg, probeInfo['xFocal'] * mmToDeg)
    print "    0.0 ! the 25m guider code handles the fiber rotation"
    print ""

def cvtPlugMap(plugFile):
    """ Write out an entire GCView block for the given plugfile. """

    ypm = YPF.YPF(plugFile)
    pm = ypm.structs['PLUGMAPOBJ'].asArray()

    gfibers = pm[numpy.where((pm.holeType == "GUIDE") & (pm.objType == "NA"))]

    cartId = ypm.vars['cartridgeId'].value
    cartInfo = getCartInfo(cartId)

    # Find the best acquisition probe. The small acquisition fibers are not marked as such,
    # so use the radius.
    acqProbes = cartInfo[(cartInfo['exists'] > 0) & (cartInfo['radius'] > 14)]
    
    # Use the closest one to the plate center. 
    dMin = 9999.9
    pMin = None
    for p in acqProbes:
        f = gfibers[gfibers['fiberId'] == p['gProbeId']]
        d = numpy.sqrt(f['xFocal']**2 + f['yFocal']**2)
        if d < dMin:
            dMin = d
            pMin = p

    if pMin == None:
        raise RuntimeError("no available acquisition fiber")

    print "! cart = %d; plate = %d" % (cartId, ypm.vars['plateId'].value)
    print
    print "PtErrProbe %d" % (pMin['gProbeId'])
    print

    for f in gfibers:
        probeInfo = cartInfo[numpy.where(cartInfo['gProbeId'] == f['fiberId'])]
        if probeInfo['exists'] > 0: 
            writeGProbe(probeInfo, f)

def getCartInfo(cartID):
    """ Return the per-cartidge guider probe info. """

    yci = YPF.YPF(os.path.join(os.environ['GUIDERACTOR_DIR'], 'etc', 'gcamFiberInfo.par'))
    allCartInfo = yci.structs['GPROBE'].asArray()
    
    cartInfo = allCartInfo[numpy.where(allCartInfo.cartridgeId == cartID)]
    return cartInfo

def main():
    plugfile = sys.argv[-1]
    
    cvtPlugMap(plugfile)


if __name__ == "__main__":
    main()

