#!/usr/bin/env ipython

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
    if len(cartInfo) != 1 or len(probeInfo) != 1:
        import pdb; pdb.set_trace()
            
    mmToDeg = 1/217.7358
    pcen = numpy.array([cartInfo.xcen, cartInfo.ycen])
    prad = cartInfo.radius * math.sqrt(2.0)
    print "! radius = %0.1f" % (cartInfo.radius)
    print "GProbe  %d T" % (probeInfo.fiberId)
    print "    %0.1f %0.1f" % (pcen[0], pcen[1])
    print "    %0.1f %0.1f" % (pcen[0]-prad, pcen[1]-prad)
    print "    %0.1f %0.1f" % (pcen[0]+prad, pcen[1]+prad)
    print "    %0.9f %0.9f" % (-probeInfo.yFocal * mmToDeg, probeInfo.xFocal * mmToDeg)
    print "    0.0 ! the 25m guider code handles the fiber rotation"
    print ""

def cvtPlugMap(plugFile):
    ypm = YPF.YPF(plugFile)
    pm = ypm.structs['PLUGMAPOBJ'].asArray()

    gfibers = pm[numpy.where((pm.holeType == "GUIDE") & (pm.objType == "NA"))]

    cartId = ypm.vars['cartridgeId'].value
    cartInfo = getCartInfo(cartId)

    # Find the best acquisition probe. CPL
    ptErrProbe = 11
    #import pdb; pdb.set_trace()
    print "PtErrProbe %d" % (ptErrProbe)
    print


    print "! cart = %d; plate = %d" % (cartId, ypm.vars['plateId'].value)
    print "! n(map probes) = %d; n(cart probes) = %d" % (len(gfibers), len(cartInfo))

    for i in range(len(gfibers)):
        f = gfibers[i:i+1]
        probeInfo = cartInfo[numpy.where(cartInfo.gProbeId == f.fiberId)]
        writeGProbe(probeInfo, f)

def getCartInfo(cartID):
    yci = YPF.YPF(os.path.join(os.environ['GUIDERACTOR_DIR'], 'etc', 'gcamFiberInfo.par'))
    allCartInfo = yci.structs['GPROBE'].asArray()
    
    cartInfo = allCartInfo[numpy.where(allCartInfo.cartridgeId == cartID)]
    return cartInfo

def main():
    plugfile = sys.argv[2]
    
    cvtPlugMap(plugfile)

instblock_template = """
!+
! SDSS Engineering Spectrograph.
! The guide probe data must be specified separately.
!-

! Instrument

! image center (x,y unbinned pix) and scale (x,y unbinned pix/deg on sky)
! center is desired loc. of zero boresight; it need not be geometric center
IIm_Ctr         384       256
IIm_Scale     23722     23722
! minimum and maximum instrument image coordinates (x,y unbinned pix)
IIm_MinXY         0         0
IIm_MaxXY       767       511

! focus (secondary piston) offset due to instrument (um)
Inst_Foc          1400                 !Changed 000329 DL
! position of the center of the rotator in instrument frame (x,y deg)
Rot_Inst_xy       0.0         0.0      !Changed 040421 DL
! angle from the instrument x axis to the rotator x axis (deg)
Rot_Inst_ang    -89.938                !Changed 060825 DL
!Restored after rotator offset corrected

! minimum and maximum rotator angle due to the instrument (deg)
! note: the actual rotator limits are the intersection of this
! with the RotLim angle limits
InstPosLim     -999       999
! Guide Camera Image

! ID: controller number (0 if no controller, <0 invalid)
! and device number (controller-specific)
! (see file config.txt on the guider Macintosh for camera
! definitions)
! 20020606 PRN Changed from device#2 to device#4 to get
!              the old 40e-/ADU inverse gain.
GCamID            1         6  !set to high gain 20021113 DL
! view number
GCViewNum         0
! maximum counts per pixel
GIm_MaxCount   65535
! default center of guide image (unbinned pix)
! (for GCAM commands; used when user does not specify center)
GIm_Ctr         384       256
! guide image scale (unbinned pix/deg on sky)
! 60.6 microns/arcsec, 9 micron pixels
GIm_Scale     23722     23722
! minimum and maximum position on guide image (unbinned pix)
GIm_MinXY         0         0
GIm_MaxXY       767       511

! binning factor (integer x,y)
GIm_BinFac        2         2

! Guide Probe(s)
! specify this data in a view file, since it varies by plug-plate

! number of probe used for pointing error measurements
PtErrProbe   16
! the minimum # of guide stars needed before rotation
! or scale will be corrected; values <2 are treated as 2
MinStarsForRot  3

! Guide Camera Mechanical
! (focus and filters)

! ID: controller number (0 if no controller, <0 invalid)
! and device number (controller-specific)
GMechID           0         0
! nominal focus setting (piston) for this camera (um)
! net focus = nominal focus + user focus + filter-specific focus correction
GCNomFocus        0
! mininum and maximum guide camera focus (piston) (um)
GCFocLim          0         0
! index to currently selected guide camera filter
! must be in range [1,GCNFilt] if gcNFilt > 0, else ignored
! warning: index is the position in the list of filters, it is not the filter number
GCCurrFiltInd               0
! number of guide camera filters available
! must be in range [0,12]; 0 if no filter, 1 if one fixed filter
GCNFilt                     0
"""

if __name__ == "__main__":
    main()

