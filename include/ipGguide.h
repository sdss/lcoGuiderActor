#ifndef GGUIDE_H
#define GGUIDE_H
/*****************************************************************************
******************************************************************************
**
** FILE:
**	ipGguide.h
**
** ABSTRACT:
**	Headers for Jim Gunn's SOP guider
**
** ENVIRONMENT:
**	ANSI C.
**
** AUTHOR:
**	Creation date:  2002-03-31
**	Jim Gunn
**      Imported into Derish: 2003-01-07
**      Eric Neilsen
** Mods
**      New Apogee camera:   2009-07-29
**      Paul Harding      
******************************************************************************
******************************************************************************
*/

/* header file for new guider package */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "shLegacy.h"


#define STATIC static
#define DOUBLE double

/*the Alta guider is 16 bit, but bitshifted in GCAM */
/* max value of guider data*/
#define DMAX 32767

/* FNORM is the normalization of the flat, so everything works with the integer 
   profile fitting code, by JG. */
/* in GMFLATTEN there is a bitshift of 12, thus FNORM must be 4096 */
#define FNORM DMAX/8			

/* PMAX is the normalization of the profile (central pixel), so everything works with the integer 
   profile fitting code, by JG. */
#define PMAX DMAX/8

/*maximum radius we have in the profile arrays seeprofile and radtable */
/*2" FWHM = sigma of 2 pixels and go to 3 sigma evaluated at 0.1 sigma steps*/  

#define MAXR 60   

/* max radius^2 we have in the profile arrays seeprofile and radtable */
/* number of pixels in profile array */
#define MAXR2 MAXR*MAXR

/* inverse gain of guider ccd */
/* photometrics 22 e-/ADU inverse gain, for low gain, prior to 20 Nov 2002 */
/* photometrics  6.5 e-/ADU inverse gain, for high gain, current as of 20 Nov 2002 */
/* Alta        1.4 e-/ADU inverse gain, current as of 29 Aug 2008*/

#define CCDGAIN 1.4

# define sigp2FwhmAs 1.0 /* conversion from sigma in pixels to fwhm in arcsec */

   
/********* VARIABLES ***********/

/*moved out of ipGguide.c*/

/* assume the CCD is square */
#define UNBINNED_CCD_SIZE 1024
#define CCD_SIZE UNBINNED_CCD_SIZE/2

/* for gextendmask(), the number of points in the fringe circle which must
 * intersect the mask in order for the center point to be set.
 */
#define MASKPTS 1

/* seeing profile table; allocated by gmakeseeprof(). Never freed. */
STATIC short int *seeprofile=NULL; 

/* radius table for indices in seeprofile table. as above, never freed */
STATIC short int *radtable=NULL;

/* mean value of real radius in extraction table */
STATIC int rvalmean[60];
/* weights for above if no edge effects */
STATIC int rvalwt[60];

/* We will find the position of the star by profile fitting on a 9-point 
 * grid which we move about to find the minimum fittin error; the structures 
 * which contain the information for these nine points are are these:
 */
/* STATIC struct gstarfit starfitarr[9]; */
/* STATIC struct gstarfit *starfit[3] = 
            {starfitarr, starfitarr + 3, starfitarr + 6}; */

STATIC struct gstarfit istarfit ;  /* initial answer */            
STATIC struct gstarfit fstarfit;  /* final answer */
STATIC int pgs2derr[9];

/* 2d fitting error matrix */
STATIC int *gs2derr[3] = {pgs2derr,pgs2derr+3,pgs2derr+6};  

/* print buffer */
STATIC char pbuf[256];

/* arrays for holding x,y vals of bad pixels we want to mask.  
 * limited to 2000 entries as we have bigger problems if more than 2%
 * of the pixels are above our threshold for badness 
 */
#define MAXBADPIX 2000


STATIC S16 badrowval[MAXBADPIX];
STATIC S16 badcolval[MAXBADPIX];
STATIC int nbadpix = 0;

/* run through histogram from the top and find peak and median 
   peak of guider fibers is ~1.5% down due to acquisition fibers
   which transmit more than guide fibers. 
   55% of the data down to find median
   70% for the pseudo bias level 

   The 2 acquisition fibers occupy 1.5%
   The 14 guide fibers occupy 1% 
   Wings of fibers bring the total fiber area up to near 5%
*/
#define PEAK_PERCENTILE 0.010
#define MEDIAN_PERCENTILE 0.55
#define   BIAS_PERCENTILE 0.70

/*jeg set bin size on photometrics to break image into 16 bins
ph maintain the bin size rather than the number of bins for Alta camera.
*/ 

#define BINSIZE	16		/*BINSIZE remains constant*/

#define BINTHRESH 500	/* can probably be increased to 500, assuming flat will have values of ~10000 at least*/



#define NBINS CCD_SIZE/BINSIZE		/* should be an integer */

#define MAXPEAKS 50

#define MAXFIB (17 + 1)                 /* Fibre IDs are sometimes 1-indexed */

#define SPOTID (MAXFIB-2)		/*fiber alignment spot lives in, zero indexed */

#define MAXFIBERS 20			/* large enough for all possible cartridges */

/* max number of iterations in refining the center; 
 * we should not need to walk futher than the range limit of 40 pix radius. */

#define MAXGFSITER 40

/* values of 100/sigma^2 used in width fitting, from sigma=10 pix to sigma=1 pix
 * Alta camera, needed 6 bins to get below sigma=1.41 which was set by max wpg=50 */

#define NCELL 31
STATIC short int wpg[NCELL] = 
    {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100};

/********* STRUCTURES **********/

typedef struct platedata{
	int *gProbeId;		/*fiber number, should iterate 1-17 */
	int *exists;		/*does this fiber exist for this plate?*/
	double *xcen;		/*ideal xcen of fibers*/
	double *ycen;		/*ideal ycen of fiber*/
	double *radius;		/*idea radius of fibers*/
} PLATEDATA;

typedef struct ghist_t{
    int ghist_medn;     /* median (55 percentile down) of whole array--essentially the bkgnd */
    int ghist_peak;     /* 1.5 th percentile, peak in guide fibers*/
    int ghist_ref;      /* halfway between medn and peak */
    int ghist_refl;     /* 1/3 way between medn and peak */
    int ghist_refgf;    /* 1/5 way between medn and peak */  
    int ghist_bias;     /* 35thpercentile--poor mans approx overscan*/
    int *ghistarr;      /* histogram array */
} GHIST;

/* doubles, python will not work properly with floats -ps */
typedef struct g_fiberdata{    
    int     g_nfibers;   /* number of useful fibers in flat frame */
    int    *g_fid;       /* id's of fibers in list below */
    double *g_xcen;      /* x centroid of fibers */
    double *g_ycen;      /* y centroid of fibers */
    double *g_fibrad;    /* radii of mask holes */
    double *g_illrad;    /* radii of illuminated portion of fiber */
    double *g_xs;        /* from psf fit, x offset of star from center */
    double *g_ys;        /* from psf fit, y offset of star from center */
    double *g_mag;       /* from psf fit, mag of star */
    double *g_fwhm;      /* from psf fit, fwhm(fullwidth at half maximum) of star */
    double *g_poserr;    /* from psf fit, estimated position error (2d) */
    double *g_fitbkgrd;  /* from psf fit  estimated background */ 
    double *g_fibercts;  /* flux through 2 arc fiber */
    double *g_skycts;    /* direct sky background measure*/   
    double *g_rmswidth;  /* rms radius of star */
    double  g_readnoise; /* readnoise in ADU */
    int     g_npixmask;  /* number of pixels in current mask */
} FIBERDATA;


/* fiber data statistics structure, to be filled in later. */
typedef struct gfiberstat{
    int ctps[16];		   /* counts per second */
    int cts_illum[16];	   /* counts illuminated */
} FIBERSTAT;

typedef struct gstarfit{
    float gsampl;       /* profile amplitude; template has ampl DMAX at org */
    float gsbkgnd;      /* mean background level */
    float gswparam;     /* width parameter--100/sigmasq */
    float gserror;      /* 1d fitting error to profile */
} GSTARFIT;

/********** FUNCTION PROTOTYPES ************/

//GHIST *ipGhistNew(void);
//RET_CODE ipGhistDel(GHIST *obj);

//FIBERDATA *ipFiberdataNew(void);
//RET_CODE ipFiberdataDel(FIBERDATA *obj);

//GSTARFIT *ipGstarfitNew(void);
//RET_CODE ipGstarfitDel(GSTARFIT *obj);

int gbin(
     REGION *inputReg,       /* the input region */
     REGION *outputReg       /* the output region */
    );

int gmakehist(
    REGION *dataReg,        /* the data region */
    struct ghist_t *gptr    /* pointer to ghist structure */     
    );

int gmakeflat(
    REGION *temReg,       /* template picture (input raw guider flat) */
    REGION *flattenReg,   /* flattener picture */
    struct ghist_t *gptr   /* pointer to hist structure */    
    );
    
int gextendmask(
    REGION *flattenReg,     /* mask (flattener) to be extended */
    MASK *maskReg,         /* new mask, extended and negated */
    int fringe              /* distance to extend (pixels) */
    );

int gfindfibers(
    REGION *flatReg,           /* the flat region */
    MASK *maskPtr,             /* pointer to the mask structure */
    struct g_fiberdata *ptr,   /* pointer to output struct */
    struct ghist_t *gptr,       /* histogram struct pointer */
    struct platedata *plptr
    );

int
gfixdark(
    REGION *darkInReg, /* input average dark frame; it has the mean
                                 of the MASKED frame removed by this routine */
    REGION *darkOutReg, /* normalized and masked dark frame */
    MASK *maskPtr,
    struct g_fiberdata *fptr,  /* fiberdata strc--this routine fills out npix*/
    float hotthresh
    );

int ginitseq(void);                      
                        /* sets flats upon a restart after flat processing
                         * data have been reread from disk
                         */
    
int 
gsubdark(
    REGION *dataReg,       /* data frame from which dark is to be subtr.*/
    REGION *maskedDarkReg, /* masked and normalized dark frame */
    REGION *darkReg,       /* normalized dark frame */
    MASK *maskPtr,         /* extended mask --used to establish background */
    struct g_fiberdata *fptr  /* fiberdata struct */    
    );

int gmflatten(
    REGION *dataReg, /* data frame to be flattened */
    REGION *flatReg  /* masked flat field as produced by gmakeflat()*/
    );

int gmakeseeprof(void);

int seeprof2(
    int r2,    /* r^2 in pixels^2 */
    int wp     /* width profile--100/sigma^2 */
    );

int seeprofr(
    int ri,    /* radial INDEX */
    int wp     /* width profile--100/sigma^2 */
    );
    
double gproffit(
    int *starprof,           /* extracted profile */
    int *starwgt,            /* number of points contributing to profile */
    int npt,                 /* number of points in extracted profile */
    int firstguess,          /* first guess at closest INDEX in wpg for
                              * width parameter; if none, use 0, which
                              * starts iteration with index 5, 2 arcsec   
                              */
    struct gstarfit *ptr     /* pointer to output struct */
    );

#if 0
int gfitparams(
    double dx,              /* x position of min in square: -.5 < dx < .5 */
    double dy,              /* y  "                "         "    dy   "  */
    struct gstarfit *ptr    /* output interpolated structure */ 
    );
#endif

int gprofext(
    short int **data,       /* picture */
    int *radprof,           /* int profile array */
    int *radcnt,            /* int weight array */
    int *r2prof,            /* int prof^2 array--should be OK for 12-b data*/
    int xc,                 /* extraction center x */
    int yc,                 /* extraction center y */
    int fxc,                /* fiber center x */
    int fyc,                /* fiber center y */
    int frad,               /* fiber mask radius */
    float *perr2d           /* pointer to 2d error wrt 1d profile */    
    );

int gfindstar(  
    short int **data,           /* data picture, dedarked and flattened */
    int fk,                     /* fiber number */
    struct g_fiberdata *ptr     /* ptr to fiberdata struct */
    );
    
int gfindstars(  
    REGION *dataReg,           /* data picture, dedarked and flattened */
    struct g_fiberdata *ptr,
    int mode     /* ptr to fiberdata struct */
    );
    
int gtranspose(
     REGION *inputReg,       /* the input region */
     REGION *outputReg       /* the output region */
    );

/* these are not, as of this moment, sent out to tcl */

int hotthresh(int xs,       /* x size (pixels) */
	      int ys,       /* y size (pixels) */
	      S16 **pic,    /* pointer to line pointer array */
	      int thresh   /* threshold */
	      );
void
hotfix( int xs,              /* x size */
        int ys,              /* y size */
        S16 **pic     /* pointer to line pointer array */
        );
    
/************************* END, GUIDSTAR.H *********************************/

   

    
    
   


   

#endif
