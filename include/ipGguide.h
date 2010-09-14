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
**      Imported into Dervish: 2003-01-07
**      Eric Neilsen
** Mods
**      New Apogee camera:   2009-07-29
**      Paul Harding      
**      Dustin Lang
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

// Error code; must match guiderImage.py:FWHM_BAD
#define FWHM_BAD 99.99

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
STATIC struct gstarfit istarfit ;  /* initial answer */            
STATIC struct gstarfit fstarfit;  /* final answer */

/* print buffer */
STATIC char pbuf[256];

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

#define MAXFIBERS 20			/* large enough for all possible cartridges */

/* max number of iterations in refining the center; 
 * we should not need to walk futher than the range limit of 40 pix radius. */

#define MAXGFSITER 40

// values of 100/sigma^2 used in width fitting, from sigma=10 pix to sigma~=0.707 pix
#define NCELL 41
#define WPG_NUMERATOR 100
STATIC short int wpg[NCELL] = 
	{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};



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
	// from PSF fit:
	// total star flux, in DN
	double *flux;
	// total sky flux, in DN ?
	double *sky;
	// full-width at half-max of the star, in (binned guidecam) pixels.
	double *fwhm;
	// error in g_xs and g_ys, in (binned guidecam) pixels.
	double *poserr;
	
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
    float sigma;        // width (sigma) in pixels
    float gserror;      /* 1d fitting error to profile */
} GSTARFIT;

/********** FUNCTION PROTOTYPES ************/

// defined in gutils.c, used in ipGguide.c
void  
grot(			   /* general array rotator; puts zeros where 
                                 data are not generated */
     float theta,              /* angle in degrees */
     short **p,              /* source array */
     short **dp,             /* destination array */
     int xsz, int ysz);              /* sizes */

void  
maskrot(
		float theta,              /* angle in degrees */
		char **p,                 /* source mask array */
		char **dp,                /* destination mask array */
		int xsz, int ysz);              /* sizes */


void rotate_region(const REGION* regin, REGION* regout, float theta);

void rotate_mask(const MASK* maskin, MASK* maskout, float theta);

FIBERDATA* fiberdata_new(int nfibers);
void fiberdata_free(FIBERDATA* f);

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

int
gfixdark(
    REGION *darkInReg, /* input average dark frame; it has the mean
                                 of the MASKED frame removed by this routine */
    REGION *darkOutReg, /* normalized and masked dark frame */
    MASK *maskPtr,
    struct g_fiberdata *fptr,  /* fiberdata strc--this routine fills out npix*/
    float hotthresh
    );

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
    
/************************* END, GUIDSTAR.H *********************************/

   

    
    
   


   

#endif
