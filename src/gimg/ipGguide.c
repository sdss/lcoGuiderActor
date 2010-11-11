/*****************************************************************************
******************************************************************************
** 
** FILE:
**	ipGguide.c
**
** ABSTRACT:
**	This file contains code for Jim Gunn's SOP guider.
**      Modified to be called by python / ctypes
**
** ENTRY POINT		SCOPE	DESCRIPTION
** -------------------------------------------------------------------------
** gmakehist            public  fill an a ghist structure based on a guider im
** gmakeflat            public  make a flattener image from a guider flat
** gextendmask          public  create a guider mask with a skirt
** gfindfibers          public  find the fiber locations in an guider im
** ginitseq             public  set up to star guiding on a new plate  ** do we need this
** gmflatten            public  flat field a guider image
** makeseeprofile       public  create a seeing profile
** seeprof2             public  look up a seeing profile point by area
** seeprofr             public  look up a seeing profile point by radius
** gproffit             public  fit a seeing profile in space
** gprofext             public  Extract radial profile as radial index array
** gfindstar            public  find the star in one fiber
** gfindstars           public  fund the stars in all fibers
**
** ENVIRONMENT:
**	ANSI C.
**
** REQUIRED PRODUCTS:
**
** AUTHORS:
**	Creation date:  2002-03-31
**	Jim Gunn
**      Mods for Python: 2009-07-
**      Phoebe Stierhoff, Paul Harding, Craig Loomis
**      Dustin Lang
******************************************************************************
******************************************************************************
*/

#define DEBUG 0			/* no debugging info */

/* ph,ps: currently, many values are hardwired in, we need to create variables for these
   in either header file, or passed in from python. */
   
/* ph,ps: need error handeling to replace dummy code */

/* ph,ps: extend profile and other arrays for large fibers, need for aquisition */

#include <stdio.h>
#include <assert.h>
#include "ipGguide.h"
#include "shLegacy.h"

/***************************************************************************
* Hardware discription
***************************************************************************
*/

/*--- Guide camera---
 * Alta camera images are: 1024x1024x13um native, 
 * Alta binned by 2-> 512x512x26um pixels -> 0.428" / pix   
 * Raw data from the Alta camera is 16 bits but is bit shifted in the gcameraICC,
 * The data regions passed to these C functions is S16 and has been overscan, dark corrected 
 * DMAX = 32767
 */
 
/*--- Cartridges ---
 * The new BOSS cartridges have 14 GUIDE probes and 2 ACQUISITION probes.
 * Six of the guide probes are in the focal plane, for measuring FWHM etc
 * Eight probes are used for focusing, 4 inside focus and 4 outside focus.
 *
 * The old SDSS cartrdges have 9 GUIDE probes 2 LargeGuide probes.  
 */ 

/*--- Fibers ---
 * The large plastic ACQUISITION fibers bundles are 1.5 mm active diameter = 57.7 (26 micron)pixels.
 * A usable size is ~ 55 pixels or ~23arcsec. Increased FRD due to squished fibers at the outer circumference
 * of the acquisition fiber bundle result a fall-off in transmission at edges. 
 *  
 * The GUIDE fibers are sumitomo IGN05-10, 450 +/- 23 microns active diameter = 17.3 (26 micron)pixels 
 * or 7-7.5 arcsec diameter. 
 *
 * The LargeGuide probes are 28 pixels diameter
 */

/* --- Fiber Layout BOSS Cartridge ---
 * Approximate locations of the fibers as seen on the guider image, +ve Y up, and +ve X right
 *
 *        1  16               
 *               13          1,8,9,16 short focus fibers
 *        8   9                
 *    7            15        
 *        5  10              2,7,5,10,15 in focus
 *  2     
 *        6  12   
 *   3            11         4, 6 12, 14 long focus fibers
 *        4  14              3,11 acquisition fibers
 *  
 *  Closest fibers 9,13 center to center  = 77 pixels,  
 *  Illuminated edge to edge min distance = 60 pixels  
 */ 

/* --- Fiber Layout SDSS Cartridge ---
 *  Approximate locations of the fibers as seen on the guider image, +ve Y up, and +ve X right
 *
 *       11
 *    10
 *   9  8  7
 *
 *   6  5  4
 *
 *   3  2  1
 * 
 */

/* All data regions passed to this package are assumed to be binned
 * Flats are taken UNBINNED, but are binned in Python
 * 
 * All the instrumental signature removal is done in python
 * Bias subtraction is a kluge because the camera wont give an overscan region
 * without crashing.
 * The flats and data images currently (Dec 2009) have a pseudo bias subtracted
 * which is derived from the histogram. 
 */
 
/*--- Algorithms ---
 * This module provides code to analyze guider flats and flatten data frames,  
 * locate the guide probe center, centroid stars in the probes, and determine
 * star and sky properties in each probe.
 */
 
/* The flat is histogramed and a series of reference levels are determined
 * which characterize the flat image 
 * A pseudbias level 
 * The MEDIAN is found and assigned to the scattered light level. 
 * The PEAK is defined near the first percentile, and the PEAK - MEDIAN 
 * difference (PMD) used to set reference levels   
 * The reference level REF: halfway between the background and the peak.
 * Low REFL LEVEL one third of PMD above median
 * Gfindfibers REFGF one fifth of PMS above the median, is needed because 
 * the larger pixels of the Alta camera result in peak height of the fibers
 * in the 16x16 binned flat can be lower than REFL, which was used previously.
 * There is too much reflected/scattered light to use the median. 
 *    
 * A mask is made which is 1 for pixels above the reference level and 0 below. 
 * This mask is used to construct a flattener
 * which is the inverse of the flat frame where the mask is 
 * nonzero and zero outside. 
 * 
 * All this stuff is stored in internal static arrays and variables.
 * 
 * Stars are found in the processed frame by running a parabolic cap 
 * corresponding to about 1.5 arcsecond seeing through the fiber images.
 * Circularly averaged profiles are then extracted from pixels in the 
 * neighborhood of this trial maximum, and a double gaussian is fitted
 * to the profile. The mean square error of this fit to the 2D image
 * is evaluated and the pixel for which it is minimized is found. The
 * 'true' minimum is found by quadratic interpolation among that point and
 * its neighbors.
 * 
 * The small fibers are about 17 pixels in diameter, the old aquisition fibers 
 * are 28 and the new aquisition fibers are 55. We use a table for the 
 * exponential (which means we can trivially use more complex functions).
 * 
 * A possible scheme is to use the first guess centers, extract a radial
 * profile, fit for A, sigma, and B in this tiny array (typically 12-20
 * elements), apply this filter to the image and find a new center, and
 * iterate. It is possible that we should force all the B's to be the same??
 * I think this is a nonissue, because the Bs are small, but we should 
 * check on moony nights.
 *
 * We use a table for the exponential (which means we can trivially use more complex functions).
 *
 */

/******************** DATA DECLARATIONS AND DEFINITIONS ******************/

STATIC int havehist = 0;

/* sum of squares of values for masked dark frame */
#if 0
STATIC float gdnorm = 0.;
#endif


/*/*static double sigma2fwhm(double sig) {
 *	// about 2.355
 *	return sig * 2.*sqrt(2.*log(2.));
 *}
 */
/*PH Use the actual conversion factor of sigma to fwhm
 *   for the double gaussian function used, which has 
 *   sigma2/sigma1=0.5 and I2/I1 = 0.1
 */
static double sigma2fwhm(double sig) {
	return sig * 2.468;
}

/*
 static double fwhm2sigma(double fwhm) {
 return fwhm / (2.*sqrt(2.*log(2.)));
 }
 */

void rotate_region(const REGION* regin, REGION* regout, float theta) {
	grot(theta, regin->rows_s16, regout->rows_s16, regin->ncol, regin->nrow);
}

void rotate_mask(const MASK* maskin, MASK* maskout, float theta) {
	maskrot(theta, (char**)maskin->rows, (char**)maskout->rows, maskin->ncol, maskin->nrow);
}

/*******************************************************************
 *
 *      FiberstatNew
 *
 ******************************************************************/
/* ph,ps: added this, will get expanded */
FIBERSTAT *
ipFiberstatNew(void)
{
   int i;
   FIBERSTAT *new_obj = (FIBERSTAT*) shMalloc(sizeof(FIBERSTAT));
   for (i = 0; i<MAXFIBERS; i++){
     new_obj->ctps[i] = 0;
     new_obj->cts_illum[i] = 0;
   }
   return(new_obj);
}



/******************* GMAKEHIST *******************************************
 * This routine constructs an integer histogram from a frame of size
 * xsize X ysize with data with maximum value DMAX.   
 * It also calculates the entries medn, peak, refh, and refl in the ghist
 * structure; these levels are defined above.
 */
 
int
gmakehist(
    REGION *dataReg,        /* the data region */
    struct ghist_t *gptr    /* pointer to ghist structure */
    ) 
{
    int * ptr2;
    S16 *ptr;
    int nc;
    int i;
    
    int ndata = dataReg->ncol * dataReg->nrow;
   
    int thresh;
    int threshm;
    int threshp;
    int threshb;
    int sum;
    S16 **data; 
	

    assert(dataReg->rows_s16 != NULL);
    data = dataReg->rows_s16;
   
    
    /* clear histogram */
    /* This sets the array pointed to by gptr->ghistarr to zero, */
    /* doing so one byte at a time. */
    memset(gptr->ghistarr,0,DMAX*sizeof(int));    

    /* make histogram */
    ptr2 = gptr->ghistarr;
    
    for(i=0;i<dataReg->nrow;i++){
    	nc = dataReg->ncol;
    	
        ptr = data[i];

		if(*ptr > DMAX-1){ 	 	
          printf ("nc ptr data[i] i %d %d %p %d /n", nc, *ptr, data[i], i);
		  shError("Data > DMAX: shorten guiderFlat exposure time");
		  return (-1);
		}
	 
		while(nc--)  {
		  if (*ptr >= 0) {
		    (*(ptr2 + (*ptr++) ))++;
		  }
		}

        /*  this does this:
            for(j=0;j<n;j++){
            (hist[ data[i][j])++;
            }   
        */
    }

    /* run through histogram from the TOP and find peak and median
    * we look for first percentile to define peak 
    * and down 55% of the data to find median, 
    * and down 65% for the pseudo bias estimate
    * The fibers contain about 10 percent of the area;
    * The scattered light is about 20 percent
    *ph values updated for alta camera*/
    
    threshp =   PEAK_PERCENTILE*ndata;  
    threshm = MEDIAN_PERCENTILE*ndata;
    threshb =   BIAS_PERCENTILE*ndata; 

    sum = 0;
    thresh = threshp;
    for(i=DMAX-1;i>=0;i--){
        sum += ptr2[i];
        if (sum > thresh){
            if(thresh == threshp){
                gptr->ghist_peak = i;
                thresh = threshm; /* look for median next */
		continue;
            }
            if(thresh == threshm){
                gptr->ghist_medn = i;
                thresh = threshb; /* look for pseudo bias level next */
		continue;
            }
            if(thresh == threshb){
                gptr->ghist_bias = i;
                break;
            }

        }
    }
    /* set the reference level halfway between median and peak      */
    gptr->ghist_ref  = (gptr->ghist_peak + gptr->ghist_medn)/2;
    /* and the 'low' reference level a third of the way up          */
    gptr->ghist_refl = (2*gptr->ghist_medn + gptr->ghist_peak)/3;
    /* and the 'gfindfiber' reference level 5th of the way up       */
    gptr->ghist_refgf =(4*gptr->ghist_medn + gptr->ghist_peak)/5;

    printf("\n peak,medn,ref,refl,refgf,bias=%d %d %d %d %d %d\n",
           gptr->ghist_peak,gptr->ghist_medn,gptr->ghist_ref,gptr->ghist_refl,gptr->ghist_refgf,gptr->ghist_bias);

    // Error handling temporarily disabled
    if(gptr->ghist_peak < FNORM/8){		//changed
    	
        shError(
    "GMAKEHIST: Flat exposure is too faint. Counts in fibers must be > 500");
	return(SH_GENERIC_ERROR);
    } 
    havehist = 1;
 
    return(SH_SUCCESS);
}


/******************* GMAKEFLAT() ******************************************* 
 * this routine takes a template image (bias and dark subtracted guider flat) 
 * and produces a multiplicative flattener picture in which all 
 * original pixels <= the threshold are zero 
 * and the pixels >= greater than the threshold are inverted about 'norm',
 * norm is usually set to ghist.ghist_peak, being representative of 
 * the high values in the flat frame. )
 * The routine returns zero if successful, -1 if the flattener overflowed
 * 12 bits--this is not possible if the processing steps here are followed
 * on a valid flat frame.
 *
 *
 * It is your responsibility that the two pictures
 * are the same size and OF size xsize X ysize.

   take image, divide by 8, make s16 -ps
 */

int
gmakeflat(
    REGION *temReg,         /* template region (guider flat) */
    REGION *flattenReg,     /* the flattener region */
    struct ghist_t *gptr   /* pointer to hist structure */
    )
{
    S16 **tem;
    S16 **flatten;
    int i, j, val;
    int thresh ;
    int norm;  
    int medn = gptr->ghist_medn;  /* background */

	
	assert(temReg->rows_s16 != NULL);
    tem = temReg->rows_s16;
    
	assert(flattenReg->rows_s16 != NULL);
    flatten = flattenReg->rows_s16;


    if ( (temReg->nrow != flattenReg->nrow) || \
	 (temReg->ncol != flattenReg->ncol) ) {
      shError("GMAKEFLAT: flat and flattener regions must be the same size");
      return(SH_GENERIC_ERROR);
    }

    if(havehist == 0) {
      shError("GMAKEFLAT: no histogram");
      return(SH_GENERIC_ERROR);
    }
  
    norm = gptr->ghist_peak - medn; /* value to set flat portion of fiber to */
    thresh = gptr->ghist_ref - medn;  /* halfway between med and peak */
    
    for(i=0;i<temReg->nrow;i++){
        for(j=0;j<temReg->ncol;j++){
            if(tem[i][j] > thresh){
                val = ((FNORM) * norm)/(tem[i][j] - medn);		//integer flat is scaled to DMAX/8
                if(val > DMAX){	
                   shError("Flattener overflowed: flat frame is not valid flat");
                   return (-1);
                }
                flatten[i][j] = val;
            }else{
                flatten[i][j] = 0;
            }
        }
    }
    return(SH_SUCCESS);
}

/******************* GEXTENDMASK ****************************************
 * This routine enlarges the nonzero areas of the input picture circularly by
 * 'fringe' pixels. It does this by stepping through the picture and
 * running around the periphery of a circle with its center at the point; 
 * if at least MASKPTS points on the periphery are a mask points, the 
 * center is assigned to the mask. In the process the mask is negated; the
 * 1s in the mask are assigned to the regions NOT in the fibers or their
 * fringes. 
 */

int
gextendmask(
    REGION *flattenReg,     /* mask (flattener) to be extended */
    MASK *maskReg,          /* new mask, extended and negated */
    int fringe              /* distance to extend (pixels) */
    )
{
    S16 **maskin;
    unsigned char **maskout; 
    int i, j, k;
    int *xbase;
    int *ybase;
    int ncirc = 0;
    int midflg;
    int ncross;
    int nmaskedin = 0;
    int nmaskedout = 0;

    assert(flattenReg->rows_s16 != NULL);
    maskin = flattenReg->rows_s16;

    if ( (flattenReg->nrow != maskReg->nrow) || \
	 (flattenReg->ncol != maskReg->ncol) ) {
      shError("GEXTENDMASK: mask not of same size as flattener region");
      return(SH_GENERIC_ERROR);
    }
    maskout = maskReg->rows;
           /* first set to all 1s */
    for(i=0;i<(maskReg->nrow);i++){
        for(j=0;j<(maskReg->ncol);j++){
         
            maskout[j][i] = 1;
	
        }
    }
    /* make arrays for circumference */
    xbase = (int *)malloc((8 * fringe + 8)*sizeof(int));
    ybase = (int *)malloc((8 * fringe + 8)*sizeof(int));
    if(xbase == NULL || ybase == NULL){
        shError("EXTENDMASK:cannot allocate periphery memory");
	return(SH_GENERIC_ERROR);
    }
    
    /*compute the first octant*/
    midflg = 0;
    for(k=0;k<=fringe;k++){
        xbase[k] = k;
        ybase[k] = sqrt((double)(fringe*fringe - k*k)) + 0.5; 
        ncirc ++;
        if(ybase[k] <= xbase[k]){
            if(ybase[k] == xbase[k]){
                midflg = 1;
                break; /* there is a midpoint on quadr, index ncirc-1,
                        * number of pts in quadrant is 2*(ncirc-1) incl origin
                        * but excluding point on x axis 
                        */
            } else {
                midflg = 0;
                ncirc --;
                break; /* overshot; no midpoint, last point ncirc-1,
                        * number of pts in quadrant is 2*ncirc-1 incl 
                        * origin but excl point on x axis 
                        */
            }
        }
    }

    /* finish the first quadrant */
    if(midflg != 0){  /* there is a midpoint */
        for(k=ncirc; k<= 2*(ncirc-1) ; k++){
            xbase[k] = ybase[2*(ncirc-1)-k];
            ybase[k] = xbase[2*(ncirc-1)-k];
        };
        ncirc = 2*(ncirc-1);
    }else{            /* there is no midpoint */
        for(k=ncirc; k<= 2*ncirc -1; k++){
            xbase[k] = ybase[2*ncirc - 1 - k];
            ybase[k] = xbase[2*ncirc - 1 - k];
        };
        ncirc = 2*ncirc -1 ;
    }
    /* finish the first half */ 
    for(k=ncirc+1; k <= 2*ncirc; k++){
        xbase[k] =  xbase[2*ncirc-k];
        ybase[k] = -ybase[2*ncirc-k];
    }
    ncirc = 2*ncirc;
    /* finish the circle */
    for(k=ncirc+1; k < 2*ncirc; k++){
        xbase[k] = -xbase[2*ncirc-k];
        ybase[k] =  ybase[2*ncirc-k];
    }
    ncirc=2*ncirc;

#ifdef FINDBUG
    printf("\nxbase,ybase,ncirc,midflg=%d %d %d %d\n",xbase,ybase,ncirc,midflg);
#endif
    /* ok, have array which describes periphery. Now apply. Note that for
       simplicity we do not consider points closer to the picture boundary
       than `fringe'. If that is not good enough, tough. */
    

    /* set up output mask -- NEGATIVE */
    for(i=fringe;i < maskReg->nrow - fringe;i++){
        for(j=fringe;j < maskReg->ncol - fringe;j++){
	    
	    if(maskin[i][j] !=0) nmaskedin++ ;
            ncross = 0;
            for(k=0;k<ncirc;k++){
            	
                if(maskin[i+ybase[k]][j+xbase[k]] !=0) ncross++ ;
            }
            if(ncross > MASKPTS){
                maskout[i][j] = 0;
                nmaskedout++;
            }
       }
    }
    //#ifdef FINDBUG
	printf("npoints masked in,out = %d %d \n",nmaskedin, nmaskedout);
    //#endif

    free(xbase);
    free(ybase); 
    return(SH_SUCCESS);
}

FIBERDATA* fiberdata_new(int N) {
	FIBERDATA* f = calloc(1, sizeof(FIBERDATA));
	f->g_nfibers = N;
	f->g_fid = malloc(N * sizeof(int));
	f->g_xcen = malloc(N * sizeof(double));
	f->g_ycen = malloc(N * sizeof(double));
	f->g_fibrad = malloc(N * sizeof(double));
	f->g_illrad = malloc(N * sizeof(double));
	f->g_xs   = malloc(N * sizeof(double));
	f->g_ys   = malloc(N * sizeof(double));
	f->flux   = malloc(N * sizeof(double));
	f->fwhm   = malloc(N * sizeof(double));
	f->poserr = malloc(N * sizeof(double));
	f->sky    = malloc(N * sizeof(double));
	return f;
}

void fiberdata_free(FIBERDATA* f) {
	free(f->g_fid);
	free(f->g_xcen);
	free(f->g_ycen);
	free(f->g_fibrad);
	free(f->g_illrad);
	free(f->g_xs);
	free(f->g_ys);
	free(f->flux);
	free(f->fwhm);
	free(f->poserr);
	free(f->sky);
	free(f);
}



/*********************** FITTING ROUTINES *******************************
 *ph
 * These routines are what all of this is about. It finds a star in a fiber,
 * returns its position, amplitude, and 'sigma'
 *
 * The guider pixels are .418 arcseconds. The little fibers are 17 pixels, 
 * 7 arcsec, diameter, and the big ones are 27 pixels, 11 arcsec, dia.
 * In 2 arcsec seeing, which is typical, the canonical double gaussian
 * has a small component with a sigma of .8 arcsec, ~2 pixels, and a
 * large component with 1.6 arcsec, 4 pixels. A 3 sigma cutoff is thus
 * at 12 pixels, larger than the small fiber radius and about the size of the
 * large fibers if the star is centered. 
 *
 * We will use the square of the radius as the variable, and evaluate
 * the filter only on even pixel centers, finding the maximum in the
 * end by quadratic interpolation. To this end, we construct a table
 * of 3600 (MAXR2) entries in the function (for photometrics ~3sig =20 pix)
 * ph: keep same for now but its bigger than it needs to be
 * 
 * exp(-s/2) + 0.1*exp(-s/8) - exp(-18) - 0.1*exp(-4.5)/(same at s=0)
 *
 * with spacings in s, which is r^2/sigma^2, of 0.01, thus evaluating
 * the function over the range r/sigma from 0 to 6, gives (3600 = 6**2/0.01)
 *  outside of which it is zero. 
 * Since r^2 is always an integer, this table look-up is
 * exact if 100*r^2/sigma^2 is an integer for any r^2, or 100/sigma^2
 * is an integer. Since sigma is typically about 2 pix or less for the alta
 * camera, this allows for a granularity in sigma^2 of about 3%, qnd 1.5% in sigma. 
 * This is good enough. One can certainly linearly interpolate if required.
 * (photometrics values were sigma ~3 pix, 10 and 5 % granualarity 
 *  for sigma**2 and sig respectivels) 
 *
 * For a small-gaussian sigma of sig and a radus of r, enter this table
 * with the value 100*r^2/sig^2. If you do not wish to interpolate,
 * with r^2 integral, 100/sig^2 must be an integer. If we use a rounded
 * value, should be OK (see above). A coarse search grid might be
 *                              alta      old photometrics
 *   100/sig^2 sig(pix)         fwhm(")   fwhm(")
 *      1      10.0            10.0       7.0     
 *      2       7.0             7.0	  4.9    
 *      3       5.7             5.7	  4.0    
 *      5       4.5             4.5	  3.1    
 *      8       3.5             3.5	  2.5    
 *     13       2.8             2.8	  1.9    
 *     20       2.2             2.2	  1.6    
 *     30       1.8             1.8	  1.3    
 *     50       1.4             1.4	  1.0    
 *    100       1.0             1.0	  0.7
 *
 *
 * dstn says: I moved all arcsec and mag units into the python: this C
 * now works only in pixels and fluxes in Data Numbers.
 *
 * 
 * Note for alta camera the conversion from sigma in pixels to fwhm in arcsec is
 *   0.428(arcec/pix) * 2.354*(sigma/fwhm) = 1.0 = sigp2FwhmAs (0.69 for the photometrics)
 *
 * We begin by finding the maximum in the frame as defined by a 
 * ~3*(sigmalarge) ~15-pixel core (was 21 for photometrics) 
 * based on 1.8 arcsecond seeing. With that center as a first guess,
 * a radial profile is extracted and a sigma is found for a profile fit,
 * which is used to refine the center. This is repeated and the center
 * interpolated quadratically if it is not on the edge.
 *
 * The half-power half-width of this is 1.234*sig = 10*1.234/sqrt(wp), where
 * wp is the width parameter 100/sig^2 so sigma=10/sqrt(wp). 
 * With 0.428 arcsec pixels, the fwhm is 1*10/sqrt(wp) in arcsec.
 *
 * The equivalent width^2 of a profile displaced by d from its center
 * is sig^2 + d^2/2, so we need to SUBTRACT d^2/2 from an offcenter
 * extracted profile to get the true width^2, or equivalently, multiply
 * by (1-d^2/2sig^2). The corrected width is the measured one multiplied
 * by (1-d^2/4sig^2) if d << sig. The relation between fwhm in arcsec and
 * sig in pixels is fwhm" = 1.0*sig_pix, so sig^2 ~ fwhm"^2, and the
 * correction is 1-d^2/4fwhm"^2)
 */

/****************** GMAKEFIBERSTAT() ************************************
 * This routine just fills out the FIBERSTAT structure. right now it
 * just calculates the counts by summing the pixals in a fiber, after
 * the image has been flattened and darks have been subtracted. currently
 * very rough, hasn't been fully tested or implemented. -ps
 * So this works byt drawing a square around each fiber using info
 * from the fiberdata structure and sums the pixals within each square.
 */
int
gmakefiberstat(
    struct g_fiberdata *fptr,
    REGION *datareg
    )
{
    int i, j, k, row0, rowf, col0, colf;
    int nrow = datareg->nrow;
    int ncol =datareg->ncol;
    int nfibers = fptr->g_nfibers;
    S16 sum[MAXFIBERS];
    double *rad = fptr->g_fibrad;
    double *xcen = fptr->g_xcen;
    double *ycen = fptr->g_ycen;
    int total = 0;
    int totaltest =0;
    for(i=0; i<MAXFIBERS; i++) sum[i]=0;
    S16 **pdata = datareg->rows_s16;
    for(k = 0; k<nfibers; k++){
	row0 = (int)(xcen[k]-rad[k]);
	col0 = (int)(ycen[k]-rad[k]);
	rowf = (int)(xcen[k]+rad[k]);
	colf = (int)(ycen[k]+rad[k]);
	for (i = row0; i<=rowf; i++){
		for (j = col0; j<=colf; j++){
			sum[k] =sum[k] + pdata[j][i]; 
		
		}
	}
    }
    for(i=0; i<MAXFIBERS; i++) {
	
	totaltest += sum[i];
    }
    for(i=0; i<ncol; i++){
	for(j = 0; j<nrow; j++){
		total += pdata[j][i];
	}
    }
    
    return(1);
}







/****************** GMAKESEEPROF() **************************************
 * This routine makes the seeing profile described above, and a `sqrt'
 * table to facilitate the calculation of radii in profile extraction. Since
 * we want to do the profile fitting in radial space and want to minimize
 * the number of points, we collapse the extracted profile to an array of
 * radial points indexed in the mirage fashion (with one extra point at
 * r^2 = 2). This may or may not be a good idea.
 */

int 
gmakeseeprof(void)
{
    int i,j;
    int norm;
    int r2,ir;
    double r;
    double x;

    if(seeprofile == NULL){
        seeprofile = (short int *)malloc((MAXR2*2) * sizeof(short int));
        if(seeprofile == NULL){
            shError("MAKEPROFILE:Cannot allocate profile memory");
            return (-1);
        }
    }
    radtable = seeprofile + MAXR2; 
    
    norm = PMAX/(1.1 - exp(-18.) - 0.1*exp(-4.5));
    
    for(i=0;i<MAXR2;i++){
        x = -.005*(float)i;
        seeprofile[i] = norm*(exp(x) + 0.1*exp(0.25*x));
    }
    seeprofile[0] = PMAX;
    
    /*ph,ps: MAXR is because r must be <= 3600 */
    for(i=0;i<MAXR;i++) rvalmean[i] = rvalwt[i] = 0 ;
    for(i=0;i<MAXR;i++){
        for(j=1;j<MAXR;j++){
            r2 = i*i + j*j;
            r = sqrt((double) r2);		//mirage convention
            if(r2 > 2) ir = r+1.5;
            else ir = r2;
            if(ir < MAXR){
                rvalwt[ir]++ ;
                rvalmean[ir] += PMAX*r + 0.5;
            }
        }
    }

    for(i=1;i<MAXR;i++) rvalmean[i] = rvalmean[i]/rvalwt[i];
    
    rvalmean[0] = 0;
    rvalwt[0] = 1;

    radtable[0] = 0;
    radtable[1] = 1;
    radtable[2] = 2;
    for(i=3;i<MAXR2;i++){ 
        radtable[i] = (int)(sqrt((double)(i)) + 1.5);
    }
    /* note that the entries in radtable are a radial INDEX, not a radius.
     * for index=3 or larger, the radius is approximately index - 1, but
     * the indices of 0 , 1, and 2 are values of r^2. This is the standard
     * mirage convention for radial profiles. The 1-d profile fits are
     * in radial profiles with this index.
     */
    return (0);
    
}

/************************* SEEPROF2() ***************************************
 * seeprof(r^2, wp) just does a table look-up, with range checking. Should
 * probably enlarge table and do a fiber size check at the beginning.
 * This could just be a macro or even just inline in that case.
 * For our chosen double gaussian, the half-power point is at 1.234*sigma
 */

int 
seeprof2(
    int r2,    /* r^2 in pixels^2 */
    int wp     /* width profile--100/sigma^2 */
    )
{
    int arg = r2*wp;
    return (arg >= MAXR2 ? 0 : seeprofile[arg]);
}


/************************* SEEPROFR() ***************************************
 * seeprofr(rindex, wp) just does a table look-up, with range checking. Should
 * probably enlarge table and do a fiber size check at the beginning.
 * This could just be a macro or even just inline in that case.
 */

int 
seeprofr(
    int ri,    /* radial INDEX */
    int wp     /* width profile--100/sigma^2 */
    )
{
    int arg;

    arg = ri<3 ? (ri*wp) : ((ri-1)*(ri-1)*wp);
    return (arg >= MAXR2 ? 0 : seeprofile[arg]);
}


/************************** GPROFFIT() ************************************
 * gproffit does the fitting in profile space to determine the width,
 * amplitude, and background. returns sum square error 
 */

/* #define GPROFDEBUG */


    
double   
gproffit(
    int *starprof,           /* extracted profile */
    int *starwgt,            /* number of points contributing to profile */
    int npt,                 /* number of points in extracted profile */
    int firstguess,          /* first guess at closest INDEX in wpg for
                              * width parameter; if none, use 0, which
                              * will start iteration with index 18 
                              * which has value of 24 which represents
                              * 2.0 arcsec [sigtoFWHM * sqrt(100/24)]
                              * for the alta camera.
                              */
    struct gstarfit *ptr     /* pointer to output struct */
    )                          
{
    int err[NCELL];
    int i;
    int sumf ;
    int sump ; 
    int sumw ;
    int sumf2 ;
    int sumfp ;
    int diff;
    int sumerr ;
    int w,f,p;
    int wpdex = (firstguess != 0 ? firstguess : 18);
    double ampl;
    double bkgnd;
    double disc; 
    int dir = 1;
    int iter;
    double amin, a, b, minerr = HUGE_VAL;
    double wpmin = HUGE_VAL;
    int wp;
    
/* defined in ipGguide.h reproduced here for clarity
NCELL = 31
wpg[NCELL] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,24,28,32,36,40,45,50,55,60,70,80,90,100} 
*/

    for(i=0;i<NCELL;i++) err[i] = 0;

    iter = 0;

    do{
        sumf = sump = sumw = sumf2 = sumfp = sumerr = 0;
		// wp == "width parameter"
        /* for current guess at wp, do linear least squares to solve for
         * amplitude, background, and evalutate rms error */
        wp = wpg[wpdex]; 
        for(i=0;i<npt;i++){
            w = starwgt[i];
            /* recall that i is a radial index which is r^2 for 0,1,2 and
             * approx r+1 for 3 or greater , so */
            f = seeprofr(i,wp);
            p = starprof[i];
            sumf  += w*f;
            sumf2 += w*f*f;
            sump  += w*p;
            sumfp += w*f*p;
            sumw  += w;

	    //#ifdef GPROFDEBUG            
            //printf("GP:f,p,w,s:f,f2,p,fp,w:%4d %2d %2d %d %d %d %d %d\n",
            //    f,p,w,sumf,sumf2,sump,sumfp,sumw);
	    //#endif
                
        }
        disc = (double)sumw*(double)sumf2 - (double)sumf*(double)sumf ;
        bkgnd = ((double)sumf2*(double)sump - (double)sumf*(double)sumfp)/disc ;
        ampl  = ((double)sumw*(double)sumfp - (double)sump*(double)sumf)/disc ;
        for(i=0;i<npt;i++){
            f = seeprofr(i,wp);
            diff = starprof[i] - ampl*f - bkgnd ;
            sumerr += starwgt[i]*diff*diff;
        }
        err[wpdex] = sumerr;        

#ifdef GPROFDEBUG        
        printf("GPROFFIT:iter,ampl,bkgnd,wp,err=%d %f %f %d %d\n",
            iter,ampl,bkgnd,wp,sumerr);
#endif

		if (iter == 0) {
			// First time through the loop: we've computed an error for our
			// initial guess.  Now we move +1 (dir = 1 above) and then
			// check whether that helped or hurt.
			wpdex += dir;
		} else if (iter == 1) {
			// We've got our initial guess from last time through the loop,
			// and initial guess + 1.  If that helped, keep going in dir=1;
			// otherwise, flip dir and backtrack.

			// We ASSUME that the errors are monotonic, as they should
			// be for reasonable first guesses.
            if (err[wpdex] > err[wpdex - dir]) {
				// got worse: backtrack and reverse direction.
                wpdex -= 2*dir;
                dir = -dir;
            } else {
				// right way, keep going
                wpdex += dir;
            }
        } else {
			// Keep going until we pass the minimum.
            if( err[wpdex] > err[wpdex - dir] ){ /* passed minimum; find it */
                wpdex -= dir;   /* minimum point */
				// interpolate.
				b = 0.5*(double)(err[wpdex + 1] - err[wpdex - 1]);
                a = -0.5*(double)(2*err[wpdex] -err[wpdex-1] - err[wpdex+1]);
                amin = -b/(2.*a);
                minerr = err[wpdex] + amin*(amin*a + b);
                /* amin is the fractional location of the minimum taking
                 * in `wpdex' space; we map this into wpg space quadratically,
                 * and extract the value of wpg at the minimum:
                 */
                wpmin = ((wpg[wpdex+1] + wpg[wpdex-1] - 2*wpg[wpdex])*amin
                        + (wpg[wpdex+1] - wpg[wpdex-1]))*amin/2. 
                        + (double) wpg[wpdex];
#ifdef GPROFDEBUG
                printf("wpdex,amin,minerr,wpmin = %d %f %f %f\n",
                    wpdex,amin,minerr, wpmin);
#endif

                break;                        
            } else {       /* just keep truckin' */
                wpdex += dir ;
            }
        }
        iter++;
		// don't walk off the edge!
    } while(wpdex > 0 && wpdex < NCELL);

    if(wpdex == 0) wpmin = wpg[0];
    if(wpdex == NCELL) wpmin = wpg[NCELL -1];  /* at end */

	assert(wpmin != HUGE_VAL);

    /* populate output struct */
    ptr->gsampl  = ampl*PMAX;
    ptr->gsbkgnd = bkgnd;
    ptr->sigma = sqrt((double)WPG_NUMERATOR / wpmin);
    ptr->gserror = minerr;
    return (minerr);
}


/*********************** GPROFEXT() ***************************************
 * Extract radial profile as radial index array; returns #of points
 */

/*#define EXTDEBUG*/

int
gprofext(
    short int **data,       /* picture */
    int *radprof,           /* int profile array--need to prevent overflow */
    int *radcnt,            /* int weight array */
    int *r2prof,            /* int prof^2 array--should be OK for 12-b data*/  
    int xc,                 /* extraction center x */
    int yc,                 /* extraction center y */
    int fxc,                /* fiber center x */
    int fyc,                /* fiber center y */
    int frad,               /* fiber mask radius */
    float *perr2d           /* pointer to 2d error wrt 1d profile */
    )
{
    int i,j,k;
    int ii,jj;
    int ir;
    int d;
    int radmax;

	
#ifdef EXTDEBUG
    printf("GPROFEXT:data,radprof,radcnt=%d %d %d\n",**data,*radprof,*radcnt);
    printf("GPROFEXT:x,y,fx,fy,frad = %d %d %d %d %d\n",xc,yc,fxc,fyc,frad);
#endif

    /* clear arrays */
    for(k=0;k<2*frad;k++){
        radprof[k] = 0;
        radcnt[k] = 0;
        r2prof[k] = 0;
    }

    radmax = 0;
    for(i= -frad; i <= frad; i++){
        ii = fyc + i;
        for(j= -frad; j <= frad; j++){
            jj = fxc + j;
            d = data[ii][jj];
            
            if(d != 0){
				int const k = (ii-yc)*(ii-yc) + (jj-xc)*(jj-xc);
				if (k >= MAXR2 - 1) {
					continue;
				}
		
                ir = radtable[k];

                radprof[ir] += d;
                r2prof[ir] += d*d;
                radcnt[ir]++;
                if(ir > radmax) radmax = ir;
            }
        }
    }
    /* normalize profile */
    *perr2d = 0.;
    for(k=0;k<=radmax;k++){
        if(radcnt[k] != 0){
        	
			// compute sample variance... ish
            r2prof[k] = r2prof[k] - 
				(double)radprof[k]*(double)radprof[k]/(double)radcnt[k];
            
            *perr2d += (float)r2prof[k];
            radprof[k] /= radcnt[k];
        }
    }
    return (radmax);
}

/*********************** GFINDSTAR() **************************************
 * This routine does the work. It first finds the star by running a 21-point
 * cap through the data, and then refines the position by extracting profiles,
 * fitting them to seeing profiles, and evaluating errors in a 9-point grid
 * around the guess. Finally, a quadratic interpolation yields the best
 * fit point and the profile parameters. This routine does one fiber; it
 * is looped over in the next routine, gfindstars(), to do all the fibers.
 */

/*#define FSTDB1*/

int 
gfindstar(  
    short int **data,           /* data picture, dedarked and flattened */
    int fk,                     /* fiber index in g_fiberdata array; NOT a fiber name */
    struct g_fiberdata *ptr     /* ptr to fiberdata struct */
    )    
{
    int i,j,k,ii,jj;
    int npt;
    int niter = 0;
    int radprof[90];
    int radcnt[90];
    int r2prof[90];
    float err2d;
    int err = 0;
    short int smcap[6];
    const int xc = (int) (ptr->g_xcen[fk] + 0.5) ;    /* fiber xcen */
    const int yc = (int) (ptr->g_ycen[fk] + 0.5) ;    /* fiber ycen */
    const int rc = (int) (ptr->g_fibrad[fk] + 0.5);   /* mask radius */
    const int ri = (int) (ptr->g_illrad[fk] + 0.5);   /* illum fiber radius */
    
    int sum, sumcap;
    int max = 0;
    int maxi;
    int maxj;
    int d;
    int maxwalk2;
    double errmin;
    double ax, bx, ay, by;
    double xs,ys;
    double dx,dy;
    float xoff;
    float yoff;
    //float fwhm0;   /* fwhm corrected for calc on integer pixels */
    float pixnoise;
    float flux;
    float sigma;
    float sigmax;
    float fudge;
    float rsq;
    float roff;

	float gs2derr[3][3];
    
    if(seeprofile == NULL){
        err = gmakeseeprof();
        if( err != 0){
            shError("GFINDSTAR: Error in gmakeseeprof()");
            return (-1);
        }
    }
    
    /* first populate smcap, the filter to be used to find maxima; this
     * corresponds to a parabolic representation of a gaussian with sigma
     * about 1.8 alta pixels, or 1.8 arcsecond seeing
     */
    //smcap[k] = 4096 - 270*k;    Original, triangular but not parabolic? k**2 *ph*
    //for now scale 4096/270   PMAX    

    for(k=0;k<6;k++) smcap[k] = PMAX - (PMAX/15.1)*k;    
    sumcap = smcap[0] + 4*smcap[1] + 4*smcap[2] + 4*smcap[4] + 8*smcap[5];

	//printf("xc,yc %g,%g\n", ptr->g_xcen[fk], ptr->g_ycen[fk]);
	//printf("rc %g, ri %g\n", ptr->g_fibrad[fk], ptr->g_illrad[fk]);
    
    /* first smooth the picture in the neighborhood of the fiber with
     * this cap and find the maximum. we do not deal properly with edge
     * effects, so if the star is closer than 2 pixels to the edge of the
     * fiber, we will not find its first-guess center properly. Tough. 
     */
    /* init maxi, maxj to fiber center xc,yc so they will always be valid */
    maxi = yc;
    maxj = xc;
    /* rc (mask radius) but only data within ri (illum radius) are non zero */
    for(i= -rc+2;i<=rc-2;i++){
        ii = yc + i;
        for(j= -rc+2;j <= rc-2;j++){
            jj = xc + j;
            d = data[ii][jj];
            if(d != 0){
              sum =  smcap[0] * (data[ii][jj] )
                 + smcap[1] * (data[ii-1][jj]   + data[ii+1][jj] 
                             + data[ii][jj-1]   + data[ii][jj+1])
                 + smcap[2] * (data[ii+1][jj+1] + data[ii-1][jj+1] 
                             + data[ii+1][jj-1] + data[ii-1][jj-1])
                 + smcap[4] * (data[ii-2][jj]   + data[ii+2][jj]
                             + data[ii][jj-2]   + data[ii][jj+2])
                 + smcap[5] * (data[ii-2][jj+1] + data[ii-1][jj+2]
                             + data[ii+1][jj+2] + data[ii+2][jj+1]
                             + data[ii+2][jj-1] + data[ii+1][jj-2]
                             + data[ii-1][jj-2] + data[ii-2][jj-1]);
              if(sum > max){
                max = sum;            
                maxi = ii;
                maxj = jj;
              }
            }
        }
    }
   
    sprintf(pbuf,"Gguide Fiber %d: xc,yc,rc = %d %d %d",
	    fk,xc,yc,rc);
    shDebug(DEBUG,pbuf);
    sprintf(pbuf,"Gguide Fiber %d: max,maxj,maxi= %d %d %d ",
	    fk,max/sumcap,maxj,maxi);
    shDebug(DEBUG,pbuf);        
   
    /* OK, use this as first guess at maximum. Extract radial profiles in
     * a 3x3 gridlet about this, and walk to find minimum fitting error */ 

   
    do{
        for(i=-1;i<=1;i++){
            ii = maxi + i;
            for(j=-1;j<=1;j++){
                jj = maxj + j;
				// Make sure we look at all pixels within the fiber, and hardly any more.
				// CPU time is expensive, dammit.
				// This could just be 2 * fiber radius.  Just sayin'.
                npt = sqrt((double)((jj-xc)*(jj-xc)+(ii-yc)*(ii-yc))) + rc + 2;
                 
				if (rc >= sizeof(radprof)/2) {
					sprintf(pbuf, "GFINDSTAR: Guider Fiber %d has radius %d, limit is %d.\n",
							fk, rc, (int)sizeof(radprof)/2);
					shError(pbuf);
					return -1;
				}
				memset(radprof, 0, sizeof(radprof));
				memset(radcnt, 0, sizeof(radcnt));
				memset(r2prof, 0, sizeof(r2prof));

                (void)gprofext(data, radprof, radcnt, r2prof, 
							   maxj + j, maxi + i, xc, yc, rc, &err2d);
                gs2derr[i+1][j+1] = err2d;
                
                (void)gproffit(radprof,radcnt,npt,0,&istarfit); 
#ifdef FSTDB1
                printf(
        "\njj=%d ii=%d err=%6.0f ampl=%5.1f bkgnd=%4.1f fwhm_pix=%4.2f, sig_pix=%4.2f\n\n",
                    jj,ii,err2d,istarfit.gsampl,istarfit.gsbkgnd,
		sigma2fwhm(istarfit.sigma), istarfit.sigma);
#endif
            }
        }
        ptr->sky[fk] = istarfit.gsbkgnd;
        
        /* have error matrix. Find minimum */
        errmin = gs2derr[1][1];
       
        ii = jj = 0;
        
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                if(gs2derr[i][j] < errmin){
                    errmin = gs2derr[i][j];
                    ii = i-1;
                    jj = j-1;
                }
            }
        }

		#ifdef FSTDB1
		printf("jj,ii,errmin:%d %d %5.0f\n",jj,ii,errmin);
		#endif


	/* do not walk too far.  we only have the profile arrays
	 * out to r^2=MAXR2 (3600) pixels, so we have to limit the 
	 * distance of the extraction center from the fiber center. 
	 * the arithmetic works out that the extraction center must at a 
	 * distance less than 40 (alta) pixels from the center of the fiber.
	 * that is 1.4 acquisition fiber radii, 4.4 guide fiber radii
         * 2.8 large guide radii. 

	 * As of now, the code walks at most one pixel in 
	 * each coord so the decision is binary. The walk was limited to match
         * the 40 pixel limit above. We should modify this for the BOSS cartridges.
	 */

        /* PH: The large fibers should let us center up quickly, 
         * so we don't need to look in the noise for stellar winds in the little fibers.
         * Thus set the max radius to be fiber dependent.         
         * Try 1.4 * illum radius, may have to make
         * if radius >= 1.4 illum rad, do not walk and bail out of this loop 
         * It turns out that letting the code walk off into never-never land is a good
         * way of flagging bad fits, so leave it as originally with 1600 */
        //maxwalk2 = 1.4*1.4*ri*ri;

	maxwalk2 = 40*40; 
	//printf("dist from starting point: %g\n", sqrt((maxi+ii-yc)*(maxi+ii-yc)+(maxj+jj-xc)*(maxj+jj-xc)));
	//printf("iter %i / %i\n", niter, MAXGFSITER);
	if((ii != 0 || jj != 0) &&
	   (niter<MAXGFSITER) &&
	   (((maxi+ii-yc)*(maxi+ii-yc)+(maxj+jj-xc)*(maxj+jj-xc))<maxwalk2)) {
		/* min at edge--walk */
            shDebug(DEBUG,"*");
       
            /* now we could save a lot of work if we just moved the good
             * stuff, marked it, and only did the new ones on the new boundary;
             * just a little bookkeeping, but skip it for now--just update
             * and go on
             */
			/* but do not walk too far.  we only have the profile arrays
			 * out to r^2=MAXR2 (3600) pixels.  if we limit the center for the 
			 * extraction to +/- 2 fiber radii, we are safe, as fibrad<30
			 * as of now, the code walks at most one pixel in 
			 * each coord, so if the resulting radius is >= 2 fibrad, go
			 * back and then bail out of this loop.  iop will do quality
			 * control and reject bad fits.   
			 */
			maxi += ii;
			maxj += jj;
			niter++;
			//printf("Walked to %i,%i\n", maxj, maxi);
        } else {   /* have minimum or at edge of range. Get out and go home */
			//printf("Quitting, at %i,%i\n", maxj, maxi);
			//printf("Middle was: %g\n", gs2derr[1][1]);
		if (rc >= sizeof(radprof)/2) {
			sprintf(pbuf, "GFINDSTAR: Guider Fiber %d has radius %d, limit is %d.\n",
					fk, rc, (int)sizeof(radprof)/2);
			shError(pbuf);
			return -1;
		}
		memset(radprof, 0, sizeof(radprof));
		memset(radcnt, 0, sizeof(radcnt));
		memset(r2prof, 0, sizeof(r2prof));

            (void)gprofext(data, radprof, radcnt, r2prof, 
                    maxj, maxi, xc, yc, rc, &err2d);
            gs2derr[1][1] = err2d;
            (void)gproffit(radprof,radcnt,npt,0,&fstarfit);

			/*
			 printf("gprofext: j=%i, i=%i, err=%g\n", maxj, maxi, err2d);
			 printf("jj=%d ii=%d err=%6.0f ampl=%5.1f bkgnd=%4.1f fwhm=%4.2f, sig=%4.2f\n",
			 maxj,maxi,err2d,fstarfit.gsampl,fstarfit.gsbkgnd,
			 fstarfit.gswparam, fstarfit.gswparam/0.69);
			 printf("Middle is: %g\n", gs2derr[1][1]);
			 */
        
            bx = 0.5*(gs2derr[1][2] - gs2derr[1][0]);
            ax = 0.5*(gs2derr[1][2] - 2.*gs2derr[1][1]
                        + gs2derr[1][0]);
            by = 0.5*(gs2derr[2][1] - gs2derr[0][1]);
            ay = 0.5*(gs2derr[2][1] - 2.*gs2derr[1][1] 
                        + gs2derr[0][1]);
            dx = -0.5*bx/ax;
            dy = -0.5*by/ay;

			//printf("dx, dy %g,%g\n", dx,dy);
           
            xs = maxj + dx;
            ys = maxi + dy; 
            xoff = xs - ptr->g_xcen[fk];  
            yoff = ys - ptr->g_ycen[fk];     
            /*modified to return star centers xs and ys rather than offsets from fiber center */
            ptr->g_xs[fk] = xs;
            ptr->g_ys[fk] = ys;
            roff = sqrt(xoff*xoff + yoff*yoff);

			sigma = fstarfit.sigma;
			flux = fstarfit.gsampl * 2.0 * M_PI * sigma*sigma;
			ptr->fwhm[fk] = sigma2fwhm(sigma);
			ptr->flux[fk] = flux;
			//fwhm0 = fwhm*(1. - (dx*dx + dy*dy)/(4*fwhm*fwhm));

            if(roff > (float)ri) {           /* peak off fiber edge */
				ptr->flux[fk] = 0;
                ptr->fwhm[fk] = FWHM_BAD;
			}

            /* error calculation--this was added later, so we have not
             * been very clever about the input quantities--we have to
             * take what we can get. We approximate the star by a gaussian
             * with the derived sigma.
             */

            /* average pixnoise for the centroid calc is the RSS of the 
             * readvariance and 2/9 the peak poisson variance--but
             * evidently wrong--simulations suggest about 1.0 */

            pixnoise = ( ptr->g_readnoise * ptr->g_readnoise + 
						 1.0 * fstarfit.gsampl/CCDGAIN );
            pixnoise = pixnoise > 0. ? sqrt(pixnoise) : 99.;

            sigmax = pixnoise * sigma * 2.*sqrt(2.0*M_PI) / flux;

            /* the error estimate is quintupled at the center and mpy'd
             * by 25 at the edge, and varies like r^6, which may be too slow.
             * chisq is more-or-less OK, but is a bit too big for bright
             * stars and too small for faint. There is some stupid error in
             * my analysis, but the form used here is essentially empirical
             * and anyway, works OK.
             */
            rsq = (xoff*xoff + yoff*yoff)/(ri*ri); 
            fudge = 5. * ( 1. + 5.*(rsq*rsq*rsq));
            ptr->poserr[fk] = sigmax * fudge;
            if (pixnoise == 99.)
				ptr->poserr[fk] = 99.99;

            /* advertise our results */
            sprintf(pbuf,"Gguide xs,ys=%5.2f %5.2f",xs,ys);
            shDebug(DEBUG,pbuf);

            sprintf(pbuf,
					"Gguide err=%6.0f ampl=%5.1f bkgnd=%4.1f",
					err2d,fstarfit.gsampl,fstarfit.gsbkgnd); 
            shDebug(DEBUG,pbuf);

            sprintf(pbuf,
					"Gguide fwhm=%4.2f sig=%4.2f",
					ptr->fwhm[fk], sigma);
            shDebug(DEBUG,pbuf);

            sprintf(pbuf,
					"Gguide flux=%4.2f poserr=%4.3f",
					ptr->flux[fk], ptr->poserr[fk]);
            shDebug(DEBUG,pbuf);

            break;
        }
    }while(1);

    return (0);
}                
               
/********************* GFINDSTARS() ****************************************
 * Does gfindstar() on each fiber in turn.
 */

int 
gfindstars(
    REGION *dataReg,            /* data picture, dedarked and flattened */
    struct g_fiberdata *ptr,	/* ptr to fiberdata struct */
    int mode				    /*data frame mode: 1 for data frames, 0 for spot frames */
    )    
{
    S16 **data;           /* data picture, dedarked and flattened */
    int nfiber = ptr->g_nfibers;
    int i;
    int err;
    
    
    assert(dataReg->rows_s16 != NULL);
    data = dataReg->rows_s16;
    
    for(i=0;i<nfiber;i++){
    	if (mode==0 && i==16){
    		err = gfindstar(data,16,ptr);
    		 if(err !=0){
    	        sprintf(pbuf,"GFINDSTARS:Fiber error, fiber %d",i);
    	        shError(pbuf);
    	        return (SH_GENERIC_ERROR);
    	  }
    	}
    	else if (mode==1){
    	    err = gfindstar(data,i,ptr);
    	 //   err = gfindstar(data,16,ptr);
    	  if(err !=0){
    	        sprintf(pbuf,"GFINDSTARS:Fiber error, fiber %d",i);
    	        shError(pbuf);
    	        return (SH_GENERIC_ERROR);
    	  }
      	   
        }
    }
    /* if mode ==0, we want to use gfindstar on i = SPOTID. else, run on nifber */
    /*  if(mode == 0){
     	   err = gfindstar(data,16,ptr);
    	   err = gfindstar(data,16,ptr);
    	   err = gfindstar(data,16,ptr);
    	   if(err !=0){
              sprintf(pbuf,"GFINDSTARS:Fiber error, fiber %d",i);
              shError(pbuf);
              return (SH_GENERIC_ERROR);
           }
	}
        else if (mode == 1){
  		
        }
     */
    return (SH_SUCCESS);
}

/********************** END MODULE, GUIDER.C *******************************/

