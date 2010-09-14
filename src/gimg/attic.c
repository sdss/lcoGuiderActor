/** Old, unused code that might be useful for reference purposes. */

/* conversion from sigma in pixels to fwhm in arcsec.
 Note here the convenient coincidence that the guide camera pixels,
 when binned by 2, are 0.428"/pixel, and that number is remarkably close
 to 1/2.35 (sigma-to-fwhm conversion factor).  This number should be
 something slightly different than 1.0
*/
# define sigp2FwhmAs 1.0

STATIC S16 badrowval[MAXBADPIX];
STATIC S16 badcolval[MAXBADPIX];
STATIC int nbadpix = 0;

/*jeg set bin size on photometrics to break image into 16 bins
ph maintain the bin size rather than the number of bins for Alta camera.
*/ 

#define BINSIZE	16		/*BINSIZE remains constant*/

#define BINTHRESH 200	/* can probably be increased to 500, assuming flat will have values of ~10000 at least*/



#define NBINS CCD_SIZE/BINSIZE		/* should be an integer */

#define MAXPEAKS 50

#define MAXFIB (17 + 1)                 /* Fibre IDs are sometimes 1-indexed */

#define SPOTID (MAXFIB-2)		/*fiber alignment spot lives in, zero indexed */

/* arrays for holding x,y vals of bad pixels we want to mask.  
 * limited to 2000 entries as we have bigger problems if more than 2%
 * of the pixels are above our threshold for badness 
 */
#define MAXBADPIX 2000

//GHIST *ipGhistNew(void);
//RET_CODE ipGhistDel(GHIST *obj);

//FIBERDATA *ipFiberdataNew(void);
//RET_CODE ipFiberdataDel(FIBERDATA *obj);

//GSTARFIT *ipGstarfitNew(void);
//RET_CODE ipGstarfitDel(GSTARFIT *obj);


int ginitseq(void);                      
                        /* sets flats upon a restart after flat processing
                         * data have been reread from disk
                         */
    
#if 0
int gfitparams(
    double dx,              /* x position of min in square: -.5 < dx < .5 */
    double dy,              /* y  "                "         "    dy   "  */
    struct gstarfit *ptr    /* output interpolated structure */ 
    );
#endif

int gtranspose(
     REGION *inputReg,       /* the input region */
     REGION *outputReg       /* the output region */
    );

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
    
int gfindfibers(
    REGION *flatReg,           /* the flat region */
    MASK *maskPtr,             /* pointer to the mask structure */
    struct g_fiberdata *ptr,   /* pointer to output struct */
    struct ghist_t *gptr,       /* histogram struct pointer */
    struct platedata *plptr
    );

   


/*********************** GMFLATTEN()  ***************************************/
/* This routine applies the masked flat to a data frame; the data frame
 * MUST have been dedarked and dezeroed by gsubdark, because it is ASSUMED
 * in this routine that the bias has been removed 
 */

int
gmflatten(
    REGION *dataReg, /* data frame to be flattened */
    REGION *flatReg  /* masked flat field as produced by gmakeflat()*/
    )
{
    S16 **data;   /* data frame to be flattened */
    S16 **flat;   /* masked flat field as produced by gmakeflat()*/
    int xsize = dataReg->ncol;     /* number of columns */ 
    int ysize = dataReg->nrow;     /* number of rows */
    int i, j;

    if(haveflatdata == 0) {
      shError("GMFLATTEN: no flat data");
      return(SH_GENERIC_ERROR);
    }
    
    assert(dataReg->rows_s16 != NULL);
    data = dataReg->rows_s16;
   
	assert(flatReg->rows_s16 != NULL);
    flat = flatReg->rows_s16;
    

    for(i=0;i<ysize;i++){
        for(j=0;j<xsize;j++){
            data[i][j] = (data[i][j] * flat[i][j]) >> 12;	//ph,ps:lets try to do this with real numbers
        }
    }
    return (SH_SUCCESS);
} 


/********************* GTRANSPOSE() ****************************************
 * Transposes the alta image to be oriented the same as Photometrics image
 * Write to new reg so can cope with possibility of non-square images 
 * ph april 2009
 */

int 
gtranspose(
	REGION *inReg,
    REGION *outReg
	)
{
	int nrow, ncol;
	int i, j;
	S16 **datain;
	S16 **dataout; 
    
#if 0
	printf("gtranspose inreg->rows_s16: 0x%p\n",inReg->rows_s16);
#endif

	assert(inReg->rows_s16 != NULL);
	datain = inReg->rows_s16;
    
	assert(outReg->rows_s16 != NULL);
	dataout = outReg->rows_s16;
    

	nrow = inReg->nrow;
	ncol = inReg->ncol;
	printf ("%d %d %d %d \n", nrow, ncol, outReg->nrow, outReg->ncol);
	printf ("%d\n", datain[0][0]);

	for(i=0; i< nrow-1; i++) {
		for(j=0; j< ncol-1; j++) {
    		dataout[j][i] = datain[i][j];
		}   
	}
	return(SH_SUCCESS);
}

#define FINDEBUG

/******************* GFINDFIBERS() ***************************************
 * This routine deals with BOSS and Marvels cartridges.
 * 512 x 512 Alta frame binned to 24 x 24 always results
 * in an image with a single maximum for each fiber or at worst two
 * adjacent maximum pixels with the same value.
 *
 * The routine therefore works by 
 * 1. binning the input frame to 32x32 
 * 2. finding maxima above the ghist reference level refgf
 * 3. rejecting adjacent maxima
 * 4. quadratically interpolate to find maxima in original picture, which
 *    we identify with fiber centers
 * 5. find median x,y offsets to canonical fiber list by Finkbeiner's method
 * 6. identify fibers by matching peak list with canonical list
 * 7. check that fibers are in zeros of skirt mask
 * 8. find radii of holes in skirt mask
 * 9. find fiber centroids
 * 10. stuff output structure
 *
 * The routine is UGLY.
 */


int
gfindfibers(
    REGION *flatReg,           /* the flat region */
    MASK *maskPtr,             /* pointer to the mask structure */
    struct g_fiberdata *ptr,   /* pointer to output struct */
    struct ghist_t *gptr,      /* histogram struct pointer */
    PLATEDATA *plptr
    )
{
    S16 **flat;                          /* flat frame */
    U8 **mask = maskPtr->rows;  /* negative 'skirt' mask from gextendmask */
    short int xsize = flatReg->ncol;     /* number of columns */ 
    short int ysize = flatReg->nrow;     /* number of rows */
    
    
    
    int thresh = gptr->ghist_refl;
    int medn = gptr->ghist_medn;
    int peak = gptr->ghist_peak;
    int xsbin = xsize/BINSIZE;  /* size of binned picture in x--normally 32 */
    int i, j, k, l;
    int ii, jj;
    int *binpic[NBINS];
    int *binarray;
    int xmax;
    int ymax;
    int npk;
    int xpk[MAXPEAKS];
    int ypk[MAXPEAKS];
    int apk[MAXPEAKS];
    int fid[MAXPEAKS];
    int firad[MAXPEAKS]; /* mask hole radius */
    int dxarr[14];
    int dyarr[14];
    int nid;
    int nposs;
    int x,y;
    int val;
    int xpkoff;
    int ypkoff;
    int xmaxpk;
    int ymaxpk;
    int dx, dy;
    int sumi, sumj, sump;
    double sum, sumx, sumy, wgt;
    double fiberx[MAXFIB + 1];		/* fiber numbers are 1-indexed */
    double fibery[MAXFIB + 1];
    int fiberid[MAXFIB + 1];
	int nfib = 0;

    /* convert input fiber arrays from 1 to zero indexed, only use fibers that exist */
    for(nfib = i = 0;i<MAXFIB;i++){
        if (plptr->exists[i]){
	
	    fiberx[nfib]=(plptr->xcen)[i];
	    fibery[nfib]=(plptr->ycen)[i];
	    fiberid[nfib]=(plptr->gProbeId)[i];
	    nfib++;
			
		}
	}
	
	printf("thresh, medn, peak: %d %d %d\n", thresh, medn, peak);
		
	//printf("%d %g %g %g %g %g %g %g %g %g\n", *ptr->g_fid, *ptr->g_xcen, *ptr->g_ycen,
	//*ptr->g_fibrad, *ptr->g_illrad, *ptr->g_xs, *ptr->g_ys,
	//*ptr->g_mag, *ptr->g_fwhm, *ptr->g_poserr);
	
    if(havemask != 2) {
      shError("GFINDFIBERS: no extended mask");
      return(SH_GENERIC_ERROR);
    }

    haveflatdata = 0;  /* incomplete until we finish this routine */
    
	assert(flatReg->rows_s16 != NULL);
    flat = flatReg->rows_s16;
    
  

    if(gptr->ghistarr == NULL){
        shError("GFINDFIBERS:No histogram. I cannot proceed");
	return(SH_GENERIC_ERROR);
    }   

    if ( (flatReg->nrow != maskPtr->nrow) || \
	 (flatReg->ncol != maskPtr->ncol) ) {
      shError("GFINDFIBERS: mask not of same size as flattener region");
      return(SH_GENERIC_ERROR);
    }
    
    
    /* allocate 32x32 picture memory (pointer to bin array) */    
    binarray = (int *)malloc(32*xsbin*sizeof(int));
    if(binarray == NULL){
        shError("GFINDFIBERS:Cannot allocate binpic memory");
	return(SH_GENERIC_ERROR);
    }
   printf ("%d %d %d %d \n", BINSIZE, xsize, ysize, xsbin);
    for(i=0;i<32;i++) binpic[i] = binarray + xsbin * i ;
    memset(binarray,0,32*xsbin*sizeof(int));     

#ifdef FINDEBUG
    printf("\nbinarray,binpic,thresh,BINSIZE,xsbin =0x%p 0x%p %d %d %d\n",
        binarray,binpic, thresh,BINSIZE,xsbin);
#endif
 

    /* bin the picture */
    ii = 0;
    ymax = BINSIZE;			/* BINSIZE */
    for(i=0;i<ysize;i++){	/* ysize=512 */
        if(i == ymax){
            ymax += BINSIZE;
	        ii++;
        }
        jj = 0;
        xmax = BINSIZE;
        for(j=0;j<xsize;j++){	/* xsize=512 */
            if(j == xmax){
                xmax += BINSIZE;
                jj++;
            }
            binpic[ii][jj] += flat[i][j];
         
        }
    }
       
    /* normalize the picture */
    for(i=0;i<32;i++){
        for(j=0;j<xsbin;j++){
            binpic[i][j] /= BINSIZE*BINSIZE;
        }
    }

    /* find peaks in binned image*/
    npk = 0;
    for(i=1;i<31;i++){				        /* maxima cannot be at edge pixals, so 1-31 */
        for(j=1;j<xsbin-1;j++){
            if((val = binpic[i][j]) >  BINTHRESH 	/* check surrounding pixals */
            && val >= binpic[i][j+1]
            && val >= binpic[i][j-1]
            && val >= binpic[i-1][j-1]
            && val >= binpic[i-1][j] 
            && val >= binpic[i-1][j+1]
            && val >= binpic[i+1][j-1]
            && val >= binpic[i+1][j]
            && val >= binpic[i+1][j+1]){
                xpk[npk] = j;
                ypk[npk] = i;
                apk[npk] = val;
                npk++;
                if(npk>MAXPEAKS){			/* mostly arbitrary, but needs to be higher
						           for the new camera, i.e. npk>25 */
		  fprintf(stderr, "nk = %d\n", npk);
                    shError("GFINDFIBERS: Invalid flat frame. Too many fibers");
                    return (-1);
                }    
            }
        }
    }
    

#ifdef FINDEBUG
    printf("\nraw xpk,ypk,npk = %d %d %d\n",*xpk,*ypk,npk);    
    for(i=0;i<npk;i++){printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){printf(" %5d",ypk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){printf(" %5d",apk[i]);}
    printf("\n");
#endif
   
    /* there should be nfib good fibers */
    /* reject adjacent maxima. replace the weaker one with the stronger... */
    /* sums values around peaks and keeps the larger sum */
    for(i=0;i<npk;i++){
        for(j=i+1; j<npk;j++){
        	
            if((abs(xpk[j]-xpk[i]) < 2) && (abs(ypk[j]-ypk[i]) < 2)){ 
                 x = xpk[i]; y = ypk[i];
                sumi = binpic[y][x]    + binpic[y+1][x]   + binpic[y-1][x]
                     + binpic[y][x-1]  + binpic[y+1][x-1] + binpic[y-1][x-1]
                     + binpic[y][x+1]  + binpic[y+1][x+1] + binpic[y-1][x+1];
                x = xpk[j]; y = ypk[j];
                sumj = binpic[y][x]    + binpic[y+1][x]   + binpic[y-1][x]
                     + binpic[y][x-1]  + binpic[y+1][x-1] + binpic[y-1][x-1]
                     + binpic[y][x+1]  + binpic[y+1][x+1] + binpic[y-1][x+1];
                if(sumj > sumi){  /* replace i with j */
                    xpk[i] = xpk[j];
                    ypk[i] = ypk[j];
                    apk[i] = apk[j];
                }     
                /* delete j */
                for(k=j+1;k<npk;k++){ 
                    xpk[k-1] = xpk[k];
                    ypk[k-1] = ypk[k];
                    apk[k-1] = apk[k];
                }
                npk--;
               
                
printf("i,j,xi,xj,yi,yj,npt= %d %d %d %d %d %d %d\n",i,j,xpk[i],xpk[j],ypk[i],
        ypk[j],npk);
                       
            }
        }
    }
    

#ifdef FINDEBUG
    printf("\nAfter culling xpk,ypk,npk = %d %d %d\n",*xpk,*ypk,npk);
    for(i=0;i<npk;i++){printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){printf(" %5d",ypk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){printf(" %5d",apk[i]);}
    printf("\n");
#endif

    if(npk > nfib){
        sprintf(pbuf,"GFINDFIBERS: Bad frame: too many fibers: %i > %i", npk, nfib);
        shError(pbuf);
		return(SH_GENERIC_ERROR);
    }
    if(npk < nfib-1){
        sprintf(pbuf,"GFINDFIBERS: Only %d fibers found",npk);
        shError(pbuf);
	/*
	  return(SH_GENERIC_ERROR);
	*/
    }


    /*location of max n0 = n + (X[n+1] - X[n-1])/(2*(2*X[n] - X[n-1] - X[n+1])) 
     *value at max X(n0) = X[n] + (2*X[n] - X[n-1] - X[n+1])/2 * (n0-n)^2 
     *inverse quardatic interpolation to fine maxima, 
     *scaled (by binsize) back to original picture coords */

    for(k=0;k<npk;k++){
        x = xpk[k];
        y = ypk[k];
        xpk[k] = x*BINSIZE + BINSIZE/2 + 
            BINSIZE*(binpic[y][x+1] - binpic[y][x-1]) /
            (2*(2*binpic[y][x] - binpic[y][x-1] - binpic[y][x+1]));
        ypk[k] = y*BINSIZE + BINSIZE/2 +
            BINSIZE*(binpic[y+1][x] - binpic[y-1][x]) /
            (2*(2*binpic[y][x] - binpic[y-1][x] - binpic[y+1][x]));
    }
    

#ifdef FINDEBUG
    printf("\nxlated xpk,ypk,npk = %d %d %d\n",*xpk,*ypk,npk);
    for(i=0;i<npk;i++){ printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",ypk[i]);}
    printf("\n");
#endif

    /* Find the shift to the canonical list: the idea here is 
     * that the coordinate differences for each
     * pair of peaks in the `real' list and the `found' list are computed,
     * and the largest cluster is identified as the true offset. 
     * We bin the data so that there is a clear peak. 
     * With MBINSIZE=8 search over a range of +/- 48 pixels in the original frame, 
     * which is +- 6 bins to avoid overlaps.
     * Bin with maximum counts gives offset with an error of +/- 8 pixels. 
'    * Good enough to uniquely identify fibers.
     */    
     /*CHANGE THIS NOW*/
    for(i=0;i<14;i++){
        dxarr[i] = dyarr[i] = 0;
    }
    for(i=0;i<npk;i++){
        for(j=0;j<nfib;j++){		
            if(abs(dx = (xpk[i] - fiberx[j])) < 48 && 		/*ph,ps: 48 also needs to change, depending on old bundles or new bundles*/
                    abs(dy = (ypk[i] - fibery[j])) < 48){
                dxarr[(dx+56)/8]++;			/*ph,ps: 56 is 7*8, size of dxarr is 14. */
                dyarr[(dy+56)/8]++;
                nposs++;		/* what is this variable for?? */
            }
        }
    }

   
    /* dxarr and dyarr are histograms of peak to fiber distances in binned pixels 
     * Smooth with a 1,2,1 filter and find the maxima.
     */
    xmaxpk=0;
    ymaxpk=0;
    xpkoff=0;
    
    ypkoff=0; 
    for(i=1;i<14-1;i++){	/*ph,ps: whole loop needs to change with length of dxarr/dyarr */
        x = 2*dxarr[i] + (dxarr[i-1] + dxarr[i+1]);
        if(x > xmaxpk){			/* wont happen until i>=7 */
            xmaxpk = x;
            xpkoff = (i-7)*8 + 4 ;	/* gives offset of star center from fiber center */
        }     
        y = 2*dyarr[i] + (dyarr[i-1] + dyarr[i+1]);
        if(y > ymaxpk){
            ymaxpk = y;
            ypkoff = (i-7)*8 + 4;
        }
    }
	
#ifdef FINDEBUG
    printf("\ndxarr,dyarr,xmaxpk,xpkoff,ymaxpk,ypkoff=%d %d %d %d %d %d\n",
	   *dxarr, *dyarr,xmaxpk,xpkoff,ymaxpk,ypkoff);    
    for(i=0;i<14;i++) printf(" %4d",dxarr[i]);  /*needs to change with length of dxarr */  
    printf("\n");
    for(i=0;i<14;i++) printf(" %4d",dyarr[i]);
    printf("\n");
#endif

   
    /* now ypkoff and xpkoff have the offsets (!?) Whew! 
     * the sense is xpk = fiberx + xpkoff 
     */
    
    /* finally, identify the fibers.*/
    nid = 0;
   
    for(i=0;i<npk;i++){
    	//if(plptr->exists[i]){
        fid[i] = 0;
        for(j=0;j<nfib;j++){		
            dx = xpk[i] - fiberx[j] - xpkoff;
            dy = ypk[i] - fibery[j] - ypkoff;
           
            /* ph,ps: below must change with CCD size and fiber spacing */
            if(dx*dx + dy*dy < 512){                /* min BOSS fiber separation: dx~69, dy~34 */
                fid[i] = fiberid[j];		    /* this makes sure the object is less than 16 */
						    /* pixels for both axes away from center of fiber? */
                nid++;				     
            }
        //}
        }
    }
    /* double check to see if any peaks are assigned to a fiber twice. checks to see which 
       fiber is closer to the center guess */
    for(i=0;i<npk;i++){
        for(j=i+1; j<npk;j++){
	    if(fid[i] == fid[j]){
	        int dx1,dx2, dy1,dy2;
		dx1 = xpk[i] - fiberx[i] - xpkoff;
		dy1 = ypk[i] - fibery[i] - ypkoff;
		dx2 = xpk[j] - fiberx[j] - xpkoff;
		dy2 = ypk[j] - fibery[j] - ypkoff;
		if (dx1*dx1+dy1*dy1>dx2*dx2+dy2*dy2){
		    xpk[i] = xpk[j];
		    ypk[i] = ypk[j];
		    fid[i] = fid[j];
		    for(k=i+1;k<npk;k++){ 
		        xpk[k-1] = xpk[k];
			ypk[k-1] = ypk[k];
			apk[k-1] = apk[k];
			fid[k-1] = fid[k];
		    }
		}
		else{
		    xpk[j] = xpk[i];
		    ypk[j] = ypk[i];
		    fid[j] = fid[i];
		    for(k=j+1;k<npk;k++){ 
		        xpk[k-1] = xpk[k];
			ypk[k-1] = ypk[k];
			apk[k-1] = apk[k];
			fid[k-1] = fid[k];
		    }
		}
		npk--;
	    }
	}
    }

#ifdef FINDEBUG    
    printf("\nid'd xpk,ypk,npk = %d %d %d\n",*xpk,*ypk,npk);        
    for(i=0;i<npk;i++){ printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",ypk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",fid[i]);}
    printf("\n");
#endif
   
    /* now the way this was constructed, no two peaks can belong to the
     * same fiber, and all fibers SHOULD be in zeros of the mask.
     * collapse the list to ids and those for which the mask is zero
     * at the peak location
     */
    for(i=0;i<npk;i++){
        if(fid[i] == 0 || mask[ypk[i]][xpk[i]] != 0){ /* throw it away */
	   for(j=i+1;j<npk;j++){
                xpk[j-1] = xpk[j];
                ypk[j-1] = ypk[j];
                fid[j-1] = fid[j];
            }
            npk--;
        }
    }
    
    /* NB!!! we should report npk at this point */ 
    /* at this point, we have a list of locations which are inside
     * holes in the mask and have valid fiber IDs--we need now to find
     * the centers and radii of the holes approximately 
     */

#ifdef FINDEBUG
    printf("\nid'd & zeromask xpk,ypk,npk = %d %d %d\n",*xpk,*ypk,npk);
    for(i=0;i<npk;i++){ printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",ypk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",fid[i]);}
    printf("\n");
#endif    

    for(l=0;l<2;l++){        /* iterate twice */
      for(k=0;k<npk;k++){    /* for each peak */
        for(i=0;i<100;){
            i++;
           
            if(mask[ypk[k]][xpk[k] + i] !=0) break;
        }
        for(j=0;j<100;){
            j++;
            if(mask[ypk[k]][xpk[k] - j] !=0) break;
        }
        xpk[k] += (i-j)/2;
        firad[k] = (i + j)/2;
        for(i=0;i<100;){
            i++;
            if(mask[ypk[k] + i][xpk[k]] !=0) break;
        }
        for(j=0;j<100;){
            j++;
            if(mask[ypk[k] - j][xpk[k]] !=0) break;
        }
        ypk[k] += (i-j)/2;
        firad[k] = (firad[k] + (i+j)/2 )/2; /* avg of x and y 'radii' */
      }
    }
	
#ifdef FINDEBUG
    printf("\nfinal fiber locations &radii, %d peaks\n",npk);        
    for(i=0;i<npk;i++){ printf(" %5d",xpk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",ypk[i]);}
    printf("\n");
    for(i=0;i<npk;i++){ printf(" %5d",firad[i]);}
    printf("\n");
#endif
	
	
    /* done--now find centroids in masked picture and fill out output struct */
   
    for(k=0;k<npk;k++){
    
        sumx = sumy = sum = 0.;
        sump = 0;
     	for(i=ypk[k] - firad[k]; i <= ypk[k] + firad[k]; i++){
            for(j=xpk[k] - firad[k]; j <= xpk[k] + firad[k]; j++){
                val = (flat[i][j] - medn);
		//                wgt  = val>0 ? sqrt(sqrt(val)) : 0;
                wgt  = val>0 ? val : 0;
                sumy += i * wgt;
                sumx += j * wgt;
                sum  += wgt;

                if(val > thresh) sump++;   /* #of pix above thresh, NEEDS to be the same as thresh for pyflat  */
                
            }
        }
        
        
        ptr->g_fid[k] = fid[k];
        ptr->g_xcen[k] = sumx/sum;
        ptr->g_ycen[k] = sumy/sum;
        ptr->g_fibrad[k] = firad[k];
        ptr->g_illrad[k] = sqrt( (float)sump/3.14159 );

        
        /*ph,ps: should be tested for each sized fiber */
        if(firad[k] >= 40){		/*ph,ps: radius size dependent on cartridge, pass in 30 value*/
            sprintf(pbuf,
              "GFINDFIBERS: shError: fiber %d has too large radius, %d pixels",
              k,firad[k]); 
            shError(pbuf);
	    return(SH_GENERIC_ERROR);
        }
    }
   
    /* need to step through loop one more time to add in location of spot fiber, if it exists */
    /* k is from the above loop, one more than the last time through */
    if (plptr->exists[SPOTID]){
    	ptr->g_fid[k] = fiberid[SPOTID];
        ptr->g_xcen[k] = fiberx[SPOTID];
        ptr->g_ycen[k] = fibery[SPOTID];
        ptr->g_fibrad[k] = 17;						/* matches size in masked flat */
        ptr->g_illrad[k] =plptr->radius[SPOTID];	
        npk++;
    }
    
    
    ptr->g_nfibers = npk;
    
// how circular are fibers?

/* generate edge detected centroid xl, yl, compare to xc, yc
   warn if uniformly different but use xl, yl for guiding reference
   add xl, yl to data structures
*/


    sprintf(pbuf,
    "Gguider found %d fibers, peak ct %d:",npk,peak-medn);  
    shDebug(DEBUG, pbuf);
    for(i=0;i<npk;i++){ 
      sprintf(pbuf,"Gguider fiber %2d xcen=%6.2f, ycen=%6.2f",
	      ptr->g_fid[i],
	      ptr->g_xcen[i],
	      ptr->g_ycen[i]); 
      shDebug(DEBUG,pbuf);
      sprintf(pbuf,"Gguider fiber %2d fibrad=%5.1f, illrad=%5.1f",
	      ptr->g_fid[i],
	      ptr->g_fibrad[i],
	      ptr->g_illrad[i]); 
      shDebug(DEBUG,pbuf);
    }
   
    haveflatdata = 1;
    free(binarray);
    return(SH_SUCCESS);
}


/**************Depreciated routines below ****************************************/




/******************* HOTLIST() *****************************************/
/* 
 * this routine finds all pixels above some threshold, and makes a list
 * of their positions. It returns the number of pixels. It does NOT
 * malloc any storage; the lists are declared static above and the
 * max size is defined as MAXBADPIX. It returns the number of bad pixels
 * found, -1 for an error.
 */

/* ph,ps:only needs to be used on hot dark pixals that dont subtract or very bright ones*/

int
hotlist(int xs,             /* x size (pixels) */
        int ys,             /* y size (pixels) */
        S16 **pic,          /* pointer to line pointer array */
        int thresh          /* threshold */
        )
{
    int i, j, n;
    
    n = 0;
    for(i=1; i < ys - 1; i++){
        for(j=1; j < xs - 1; j++){
            if(pic[i][j] > thresh){
                badcolval[n] = j;
                badrowval[n] = i;
                n++;
                if(n >= MAXBADPIX){
		  return(-1);
                }
            }
        }
    }
    return(n);
}

/* helper for hotlist */
/***************** SSHINSORT() *****************************************/
/* sort a short array in place using straight insertion.
 * From Numerical recipes.
 */
 
static void
sshinsort(S16 *arr,int n)
{
    int i,j;
    S16 a;
    
    for(j=1;j<n;j++){
        a = arr[j];
        i = j-1;
        while(i>=0 && arr[i] > a){
            arr[i+1] = arr[i];
            i--;
        }
        arr[i+1] = a;
    }
}        


/******************** HOTFIX() ******************************************/
/* 
 * This routine replaces a list of hot pixels by the medians of the 
 * surrounding pixels.  the row,col coordinates of the bad pixels are
 * stored in the static arrays badrowval, badcolval, and there are 
 * nbadpix<MAXBADPIX of them. 
 */
  
void
hotfix( int xs,              /* x size */
        int ys,              /* y size */
        S16 **pic     /* pointer to line pointer array */
        )
{
    int i, j, n;
    S16 nbrlist[8];   
    int med;
   
    for(n=0; n<nbadpix; n++){
        i = badrowval[n];
        j = badcolval[n];
        if(i > 0 && i < ys - 1 && j > 0 && j < xs - 1){
            nbrlist[0] = pic[i][j-1];
            nbrlist[1] = pic[i][j+1];
            nbrlist[2] = pic[i+1][j-1];
            nbrlist[3] = pic[i+1][j];
            nbrlist[4] = pic[i+1][j+1];
            nbrlist[5] = pic[i-1][j-1];
            nbrlist[6] = pic[i-1][j];
            nbrlist[7] = pic[i-1][j+1];
            sshinsort(nbrlist,8);
            med = (nbrlist[3] + nbrlist[4])/2;
            pic[i][j] = med;
        }
    }
}


