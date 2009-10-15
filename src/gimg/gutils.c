#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* use ours */
#ifdef u_short
#undef u_short
#endif
#define u_short unsigned short int

#ifdef PI
#undef PI
#endif
#define PI 3.1415926

#define ERROR puts /* you may want to change this to some other error
                        logging routine */

/*
 * This routine extracts a u_short square patch from a u_short parent 
 * picture. Data are extracted from a circular disk of radius r which
 * must be no greater than half the (odd) patch size if you want the
 * results to be nice, though no check is made. NB!!!! NO CHECK IS
 * MADE TO BE SURE THE PATCH FITS IN THE PARENT PICTURE--the routine
 * does not even know the size of the parent picture. You must ensure
 * that this is so; the code will segfault if it is not. Likewise you
 * must ensure that the patch array is big enough. We insist for
 * nicety that the patch edge length is odd, though again no check is
 * made.
 */

void 
cextract(p,cp,xc,yc,r,cpsiz)
    u_short **p;      /* parent picture */
    u_short **cp;     /* patch (square) */
    int cpsiz;        /* edge length of patch (odd) */
    int xc;                 
    int yc;           /* extraction center of patch */
    int r;            /* extraction radius; <= cpsize/2 */
{
    int i,j;
    int cphs = cpsiz/2;
    int r2;
    int R2 = r*r;
    
    for( i= -cphs ; i <= cphs ; i++ ){
        for( j= -cphs; j <= cphs; j++ ){
           r2 = i*i + j*j;
           if( r2 <= R2 ){
                cp[j + cphs ][i + cphs ] = p[ j + yc ][ i + xc ] ;
            } else {
                cp[ j + cphs ][ i + cphs ] = 0;
            }
        }
    }
}



void 
cmextract(p,cp,xc,yc,r,cpsiz)
    char **p;         /* parent mask */
    char **cp;        /* patch mask (square) */
    int cpsiz;        /* edge length of patch (odd) */
    int xc;                 
    int yc;           /* extraction center of patch */
    int r;            /* extraction radius; <= cpsize/2 */
{
    int i,j;
    int cphs = cpsiz/2;
    int r2;
    int R2 = r*r;
    
    for( i= -cphs ; i <= cphs ; i++ ){
        for( j= -cphs; j <= cphs; j++ ){
           r2 = i*i + j*j;
           if( r2 <= R2 ){
                cp[j + cphs ][i + cphs ] = p[ j + yc ][ i + xc ] ;
            } else {
                cp[ j + cphs ][ i + cphs ] = 0;
            }
        }
    }
}



/* 
 * this rotates the IMAGE with respect to the PIXELS counterclockwise
 * by theta (deg). Uses bilinear linear interpolation.
 */
 
void  
grot(			   /* general array rotator; puts zeros where 
                                 data are not generated */
     float theta,              /* angle in degrees */
     u_short **p,              /* source array */
     u_short **dp,             /* destination array */
     int xsz, int ysz)              /* sizes */
{
    int sthet;
    int cthet;
    int i,j;
    int xrf,yrf,xr,yr;
    int x2 = xsz/2;
    int y2 = ysz/2;
    int xroff = x2*8192;
    int yroff = y2*8192;
    int xrm = (xsz-1)*8192;
    int yrm = (ysz-1)*8192;
    int xrem,yrem,xrem1,yrem1;
    double thetrad = PI * theta / 180.;
    
    cthet = 8192.*cos(thetrad);
    sthet = 8192.*sin(thetrad);
    
    for(i=0;i<ysz;i++){
        for(j=0;j<xsz;j++){
            xrf = (j-x2)*cthet + (i-y2)*sthet + xroff;
            yrf = (i-y2)*cthet - (j-x2)*sthet + yroff;
            if(xrf < 0 || xrf >= xrm || yrf < 0 || yrf >= yrm){
                dp[i][j] = 0;
                continue;
            }
            /* have pixel in interior; now interpolate */
            xr = xrf>>13;
            yr = yrf>>13;
            xrem = (xrf & 8191)>>5;  /* xrem <= 255  */
            yrem = (yrf & 8191)>>5;  /* yrem <= 255  */
            xrem1 = 256-xrem;
            yrem1 = 256-yrem;
            
            dp[i][j] = (  p[yr][xr]*xrem1*yrem1
                        + p[yr+1][xr]*yrem*xrem1
                        + p[yr][xr+1]*yrem1*xrem
                        + p[yr+1][xr+1]*xrem*yrem)>>16;
        }
    }
}

void  
maskrot(theta,p,dp,xsz,ysz)   /* general byte-sized mask rotator puts zeros 
                               * where data are not generated */
    float theta;              /* angle in degrees */
    char **p;                 /* source mask array */
    char **dp;                /* destination mask array */
    int xsz,ysz;              /* sizes */
{
    int sthet;
    int cthet;
    int i,j;
    int xrf,yrf,xr,yr;
    int x2 = xsz/2;
    int y2 = ysz/2;
    int xroff = x2*256;
    int yroff = y2*256;
    int xrm = (xsz-1)*256;
    int yrm = (ysz-1)*256;
    double thetrad = PI * theta / 180.;
    
    cthet = 256.*cos(thetrad);
    sthet = 256.*sin(thetrad);
    
    for(i=0;i<ysz;i++){
        for(j=0;j<xsz;j++){
            xrf = (j-x2)*cthet + (i-y2)*sthet + xroff;
            yrf = (i-y2)*cthet - (j-x2)*sthet + yroff;
            /* scaled coordinates in rotated frame */
            if(xrf < 0 || xrf >= xrm || yrf < 0 || yrf >= yrm){
                dp[i][j] = 0;
                continue;
            }
            /* have pixel in interior; now find nearest pixel in dest. */
            /* first shift by 7 bits, so have twice coordinates; if
             * twice is odd, round up, if even, round down */
            xr = xrf>>7;
            xr = ((xr % 2) == 1) ? xr/2 + 1 : xr/2 ;
            yr = yrf>>7;
            yr = ((yr % 2) == 1) ? yr/2 + 1 : yr/2 ;
            
            dp[i][j] = p[yr][xr];
        }
    }
}



int 
cexrot(p, cp, cpsiz, xc, yc, r, theta)
    unsigned short **p;      /* parent picture */
    unsigned short **cp;     /* patch (square) */
    int cpsiz;               /* edge length of patch (odd) */
    int xc;                  
    int yc;                  /* extraction center of patch */
    int r;                   /* extraction radius; <= cpsize/2 */
    double theta;            /* rotation angle in degrees */
{
    u_short **p1;     /* intermediate array */
    int i;
    
    p1 = (u_short **)malloc(cpsiz * (sizeof(u_short*) + cpsiz*sizeof(u_short)));
    if(p1 == (u_short **)NULL){
        ERROR("\ncexrot:Cannot allocate memory for intermediate patch storage\n");
        return (1);
    }


    
    /* make line pointers */
    for( i=0 ; i < cpsiz; i++){
        p1[i] = (u_short *)((char *)p1 + cpsiz*(sizeof(u_short *) + i*sizeof(u_short)));
    }
    /*printf("\np1=%u,p1[0]=%u,p1[cpsiz-1]=%u\n",p1,p1[0],p1[cpsiz-1]); debug*/
    cextract(p,p1,xc,yc,r,cpsiz);
    grot(theta,p1,cp,cpsiz,cpsiz); 
    free(p1);
    return (0);
}


int 
cmexrot(p, cp, cpsiz, xc, yc, r, theta)
    char **p;                /* parent mask */
    char **cp;               /* mask for patch (square) */
    int cpsiz;               /* edge length of patch (odd) */
    int xc;                  
    int yc;                  /* extraction center of patch */
    int r;                   /* extraction radius; <= cpsize/2 */
    double theta;            /* rotation angle in degrees */
{
    char **p1;               /* intermediate array */
    int i;
    
    p1 = (char **)malloc(cpsiz * (sizeof(char*) + cpsiz*sizeof(char)));
    if(p1 == (char **)NULL){
        ERROR("\nmaskexrot:Cannot allocate memory for intermediate patch storage\n");
        return (1);
    }

    
    /* make line pointers */
    for( i=0 ; i < cpsiz; i++){
        p1[i] = (char *)((char *)p1 + cpsiz*(sizeof(char *) + i*sizeof(char)));
    }
    /*printf("\np1=%u,p1[0]=%u,p1[cpsiz-1]=%u\n",p1,p1[0],p1[cpsiz-1]); debug*/
    cmextract(p,p1,xc,yc,r,cpsiz);
    maskrot(theta,p1,cp,cpsiz,cpsiz); 
    free(p1);
    return (0);
}


