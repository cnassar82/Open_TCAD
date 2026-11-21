/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *************************************************************************/
/*   misc.c                Version 5.1     */
/*   Last Modification : 7/3/91  08:34:05 */

#include <math.h>
#include <float.h>
#include "global.h"
#include "constant.h"
#include "geom.h"

#if FLT_RADIX != 2 || DBL_MANT_DIG < 2*FLT_MANT_DIG || \
    DBL_MIN_EXP > 2*(FLT_MIN_EXP-FLT_MANT_DIG) || DBL_MAX_EXP < 2*FLT_MAX_EXP
#error "Exact stored-coordinate predicates require binary32 products representable in binary64"
#endif
#ifdef __FAST_MATH__
#error "Exact stored-coordinate predicates require error-free floating-point sums"
#endif

/************************************************************************
 *									*
 *	dist( a     ,  b     ) - this function returns the distance 	*
 *  between two points. if b is null, just return magnitude.            *
 *									*
 ************************************************************************/
double dist ( a, b)
float a[MAXDIM], b[MAXDIM];
{
    float tx, ty;

    tx = a[0];
    ty = a[1]; 
    if (b) {
	tx -= b[0];
	ty -= b[1]; 
    } 
    switch( mode ) {
    case ONED :
	return( fabs(tx) );
    case TWOD :
	return ( sqrt( tx*tx + ty*ty ) );
    }
    return(0.0);
}
    

/************************************************************************
 *									*
 *	vol_ele(t) - this routine returns the volume of the element	*
 *									*
 ************************************************************************/
float vol_ele(t)
int t;
{
    switch(mode) {
    case ONED :
	return( fabs( cordinate(pt_nd(vert_tri(t,0)),0) - 
		      cordinate(pt_nd(vert_tri(t,1)),0) ) );
	break;
    case TWOD :
	return( area_tri(
		      cordinate(pt_nd(vert_tri(t,0)),0),
		      cordinate(pt_nd(vert_tri(t,0)),1),
		      cordinate(pt_nd(vert_tri(t,1)),0),
		      cordinate(pt_nd(vert_tri(t,1)),1),
		      cordinate(pt_nd(vert_tri(t,2)),0),
		      cordinate(pt_nd(vert_tri(t,2)),1) ) );
	break;
    }
    return(0.0);
}

/************************************************************************
 *									*
 *	area_tri(a1,a2,b1,b2,c1,c2) - this routine returns the area of 	*
 *  the triangle formed by the coordintes passed.			*
 *									*
 ************************************************************************/
float area_tri(a1,a2,b1,b2,c1,c2)
float a1,a2,b1,b2,c1,c2;
{
    /* Translate before multiplying and promote before subtraction.
     * The expanded float determinant loses small triangle areas at large
     * lateral coordinates, even rejecting points on an actual edge. */
    return 0.5*(((double)b1-a1)*((double)c2-a2)-
                ((double)b2-a2)*((double)c1-a1));
}



/************************************************************************
 *									*
 *	pt_in_tri(x,y,tr) - this routine determins if the pt x,y is in  *
 *  element number tr.							*
 *									*
 ************************************************************************/
pt_in_tri(c,tr)
float c[3];
int tr;
{
    float tx1, tx2, tx3, ty1, ty2, ty3;
    float abtol, artol, tsign;
    float x = c[0];
    float y = c[1];

    switch( mode ) {
    case ONED :
	if ( (c[0] < cordinate(pt_nd(vert_tri(tr,0)),0)) != 
	     (c[0] <= cordinate(pt_nd(vert_tri(tr,1)),0)) )
	    return(TRUE);
	else
	    return(FALSE);
	break;
    case TWOD :
	/*get triangle coordinates*/
	tx1 = cordinate( pt_nd( vert_tri(tr,0) ), 0);
	ty1 = cordinate( pt_nd( vert_tri(tr,0) ), 1);
	tx2 = cordinate( pt_nd( vert_tri(tr,1) ), 0);
	ty2 = cordinate( pt_nd( vert_tri(tr,1) ), 1);
	tx3 = cordinate( pt_nd( vert_tri(tr,2) ), 0);
	ty3 = cordinate( pt_nd( vert_tri(tr,2) ), 1);

	/*calculate tolerances*/
	abtol = 1e-04;
	artol = area_tri(tx1,ty1, tx2,ty2, tx3,ty3);
	tsign = artol/fabs(artol);
	artol = fabs(artol);
	artol = artol * abtol;

	/*begin the testing*/
	if ( (tsign * area_tri( tx1,ty1, tx2,ty2, x,y )) <= -artol)
	    return( FALSE );
	if ( (tsign * area_tri( tx2,ty2, tx3,ty3, x,y )) <= -artol)
	    return( FALSE );
	if ( (tsign * area_tri( tx3,ty3, tx1,ty1, x,y )) <= -artol)
	    return( FALSE );
    
	/*if all that failed, retun true*/
	return( TRUE );
    }
    return(-1);
}


/* Exact orientation sign for stored float coordinates. Each product of two
 * IEEE binary32 coordinates is exact in binary64 (at most 48 significant
 * bits, with exponent range safely inside binary64). Error-free TwoSum
 * accumulates the six products without losing small terms near an origin. */
int orient2d_float_sign(const float *a, const float *b, const float *c)
{
    double terms[6], expansion[6], next[6], q, sum, bv, error;
    int i,j,n=0,m;
    terms[0]=(double)a[0]*b[1]; terms[1]=-(double)a[0]*c[1];
    terms[2]=(double)b[0]*c[1]; terms[3]=-(double)b[0]*a[1];
    terms[4]=(double)c[0]*a[1]; terms[5]=-(double)c[0]*b[1];
    for(i=0;i<6;i++) {
        q=terms[i];m=0;
        for(j=0;j<n;j++) {
            sum=q+expansion[j]; bv=sum-q;
            error=(q-(sum-bv))+(expansion[j]-bv);
            if(error!=0)next[m++]=error;
            q=sum;
        }
        if(q!=0 || m==0)next[m++]=q;
        n=m;
        for(j=0;j<n;j++)expansion[j]=next[j];
    }
    return expansion[n-1]>0 ? 1 : expansion[n-1]<0 ? -1 : 0;
}

