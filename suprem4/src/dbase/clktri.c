/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *     Copyright c 1991 The board of regents of the University of 	 *
 *                      Florida.  All rights reserved.			 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of the University of 	 *
 *     Florida.								 *
 *									 *
 *									 *
 *************************************************************************/
/*   clktri.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:20:37 */

#include <stdio.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"



/************************************************************************
 *									*
 *	clock_tri() - This routine makes sure the triangles are clock	*
 *  wise.								*
 *									*
 *      Original     	MEL	90	Based on meshgen		*
 *									*
 ************************************************************************/
clock_tri()
{
    register int ie;

    /*preform this check for everybody...*/
    for(ie = 0; ie < ne; ie++) {
	if ( ! ask(tri[ie], CLKWS) ) sclk_tri(ie);
    }
}


sclk_tri(ie)
int ie;
{
    register int i, j, k, tt;
    float xi, xj;

    set(tri[ie], CLKWS);
    switch( mode ) {
    case ONED :
	/*get the node numbers*/
	i = tri[ie]->nd[0];
	j = tri[ie]->nd[1];

	/*get the coordinate terms*/
	xi = pt[ nd[i]->pt ]->cord[0];
	xj = pt[ nd[j]->pt ]->cord[0];

	if ((xi - xj) < 0.0) {
	    tri[ie]->nd[0] = j;
	    tri[ie]->nd[1] = i;
	    tt = tri[ie]->nb[1];
	    tri[ie]->nb[1] = tri[ie]->nb[0];
	    tri[ie]->nb[0] = tt;
	}
	break;
    
    case TWOD :
	/*get the point numbers*/
	i = tri[ie]->nd[0];
	j = tri[ie]->nd[1];
	k = tri[ie]->nd[2];

	/* Orientation is decided by the exact predicate on the stored float
	 * coordinates. A float determinant loses the sign of valid slivers
	 * (e.g. two vertices one ULP apart), leaving an order that tri_geom's
	 * more accurate evaluation then rejects as corrupt. Positive exact
	 * orient2d(i,j,k) is this database's clockwise (y-down) convention. */
	if (orient2d_float_sign(pt[ nd[i]->pt ]->cord, pt[ nd[j]->pt ]->cord,
	                        pt[ nd[k]->pt ]->cord) < 0) {
	    tri[ie]->nd[1] = k;
	    tri[ie]->nd[2] = j;
	    /*maintain neighbor and edge ordering*/
	    tt = tri[ie]->nb[1];
	    tri[ie]->nb[1] = tri[ie]->nb[2];
	    tri[ie]->nb[2] = tt;
	    tt = tri[ie]->edg[1];
	    tri[ie]->edg[1] = tri[ie]->edg[2];
	    tri[ie]->edg[2] = tt;
	}
    }
}
