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
/*   rect_tri.c                Version 5.1     */
/*   Last Modification : 7/3/91 15:40:49 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "dbaccess.h"
#include "skel.h"


/*-----------------DIVIDE-----------------------------------------------
 * Choose good ways to cut a region in half.
 * Returns pointers to the edge following each end of the cut.
 *----------------------------------------------------------------------*/
int rect_div(r, lep1, lep2)
struct sreg_str *r;
struct LLedge **lep1,**lep2;
{
    struct LLedge *bp1,*bp2,*bnd, *t1, *t2, *bst1, *bst2;
    double piby2, ang1, ang2, ang3, ang4, intang(), sk_cross();
    double eps = 1.0e-2;
    double eps1, eps2, eps3, eps4, beps;
    int f1;

    piby2 = PI/2;
    bst1=NULL;
    bst2=NULL;
    beps = MAXFLOAT;

  /*...Walk the whole region boundary */
    bnd = r->bnd;
    for (f1=1, bp1=bnd; ((bp1 != bnd) || f1); bp1=bp1->next, f1=0) {
	if ( (eps1 = dabs(bp1->ang - piby2)) < eps ) {

	    /*walk around all the following angles*/
	    for(bp2 = bp1->next; (bp2 != bp1); bp2=bp2->next) {

		if ( (eps2 = dabs(bp2->ang - piby2)) < eps ) {
		    /*test this case to see if we found a cut*/
		    t1 = bp1->prev;
		    t2 = bp2->next;

		    if ( t2->next != t1 ) {
			ang1 = intang(nB(bp2), nF(bp2), nB(t1));
			ang3 = intang(nB(t1), nF(t1), nB(bp2));
			eps1 = fabs(cos(ang1) * ndist(nF(bp2), nB(t1)));
			eps3 = fabs(cos(ang3) * ndist(nF(t1), nB(bp2)));

			if ((eps1>1.0e-8)||(eps3>1.0e-8)) {
			    break;
			}
			else if (sk_cross(r, t1, t2, TRUE) < 0.5) {
			    break;
			}
			else {
			    if (skel_split_valid(t1,t2) && (eps1+eps3) < beps) {
				bst1 = t1;
				bst2 = t2;
				beps = eps1+eps3;
			    }
			}
		    }
		}

		else if ( dabs(bp2->ang - PI) > eps ) break;
	    }
	}
    }

    if ( bst1 != NULL ) {
	*lep1 = bst1;
	*lep2 = bst2;
	return(TRUE);
    }
    else return(FALSE);
}


/*-----------------DVEDGE-----------------------------------------------
 * Divide an edge, using the spacing parameters calculated in dvpram.
 * Called from bfill, and also during the triangulation (triang.c).
 *----------------------------------------------------------------------*/
dvrecedg(r1, r2, ie, t1 )
struct sreg_str *r1, *r2;
int ie;
struct LLedge *t1;
{
    struct LLedge *bp1,*bnd;
    float *p;
    int nump, n;
    int ip, ied, fnd;
    struct sreg_str *r;
    int ret;
    int nbase;
    int axis, first, next;
    double direction;
    float value[MAXDIM];

    n = (r1->len > r2->len)?r1->len:r2->len;
    p = scalloc(float, MAXDIM * n);

    /*which region has the passed edge?*/
    fnd = FALSE;
    for(n=1, bp1 = r1->bnd; n || bp1 != r1->bnd; n = 0, bp1 = bp1->next)
	if ( bp1 == t1 ) fnd = TRUE;

    if ( fnd ) {
	r = r1;
	ret = TRUE;
    }
    else {
	r = r2;
	ret = FALSE;
    }

    /*compute the number of points to add*/
    nump = rect_cnt( r, p, ie );

    /* Opposite-boundary traversal is not necessarily monotone in its
     * projection onto this diagonal. Insert along the original edge, so
     * each split lies in the remaining piece instead of doubling back. */
    axis=fabs((double)cordinate(pt_edg(ie,1),0)-cordinate(pt_edg(ie,0),0)) >=
         fabs((double)cordinate(pt_edg(ie,1),1)-cordinate(pt_edg(ie,0),1)) ? 0:1;
    direction=cordinate(pt_edg(ie,1),axis)>=cordinate(pt_edg(ie,0),axis)?1:-1;
    for(first=1;first<nump;first++) {
        for(n=0;n<MAXDIM;n++)value[n]=p[first*MAXDIM+n];
        for(next=first;next>0 &&
            direction*(value[axis]-p[(next-1)*MAXDIM+axis])<0;next--)
            for(n=0;n<MAXDIM;n++)p[next*MAXDIM+n]=p[(next-1)*MAXDIM+n];
        for(n=0;n<MAXDIM;n++)p[next*MAXDIM+n]=value[n];
    }

    /*if this infinite loops, big trouble*/
    for( bnd = r->bnd; bnd->edge != ie; bnd=bnd->next );

#define min2(A,B) ((A)<(B))?(A):(B)
    if (nump != 0) {
	nbase = nd_edg(ie,1);

	for(n = 0; n < nump; n++) {
	    ied = sp_edge(ie,&(p[n*MAXDIM]),&ip,FALSE);
	    if (ied<0)continue;
	    if ( (nd_edg(ie,0) != nbase) && (nd_edg(ie,1) != nbase) )
		ie = ied;
	}
    }

    free(p);
    return(ret);
}



rect_cnt(r, p, ie)
struct sreg_str *r;
float *p;
int ie;
{
    struct LLedge *bp1, *bnd;
    float x[MAXDIM], dx[MAXDIM], q[MAXDIM], dq[MAXDIM], alp[MAXDIM];
    int nump, n;

    /*if this infinite loops, big trouble*/
    for( bnd = r->bnd; bnd->edge != ie; bnd=bnd->next );

    x[0]  = cordinate( pt_edg(ie,0), 0);
    x[1]  = cordinate( pt_edg(ie,0), 1);
    dx[0] = cordinate( pt_edg(ie,1), 0) - x[0];
    dx[1] = cordinate( pt_edg(ie,1), 1) - x[1];
    dq[0] = -dx[1];
    dq[1] = dx[0];

    /*walk around the rest of the edges*/
    nump = 0;
    for( bp1=bnd->next; (bp1->next != bnd); bp1 = bp1->next ) {
	if ( dabs(bp1->ang - PI) < 1.0e-3 ) {
	    n = pt_nd(nB(bp1));

	    q[0] = cordinate( n, 0 );
	    q[1] = cordinate( n, 1 );
	    if (lil(x, dx, q, dq, alp)) panic("rect_cnt: parallel lines?!?");
	    if ((alp[0]<0.0)||(alp[0]>1.0)) panic("new node not on old edge");
	    q[0] = alp[0] * dx[0] + x[0];
	    q[1] = alp[0] * dx[1] + x[1];

	    /*check to make sure this point isn't on the region boundary*/
	    if ( chk_pt( q, r ) ) {
		p[nump*MAXDIM]   = q[0];
		p[nump*MAXDIM+1] = q[1];
		nump++;
	    }
	}
    }

    return(nump);
}


chk_pt(q, r)
struct sreg_str *r;
float *q;
{
    struct LLedge *bp1;
    int p, f;

    for(bp1 = r->bnd, f=1; (bp1 != r->bnd) || f; bp1 = bp1->next, f=0) {
	p = pt_nd(nB( bp1 ));
	if ( dist(cord_arr(p), q) < EPS ) return(0);
    }
    return(1);
}




rect_tri( r, rbn )
int r;
int rbn;
{
    struct LLedge *bp, *quadsplit();

    while ( sreg[r]->len > 4 ) {

        /* Angle-only clipping can remove the sole noncollinear vertex,
         * leaving a flat remainder. Keep every boundary vertex and require
         * an interior ear with no other vertex on or inside its triangle. */
        bp = geometric_ear(sreg[r]);
        if (!bp) panic("No interior ear in rectangular triangulation");
        (void)sp_reg(bp->prev, bp->next, NULL);
        if (sreg[nsreg-1]->len != 3) panic("Invalid rectangular ear split");
        (void)cr_tri(rbn, nsreg-1);
        free_skel(nsreg-1);
    }

    if ( sreg[r]->len != 4 ) panic("supposed to quit with quads!!");
    bp = quadsplit(sreg[r]);
    if (!bp) panic("No interior diagonal in rectangular triangulation");
    (void)sp_reg(bp->prev, bp->next, NULL);
    (void)cr_tri(rbn, nsreg-1);
    free_skel(nsreg-1);
    (void)cr_tri(rbn, r);
}
