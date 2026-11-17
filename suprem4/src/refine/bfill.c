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
/*   bfill.c                Version 5.1     */
/*   Last Modification : 7/3/91 15:40:34 */

#include <stdio.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "dbaccess.h"
#include "skel.h"

static int refinement_preserves_boundary(int edge,float *c);

/*-----------------DVEDGE-----------------------------------------------
 * Divide an edge
 *----------------------------------------------------------------------*/
div_edg(ie)
int ie;
{
    double lej, ha, hb, rr, ff;
    int nn, inew, i;
    double updn = 10.0;
    float c[MAXDIM];

    lej = len_edg(ie);
    ha  = get_space(pt_edg(ie,0));
    hb  = get_space(pt_edg(ie,1));

    if ( (lej / ha > updn) && (lej / hb > updn) ) {
	/*split the edge in half, than treat the pieces...*/
	for(i = 0; i < mode; i++) 
    c[i] = 0.5 * ((double)cordinate(pt_edg(ie,0),i)+cordinate(pt_edg(ie,1),i));
	if(mode==TWOD && !refinement_preserves_boundary(ie,c)) return;
	inew = sp_edge(ie, c, &nn, FALSE);
	/* A negative result identifies an existing endpoint, not an edge. */
	if (inew < 0) return;
	ha  = get_space(pt_edg(ie,0));
	hb  = get_space(pt_edg(ie,1));
	dvpram(ha, hb, len_edg(ie), &nn, &rr, &ff);
	dvedge(ie, len_edg(ie), nn, rr, ff);

	ha  = get_space(pt_edg(inew,0));
	hb  = get_space(pt_edg(inew,1));
	dvpram(ha, hb, len_edg(inew), &nn, &rr, &ff);
	dvedge(inew, len_edg(inew), nn, rr, ff);
    }
    else {
	dvpram(ha, hb, lej, &nn, &rr, &ff);
	dvedge(ie, lej, nn, rr, ff);
    }
}




/*-----------------DVEDGE-----------------------------------------------
 * Divide an edge, using the spacing parameters calculated in dvpram.
 * Called from bfill, and also during the triangulation (triang.c).
 *----------------------------------------------------------------------*/
static double refinement_side(double ax,double ay,double bx,double by,
                              double cx,double cy)
{
    return (bx-ax)*(cy-ay)-(by-ay)*(cx-ax);
}

static int refinement_crosses(float *a,float *b,int p,int q)
{
    double px=cordinate(p,0),py=cordinate(p,1);
    double qx=cordinate(q,0),qy=cordinate(q,1);
    double s1=refinement_side(a[0],a[1],b[0],b[1],px,py);
    double s2=refinement_side(a[0],a[1],b[0],b[1],qx,qy);
    double s3=refinement_side(px,py,qx,qy,a[0],a[1]);
    double s4=refinement_side(px,py,qx,qy,b[0],b[1]);
    return ((s1<0 && s2>0)||(s1>0 && s2<0)) &&
           ((s3<0 && s4>0)||(s3>0 && s4<0));
}

/* Refinement may only insert representable points which preserve the
 * existing boundary. Even a correctly rounded float on a mathematical
 * segment can cross the other wall of a sub-ulp-wide polygon. In that
 * case leave the original edge intact; never change its topology merely
 * to satisfy a requested spacing. No physical length cutoff is used. */
static int refinement_preserves_boundary(int edge,float *c)
{
    int i,sr;
    struct LLedge *b,*start;
    float *a=cord_arr(pt_edg(edge,0)),*d=cord_arr(pt_edg(edge,1));
    for(i=0;i<num_skel_edg(edge);i++) {
        sr=skel_edg(edge,i);start=b=sreg[sr]->bnd;
        do {
            if(b->edge!=edge &&
               (refinement_crosses(a,c,pt_nd(nB(b)),pt_nd(nF(b))) ||
                refinement_crosses(c,d,pt_nd(nB(b)),pt_nd(nF(b)))))
                return FALSE;
            b=b->next;
        }while(b!=start);
    }
    return TRUE;
}

int dvedge (ie, lej, nnew, ratio, first)
int ie;			/* Edge to divide. */
int nnew;		/* Number of new nodes to add */
double lej;		/* Length of edge. */
double ratio;		/* Ratio of each interval to next. */
double first;		/* First interval as a fraction of lej. */
{
    int inew, ip, ied;
    int sp_edge();
    double dx,dy,h,frac,x0,y0,position;
    float c[MAXDIM];
 
    /*if we don't need to add any, we don't need to add any*/
    if (nnew==0) return;

    x0 = cordinate( pt_edg(ie,0), 0);
    y0 = cordinate( pt_edg(ie,0), 1);
    dx = (double)cordinate( pt_edg(ie,1), 0) - x0;
    dy = (double)cordinate( pt_edg(ie,1), 1) - y0;
    position = 0;
    frac = first;
    h =  frac*lej;
    if (ratio < 1) h *= ratio;

    for (inew=1; inew<=nnew; inew++) {

	/* Compute each point from the original segment in double precision.
	 * Accumulating rounded float increments bends a straight boundary and
	 * can cross the opposite wall of a thin retained region. */
	position += frac;
	c[0] = x0 + position*dx;
	c[1] = y0 + position*dy;
	if(mode==TWOD && !refinement_preserves_boundary(ie,c)) {
	    h *= ratio;
	    frac *= ratio;
	    continue;
	}
	
	ied = sp_edge(ie,c,&ip,FALSE); 
	set_space(ip, h);
	/* Splits below MIN_GRID_SPACE can return an endpoint.  Keep
	 * walking the current edge; never use that status as an edge ID. */
	if (ied == -1) break;
	if (ied < 0) {
	    h *= ratio;
	    frac *= ratio;
	    continue;
	}

	if (pt_edg(ie,0) != ip) ie = ied;

	h *= ratio;
	frac *= ratio;
    }

    return;
}
