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
/*   flip.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:20:46 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <math.h>

#include "global.h"
#include "constant.h"
#include "geom.h"
#include "conservative_transfer.h"
#include "mesh_trial.h"

static double flip_area(int a,int b,int c)
{
    float *p=cord_arr(pt_nd(a)),*q=cord_arr(pt_nd(b)),*r=cord_arr(pt_nd(c));
    return .5*(((double)q[0]-p[0])*((double)r[1]-p[1])-
               ((double)q[1]-p[1])*((double)r[0]-p[0]));
}

/* A diagonal swap changes P1 integration weights even with fixed vertices.
 * Prepare a bounded conservative transfer before changing connectivity.
 * Include each affected node's entire support, not just this quadrilateral:
 * its single concentration also contributes to unchanged neighboring cells.
 */
static int flip_fields(int *nodes,int tp,int tq,double out[4][MAXIMP])
{
    double oldw[4]={0},neww[4],values[4],fraction[4][MAXIMP];
    /* Edge endpoint order is arbitrary; sclk_tri orders committed children. */
    double a=fabs(flip_area(nodes[0],nodes[1],nodes[2]))/3;
    double b=fabs(flip_area(nodes[0],nodes[2],nodes[3]))/3;
    int i,j,k,t,parent,stream,material=mat_nd(nodes[0]);
    if(!(a>0 && b>0))return 0;
    for(i=0;i<4;i++) {
        if(mat_nd(nodes[i])!=material)return 0;
        neww[i]=0;
        for(k=0;k<num_tri_nd(nodes[i]);k++) {
            t=tri_nd(nodes[i],k);
            if(dead_tri(t) || !leaf(tri[t]))continue;
            oldw[i]+=flip_area(vert_tri(t,0),vert_tri(t,1),vert_tri(t,2))/3;
            if(t!=tp && t!=tq)neww[i]+=flip_area(vert_tri(t,0),vert_tri(t,1),vert_tri(t,2))/3;
        }
        neww[i]+=i==0 || i==2 ? a+b:i==1 ? a:b;
        for(j=0;j<n_imp;j++) {
            out[i][j]=sol_nd(nodes[i],j);
            if(!isfinite(out[i][j]))return 0;
        }
    }
    for(j=0;j<n_imp;j++) {
        parent=conserved_parent(soltoimp[j]);
        if(parent<0)continue;
        stream=imptosol[parent];if(stream<0)return 0;
        for(i=0;i<4;i++) {
            double total=out[i][stream],active=out[i][j];
            if(total<0 || active<0 || active>total)return 0;
            fraction[i][j]=total>0?active/total:0;
        }
    }
    for(j=0;j<n_imp;j++)if(conserved_density(soltoimp[j],material)) {
        for(i=0;i<4;i++)values[i]=out[i][j];
        if(!conserved_transfer(4,oldw,neww,values))return 0;
        for(i=0;i<4;i++)out[i][j]=values[i];
    }
    for(j=0;j<n_imp;j++) {
        parent=conserved_parent(soltoimp[j]);if(parent<0)continue;
        stream=imptosol[parent];
        for(i=0;i<4;i++)out[i][j]=fraction[i][j]*out[i][stream];
    }
    return 1;
}


/*-----------------FLIP-------------------------------------------------
 * Triangle flipper - makes triangulation Delauney.
 * Returns number of triangles flipped or -1 for error.
 *----------------------------------------------------------------------*/
int flip()
{
    int it,j,count,flips=0,all_done,*todo,*done;
    struct tms before, after;
    char buf[100];
    int yes;

    /*no flipping required in one dimension*/
    if ( mode != TWOD ) return(0);

    times(&before);

    /*...Need to remember which points have been done, which to do. */
    if (!(
	(todo = scalloc(int, 1+ned)) &&
	(done = scalloc(int, 1+ned))
       )) return(-1);

    for (it=0; it < ned; it++) todo[it] = TRUE;

  /*...Smoothing loop. */
    count = 0;
    do {
	all_done = TRUE;
	for (it=0; it < ned; it++) {
	    done[it] = !todo[it];
	    todo[it] = FALSE;
	    }

	for(it = 0; it < ned; it++) {
	    /*if we've checked it*/
	    if (done[it]) continue;

	    /*if it belongs only to one element*/
	    if (num_tri_edg(it) != 2) continue;

	    /*make sure the two elements are in the same region*/
	    if ( reg_tri(tri_edg(it,0)) != reg_tri(tri_edg(it,1)) ) continue;

	    if ( gimme_ehed(it) < 0 ) {
		if ((yes = ck_flip(it)) != 1) continue;

		if(!tri_flip( it, tri_edg(it,0), tri_edg(it,1) ))continue;
		flips++;
		all_done = FALSE;

		/*mark all neighbor edges as needing checking*/
		for(j = 0; j < num_edge_nd( nd_edg(it,0) ); j++)
		    todo[ edge_nd( nd_edg(it,0), j) ] = TRUE;
		for(j = 0; j < num_edge_nd( nd_edg(it,1) ); j++)
		    todo[ edge_nd( nd_edg(it,1), j) ] = TRUE;
	    }
	}
    } while (!all_done && ++count < 30);

    free(todo);	free(done);

    times(&after);
    sprintf(buf, "obtuse triangle flip (%d)", flips);
    print_time( buf, &before, &after);
    return(flips);
}




/*-----------------DO_FLIP----------------------------------------------
 * Mechanics of flip.
 *----------------------------------------------------------------------*/
tri_flip( it, tp, tq )
int it, tp, tq;
{
    int i, n1, n2, jp, jq, tr;
    int sp[4], sn[6], se[6], fixme;
    double fields[4][MAXIMP];

    /*...Find which points are unique. */
    jp = twhich (tp, tq);
    jq = twhich (tq, tp);

    /*get the four nodes the form the quad*/
    sp[0] = vert_tri(tp, jp);
    sp[1] = nd_edg(it,0);
    sp[2] = vert_tri(tq, jq);
    sp[3] = nd_edg(it,1);

    /* In an implicit moving-mesh solve, old/mid/new concentrations and
     * their respective integration weights carry the inventory. nd[].sol
     * remains the committed input until that solve succeeds. Rezone its
     * fields only outside that time-integration transaction. */
    if(n_imp && !mesh_trial_active() && !flip_fields(sp,tp,tq,fields))return 0;

    /*change the edge in question*/
    sub_edge_nd(nd_edg(it,0), it );
    sub_edge_nd(nd_edg(it,1), it );
    n1 = vert_tri(tp, jp);
    edg[it]->nd[0] = n1;
    add_edge_nd(n1, it);
    n2 = vert_tri(tq, jq);
    edg[it]->nd[1] = n2;
    add_edge_nd(n2, it);

    /*...Copy points out - order is guaranteed to remain c-clockwise. */
    sn[0] = tri[tq]->nb[nwhich(tq,sp[3])];
    sn[1] = tq;
    sn[2] = tri[tp]->nb[nwhich(tp,sp[3])];
    sn[3] = tri[tp]->nb[nwhich(tp,sp[1])];
    sn[4] = tp;
    sn[5] = tri[tq]->nb[nwhich(tq,sp[1])];

    se[0] = tri[tq]->edg[nwhich(tq,sp[3])];
    se[1] = it;
    se[2] = tri[tp]->edg[nwhich(tp,sp[3])];
    se[3] = tri[tp]->edg[nwhich(tp,sp[1])];
    se[4] = it;
    se[5] = tri[tq]->edg[nwhich(tq,sp[1])];

    /*...and back in. */
    for (i=0; i<3; i++) {
	tri[tp]->nd[i] = sp[i];
	tri[tq]->nd[i] = sp[(i+2)%4];
	tri[tp]->nb[i] = sn[i];
	tri[tq]->nb[i] = sn[i+3];
	tri[tp]->edg[i] = se[i];
	tri[tq]->edg[i] = se[i+3];
    }

    /*update the edge to element pointers*/
    sub_ele_edg(se[3], tp);
    sub_ele_edg(se[0], tq);
    add_ele_edg(se[0], tp);
    add_ele_edg(se[3], tq);


    /*update the neighbors of the triangles*/
    for (tr=tp; tr==tp || tr==tq; tr += tq-tp) { /* for each of tp, tq */
	if ((fixme = tri[tr]->nb[0]) < 0) continue;
	for (i = 0; i < 3; i++)
	    if (tri[fixme]->nb[i] == tp+tq-tr) tri[fixme]->nb[i] = tr;
    }

    /*now update the nd to tri list*/
    add_tri_nd(n1, tq);
    add_tri_nd(n2, tp);
    sub_tri_nd(vert_tri(tp, twhich(tp,tq)), tq);
    sub_tri_nd(vert_tri(tq, twhich(tq,tp)), tp);

    /*repair the geometry*/
    sclk_tri(tp);
    do_geom(tp);
    sclk_tri(tq);
    do_geom(tq);
    for(i=0;i<4;i++) {
        int j;
        if(!mesh_trial_active())for(j=0;j<n_imp;j++)nd[sp[i]]->sol[j]=fields[i][j];
    }
    return 1;
}


/*-----------------CK_FLIP----------------------------------------------
* Check the right situation for flip.
* Returns 1 for flipped or -1 for error.
*----------------------------------------------------------------------*/
int ck_flip(we)
int we;
{
    int t1, t2, ix1, ix2;
    int n1, n2, n3, n4, p1, p2, p3, p4;
    int s1, s2, s3, s4;

    t1  = tri_edg(we,0);
    t2  = tri_edg(we,1);
    ix1 = ewhich(t1,we);
    ix2 = ewhich(t2,we);

    n1 = vert_tri(t1,ix1);
    n2 = vert_tri(t2,ix2);
    n3 = nd_edg(we,0);
    n4 = nd_edg(we,1);

    p1 = pt_nd(n1);
    p2 = pt_nd(n2);
    p3 = pt_nd(n3);
    p4 = pt_nd(n4);

    /* Both diagonals must cross strictly in the stored coordinates.
     * Slope/intercept arithmetic loses the side of a near-collinear point:
     * even its products were evaluated as floats before promotion. A false
     * refusal leaves an avoidable sliver; a false acceptance changes domain.
     * Use the exact predicate, including horizontal/vertical diagonals. */
    s1 = orient2d_float_sign(pt[p1]->cord, pt[p2]->cord, pt[p3]->cord);
    s2 = orient2d_float_sign(pt[p1]->cord, pt[p2]->cord, pt[p4]->cord);
    s3 = orient2d_float_sign(pt[p3]->cord, pt[p4]->cord, pt[p1]->cord);
    s4 = orient2d_float_sign(pt[p3]->cord, pt[p4]->cord, pt[p2]->cord);
    return s1*s2 < 0 && s3*s4 < 0;
}




