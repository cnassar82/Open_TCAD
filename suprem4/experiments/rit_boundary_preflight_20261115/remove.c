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
/*   remove.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:21:13 */

#include <stdio.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "dbaccess.h"
#include "material.h"
#include "impurity.h"

/************************************************************************
 *									*
 *	rem_nd( n ) - This routine removes a node and eliminates	*
 * all the data associated with it. THIS IS NOT A SAFE ENTRY POINT!!!	*
 *									*
 * Original:	MEL	8/90						*
 *									*
 ************************************************************************/
rem_1nd( n )
int n;
{
    register int p, n1;

    if ( dead_nd(n) ) return;

    /*remove the node from the point list*/
    p = pt_nd(n);
    sub_nd_pt(n, p);

    /*discard the node*/
    fing_nd(n);

    /*remove the parent node if this is the last node...*/
    if ( num_nd(p) == 0 ) fing_pt(p);

    /*if only a gas node left, kill the node and point*/
    if ( num_nd(p) == 1 ) {
	n1 = nd_pt(p, 0);
	if ( mat_nd(n1) == GAS ) rem_1nd(n1);
    }
}
	




/************************************************************************
 *									*
 *	rem_ele( t ) - This routine removes an element and eliminates	*
 * all the data associated with it. This does not maintain region	*
 * information.								*
 *									*
 * Original:	MEL	8/90						*
 *									*
 ************************************************************************/
rem_ele( t, expos )
int t;
{
    register int i, n, j, k, n1;

    /*start off with the nodes*/
    for(i = 0; i < num_vert(t); i++) {
	n = vert_tri(t,i);
	sub_tri_nd(n, t);
	if (num_tri_nd(n) == 0) rem_1nd(n);
    }

    /*fix up the neighbor triangles*/
    for(i = 0; i < num_face(t); i++) {
	n = neigh_fc(t, i);

	/*if there is a neighbor to fix*/
	if ( n >= 0 ) {
	    j = twhich(n, t);

	    /*fix the boundary code of this location*/
	    if ( expos ) {
		set_face_exp( n, j );

		for(k = 0; k < num_nd_fc(n, j); k++) {
		    n1 = nd_face_ele(n,j,k);
		    set_surf( pt_nd(n1) );
		    /*add a gas node if one doesn't exist*/
		    if( node_mat(n1, GAS) == -1 ) (void)mk_nd( pt_nd(n1), GAS );
		}
	    }

	    /*fix the boundary code of this location*/
	    tri[n]->nb[j] = BC_OFFSET;
	    clr( tri[n], NEIGHB );
	    neigh_dirty = TRUE;
	}
    }

    /*finally, discard the triangle*/
    fing_tri(t);
    sub_tri_reg(reg_tri(t),t);
    if (num_tri_reg(reg_tri(t)) == 0) fing_reg(reg_tri(t));
    edge_dirty = TRUE;
}
	


/************************************************************************
 *									*
 *	rem_ele( t ) - This routine removes an element and eliminates	*
 * all the data associated with it. This does not maintain region	*
 * information.								*
 *									*
 * Original:	MEL	8/90						*
 *									*
 ************************************************************************/
rem_reg( r )
int r;
{
    register int i;

    /*eliminate all the triangles*/
    while(num_tri_reg(r) > 0) rem_ele( tri_reg(r, 0), TRUE );
    fing_reg(r);
}


#define mag(A) ((mode==ONED)?(sqrt(A[0]*A[0])):sqrt(A[0]*A[0]+A[1]*A[1]))

/************************************************************************
 *									*
 *	rem_edg( ie ) - This routine removes an edge and eliminates	*
 * all the data associated with it. This does not maintain region	*
 * information.								*
 *									*
 * Original:	MEL	8/90						*
 *									*
 ************************************************************************/
rem_edg( ie )
int ie;
{
    int nr, i;
    int p1, p2, n1, n2;
    int newnd[20], l;
    double v1, v2;
    
    /*compute the magnitude of the point velocities*/
    v1 = mag( pt[(pt_edg(ie,0))]->vel );
    v2 = mag( pt[(pt_edg(ie,1))]->vel );

    /*pick the node moving slower*/
    if ( v1 > v2 ) nr = 1;
    else nr = 0;

    p1 = pt_edg(ie, nr);
    p2 = pt_edg(ie,(nr+1)%2);
    n1 = nd_edg(ie, nr);
    n2 = nd_edg(ie,(nr+1)%2);

    /*if either end is dead...*/
    if ( dead_pt(p1) || dead_pt(p2) ) return;
    if ( dead_nd(n1) || dead_nd(n2) ) return;

    /*figure out if joined node will be a triple*/
    for(l = i = 0; i < num_nd(p1); i++) {
	n1 = nd_pt(p1,i);
	if ( (n2 = node_mat( nd_pt(p2,0), mat_nd(n1))) == -1 ) {
	    newnd[l++] = n1;
	}
    }

    /*just remove the point, unless material goes away....*/
    if ( l + num_nd(p2) > 2 ) {
	for(i = 0; i < num_ele_edg(ie); i++)
	    if (dead_tri(ele_edg(ie,i))) return;

	for(i = 0; i < l; i++) {
	    add_nd_pt(p2, newnd[i]);
	    sub_nd_pt(newnd[i], p1);
	    pt_nd(newnd[i]) = p2;
	}

	rem_edg_end(ie, nr);
    }
    else {
	rem_pt(p1);
    }
}




rem_edg_end(ie, nr)
int ie, nr;
{
    int i, j, k, t, doit;
    int p1, p2, n1, n2, dl[20], l=0;

    n1 = nd_edg(ie,nr);
    n2 = nd_edg(ie,(nr+1)%2);
    p1 = pt_edg(ie,nr);
    p2 = pt_edg(ie,(nr+1)%2);
    if ( is_surf(p1) ) set_surf(p2);
    if ( is_back(p1) ) set_back(p2);

    /*update all neighbor triangles to have the new node*/
    for( i = 0; i < num_tri_nd(n1); i++ ) {
	t = tri_nd(n1,i);
	doit = TRUE;
	for(k = 0; k < num_ele_edg(ie); k++) 
	    if (ele_edg(ie,k) == t) doit = FALSE;

	/*add these triangles to the other node*/
	if ( doit ) {
	    for(j = 0; j < num_vert(t); j++) {
		if ( vert_tri(t,j) == n1 ) {
		    vert_tri(t,j) = n2;
		    add_tri_nd(n2, t);
		}
	    }

	    if ( vol_ele(t) == 0.0 ) {
		dl[l++] = t;
	    }
	    else {
		clr(tri[t], GEOMDN);
		clr(tri[t], NEIGHB);
		clr(tri[t], CLKWS);
		geom_dirty = clkws_dirty = neigh_dirty = TRUE;
	    }
	}
	else {
	    dl[l++] = t;
	}
    }

    if (l > 2) panic("too many triangles depend on edge");
    for(i = 0; i < l; i++) rem_ele(dl[i], FALSE);
    if ( !dead_nd(n1) ) rem_1nd(n1);

    edge_dirty = TRUE;
}



/* Experimental exact-domain preflight. An open fan may be removed without
 * changing its region boundary only if its two boundary rays are collinear.
 * Inspect every material at the physical point before deleting any node. */
extern int orient2d_float_sign(float *, float *, float *);
static int removal_preserves_boundary(int ip)
{
    int i, n, nt, j, k, t, v, other, region, boundaries, endpoint[2];
    int *ray, *count, nr, cap, ok;
    for (i = 0; i < num_nd(ip); i++) {
        n = nd_pt(ip,i);
        nt = num_tri_nd(n);
        if (nt == 0) continue;
        cap = 2 * nt;
        ray = malloc((size_t)cap * sizeof(*ray));
        count = malloc((size_t)cap * sizeof(*count));
        if (!ray || !count) {
            free(ray); free(count);
            panic("Out of memory in removal preflight");
            return FALSE;
        }
        nr = 0; region = -1; ok = TRUE;
        for (j = 0; j < nt && ok; j++) {
            t = tri_nd(n,j);
            if (dead_tri(t)) {ok = FALSE; break;}
            if (region < 0) region = reg_tri(t);
            if (reg_tri(t) != region) {ok = FALSE; break;}
            for (v = 0; v < num_vert(t); v++) {
                other = vert_tri(t,v);
                if (other == n) continue;
                if (pt_nd(other) == ip) {ok = FALSE; break;}
                for (k = 0; k < nr && ray[k] != other; k++);
                if (k == nr) {ray[nr] = other; count[nr++] = 1;}
                else count[k]++;
            }
        }
        boundaries = 0;
        for (k = 0; k < nr; k++) {
            if (count[k] == 1 && boundaries < 2)
                endpoint[boundaries++] = pt_nd(ray[k]);
            else if (count[k] != 2) ok = FALSE;
        }
        if (boundaries != 0 && boundaries != 2) ok = FALSE;
        if (boundaries == 2 &&
            orient2d_float_sign(pt[endpoint[0]]->cord, pt[ip]->cord,
                                pt[endpoint[1]]->cord) != 0) ok = FALSE;
        free(ray); free(count);
        if (!ok) return FALSE;
    }
    return TRUE;
}

rem_pt(ip)
int ip;
{
    int i;
    int loc[MAXMAT], maxn;

    if ( dead_pt(ip) ) return;
    if (mode == TWOD && getenv("SUPREM4_EXACT_DOMAIN_REMOVAL") &&
        !removal_preserves_boundary(ip)) {
        if (getenv("SUPREM4_REMOVAL_AUDIT"))
            fprintf(stderr, "[removal-audit] refused point=%d: exact domain not preserved\n", ip);
        clr(pt[ip], KILL_LATER);
        return;
    }

    if ( num_nd(ip) > 2 )
	panic("can't remove triple point");
    else {
	/*find the nodes at this point*/
	for(i = 0; i < num_nd(ip); i++) loc[i] = nd_pt(ip,i);
	maxn = num_nd(ip);
        if (getenv("SUPREM4_REMOVAL_TRACE_POINT") &&
            atoi(getenv("SUPREM4_REMOVAL_TRACE_POINT")) == ip) {
            fprintf(stderr, "[rem-pt-trace] point=%d surface=%d back=%d maxn=%d\n",
                    ip, is_surf(ip), is_back(ip), maxn);
            for(i=0;i<maxn;i++)
                fprintf(stderr, "[rem-pt-trace] step=%d node=%d material=%d triangles=%d\n",
                        i, loc[i], mat_nd(loc[i]), num_tri_nd(loc[i]));
        }

	/*for each node at this point*/
	for(i = 0; i < maxn; i++)  {
	    if (num_tri_nd(loc[i]) == 0 )
		rem_1nd(loc[i]);
	    else
		kill_node(loc[i]);
	}

    }
}



