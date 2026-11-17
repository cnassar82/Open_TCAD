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
/*   kill_node.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:20:57 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "constant.h"
#include "global.h"
#include "dbaccess.h"
#include "material.h"
#include "skel.h"
extern int orient2d_float_sign(float *, float *, float *);

static int kn_diag_on()
{
    char *v = (char *)getenv("SUPREM4_MESH_DIAG");
    return(v && *v && *v != '0');
}

/*-----------------KILL_NODE--------------------------------------------
 * Remove a node from the grid and fix up the triangulation.
 * Algorithm:
 *     if edge node, internalize the node to be removed
 *     call kill_internal_node to do the dirty work
 * Node is passed in as a triangle/vertex pair to avoid O(N) search every
 * time we do this.
 *----------------------------------------------------------------------*/
kill_node(in)
int in; 		/* Number of node to operate on */
{
    register int i, j;
    int maxn, *nb, tn, e, sr;
    float vec[MAXDIM], *ang, tf;
    int ireg,reg1,reg2;
    int ne_before;

    /*get all of the neighbors pulled together*/
    maxn = num_neigh_nd(in);
    /* Connectivity changes can produce valences beyond the old stack
     * buffers (100 neighbors, but only 50 angles). Size both alike. */
    nb = malloc((size_t)(maxn > 0 ? maxn : 1) * sizeof(*nb));
    ang = malloc((size_t)(maxn > 0 ? maxn : 1) * sizeof(*ang));
    if (!nb || !ang) {
	free(nb);
	free(ang);
	panic("Out of storage in kill_node");
	return;
    }
    for(i = 0; i < maxn; i++) nb[i] = neigh_nd(in,i);

    /*eliminate any neighbors that share the point*/
    for(i = 0; i < maxn; i++) {
	if ( pt_nd(nb[i]) == pt_nd(in) ) {
	    maxn--;
	    for(j = i; j < maxn; j++) nb[j] = nb[j+1];
	    i--;
	}
    }

    /*error checks on the number of neighbors*/
    if ( mode == ONED ) {
	if (maxn != 2) panic("can't do node removal in one dimension");
    }
    else if ( mode == TWOD ) {
	if ( maxn < 3 ) {
	    while ( num_tri_nd(in) > 0 ) rem_ele( tri_nd(in, 0), FALSE );
	    free(nb);
	    free(ang);
	    return;
	}

	/*compute the angular dependence*/
	for(i = 0; i < maxn; i++) {
	    vec[0] = cordinate(pt_nd(nb[i]),0) - cordinate(pt_nd(in),0);
	    vec[1] = cordinate(pt_nd(nb[i]),1) - cordinate(pt_nd(in),1);
	    ang[i] = atan2( vec[1], vec[0] );
	}

	/*bubble sort them into order*/
	for(i = 0; i < maxn; i++) {
	    for(j = i+1; j < maxn; j++) {
		if (ang[j] < ang[i]) {
		    tf = ang[i];
		    ang[i] = ang[j];
		    ang[j] = tf;
		    tn = nb[i];
		    nb[i] = nb[j];
		    nb[j] = tn;
		}
	    }
	}
    }

    /*build a skeleton from the list of nodes*/
    sr = cr_sreg(mat_nd(in));
    for(i = 0; i < maxn-1; i++) {
	e = mk_edg( nb[i], nb[i+1] );
	ad_edge(sr, e, sreg[sr]->bnd, MAYBE, BEFORE);
    }

    /*don't loop back around in one dimension*/
    if ( mode != ONED ) {
	e = mk_edg( nb[maxn-1], nb[0] );
	ad_edge(sr, e, sreg[sr]->bnd, MAYBE, BEFORE);
    }
    hinit(sr, MAXFLOAT);

    /*retriangulate the region*/
    ireg = reg_tri(tri_nd(in,0));
    if (pl_interface(in,&reg1,&reg2) > 1)
       if (reg1 > reg2) ireg = reg1;
       else ireg = reg2;

    if (kn_diag_on()) {
	fprintf(stderr, "[kn] kill in=%d pt=%d mat=%d ireg=%d ne=%d maxn=%d ntri=%d\n",
		in, pt_nd(in), mat_nd(in), ireg, ne, maxn, num_tri_nd(in));
	for (i = 0; i < maxn; i++)
	    fprintf(stderr, "[kn]   nb[%d]=%d pt=%d\n", i, nb[i], pt_nd(nb[i]));
    }

    if (getenv("SUPREM4_REMOVAL_TRACE_POINT") &&
        atoi(getenv("SUPREM4_REMOVAL_TRACE_POINT")) == pt_nd(in)) {
        int q, v, ti, node, point;
        if (getenv("SUPREM4_DIAG_LIVE_TRIANGLES"))
            ig2_write("removal_before.str", 0, 1.0);
        fprintf(stderr, "[removal-trace] point=%d node=%d material=%d region=%d ntri=%d ne=%d\n",
                pt_nd(in), in, mat_nd(in), ireg, num_tri_nd(in), ne);
        fprintf(stderr, "[removal-trace] center point=%d surface=%d back=%d nodes=", pt_nd(in),
                is_surf(pt_nd(in)), is_back(pt_nd(in)));
        for (q = 0; q < num_nd(pt_nd(in)); q++) {
            node = nd_pt(pt_nd(in),q);
            fprintf(stderr, "%s%d:%d:%d", q ? "," : "", node,
                    mat_nd(node), dead_nd(node));
        }
        fprintf(stderr, "\n");
        for (q = 0; q < maxn; q++) {
            point = pt_nd(nb[q]);
            fprintf(stderr, "[removal-trace] polygon[%d] node=%d point=%d x=%a y=%a\n",
                    q, nb[q], point, (double)xcord(point), (double)ycord(point));
            { int k;
              fprintf(stderr, "[removal-trace] owner point=%d surface=%d back=%d nodes=", point,
                      is_surf(point), is_back(point));
              for (k = 0; k < num_nd(point); k++) {
                  node = nd_pt(point,k);
                  fprintf(stderr, "%s%d:%d:%d", k ? "," : "", node,
                          mat_nd(node), dead_nd(node));
              }
              fprintf(stderr, "\n");
            }
        }
        for (q = 0; q < num_tri_nd(in); q++) {
            ti = tri_nd(in,q);
            fprintf(stderr, "[removal-trace] oldtri[%d] id=%d region=%d", q, ti, reg_tri(ti));
            for (v = 0; v < num_vert(ti); v++) {
                node = vert_tri(ti,v); point = pt_nd(node);
                fprintf(stderr, " v%d=node%d/point%d/%a,%a", v, node, point,
                        (double)xcord(point), (double)ycord(point));
            }
            fprintf(stderr, "\n");
        }
    }
    if (getenv("SUPREM4_REMOVAL_AUDIT") && mode == TWOD) {
        int a, b, pa, pb, pc, pd;
        for (a = 0; a < maxn; a++) {
            pa = pt_nd(nb[a]); pb = pt_nd(nb[(a+1)%maxn]);
            for (b = a+2; b < maxn; b++) {
                if (a == 0 && b == maxn-1) continue;
                pc = pt_nd(nb[b]); pd = pt_nd(nb[(b+1)%maxn]);
                if (orient2d_float_sign(pt[pa]->cord, pt[pb]->cord, pt[pc]->cord) *
                    orient2d_float_sign(pt[pa]->cord, pt[pb]->cord, pt[pd]->cord) < 0 &&
                    orient2d_float_sign(pt[pc]->cord, pt[pd]->cord, pt[pa]->cord) *
                    orient2d_float_sign(pt[pc]->cord, pt[pd]->cord, pt[pb]->cord) < 0)
                    fprintf(stderr, "[removal-audit] crossing proposal point=%d node=%d region=%d edges=(%d,%d),(%d,%d) ntri=%d\n",
                            pt_nd(in), in, ireg, pa, pb, pc, pd, num_tri_nd(in));
            }
        }
    }
    ne_before = ne;
    grid( sr, ireg, FALSE, FALSE );

    /*
     * If grid() created a triangle that duplicates a pre-existing one —
     * possible when the skeleton polygon overlaps an external triangle —
     * remove only the new copy.  The pre-existing triangle is ground truth.
     */
    if (mode == TWOD) {
	int t, t2, s0, s1, s2, q0, q1, q2, tmp, ndup = 0;
	for (t = ne_before; t < ne; t++) {
	    if (dead_tri(t)) continue;
	    if (reg_tri(t) != ireg) continue;
	    s0 = pt_nd(tri[t]->nd[0]); s1 = pt_nd(tri[t]->nd[1]); s2 = pt_nd(tri[t]->nd[2]);
	    if(s0>s1){tmp=s0;s0=s1;s1=tmp;} if(s1>s2){tmp=s1;s1=s2;s2=tmp;} if(s0>s1){tmp=s0;s0=s1;s1=tmp;}
	    for (t2 = 0; t2 < ne_before; t2++) {
		if (dead_tri(t2)) continue;
		q0 = pt_nd(tri[t2]->nd[0]); q1 = pt_nd(tri[t2]->nd[1]); q2 = pt_nd(tri[t2]->nd[2]);
		if(q0>q1){tmp=q0;q0=q1;q1=tmp;} if(q1>q2){tmp=q1;q1=q2;q2=tmp;} if(q0>q1){tmp=q0;q0=q1;q1=tmp;}
		if (s0==q0 && s1==q1 && s2==q2) {
		    if (kn_diag_on())
			fprintf(stderr, "[kn] dup new=t%d old=t%d spts=(%d,%d,%d) removing new\n",
				t, t2, s0, s1, s2);
		    rem_ele(t, FALSE);
		    ndup++;
		    break;
		}
	    }
	}
	if (ndup && kn_diag_on())
	    fprintf(stderr, "[kn] removed %d dup(s) in=%d ireg=%d\n", ndup, in, ireg);
    }

    /*
     * If a degenerate skeleton boundary has coincident physical points,
     * grid() can create two new triangles that both claim the same physical
     * boundary edge.  Later neighbor construction can then produce
     * asymmetric triangle neighbors and stale edge slots.  Keep the first
     * new triangle on each skeleton boundary edge and remove later claimants.
     */
    if (mode == TWOD) {
	int bpa, bpb, f, cnt, ndup = 0;
	int t, p0, p1;

	for (i = 0; i < maxn; i++) {
	    bpa = pt_nd(nb[i]);
	    bpb = pt_nd(nb[(i+1) % maxn]);
	    cnt = 0;

	    for (t = ne_before; t < ne; t++) {
		if (dead_tri(t)) continue;
		if (reg_tri(t) != ireg) continue;

		for (f = 0; f < 3; f++) {
		    p0 = pt_nd(nd_face_ele(t, f, 0));
		    p1 = pt_nd(nd_face_ele(t, f, 1));
		    if (!((p0 == bpa && p1 == bpb) ||
			  (p0 == bpb && p1 == bpa))) continue;

		    cnt++;
		    if (cnt > 1) {
			if (kn_diag_on())
			    fprintf(stderr,
				    "[kn] bedge-dup new=t%d pts=(%d,%d) removing\n",
				    t, bpa, bpb);
			rem_ele(t, FALSE);
			ndup++;
		    }
		    break;
		}
	    }
	}
	if (ndup && kn_diag_on())
	    fprintf(stderr, "[kn] removed %d boundary dup(s) in=%d ireg=%d\n",
		    ndup, in, ireg);
    }

    /*kill triangles until none are left*/
    while ( num_tri_nd(in) > 0 ) rem_ele( tri_nd(in, 0), FALSE );
    free_skel(sr);
    free(nb);
    free(ang);
}

	

/*-----------------TROTATE----------------------------------------------
 * Tricky little routine to walk around a boundary edge looking
 * for the next node.
 *----------------------------------------------------------------------*/
trotate( ie, iv, clock, Pje, Poje, Pj, Pjj)
     int ie;			/* The triangle */
     int iv;			/* The index of the boundary node */
     int clock;			/* 1 for clockwise, 0 for cc */
     int *Pje;			/* The outside neighbor (loop-breaker) */
     int *Poje;			/* The inside neighbor */
     int *Pj;			/* index of rotation node in oje */
     int *Pjj;			/* index of leading node or trailing triangle */
{
    int j, jj, je, oje, delta, ir=tri[ie]->regnum, in=tri[ie]->nd[iv];
    
    /*
     * Rotate around given node (in) until we
     *  a: come back to where we started    -or-
     *  b: hit an interface.
     *  Loop invariants: tri[oje]->nd[j] is the pivot (in)
     *  		 tri[oje]->nb[jj] precedes oje.
     *                   je succeeds oje
     *  The neighbors must be correct for this to work.
     */

    delta = clock? 1 : 2;
    for (oje=ie, j = iv, jj = (j+delta)%3; ; ) {
	 
	 je = tri[ oje]->nb[ 3-j-jj];
	 if ((je < 0) || (je == ie) || (tri[ je]->regnum != ir)) break;

	 jj = twhich (je, oje);
	 j  = nwhich (je, in);
	 oje = je;
     }

    *Pjj = jj;
    *Pj = j;
    *Poje = oje;
    *Pje = je;
}
          
