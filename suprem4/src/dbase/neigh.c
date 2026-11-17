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
/*   neigh.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:21:04 */

#include "global.h"
#include "constant.h"
#include "dbaccess.h"

static struct list_str node_neighbors;
static int work_nd = -1;
static struct list_str point_neighbors;
static int work_pt = -1;


/*
 *	num_neigh_nd() - This routine computes a nodes neighbors,
 * and stores them up for fast access.  It returns the number of
 * neighbors.
 */
num_neigh_nd( n )
int n;
{
    int i, j, p, t;

    node_neighbors.num = 0;
    work_nd = n;

    /*there is always coupling between the different solutions at a node*/
    p = pt_nd(n);
    for(i = 0; i < num_nd(p); i++) {
	add_list(&node_neighbors, nd_pt(p,i));
    }

    /*step through all the triangles and mark their spots true*/
    for(i = 0; i < num_tri_nd(n); i++) {
	t = tri_nd(n, i);
	for(j = 0; j < num_vert(t); j++) {
	    add_list(&node_neighbors, vert_tri(t,j));
	}
    }

    elim_dups(node_neighbors.list, &node_neighbors.num);

    return(node_neighbors.num);
}



/*
 *	neigh_nd() - This routine computes a nodes neighbors,
 * and stores them up for fast access.  It returns the number of
 * neighbors.
 */
int neigh_nd( n, i )
int n, i;
{
    if ( n != work_nd ) (void)num_neigh_nd(n);
    return( node_neighbors.list[i] );
}


/*
 *	num_neigh_pt() - This routine computes a nodes neighbors,
 * and stores them up for fast access.  It returns the number of
 * neighbors.
 */
num_neigh_pt( p )
int p;
{
    int i, j, k, n, t;

    point_neighbors.num = 0;
    work_pt = p;

    /*there is always coupling between the different solutions at a node*/
    for(i = 0; i < num_nd(p); i++) {
	n = nd_pt(p,i);
	add_list(&point_neighbors, pt_nd(n));

	/*step through all the triangles and mark their spots true*/
	for(k = 0; k < num_tri_nd(n); k++) {
	    t = tri_nd(n, k);
	    for(j = 0; j < num_vert(t); j++) {
		add_list(&point_neighbors, pt_nd(vert_tri(t,j)));
	    }
	}
    }

    elim_dups(point_neighbors.list, &point_neighbors.num);

    return(point_neighbors.num);
}


/*
 *	neigh_nd() - This routine computes a nodes neighbors,
 * and stores them up for fast access.  It returns the number of
 * neighbors.
 */
int neigh_pt( p, i )
int p, i;
{
    if ( p != work_pt ) (void)num_neigh_pt(p);
    return( point_neighbors.list[i] );
}


/************************************************************************
 *									*
 *	clean_list( taken, num) - This routine sorts and removes 	*
 *  duplicate entries from a list of nodes.				*
 *									*
 ************************************************************************/

elim_dups(taken, num)
int *taken;
int *num;
{
    int i, j, swap, ndup;

    ndup = *num;
    if (ndup == 0) return;

    /*sort the neighbor list (bubble sort since lists are short)*/
    for(i = 0; i < ndup; i++) {
	for(j = i + 1; j < ndup; j++) {
	    if ( taken[i] > taken[j] ) {
		swap = taken[i];
		taken[i] = taken[j];
		taken[j] = swap;
	    }
	}
    }

    /*eliminate duplicates*/
    *num=1;
    for(i = 1; i < ndup; i++)
	if (taken[i] != taken[*num -1])
	    taken[ (*num)++] = taken[ i];

}
