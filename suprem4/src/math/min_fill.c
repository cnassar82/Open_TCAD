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
/*   min_fill.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:44:15 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#ifdef SUPREM_USE_AMD
#include <limits.h>
#include <suitesparse/amd.h>
#endif
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
extern o_compar();

/*#define salloc(S,N) (S *) malloc( (N)*sizeof(S))
    /*printf("allocated N S's = %d bytes\n", (N)*sizeof(S))*/

/*#define sralloc(Old,Type,Num) (Type *) realloc( Old, Num * sizeof(Type))
    /*printf("reallocated Num Type's = %d bytes\n", (Num)*sizeof(Type))*/

struct d_sort {
       float x;
       int  nbin, mat;
       } ;

/* Blocks never move: active adjacency lists contain pointers into them. */
typedef struct _link {struct _link *n; int v;} link;
typedef struct link_block {
    struct link_block *next;
    link *items;
} link_block;
static link_block *storage;
static link *freelist;
static size_t storage_size;
static size_t su, hwm;		/* Storage used, high water mark */
static void new_storage(void);
static void free_storage(void);

#define takefrom( F, O) {link *t = (O); (O) = (F); (F) = (F)->n; (O)->n = t;}

#define new_link(A) {if (!freelist) new_storage(); takefrom( freelist, A);\
		     if (++su > hwm) hwm = su;}
#define free_link(A) {takefrom((A)->n, freelist); su--;}

#ifdef SUPREM_USE_AMD
static int amd_fill(int *ia, int *reorder, int talk)
{
    int n = ia[0] - 1, base = ia[0], i, status;
    int *columns, *permutation;
    double info[AMD_INFO];
    columns = salloc(int, (size_t)n + 1);
    permutation = salloc(int, n);
    if (!columns || !permutation) panic("Out of memory");
    for (i = 0; i <= n; i++) columns[i] = ia[i] - base;
    /* The adjacency map is CSR. AMD orders the pattern of A+A', so its
     * transpose (CSC interpretation) gives the same symmetric graph. */
    status = amd_order(n, columns, ia + base, permutation, NULL, info);
    free(columns);
    if (status != AMD_OK && status != AMD_OK_BUT_JUMBLED)
        panic("AMD sparse ordering failed");
    for (i = 0; i < n; i++) reorder[i] = -1;
    for (i = 0; i < n; i++) {
        int old = permutation[i];
        if (old < 0 || old >= n || reorder[old] != -1)
            panic("AMD returned an invalid permutation");
        /* SUPREM wants old->new; AMD returns new->old. */
        reorder[old] = i;
    }
    free(permutation);
    if (!isfinite(info[AMD_LNZ]) || info[AMD_LNZ] < 0 || info[AMD_LNZ] > INT_MAX)
        panic("AMD fill estimate exceeds integer capacity");
    if (talk) printf("AMD ordering: nodes=%d, estimated strict L entries=%.0f\n",
                     n, info[AMD_LNZ]);
    /* FEsymb uses this only as an initial symbolic storage estimate. */
    return (int)ceil(info[AMD_LNZ]);
}
#endif


/************************************************************************
 *									*
 *	min_ia_fill( ia, reorder, talk) - Find a minimum fill ordering. *
 *									*
 *  Original:	MEL	11/84						*
 *  Modified:   CSR     09/86 No inkludge files, thanks			*
 *  Rewritten:	CSR	02/86 Don't give up until malloc does		*
 *	        LCC	06/90 Add 1-D mode
 *									*
 * the idea here is to do a mock gauss elimantion, updating the pt to pt*
 * list at each step.  At each step, a decision is made about the	*
 * optimal node to work on next. Returns size of strict upper triangle.	*
 ************************************************************************/
min_ia_fill(ia, reorder, talk)
     int *ia;			/* Matrix map including l and u */
     int *reorder;		/* Result vector */
     int talk;
{
    link *nbrs;			/* Base of each nodes nbr list: .v=length */
    int in, nn, j;

    struct d_sort *order;
    int count=0;
    float tmp;

    int i, best, besti, bestn;	/* Best is node being eliminated */
    int *remaining, nremain;	/* Array of uneliminated nodes */
    int *mask;			/* Tmp for computing union */
    int work;			/* neighbor of best being worked on */
    link *p, *workl, *bestl;
    register link *q, *lastq;
    int  nnb, step, Obestn = -1, total = 0, e;

#ifdef SUPREM_USE_AMD
    if (mode == TWOD) return amd_fill(ia, reorder, talk);
#endif

    /* Initialize storage - make a guess of l ~ 5*a, update later.*/
    nn = ia[0] - 1;

    order = salloc( struct d_sort, nn );

    if ( mode == TWOD ) {
        mask = salloc( int, nn);
        remaining = salloc( int, nn);
    }

    nbrs = salloc( link, nn);
    storage_size = suprem_array_bytes((size_t)ia[nn], 10);
    su = hwm = 0;
    storage = 0;
    freelist = 0;
    new_storage();

    /* Initialize link table */
    for (in = 0; in < nn; in++) {

	/* 1-D we store data in ORDER and qsort it ! */
	if ( mode == ONED ) {
	    order[in].x = pt[ nd[in]->pt ]->cord[0];
	    order[in].nbin = in;
	    order[in].mat  = nd[in]->mater;
	    count++;
        }
        else {
       	    mask[ in] = 0;
            remaining[ in] = in;
	}

	nbrs[ in].v = ia[ in+1] - ia[ in];
	for (nbrs[ in].n = 0, j = ia[ in]; j < ia[in+1]; j++) {
	    new_link( nbrs[in].n);
	    nbrs[ in].n->v = ia[j];
	}
    }


    /*Switch on the mode*/
    switch ( mode ) {
    case  ONED :

	/*do a linear search to find an exposed surface node*/
	for(best = -1, i = 0; i<ne && best < 0; i++) {
	    if (tri[i]->nb[0] == EXPOSED) best = tri[i]->nd[1];
	    if (tri[i]->nb[1] == EXPOSED) best = tri[i]->nd[0];
	}

	/*get the gas node!*/
	best = node_mat(best, GAS);
	step = 0;
	besti = -1;
	reorder[best] = step++;

	while(step < nn) {
	    /*if this is a multiple node point, go to the other one*/
	    if ( num_nd(pt_nd(best)) == 2 ) {
		besti = best;
		if ( nd_pt(pt_nd(best),0) != best ) best = nd_pt(pt_nd(best),0);
		else best = nd_pt(pt_nd(best),1);
		reorder[best] = step++;
	    }

	    /*walk the edges at this node*/
	    for(work = -1, i = 0; work == -1 && i < num_edge_nd(best); i++) {
		e = edge_nd(best, i);
		if ( (nd_edg(e,0) == best) ) {
		    if ( nd_edg(e,1) != besti ) {
			work = besti = best;
			best = nd_edg(e,1);
			reorder[best] = step++;
		    }
		}
		else {
		    if ( nd_edg(e,0) != besti ) {
			work = besti = best;
			best = nd_edg(e,0);
			reorder[best] = step++;
		    }
		}
	    }
	}

        break;


    case TWOD :

    for (nremain = nn, step = 0; step < nn; step++) {

	/* Find the node most wanting elimination */
	for (bestn = nn+1, i = 0; i < nremain; i++) {
	    if ((nnb = nbrs[ remaining[ i]].v) < bestn) {
		bestn = nnb;
		besti = i;
		if (bestn <= Obestn) break;
	    }
	}
	best = remaining[ besti];
	Obestn = bestn;
	bestl = &nbrs[ best];

	/* Take it off the remaining list and add it to the reorder list*/
	remaining[ besti] = remaining[ --nremain];
	reorder[ best] = step;

	/* Eliminate node best: take best out of each of its neighbors
	   adjacency lists, and union its adjacency list in instead */
	for (p = bestl->n; p; p = p->n) {
	    work = p->v;
	    if (work == best) continue;
	    workl = &nbrs[ work];

	    /* To speed union, generate membership mask */
	    for (q = workl->n; q; q = q->n) mask[q->v] = 1;

	    /* Take best out*/
	    for (lastq = workl, q = lastq->n; q; lastq = q, q=q->n)
		if (q->v == best)
		    { free_link( lastq); break; }
	    workl->v--;

	    /* Form union by stuffing new elements in front */
	    for (q = bestl->n; q; q = q->n)
		if (!mask[ q->v]) {
		    new_link( workl->n);
		    workl->n->v = q->v;
		    workl->v++;
		}

	    /* Reset mask */
	    mask[best] = 0;
	    for (q = workl->n; q; q = q->n) mask[ q->v] = 0;

	    /* Might have a new lower count [never happens?] */
	    if (workl->v < Obestn) Obestn = workl->v;
	}

	/* Mark best as dead and gone */
	while( bestl->n)
	    { free_link( bestl); total++;}
    }

    /* Self protection */
    for (in = 0; in < nn; in++) mask[ reorder[ in]]++;
    for (in = 0; in < nn; in++) assert(mask[ in] == 1);
    free(remaining);
    free(mask);
    break;
   }

    if (talk)
	 printf("high water mark=%zu [%f *asize]  loff=%d\n",
		hwm, (float) hwm / ia[ nn], total);

    free(nbrs);
    free_storage();
    free(order);
    return( total - nn);
}

o_compar(f1, f2)
struct d_sort  *f1, *f2;
{
    if ( f1->x > f2->x )
	return( 1 );
    else if ( f1->x < f2->x )
	return( -1 );
    else if ( f1->mat > f2->mat )
	return( 1 );
    else if ( f1->mat < f2->mat )
	return( -1 );
    else
	return( 0 );
}

static void new_storage(void)
{
    link *start, *end, *p;
    link_block *block;

    assert(!freelist && storage_size > 0);
    block = salloc(link_block, 1);
    if (!block) panic("Out of memory");
    block->items = salloc(link, storage_size);
    if (!block->items) {
        free(block);
        panic("Out of memory");
    }
    block->next = storage;
    storage = block;
    start = block->items;
    end = start + storage_size - 1;

    for (p = start; p < end; p++) p->n = p+1;
    end->n = 0;
    freelist = start;
}

static void free_storage(void)
{
    link_block *block;
    while (storage) {
        block = storage;
        storage = block->next;
        free(block->items);
        free(block);
    }
    freelist = 0;
    su = 0;
}

#ifdef notdef
run_length(nn,nbrs)
     int nn;
     link *nbrs;
{
    int in, last, nrlen, rlen, sumrlen;
    float averlen;
    link *j;

    sumrlen = nrlen = 0;
    for (in = 0; in < nn; in++) {
	rlen = 1;
	for (j = nbrs[in].n->n; last = j->v, j->n; j = j->n) {
	    if (j->n->v == last-1)
		rlen++;
	    else {
		sumrlen += rlen;
		nrlen ++;
		rlen = 1;
	    }
	}
    }
    averlen = sumrlen / nrlen;
    printf("average run length: %f\n", averlen);
}


#endif
