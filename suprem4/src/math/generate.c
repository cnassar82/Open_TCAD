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
/*   generate.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:44:13 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "global.h"

/* Grow before writing; keep all array offsets representable by int. */
static int reserve_map(int **map, int *capacity, int needed)
{
    int next;
    int *grown;
    if (needed <= *capacity) return(0);
    next = *capacity > INT_MAX / 2 ? INT_MAX : 2 * *capacity;
    if (next < needed) next = needed;
    grown = sralloc(int, next, *map);
    if (!grown) return(-1);
    *map = grown;
    *capacity = next;
    return(0);
}

/************************************************************************
 *									*
 *	generate - This routine generates the ia, ja list		*
 *									*
 ************************************************************************/

generate( nn, symb_var, symm, get_connect, pia, aoff)
    int nn;			/* Number of nodes */
    int symb_var;		/* Numbers of degrees of freedom per node */
    int symm;			/* Do upper *and* lower halves? */
    int (*get_connect)(int, int *, int *); /* Returns full-width node IDs */
    int **pia;			/* Symbolic matrix map */
    int *aoff;			/* Size of upper triangle (redundant) */
{
    int i, j;			/*indices*/
    int cnt;			/*first empty space in ia*/
    int *neigh, *grown;		/*neighbors of the current point*/
    int len;
    int num;			/*#of entries in neigh*/
    int ih, il, row, col, iamax = *aoff, *ia = *pia;


    len = 100;
    neigh = salloc(int, len);
    if (!neigh) return(-1);

    /*symbolic entries for each node begin just beyond the head pointers */
    /*or to put it another way, ja[0] == ia[neq]+1*/
    if (nn < 0 || symb_var <= 0 || nn > (INT_MAX - 1) / symb_var)
        panic("Sparse matrix dimension exceeds integer capacity");
    if (reserve_map(pia, &iamax, symb_var * nn + 1) < 0) {
        free(neigh);
        return(-1);
    }
    ia = *pia;
    ia[ 0]= symb_var * nn + 1;
    cnt = ia[ 0];

    /*foreach node in the mesh*/
    for(i = 0; i < nn; i++) {

	/*call the connection routine for a list of the node's neighbors */
	num = len;
	while ( ! (0== (*get_connect)( i, neigh, &num)) ) {
	    if (len > INT_MAX / 2)
                panic("Sparse neighbor count exceeds integer capacity");
	    len *= 2;
	    num = len;
	    grown = sralloc(int, len, neigh);
            if (!grown) { free(neigh); return(-1); }
            neigh = grown;
	}
        if (num < 0 || num > len) panic("Invalid sparse neighbor count");
        for (j = 0; j < num; j++)
            if (neigh[j] < 0 || neigh[j] >= nn)
                panic("Sparse neighbor is outside the mesh");

	/*put this fine group into the ia, ja list*/
	/*remember each node may have multiple variables*/
	for( ih = 0; ih < symb_var; ih++) {
	    row = symb_var * i + ih;

	    for( j = 0; j < num; j++ )

		for( il = 0; il < symb_var; il++) {
		    col = symb_var * neigh[ j] + il;

		    /*Symmetry: just keep the upper triangle */
		    if (symm || col > row) {
                        if (cnt == INT_MAX)
                            panic("Sparse matrix entries exceed integer capacity");
                        if (reserve_map(pia, &iamax, cnt + 1) < 0) {
                            free(neigh);
                            return(-1);
                        }
                        ia = *pia;
			ia[ cnt++ ] = col;
		    }

		}
	    ia[row + 1] = cnt;
	}
    }

    *aoff = ia[symb_var*nn] - ia[ 0];
    free(neigh);
    return(0);

#ifdef DEBUG
    check_db();
#endif

}

/************************************************************************
 *									*
 *	iaja()	- Find the matrix location given i,j			*
 *									*
 ************************************************************************/
iaja( ia, aoff, row, col)
    int *ia, aoff;
    int row;
    int col;
{
    int i, off;

    /* Diagonal comes first */
    if (row == col)
	return (row);

    /* Off-diagonal entries: only store half the connectivity*/
    if (row < col)
	off = aoff;
    else {
	i = col; col = row; row = i;
	off = 0;
    }

    /* Your basic linear search of the elements in this row */
    for (i = ia[ row]; i < ia[ row+1]; i++)
	if (ia[ i] == col) return( i + off);

    /* Trouble right here in river city */
    return(-32767);
}

#ifdef DEBUG
/************************************************************************
 *									*
 *	debug routines							*
 *	useful for looking at variables dbx won't let us see (snarl)	*
 *									*
 ************************************************************************/
#include <stdio.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "matrix.h"

dump_ij(ia)
     int *ia;
{
    int i, j;
    for (i = 0; i < ia[ 0]; i++)
	printf ("ia[ %d] =  %d\n", i, ia[i]);

    for (j = 0; j < ia[ 0]-1; j++) {
	printf("row %d of ja: ", j);
	for (i = ia[j]; i < ia[j+1]; i++)
	    printf("%d ", ia[i]);
	printf("\n");
    }
}

dump_a( ia, a, rhs)
     int *ia;
     double *a;
     double *rhs;
{
    int i, j;
    float yo;
    for (i = 0; i < ia[ 0]; i++)
	printf ("a[ %d] =  %g \twrk[ %d] = %g\n", i, a[i], i, rhs[i]);

    for (j = 0; j < ia[ 0]-1; j++) {
	printf( "col %d of a: ", j);
	for (i = ia[ j]; i < ia[ j+1]; i++) {
	    yo = a[ i];
	    if (yo)
		printf("a[%d,%d] = %9.2g ", ia[i], j, yo);
	    else
		printf("a[%d,%d] ", ia[i], j);
	}
	printf("\n");
	printf( "row %d of a: ", j);
	for (i = ia[ j]; i < ia[ j+1]; i++) {
	    yo = a[ i+aoff];
	    if (yo)
		printf("a[%d,%d] = %9.2g ", j, ia[i], yo);
	    else
		printf("a[%d,%d] ", j, ia[i]);
	}
	printf("\n");
    }
}

#endif
