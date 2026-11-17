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
/*   list.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:20:59 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <global.h>
#include <constant.h>
#include <dbaccess.h>

/*
 *	Add to a list.
 */
add_list( l, item )
struct list_str *l;
int item;
{
    int capacity;
    int *items;

    /*check for space consideration*/
    if ( l->num >= l->all ) {
	if (l->num == INT_MAX) {
	    panic("Connectivity list exceeds integer capacity");
	    return;
	}
	capacity = l->all == 0 ? 10 :
	    (l->all > INT_MAX / 2 ? INT_MAX : 2 * l->all);
	if ((size_t)capacity > (size_t)-1 / sizeof(*items)) {
	    panic("Connectivity list allocation size overflow");
	    return;
	}
	items = realloc(l->list, (size_t)capacity * sizeof(*items));
	if (!items) {
	    panic("Out of storage growing connectivity list");
	    return;
	}
	l->list = items;
	l->all = capacity;
    }

    /*store the triangle number away*/
    l->list[ l->num++] = item;
}


dis_list(l)
struct list_str *l;
{
    if (l->list) free(l->list);
    l->list = NULL;
    l->num = 0;
    l->all = 0;
}

pr_list(l)
struct list_str *l;
{
    int i;

    printf("    ");
    for(i = 0; i < l->num; i++) printf("%d ", l->list[i]);
    printf("\n");
}



sub_list(l, item)
struct list_str *l;
int item;
{
    int i, j;

    for(i = 0; (i < l->num) && (l->list[i] != item); i++);

    if ( i == l->num )
	panic("subtracting nonexistant value from list");

    for(j = i; j < l->num-1; j++) l->list[j] = l->list[j+1];

    l->num--;
}

sub_nd_pt(n, p)
int n, p;
{
    int i, j;

    for(i = 0; (i < pt[p]->nn) && (pt[p]->nd[i] != n); i++);

    if ( i == pt[p]->nn )
	panic("subtracting nonexistant node from point");

    for(j = i; j < pt[p]->nn; j++) pt[p]->nd[j] = pt[p]->nd[j+1];

    pt[p]->nn--;
}

