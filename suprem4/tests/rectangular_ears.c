/* A right-angle ear can contain all remaining vertices on its diagonal.
 * The former rectangle shortcut clipped it and left a flat quadrilateral. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "skel.h"
#include "material.h"
#include "impurity.h"
#include "diffuse.h"
#include "regrid.h"

void panic(char *message) { fprintf(stderr, "%s\n", message); abort(); }
int edge_pl(int edge) { assert(!"debug plotting is disabled"); return 0; }

static void check(int rotation)
{
    float coordinates[5][3] = {{0,0,0},{4,0,0},{3,1,0},{2,2,0},{0,4,0}};
    int nodes[5], boundary[5], seen[5] = {0}, r, i, j, t, start = ne;
    double total = 0;
    r = cr_sreg(SiO2);
    for (i = 0; i < 5; i++) nodes[i] = mk_nd(mk_pt(2, coordinates[(i+rotation)%5]), SiO2);
    for (i = 0; i < 5; i++) {
        boundary[i] = mk_edg(nodes[i], nodes[(i+1)%5]);
        ad_edge(r, boundary[i], sreg[r]->bnd, MAYBE, BEFORE);
    }
    ck_clock(r, TRUE);
    rect_tri(r, 0);
    free_skel(r);
    assert(ne-start == 3);
    for (t = start; t < ne; t++) {
        float *a = pt[pt_nd(tri[t]->nd[0])]->cord;
        float *b = pt[pt_nd(tri[t]->nd[1])]->cord;
        float *c = pt[pt_nd(tri[t]->nd[2])]->cord;
        double area = ((double)b[0]-a[0])*((double)c[1]-a[1]) -
                      ((double)b[1]-a[1])*((double)c[0]-a[0]);
        assert(area > 0);
        total += area / 2;
        for (i = 0; i < 3; i++) for (j = 0; j < 5; j++)
            if (tri[t]->nd[i] == nodes[j]) seen[j]++;
    }
    assert(fabs(total-8) < 1e-12);
    for (i = 0; i < 5; i++) {
        assert(seen[i] > 0);
        assert(edg[boundary[i]]->ele.num == 1);
    }
}

int main(void)
{
    int i;
    mode = TWOD; nvrt = nedg = 3;
    for (i = 0; i < 5; i++) check(i);
    for (i = 0; i < ne; i++) free(tri[i]);
    for (i = 0; i < ned; i++) { dis_list(&edg[i]->ele); dis_list(&edg[i]->skel); free(edg[i]); }
    for (i = 0; i < nn; i++) { dis_list(&nd[i]->edg); dis_list(&nd[i]->ltri); free(nd[i]); }
    for (i = 0; i < np; i++) free(pt[i]);
    free(tri); free(edg); free(nd); free(pt);
    puts("PASS: rectangular ears preserve all boundary vertices and positive area");
    return 0;
}
