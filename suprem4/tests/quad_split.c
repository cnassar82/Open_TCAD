/* Regression for the concave oxidation cavity that produced an exterior
 * diagonal. Exercise the production splitter, including both windings. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "dbaccess.h"
#include "skel.h"

extern struct LLedge *quadsplit();

void panic(char *message)
{
    fprintf(stderr, "Unexpected panic: %s\n", message);
    abort();
}

static double area2(int a, int b, int c)
{
    double ax = pt[a]->cord[0], ay = pt[a]->cord[1];
    double bx = pt[b]->cord[0], by = pt[b]->cord[1];
    double cx = pt[c]->cord[0], cy = pt[c]->cord[1];
    return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
}

static void check(const float coordinates[4][2], int valid)
{
    struct pt_str points[4], *point_ptrs[4];
    struct nd_str nodes[4], *node_ptrs[4];
    struct edg_str edges[4], *edge_ptrs[4];
    struct LLedge links[4], *split;
    struct sreg_str region;
    int rotation, reverse, i, a, b, c, d;
    double polygon, triangles, side1, side2;
    memset(&region, 0, sizeof(region));
    memset(points, 0, sizeof(points));
    memset(nodes, 0, sizeof(nodes));
    memset(edges, 0, sizeof(edges));
    memset(links, 0, sizeof(links));
    pt = point_ptrs; nd = node_ptrs; edg = edge_ptrs;
    mode = TWOD;
    region.len = 4;
    for (reverse = 0; reverse < 2; reverse++) {
        for (i = 0; i < 4; i++) {
            int source = reverse ? 3-i : i;
            pt[i] = &points[i]; nd[i] = &nodes[i]; edg[i] = &edges[i];
            points[i].cord[0] = coordinates[source][0];
            points[i].cord[1] = coordinates[source][1];
            nodes[i].pt = i;
            edges[i].nd[0] = i; edges[i].nd[1] = (i+1)%4;
            links[i].edge = i; links[i].iscc = TRUE;
            links[i].prev = &links[(i+3)%4];
            links[i].next = &links[(i+1)%4];
        }
        for (rotation = 0; rotation < 4; rotation++) {
            region.bnd = &links[rotation];
            split = quadsplit(&region);
            if (!valid) {
                assert(split == NULL);
                assert(!skel_split_valid(&links[0], &links[2]));
                assert(!skel_split_valid(&links[1], &links[3]));
                continue;
            }
            assert(split != NULL);
            a = nB(split->prev); b = nB(split); c = nB(split->next);
            d = nB(split->next->next);
            side1 = area2(a,c,b); side2 = area2(a,c,d);
            /* An interior diagonal separates the two remaining vertices. */
            assert((side1 < 0 && side2 > 0) || (side1 > 0 && side2 < 0));
            polygon = fabs(area2(0,1,2) + area2(0,2,3));
            assert(skel_split_valid(split->prev, split->next) ==
                   (area2(0,1,2) + area2(0,2,3) > 0));
            triangles = fabs(side1) + fabs(side2);
            assert(fabs(triangles-polygon) <= polygon * 1e-12);
        }
    }
}

int main(void)
{
    const float oxidation[4][2] = {
        {0.012950016185641289, 4.3748914322350174e-5},
        {0.0129757, 3.59436e-5},
        {0.012979512102901936, 4.4311542296782136e-5},
        {0.013000000268220901, 4.7856512537691742e-5}
    };
    const float convex[4][2] = {{0,0}, {3,0}, {2,2}, {0,1}};
    const float concave[4][2] = {{0,0}, {4,0}, {1,1}, {0,4}};
    const float collinear_corner[4][2] = {{0,0}, {1,0}, {2,0}, {0,2}};
    const float crossed[4][2] = {{0,0}, {2,2}, {0,2}, {2,0}};
    const float shifted_mask_crossed[4][2] = {
        {.012146644294261932,2.980830322485417e-5},
        {.012171067297458649,3.229017602279782e-5},
        {.012163200415670872,3.149074109387584e-5},
        {.012155062519013882,3.0663766665384173e-5}
    };
    check(oxidation, 1);
    check(convex, 1);
    check(concave, 1);
    check(collinear_corner, 1);
    check(crossed, 0);
    check(shifted_mask_crossed, 0);
    puts("PASS: internal quad diagonals and area conservation for both windings");
    return 0;
}
