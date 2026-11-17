/* Exercise new_layer's edge-selection rule at an oxidation mask corner. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "material.h"

static int splits;
int sp_edge(int edge, float *coordinate, int *point, int flag)
{
    assert(edge==0 && flag==FALSE);
    assert(fabs(coordinate[0]-5e-8)<1e-13);
    assert(coordinate[1]==0);
    splits++;
    *point=0;
    return 0;
}

int main(void)
{
    struct nd_str nodes[3] = {0};
    struct pt_str points[2] = {0};
    struct edg_str edge = {0};
    struct nd_str *node_ptrs[3] = {&nodes[0],&nodes[1],&nodes[2]};
    struct pt_str *point_ptrs[2] = {&points[0],&points[1]};
    struct edg_str *edge_ptrs[1] = {&edge};
    int adjacent[1] = {0};
    nd=node_ptrs; pt=point_ptrs; edg=edge_ptrs;
    mode=TWOD; nn=3; np=2; ned=1;
    nodes[0].mater=nodes[2].mater=SiO2; nodes[1].mater=Si;
    nodes[0].pt=nodes[1].pt=0; nodes[2].pt=1;
    nodes[0].edg.list=adjacent; nodes[0].edg.num=1;
    points[0].nn=2; points[0].nd[0]=0; points[0].nd[1]=1;
    points[1].nn=1; points[1].nd[0]=2;
    points[1].cord[0]=1e-4;
    edge.nd[0]=0; edge.nd[1]=2; edge.len=1e-4;
    /* Neighbor is ahead of the moving interface: no eligible edge. */
    points[0].vel[0]=1e-7;
    assert(new_layer(Si,SiO2,1e-5,1000.0)==FALSE && splits==0);
    /* Tangential and stationary motion must also leave it alone. */
    points[0].vel[0]=0; points[0].vel[1]=1e-7;
    assert(new_layer(Si,SiO2,1e-5,1000.0)==FALSE && splits==0);
    points[0].vel[1]=0;
    assert(new_layer(Si,SiO2,1e-5,1000.0)==FALSE && splits==0);
    /* Opposing edge, with sufficient spacing and motion: do split. */
    points[0].vel[0]=-1e-7;
    assert(new_layer(Si,SiO2,1e-5,1000.0)==TRUE && splits==1);
    puts("PASS: oxidation layer insertion handles mask-corner edge selection");
    return 0;
}
