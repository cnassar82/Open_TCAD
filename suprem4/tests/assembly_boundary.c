/* Assembly endpoint sorting must preserve both sides' boundary traversal. */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "skel.h"
#include "impurity.h"
#include "material.h"
#include "matrix.h"

void panic(char *message) { fprintf(stderr,"Unexpected panic: %s\n",message); abort(); }

int main(void)
{
    struct nd_str nodes[2]={0},*node_ptrs[2]={&nodes[0],&nodes[1]};
    struct edg_str edge={0},*edges[1]={&edge};
    int adjacent[1]={0}, membership[2]={0,1};
    int map[4]={3,4,4,1};
    struct LLedge forward={0},backward={0};
    struct sreg_str a={0},b={0};
    nd=node_ptrs;edg=edges;nn=2;ned=1;nsreg=2;ia=map;
    nodes[0].pt=0;nodes[1].pt=1;
    nodes[0].edg.num=nodes[1].edg.num=1;
    nodes[0].edg.list=nodes[1].edg.list=adjacent;
    edge.nd[0]=1;edge.nd[1]=0;edge.skel.num=2;edge.skel.list=membership;
    sreg[0]=&a;sreg[1]=&b;a.bnd=&forward;b.bnd=&backward;
    forward.next=forward.prev=&forward;backward.next=backward.prev=&backward;
    forward.iscc=1;backward.iscc=0;
    assert(nB((&forward))==1 && nF((&forward))==0);
    assert(nB((&backward))==0 && nF((&backward))==1);
    assem_map();
    assert(edge.nd[0]==0 && edge.nd[1]==1);
    assert(nB((&forward))==1 && nF((&forward))==0);
    assert(nB((&backward))==0 && nF((&backward))==1);
    assert(edge.l==3 && edge.u==4);
    assem_map(); /* Repeated mapping must not reverse traversal again. */
    assert(nB((&forward))==1 && nF((&forward))==0);
    assert(nB((&backward))==0 && nF((&backward))==1);
    free(bval);
    puts("PASS: matrix endpoint sorting preserves shared skeleton directions");
    return 0;
}
