/* Projected refinement must split the remaining edge monotonically. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "skel.h"
static int calls,next_edge=6,next_node=6;
double dabs(double x){return fabs(x);}
void panic(char *s){fprintf(stderr,"%s\n",s);abort();}
int sp_edge(int edge,float *c,int *point,int move) {
 const float expected[3]={.25,.5,.75};int node,newedge;
 assert(calls<3 && c[0]==expected[calls]);calls++;
 assert(c[0]>=cordinate(pt_edg(edge,0),0) && c[0]<=cordinate(pt_edg(edge,1),0));
 if(calls==2){*point=pt_edg(edge,0);return -2;}
 node=next_node++;newedge=next_edge++;
 nd[node]->pt=node;pt[node]->cord[0]=c[0];pt[node]->cord[1]=c[1];*point=node;
 edg[newedge]->nd[0]=node;edg[newedge]->nd[1]=edg[edge]->nd[1];edg[edge]->nd[1]=node;
 return newedge;
}
int main(void) {
 struct pt_str points[16]={0},*pp[16];struct nd_str nodes[16]={0},*npointers[16];
 struct edg_str edges[16]={0},*ep[16];struct LLedge links[6]={0};struct sreg_str r1={0},r2={0};
 float xy[6][2]={{0,0},{1,0},{.75,.2},{.25,.4},{.5,.6},{0,1}};int i;
 pt=pp;nd=npointers;edg=ep;mode=TWOD;
 for(i=0;i<16;i++){pp[i]=&points[i];npointers[i]=&nodes[i];ep[i]=&edges[i];nodes[i].pt=i;}
 for(i=0;i<6;i++){
  points[i].cord[0]=xy[i][0];points[i].cord[1]=xy[i][1];edges[i].nd[0]=i;edges[i].nd[1]=(i+1)%6;
  links[i].edge=i;links[i].iscc=TRUE;links[i].next=&links[(i+1)%6];links[i].ang=PI;
 }
 r1.len=6;r1.bnd=&links[0];r2.bnd=&links[0];
 assert(dvrecedg(&r1,&r2,0,&links[1]));assert(calls==3);
 puts("PASS rectangular refinement: ordered projections and existing-endpoint return");return 0;
}
