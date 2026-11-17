/* A reflecting face constrains the normal, before it creates mesh points. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
int nwhich(int t,int n) {int i;for(i=0;i<3;i++)if(tri[t]->nd[i]==n)return i;return -1;}
int local_normal();
static void check(float length,float slope,int turn)
{
    struct pt_str points[3]={0},*pp[3];
    struct nd_str nodes[3]={0},*npointers[3];
    struct tri_str triangle={0},*tp[1]={&triangle};
    int adjacent=0,i,j;float normal[MAXDIM],x,y;
    pt=pp;nd=npointers;tri=tp;mode=TWOD;
    points[1].cord[0]=length;points[1].cord[1]=length*slope;
    points[2].cord[1]=length*.731f;
    for(i=0;i<3;i++) {
        pp[i]=&points[i];npointers[i]=&nodes[i];nodes[i].pt=i;
        triangle.nd[i]=i;triangle.nb[i]=EXPOSED;
        for(j=0;j<turn;j++) {
            x=points[i].cord[0];y=points[i].cord[1];
            points[i].cord[0]=-y;points[i].cord[1]=x;
        }
    }
    triangle.nb[1]=BC_OFFSET;
    nodes[0].ltri.num=1;nodes[0].ltri.list=&adjacent;
    assert(local_normal(0,GAS,normal)==0);
    assert(isfinite(normal[0]) && isfinite(normal[1]));
    /* The constrained component must be exactly zero: it seeds geometry. */
    assert(normal[turn%2]==0.0f);
    assert(fabs(hypot(normal[0],normal[1])-1)<1e-6);
}
int main(void)
{
    int i,r,s;
    for(norm_style=0;norm_style<=1;norm_style++)
    for(s=-30;s<=30;s+=10)for(r=0;r<4;r++)for(i=1;i<=100;i++)
        check(ldexpf(1.0f,s),i*.0000137f,r);
    puts("PASS reflecting normal: 5600 rotated/scaled constraints");return 0;
}
