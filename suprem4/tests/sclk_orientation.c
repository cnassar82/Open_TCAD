/* Vertex ordering (sclk_tri) and the geometry kernel (tri_geom) must agree
 * on orientation for every valid triangle, however thin. A float
 * determinant in sclk_tri once left a captured RIT NMOS sliver (two vertices
 * one ULP apart) reversed; tri_geom's double evaluation then panicked with
 * "triangles are not clock wise". Both now use the exact predicate. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
void panic(char *s){fprintf(stderr,"panic: %s\n",s);abort();}
int sclk_tri();
static void check(const float xy[3][2])
{
    struct pt_str points[3]={0},*pp[3];
    struct nd_str nodes[3]={0},*npointers[3];
    struct tri_str triangle={0},*tp[1]={&triangle};
    int i,sign;
    pt=pp;nd=npointers;tri=tp;mode=TWOD;
    for(i=0;i<3;i++) {
        pp[i]=&points[i];npointers[i]=&nodes[i];nodes[i].pt=i;
        memcpy(points[i].cord,xy[i],sizeof(xy[i]));triangle.nd[i]=i;
        triangle.nb[i]=10+i;triangle.edg[i]=20+i;
    }
    sign=orient2d_float_sign(points[0].cord,points[1].cord,points[2].cord);
    assert(sign!=0);
    sclk_tri(0);
    /* After ordering, the stored order is positive by the exact predicate,
     * and neighbour/edge slots were permuted with their vertices. */
    assert(orient2d_float_sign(points[triangle.nd[0]].cord,points[triangle.nd[1]].cord,
                               points[triangle.nd[2]].cord)>0);
    for(i=0;i<3;i++) {
        assert(triangle.nb[i]==10+triangle.nd[i]);
        assert(triangle.edg[i]==20+triangle.nd[i]);
    }
    tri_geom(0);
    for(i=0;i<3;i++)assert(isfinite(triangle.earea[i]) && triangle.earea[i]>=0);
}
int main(void)
{
    /* The captured flipped triangle: tri_flip(it=95133), NMOS photoresist
     * deposit after the spacer etch. */
    const float sliver[3][2]={
        {0.00063261966f,-4.26869155e-05f},
        {0.000632619718f,-4.26869155e-05f},
        {0.00013276501f,-7.40009054e-05f}};
    static const int perm[6][3]={{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    float t[3][2];
    int p,i,j,k,s;
    for(p=0;p<6;p++)for(s=-20;s<=20;s+=5)for(k=0;k<3;k++) {
        /* Permutations, binary scalings and translations. */
        for(i=0;i<3;i++)for(j=0;j<2;j++)
            t[i][j]=ldexpf(sliver[perm[p][i]][j],s)+(k==1?ldexpf(1.0f,s-8):k==2?-ldexpf(1.0f,s-9):0);
        if(orient2d_float_sign(t[0],t[1],t[2]))check(t);
    }
    /* Random near-collinear slivers: a third vertex one to a few ULPs off a
     * long edge. The float determinant misorders a share of these. */
    srand(12345);
    for(i=0;i<20000;i++) {
        float a[2],b[2],c[2];
        double u=(double)rand()/RAND_MAX;
        a[0]=(float)((double)rand()/RAND_MAX*1e-3);a[1]=(float)((double)rand()/RAND_MAX*1e-4-5e-5);
        b[0]=(float)((double)rand()/RAND_MAX*1e-3);b[1]=(float)((double)rand()/RAND_MAX*1e-4-5e-5);
        c[0]=(float)(a[0]+u*((double)b[0]-a[0]));c[1]=(float)(a[1]+u*((double)b[1]-a[1]));
        c[rand()%2]=nextafterf(c[rand()%2],rand()%2?1.0f:-1.0f);
        memcpy(t[0],a,sizeof a);memcpy(t[1],b,sizeof b);memcpy(t[2],c,sizeof c);
        if(orient2d_float_sign(t[0],t[1],t[2]) && !(a[0]==b[0] && a[1]==b[1]))check(t);
    }
    puts("PASS sclk_tri and tri_geom agree on exact orientation");return 0;
}
