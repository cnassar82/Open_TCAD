/* A Delaunay patch can contain an obtuse triangle. Its mass must still
 * equal its physical area; positive net edge couplings do not repair mass. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
void panic(char *s){fprintf(stderr,"%s\n",s);abort();}
static void check(int scale,int turn)
{
    struct pt_str points[4]={0},*pp[4];
    struct nd_str nodes[4]={0},*npointers[4];
    struct tri_str triangles[2]={0},*tp[2];
    const float base[4][2]={{0,0},{2,0},{1,.25},{1,-5}};
    const int vertices[2][3]={{0,1,2},{0,3,1}};
    double mass=0,linear_mass=0,expected=0;
    float x,y;int i,j,r;
    pt=pp;nd=npointers;tri=tp;
    for(i=0;i<4;i++) {
        pp[i]=&points[i];npointers[i]=&nodes[i];nodes[i].pt=i;
        x=ldexpf(base[i][0],scale);y=ldexpf(base[i][1],scale);
        for(r=0;r<turn;r++){float old=x;x=-y;y=old;}
        points[i].cord[0]=x;points[i].cord[1]=y;
    }
    for(i=0;i<2;i++) {
        double area,mean=0;float *a,*b,*c;
        tp[i]=&triangles[i];
        for(j=0;j<3;j++)triangles[i].nd[j]=vertices[i][j];
        tri_geom(i);
        a=points[vertices[i][0]].cord;b=points[vertices[i][1]].cord;c=points[vertices[i][2]].cord;
        area=.5*(((double)b[0]-a[0])*((double)c[1]-a[1])-((double)b[1]-a[1])*((double)c[0]-a[0]));
        for(j=0;j<3;j++) {
            double value=2+ldexp((double)points[vertices[i][j]].cord[0],-scale);
            assert(triangles[i].earea[j]>0 && isfinite(triangles[i].earea[j]));
            mass+=triangles[i].earea[j];linear_mass+=value*triangles[i].earea[j];mean+=value/3;
        }
        expected+=area*mean;
    }
    /* All boundary couplings and the shared-edge sum are nonnegative:
     * geom's legacy negative-edge repair will leave this patch unchanged. */
    assert(triangles[0].ehed[0]>0 && triangles[0].ehed[1]>0);
    assert(triangles[1].ehed[0]>0 && triangles[1].ehed[2]>0);
    assert(triangles[0].ehed[2]+triangles[1].ehed[1]>0);
    assert(fabs(mass-ldexp(5.25,2*scale)) < ldexp(5.25e-7,2*scale));
    assert(fabs(linear_mass-expected) < 1e-6*fabs(expected)+ldexp(1e-7,2*scale));
    /* Constrained obtuse treatment may change stiffness, never P1 mass. */
    chp_fix(0);
    for(j=0;j<3;j++)assert(fabs(triangles[0].earea[j]-ldexp(.25/3,2*scale))<ldexp(1e-8,2*scale));
}
int main(void)
{int s,r;setenv("SUPREM4_MASS_LUMPING","p1",1);for(s=-24;s<=12;s+=6)for(r=0;r<4;r++)check(s,r);puts("PASS P1 mass: area and linear-field integration on 28 Delaunay obtuse patches");return 0;}
