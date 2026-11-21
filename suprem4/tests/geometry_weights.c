/* Valid triangles must yield finite, conservative dual-cell areas. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
void panic(char *s){fprintf(stderr,"%s\n",s);abort();}
static void check(const float xy[3][2],int obtuse)
{
    struct pt_str points[3]={0},*pp[3];
    struct nd_str nodes[3]={0},*npointers[3];
    struct tri_str triangle={0},*tp[1]={&triangle};
    double area,sum=0;
    int i;
    pt=pp;nd=npointers;tri=tp;
    for(i=0;i<3;i++) {
        pp[i]=&points[i];npointers[i]=&nodes[i];nodes[i].pt=i;
        memcpy(points[i].cord,xy[i],sizeof(xy[i]));triangle.nd[i]=i;
    }
    area=.5*(((double)xy[1][0]-xy[0][0])*((double)xy[2][1]-xy[0][1])-
             ((double)xy[1][1]-xy[0][1])*((double)xy[2][0]-xy[0][0]));
    assert(area>0);
    /* Before any CHP repair: an obtuse triangle whose edge couplings sum
     * positive keeps its own ehed, but its control volumes must still be
     * positive and sum to its area (the legacy |cot| pieces did not). */
    tri_geom(0);
    for(i=0;i<3;i++) {
        assert(isfinite(triangle.earea[i]) && triangle.earea[i]>=0);
        sum+=triangle.earea[i];
    }
    assert(fabs(sum-area)<2e-6*area);
    if(obtuse)chp_fix(0);
    for(sum=0,i=0;i<3;i++) {
        assert(isfinite(triangle.earea[i]) && triangle.earea[i]>=0);
        assert(isfinite(triangle.ehed[i]) && triangle.ehed[i]>=0);
        sum+=triangle.earea[i];
    }
    assert(fabs(sum-area)<2e-6*area);
}
int main(void)
{
    /* Captured valid PMOS triangles, failed-step checkpoint at 1716.12 s.
     * 6448's acute-angle denominators round to zero in the old CHP formula. */
    const float narrow[3][2]={
        {0.0006372862262651324f,2.293843135703355e-05f},
        {0.0006372864590957761f,-2.7061571017839015e-05f},
        {0.000637286517303437f,-2.7061571017839015e-05f}};
    const float right[3][2]={
        {0.0007891380810178816f,2.5165551051031798e-05f},
        {0.000782337156124413f,2.5165551051031798e-05f},
        {0.0007891380810178816f,2.5164283215417527e-05f}};
    const float acute[3][2]={{0,0},{2,0},{1,1.5}};
    const float obtuse[3][2]={{0,0},{2,0},{.25,.25}};
    float shifted[3][2];int i,j,p,scale;
    /* Flat fan element: a refined interface node over a distant bulk
     * node (thin-film barrier deck); legacy volume was ~10^4 x its area. */
    const float fan[3][2]={{.25f,0},{.265625f,0},{0,.0049f}};
    check(narrow,1);check(right,0);check(acute,0);check(obtuse,1);check(fan,1);
    /* Cyclic ordering and scale changes exercise the formula, not one ID. */
    for(p=0;p<3;p++)for(scale=-12;scale<=12;scale+=6) {
        for(i=0;i<3;i++)for(j=0;j<2;j++)shifted[i][j]=ldexpf(narrow[(i+p)%3][j],scale);
        check(shifted,1);
    }
    for(i=1;i<=200;i++) {
        double a,b,c;
        shifted[0][0]=shifted[0][1]=shifted[1][1]=0;
        shifted[1][0]=2;
        shifted[2][0]=(i%31)*.125f-.875f;
        shifted[2][1]=(i%23+1)*.0625f;
        a=4;b=shifted[2][0]*shifted[2][0]+shifted[2][1]*shifted[2][1];
        c=(shifted[2][0]-2)*(shifted[2][0]-2)+shifted[2][1]*shifted[2][1];
        check(shifted,a>b+c || b>a+c || c>a+b);
    }
    puts("PASS finite positive geometry weights and area conservation");return 0;
}
