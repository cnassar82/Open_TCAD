#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "global.h"
#include "constant.h"
#include "implant.h"
int mode=TWOD;
int imp_model;
double evaluate(double *,double *,double *,struct imp_info *);
double prson(double x,struct imp_info *d) { return exp(-x*x); }
int main(void)
{
    struct imp_info d;
    double l[3]={0,0,0},a[3]={0,0,0},b[3]={0,0,0},v,w;
    memset(&d,0,sizeof(d));imp_model=GAUSS;
    d.maxz=100;d.vert[GUS_RP]=0;d.vert[GUS_SIG]=1;
    d.lat[LATSIG]=.01;d.lat[LATAREA]=sqrt(2*PI)*.01;
    /* A narrow Gaussian far from every Simpson probe on a broad source.
     * The integral is unity, even when a sampled error estimate says zero. */
    a[0]=-1e-2;b[0]=1e-2;l[1]=.003123;
    v=evaluate(l,a,b,&d);assert(fabs(v-1)<1e-14);
    a[0]=0;b[0]=1e-3;l[1]=0;
    assert(fabs(evaluate(l,a,b,&d)-.5)<1e-14);
    /* Stable tails and invariance under subdivision and translation. */
    a[0]=0;b[0]=1e-6;l[1]=-1e-5;
    v=evaluate(l,a,b,&d);
    assert(v>0 && fabs(v/(.5*(erfc(10/sqrt(2.0))-erfc(11/sqrt(2.0))))-1)<1e-12);
    a[0]=0;b[0]=.000017;l[1]=.0000081;v=evaluate(l,a,b,&d);
    b[0]=.000006;w=evaluate(l,a,b,&d);a[0]=b[0];b[0]=.000017;w+=evaluate(l,a,b,&d);
    assert(fabs(v-w)<1e-14);
    w=evaluate(l,a,b,&d);
    a[0]+=.0123;b[0]+=.0123;l[1]+=.0123;
    assert(fabs(evaluate(l,a,b,&d)-w)<1e-11);
    /* Tilted source: compare with the closed-form convolution of two
     * Gaussians, with both peaks between the old quadrature probes. */
    d.vert[GUS_SIG]=.001;
    a[0]=-.01;b[0]=.01;a[1]=-.005;b[1]=.005;
    l[1]=.003123;l[0]=.5*l[1];
    v=evaluate(l,a,b,&d);
    assert(fabs(v-1/sqrt(1+25.0))<1e-10);
    puts("PASS: exact flat-source implant integral, narrow peaks, tails and subdivisions");
    return 0;
}
