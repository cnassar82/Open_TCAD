#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "impurity.h"
#include "matrix.h"
#include "diffuse.h"
#include "charged_flux.h"
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int chrg_block_set();
static double blocks[4][5],right_hand[2];
static void evaluate(double x[2][4],double sign)
{
    double df[2],dp[2],val[2],lv[2],sv[2],dv[4][2],*derivatives[MAXIMP]={0},pot[2],e,f;
    int i,k;
    memset(blocks,0,sizeof(blocks));memset(right_hand,0,sizeof(right_hand));
    for(k=0;k<4;k++)derivatives[k]=dv[k];
    for(i=0;i<2;i++) {
        pot[i]=x[i][3];e=exp(.2*pot[i]);f=1+.1*x[i][1]+.03*x[i][2];
        val[i]=x[i][0]*f*e;dv[0][i]=f*e;dv[1][i]=.1*x[i][0]*e;
        dv[2][i]=.03*x[i][0]*e;dv[3][i]=.2*val[i];
        lv[i]=log(val[i]);sv[i]=sqrt(val[i]);df[i]=1.7*exp(.4*pot[i]);dp[i]=.4*df[i];
    }
    chrg_block_set(0,0,df,dp,val,lv,sv,derivatives,pot,sign);
}
int main(void)
{
    struct edg_str edge={0},*edges[1]={&edge};
    double samples[]={0,1e-12,-1e-12,1e-6,-1e-6,.009,-.009,.01,-.01,.1,-.1,1,-1,50,-50,700,-700,1000,-1000};
    double x[2][4]={{.8,.2,.3,.35},{2.1,.7,.4,-.15}};
    double jac[4][5],plus[2],minus[2],h,analytic,numeric,old,b,bp,d,dp,p,scale;
    double df[2]={1,1},dprt[2]={0},val[2],lv[2],sv[2],dv[4][2]={{1,1},{0},{0},{0}},*derivatives[MAXIMP]={0},pot[2];
    int i,j,k,s,node,row;
    for(i=0;i<sizeof(samples)/sizeof(*samples);i++) {
        long double z=samples[i],reference=z==0?1:z/expm1l(z);
        charged_bernoulli(samples[i],&b,&d);charged_bernoulli(-samples[i],&bp,&dp);
        assert(isfinite(b) && isfinite(d) && b>=0 && d<=0);
        assert(fabs(b-(double)reference)<=32*DBL_EPSILON*fmax(fabs((double)reference),1e-300));
        assert(fabs(bp-b-samples[i])<=32*DBL_EPSILON*fmax(1,fabs(samples[i])));
        assert(fabs(d+dp+1)<1e-12);
    }
    nn=2;ned=1;edg=edges;edge.nd[0]=0;edge.nd[1]=1;edge.cpl=.8;edge.l=3;edge.u=4;
    imptosol[I]=1;imptosol[V]=2;imptosol[Psi]=3;
    for(k=0;k<4;k++){a[0][k]=blocks[k];derivatives[k]=dv[k];}newb[0]=right_hand;
    for(s=-1;s<=1;s+=2) {
        evaluate(x,s);memcpy(jac,blocks,sizeof(jac));assert(right_hand[0]+right_hand[1]==0);
        for(node=0;node<2;node++)for(k=0;k<4;k++) {
            old=x[node][k];h=1e-6*(1+fabs(old));x[node][k]=old+h;evaluate(x,s);memcpy(plus,right_hand,sizeof(plus));
            x[node][k]=old-h;evaluate(x,s);memcpy(minus,right_hand,sizeof(minus));x[node][k]=old;
            for(row=0;row<2;row++) {
                analytic=jac[k][node==row?row:row==0?4:3];numeric=-(plus[row]-minus[row])/(2*h);
                assert(fabs(analytic-numeric)<2e-7*(1+fabs(numeric)));
            }
        }
    }
    /* The captured implant-tail contrast makes the historical geometric-
     * mean drift Jacobian anti-monotone. A fitted edge has nonnegative
     * diagonal and nonpositive off-diagonal concentration coefficients. */
    val[0]=1e5;val[1]=1e20;pot[0]=0;pot[1]=4.3;
    for(i=0;i<2;i++){lv[i]=log(val[i]);sv[i]=sqrt(val[i]);}
    memset(blocks,0,sizeof(blocks));memset(right_hand,0,sizeof(right_hand));
    chrg_block_set(0,0,df,dprt,val,lv,sv,derivatives,pot,1.);
    assert(blocks[0][0]>=0 && blocks[0][1]>=0 && blocks[0][3]<=0 && blocks[0][4]<=0);
    assert(right_hand[0]+right_hand[1]==0);
    for(j=-6;j<=6;j+=6) {
        double dt=pow(10.,j),m0=.7,m1=1.3,det=m0*m1+dt*(m1*blocks[0][0]+m0*blocks[0][1]);
        double c0=(m0*val[0]*(m1+dt*blocks[0][1])-dt*blocks[0][4]*m1*val[1])/det;
        double c1=(m1*val[1]*(m0+dt*blocks[0][0])-dt*blocks[0][3]*m0*val[0])/det;
        assert(c0>=0 && c1>=0);
        assert(fabs((m0*c0+m1*c1)/(m0*val[0]+m1*val[1])-1)<64*DBL_EPSILON);
    }
    val[0]=.3;val[1]=1.7;pot[0]=pot[1]=0;
    for(i=0;i<2;i++){lv[i]=log(val[i]);sv[i]=sqrt(val[i]);}
    memset(blocks,0,sizeof(blocks));memset(right_hand,0,sizeof(right_hand));
    chrg_block_set(0,0,df,dprt,val,lv,sv,derivatives,pot,1.);
    assert(fabs(right_hand[1]-edge.cpl*(val[0]-val[1]))<4*DBL_EPSILON);
    assert(blocks[0][0]==edge.cpl && blocks[0][1]==edge.cpl);
    for(j=-100;j<=100;j+=5) {
        p=j;pot[0]=p;pot[1]=0;
        val[0]=p>0?exp(-p):1;val[1]=p>0?1:exp(p);
        for(i=0;i<2;i++){lv[i]=log(val[i]);sv[i]=sqrt(val[i]);}
        memset(blocks,0,sizeof(blocks));memset(right_hand,0,sizeof(right_hand));
        chrg_block_set(0,0,df,dprt,val,lv,sv,derivatives,pot,1.);
        charged_bernoulli(p,&b,&d);charged_bernoulli(-p,&bp,&dp);scale=.8*(bp*val[0]+b*val[1]);
        assert(fabs(right_hand[0])<=64*DBL_EPSILON*fmax(scale,1e-300));
    }
    puts("PASS: stable Bernoulli limits, actual charged-edge Jacobians, conservation, monotone tail transport and detailed balance");
    return 0;
}
