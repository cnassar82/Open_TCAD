#ifndef SUPREM_DOGLEG_H
#define SUPREM_DOGLEG_H
#include <math.h>
/* Classical two-segment trust-region dogleg in scaled unknowns. The caller
 * supplies the Newton and Cauchy points of the SAME residual model. */
static double suprem_vector_norm(int n,const double *x)
{
    double norm=0;int i;for(i=0;i<n;i++)norm=hypot(norm,x[i]);return norm;
}
static int suprem_dogleg(int n,const double *newton,const double *cauchy,
                         double radius,double *step)
{
    double nn=suprem_vector_norm(n,newton),nc=suprem_vector_norm(n,cauchy),tau;
    long double a=0,b=0,d,root;int i;
    if(!(radius>0) || !isfinite(radius) || !isfinite(nn) || !isfinite(nc))return 0;
    if(nn<=radius){for(i=0;i<n;i++)step[i]=newton[i];return 1;}
    if(nc>=radius) {
        if(nc==0)return 0;
        for(i=0;i<n;i++)step[i]=radius*(cauchy[i]/nc);
        return 2;
    }
    for(i=0;i<n;i++) {
        long double v=(long double)newton[i]-cauchy[i];
        a+=v*v;b+=2*(long double)cauchy[i]*v;
    }
    if(a==0)return 0;
    d=(long double)radius*radius-(long double)nc*nc;
    root=sqrtl(b*b+4*a*d);
    tau=b>=0?(double)(2*d/(b+root)):(double)((root-b)/(2*a));
    if(tau<0)tau=0;if(tau>1)tau=1;
    for(i=0;i<n;i++)step[i]=cauchy[i]+tau*(newton[i]-cauchy[i]);
    return 3;
}
#endif
