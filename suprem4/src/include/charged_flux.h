#ifndef SUPREM_CHARGED_FLUX_H
#define SUPREM_CHARGED_FLUX_H
#include <math.h>
/* B(x)=x/(exp(x)-1), including its analytic derivative. Stable limits are
 * numerical evaluations of the same function, not coefficient floors. */
static void charged_bernoulli(double x,double *b,double *derivative)
{
    double e,d,x2;
    if(fabs(x)<.01) {
        x2=x*x;
        *b=1-.5*x+x2*(1./12+x2*(-1./720+x2/30240));
        *derivative=-.5+x*(1./6+x2*(-1./180+x2/5040));
    } else if(x>50) {
        e=exp(-x);d=1-e;
        *b=x*e/d;*derivative=e*(1-x-e)/(d*d);
    } else if(x< -50) {
        e=exp(x);d=1-e;
        *b=-x/d;*derivative=(-1+e-x*e)/(d*d);
    } else {
        d=expm1(x);e=d+1;
        *b=x/d;*derivative=(d-x*e)/(d*d);
    }
}
#endif
