#include <assert.h>
#include <math.h>
#include <stdio.h>
int exp_l(double,double *,double *,double);
int main(void)
{
 double x,l,r,dl,dr,v,d,vp,vm,unused,h=1e-7;
 int i;
 exp_l(-4.6-h,&l,&dl,1);exp_l(-4.6+h,&r,&dr,1);
 /* A constitutive continuation must agree in both value and derivative.
  * A discontinuity here makes Newton's Jacobian inconsistent at the seam. */
 assert(fabs(l/r-1)<3*h);assert(fabs(dl/dr-1)<5*h);
 for(i=0;i<100;i++) {
  x=-20+i*.2;exp_l(x,&v,&d,1);exp_l(x+h,&vp,&unused,1);exp_l(x-h,&vm,&unused,1);
  assert(v>0 && d>0 && isfinite(v) && isfinite(d));
  assert(fabs((vp-vm)/(2*h)/d-1)<1e-6);
 }
 exp_l(0,&v,&d,1);assert(v==1 && d==1);
 puts("PASS oxidation exponential: C1 seam, positive rates and consistent derivatives");return 0;
}
