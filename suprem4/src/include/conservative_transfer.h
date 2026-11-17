/* Bounded conservative nodal transfer for local connectivity changes.
 * Caller supplies full material-node P1 integration supports.
 * Identity interpolation preserves affine fields when domain is unchanged;
 * a bounded scalar projection enforces the inventory constraint otherwise.
 */
#ifndef SUPREM_CONSERVATIVE_TRANSFER_H
#define SUPREM_CONSERVATIVE_TRANSFER_H
#include "impurity.h"
#include "material.h"
static int conserved_density(int imp,int material)
{
    /* do_setup solves STEADY and PSEUDO variables algebraically, without
     * a concentration-storage term. This includes the default oxidants. */
    if(IS_PSEUDO(imp) || IS_STEADY(imp))return 0;
    /* Insulator defect entries are equilibrium placeholders (one), not
     * transported point-defect inventories. T is a trap occupancy fraction.
     * Preserve algebraic/fractional fields by interpolation. */
    if(imp==V || imp==I)return !IS_PSEUDO(imp) && (material==Si || material==Poly);
    return imp==As || imp==P || imp==Sb || imp==B || imp==O2 || imp==H2O ||
           imp==Au || imp==Cs || imp==Ga ||
           (imp>=iBe && imp<=iG && (imp&1));
}
static int conserved_parent(int imp)
{
    switch(imp){case Asa:return As;case Pa:return P;case Sba:return Sb;case Ba:return B;}
    return imp>=iBea && imp<=iGa && !(imp&1) ? imp-1:-1;
}
/* Bounded least-change correction, using the complete supports of affected
 * nodes. All preparation precedes the connectivity change. Roundoff-sized residuals
 * are accepted only at the same 1e-12 relative inventory criterion as the
 * existing constrained point-removal transfer. */
static int conserved_transfer(int n,const double *oldw,const double *neww,double *value)
{
    int i,k;double target=0,current=0,lo=HUGE_VAL,hi=-HUGE_VAL,l,h,m=0,sum;
    for(i=0;i<n;i++) {
        if(!isfinite(value[i]) || value[i]<0 || neww[i]<=0)return 0;
        target+=oldw[i]*value[i];current+=neww[i]*value[i];
        lo=fmin(lo,value[i]);hi=fmax(hi,value[i]);
    }
    if(fabs(current-target)<=1e-12*fmax(fabs(target),1e-300))return 1;
    l=lo-hi;h=hi-lo;
    for(k=0;k<100;k++) {
        m=l+(h-l)*.5;sum=0;
        for(i=0;i<n;i++)sum+=neww[i]*fmax(lo,fmin(hi,value[i]+m));
        if(sum<target)l=m;else h=m;
    }
    sum=0;m=l+(h-l)*.5;
    for(i=0;i<n;i++){value[i]=fmax(lo,fmin(hi,value[i]+m));sum+=neww[i]*value[i];}
    return isfinite(sum) && fabs(sum-target)<=1e-12*fmax(fabs(target),1e-300);
}

#endif
