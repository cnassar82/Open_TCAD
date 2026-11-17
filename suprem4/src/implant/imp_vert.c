/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *     Copyright c 1991 The board of regents of the University of 	 *
 *                      Florida.  All rights reserved.			 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of the University of 	 *
 *     Florida.								 *
 *									 *
 *									 *
 *************************************************************************/
/*   imp_vert.c                Version 5.1     */
/*   Last Modification : 7/3/91 12:08:45 */

#include <stdio.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "implant.h"

/************************************************************************
 *									*
 *	imp_vert() - This function returns a value between 0 and 1 	*
 *  corresponding to the value of the distribution function at the 	*
 *  given depth.							*
 *									*
 *  Original:	MEL	1/87						*
 *									*
 ************************************************************************/
double imp_vert(depth, dat)
double	depth ;
struct imp_info	*dat;
{
    double val;

    if ( imp_model == PEARS ) {
	if ( depth < dat->maxz ) {
	    val = prson(depth - dat->vert[PRS_RP], dat);
	}
	else 
	    val = 0.0;
    }
    else if ( imp_model == GAUSS ) {
	if ( depth < dat->maxz ) {
	    val = (depth - dat->vert[GUS_RP]) / dat->vert[GUS_SIG];
	    val = exp( - val * val * 0.5 );
	}
	else val = 0.0;
    }
    return( val );
}


double func(l, pt, dat)
double l[MAXDIM];
double pt[MAXDIM];
struct imp_info	*dat;
{
    return( imp_vert((l[0] - pt[0])*1.0e4, dat) * imp_latr((l[1] - pt[1])*1.0e4, dat) );
}



double simp_int( l, a, b, eps, dat )
double l[MAXDIM];
double a[MAXDIM];        
double b[MAXDIM];        
double eps;
struct imp_info *dat;
{
    int i;
    double s, sbar;
    double ab[MAXDIM], ab3[MAXDIM], a3b[MAXDIM];
    double h = b[1] - a[1];
    double fa, fab, fb, fa3, fb3;
    double integ = 0.0;

    for(i = 0; i < mode; i++) {
	ab[i] = 0.5 * (a[i] + b[i]);
	ab3[i] = 0.75*a[i]+0.25*b[i];
	a3b[i] = 0.75*b[i]+0.25*a[i];
    }

    /*Evaluate using the adaptive Simpson's rule*/
    fa = func(l, a, dat);
    fab = func(l, ab, dat);
    fb = func(l, b, dat);
    fa3 = func(l, ab3, dat);
    fb3 = func(l, a3b, dat);

    s = h * (fa + 4.0 * fab + fb) / 6.0;
    sbar = h * (fa + 4.0*fa3 + 2*fab + 4*fb3 + fb) / 12.0;

    if ( fabs(s - sbar) > 15.0 * eps ) {
	integ = simp_int(l, a, ab, eps, dat) + simp_int(l, ab, b, eps, dat);
    }
    else
	integ = sbar;

    return(integ);
}


/* Seed quadrature at both physical peaks.  An unseeded adaptive rule can
 * report zero error when all its probes miss a narrow positive integrand.
 * Gaussian lateral width and vertical peak curvature provide scales even
 * for a tilted source segment.  Integrate the complete segment: no tail
 * clipping or renormalization of the requested dose is used. */
static double source_integral(double *l,double *a,double *b,struct imp_info *dat)
{
    static const double offsets[]={-32,-16,-8,-4,-2,-1,0,1,2,4,8,16,32};
    double cuts[28],p[3]={0},q[3]={0},slope,peak,width,x,tmp,value=0;
    int n=2,i,j,k;
    cuts[0]=a[1];cuts[1]=b[1];slope=(b[0]-a[0])/(b[1]-a[1]);
    for(i=0;i<13;i++) {
        x=l[1]+offsets[i]*dat->lat[LATSIG]*1e-4;
        if(x>a[1] && x<b[1])cuts[n++]=x;
    }
    if(imp_model==GAUSS) {
        peak=dat->vert[GUS_RP];width=dat->vert[GUS_SIG];
    }
    else {
        double z=dat->vert[PRS_A0];
        peak=dat->vert[PRS_RP]+z;
        width=sqrt(fabs(dat->vert[PRS_B0]+z*z*(1+dat->vert[PRS_B2])));
    }
    if(slope!=0) {
        peak=a[1]+(l[0]-a[0]-peak*1e-4)/slope;
        width=width*1e-4/fabs(slope);
        for(i=0;i<13;i++) {
            x=peak+offsets[i]*width;
            if(x>a[1] && x<b[1])cuts[n++]=x;
        }
    }
    for(i=1;i<n;i++) {
        tmp=cuts[i];j=i;
        while(j>0 && cuts[j-1]>tmp) {cuts[j]=cuts[j-1];j--;}
        cuts[j]=tmp;
    }
    for(k=1;k<n;k++) {
        if(cuts[k]==cuts[k-1])continue;
        p[1]=cuts[k-1];q[1]=cuts[k];
        p[0]=a[0]+(p[1]-a[1])*slope;
        q[0]=a[0]+(q[1]-a[1])*slope;
        value+=simp_int(l,p,q,1e-12/(n-1),dat);
    }
    return value;
}

double evaluate(l, lft, rht, dat)
double l[MAXDIM];  
double lft[MAXDIM];        
double rht[MAXDIM];        
struct imp_info *dat;
{
    double a[MAXDIM], b[MAXDIM];
    double val;

    /*This is confusing due to a mismatch - parallel to beam is
      cooridinate 0, perp is cord 0*/
    a[0] = lft[1] - 0.5 * (lft[1] + rht[1]);
    a[1] = lft[0];
    b[0] = rht[1] - 0.5 * (lft[1] + rht[1]);
    b[1] = rht[0];

    if ( mode == ONED ) {
	val = imp_vert(l[0]*1.0e4, dat);
    }
    else if ( mode == TWOD ) {
        if(lft[1]==rht[1]) {
            /* A flat source segment separates exactly into a vertical
             * profile and a Gaussian lateral integral.  Point-sampled
             * Simpson quadrature can completely miss a narrow lateral
             * peak on a long segment, regardless of its error estimate. */
            double scale=sqrt(2.0)*dat->lat[LATSIG]*1e-4;
            double left=(l[1]-lft[0])/scale;
            double right=(l[1]-rht[0])/scale;
            double lateral;
            if(right>=0)lateral=0.5*(erfc(right)-erfc(left));
            else if(left<=0)lateral=0.5*(erfc(-left)-erfc(-right));
            else lateral=0.5*(erf(left)-erf(right));
            val=imp_vert(l[0]*1e4,dat)*lateral;
        }
        else val = source_integral(l,a,b,dat);
    }
    return(val);
}
