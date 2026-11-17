/* date: 05 jul 85 (mje)
 *
 * "zeqv" : Find the equivalent dose through a given depth.
 *
 * calling sequence:
 *	double
 *	zeqv(dose, dx, kp)
 *
 * where:
 *	<return>	- (double) depth over which Dose extends.
 *	dose		- (double) Dose over which to find an 'zeqiv'
 *	dx		- (double) Delta x to use in calculations.
 *	kp		- (double []) Pearson constants.
 *
 * notes:
 *	> What happens if the target dose is too large for the implantation
 *	  parameters?
 *
 * written: Michael Eldredge (apr 85)
 *
 * imp_zeqv.c	5.1	7/3/91	12:08:46
 */

#include <stdlib.h>
#include <math.h>
#define  abs(V) ((V) < 0.? -(V) : (V))

/* Offsets to the Pearson-IV constants array */
#include "global.h"
#include "constant.h"
#include "implant.h"

double
zeqv(dose, qtot, dx, bt, data)
	double	 dose ;
	double	 qtot ;
	double	 dx ;
	double bt;
	struct imp_info	 *data ;
	{

    double x=0, next, val, ov, q=0, target, segment, remaining;
    double width, slope, discriminant, fraction, denominator;

    if (!(dx>0) || !(qtot>0) || !(data->area>0) || !isfinite(dose))
        return NAN;
    if(dose<=0)return 0;
    if(dose>=qtot)return data->maxz;
    target=(dose/qtot)*data->area*1e4;
    ov=imp_vert(0.0,data);
    while(x<data->maxz) {
        next=fmin(x+dx,data->maxz);
        if(next<=x)return NAN;
        width=next-x;
        val=imp_vert(next,data);
        segment=0.5*(val+ov)*width;
        if(q+segment>=target) {
            /* Invert this trapezoid, keeping the solution inside its
             * bracket. Dividing the excess by the right endpoint value
             * is unbounded when the profile falls to zero. */
            remaining=(target-q)/width;
            slope=val-ov;
            discriminant=fmax(0.0,ov*ov+2*slope*remaining);
            denominator=ov+sqrt(discriminant);
            fraction=denominator>0 ? 2*remaining/denominator:0;
            fraction=fmin(1.0,fmax(0.0,fraction));
            return x+width*fraction;
        }
        q+=segment;
        x=next; ov=val;
    }
    return data->maxz;
}
