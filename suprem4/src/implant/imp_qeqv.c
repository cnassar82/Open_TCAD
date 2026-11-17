/* date: 05 jul 85 (mje)
 *
 * "qeqv" : Find the equivalent dose through a given depth.
 *
 * calling sequence:
 *	double
 *	qeqv(xdep, dx, kp)
 *
 * where:
 *	<return>	- (double) Dose through 'xdep'.
 *	xdep		- (double) Depth through which to calculate dose.
 *	dx		- (double) Delta x to use in calculations.
 *	kp		- (double []) Pearson constants.
 *
 * written: Michael Eldredge (apr 85)
 *
 * imp_qeqv.c	5.1	7/3/91	12:08:40
 */

#include <stdlib.h>
#include <math.h>
#define  abs(V) ((V) < 0.? -(V) : (V))

/* Offsets to the Pearson-IV constants array */
#include "global.h"
#include "constant.h"
#include "implant.h"

double
qeqv(xst, xen, qtot, dx, data)
double	 xst, xen ;
double	 qtot ;		/* total possible dose */
double	 dx ;
struct imp_info	 *data ;
{

    double x, end, next, val, ov, q=0;

    if (!(dx>0) || !(data->area>0) || !isfinite(xst) || !isfinite(xen))
        return NAN;
    x=fmax(0.0,xst);
    end=fmin(xen,data->maxz);
    if(end<=x)return 0;
    ov=imp_vert(x,data);
    while(x<end) {
        next=fmin(x+dx,end);
        if(next<=x)return NAN;
        val=imp_vert(next,data);
        q+=0.5*(val+ov)*(next-x);
        x=next; ov=val;
    }
    return q*1e-4*qtot/data->area;
}
