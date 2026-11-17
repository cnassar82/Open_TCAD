#ifndef SUPREM_TIME_CONTROL_H
#define SUPREM_TIME_CONTROL_H
#include <math.h>
#include <float.h>
/* Startup bound for d(M C)/dt = F. Limit the explicit depletion estimate
 * to half the current inventory. Sources impose no positivity bound; the
 * time integrator continues to control their temporal error. */
static double diffusion_loss_step(double concentration,double mass,double flux)
{
    long double step;
    if(!(concentration>0) || !(mass>0) || !(flux<0) ||
       !isfinite(concentration) || !isfinite(mass) || !isfinite(flux))return DBL_MAX;
    step=.5L*(long double)concentration*(long double)mass/-(long double)flux;
    if(step>DBL_MAX)return DBL_MAX;
    return (double)step;
}
#endif
