#include <assert.h>
#include <stdio.h>
#include "time_control.h"

int main(void)
{
    int i;
    for(i=-12;i<=12;i++) {
        double mass=pow(10.,i),rate=3e9,concentration=7e18;
        double dt=diffusion_loss_step(concentration,mass,-rate*mass*concentration);
        /* For scalar decay the TR stage amplification remains positive.
         * The startup estimate must not depend on mesh-volume units. */
        double z=(2-sqrt(2.))*dt*rate;
        assert(fabs(dt-.5/rate)<1e-15/rate);
        assert((1-z/2)/(1+z/2)>0);
    }
    assert(diffusion_loss_step(0,1,-1)==DBL_MAX);
    assert(diffusion_loss_step(1,0,-1)==DBL_MAX);
    assert(diffusion_loss_step(1,1,1)==DBL_MAX);
    assert(diffusion_loss_step(1,1,0)==DBL_MAX);
    assert(diffusion_loss_step(INFINITY,1,-1)==DBL_MAX);
    puts("PASS inventory-based startup: stiff decay positivity, volume scaling, zero/algebraic rows");
    return 0;
}
