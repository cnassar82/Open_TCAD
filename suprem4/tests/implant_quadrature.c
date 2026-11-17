#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "global.h"
#include "constant.h"
#include "implant.h"

/* Analytic profiles isolate integration/inversion from moment fitting. */
static int triangular;
double imp_vert(double x,struct imp_info *data)
{
    assert(x>=0 && x<=data->maxz);
    return triangular ? 2-x:1;
}

int main(void)
{
    struct imp_info data;
    double dose=2e15,z,f;
    int i;
    memset(&data,0,sizeof(data));
    data.maxz=2;data.area=2e-4;
    /* Sub-step layers used to receive a negative dose from the final
     * interval: xen-x was measured from the next, unused grid point. */
    assert(fabs(qeqv(.1,.10001,dose,.001,&data)/dose-5e-6)<1e-15);
    assert(fabs(qeqv(-1,3,dose,.003,&data)/dose-1)<1e-12);
    assert(qeqv(2,1,dose,.001,&data)==0);
    for(triangular=0;triangular<=1;triangular++) {
        for(i=0;i<=100;i++) {
            f=i/100.0;
            z=zeqv(f*dose,dose,.003,1,&data);
            assert(isfinite(z) && z>=0 && z<=2);
            assert(fabs(z-(triangular ? 2-2*sqrt(1-f):2*f))<1e-10);
            assert(fabs(qeqv(0,z,dose,.003,&data)/dose-f)<1e-11);
        }
    }
    /* A nearly exhausted triangular profile ends at zero. Its inverse
     * must remain in the last bracket rather than jumping negative. */
    z=zeqv((1-1e-12)*dose,dose,.003,1,&data);
    assert(isfinite(z) && z>1.99999 && z<=2);
    puts("PASS: implant quadrature and inverse preserve dose and depth bounds");
    return 0;
}
