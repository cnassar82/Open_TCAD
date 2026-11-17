/* Conservative transfer contract for inserted shared-interface points. */
#include <assert.h>
#define STATIC_ALLOCATION_TIME
#include "../src/dbase/vertical_constraints.c"
int main(void)
{
    double oldw[4]={1,2,3,0},neww[4]={.8,1.9,2.8,.5};
    double c[4]={1,3,5,2},constant[4]={7,7,7,7},impossible[4]={1,1,1,1};
    double expanded[4]={1,2,3,1},target=0,result=0;int i;
    for(i=0;i<4;i++)target+=oldw[i]*c[i];
    assert(vc_transfer(4,oldw,neww,c));
    for(i=0;i<4;i++){assert(c[i]>=1 && c[i]<=5);result+=neww[i]*c[i];}
    assert(fabs(result-target)<1e-12*target);
    assert(vc_transfer(4,oldw,neww,constant));
    for(i=0;i<4;i++)assert(constant[i]==7);
    assert(!vc_transfer(4,oldw,expanded,impossible));
    assert(vc_density(B,SiO2) && vc_density(P,Poly));
    assert(!vc_density(V,SiO2) && !vc_density(I,SiO2) && !vc_density(T,Poly));
    impur[V].flags|=PSEUDO;assert(!vc_density(V,Poly));
    impur[V].flags&=~PSEUDO;assert(vc_density(V,Poly));
    impur[O2].flags=STEADY|LOCKSTEP;assert(!vc_density(O2,SiO2));
    impur[O2].flags=0;assert(vc_density(O2,SiO2));
    puts("PASS constraint transfer: inventory, positivity, constant fields, infeasible bounds, algebraic defect semantics");
    return 0;
}
