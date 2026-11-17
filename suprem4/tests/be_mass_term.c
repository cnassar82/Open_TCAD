/* The real BE time term balances flux against material inventory. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "matrix.h"
int time_val();
int main(void)
{
 struct call_str call={0};
 double old_values[2]={10,2},new_values[2],old_area[2]={1,3},new_area[2]={1,3};
 double flux[2],residual[2],diagonal[2],*old_pointer[1]={old_values},*new_pointer[1]={new_values};
 double norm=0,dt=5,k=2,det;
 nn=2;a[0][0]=diagonal;rhs[0]=residual;newb[0]=flux;
 call.old=old_pointer;call.new=new_pointer;call.oco=old_area;call.nco=new_area;call.new_del=dt;call.type=BE;
 det=(old_area[0]+dt*k)*(old_area[1]+dt*k)-dt*dt*k*k;
 new_values[0]=((old_area[1]+dt*k)*old_area[0]*old_values[0]+dt*k*old_area[1]*old_values[1])/det;
 new_values[1]=(dt*k*old_area[0]*old_values[0]+(old_area[0]+dt*k)*old_area[1]*old_values[1])/det;
 diagonal[0]=diagonal[1]=k;
 flux[0]=k*(new_values[1]-new_values[0]);flux[1]=-flux[0];
 time_val(0,&norm,&call);
 assert(fabs(residual[0])+fabs(residual[1])<1e-13 && norm<1e-24);
 assert(fabs(new_values[0]+3*new_values[1]-16)<1e-13);
 assert(diagonal[0]>0 && diagonal[1]>0);
 /* Changing volumes alone must preserve each node's stored inventory. */
 new_area[0]=.25;new_area[1]=7;
 new_values[0]=old_area[0]*old_values[0]/new_area[0];
 new_values[1]=old_area[1]*old_values[1]/new_area[1];
 flux[0]=flux[1]=diagonal[0]=diagonal[1]=norm=0;
 time_val(0,&norm,&call);
 assert(fabs(residual[0])+fabs(residual[1])<1e-13);
 puts("PASS BE mass term: analytic two-volume diffusion and changing-volume inventory");return 0;
}
