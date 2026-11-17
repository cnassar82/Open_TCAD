/* Production BE controller, analytic conservative two-volume diffusion. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "matrix.h"
#include "mesh_trial.h"
static int calls,fail_call;
static double rate=1;
void panic(char *s){fprintf(stderr,"%s\n",s);abort();}
int time_val(){return 0;}
int upd_diff(double dt,double *area,int doit,int min_ok)
{area[0]=area[1]=1;return 0;}
int init_rhs(float t,double **o,double **oa,double *area){return 0;}
int soldif_be(double dt,float temp,double **o,double **n,double **oa,
              double **na,double *a,double *b)
{
 double sum=o[0][0]+o[0][1],difference=(o[0][0]-o[0][1])/(1+2*rate*dt);
 n[0][0]=(sum+difference)*.5;n[0][1]=(sum-difference)*.5;
 return ++calls==fail_call;
}
int soldif_tr(){abort();return 0;}
int soldif_bdf(){abort();return 0;}
struct mesh_trial *mesh_trial_save(void){abort();return NULL;}
void mesh_trial_restore(const struct mesh_trial *t){assert(!t);}
void mesh_trial_discard(struct mesh_trial *t){assert(!t);}
int moving_bound(){abort();return 0;}
int tot_areas(){abort();return 0;}
extern double timed_trbdf(double,float,int);
int main(void)
{
 struct nd_str nodes[2]={0},*p[2]={&nodes[0],&nodes[1]};
 double buffers[8][2]={{0}},kt,e1,e2;
 int i,fail;
 nd=p;nn=2;ned=0;n_imp=tm_nsol=1;tm_sol[0]=0;soltoimp[0]=B;
 impur[B].flags=DIFFUSING;impur[B].time_val=time_val;
 for(i=0;i<MAXIMP;i++)imptosol[i]=-1;
 imptosol[B]=0;old[0]=buffers[0];mid[0]=buffers[1];new[0]=buffers[2];
 olda[0]=buffers[3];mida[0]=buffers[4];newa[0]=buffers[5];
 oldb[0]=buffers[6];newb[0]=buffers[7];LTE[B]=.01;ABE[B]=1e-6;
 setenv("SUPREM4_TIME_INTEGRATOR","be",1);total=17;
 for(fail=1;fail<=3;fail++) {
  nodes[0].sol[0]=1;nodes[1].sol[0]=0;calls=0;fail_call=fail;
  kt=timed_trbdf(.001,1000,1);
  assert(kt<0 && nodes[0].sol[0]==1 && nodes[1].sol[0]==0 && total==17);
 }
 fail_call=0;nodes[0].sol[0]=1;nodes[1].sol[0]=0;
 /* Maximum local LTE rejects a large inaccurate trial without publishing. */
 assert(timed_trbdf(1,1000,0)<0 && nodes[0].sol[0]==1 && total==17);
 LTE[B]=10;
 assert(timed_trbdf(.01,1000,1)>0);
 e1=fabs(nodes[0].sol[0]-(1+exp(-.02))*.5);
 nodes[0].sol[0]=1;nodes[1].sol[0]=0;
 assert(timed_trbdf(.005,1000,1)>0);
 e2=fabs(nodes[0].sol[0]-(1+exp(-.01))*.5);
 assert(e1/e2>3.8 && e1/e2<4.1);
 /* Stiff mode is damped with positive concentrations and exact inventory. */
 rate=1e8;nodes[0].sol[0]=1;nodes[1].sol[0]=0;
 assert(timed_trbdf(1,1000,1)>0);
 assert(nodes[0].sol[0]>0 && nodes[1].sol[0]>0);
 assert(fabs(nodes[0].sol[0]+nodes[1].sol[0]-1)<1e-14);
 assert(fabs(nodes[0].sol[0]-.5)<1e-12 && total==17);
 puts("PASS BE: stage rollback, LTE rejection, convergence, stiff positivity and conserved inventory");return 0;
}
