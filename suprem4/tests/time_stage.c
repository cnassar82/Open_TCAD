/* Production TR/BDF stage acceptance, with deterministic nonlinear solves. */
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
static int scenario,tr_calls,restores;
int time_val(){abort();return 0;}
int soldif_be(double dt,float temp,double **o,double **n,double **oa,
              double **na,double *a,double *b){abort();return 0;}
void panic(char *message){fprintf(stderr,"%s\n",message);abort();}
int upd_diff(double dt,double *area,int doit,int min_ok)
{area[0]=1;return 0;}
int init_rhs(float temp,double **o,double **oa,double *area)
{assert(o[0][0]==10);return 0;}
int soldif_tr(double dt,float temp,double **o,double **n,double **oa,
              double **na,double *a,double *b,int init)
{
    tr_calls++;
    n[0][0]=tr_calls==1?9:8.5;
    if(scenario==4 && tr_calls==2)n[0][0]=NAN;
    return scenario==1 || (scenario==3 && tr_calls==2);
}
int soldif_bdf(double dt,double olddt,float temp,double **n,double **m,
    double **o,double **na,double **ma,double **oa,double *a,double *b,
    double *c,int init)
{n[0][0]=8;return scenario==2;}
struct mesh_trial *mesh_trial_save(void){abort();return NULL;}
void mesh_trial_restore(const struct mesh_trial *trial)
{assert(!trial);restores++;}
void mesh_trial_discard(struct mesh_trial *trial){assert(!trial);}
int moving_bound(double dt,double gam){abort();return 0;}
int tot_areas(double *a,double *b){abort();return 0;}
extern double timed_trbdf(double,float,int);
int main(int argc,char **argv)
{
    struct nd_str node={0},*nodes[1]={&node};
    double buffers[8][1]={{0}},kt;
    int i;
    scenario=argc>1?atoi(argv[1]):0;
    nd=nodes;nn=1;ned=0;n_imp=tm_nsol=1;tm_sol[0]=0;
    soltoimp[0]=B;impur[B].flags=DIFFUSING;LTE[B]=.01;ABE[B]=1;
    for(i=0;i<MAXIMP;i++)imptosol[i]=-1;
    imptosol[B]=0;
    old[0]=buffers[0];mid[0]=buffers[1];new[0]=buffers[2];
    olda[0]=buffers[3];mida[0]=buffers[4];newa[0]=buffers[5];
    oldb[0]=buffers[6];newb[0]=buffers[7];
    node.sol[0]=10;total=17;verbose=0;
    kt=timed_trbdf(1,1000,1);
    assert(total==17);
    if(scenario) {
        assert(kt<0 && node.sol[0]==10 && restores==1);
        /* A rejected comparison must not poison a fresh retry. */
        scenario=0;tr_calls=0;
        kt=timed_trbdf(.5,1000,1);
    }
    assert(isfinite(kt) && kt>0 && node.sol[0]==8 && total==17);
    puts("PASS primary/comparison stage rejection and fresh retry");return 0;
}
