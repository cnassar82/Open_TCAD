/* Exercise the production retry controller with deterministic solver results. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "impurity.h"
#include "matrix.h"
static int calls,scenario;
void panic(char *s){fprintf(stderr,"%s\n",s);abort();}
double timed_trbdf(double dt,float temp,int compute)
{
    int i;
    calls++;
    if(scenario==3) {
        struct nd_str *tmp;
        for(i=0;i<nn;i++) {
            double start=2+8*nd[i]->pt+4*total,end=start+4*dt;
            assert(fabs(nd[i]->sol[MAXIMP-5]-start)<1e-12);
            assert(fabs(nd[i]->sol[0]-end)<1e-12);
            old[0][i]=start;new[0][i]=end;mid[0][i]=.5*(start+end);
        }
        tmp=nd[0];nd[0]=nd[1];nd[1]=tmp;
        return dt>.25?-1:1;
    }
    if(scenario==0)return calls==1?-1:1;
    if(scenario==1)return -1;
    return dt>1.0/128?-1:1;
}
int ig2_write(char *name,int flip,double scale)
{
    assert(!strcmp(name,"step_failure.str"));
    printf("checkpoint calls=%d time=%.17g\n",calls,total);
    return 0;
}
extern double solstep(double,float,int,int);
int main(int argc,char **argv)
{
    double step;
    double arrays[3][3][2]={{{0}}};
    struct nd_str nodes[2]={0},*p[2]={&nodes[0],&nodes[1]};int i;
    scenario=argc>1?atoi(argv[1]):0;n_imp=1;nn=0;tm_nsol=0;total=0;
    imptosol[XVEL]=imptosol[YVEL]=imptosol[DELA]=-1;
    if(scenario==3) {
        nn=2;nd=p;imptosol[XVEL]=0;imptosol[YVEL]=1;imptosol[DELA]=2;
        for(i=0;i<3;i++){old[i]=arrays[0][i];mid[i]=arrays[1][i];new[i]=arrays[2][i];}
        for(i=0;i<2;i++) {
            nodes[i].pt=i;nodes[i].sol[MAXIMP-5]=2+8*i;nodes[i].sol[0]=6+8*i;
        }
        step=solstep(1,1000,TRBDF,1);
        assert(calls==7 && step>0 && total==0);
        for(i=0;i<2;i++)assert(nodes[i].sol[MAXIMP-5]==2+8*i && nodes[i].sol[0]==6+8*i);
        puts("PASS velocity history: nested half steps and node reordering");return 0;
    }
    step=solstep(1,1000,TRBDF,1);
    assert(scenario==0 && calls==3 && step>0 && total==0);
    /* A fresh top-level step must have a new retry budget/depth. */
    calls=0;step=solstep(1,1000,TRBDF,1);
    assert(calls==3 && step>0 && total==0);
    puts("PASS split step succeeds and resets retry accounting");return 0;
}
