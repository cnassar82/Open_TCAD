/* Flux assembly must be a read-only operation on geometry. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "diffuse.h"
#include "material.h"
#include "impurity.h"

int main(void)
{
    struct pt_str points[4]={0}, *pp[4];
    struct nd_str nodes[6]={0}, *npointers[6];
    struct tri_str triangles[2]={0}, *tp[2];
    struct reg_str regions[2]={0};
    double ov[3][6]={{0}},mv[3][6]={{0}},nv[3][6]={{0}};
    float saved[4][MAXDIM];
    int incident=0, i,j,k,repeat;
    const double dt=21.9393,scale=2-sqrt(2.0);
    double expected[3][2],before_flux[3][6],*stages[3][3];
    pt=pp;nd=npointers;tri=tp;np=4;nn=6;ne=2;nreg=2;
    mode=TWOD;nvrt=nedg=3;
    for(i=0;i<MAXIMP;i++)imptosol[i]=-1;
    imptosol[XVEL]=0;imptosol[YVEL]=1;imptosol[DELA]=2;
    for(i=0;i<3;i++) {
        old[i]=ov[i];mid[i]=mv[i];new[i]=nv[i];
        stages[0][i]=ov[i];stages[1][i]=mv[i];stages[2][i]=nv[i];
    }
    for(i=0;i<4;i++){pt[i]=&points[i];pt[i]->nn=i<2?2:1;}
    pt[0]->nd[0]=0;pt[0]->nd[1]=1;pt[1]->nd[0]=2;pt[1]->nd[1]=3;
    pt[2]->nd[0]=4;pt[3]->nd[0]=5;
    pt[0]->cord[0]=0.0004;pt[0]->cord[1]=2.3e-5;
    pt[1]->cord[0]=0.000402;pt[1]->cord[1]=2.3001e-5;
    for(i=0;i<4;i++)memcpy(saved[i],pt[i]->cord,sizeof saved[i]);
    for(i=0;i<6;i++) {
        nd[i]=&nodes[i];nd[i]->mater=i%2?SiO2:Si;
        nd[i]->pt=i<4?i/2:i-2;
        for(j=0;j<2;j++) {
            nv[j][i]=(i%2?1:-1)*(i+1)*(j+1)*1.37e-11;
            ov[j][i]=nv[j][i]*.75;
            mv[j][i]=ov[j][i]*(1-scale)+nv[j][i]*scale;
            nd[i]->sol[j]=nv[j][i];
        }
    }
    for(i=0;i<2;i++) {
        tri[i]=&triangles[i];reg[i]=&regions[i];reg[i]->mater=i?SiO2:Si;
        tri[i]->regnum=i;tri[i]->nb[0]=1-i;tri[i]->nb[1]=tri[i]->nb[2]=-1;
    }
    tri[0]->nd[0]=4;tri[0]->nd[1]=0;tri[0]->nd[2]=2;
    tri[1]->nd[0]=5;tri[1]->nd[1]=3;tri[1]->nd[2]=1;
    nd[0]->ltri.num=nd[2]->ltri.num=1;
    nd[0]->ltri.list=nd[2]->ltri.list=&incident;
    /* Independently form the virtual endpoint positions at each stage. */
    for(k=0;k<3;k++) {
        double back=k==0?1:k==1?1-scale:0;
        double hx=(double)saved[0][0]-back*dt*nv[0][0];
        double hy=(double)saved[0][1]-back*dt*nv[1][0];
        double tx=(double)saved[1][0]-back*dt*nv[0][2];
        double ty=(double)saved[1][1]-back*dt*nv[1][2];
        for(i=0;i<2;i++)expected[k][i]=.5*(stages[k][0][2*i]*(ty-hy)+stages[k][1][2*i]*(hx-tx));
    }
    for(repeat=0;repeat<100;repeat++) {
        moving_bound(dt,scale);
        for(i=0;i<4;i++)assert(memcmp(saved[i],pt[i]->cord,sizeof saved[i])==0);
        for(k=0;k<3;k++) {
            for(i=0;i<2;i++)assert(fabs(stages[k][2][2*i]-expected[k][i])<1e-10*fabs(expected[k][i]));
            for(i=0;i<6;i++) {
                if(repeat)assert(stages[k][2][i]==before_flux[k][i]);
                before_flux[k][i]=stages[k][2][i];
            }
        }
    }
    /* A separate interface across the opposite face cannot sweep volume
     * at a node which is not on that face. Give that node a coincident
     * oxide node so the production pair traversal actually exercises it. */
    pt[2]->nn=2;pt[2]->nd[1]=5;nd[5]->pt=2;
    nd[4]->ltri.num=1;nd[4]->ltri.list=&incident;
    moving_bound(dt,scale);
    for(k=0;k<3;k++)assert(stages[k][2][4]==0.0);
    puts("PASS: analytic moving-interface flux, repeatable values, bitwise unchanged coordinates");
    return 0;
}
