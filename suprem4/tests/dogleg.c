#include <assert.h>
#include <stdio.h>
#include "dogleg.h"
static void residual(double *x,double *f,double j[2][2],int kind)
{
    if(!kind) {
        f[0]=10*(x[1]-x[0]*x[0]);f[1]=1-x[0];
        j[0][0]=-20*x[0];j[0][1]=10;j[1][0]=-1;j[1][1]=0;
    } else {
        f[0]=exp(x[0])+x[0]+.2*x[1]-1;f[1]=.2*x[0]+2*exp(x[1])+x[1]-2;
        j[0][0]=exp(x[0])+1;j[0][1]=.2;j[1][0]=.2;j[1][1]=2*exp(x[1])+1;
    }
}
int main(void)
{
    double newton[2]={2,1},cauchy[2]={.2,.4},step[2],x[2],f[2],j[2][2],trial[2],ft[2],jt[2][2];
    double g[2],ag[2],det,radius,before,after,predicted,rho,gg,aa;
    int kind,i,k,it,search,accepted;
    assert(suprem_dogleg(2,newton,cauchy,3,step)==1);
    assert(step[0]==2 && step[1]==1);
    assert(suprem_dogleg(2,newton,cauchy,.1,step)==2);
    assert(fabs(suprem_vector_norm(2,step)-.1)<1e-15);
    assert(suprem_dogleg(2,newton,cauchy,1,step)==3);
    assert(fabs(suprem_vector_norm(2,step)-1)<1e-15);
    assert(!suprem_dogleg(2,newton,cauchy,NAN,step));
    for(kind=0;kind<2;kind++) {
        x[0]=kind?-2:-1.2;x[1]=kind?1.5:1;radius=1;
        for(it=0;it<150;it++) {
            residual(x,f,j,kind);before=suprem_vector_norm(2,f);
            if(before<1e-10)break;
            det=j[0][0]*j[1][1]-j[0][1]*j[1][0];assert(det!=0);
            newton[0]=(-j[1][1]*f[0]+j[0][1]*f[1])/det;
            newton[1]=(j[1][0]*f[0]-j[0][0]*f[1])/det;
            for(i=0;i<2;i++){g[i]=0;for(k=0;k<2;k++)g[i]-=j[k][i]*f[k];}
            for(i=0;i<2;i++)ag[i]=j[i][0]*g[0]+j[i][1]*g[1];
            gg=g[0]*g[0]+g[1]*g[1];aa=ag[0]*ag[0]+ag[1]*ag[1];assert(aa>0);
            for(i=0;i<2;i++)cauchy[i]=g[i]*gg/aa;
            accepted=0;
            for(search=0;search<60;search++) {
                assert(suprem_dogleg(2,newton,cauchy,radius,step));
                for(i=0;i<2;i++){trial[i]=x[i]+step[i];ag[i]=f[i]+j[i][0]*step[0]+j[i][1]*step[1];}
                predicted=before*before-ag[0]*ag[0]-ag[1]*ag[1];
                residual(trial,ft,jt,kind);after=suprem_vector_norm(2,ft);
                rho=(before-after)*(before+after)/predicted;
                if(rho<.25 || !isfinite(rho))radius=.25*suprem_vector_norm(2,step);
                else if(rho>.75 && suprem_vector_norm(2,step)>.9*radius)radius*=2;
                if(predicted>0 && isfinite(after) && after<before && rho>.001) {
                    x[0]=trial[0];x[1]=trial[1];accepted=1;break;
                }
            }
            assert(accepted);
        }
        assert(it<150);assert(fabs(x[0]-(kind?0:1))<1e-8);assert(fabs(x[1]-(kind?0:1))<1e-8);
        printf("PASS trust region: nonlinear control %d reached the known root in %d iterations\n",kind,it);
    }
    return 0;
}
