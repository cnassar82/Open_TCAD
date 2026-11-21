#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "triangle7.h"
int noderiv=0;
/* The historical default modifies obtuse geometry. The patch test must
 * work with this value too; it must not need a model-card workaround. */
float oxide_obfix=0;
int dcopy(double *a,double *b,int n){memcpy(b,a,n*sizeof(double));return 0;}
int tri7_stiff();int tri7_nodal_stress();int triox_stiff();int triox_nodal_stress();
static double bulk(void){return getenv("SUPREM_TEST_BULK")?atof(getenv("SUPREM_TEST_BULK")):1;}
int oxide(double *m,int mat,double *e,double *s) {
 double lam=bulk();memset(m,0,40*sizeof(double));m[0]=lam+2;m[1]=lam;m[2]=1;
 m[5]=1;m[10]=lam+2;m[11]=lam;m[13]=lam;m[14]=lam+2;m[18]=1;m[20]=1;
 s[0]=(lam+2)*e[0]+lam*e[1];s[1]=lam*e[0]+(lam+2)*e[1];s[2]=e[2];return 0;
}
int elastic(double *m,int mat,double *e,double *s) {
 double lam=bulk();memset(m,0,20*sizeof(double));m[0]=lam+2;m[1]=lam;m[2]=1;
 s[0]=(lam+2)*e[0]+lam*e[1];s[1]=lam*e[0]+(lam+2)*e[1];s[2]=e[2];return 0;
}
static double factorial(int i){return i<2?1:i*factorial(i-1);}
int main(void){
 double vertices[4][3][2]={{{0,0},{1,0},{0,1}},{{0,0},{2,0},{.2,.1}},
                         {{0,0},{1,0},{.01,1e-5}},
                         {{.0006,.0002},{.0007,.0002},{.000601,.000200001}}};
 double reference[7][2]={{0,0},{1,0},{0,1},{.5,.5},{0,.5},{.5,0},{1./3,1./3}};
 double x[7][2],u[7][2],k[14][14],r[7][2],v[7][3],a[21][21],b[7][3],shp[7][3];
 double absolute_sum,absolute_oxide_sum;
 double sum,exact,area,energy,expect,epsx=.7,epsy=-.2,shear=.4,oxenergy,px,py;
 int i,j,q,powerx,powery,caseid;
 for(powerx=0;powerx<=5;powerx++)for(powery=0;powery<=5-powerx;powery++) {
  sum=0;for(q=0;q<TRIANGLE7_NQ;q++)sum+=triangle7_q[q][2]*pow(triangle7_q[q][0],powerx)*pow(triangle7_q[q][1],powery);
  exact=factorial(powerx)*factorial(powery)/factorial(powerx+powery+2);
  assert(fabs(sum-exact)<1e-14);
 }
 for(q=0;q<TRIANGLE7_NQ;q++) {
  triangle7_basis(triangle7_q[q][0],triangle7_q[q][1],shp);
  for(j=0;j<3;j++) {
   sum=0;px=0;py=0;for(i=0;i<7;i++){sum+=shp[i][j];px+=reference[i][0]*shp[i][j];py+=reference[i][1]*shp[i][j];}
   assert(fabs(sum-(j==2?1:0))<1e-14);
   assert(fabs(px-(j==0?1:j==2?triangle7_q[q][0]:0))<1e-14);
   assert(fabs(py-(j==1?1:j==2?triangle7_q[q][1]:0))<1e-14);
  }
 }
 /* An interior bubble has zero corner strain. The historical linear
  * extrapolation produced a false corner strain and reaction feedback. */
 {
  double eps[7][3],sigma[7][3],elastic_eps[7][3],jac;
  memset(v,0,sizeof(v));memset(u,0,sizeof(u));v[6][0]=u[6][0]=1;
  triox_nodal_stress(reference,v,eps,sigma,oxide,1,&jac);
  tri7_nodal_stress(reference,u,elastic_eps,sigma,elastic,1,&jac);
  for(i=0;i<7;i++) {
   double s=reference[i][0],t=reference[i][1],l0=1-s-t;
   double dx=9*t*(l0-s),projected=-1.8*(s-l0);
   assert(fabs(eps[i][0]-(getenv("SUPREM4_FE_BBAR_P1")?.5*(dx+projected):dx))<1e-13);
   assert(fabs(eps[i][1]-(getenv("SUPREM4_FE_BBAR_P1")?.5*(projected-dx):0))<1e-13);
   assert(fabs(eps[i][2]-9*s*(l0-t))<1e-13);
   assert(fabs(elastic_eps[i][0]-dx)<1e-13);
   assert(fabs(elastic_eps[i][1])<1e-13);
   assert(fabs(elastic_eps[i][2]-9*s*(l0-t))<1e-13);
  }
 }
 for(caseid=0;caseid<4;caseid++) {
  memset(k,0,sizeof(k));memset(a,0,sizeof(a));memset(r,0,sizeof(r));memset(b,0,sizeof(b));
  for(i=0;i<3;i++)memcpy(x[i],vertices[caseid][i],sizeof(x[i]));
  for(i=3;i<6;i++)for(j=0;j<2;j++)x[i][j]=.5*(x[(i+1)%3][j]+x[(i+2)%3][j]);
  for(j=0;j<2;j++)x[6][j]=(x[0][j]+x[1][j]+x[2][j])/3;
  for(i=0;i<7;i++){u[i][0]=epsx*(x[i][0]-x[0][0])+shear*(x[i][1]-x[0][1]);u[i][1]=epsy*(x[i][1]-x[0][1]);v[i][0]=u[i][0];v[i][1]=u[i][1];v[i][2]=0;}
  tri7_stiff(k,r,x,u,elastic,1);triox_stiff(a,b,x,v,oxide,1);
  area=((x[1][0]-x[0][0])*(x[2][1]-x[0][1])-(x[1][1]-x[0][1])*(x[2][0]-x[0][0]))/2;
  expect=area*(bulk()*(epsx+epsy)*(epsx+epsy)+2*(epsx*epsx+epsy*epsy)+shear*shear);
  energy=oxenergy=absolute_sum=absolute_oxide_sum=0;
  for(i=0;i<14;i++)for(j=0;j<14;j++){double term=u[i/2][i%2]*k[i][j]*u[j/2][j%2];
   double oxterm=u[i/2][i%2]*a[3*(i/2)+i%2][3*(j/2)+j%2]*u[j/2][j%2];
   energy+=term;oxenergy+=oxterm;absolute_sum+=fabs(term);absolute_oxide_sum+=fabs(oxterm);}
  /* Use an absolute rounding bound for cancellation in u^T K u on
   * extremely thin elements, not a widened physical relative tolerance. */
  assert(fabs(energy-expect)<128*DBL_EPSILON*(absolute_sum+fabs(expect)));
  assert(fabs(oxenergy-expect)<128*DBL_EPSILON*(absolute_oxide_sum+fabs(expect)));
  if(caseid==0) {
   for(i=0;i<14;i++)for(j=0;j<14;j++)printf("%.17g%c",k[i][j],j==13?'\n':' ');
   for(i=0;i<21;i++)for(j=0;j<21;j++)printf("%.17g%c",a[i][j],j==20?'\n':' ');
  }
 }
 return 0;
}
