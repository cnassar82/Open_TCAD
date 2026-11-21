/* Differentiate the actual nonlinear oxide law, including the residual-only
 * path used in line searches. The Deal-Grove transport coefficient is held
 * fixed so this test isolates stress/strain and its diffusion multiplier. */
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "material.h"
#include "impurity.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
float proc_temp;
int noderiv;
double lambda;
int FE_oxbulk();
void oxide_staged_homotopy(int);
int do_deal_grove(float temp,int species){B_DealGrove=2;return 0;}
int main(void)
{
 double e[3],direction[3],trial[3],m[40],mp[40],mm[40],s[3],sp[3],sm[3];
 double scale,h=1e-6,analytic,numeric,error=0,denominator;
 int temperature,path,stage,level,trace,nearly_isochoric,j,k;
 mu0[SiO2]=2413.9;mue[SiO2]=2.525;nux[SiO2]=.499;
 stress_dep=1;Vc=300e-24;Vd=25e-24;Dlim=1;pressure=1;gas_type=O2;
 impur[O2].constant[SiO2][0][0]=1;impur[O2].constant[SiO2][1][0]=1;
 for(temperature=0;temperature<2;temperature++)for(path=0;path<2;path++)for(stage=0;stage<5;stage++)for(level=0;level<5;level++)for(trace=-1;trace<=1;trace+=2)for(nearly_isochoric=0;nearly_isochoric<2;nearly_isochoric++) {
  proc_temp=temperature?1373:1173;oxide_staged_homotopy(path);lambda=.25*stage;
  scale=pow(10.,-10+2*level);e[0]=trace*.7*scale;e[1]=trace*(nearly_isochoric?-.699:.4)*scale;e[2]=.3*scale;
  for(k=0;k<3;k++)direction[k]=scale*sin(k+1.);
  noderiv=0;memset(m,0,sizeof(m));FE_oxbulk(m,SiO2,e,s);
  for(k=0;k<3;k++)trial[k]=e[k]+h*direction[k];
  noderiv=1;memset(mp,0,sizeof(mp));FE_oxbulk(mp,SiO2,trial,sp);
  for(k=0;k<3;k++)trial[k]=e[k]-h*direction[k];
  memset(mm,0,sizeof(mm));FE_oxbulk(mm,SiO2,trial,sm);
  for(j=0;j<4;j++) {
   analytic=0;
   for(k=0;k<3;k++)analytic+=(j<3?m[10+3*j+k]:m[21+k])*direction[k];
   numeric=j<3?(sp[j]-sm[j])/(2*h):(mp[20]-mm[20])/(2*h);
   denominator=fmax(fabs(analytic)+fabs(numeric),1e-20);
   if(fabs(analytic-numeric)/denominator>error)error=fabs(analytic-numeric)/denominator;
   assert(isfinite(analytic) && isfinite(numeric));
   if(fabs(analytic-numeric)>=2e-6*denominator)
    fprintf(stderr,"material tangent: T=%g path=%d lambda=%g scale=%g trace=%d nearly_isochoric=%d component=%d analytic=%.17g finite_difference=%.17g relative=%g\n",proc_temp,path,lambda,scale,trace,nearly_isochoric,j,analytic,numeric,fabs(analytic-numeric)/denominator);
   assert(fabs(analytic-numeric)<2e-6*denominator);
  }
  memset(mp,0,sizeof(mp));FE_oxbulk(mp,SiO2,e,sp);
  for(j=0;j<3;j++)assert(s[j]==sp[j]);
  assert(m[20]==mp[20]);
 }
 printf("PASS actual oxide material Jacobian and residual-only consistency, max relative difference %.3g\n",error);return 0;
}
