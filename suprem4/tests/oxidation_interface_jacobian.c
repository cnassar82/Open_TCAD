/* Differentiate the assembled residual at an oxide/nitride/Si junction.
 * A perturbation confined to the nitride still changes the averaged strain
 * used by the oxide reaction law. Test that dependency explicitly. */
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "FEgeom.h"
#include "material.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int noderiv;
int triox_stiff(),triox_bc(),triox_nodal_stress(),FEassmb();
int dzero(double *a,int n){memset(a,0,n*sizeof(*a));return 0;}
int dcopy(double *a,double *b,int n){memcpy(b,a,n*sizeof(*a));return 0;}
int iaja(int *ia,int off,int row,int col) {
 int i,swap;
 if(row==col)return row;
 if(row<col){swap=row;row=col;col=swap;}else off=0;
 for(i=ia[row];i<ia[row+1];i++)if(ia[i]==col)return i+off;
 abort();
}
int FE_oxbulk(double *m,int mat,double *e,double *s) {
 memset(m,0,40*sizeof(*m));m[0]=3;m[1]=1;m[2]=1;m[5]=1;
 m[10]=3;m[11]=1;m[13]=1;m[14]=3;m[18]=1;m[20]=1;
 s[0]=3*e[0]+e[1];s[1]=e[0]+3*e[1];s[2]=e[2];return 0;
}
int ksubs(double *k,double *s,double *normal,double *extra) {
 memset(k,0,11*sizeof(*k));k[0]=.1*exp(.1*s[0]);
 k[3]=1;k[4]=.5;k[5]=exp(.1*s[0]);return 0;
}
static int node(double x,double y) {
 int i;for(i=0;i<FEnn;i++)if(FEnd[i]->cord[0]==x && FEnd[i]->cord[1]==y)return i;
 i=FEnn++;FEnd[i]=calloc(1,sizeof(**FEnd));FEnd[i]->cord[0]=x;FEnd[i]->cord[1]=y;
 FEnd[i]->sol[0]=.01*x*x;FEnd[i]->sol[1]=.01*y*y;
 FEnd[i]->sol[2]=2+.1*x+.2*y+.01*x*x;
 return i;
}
int main(int argc,char **argv) {
 double vertices[4][3][2]={{{0,0},{1,0},{0,1}},{{0,0},{0,1},{-1,0}},{{0,0},{0,-1},{1,0}},{{0,0},{-1,0},{0,-1}}};
 double *a,*rhs,*plus,*minus,h=1e-5,analytic,numeric,error=0,ref=0;
 int *ia,i,j,e,n,off,next,center,perturb,row,col,fail=0;
 FEnd=calloc(30,sizeof(*FEnd));FEelt=calloc(4,sizeof(*FEelt));FEne=argc>1?4:3;FEdf=3;FEdm=2;FEnsc=3;
 FEdesc[0].nel=7;FEdesc[0].stiff=triox_stiff;FEdesc[0].coeff=FE_oxbulk;
 FEdesc[0].bound=triox_bc;FEdesc[0].nodal_stress=triox_nodal_stress;FEdesc[1].nel=7;
 for(e=0;e<FEne;e++) {
  FEelt[e]=calloc(1,sizeof(**FEelt));FEelt[e]->desc=e==2?1:0;
  FEelt[e]->mat=e==1?SiNi:e==2?Si:SiO2;
  for(i=0;i<3;i++){FEelt[e]->nd[i]=node(vertices[e][i][0],vertices[e][i][1]);FEelt[e]->face[i]=-1;}
  for(i=3;i<6;i++)FEelt[e]->nd[i]=node(.5*(vertices[e][(i+1)%3][0]+vertices[e][(i+2)%3][0]),.5*(vertices[e][(i+1)%3][1]+vertices[e][(i+2)%3][1]));
  FEelt[e]->nd[6]=node((vertices[e][0][0]+vertices[e][1][0]+vertices[e][2][0])/3,(vertices[e][0][1]+vertices[e][1][1]+vertices[e][2][1])/3);
 }
 FEelt[0]->face[2]=2;if(FEne==4)FEelt[3]->face[1]=2;
 center=node(0,0);perturb=node(-1,0);
 FEnd[center]->fixity[0]=FEnd[center]->fixity[1]=4;FEnd[center]->garbage[0]=1;
 if(getenv("SUPREM4_FE_LINEAR_TRACE")) {
  int l,in;
  for(e=0;e<FEne;e++)if(FEelt[e]->mat==SiO2)for(i=0;i<3;i++)if(FEelt[e]->face[i]==2) {
   int local[3];local[0]=(i+1)%3;local[1]=(i+2)%3;local[2]=i+3;
   for(l=0;l<3;l++){in=FEelt[e]->nd[local[l]];FEnd[in]->fixity[0]=FEnd[in]->fixity[1]=4;FEnd[in]->garbage[0]=1;}
  }
 }
 n=FEnn*3;off=n*(n-1)/2;ia=calloc(n+1+off,sizeof(*ia));next=n+1;
 for(i=0;i<n;i++){ia[i]=next;for(j=0;j<i;j++)ia[next++]=j;}ia[n]=next;
 a=calloc(next+off,sizeof(*a));rhs=calloc(n,sizeof(*rhs));plus=calloc(n,sizeof(*plus));minus=calloc(n,sizeof(*minus));
 noderiv=0;FEassmb(ia,off,a,rhs);
 /* Probe every active velocity column, not just the nitride-only one.
  * The mixed-order experiment also checks scalar columns and all rows,
  * including interpolation constraints away from their solution manifold. */
 for(col=0;col<n;col++)if(col%3!=2 || getenv("SUPREM4_FE_OXIDANT_P1")) {
  FEnd[col/3]->sol[col%3]+=h;noderiv=1;FEassmb(ia,off,a,plus);
  FEnd[col/3]->sol[col%3]-=2*h;FEassmb(ia,off,a,minus);
  FEnd[col/3]->sol[col%3]+=h;
  for(j=0;j<(getenv("SUPREM4_FE_OXIDANT_P1")?n:3);j++) {
   row=getenv("SUPREM4_FE_OXIDANT_P1")?j:3*center+j;i=row-3*center;
   analytic=a[iaja(ia,off,row,col)];numeric=-(plus[row]-minus[row])/(2*h);
   error+=(analytic-numeric)*(analytic-numeric);ref+=numeric*numeric;
   if(col==3*perturb && (i==0 || i==2)) {
    printf("nitride-only perturbation: row=%d analytic=%.12g finite_difference=%.12g\n",i,analytic,numeric);
    if(fabs(analytic-numeric)>1e-8*(1+fabs(numeric)))fail=1;
   }
  }
 }
 printf("interface Jacobian relative error: %.12g\n",sqrt(error/ref));
 if(sqrt(error/ref)>1e-8)fail=1;
 for(e=0;e<FEne;e++)free(FEelt[e]);for(i=0;i<FEnn;i++)free(FEnd[i]);
 free(FEelt);free(FEnd);free(ia);free(a);free(rhs);free(plus);free(minus);
 return fail;
}
