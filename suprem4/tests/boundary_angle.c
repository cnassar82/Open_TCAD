#include <assert.h>
#include <math.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
extern double intang(int,int,int);
int main(void) {
 struct pt_str points[3]={0},*pp[3];struct nd_str nodes[3]={0},*npointers[3];
 int i,turn,s,e,sign;pt=pp;nd=npointers;
 for(i=0;i<3;i++){pp[i]=&points[i];npointers[i]=&nodes[i];nodes[i].pt=i;}
 for(turn=0;turn<4;turn++)for(s=-20;s<=20;s+=10)for(e=12;e<=36;e+=4)for(sign=-1;sign<=1;sign+=2){
  double expected=sign*atan2(ldexp(1.,-e),1.);int r;
  for(i=0;i<3;i++) {
   float x=ldexpf(i,s),y=i==2?sign*ldexpf(1.f,s-e):0;
   for(r=0;r<turn;r++){float old=x;x=-y;y=old;}
   points[i].cord[0]=x;points[i].cord[1]=y;
  }
  assert(fabs((PI-intang(0,1,2))-expected)<5e-16);
 }
 puts("PASS boundary angles: signed near-collinear turns under rotation and scaling");return 0;
}
