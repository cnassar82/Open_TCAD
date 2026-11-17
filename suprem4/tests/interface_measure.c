/* Interface quadrature must belong to the face's endpoints. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "matrix.h"
#include "impurity.h"
int sum_coeff();
int main(void)
{
 struct nd_str nodes[6]={0},*npointers[6];
 struct tri_str triangles[2]={0},*tp[2];
 struct reg_str regions[2]={0};
 struct bound_str bounds[3]={0};
 int i,j,one=0;double area[6];
 nd=npointers;tri=tp;bval=bounds;nn=6;ne=2;bcnum=3;
 mode=TWOD;nvrt=nedg=3;
 for(i=0;i<6;i++){npointers[i]=&nodes[i];nodes[i].mater=i<3?Si:SiO2;}
 for(i=0;i<2;i++){
  tp[i]=&triangles[i];reg[i]=&regions[i];regions[i].mater=i?SiO2:Si;
  triangles[i].regnum=i;triangles[i].son=NO_TRIANGLE;
  for(j=0;j<3;j++){triangles[i].nd[j]=3*i+j;triangles[i].nb[j]=BC_OFFSET;triangles[i].earea[j]=1;}
 }
 /* Only face opposite node 0 contacts oxide; its half length is 2. */
 triangles[0].nb[0]=1;triangles[0].d[0]=2;
 for(i=0;i<3;i++) {
  nodes[i].ltri.num=1;nodes[i].ltri.list=&one;
  bounds[i].nx[0]=i;bounds[i].nx[1]=i+3;
 }
 sum_coeff(area);
 assert(bounds[0].cpl==0);
 assert(bounds[1].cpl==2 && bounds[2].cpl==2);
 assert(bounds[0].cpl+bounds[1].cpl+bounds[2].cpl==4);
 /* A free surface obeys the same incidence rule. */
 triangles[0].nb[0]=EXPOSED;
 for(i=3;i<6;i++)nodes[i].mater=GAS;
 sum_coeff(area);
 assert(bounds[0].cpl==0 && bounds[1].cpl==2 && bounds[2].cpl==2);
 /* In 1D face s is the vertex (s+1)%2, by nd_face_ele. */
 mode=ONED;nvrt=nedg=2;bcnum=2;
 triangles[0].nb[0]=EXPOSED;triangles[0].nb[1]=BC_OFFSET;
 triangles[0].d[0]=1;
 sum_coeff(area);
 assert(bounds[0].cpl==0 && bounds[1].cpl==1);
 puts("PASS interface measure: nonincident faces excluded, face length conserved");return 0;
}
