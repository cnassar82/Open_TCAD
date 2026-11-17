#include <assert.h>
#include <stdio.h>
#include "pivoted_solver.h"
int main(void)
{
    struct suprem_pivoted p={0};
    int ia[]={4,6,7,7,1,2,2},i,status;
    double a[10]={0},c[]={1e-9,1e16,3},row[]={1e-20,1e22,1e-6};
    double b[3],base[3][3]={{0,2,3},{4,5,6},{7,8,10}};
    for(i=0;i<3;i++)a[i]=row[i]*base[i][i]/c[i];
    a[4]=row[1]*4/c[0];a[5]=row[2]*7/c[0];a[6]=row[2]*8/c[1];
    a[7]=row[0]*2/c[1];a[8]=row[0]*3/c[2];a[9]=row[1]*6/c[2];
    status=suprem_pivoted_factor(&p,3,ia,3,a,c);assert(status==UMFPACK_OK);
    b[0]=13*row[0];b[1]=32*row[1];b[2]=53*row[2];
    assert(suprem_pivoted_solve(&p,b)==UMFPACK_OK);
    for(i=0;i<3;i++)assert(fabs(b[i]/(c[i]*(i+1))-1)<1e-12);
    b[0]=26*row[0];b[1]=64*row[1];b[2]=106*row[2];
    assert(suprem_pivoted_solve(&p,b)==UMFPACK_OK);
    for(i=0;i<3;i++)assert(fabs(b[i]/(2*c[i]*(i+1))-1)<1e-12);
    suprem_pivoted_free(&p);
    /* Singular systems must be reported, never padded or regularized. */
    for(i=0;i<10;i++)a[i]=1;
    assert(suprem_pivoted_factor(&p,3,ia,3,a,c)!=UMFPACK_OK);
    assert(p.numeric==NULL);
    a[0]=NAN;
    assert(suprem_pivoted_factor(&p,3,ia,3,a,c)!=UMFPACK_OK);
    puts("PASS: pivoting, 42-order equation scales, cached RHS, singular/nonfinite rejection");
    return 0;
}
