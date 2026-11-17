#ifndef SUPREM_PIVOTED_SOLVER_H
#define SUPREM_PIVOTED_SOLVER_H
/* Scaled, pivoted sparse LU for the paired-triangle SUPREM matrix format.
 * Scaling changes units only: A x=b becomes R A C y=R b, x=C y.
 * No diagonal replacement, pivot floor, or material parameter is introduced.
 * UMFPACK owns symbolic ordering, numerical pivoting and iterative refinement.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <suitesparse/umfpack.h>
struct suprem_pivoted {
    int n,*ap,*ai;
    double *ax,*row_scale,*column_scale,*rhs,*solution;
    void *numeric;
};
static void suprem_pivoted_free(struct suprem_pivoted *p)
{
    if(p->numeric)umfpack_di_free_numeric(&p->numeric);
    free(p->ap);free(p->ai);free(p->ax);free(p->row_scale);
    free(p->column_scale);free(p->rhs);free(p->solution);
    memset(p,0,sizeof(*p));
}
static int suprem_pivoted_factor(struct suprem_pivoted *p,int n,int *ia,
                                int aoff,double *a,double *scale)
{
    int *ti=NULL,*tj=NULL,nz,k,j,z=0,status=UMFPACK_ERROR_invalid_matrix;
    double *tx=NULL,*maximum=NULL,control[UMFPACK_CONTROL];
    void *symbolic=NULL;
    suprem_pivoted_free(p);
    if(n<=0 || ia[0]!=n+1 || ia[n]<n+1 ||
       ia[n]-n-1>(INT_MAX-n)/2)return status;
    nz=n+2*(ia[n]-n-1);p->n=n;
    ti=malloc((size_t)nz*sizeof(int));tj=malloc((size_t)nz*sizeof(int));
    tx=malloc((size_t)nz*sizeof(double));maximum=calloc(n,sizeof(double));
    p->ap=malloc(((size_t)n+1)*sizeof(int));p->ai=malloc((size_t)nz*sizeof(int));
    p->ax=malloc((size_t)nz*sizeof(double));p->row_scale=malloc((size_t)n*sizeof(double));
    p->column_scale=malloc((size_t)n*sizeof(double));p->rhs=malloc((size_t)n*sizeof(double));
    p->solution=malloc((size_t)n*sizeof(double));
    if(!ti || !tj || !tx || !maximum || !p->ap || !p->ai || !p->ax ||
       !p->row_scale || !p->column_scale || !p->rhs || !p->solution) {
        status=UMFPACK_ERROR_out_of_memory;goto done;
    }
    for(k=0;k<n;k++) {
        if(!(scale[k]>0) || !isfinite(scale[k]))goto done;
        p->column_scale[k]=scale[k];
        ti[z]=tj[z]=k;tx[z++]=a[k];
        if(ia[k]>ia[k+1])goto done;
        for(j=ia[k];j<ia[k+1];j++) {
            if(ia[j]<=k || ia[j]>=n)goto done;
            ti[z]=ia[j];tj[z]=k;tx[z++]=a[j];
            ti[z]=k;tj[z]=ia[j];tx[z++]=a[j+aoff];
        }
    }
    if(z!=nz)goto done;
    for(z=0;z<nz;z++) {
        tx[z]*=p->column_scale[tj[z]];
        if(!isfinite(tx[z]))goto done;
        if(fabs(tx[z])>maximum[ti[z]])maximum[ti[z]]=fabs(tx[z]);
    }
    for(k=0;k<n;k++) {
        if(!(maximum[k]>0)) {status=UMFPACK_WARNING_singular_matrix;goto done;}
        p->row_scale[k]=1/maximum[k];
        if(!isfinite(p->row_scale[k]))goto done;
    }
    for(z=0;z<nz;z++)tx[z]*=p->row_scale[ti[z]];
    status=umfpack_di_triplet_to_col(n,n,nz,ti,tj,tx,p->ap,p->ai,p->ax,NULL);
    if(status!=UMFPACK_OK)goto done;
    umfpack_di_defaults(control);control[UMFPACK_PRL]=0;
    status=umfpack_di_symbolic(n,n,p->ap,p->ai,p->ax,&symbolic,control,NULL);
    if(status!=UMFPACK_OK)goto done;
    status=umfpack_di_numeric(p->ap,p->ai,p->ax,symbolic,&p->numeric,control,NULL);
 done:
    if(symbolic)umfpack_di_free_symbolic(&symbolic);
    free(ti);free(tj);free(tx);free(maximum);
    if(status!=UMFPACK_OK)suprem_pivoted_free(p);
    return status;
}
static int suprem_pivoted_solve(struct suprem_pivoted *p,double *rhs)
{
    int k,status;
    double control[UMFPACK_CONTROL];
    if(!p->numeric)return UMFPACK_ERROR_invalid_Numeric_object;
    for(k=0;k<p->n;k++) {
        p->rhs[k]=rhs[k]*p->row_scale[k];
        if(!isfinite(p->rhs[k]))return UMFPACK_ERROR_invalid_matrix;
    }
    umfpack_di_defaults(control);control[UMFPACK_PRL]=0;
    status=umfpack_di_solve(UMFPACK_A,p->ap,p->ai,p->ax,p->solution,p->rhs,
                          p->numeric,control,NULL);
    if(status!=UMFPACK_OK)return status;
    for(k=0;k<p->n;k++) {
        rhs[k]=p->solution[k]*p->column_scale[k];
        if(!isfinite(rhs[k]))return UMFPACK_ERROR_invalid_matrix;
    }
    return UMFPACK_OK;
}
#endif
