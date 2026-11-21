/* Read-only, region-constrained cavity proposal. A boundary point may be
 * removed only when it subdivides one exactly straight boundary segment.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mesh_patch.h"
extern int orient2d_float_sign(const float *, const float *, const float *);

static int sign3(const struct mesh_patch_point *p, int a, int b, int c)
{
    float x[2], y[2], z[2];
    memcpy(x,p[a].xy,sizeof x); memcpy(y,p[b].xy,sizeof y); memcpy(z,p[c].xy,sizeof z);
    return orient2d_float_sign(x,y,z);
}
static double quality(const struct mesh_patch_point *p, int a, int b, int c)
{
    double x=(double)p[b].xy[0]-p[a].xy[0], y=(double)p[b].xy[1]-p[a].xy[1];
    double u=(double)p[c].xy[0]-p[a].xy[0], v=(double)p[c].xy[1]-p[a].xy[1];
    double denom=x*x+y*y+u*u+v*v+(x-u)*(x-u)+(y-v)*(y-v);
    return denom>0 ? 2.0*sqrt(3.0)*(x*v-y*u)/denom : 0;
}
static int overlap(const struct mesh_patch_point *p, const int *a, const int *b)
{
    int side, i, j, separated;
    for(side=0;side<2;side++) {
        const int *u=side ? b:a, *v=side ? a:b;
        for(i=0;i<3;i++) {
            separated=1;
            for(j=0;j<3;j++) if(sign3(p,u[i],u[(i+1)%3],v[j])>0) separated=0;
            if(separated)return 0;
        }
    }
    return 1;
}
static int segment_bad(const struct mesh_patch_point *p,int a,int b,int c,int d)
{
    int ac=sign3(p,a,b,c),ad=sign3(p,a,b,d);
    int ca=sign3(p,c,d,a),cb=sign3(p,c,d,b), axis, i;
    int end[4]={a,b,c,d};
    if(ac*ad<0 && ca*cb<0)return 1;
    if(ac==0 && ad==0 && ca==0 && cb==0) {
        double lo,hi;
        axis=p[a].xy[0]!=p[b].xy[0] ? 0:1;
        lo=fmax(fmin(p[a].xy[axis],p[b].xy[axis]),fmin(p[c].xy[axis],p[d].xy[axis]));
        hi=fmin(fmax(p[a].xy[axis],p[b].xy[axis]),fmax(p[c].xy[axis],p[d].xy[axis]));
        if(hi>lo)return 1;
    }
    for(i=0;i<4;i++) {
        int e=end[i],u=i<2?c:a,v=i<2?d:b;
        if(e==u || e==v)continue;
        if(sign3(p,u,v,e)==0 &&
           p[e].xy[0]>=fmin(p[u].xy[0],p[v].xy[0]) &&
           p[e].xy[0]<=fmax(p[u].xy[0],p[v].xy[0]) &&
           p[e].xy[1]>=fmin(p[u].xy[1],p[v].xy[1]) &&
           p[e].xy[1]<=fmax(p[u].xy[1],p[v].xy[1]))return 1;
    }
    return 0;
}
/* Visibility tests use only exact signs on the stored coordinates. */
static int patch_in_cone(const struct mesh_patch_point *p,const int *order,
                         int n,int i,int j)
{
    int a=order[i],b=order[(i+1)%n],c=order[(i+n-1)%n],v=order[j];
    if(sign3(p,a,b,c)>=0)return sign3(p,a,v,c)>0 && sign3(p,v,a,b)>0;
    return !(sign3(p,a,v,b)>=0 && sign3(p,v,a,c)>=0);
}
/* Bottleneck-quality dynamic programming over all valid polygon diagonals.
 * Q(i,j)=max_k min(Q(i,k),Q(k,j),quality(i,k,j)). Unlike choosing the best
 * current ear, this optimizes the worst triangle of the complete cavity.
 * Coordinates, boundary constraints and the acceptance score are unchanged. */
static int patch_optimal_triangulation(const struct mesh_patch_point *p,
    const int *order,int n,int region,struct mesh_patch_result *r)
{
    double *score=NULL,q,candidate;int *split=NULL,*tasks=NULL;
    unsigned char *visible=NULL;size_t cells;int i,j,k,a,b,gap,pending=0;
    int status=MESH_PATCH_INVALID;
    if(n<3 || (size_t)n>((size_t)-1)/(size_t)n)return status;
    cells=(size_t)n*n;
    if(cells>((size_t)-1)/sizeof(double) || cells>((size_t)-1)/sizeof(int) ||
       (size_t)n>((size_t)-1)/(2*sizeof(int)))return status;
    score=malloc(cells*sizeof(double));split=malloc(cells*sizeof(int));
    visible=calloc(cells,sizeof(*visible));tasks=malloc((size_t)2*n*sizeof(int));
    if(!score || !split || !visible || !tasks){status=MESH_PATCH_NOMEM;goto out;}
    for(i=0;i<n;i++)for(j=0;j<n;j++){score[(size_t)i*n+j]=-1;split[(size_t)i*n+j]=-1;}
    for(i=0;i<n;i++)for(j=i+1;j<n;j++) {
        int valid=(j==i+1 || (i==0 && j==n-1));
        if(!valid && patch_in_cone(p,order,n,i,j) && patch_in_cone(p,order,n,j,i)) {
            for(k=0;k<n;k++)if(segment_bad(p,order[i],order[j],order[k],order[(k+1)%n]))break;
            valid=k==n;
        }
        visible[(size_t)i*n+j]=valid;
        if(j==i+1)score[(size_t)i*n+j]=HUGE_VAL;
    }
    for(gap=2;gap<n;gap++)for(i=0;i+gap<n;i++) {
        j=i+gap;if(!visible[(size_t)i*n+j])continue;
        for(k=i+1;k<j;k++) {
            if(score[(size_t)i*n+k]<0 || score[(size_t)k*n+j]<0 ||
               sign3(p,order[i],order[k],order[j])<=0)continue;
            q=quality(p,order[i],order[k],order[j]);
            candidate=fmin(q,fmin(score[(size_t)i*n+k],score[(size_t)k*n+j]));
            if(candidate>score[(size_t)i*n+j]) {
                score[(size_t)i*n+j]=candidate;split[(size_t)i*n+j]=k;
            }
        }
    }
    if(score[n-1]<0)goto out;
    tasks[0]=0;tasks[1]=n-1;pending=1;r->count=0;r->new_min_quality=HUGE_VAL;
    while(pending) {
        pending--;i=tasks[2*pending];j=tasks[2*pending+1];k=split[(size_t)i*n+j];
        if(k<=i || k>=j || r->count>=n-2)goto out;
        a=r->count++;
        r->triangle[a].vertex[0]=order[i];r->triangle[a].vertex[1]=order[k];
        r->triangle[a].vertex[2]=order[j];r->triangle[a].region=region;
        q=quality(p,order[i],order[k],order[j]);r->new_min_quality=fmin(r->new_min_quality,q);
        for(b=0;b<2;b++) {
            int first=b?k:i,last=b?j:k;
            if(last==first+1)continue;
            if(pending>=n)goto out;
            tasks[2*pending]=first;tasks[2*pending+1]=last;pending++;
        }
    }
    if(r->count==n-2)status=MESH_PATCH_OK;
out:
    free(score);free(split);free(visible);free(tasks);return status;
}

void mesh_patch_release(struct mesh_patch_result *r)
{
    free(r->triangle);r->triangle=NULL;r->count=0;
}
int mesh_patch_propose(int center,int np,const struct mesh_patch_point *p,int nt,
    const struct mesh_patch_triangle *t,struct mesh_patch_result *r)
{
    int *next=NULL,*incoming=NULL,*order=NULL,*seen=NULL;
    int i,j,k,a,b,c,n,start=-1,status=MESH_PATCH_INVALID,region;
    int open_start=-1,open_end=-1,boundary=0;
    double q,minq=1.0;
    memset(r,0,sizeof *r);r->reason="invalid input";
    if(np<4 || nt<1 || center<0 || center>=np || !p || !t)return status;
    if((size_t)np>(size_t)-1/sizeof(int) || (size_t)nt>(size_t)-1/sizeof(*t))return status;
    for(i=0;i<np;i++)if(!isfinite(p[i].xy[0]) || !isfinite(p[i].xy[1]))return status;
    next=malloc((size_t)np*sizeof(int));incoming=calloc(np,sizeof(int));
    order=malloc((size_t)np*sizeof(int));seen=calloc(np,sizeof(int));
    if(!next || !incoming || !order || !seen){status=MESH_PATCH_NOMEM;r->reason="proposal allocation";goto out;}
    for(i=0;i<np;i++)next[i]=-1;
    region=t[0].region;
    for(i=0;i<nt;i++) {
        int at=-1;
        for(j=0;j<3;j++) {
            a=t[i].vertex[j];
            if(a<0 || a>=np)goto out;
            if(a==center)at=j;
        }
        if(at<0 || sign3(p,t[i].vertex[0],t[i].vertex[1],t[i].vertex[2])<=0)goto out;
        if(t[i].region!=region){status=MESH_PATCH_PINNED;r->reason="region constraint";goto out;}
        a=t[i].vertex[(at+1)%3];b=t[i].vertex[(at+2)%3];
        if(next[a]>=0 || incoming[b]++)goto out;
        next[a]=b;
        if(start<0 || a<start)start=a;
        q=quality(p,t[i].vertex[0],t[i].vertex[1],t[i].vertex[2]);if(q<minq)minq=q;
    }
    for(i=0;i<nt;i++)for(j=i+1;j<nt;j++)if(overlap(p,t[i].vertex,t[j].vertex))goto out;
    for(i=0;i<np;i++)if((next[i]>=0)!=(incoming[i]>0)) {
        if(next[i]>=0) {
            if(open_start>=0)goto out;
            open_start=i;
        }else {
            if(open_end>=0)goto out;
            open_end=i;
        }
    }
    if(open_start>=0 || open_end>=0) {
        if(open_start<0 || open_end<0)goto out;
        if(sign3(p,open_end,center,open_start)!=0) {
            status=MESH_PATCH_PINNED;r->reason="boundary constraint";goto out;
        }
        for(i=0;i<2;i++)if(p[center].xy[i]<fmin(p[open_start].xy[i],p[open_end].xy[i]) ||
            p[center].xy[i]>fmax(p[open_start].xy[i],p[open_end].xy[i]))goto out;
        /* Replace two collinear subdivisions by the same directed segment.
         * No corner or material interface is approximated. */
        next[open_end]=open_start;incoming[open_start]=1;boundary=1;
    }
    n=0;a=start;
    do {
        if(a<0 || a>=np || seen[a] || n>=np)goto out;
        seen[a]=1;order[n++]=a;a=next[a];
    }while(a!=start);
    if(n!=nt+boundary || n<3)goto out;
    for(i=0;i<np;i++)if(next[i]>=0 && !seen[i])goto out;
    for(i=0;i<n;i++)for(j=i+2;j<n;j++) {
        if(i==0 && j==n-1)continue;
        if(segment_bad(p,order[i],order[(i+1)%n],order[j],order[(j+1)%n]))goto out;
    }
    r->triangle=malloc((size_t)(n-2)*sizeof(*r->triangle));
    if(!r->triangle){status=MESH_PATCH_NOMEM;r->reason="triangle allocation";goto out;}
    r->old_min_quality=minq;r->new_min_quality=1.0;
    {
        int proposed=patch_optimal_triangulation(p,order,n,region,r);
        if(proposed!=MESH_PATCH_OK) {
            status=proposed;r->reason=proposed==MESH_PATCH_NOMEM?"optimal triangulation allocation":"no valid constrained triangulation";
            goto out;
        }
    }
    for(i=0;i<r->count;i++)for(j=i+1;j<r->count;j++)
        if(overlap(p,r->triangle[i].vertex,r->triangle[j].vertex))goto out;
    /* Verify directed attachment faces and opposite internal owners. */
    for(i=0;i<r->count;i++)for(j=0;j<3;j++) {
        int same=0,reverse=0,u,v,ti,f;
        u=r->triangle[i].vertex[j];v=r->triangle[i].vertex[(j+1)%3];
        for(ti=0;ti<r->count;ti++)for(f=0;f<3;f++) {
            a=r->triangle[ti].vertex[f];b=r->triangle[ti].vertex[(f+1)%3];
            if(a==u && b==v)same++;
            if(a==v && b==u)reverse++;
        }
        if(same!=1 || (next[u]==v ? reverse!=0 : reverse!=1))goto out;
    }
    for(i=0;i<np;i++)if(next[i]>=0) {
        int count=0;
        for(j=0;j<r->count;j++)for(k=0;k<3;k++)
            if(r->triangle[j].vertex[k]==i &&
               r->triangle[j].vertex[(k+1)%3]==next[i])count++;
        if(count!=1)goto out;
    }
    if(r->new_min_quality<r->old_min_quality*(1.0-1e-12)) {
        status=MESH_PATCH_PINNED;r->reason="quality would deteriorate";goto out;
    }
    status=MESH_PATCH_OK;r->reason=boundary ? "validated straight-boundary proposal" : "validated interior proposal";
out:
    free(next);free(incoming);free(order);free(seen);
    if(status!=MESH_PATCH_OK)mesh_patch_release(r);
    return status;
}
