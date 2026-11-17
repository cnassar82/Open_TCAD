/* Prepare a complete single-region removal before publishing it. Straight
 * side-boundary subdivisions preserve the boundary code and exact segment.
 * Exposed, backside and multi-material points remain pinned.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "impurity.h"
#include "mesh_patch.h"

static double triangle_weight(int a,int b,int c)
{
    double x=(double)xcord(pt_nd(b))-xcord(pt_nd(a));
    double y=(double)ycord(pt_nd(b))-ycord(pt_nd(a));
    double u=(double)xcord(pt_nd(c))-xcord(pt_nd(a));
    double v=(double)ycord(pt_nd(c))-ycord(pt_nd(a));
    return (x*v-y*u)/6.0;
}
static int in_set(int v,const int *set,int n)
{
    int i;for(i=0;i<n;i++)if(set[i]==v)return 1;return 0;
}
static int chemical(int imp)
{
    return imp==As || imp==P || imp==Sb || imp==B || imp==O2 || imp==H2O ||
           imp==Au || imp==Cs || imp==Ga || (imp>=iBe && imp<=iG && (imp&1));
}
static int active_parent(int imp)
{
    switch(imp) {
    case Asa:return As;case Pa:return P;case Sba:return Sb;case Ba:return B;
    default:if(imp>=iBea && imp<=iGa && !(imp&1))return imp-1;
    }
    return -1;
}
/* Bounded least-change nodal correction. Weights include the entire support
 * of modified nodes, so changing their fields also accounts for external
 * triangles. All supports must belong to the same region in this version. */
static int transfer_density(int count,const double *oldw,const double *neww,
    const int *nodes,int stream,double *values)
{
    int i,iteration;
    double target=0,total=0,lo=HUGE_VAL,hi=-HUGE_VAL,shiftlo,shifthi,shift,sum;
    for(i=0;i<count;i++) {
        double c=sol_nd(nodes[i],stream);
        if(!isfinite(c) || c<0 || !isfinite(oldw[i]) || oldw[i]<=0)return 0;
        target+=oldw[i]*c;
        if(c<lo)lo=c;
        if(c>hi)hi=c;
        if(i>0)total+=neww[i];
    }
    if(!isfinite(target) || total<=0)return 0;
    /* A constant remains constant; exact-domain weight equality is checked
     * by the caller. Floating summation residuals are bounded below. */
    shiftlo=lo-hi;shifthi=hi-lo;
    for(iteration=0;iteration<100;iteration++) {
        shift=shiftlo+(shifthi-shiftlo)*0.5;sum=0;
        for(i=1;i<count;i++) {
            double c=sol_nd(nodes[i],stream)+shift;
            c=fmax(lo,fmin(hi,c));sum+=neww[i]*c;
        }
        if(sum<target)shiftlo=shift;else shifthi=shift;
    }
    shift=shiftlo+(shifthi-shiftlo)*0.5;sum=0;
    for(i=1;i<count;i++) {
        values[i]=fmax(lo,fmin(hi,sol_nd(nodes[i],stream)+shift));
        sum+=neww[i]*values[i];
    }
    return isfinite(sum) && fabs(sum-target)<=1e-12*fmax(fabs(target),1e-300);
}

int mesh_patch_remove_point(int ip,const char **reason)
{
    struct mesh_patch_point *points=NULL;
    struct mesh_patch_triangle *input=NULL;
    struct mesh_patch_result result;
    int *nodes=NULL,*old=NULL,**lists=NULL,*counts=NULL,*capacity=NULL;
    int *newnb=NULL,*external_t=NULL,*external_f=NULL,*newregion=NULL;
    double *oldw=NULL,*neww=NULL,*fields=NULL,*values=NULL;
    int node,nt,cap,n=1,i,j,k,t,v,r,status=MESH_PATCH_INVALID,newcount=0;
    int a,b,at,nb,regioncount=0;
    int boundary=0,boundary_count=0,boundary_code=0,boundary_end[2];
    double oldsum=0,newsum=0;
    memset(&result,0,sizeof result);*reason="invalid mesh patch";
    if(mode!=TWOD || n_imp<0 || n_imp>MAXIMP || ip<0 || ip>=np || dead_pt(ip))goto out;
    if(is_surf(ip) || is_back(ip) || num_nd(ip)!=1) {
        status=MESH_PATCH_PINNED;*reason="physical boundary/material constraint";goto out;
    }
    node=nd_pt(ip,0);
    if(node<0 || node>=nn || dead_nd(node) || mat_nd(node)==GAS)goto out;
    nt=num_tri_nd(node);
    if(nt<2 || nt>INT_MAX/2-1)goto out;
    cap=2*nt+1;
    nodes=malloc((size_t)cap*sizeof(int));old=malloc((size_t)nt*sizeof(int));
    points=malloc((size_t)cap*sizeof(*points));input=malloc((size_t)nt*sizeof(*input));
    if(!nodes || !old || !points || !input){status=MESH_PATCH_NOMEM;goto out;}
    nodes[0]=node;r=-1;
    for(i=0;i<nt;i++) {
        t=tri_nd(node,i);old[i]=t;
        if(t<0 || t>=ne || dead_tri(t))goto out;
        if(!leaf(tri[t]) || !root(tri[t])) {
            status=MESH_PATCH_PINNED;*reason="refinement hierarchy constraint";goto out;
        }
        if(i && in_set(t,old,i))goto out;
        if(reg_tri(t)<0 || reg_tri(t)>=nreg || !reg[reg_tri(t)])goto out;
        if(r<0)r=reg_tri(t);
        input[i].region=reg_tri(t);
        for(j=0;j<3;j++) {
            v=vert_tri(t,j);
            if(v<0 || v>=nn || dead_nd(v) || mat_nd(v)!=mat_nd(node))goto out;
            for(k=0;k<n && nodes[k]!=v;k++);
            if(k==n)nodes[n++]=v;
            input[i].vertex[j]=k;
        }
    }
    for(i=0;i<n;i++) {
        points[i].xy[0]=xcord(pt_nd(nodes[i]));points[i].xy[1]=ycord(pt_nd(nodes[i]));
        for(j=0;j<i;j++)if(pt_nd(nodes[i])==pt_nd(nodes[j]))goto out;
    }
    status=mesh_patch_propose(0,n,points,nt,input,&result);*reason=result.reason;
    if(status!=MESH_PATCH_OK)goto out;
    /* From here every unsuccessful preparation must still leave globals intact. */
    status=MESH_PATCH_INVALID;*reason="patch attachment or support invalid";
    newcount=result.count;
    boundary=newcount==nt-1;
    if(!boundary && newcount!=nt-2)goto out;
    /* Validate boundary-condition identity before allocations or publication.
     * The geometry kernel has proved that the two segments are collinear. */
    for(i=0;i<nt;i++)for(j=0;j<3;j++) {
        a=input[i].vertex[(j+1)%3];b=input[i].vertex[(j+2)%3];
        nb=tri[old[i]]->nb[j];
        if((a!=0 && b!=0) || nb>=0)continue;
        if(!boundary || boundary_count>=2)goto out;
        if(nb==EXPOSED || (boundary_count && nb!=boundary_code)) {
            status=MESH_PATCH_PINNED;*reason="boundary condition constraint";goto out;
        }
        boundary_code=nb;boundary_end[boundary_count++]=a==0?b:a;
    }
    if(boundary && (boundary_count!=2 || boundary_end[0]==boundary_end[1]))goto out;
    oldw=calloc(n,sizeof(double));neww=calloc(n,sizeof(double));
    fields=malloc((size_t)n*MAXIMP*sizeof(double));values=malloc((size_t)n*sizeof(double));
    lists=calloc(n,sizeof(int *));counts=calloc(n,sizeof(int));capacity=calloc(n,sizeof(int));
    newnb=malloc((size_t)newcount*3*sizeof(int));
    external_t=malloc((size_t)nt*sizeof(int));external_f=malloc((size_t)nt*sizeof(int));
    newregion=malloc((size_t)reg[r]->tri.num*sizeof(int));
    if(!oldw || !neww || !fields || !values || !lists || !counts || !capacity ||
       !newnb || !external_t || !external_f || !newregion){status=MESH_PATCH_NOMEM;goto out;}
    for(i=0;i<n;i++) {
        int ndx=nodes[i],outside=0;
        for(j=0;j<num_tri_nd(ndx);j++) {
            t=tri_nd(ndx,j);
            if(t<0 || t>=ne || dead_tri(t) || !leaf(tri[t]))goto out;
            if(reg_tri(t)!=r){status=MESH_PATCH_PINNED;*reason="field support spans region boundary";goto out;}
            for(k=0;k<j;k++)if(tri_nd(ndx,k)==t)goto out;
            if(!in_set(ndx,tri[t]->nd,3))goto out;
            oldw[i]+=triangle_weight(tri[t]->nd[0],tri[t]->nd[1],tri[t]->nd[2]);
            if(!in_set(t,old,nt))outside++;
        }
        if(oldw[i]<=0 || !isfinite(oldw[i]))goto out;
        neww[i]=oldw[i];
        capacity[i]=outside+newcount;
        if(capacity[i]>0)lists[i]=malloc((size_t)capacity[i]*sizeof(int));
        if(capacity[i]>0 && !lists[i]){status=MESH_PATCH_NOMEM;goto out;}
        for(j=0;j<num_tri_nd(ndx);j++)if(!in_set(tri_nd(ndx,j),old,nt))
            lists[i][counts[i]++]=tri_nd(ndx,j);
        for(j=0;j<n_imp;j++) {
            if(!isfinite(sol_nd(ndx,j)))goto out;
            fields[i*MAXIMP+j]=sol_nd(ndx,j);
        }
    }
    for(i=0;i<nt;i++) {
        double w=triangle_weight(tri[old[i]]->nd[0],tri[old[i]]->nd[1],tri[old[i]]->nd[2]);
        for(j=0;j<3;j++)neww[input[i].vertex[j]]-=w;
    }
    for(i=0;i<newcount;i++) {
        const int *vv=result.triangle[i].vertex;
        double w=triangle_weight(nodes[vv[0]],nodes[vv[1]],nodes[vv[2]]);
        if(w<=0)goto out;
        for(j=0;j<3;j++) {
            v=vv[j];neww[v]+=w;lists[v][counts[v]++]=old[i];
        }
    }
    if(counts[0]!=0 || fabs(neww[0])>1e-12*oldw[0])goto out;
    for(i=0;i<n;i++)oldsum+=oldw[i];
    for(i=1;i<n;i++) {
        if(neww[i]<=0 || !isfinite(neww[i]))goto out;
        newsum+=neww[i];
    }
    if(fabs(oldsum-newsum)>1e-12*oldsum)goto out;
    for(j=0;j<n_imp;j++) {
        int imp=soltoimp[j],parent=active_parent(imp);
        if(chemical(imp) || imp==V || imp==I || imp==T) {
            if(!transfer_density(n,oldw,neww,nodes,j,values)) {
                status=MESH_PATCH_PINNED;*reason="bounded conservative transfer unavailable";goto out;
            }
            for(i=1;i<n;i++)fields[i*MAXIMP+j]=values[i];
        }else if(parent<0 && imp!=N && imp!=H && imp!=Psi && imp!=XVEL && imp!=YVEL &&
                 imp!=GRN && imp!=DELA && imp!=Sxx && imp!=Syy && imp!=Sxy) {
            status=MESH_PATCH_PINNED;*reason="unsupported field";goto out;
        }
    }
    /* Active species retain their local active fraction of the remapped total;
     * they are subsequently recomputed by the physical activation model. */
    for(j=0;j<n_imp;j++) {
        int parent=active_parent(soltoimp[j]),stream;
        if(parent<0)continue;
        stream=imptosol[parent];
        if(stream<0 || stream>=n_imp)goto out;
        for(i=1;i<n;i++) {
            double total=sol_nd(nodes[i],stream),active=sol_nd(nodes[i],j),fraction;
            if(!isfinite(total) || !isfinite(active) || total<0 || active<0 || active>total)goto out;
            fraction=total>0 ? active/total:0;
            fields[i*MAXIMP+j]=fraction*fields[i*MAXIMP+stream];
        }
    }
    /* Prepare all neighbor attachments, including reciprocal outside slots. */
    for(i=0;i<nt;i++) {
        at=-1;
        for(j=0;j<3;j++)if(input[i].vertex[j]==0)at=j;
        if(at<0)goto out;
        a=tri[old[i]]->nd[(at+1)%3];b=tri[old[i]]->nd[(at+2)%3];
        nb=tri[old[i]]->nb[at];external_t[i]=nb;external_f[i]=-1;
        if(nb>=0) {
            if(nb>=ne || dead_tri(nb) || !leaf(tri[nb]) || in_set(nb,old,nt))goto out;
            for(j=0;j<3;j++) {
                int x=pt_nd(tri[nb]->nd[(j+1)%3]),y=pt_nd(tri[nb]->nd[(j+2)%3]);
                if((x==pt_nd(a) && y==pt_nd(b)) || (x==pt_nd(b) && y==pt_nd(a)))external_f[i]=j;
            }
            if(external_f[i]<0 || tri[nb]->nb[external_f[i]]!=old[i])goto out;
        }
    }
    for(i=0;i<newcount;i++)for(j=0;j<3;j++) {
        int found=0;
        a=result.triangle[i].vertex[(j+1)%3];b=result.triangle[i].vertex[(j+2)%3];
        for(k=0;k<newcount;k++)if(k!=i)for(v=0;v<3;v++)
            if(result.triangle[k].vertex[(v+1)%3]==b && result.triangle[k].vertex[(v+2)%3]==a) {
                newnb[i*3+j]=old[k];found++;
            }
        if(!found)for(k=0;k<nt;k++) {
            int at0=-1;
            for(v=0;v<3;v++)if(input[k].vertex[v]==0)at0=v;
            if(input[k].vertex[(at0+1)%3]==a && input[k].vertex[(at0+2)%3]==b) {
                newnb[i*3+j]=external_t[k];found++;
            }
        }
        if(!found && boundary &&
            ((a==boundary_end[0] && b==boundary_end[1]) ||
             (b==boundary_end[0] && a==boundary_end[1]))) {
            newnb[i*3+j]=boundary_code;found=1;
        }
        if(found!=1)goto out;
    }
    for(i=0;i<reg[r]->tri.num;i++) {
        t=reg[r]->tri.list[i];
        if(!in_set(t,old,nt))newregion[regioncount++]=t;
    }
    if(regioncount!=reg[r]->tri.num-nt)goto out;
    for(i=0;i<newcount;i++)newregion[regioncount++]=old[i];
    /* Commit: only assignments, retirement flags and freeing old list storage.
     * Every allocation and validation has completed above. */
    for(i=0;i<newcount;i++) {
        for(j=0;j<3;j++) {
            tri[old[i]]->nd[j]=nodes[result.triangle[i].vertex[j]];
            tri[old[i]]->nb[j]=newnb[i*3+j];
        }
        clr(tri[old[i]],GEOMDN|CLKWS);set(tri[old[i]],NEIGHB);
    }
    for(i=0;i<newcount;i++)for(j=0;j<3;j++) {
        nb=newnb[i*3+j];if(nb<0 || in_set(nb,old,nt))continue;
        a=pt_nd(tri[old[i]]->nd[(j+1)%3]);b=pt_nd(tri[old[i]]->nd[(j+2)%3]);
        for(k=0;k<3;k++) {
            int x=pt_nd(tri[nb]->nd[(k+1)%3]),y=pt_nd(tri[nb]->nd[(k+2)%3]);
            if((x==a && y==b) || (x==b && y==a))tri[nb]->nb[k]=old[i];
        }
    }
    for(i=newcount;i<nt;i++)set(tri[old[i]],KILL_ME);
    for(i=0;i<n;i++) {
        free(nd[nodes[i]]->ltri.list);
        nd[nodes[i]]->ltri.list=lists[i];lists[i]=NULL;
        nd[nodes[i]]->ltri.num=counts[i];nd[nodes[i]]->ltri.all=capacity[i];
        for(j=0;j<n_imp;j++)sol_nd(nodes[i],j)=fields[i*MAXIMP+j];
    }
    free(reg[r]->tri.list);reg[r]->tri.list=newregion;newregion=NULL;
    reg[r]->tri.num=regioncount;reg[r]->tri.all=regioncount;
    set(nd[node],KILL_ME);set(pt[ip],KILL_ME);pt[ip]->nn=0;
    geom_dirty=clkws_dirty=neigh_dirty=edge_dirty=nd2edg_dirty=reg_dirty=need_waste=TRUE;
    status=MESH_PATCH_OK;*reason="committed conservative interior removal";
out:
    if(lists)for(i=0;i<n;i++)free(lists[i]);
    free(lists);free(counts);free(capacity);free(nodes);free(old);free(points);free(input);
    free(oldw);free(neww);free(fields);free(values);free(newnb);free(external_t);free(external_f);free(newregion);
    mesh_patch_release(&result);
    return status;
}
