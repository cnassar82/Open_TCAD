/* Quotient physically equivalent regions using shared oriented incidence.
 * The same-domain operation preserves every triangle, coordinate and field;
 * it removes region constraints only when the merged legacy face is a disk.
 */
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "material_domain.h"
struct domain_edge {int a,b,from,to,region;};
struct domain_pair {int a,b;};
static int edge_order(const void *aa,const void *bb)
{
    const struct domain_edge *a=aa,*b=bb;
    if(a->a!=b->a)return a->a<b->a?-1:1;
    return a->b<b->b?-1:a->b>b->b;
}
static int pair_order(const void *aa,const void *bb)
{
    const struct domain_pair *a=aa,*b=bb;
    if(a->a!=b->a)return a->a<b->a?-1:1;
    return a->b<b->b?-1:a->b>b->b;
}
static int root_region(const int *parent,int r)
{while(parent[r]!=r)r=parent[r];return r;}
static int inside(const int *parent,int r,int a,int b)
{r=root_region(parent,r);return r==a || r==b;}
static int disk(int np,int count,const struct domain_edge *edges,const int *parent,
                int a,int b,int *next,int *incoming)
{
    int i,j,k,n=0,start=-1,v,steps,owner;
    for(i=0;i<np;i++)next[i]=-1,incoming[i]=0;
    for(i=0;i<count;i=j) {
        for(j=i+1;j<count && edges[j].a==edges[i].a && edges[j].b==edges[i].b;j++);
        owner=-1;k=0;
        if(inside(parent,edges[i].region,a,b)){owner=i;k++;}
        if(j==i+2 && inside(parent,edges[i+1].region,a,b)){owner=i+1;k++;}
        if(k!=1)continue;
        v=edges[owner].from;
        if(next[v]>=0 || incoming[edges[owner].to]++)return 0;
        next[v]=edges[owner].to;start=v;n++;
    }
    if(n<3)return 0;
    for(i=0;i<np;i++)if((next[i]>=0)!=(incoming[i]>0))return 0;
    v=start;
    for(steps=0;steps<n;steps++) {
        if(v<0 || next[v]<0)return 0;
        v=next[v];
        if(v==start)return steps+1==n;
    }
    return 0;
}
int material_domain_plan(int np,int nr,const struct domain_region *regions,
    int nt,const struct domain_triangle *triangles,int *map,int *domains)
{
    struct domain_edge *edges=NULL;
    struct domain_pair *pairs=NULL;
    int *parent=NULL,*next=NULL,*incoming=NULL,*result=NULL;
    int i,j,k,a,b,n=0,count,changed,nd=0,status=-1;
    if(np<3 || nr<1 || nt<1 || nt>INT_MAX/3 || !regions || !triangles || !map || !domains)return -1;
    count=3*nt;
    if((size_t)count>(size_t)-1/sizeof(*edges))return -1;
    edges=malloc((size_t)count*sizeof(*edges));pairs=malloc((size_t)count*sizeof(*pairs));
    parent=malloc((size_t)nr*sizeof(int));result=malloc((size_t)nr*sizeof(int));
    next=malloc((size_t)np*sizeof(int));incoming=malloc((size_t)np*sizeof(int));
    if(!edges || !pairs || !parent || !result || !next || !incoming){status=-2;goto out;}
    for(i=0;i<nr;i++)parent[i]=i,result[i]=-1;
    for(i=0;i<nt;i++) {
        const struct domain_triangle *t=&triangles[i];
        if(t->region<0 || t->region>=nr)goto out;
        for(j=0;j<3;j++) {
            a=t->vertex[j];b=t->vertex[(j+1)%3];
            if(a<0 || a>=np || b<0 || b>=np || a==b)goto out;
            edges[3*i+j].a=a<b?a:b;edges[3*i+j].b=a<b?b:a;
            edges[3*i+j].from=a;edges[3*i+j].to=b;edges[3*i+j].region=t->region;
        }
    }
    qsort(edges,count,sizeof(*edges),edge_order);
    for(i=0;i<count;i=j) {
        for(j=i+1;j<count && edges[j].a==edges[i].a && edges[j].b==edges[i].b;j++);
        if(j-i>2)goto out;
        if(j-i==2) {
            if(edges[i].from!=edges[i+1].to || edges[i].to!=edges[i+1].from)goto out;
            a=edges[i].region;b=edges[i+1].region;
            if(a==b || regions[a].protected_region || regions[b].protected_region ||
               regions[a].material!=regions[b].material || regions[a].model_key!=regions[b].model_key)continue;
            pairs[n].a=a<b?a:b;pairs[n++].b=a<b?b:a;
        }
    }
    qsort(pairs,n,sizeof(*pairs),pair_order);
    for(i=0,k=0;i<n;i++)if(k==0 || pair_order(&pairs[i],&pairs[k-1]))pairs[k++]=pairs[i];
    n=k;
    do {
        changed=0;
        for(i=0;i<n;i++) {
            a=root_region(parent,pairs[i].a);b=root_region(parent,pairs[i].b);
            if(a==b)continue;
            if(disk(np,count,edges,parent,a,b,next,incoming)) {
                if(a<b)parent[b]=a;else parent[a]=b;
                changed=1;
            }
        }
    }while(changed);
    for(i=0;i<nr;i++) {
        a=root_region(parent,i);
        if(result[a]<0)result[a]=nd++;
        result[i]=result[a];
    }
    memcpy(map,result,(size_t)nr*sizeof(int));*domains=nd;status=0;
out:
    free(edges);free(pairs);free(parent);free(result);free(next);free(incoming);
    return status;
}
