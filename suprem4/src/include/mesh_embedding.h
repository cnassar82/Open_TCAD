#ifndef SUPREM_MESH_EMBEDDING_H
#define SUPREM_MESH_EMBEDDING_H
/* Broad-phase BVH plus exact stored-coordinate predicates. Positive local
 * Jacobians do not exclude overlap between non-neighboring elements. */
#include <stdlib.h>
#include <math.h>
#include <limits.h>
struct embedding_triangle { int vertex[3], boundary, identity; };
struct embedding_box { double lo[2],hi[2],center[2]; int triangle; };
struct embedding_tree { double lo[2],hi[2]; int first,count,left,right; };
struct embedding_work {
    const float (*xy)[3];
    const struct embedding_triangle *triangle;
    struct embedding_box *box;
    struct embedding_tree *tree;
    int used,bad;
};
extern int orient2d_float_sign(const float *,const float *,const float *);
static int embedding_sort_axis;
static int embedding_compare(const void *aa,const void *bb)
{
    const struct embedding_box *a=aa,*b=bb;
    double d=a->center[embedding_sort_axis]-b->center[embedding_sort_axis];
    if(d)return d<0?-1:1;
    return (a->triangle>b->triangle)-(a->triangle<b->triangle);
}
static int embedding_build(struct embedding_work *w,int first,int count)
{
    int node=w->used++,i,k,half;
    struct embedding_tree *t=&w->tree[node];
    t->first=first;t->count=count;t->left=t->right=-1;
    for(k=0;k<2;k++){t->lo[k]=HUGE_VAL;t->hi[k]=-HUGE_VAL;}
    for(i=first;i<first+count;i++)for(k=0;k<2;k++) {
        if(w->box[i].lo[k]<t->lo[k])t->lo[k]=w->box[i].lo[k];
        if(w->box[i].hi[k]>t->hi[k])t->hi[k]=w->box[i].hi[k];
    }
    if(count>1) {
        embedding_sort_axis=t->hi[1]-t->lo[1]>t->hi[0]-t->lo[0];
        qsort(w->box+first,count,sizeof(*w->box),embedding_compare);
        half=count/2;
        t->left=embedding_build(w,first,half);
        t->right=embedding_build(w,first+half,count-half);
    }
    return node;
}
static int embedding_separated(struct embedding_work *w,const int *a,const int *b)
{
    int i,j;
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++)if(orient2d_float_sign(w->xy[a[i]],w->xy[a[(i+1)%3]],w->xy[b[j]])>0)break;
        if(j==3)return 1;
    }
    return 0;
}
static int embedding_on_segment(struct embedding_work *w,int p,int a,int b)
{
    int k;
    if(orient2d_float_sign(w->xy[a],w->xy[b],w->xy[p]))return 0;
    for(k=0;k<2;k++)if(w->xy[p][k]<fmin(w->xy[a][k],w->xy[b][k]) ||
                         w->xy[p][k]>fmax(w->xy[a][k],w->xy[b][k]))return 0;
    return 1;
}
static int embedding_boundary_contact(struct embedding_work *w,int a,int b,int c,int d)
{
    int ac=orient2d_float_sign(w->xy[a],w->xy[b],w->xy[c]);
    int ad=orient2d_float_sign(w->xy[a],w->xy[b],w->xy[d]);
    int ca=orient2d_float_sign(w->xy[c],w->xy[d],w->xy[a]);
    int cb=orient2d_float_sign(w->xy[c],w->xy[d],w->xy[b]);
    int axis;
    if(ac*ad<0 && ca*cb<0)return 1;
    if(!ac && !ad && !ca && !cb) {
        axis=w->xy[a][0]!=w->xy[b][0]?0:1;
        if(fmin(fmax(w->xy[a][axis],w->xy[b][axis]),fmax(w->xy[c][axis],w->xy[d][axis]))>
           fmax(fmin(w->xy[a][axis],w->xy[b][axis]),fmin(w->xy[c][axis],w->xy[d][axis])))return 1;
    }
    return (a!=c && a!=d && embedding_on_segment(w,a,c,d)) ||
           (b!=c && b!=d && embedding_on_segment(w,b,c,d)) ||
           (c!=a && c!=b && embedding_on_segment(w,c,a,b)) ||
           (d!=a && d!=b && embedding_on_segment(w,d,a,b));
}
static int embedding_pair(struct embedding_work *w,int aa,int bb)
{
    const struct embedding_triangle *a=&w->triangle[aa],*b=&w->triangle[bb];
    int i,j;
    if(!embedding_separated(w,a->vertex,b->vertex) &&
       !embedding_separated(w,b->vertex,a->vertex))goto invalid;
    for(i=0;i<3;i++)for(j=0;j<3;j++)if((a->boundary&(1<<i)) || (b->boundary&(1<<j))) {
        if(embedding_boundary_contact(w,a->vertex[i],a->vertex[(i+1)%3],
                                      b->vertex[j],b->vertex[(j+1)%3]))goto invalid;
    }
    return 1;
invalid:
    w->bad=a->identity;return 0;
}
static int embedding_visit(struct embedding_work *w,int a,int b)
{
    struct embedding_tree *x=&w->tree[a],*y=&w->tree[b];
    int k;
    if(a==b) {
        if(x->count==1)return 1;
        return embedding_visit(w,x->left,x->left) &&
               embedding_visit(w,x->left,x->right) &&
               embedding_visit(w,x->right,x->right);
    }
    for(k=0;k<2;k++)if(x->hi[k]<y->lo[k] || y->hi[k]<x->lo[k])return 1;
    if(x->count==1 && y->count==1)
        return embedding_pair(w,w->box[x->first].triangle,w->box[y->first].triangle);
    if(x->count>=y->count)
        return embedding_visit(w,x->left,b) && embedding_visit(w,x->right,b);
    return embedding_visit(w,a,y->left) && embedding_visit(w,a,y->right);
}
static int mesh_embedding_check(const float (*xy)[3],int count,
                              const struct embedding_triangle *triangles,int *bad)
{
    struct embedding_work w;
    int i,j,k,p,valid;
    *bad=-1;
    if(count<0 || count>INT_MAX/2)return 0;
    if(!count)return 1;
    w.xy=xy;w.triangle=triangles;w.used=0;w.bad=-1;
    w.box=calloc((size_t)count,sizeof(*w.box));
    w.tree=calloc((size_t)count*2,sizeof(*w.tree));
    if(!w.box || !w.tree){free(w.box);free(w.tree);return 0;}
    for(i=0;i<count;i++) {
        if(orient2d_float_sign(xy[triangles[i].vertex[0]],xy[triangles[i].vertex[1]],
                              xy[triangles[i].vertex[2]])<=0) {
            *bad=triangles[i].identity;free(w.box);free(w.tree);return 0;
        }
        w.box[i].triangle=i;
        for(k=0;k<2;k++){w.box[i].lo[k]=HUGE_VAL;w.box[i].hi[k]=-HUGE_VAL;}
        for(j=0;j<3;j++) {
            p=triangles[i].vertex[j];
            for(k=0;k<2;k++) {
                if(xy[p][k]<w.box[i].lo[k])w.box[i].lo[k]=xy[p][k];
                if(xy[p][k]>w.box[i].hi[k])w.box[i].hi[k]=xy[p][k];
            }
        }
        for(k=0;k<2;k++)w.box[i].center[k]=.5*(w.box[i].lo[k]+w.box[i].hi[k]);
    }
    embedding_build(&w,0,count);valid=embedding_visit(&w,0,0);*bad=w.bad;
    free(w.box);free(w.tree);return valid;
}
#endif
