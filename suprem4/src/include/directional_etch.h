#ifndef SUPREM_DIRECTIONAL_ETCH_H
#define SUPREM_DIRECTIONAL_ETCH_H
#include <stdlib.h>
#include <math.h>
/* Vertical rays travel toward increasing y.  The front is a piecewise-linear
 * lower envelope of the depth-limited free surface and every opaque material.
 * Separate segments retain jumps at mask edges; no mesh-size tolerance or
 * sampled rays can erase a thin mask or a buried film. */
struct etch_ray_segment { double x0,y0,x1,y1; };
struct etch_ray_front { struct etch_ray_segment *s; int n; };
static double etch_ray_y(struct etch_ray_segment s,double x)
{
    if(x==s.x0)return s.y0;
    if(x==s.x1)return s.y1;
    return s.y0+(s.y1-s.y0)*((x-s.x0)/(s.x1-s.x0));
}
/* Overlay one finite blocker (or depth-limit segment) on the current front.
 * Split at both projected endpoints and at an interior line intersection. */
static int etch_ray_limit(struct etch_ray_front *front,
                         double x0,double y0,double x1,double y1)
{
    struct etch_ray_segment b,*out,a,c;
    int i,j,k,m,n=0;
    double cuts[5],v,t,l,r,dl,dr,yl,yr;
    if(x0==x1)return 0; /* zero projected area */
    if(x0>x1){t=x0;x0=x1;x1=t;t=y0;y0=y1;y1=t;}
    b.x0=x0;b.y0=y0;b.x1=x1;b.y1=y1;
    out=malloc((size_t)(front->n*4)*sizeof(*out));
    if(!out)return -1;
    for(i=0;i<front->n;i++) {
        a=front->s[i];k=0;cuts[k++]=a.x0;cuts[k++]=a.x1;
        l=fmax(a.x0,x0);r=fmin(a.x1,x1);
        if(l>=r){out[n++]=a;continue;}
        dl=etch_ray_y(a,l)-etch_ray_y(b,l);
        dr=etch_ray_y(a,r)-etch_ray_y(b,r);
        if(dl<=0 && dr<=0){out[n++]=a;continue;}
        if(l<r) {
            if(l>a.x0)cuts[k++]=l;
            if(r<a.x1)cuts[k++]=r;
            dl=etch_ray_y(a,l)-etch_ray_y(b,l);
            dr=etch_ray_y(a,r)-etch_ray_y(b,r);
            if((dl<0 && dr>0)||(dl>0 && dr<0))
                cuts[k++]=l+(r-l)*(dl/(dl-dr));
        }
        for(j=1;j<k;j++){v=cuts[j];for(m=j;m>0 && cuts[m-1]>v;m--)cuts[m]=cuts[m-1];cuts[m]=v;}
        for(j=1;j<k;j++) {
            if(cuts[j]<=cuts[j-1])continue;
            c.x0=cuts[j-1];c.x1=cuts[j];
            yl=etch_ray_y(a,c.x0);yr=etch_ray_y(a,c.x1);
            if(c.x0>=x0 && c.x1<=x1) {
                yl=fmin(yl,etch_ray_y(b,c.x0));
                yr=fmin(yr,etch_ray_y(b,c.x1));
            }
            c.y0=yl;c.y1=yr;out[n++]=c;
        }
    }
    free(front->s);front->s=out;front->n=n;
    return 0;
}
/* Terminate rays strictly inside an opaque region wherever its vertical chord
 * has positive length. Entry and exit shadows remove the same target volume,
 * but putting the cutter on either material interface needlessly creates
 * coincident Boolean edges. Mid-chords are an exact alternative, without an
 * arbitrary geometric displacement. Use the polygon's vertical decomposition,
 * not individual mesh triangles: internal vertices must not imprint zero-width
 * contacts with the target interface on an otherwise smooth cutter. */
struct etch_ray_crossing { double left,right,middle; };
static int etch_ray_double_order(const void *a,const void *b)
{double x=*(const double *)a,y=*(const double *)b;return (x>y)-(x<y);}
static int etch_ray_crossing_order(const void *a,const void *b)
{double x=((const struct etch_ray_crossing *)a)->middle,y=((const struct etch_ray_crossing *)b)->middle;return (x>y)-(x<y);}
static int etch_ray_opaque_polygon(struct etch_ray_front *front,
    const struct etch_ray_segment *boundary,int count)
{
    double *x,l,r,m;
    struct etch_ray_crossing *cross;
    int i,j,n,k,rc=0;
    x=malloc((size_t)(2*count)*sizeof(*x));
    cross=malloc((size_t)count*sizeof(*cross));
    if(!x || !cross){free(x);free(cross);return -1;}
    for(i=0;i<count;i++){x[2*i]=boundary[i].x0;x[2*i+1]=boundary[i].x1;}
    qsort(x,2*count,sizeof(*x),etch_ray_double_order);
    n=0;for(i=0;i<2*count;i++)if(!n || x[i]!=x[n-1])x[n++]=x[i];
    for(i=1;i<n;i++) {
        l=x[i-1];r=x[i];m=l+(r-l)*.5;k=0;
        for(j=0;j<count;j++) {
            struct etch_ray_segment s=boundary[j];
            if(s.x0==s.x1)continue;
            if(fmin(s.x0,s.x1)>l || fmax(s.x0,s.x1)<r)continue;
            cross[k].left=etch_ray_y(s,l);cross[k].right=etch_ray_y(s,r);
            cross[k++].middle=etch_ray_y(s,m);
        }
        if(k%2){rc=-2;break;} /* Invalid polygon, never guess its interior. */
        qsort(cross,k,sizeof(*cross),etch_ray_crossing_order);
        for(j=0;j<k;j+=2)
            if(etch_ray_limit(front,l,(cross[j].left+cross[j+1].left)*.5,
                                   r,(cross[j].right+cross[j+1].right)*.5)) {rc=-1;break;}
        if(rc)break;
    }
    free(x);free(cross);return rc;
}
#endif
