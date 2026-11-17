#include <assert.h>
#include <stdio.h>
#include "directional_etch.h"
static struct etch_ray_front flat(double depth)
{
    struct etch_ray_front f;
    f.n=1;f.s=malloc(sizeof(*f.s));
    f.s[0]=(struct etch_ray_segment){0,depth,10,depth};return f;
}
static double value(struct etch_ray_front *f,double x)
{
    int i;for(i=0;i<f->n;i++)if(x>f->s[i].x0 && x<f->s[i].x1)
        return etch_ray_y(f->s[i],x);
    assert(0);return 0;
}
static int triangle(struct etch_ray_front *f,
    double x0,double y0,double x1,double y1,double x2,double y2)
{
    struct etch_ray_segment s[3]={{x0,y0,x1,y1},{x1,y1,x2,y2},{x2,y2,x0,y0}};
    return etch_ray_opaque_polygon(f,s,3);
}
int main(void)
{
    struct etch_ray_front f=flat(10),g=flat(10);
    int i;
    /* A suspended gate stops rays even when oxide is connected around it.
     * Open windows to either side retain their full etch budget. */
    assert(!etch_ray_limit(&f,4,2,6,2));
    assert(value(&f,5)==2 && value(&f,3)==10 && value(&f,7)==10);
    /* A sloping second mask intersects the first: resolve the intersection,
     * rather than sampling or extending the mask across its endpoints. */
    assert(!etch_ray_limit(&f,3,0,7,8));
    assert(value(&f,3.5)==1 && value(&f,4.5)==2 && value(&f,6.5)==7);
    /* Order-independent lower envelope and endpoint reversal. */
    assert(!etch_ray_limit(&g,7,8,3,0));
    assert(!etch_ray_limit(&g,6,2,4,2));
    for(i=0;i<100;i++) {
        double x=(i+.37)*.1;
        assert(fabs(value(&f,x)-value(&g,x))<1e-13);
    }
    free(f.s);free(g.s);
    /* A sub-grid-width mask must not disappear; a vertical edge has no
     * projected width and cannot cast a finite shadow. */
    f=flat(10);
    assert(!etch_ray_limit(&f,5,1,5+1e-10,1));
    assert(value(&f,5+5e-11)==1);
    assert(!etch_ray_limit(&f,2,-10,2,10));
    assert(value(&f,2.01)==10);
    free(f.s);
    /* Different triangulations of the same opaque rectangle may terminate
     * rays at different interior points. Their effect on target material
     * above and below the rectangle must be identical. */
    f=flat(10);g=flat(10);
    assert(!triangle(&f,4,0,6,0,6,2));
    assert(!triangle(&f,4,0,6,2,4,2));
    assert(!triangle(&g,4,0,6,0,4,2));
    assert(!triangle(&g,6,0,6,2,4,2));
    for(i=0;i<100;i++) {
        double x=4+(i+.37)*.02;
        assert(value(&f,x)>0 && value(&f,x)<2);
        assert(value(&g,x)>0 && value(&g,x)<2);
    }
    free(f.s);free(g.s);
    /* Subdividing a straight polygon boundary must not introduce contacts
     * between the cutter and the solid's entry interface. */
    {
        struct etch_ray_segment s[6]={{4,0,5,0},{5,0,6,0},{6,0,6,2},
                                     {6,2,5,2},{5,2,4,2},{4,2,4,0}};
        f=flat(10);assert(!etch_ray_opaque_polygon(&f,s,6));
        assert(value(&f,4.5)==1 && value(&f,5.5)==1);free(f.s);
    }
    /* A nonconvex region has two occupied intervals on a ray. Its first
     * arm must shield target material in the gap between the two arms. */
    {
        struct etch_ray_segment s[8]={{2,0,8,0},{8,0,8,1},{8,1,3,1},
            {3,1,3,3},{3,3,8,3},{8,3,8,4},{8,4,2,4},{2,4,2,0}};
        f=flat(10);assert(!etch_ray_opaque_polygon(&f,s,8));
        assert(value(&f,5)==.5 && value(&f,2.5)==2 && value(&f,1)==10);
        free(f.s);
    }
    puts("PASS directional ray envelope and opaque-cell triangulation invariance");return 0;
}
