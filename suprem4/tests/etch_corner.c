/* Local oxide edges captured from the doubled-density PMOS failure,
 * clipped by the supporting line of its first etch segment.
 * Use real point/node/edge insertion, material-interface splitting, and
 * skeleton validation. No clipping or topology routines are mocked. */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "skel.h"
#include "material.h"
#include "impurity.h"
#include "diffuse.h"
#include "matrix.h"
#include "regrid.h"

void panic(char *message) { fprintf(stderr,"Unexpected panic: %s\n",message); abort(); }
/* The junction-collapse fixture must use rem_edg_end, not vertex remeshing. */
int kill_node(int node) { assert(0); return -1; }

static double area(float *a, float *b, float *c)
{
    return fabs(((double)b[0]-a[0])*((double)c[1]-a[1])-
                ((double)b[1]-a[1])*((double)c[0]-a[0]))/2;
}

static void crossing(int edge, float *c)
{
    float p[3], dp[3], q[3]={.00798112992,4.37345952e-5,0};
    float end[3]={.00800000038,4.21760960e-5,0}, dq[3], alpha[3];
    int i;
    for(i=0;i<2;i++) {
        p[i]=cordinate(pt_edg(edge,0),i);
        dp[i]=cordinate(pt_edg(edge,1),i)-p[i];
        dq[i]=end[i]-q[i];
    }
    assert(lil(p,dp,q,dq,alpha)==0);
    assert(alpha[0]>=0 && alpha[0]<=1);
    for(i=0;i<2;i++) c[i]=p[i]+alpha[0]*dp[i];
    c[2]=0;
}

static void check(int reverse)
{
    float coordinates[3][3]={
        {.00799527951,4.93631014e-5,0},
        {.00800000038,8.85594000e-6,0},
        {.00792452972,3.65687010e-5,0}};
    float cuts[2][3], legacy[3], current[3];
    int points[3], nodes[3], silicon[3], edges[2], kept[2];
    int i,j,p,other,split,region,closure,endpoint;
    double expected;
    for(i=0;i<3;i++) {
        points[i]=mk_pt(2,coordinates[i]);
        silicon[i]=mk_nd(points[i],Si);
        nodes[i]=mk_nd(points[i],SiO2);
        nd[nodes[i]]->sol[0]=nd[silicon[i]]->sol[0]=coordinates[i][0]*1e6;
    }
    for(i=0;i<2;i++) {
        edges[i]=mk_edg(nodes[reverse ? i+1 : 0],nodes[reverse ? 0 : i+1]);
        mk_edg(silicon[reverse ? i+1 : 0],silicon[reverse ? 0 : i+1]);
        crossing(edges[i],cuts[i]);
    }
    expected=area(coordinates[0],cuts[0],cuts[1]);
    assert(expected>0);
    /* The legacy operation moves the OUT corner onto the first cut. */
    memcpy(legacy,pt[points[0]]->cord,sizeof(legacy));
    assert(sp_edge(edges[0],cuts[0],&p,TRUE)<0 && p==points[0]);
    assert(memcmp(pt[p]->cord,legacy,2*sizeof(float))!=0);
    memcpy(pt[p]->cord,legacy,sizeof(legacy));
    edg[edges[0]]->len=l_edge(edges[0]);

    for(i=0;i<2;i++) {
        crossing(edges[i],current);
        assert(memcmp(current,cuts[i],2*sizeof(float))==0);
        split=sp_edge_keep(edges[i],current,&p,TRUE,nodes[0]);
        assert(split>=0 && p!=points[0]);
        assert(memcmp(pt[points[0]]->cord,coordinates[0],2*sizeof(float))==0);
        assert(num_nd(p)==2); /* Both material edges share the new point. */
        for(j=0;j<num_nd(p);j++)
            assert(fabs(nd[nd_pt(p,j)]->sol[0]-pt[p]->cord[0]*1e6)<.001);
        kept[i]=(nd_edg(edges[i],0)==nodes[0] || nd_edg(edges[i],1)==nodes[0]) ? edges[i] : split;
        endpoint=nd_edg(kept[i],0)==nodes[0] ? 1 : 0;
        assert(pt_edg(kept[i],endpoint)==p);
    }
    endpoint=nd_edg(kept[0],0)==nodes[0] ? 1 : 0;
    p=pt_edg(kept[0],endpoint);
    other=nd_edg(kept[0],endpoint);
    endpoint=nd_edg(kept[1],0)==nodes[0] ? 1 : 0;
    assert(fabs(area(pt[points[0]]->cord,pt[p]->cord,
                    pt[pt_edg(kept[1],endpoint)]->cord)-expected)<expected*1e-12);
    closure=mk_edg(other,nd_edg(kept[1],endpoint));
    region=cr_sreg(SiO2);
    ad_edge(region,kept[0],sreg[region]->bnd,MAYBE,BEFORE);
    ad_edge(region,closure,sreg[region]->bnd,MAYBE,BEFORE);
    ad_edge(region,kept[1],sreg[region]->bnd,MAYBE,BEFORE);
    ck_clock(region,TRUE);
    assert(sreg[region]->len==3);
    free_skel(region);
}

static void check_minimum_spacing(void)
{
    int i, a,b,e,p,before;
    float left[3]={0,0,0}, right[3]={1e-4,0,0}, cut[3]={1e-9,0,0};
    for(i=0;i<2;i++) {
        a=mk_nd(mk_pt(2,left),SiO2); b=mk_nd(mk_pt(2,right),SiO2);
        e=mk_edg(a,b);
        if(i) set(edg[e],REGS);
        before=np;
        assert(sp_edge_keep(e,left,&p,FALSE,a)==-2);
        assert(np==before && p==pt_nd(a));
        /* Even a cut below MIN_GRID_SPACE must not move a retained point. */
        assert(sp_edge_keep(e,cut,&p,FALSE,a)>=0);
        assert(np==before+1 && p!=pt_nd(a));
        assert(pt[pt_nd(a)]->cord[0]==0 && pt[p]->cord[0]==cut[0]);
    }
}

static void check_shared_crossing_endpoint(void)
{
    float a[3]={0,0,0}, b[3]={1e-4,0,0}, c[3]={9.9e-5,0,0};
    int na=mk_nd(mk_pt(2,a),SiO2), nb=mk_nd(mk_pt(2,b),SiO2);
    int e=mk_edg(na,nb), p;
    /* The opposite endpoint can also be a previously recorded ON
     * junction. Neither its coordinate nor its node identity may move. */
    assert(sp_edge_keep(e,c,&p,TRUE,na)>=0);
    assert(pt[pt_nd(na)]->cord[0]==a[0]);
    assert(pt[pt_nd(nb)]->cord[0]==b[0]);
    assert(p!=pt_nd(na) && p!=pt_nd(nb));
    assert(pt[p]->cord[0]==c[0]);
}

static void check_refinement_endpoint_status(void)
{
    float a[3]={0,2e-4,0},b[3]={3.3186911e-8,2e-4,0};
    float before_a[3],before_b[3];
    int pa=mk_pt(2,a),pb=mk_pt(2,b),edge;
    edge=mk_edg(mk_nd(pa,SiO2),mk_nd(pb,SiO2));
    memcpy(before_a,pt[pa]->cord,sizeof(before_a));
    memcpy(before_b,pt[pb]->cord,sizeof(before_b));
    /* Requests smaller than the minimum spacing return -1/-2.  They
     * must neither move the boundary nor become the next edge index. */
    dvedge(edge,l_edge(edge),4,.688112897,.5);
    assert(memcmp(before_a,pt[pa]->cord,sizeof(before_a))==0);
    assert(memcmp(before_b,pt[pb]->cord,sizeof(before_b))==0);
}

static void check_straight_refinement(void)
{
    float a[3]={7.6e-4,7.34e-6,0}, b[3]={7.9e-4,7.34001e-6,0};
    int first=np, i, e=mk_edg(mk_nd(mk_pt(2,a),SiO2),
                            mk_nd(mk_pt(2,b),SiO2));
    double dx=(double)b[0]-a[0],dy=(double)b[1]-a[1];
    double x,y,error,bound;
    dvedge(e,l_edge(e),50,1.0,1.0/51.0);
    assert(np==first+52);
    for(i=first;i<np;i++) {
        x=pt[i]->cord[0];y=pt[i]->cord[1];
        error=fabs((x-a[0])*dy-(y-a[1])*dx);
        bound=FLT_EPSILON*(fabs(x*dy)+fabs(y*dx));
        assert(error<=bound);
    }
}

static int polygon(float (*coordinates)[3], int count)
{
    int region=cr_sreg(SiO2), nodes[32],i;
    assert(count<=32);
    for(i=0;i<count;i++) nodes[i]=mk_nd(mk_pt(2,coordinates[i]),SiO2);
    for(i=0;i<count;i++)
        ad_edge(region,mk_edg(nodes[i],nodes[(i+1)%count]),
                sreg[region]->bnd,MAYBE,BEFORE);
    ck_clock(region,TRUE);
    return region;
}

static double signed_side(float *a,float *b,float *c)
{
    return ((double)b[0]-a[0])*((double)c[1]-a[1])-
           ((double)b[1]-a[1])*((double)c[0]-a[0]);
}

static void check_thin_boundary_refinement(void)
{
    /* Valid retained NMOS oxide fragment. Rounded refinement positions
     * formerly crossed its opposite boundary. Keep the finite fragment. */
    float shape[9][3]={
        {0.0007891380810178816,7.343683137150947e-6,0},
        {0.000782337156124413,7.342945536947809e-6,0},
        {0.000775605789385736,7.3426035669399425e-6,0},
        {0.0007689432241022587,7.342403932852903e-6,0},
        {0.0007623487035743892,7.342317985603586e-6,0},
        {0.0007623490528203547,7.342317985603586e-6,0},
        {0.0007689438643865287,7.342401659116149e-6,0},
        {0.0007756073027849197,7.342601747950539e-6,0},
        {0.0007823397172614932,7.342943717958406e-6,0}};
    int sr=polygon(shape,9),edges[9],i;
    struct LLedge *b=sreg[sr]->bnd,*q;
    for(i=0;i<9;i++){edges[i]=b->edge;b=b->next;}
    for(i=0;i<9;i++)dvedge(edges[i],l_edge(edges[i]),50,1,1.0/51);
    b=sreg[sr]->bnd;
    do {
        q=b->next;
        while(q!=b) {
            float *a=cord_arr(pt_nd(nB(b))),*d=cord_arr(pt_nd(nF(b)));
            float *c=cord_arr(pt_nd(nB(q))),*e=cord_arr(pt_nd(nF(q)));
            assert(!(signed_side(a,d,c)*signed_side(a,d,e)<0 &&
                     signed_side(c,e,a)*signed_side(c,e,d)<0));
            q=q->next;
        }
        b=b->next;
    }while(b!=sreg[sr]->bnd);
    ck_clock(sr,TRUE);
    free_skel(sr);
}

static void check_subtraction(void)
{
    float material[3][3]={{0,1e-4,0},{-1e-4,0,0},{1e-4,0,0}};
    float mask[4][3]={{-2e-4,-1e-4,0},{2e-4,-1e-4,0},
                      {2e-4,8e-5,0},{-2e-4,8e-5,0}};
    int region=polygon(material,3), etch=polygon(mask,4), results[MAXREG],count;
    int corner=pt_nd(nB(sreg[region]->bnd));
    struct LLedge *b;
    double height=(double)material[0][1]-mask[2][1],expected=height*height;
    count=sub_skel(region,etch,results);
    assert(count==1 && sreg[results[0]]->len==3);
    assert(pt[corner]->cord[0]==material[0][0] && pt[corner]->cord[1]==material[0][1]);
    b=sreg[results[0]]->bnd;
    assert(fabs(area(pt[pt_nd(nB(b))]->cord,pt[pt_nd(nB(b->next))]->cord,
                     pt[pt_nd(nB(b->prev))]->cord)-expected)<expected*1e-5);
    while(nsreg) free_skel(nsreg-1);
}

/* Four crossings on one clipping edge: closing each retained arc on its
 * own fills the U's void and overlaps the actual material component. */
static void check_concave_subtraction(int rotation, int cut)
{
    float shape[12][3]={{0,0,0},{3,0,0},{3,2,0},{3,3,0},
                       {2,3,0},{2,2,0},{2,1,0},{1,1,0},
                       {1,2,0},{1,3,0},{0,3,0},{0,2,0}};
    float mask[4][3]={{-1,-1,0},{4,-1,0},{4,2,0},{-1,2,0}};
    double expected, sum=0, x,y;
    int i,j,r,clip,count,results[MAXREG];
    struct LLedge *b,*start;
    if (cut==0) { /* remove upper arms, retain the connected base */
        mask[0][1]=mask[1][1]=2;
        mask[2][1]=mask[3][1]=4;
        expected=5;
    } else if (cut==1) { /* retain two separate arm tips */
        expected=2;
    } else { /* middle slab: retain three separate components */
        mask[0][1]=mask[1][1]=1.5;
        mask[2][1]=mask[3][1]=2.5;
        expected=5;
    }
    for(j=0;j<2;j++) for(i=0;i<(j?4:12);i++) {
        float *p=j?mask[i]:shape[i];
        x=p[0]*1e-4; y=p[1]*1e-4;
        if(rotation==0) {p[0]=x;p[1]=y;}
        if(rotation==1) {p[0]=-y;p[1]=x;}
        if(rotation==2) {p[0]=-x;p[1]=-y;}
        if(rotation==3) {p[0]=y;p[1]=-x;}
    }
    r=polygon(shape,12); clip=polygon(mask,4);
    count=sub_skel(r,clip,results);
    if(count!=cut+1) fprintf(stderr,"rotation=%d cut=%d components=%d expected=%d\n",
                            rotation,cut,count,cut+1);
    assert(count==cut+1);
    for(i=0;i<count;i++) {
        double signed_area=0;
        start=b=sreg[results[i]]->bnd;
        do {
            float *p=pt[pt_nd(nB(b))]->cord;
            float *q=pt[pt_nd(nF(b))]->cord;
            signed_area+=(double)p[0]*q[1]-(double)p[1]*q[0];
            b=b->next;
        } while(b!=start);
        assert(signed_area>0);
        sum+=signed_area/2;
    }
    assert(fabs(sum-expected*1e-8)<1e-13);
    while(nsreg) free_skel(nsreg-1);
}

static void check_near_corner_intersection(void)
{
    /* A dry-etch boundary captured from PMOS. The first short segment
     * ends just left of the subject edge. A tolerant hit there must not
     * hide the true intersection with the subsequent horizontal segment. */
    float mask[6][3]={{-1e-4,-3.65822707e-5,0},
                     {-3.63797881e-12,-3.65822707e-5,0},
                     {-3.63797881e-12,-3.65822561e-5,0},
                     {11e-4,-3.65822561e-5,0},
                     {11e-4,-2e-4,0},{-1e-4,-2e-4,0}};
    float a[3]={0,-5e-5,0},b[3]={0,3.0401854e-5,0},c[3];
    int clip=polygon(mask,6),edge;
    struct LLedge *hit,*edg_crs();
    edge=mk_edg(mk_nd(mk_pt(2,a),SiO2),mk_nd(mk_pt(2,b),SiO2));
    hit=edg_crs(edge,clip,c);
    assert(hit && c[0]==0 && c[1]==mask[2][1]);
    assert(cordinate(pt_nd(nB(hit)),1)==mask[2][1]);
    assert(cordinate(pt_nd(nF(hit)),1)==mask[2][1]);
    while(nsreg) free_skel(nsreg-1);
}

static void check_coincident_clip_corner(void)
{
    float shape[4][3]={{0,0,0},{1e-4,0,0},{1e-4,1e-4,0},{0,1e-4,0}};
    float mask[4][3]={{0,-1e-4,0},{2e-4,-1e-4,0},
                     {2e-4,5e-5,0},{0,5e-5,0}};
    int r=polygon(shape,4), clip=polygon(mask,4),results[MAXREG],count;
    struct LLedge *b,*start;
    double sum=0;
    count=sub_skel(r,clip,results);
    assert(count==1 && sreg[results[0]]->len==4);
    start=b=sreg[results[0]]->bnd;
    do {
        float *p=pt[pt_nd(nB(b))]->cord,*q=pt[pt_nd(nF(b))]->cord;
        assert(p[0]!=q[0] || p[1]!=q[1]);
        sum+=(double)p[0]*q[1]-(double)p[1]*q[0];
        b=b->next;
    } while(b!=start);
    assert(fabs(sum/2-5e-9)<1e-15);
    while(nsreg) free_skel(nsreg-1);
}

static void check_tangent_components(void)
{
    float shape[7][3]={{0,-1e-4,0},{4e-4,-1e-4,0},{4e-4,2e-4,0},
                      {3e-4,2e-4,0},{2e-4,0,0},{1e-4,2e-4,0},{0,2e-4,0}};
    float mask[4][3]={{-1e-4,-2e-4,0},{5e-4,-2e-4,0},
                     {5e-4,0,0},{-1e-4,0,0}};
    int r=polygon(shape,7),clip=polygon(mask,4),results[MAXREG],count,i;
    double sum=0;
    struct LLedge *b,*start;
    count=sub_skel(r,clip,results);
    assert(count==2);
    for(i=0;i<count;i++) {
        start=b=sreg[results[i]]->bnd;
        do {
            float *p=pt[pt_nd(nB(b))]->cord,*q=pt[pt_nd(nF(b))]->cord;
            sum+=(double)p[0]*q[1]-(double)p[1]*q[0];
            b=b->next;
        } while(b!=start);
    }
    assert(fabs(sum/2-6e-8)<1e-14);
    while(nsreg) free_skel(nsreg-1);
}

static void check_roundoff_clip_corner(void)
{
    /* Oxidation can leave an axis boundary infinitesimally below x=0.
     * Do not insert a second, reversing corner at the same float position. */
    float shape[4][3]={{-1e-23,-1e-4,0},{1e-4,-1e-4,0},
                      {1e-4,1e-4,0},{-1e-23,1e-4,0}};
    float mask[4][3]={{-1e-4,-2e-4,0},{2e-4,-2e-4,0},
                     {1e-4,-5e-5,0},{-4e-23,-5e-5,0}};
    int r=polygon(shape,4),clip=polygon(mask,4),results[MAXREG],count;
    double sum=0;
    struct LLedge *b,*start;
    count=sub_skel(r,clip,results);
    assert(count==1 && sreg[results[0]]->len==4);
    start=b=sreg[results[0]]->bnd;
    do {
        float *p=pt[pt_nd(nB(b))]->cord,*q=pt[pt_nd(nF(b))]->cord;
        sum+=(double)p[0]*q[1]-(double)p[1]*q[0];
        b=b->next;
    } while(b!=start);
    assert(fabs(sum/2-1.5e-8)<1e-14);
    while(nsreg) free_skel(nsreg-1);
}

static void check_consistent_boundary_predicate(void)
{
    float mask[4][3]={{0,0,0},{1e-4,0,0},{1e-4,5e-5,0},{0,5e-5,0}};
    float a[3]={3e-5,5.0002e-5,0},b[3]={7e-5,5.0002e-5,0},c[3];
    int clip=polygon(mask,4),pa=mk_pt(2,a),pb=mk_pt(2,b),edge;
    struct LLedge *hit,*edg_crs();
    /* Nearby parallel edges are outside, not falsely classified ON. */
    assert(!on_bound(pa,clip) && !on_bound(pb,clip));
    edge=mk_edg(mk_nd(pa,SiO2),mk_nd(pb,SiO2));
    assert(edg_crs(edge,clip,c)==NULL);
    /* Collinear ON endpoints have a defined intersection too. */
    pt[pa]->cord[1]=pt[pb]->cord[1]=5e-5;
    edg[edge]->len=l_edge(edge);
    assert(on_bound(pa,clip) && on_bound(pb,clip));
    hit=edg_crs(edge,clip,c);
    assert(hit && c[0]==pt[pa]->cord[0] && c[1]==pt[pa]->cord[1]);
    while(nsreg) free_skel(nsreg-1);
}

static void check_multiple_crossings(void)
{
    /* The notch crosses one unsplit subject edge twice; all four original
     * subject vertices are outside. Endpoint-only clipping misses it. */
    float shape[4][3]={{0,0,0},{4e-4,0,0},{4e-4,2e-4,0},{0,2e-4,0}};
    float mask[8][3]={{-1e-4,-1e-4,0},{5e-4,-1e-4,0},
                     {5e-4,-5e-5,0},{3e-4,-5e-5,0},
                     {3e-4,1e-4,0},{1e-4,1e-4,0},
                     {1e-4,-5e-5,0},{-1e-4,-5e-5,0}};
    int r=polygon(shape,4),clip=polygon(mask,8),results[MAXREG],count;
    int corner1=0,corner2=0;
    double sum=0;
    struct LLedge *b,*start;
    count=sub_skel(r,clip,results);
    assert(count==1);
    start=b=sreg[results[0]]->bnd;
    do {
        float *p=pt[pt_nd(nB(b))]->cord,*q=pt[pt_nd(nF(b))]->cord;
        sum+=(double)p[0]*q[1]-(double)p[1]*q[0];
        if(p[0]==mask[4][0] && p[1]==mask[4][1])corner1=1;
        if(p[0]==mask[5][0] && p[1]==mask[5][1])corner2=1;
        b=b->next;
    } while(b!=start);
    assert(corner1 && corner2 && fabs(sum/2-6e-8)<1e-14);
    while(nsreg) free_skel(nsreg-1);
}

static void check_remove_material_junction(void)
{
    float a[3]={0,0,0},b[3]={1e-7,0,0};
    int pa=mk_pt(2,a),pb=mk_pt(2,b),na,nb,oxide,edge;
    na=mk_nd(pa,Si);oxide=mk_nd(pa,SiO2);
    nb=mk_nd(pb,Si);mk_nd(pb,GAS);
    edge=mk_edg(na,nb);
    /* The material search ends with node_mat(...,SiO2)==-1. That
     * lookup result must never be passed to dead_nd(). */
    rem_edg(edge);
    assert(pt_nd(oxide)==pb && node_mat(nb,SiO2)==oxide);
    assert(dead_nd(na) && !dead_nd(nb));
}

/* A rounded cut shared with a neighboring oxide region is ON under the
 * classifier. Re-intersecting the nearly parallel edge shifts that junction
 * by 28 float ULPs. The retained arc must use the existing point identity. */
static void check_rounded_shared_cut(int reverse, int start, int variant)
{
    float shape[7][3]={
        {0.00045063867582939565,7.3422506829956546e-6,0},
        {0.00045703299110755324,7.3422161221969873e-6,0},
        {0.00045703334035351872,7.3422161221969873e-6,0},
        {0.00046336193918250501,7.3422224886598997e-6,0},
        {0.00046336193918250501,9.0098219516221434e-6,0},
        {0.00045063867582939565,9.0098501459578983e-6,0},
        {0.000449,7.0e-6,0}};
    float mask[7][3]={
        {0.00044,-1e-4,0},{0.00047,-1e-4,0},
        {0.00047,7.342225e-6,0},
        {0.00046336199739016593,7.3422197601757944e-6,0},
        {0.00045703272917307913,7.3422161221969873e-6,0},
        {0.00045063853031024337,7.3422488640062511e-6,0},
        {0.000448,7.34226e-6,0}};
    float ordered[7][3];
    int i,j,r,clip,results[MAXREG],count,k,point=-1,found=0;
    double x,y,scale=variant==1 ? .5 : variant==2 ? 2 : 1;
    struct LLedge *bp,*first;
    for(j=0;j<2;j++)for(i=0;i<7;i++) {
        float *p=j ? mask[i]:shape[i];
        x=(double)p[0]*scale; y=(double)p[1]*scale;
        if(variant==3)x-=(double)(float)0.0004;
        if(variant==4) {p[0]=-y;p[1]=x;}
        else {p[0]=x;p[1]=y;}
    }
    for(i=0;i<7;i++)memcpy(ordered[i],shape[(start+(reverse ? 7-i:i))%7],sizeof(ordered[i]));
    r=polygon(ordered,7);clip=polygon(mask,7);
    first=bp=sreg[r]->bnd;
    do {
        i=pt_nd(nB(bp));
        if(pt[i]->cord[0]==shape[2][0] && pt[i]->cord[1]==shape[2][1])point=i;
        bp=bp->next;
    }while(bp!=first);
    assert(point>=0 && on_bound(point,clip));
    count=sub_skel(r,clip,results);
    if(count<1)fprintf(stderr,"rounded cut reverse=%d start=%d count=%d\n",reverse,start,count);
    assert(count>=1);
    for(k=0;k<count;k++) {
        first=bp=sreg[results[k]]->bnd;
        do { if(pt_nd(nB(bp))==point)found=1;bp=bp->next; }while(bp!=first);
    }
    if(!found) {
        fprintf(stderr,"rounded shared cut reverse=%d start=%d lost point %.17g,%.17g\n",reverse,start,
                (double)pt[point]->cord[0],(double)pt[point]->cord[1]);
        first=bp=sreg[results[0]]->bnd;
        do { i=pt_nd(nB(bp));fprintf(stderr,"retained %.17g,%.17g\n",(double)pt[i]->cord[0],(double)pt[i]->cord[1]);bp=bp->next; }while(bp!=first);
    }
    assert(found); /* No coordinate-based substitution for the shared junction. */
    while(nsreg)free_skel(nsreg-1);
}

/* Captured spacer-mask endpoint: line roundoff must not classify a point
 * beyond the finite segment as ON. Exercise rotation and scale as well. */
static void check_finite_clip_endpoints(void)
{
    int rot, scale, i, r, p;
    for (rot=0;rot<4;rot++) for(scale=0;scale<2;scale++) {
        float mask[3][3]={{5.684342356405358e-21,-1.0806987120304257e-5,0},
                         {-.001,-.00095409882487729192,0},
                         {.002,-.00095409882487729192,0}};
        float probe[3]={5.4492412069028113e-27,-1.0806985301314853e-5,0};
        for(i=0;i<4;i++) {
            float *v=i==3 ? probe:mask[i];
            float x=v[0]*(scale ? 2.0f:1.0f),y=v[1]*(scale ? 2.0f:1.0f);
            v[0]=rot==0 ? x:rot==1 ? -y:rot==2 ? -x:y;
            v[1]=rot==0 ? y:rot==1 ? x:rot==2 ? -y:-x;
        }
        r=polygon(mask,3);p=mk_pt(2,probe);
        assert(!on_bound(p,r));
        p=mk_pt(2,mask[0]);assert(on_bound(p,r));
        while(nsreg)free_skel(nsreg-1);
    }
}

/* The NMOS spacer cut meets the right side at a clipping vertex. A nearby
 * original vertex is ON under storage-roundoff classification, but is below
 * the actual cut. Replacing the intersection by that vertex creates a
 * collinear reversal when the clip walker subsequently visits the corner. */
static void check_constructed_corner(int variant)
{
    float shape[5][3]={{.00099,-.00004,0},{.001,-.00004,0},
        {.001,-2.5622492103138939e-5,0},{.001,-.00002,0},{.00099,-.00002,0}};
    float mask[5][3]={{0,-.001,0},{.002,-.001,0},
        {.002,-.00098533881828188896,0},{.001,-2.5622462999308482e-5,0},
        {0,-2.5622462999308482e-5,0}};
    int i,j,subject,clip,results[MAXREG],count;
    double scale=variant==1?.5:variant==2?2:1,actual=0,expected;
    struct LLedge *bp;
    for(j=0;j<2;j++)for(i=0;i<5;i++) {
        float *p=j?mask[i]:shape[i];
        double x=p[0]*scale,y=p[1]*scale;
        if(variant==3){p[0]=-y;p[1]=x;}else{p[0]=x;p[1]=y;}
    }
    expected=area(shape[4],shape[3],mask[3])*2;
    subject=polygon(shape,5);clip=polygon(mask,5);
    count=sub_skel(subject,clip,results);
    assert(count==1);
    bp=sreg[results[0]]->bnd;
    do {
        float *a=cord_arr(pt_nd(nB(bp))),*b=cord_arr(pt_nd(nF(bp)));
        float *c=cord_arr(pt_nd(nF(bp->next)));
        double ux=(double)b[0]-a[0],uy=(double)b[1]-a[1];
        double vx=(double)c[0]-b[0],vy=(double)c[1]-b[1];
        assert(signed_side(a,b,c)!=0 || ux*vx+uy*vy>=0);
        actual+=(double)a[0]*b[1]-(double)a[1]*b[0];
        bp=bp->next;
    }while(bp!=sreg[results[0]]->bnd);
    assert(fabs(actual/2-expected)<expected*1e-10);
    while(nsreg)free_skel(nsreg-1);
}

/* A source edge starts at a clip corner, overlaps a straight clip piece,
 * then exits at an interior corner. Both contacts must be noded; choosing
 * the first edge hit at the starting endpoint makes the closure reverse. */
static void check_overlap_endpoint(int rotation)
{
    float shape[6][3]={{0,-1,0},{4,-1,0},{4,0,0},
                       {3.25,0,0},{2.75,0,0},{0,0,0}};
    float mask[9][3]={{0,-2,0},{4,-2,0},{4,-.25,0},
        {3.25,0,0},{3,0,0},{3,-.5,0},{1,-.5,0},{1,0,0},{0,0,0}};
    int i,j,subject,clip,count,results[MAXREG];
    double actual=0,expected;
    struct LLedge *bp;
    for(j=0;j<2;j++)for(i=0;i<(j?9:6);i++) {
        float *p=j?mask[i]:shape[i];
        float x=p[0]*1e-4,y=p[1]*1e-4;
        p[0]=rotation ? -y:x;p[1]=rotation ? x:y;
    }
    expected=area(shape[2],mask[3],mask[2])+2*area(mask[4],mask[5],mask[6]);
    subject=polygon(shape,6);clip=polygon(mask,9);
    count=sub_skel(subject,clip,results);
    assert(count==2);
    for(i=0;i<count;i++) {
        bp=sreg[results[i]]->bnd;
        do {
            float *a=cord_arr(pt_nd(nB(bp))),*b=cord_arr(pt_nd(nF(bp)));
            float *c=cord_arr(pt_nd(nF(bp->next)));
            double ux=(double)b[0]-a[0],uy=(double)b[1]-a[1];
            double vx=(double)c[0]-b[0],vy=(double)c[1]-b[1];
            assert(signed_side(a,b,c)!=0 || ux*vx+uy*vy>=0);
            actual+=(double)a[0]*b[1]-(double)a[1]*b[0];
            bp=bp->next;
        }while(bp!=sreg[results[i]]->bnd);
    }
    assert(fabs(actual/2-expected)<expected*1e-10);
    while(nsreg)free_skel(nsreg-1);
}

/* Both endpoints of the retained upper edge lie on different clip walls.
 * Each lower edge has only one proper crossing, so the old multiple-crossing
 * prepass returned early. The arc walker then replaced the upper subject edge
 * with the recessed clip boundary, inventing material inside the opaque gate. */
static void check_on_on_retained_arc(int rotation)
{
    float shape[7][3]={{0,0,0},{4,0,0},{6,0,0},{10,0,0},
                       {10,.02,0},{5,.02,0},{0,.02,0}};
    float mask[8][3]={{-1,-2,0},{11,-2,0},{11,1,0},{6,1,0},
                      {6,-.25,0},{4,-.25,0},{4,1,0},{-1,1,0}};
    int i,j,k,subject,clip,count,results[MAXREG];
    double actual=0,expected;
    struct LLedge *bp;
    expected=2e-4*(double)(float)(.02e-4);
    for(j=0;j<2;j++)for(i=0;i<(j?8:7);i++) {
        float *p=j?mask[i]:shape[i];
        float x=p[0]*1e-4,y=p[1]*1e-4,t;
        for(k=0;k<rotation;k++){t=x;x=-y;y=t;}
        p[0]=x;p[1]=y;
    }
    subject=polygon(shape,7);clip=polygon(mask,8);
    count=sub_skel(subject,clip,results);assert(count==1);
    bp=sreg[results[0]]->bnd;
    do {
        float *a=cord_arr(pt_nd(nB(bp))),*b=cord_arr(pt_nd(nF(bp)));
        actual+=(double)a[0]*b[1]-(double)a[1]*b[0];bp=bp->next;
    }while(bp!=sreg[results[0]]->bnd);
    assert(fabs(actual/2-expected)<expected*1e-6);
    while(nsreg)free_skel(nsreg-1);
}

int main(void)
{
    int i,j,k;
    mode=TWOD; nvrt=nedg=3; n_imp=1; soltoimp[0]=B;
    check(0); check(1); check_minimum_spacing();
    check_shared_crossing_endpoint(); check_straight_refinement();
    check_thin_boundary_refinement(); check_subtraction();
    for(i=0;i<4;i++) for(j=0;j<3;j++) check_concave_subtraction(i,j);
    check_near_corner_intersection();
    check_coincident_clip_corner();
    check_tangent_components();
    check_roundoff_clip_corner();
    check_consistent_boundary_predicate();
    check_multiple_crossings();
    check_refinement_endpoint_status();
    check_remove_material_junction();
    check_finite_clip_endpoints();
    for(i=0;i<4;i++)check_constructed_corner(i);
    for(i=0;i<2;i++)check_overlap_endpoint(i);
    for(i=0;i<4;i++)check_on_on_retained_arc(i);
    for(k=0;k<5;k++)for(i=0;i<2;i++)for(j=0;j<7;j++)check_rounded_shared_cut(i,j,k);
    for(i=0;i<ned;i++) { dis_list(&edg[i]->ele); dis_list(&edg[i]->skel); free(edg[i]); }
    for(i=0;i<nn;i++) { dis_list(&nd[i]->edg); dis_list(&nd[i]->ltri); free(nd[i]); }
    for(i=0;i<np;i++) free(pt[i]);
    free(edg); free(nd); free(pt);
    puts("PASS: etch preserves corners, interfaces, and concave material components");
    return 0;
}
