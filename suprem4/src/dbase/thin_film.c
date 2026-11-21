/*************************************************************************
 * thin_film.c - through-thickness resolution of dielectric films.
 *
 * Dopant transport through a thin dielectric (a gate oxide between n+ poly
 * and the channel) is a barrier problem: the diffusion length in the film
 * is a fraction of its thickness. A film meshed one element thick couples
 * its two interfaces directly, and the discrete transport through it is
 * grossly wrong. For the RIT gate stack (16.9 nm oxide, P-doped poly, 117 min
 * activation) the phosphorus dose entering the channel is
 *     1, 2, 4, 8, 16 element layers:  7.7e10, 2.5e10, 4.3e8, 1.4e8, 1.0e8 cm-2
 * so a one-layer oxide inverts the NMOS channel. Physics and parameters are
 * identical in each case; only the discretization differs.
 *
 * thin_film_conform(L) enforces at least L element layers across every
 * dielectric film that separates two distinct semiconductor regions
 * (see tf_face_class). For each neighbour class c
 * a film coordinate s_c = d_c/(d_c + d_other) is evaluated at every dielectric
 * node, where d_c is the Euclidean distance to the nearest interface segment
 * of class c and d_other the distance to the nearest interface segment of any
 * other class (reflecting sides are not interfaces). s_c depends only on the
 * film geometry, never on the mesh inside it, so the levels are fixed and the
 * construction is idempotent: a film that already has the layers is left alone.
 *
 * An edge of the dielectric (interior or boundary) whose jump |ds_c| exceeds
 * 1/L for some c is split at its crossings of the levels j/M, M = ceil(3L/2),
 * at the exact parametric positions on that edge. Inserted points lie within
 * about one level of their level set, so edges between neighbouring levels
 * jump by about 1/M < 1/L and are not split again; any path between two
 * distinct interfaces therefore crosses at least L edges. Points are added by
 * conforming bisection (sp_gredg), which gives each new node the linear
 * interpolant of every field: the piecewise-linear representation of every
 * dopant is unchanged exactly. Refinement is anisotropic: edges along a film
 * have ds ~ 0 and are never split. A split is made only if every child
 * triangle stays positively oriented (exact predicate).
 *************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "thin_film.h"
extern int sp_gredg();

/* Dielectrics in which dopants diffuse and whose thickness is a barrier. */
static int thin_film_dielectric(int material)
{
    return material==SiO2 || material==OxNi;
}

#define TF_CLASSES MAXREG
struct tf_segment {double x0,y0,x1,y1; int cls;};

static double tf_segment_distance(double px,double py,const struct tf_segment *g)
{
    double dx=g->x1-g->x0,dy=g->y1-g->y0,l2=dx*dx+dy*dy,t=0;
    if(l2>0) {
        t=((px-g->x0)*dx+(py-g->y0)*dy)/l2;
        if(t<0)t=0;
        if(t>1)t=1;
    }
    return hypot(px-(g->x0+t*dx),py-(g->y0+t*dy));
}

/* Interface class of face f of dielectric triangle t: the region across it
 * when that region is a semiconductor (a dopant reservoir), otherwise -1.
 * The barrier problem is transport between two reservoirs through the film:
 * gate oxide between poly and channel, interpoly dielectric, buried oxide.
 * Gas carries no impurities and its surface is a no-flux boundary, metal,
 * nitride and resist are not dopant reservoirs, and a dielectric neighbour,
 * reflecting side or backside is not an interface; a film facing only one
 * reservoir cannot leak and is not refined. Classes are regions, not
 * materials, so two distinct regions of one semiconductor (SOI device layer
 * and substrate) are distinct reservoirs. Edges join material-specific nodes,
 * so the medium across a face is found from the face neighbour. */
static int tf_face_class(int t,int f)
{
    int nb=neigh_fc(t,f),r;
    if(nb<0 || dead_tri(nb))return -1;
    r=reg_tri(nb);
    if(r<0 || r>=TF_CLASSES || INSULATE(mat_reg(r)))return -1;
    return r;
}

/* An edge of at least one live dielectric triangle, with dielectric nodes.
 * Film boundary edges count: an edge along a reflecting side, an end face or
 * a gas surface couples the two film interfaces just as an interior chord
 * does. sp_gredg splits a boundary edge conformingly in every material. */
static int thin_film_interior(int e)
{
    int k,t,dielectric=0;
    if(dead_edg(e) || num_tri_edg(e)<1)return 0;
    for(k=0;k<num_tri_edg(e);k++) {
        t=tri_edg(e,k);
        if(dead_tri(t) || !leaf(tri[t]))return 0;
        if(thin_film_dielectric(mat_reg(reg_tri(t))))dielectric=1;
    }
    return dielectric && thin_film_dielectric(mat_nd(nd_edg(e,0))) &&
           thin_film_dielectric(mat_nd(nd_edg(e,1)));
}

/* Would splitting edge e at c leave every child triangle positive? */
static int thin_film_children_valid(int e,const float *c)
{
    int k,j,t,a,b,child;
    float v[3][MAXDIM];
    for(k=0;k<num_tri_edg(e);k++) {
        t=tri_edg(e,k);a=b=-1;
        for(j=0;j<3;j++) {
            if(vert_tri(t,j)==nd_edg(e,0))a=j;
            if(vert_tri(t,j)==nd_edg(e,1))b=j;
            memcpy(v[j],cord_arr(pt_nd(vert_tri(t,j))),sizeof(v[j]));
        }
        if(a<0 || b<0)return 0;
        for(child=0;child<2;child++) {
            float saved[MAXDIM];
            j=child?b:a;memcpy(saved,v[j],sizeof(saved));memcpy(v[j],c,sizeof(saved));
            if(orient2d_float_sign(v[0],v[1],v[2])<=0)return 0;
            memcpy(v[j],saved,sizeof(saved));
        }
    }
    return 1;
}

static int thin_film_round(int layers)
{
    int *index=NULL,count=0,i,j,k,e,c,classes[TF_CLASSES],nclass=0,pass=0,splits=0,refused=0;
    int nseg=0,segcap=0,levels_m=(3*layers+1)/2;
    struct tf_segment *seg=NULL;
    double *s=NULL; /* s[class*count + dielectric index], -1 = undefined */
    if(mode!=TWOD || layers<2)return 0;
    /* Interface segments of every class bounding the dielectric. */
    memset(classes,0,sizeof(classes));
    for(i=0;i<ne;i++) {
        int f;
        if(dead_tri(i) || !leaf(tri[i]) || !thin_film_dielectric(mat_reg(reg_tri(i))))continue;
        for(f=0;f<num_face(i);f++) {
            int cls=tf_face_class(i,f),pa,pb;
            if(cls<0)continue;
            if(nseg==segcap) {
                segcap=segcap?2*segcap:256;
                seg=realloc(seg,(size_t)segcap*sizeof(*seg));
                if(!seg)panic("thin film conformity: out of memory");
            }
            pa=pt_nd(nd_face_ele(i,0,f));pb=pt_nd(nd_face_ele(i,1,f));
            seg[nseg].x0=xcord(pa);seg[nseg].y0=ycord(pa);
            seg[nseg].x1=xcord(pb);seg[nseg].y1=ycord(pb);
            seg[nseg].cls=cls;classes[cls]=1;nseg++;
        }
    }
    for(c=0;c<TF_CLASSES;c++)nclass+=classes[c];
    if(nclass<2){free(seg);return 0;}
    /* Film coordinates at every dielectric node. */
    index=malloc((size_t)nn*sizeof(int));
    if(!index)panic("thin film conformity: out of memory");
    for(i=0;i<nn;i++) {
        index[i]=-1;
        if(dead_nd(i) || !thin_film_dielectric(mat_nd(i)) || num_tri_nd(i)==0)continue;
        index[i]=count++;
    }
    s=malloc((size_t)TF_CLASSES*(count?count:1)*sizeof(double));
    if(!s)panic("thin film conformity: out of memory");
    for(i=0;i<nn;i++) if(index[i]>=0) {
        double d[TF_CLASSES],px=xcord(pt_nd(i)),py=ycord(pt_nd(i));
        for(c=0;c<TF_CLASSES;c++)d[c]=INFINITY;
        for(k=0;k<nseg;k++) {
            double dist=tf_segment_distance(px,py,&seg[k]);
            if(dist<d[seg[k].cls])d[seg[k].cls]=dist;
        }
        /* Nearest and second-nearest class: the nearest other class of c. */
        double d1=INFINITY,d2=INFINITY;
        int c1=-1;
        for(c=0;c<TF_CLASSES;c++) if(classes[c]) {
            if(d[c]<d1){d2=d1;d1=d[c];c1=c;}
            else if(d[c]<d2)d2=d[c];
        }
        for(c=0;c<TF_CLASSES;c++) {
            double other=c==c1?d2:d1,v=-1;
            if(classes[c]) {
                if(d[c]==0)v=0;
                else if(other==0)v=1;
                else if(isfinite(d[c]) && isfinite(other))v=d[c]/(d[c]+other);
            }
            s[(size_t)c*count+index[i]]=v;
        }
    }
    free(seg);
    /* Insert the iso-levels j/K of every coordinate on each original interior
     * dielectric edge that crosses them, at their exact parametric positions.
     * Each original edge receives at most K-1 points per coordinate, so the
     * construction is bounded; the Delaunay pass in bd_connect then connects
     * the inserted points into layers. */
    {
        int limit=ned,nlev,lev;
        double *levels=malloc((size_t)TF_CLASSES*levels_m*sizeof(double));
        if(!levels)panic("thin film conformity: out of memory");
        for(e=0;e<limit;e++) {
            int a,b,q;
            double ax,ay,bx,by;
            if(!thin_film_interior(e))continue;
            a=nd_edg(e,0);b=nd_edg(e,1);
            if(index[a]<0 || index[b]<0)continue;
            /* Parametric positions (from a) of all level crossings. */
            nlev=0;
            for(c=0;c<TF_CLASSES;c++) {
                double sa=s[(size_t)c*count+index[a]],sb=s[(size_t)c*count+index[b]];
                if(!classes[c] || sa<0 || sb<0 || fabs(sb-sa)*layers<=1+1e-9)continue;
                for(lev=1;lev<levels_m;lev++) {
                    double l=(double)lev/levels_m,t=(l-sa)/(sb-sa);
                    if(t>1e-6 && t<1-1e-6)levels[nlev++]=t;
                }
            }
            if(!nlev)continue;
            /* Split from the far end so the remaining edge keeps endpoint a. */
            for(i=1;i<nlev;i++)for(j=i;j>0 && levels[j]>levels[j-1];j--) {
                double tmp=levels[j];levels[j]=levels[j-1];levels[j-1]=tmp;
            }
            ax=xcord(pt_nd(a));ay=ycord(pt_nd(a));bx=xcord(pt_nd(b));by=ycord(pt_nd(b));
            for(q=0;q<nlev;q++) {
                float point_c[MAXDIM]={0};
                int point,node,end=nd_edg(e,1);
                if(q && levels[q]>levels[q-1]-1e-6)continue; /* coincident levels */
                point_c[0]=ax+levels[q]*(bx-ax);point_c[1]=ay+levels[q]*(by-ay);
                if((point_c[0]==xcord(pt_nd(a)) && point_c[1]==ycord(pt_nd(a))) ||
                   (point_c[0]==xcord(pt_nd(end)) && point_c[1]==ycord(pt_nd(end))) ||
                   !thin_film_children_valid(e,point_c)) {refused++;continue;}
                sp_gredg(e,point_c,&point);
                node=node_mat(nd_pt(point,0),mat_nd(a));
                if(node<0)panic("thin film conformity: split lost its dielectric node");
                splits++;
            }
        }
        free(levels);
        pass=0;
    }
    free(s);free(index);
    if(splits) {
        bd_connect("thin film conformity");
        if(getenv("SUPREM4_MESH_DIAG") || getenv("SUPREM4_THIN_FILM_DIAG"))
            fprintf(stderr,"[thin-film] layers=%d splits=%d refused=%d passes=%d\n",layers,splits,refused,pass+1);
    }
    return splits;
}

/* Delaunay flips after insertion can reconnect two levels directly. Recompute
 * the coordinates on the refined mesh and repeat; each round inserts at most
 * K-1 points per offending edge, and rounds are capped. */
int thin_film_conform(int layers)
{
    int round,splits,total=0;
    for(round=0;round<3;round++) {
        splits=thin_film_round(layers);
        total+=splits;
        if(!splits)break;
    }
    return total;
}
