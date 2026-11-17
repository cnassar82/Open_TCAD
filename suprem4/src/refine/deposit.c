/*************************************************************************
 *									 *
 *     Copyright c 1984 The board of trustees of the Leland Stanford 	 *
 *                      Junior University. All rights reserved.		 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of Stanford University. *
 *									 *
 *     Copyright c 1991 The board of regents of the University of 	 *
 *                      Florida.  All rights reserved.			 *
 *     This subroutine may not be used outside of the SUPREM4 computer	 *
 *     program without the prior written consent of the University of 	 *
 *     Florida.								 *
 *									 *
 *									 *
 *************************************************************************/
/*   deposit.c                Version 5.1     */
/*   Last Modification : 7/3/91 15:40:39 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <float.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "impurity.h"
#include "regrid.h"
#include "refine.h"
#include "poly.h"
#include "process_surface.h"
#include "matrix.h"

/* Boundary subdivisions alone do not constrain interior elements.  A size
 * queue drives conforming longest-edge bisection of the new film, with an
 * h*sqrt(2) diagonal bound.  Lazy keys discard entries invalidated by a split;
 * no O(N^2) repeated whole-mesh scan is needed for fine films. */
struct deposition_edge { int edge; double length2; };
struct deposition_heap {
    struct deposition_edge *entry;
    size_t count,capacity;
};

static int deposition_region_edge(int e,int region)
{
    int k,t;
    if(dead_edg(e))return FALSE;
    for(k=0;k<num_tri_edg(e);k++) {
        t=tri_edg(e,k);
        if(!dead_tri(t) && leaf(tri[t]) && reg_tri(t)==region)return TRUE;
    }
    return FALSE;
}

static void deposition_push(struct deposition_heap *q,int e,int region,double limit)
{
    size_t i,parent;
    struct deposition_edge value,*grown;
    if(!deposition_region_edge(e,region))return;
    value.edge=e;value.length2=len_edg(e)*len_edg(e);
    if(value.length2<=limit)return;
    if(q->count==q->capacity) {
        size_t capacity=q->capacity ? 2*q->capacity:256;
        grown=realloc(q->entry,capacity*sizeof(*grown));
        if(!grown)panic("cannot allocate deposition size queue");
        q->entry=grown;q->capacity=capacity;
    }
    i=q->count++;
    while(i) {
        parent=(i-1)/2;
        if(q->entry[parent].length2>=value.length2)break;
        q->entry[i]=q->entry[parent];i=parent;
    }
    q->entry[i]=value;
}

static struct deposition_edge deposition_pop(struct deposition_heap *q)
{
    struct deposition_edge result=q->entry[0],last=q->entry[--q->count];
    size_t i=0,child;
    while(2*i+1<q->count) {
        child=2*i+1;
        if(child+1<q->count && q->entry[child+1].length2>q->entry[child].length2)child++;
        if(last.length2>=q->entry[child].length2)break;
        q->entry[i]=q->entry[child];i=child;
    }
    if(q->count)q->entry[i]=last;
    return result;
}

static int deposition_incircle(int e)
{
    int t=tri_edg(e,0),u=tri_edg(e,1);
    float *a=cord_arr(pt_nd(vert_tri(t,ewhich(t,e))));
    float *b=cord_arr(pt_edg(e,0)),*c=cord_arr(pt_edg(e,1));
    float *d=cord_arr(pt_nd(vert_tri(u,ewhich(u,e))));
    long double ax=(long double)a[0]-d[0],ay=(long double)a[1]-d[1];
    long double bx=(long double)b[0]-d[0],by=(long double)b[1]-d[1];
    long double cx=(long double)c[0]-d[0],cy=(long double)c[1]-d[1];
    long double al=ax*ax+ay*ay,bl=bx*bx+by*by,cl=cx*cx+cy*cy;
    long double determinant=al*(bx*cy-by*cx)+bl*(cx*ay-cy*ax)+cl*(ax*by-ay*bx);
    long double permanent=al*(fabsl(bx*cy)+fabsl(by*cx))+
                          bl*(fabsl(cx*ay)+fabsl(cy*ax))+
                          cl*(fabsl(ax*by)+fabsl(ay*bx));
    /* A filtered strict decision: leave numerically cocircular edges to
     * the global geometry builder, rather than oscillating a local queue. */
    return orient2d_float_sign(a,b,c)*determinant>64*LDBL_EPSILON*permanent;
}

/* Legalize as points are inserted.  Deferring all flips until the end
 * needlessly refines the long fans produced by a sparse offset boundary.
 * Only the new film is rezoned; every flip uses the existing conservative
 * field transfer and strict convexity check. */
static void deposition_legalize(int point,int region,
                               struct deposition_heap *sizes,double limit)
{
    int *work=NULL,e,t,u,i,j,node,count=0,capacity=0;
    int k,l,m,n,other,side,opposite;
    /* sp_edge updates edge ownership immediately, but deliberately defers
     * triangle neighbor rebuilding. Reconnect this complete point star on
     * all material sides before using the diagonal-swap API. */
    for(i=0;i<num_nd(point);i++)for(j=0;j<num_tri_nd(nd_pt(point,i));j++) {
        t=tri_nd(nd_pt(point,i),j);
        if(dead_tri(t) || !leaf(tri[t]))continue;
        for(k=0;k<3;k++)tri[t]->nb[k]=BC_OFFSET;
    }
    for(i=0;i<num_nd(point);i++)for(j=0;j<num_tri_nd(nd_pt(point,i));j++) {
        t=tri_nd(nd_pt(point,i),j);
        if(dead_tri(t) || !leaf(tri[t]))continue;
        for(k=0;k<3;k++) {
            int p=pt_nd(vert_tri(t,k));
            for(l=0;l<num_nd(p);l++) {
                n=nd_pt(p,l);
                for(m=0;m<num_tri_nd(n);m++) {
                    other=tri_nd(n,m);
                    if(other==t || dead_tri(other) || !leaf(tri[other]))continue;
                    if(tnabor(t,other,&side,&opposite)) {
                        tri[t]->nb[side]=other;tri[other]->nb[opposite]=t;
                    }
                }
            }
        }
    }
#define DEPOSITION_CHECK_EDGE(E) do { \
    if(count==capacity) { \
        int *grown;capacity=capacity?2*capacity:64; \
        grown=realloc(work,capacity*sizeof(*work)); \
        if(!grown)panic("cannot allocate deposition legalization queue"); \
        work=grown; \
    } \
    work[count++]=(E); \
} while(0)
    for(i=0;i<num_nd(point);i++) {
        node=nd_pt(point,i);
        for(j=0;j<num_tri_nd(node);j++) {
            t=tri_nd(node,j);
            if(dead_tri(t) || !leaf(tri[t]) || reg_tri(t)!=region)continue;
            DEPOSITION_CHECK_EDGE(edg_ele(t,nwhich(t,node)));
        }
    }
    while(count) {
        e=work[--count];
        if(dead_edg(e) || num_tri_edg(e)!=2)continue;
        t=tri_edg(e,0);u=tri_edg(e,1);
        if(reg_tri(t)!=region || reg_tri(u)!=region)continue;
        if(!deposition_incircle(e) || ck_flip(e)!=1)continue;
        if(!tri_flip(e,t,u))continue;
        edg[e]->len=l_edge(e);
        for(i=0;i<3;i++) {
            int f=edg_ele(t,i),g=edg_ele(u,i);
            deposition_push(sizes,f,region,limit);
            deposition_push(sizes,g,region,limit);
            if(f!=e)DEPOSITION_CHECK_EDGE(f);
            if(g!=e)DEPOSITION_CHECK_EDGE(g);
        }
    }
    free(work);
#undef DEPOSITION_CHECK_EDGE
}

static int deposition_spacing(int region,double h)
{
    int i,j,k,t,e,point,a,b,splits=0;
    double limit=2*h*h;
    float c[MAXDIM],v[3][MAXDIM];
    struct deposition_heap queue={0};
    struct deposition_edge next;
    if(mode!=TWOD)return 0;
    for(;;) {
      for(e=0;e<ned;e++)deposition_push(&queue,e,region,limit);
      if(!queue.count)break;
      while(queue.count) {
        next=deposition_pop(&queue);e=next.edge;
        if(!deposition_region_edge(e,region) ||
           next.length2!=len_edg(e)*len_edg(e))continue;
        for(j=0;j<MAXDIM;j++)c[j]=0;
        for(j=0;j<2;j++)c[j]=0.5*((double)cordinate(pt_edg(e,0),j)+
                                               cordinate(pt_edg(e,1),j));
        /* Reject an unrepresentable split instead of silently violating h.
         * Check both children on both material sides before any mutation. */
        for(i=0;i<num_tri_edg(e);i++) {
            t=tri_edg(e,i);
            if(dead_tri(t) || !leaf(tri[t]))continue;
            for(a=0;a<3;a++)if(vert_tri(t,a)==nd_edg(e,0))break;
            for(b=0;b<3;b++)if(vert_tri(t,b)==nd_edg(e,1))break;
            if(a==3 || b==3)panic("deposition spacing: inconsistent edge");
            for(k=0;k<3;k++)memcpy(v[k],cord_arr(pt_nd(vert_tri(t,k))),sizeof(v[k]));
            memcpy(v[a],c,sizeof(c));
            if(orient2d_float_sign(v[0],v[1],v[2])<=0)goto unresolved;
            memcpy(v[a],cord_arr(pt_nd(vert_tri(t,a))),sizeof(c));
            memcpy(v[b],c,sizeof(c));
            if(orient2d_float_sign(v[0],v[1],v[2])<=0)goto unresolved;
        }
        /* Preserve endpoints and prolong all fields on each material side. */
        if(sp_edge_keep(e,c,&point,FALSE,nd_edg(e,0))<0)goto unresolved;
        splits++;
        for(i=0;i<num_nd(point);i++) {
            int node=nd_pt(point,i);
            for(j=0;j<num_edge_nd(node);j++)
                deposition_push(&queue,edge_nd(node,j),region,limit);
        }
        deposition_legalize(point,region,&queue,limit);
      }
      /* Geometry rebuilding legalizes diagonals.  A legal flip can create
       * a longer edge, so size and Delaunay constraints must be satisfied
       * together, not checked on an intermediate triangulation. */
      bd_connect("deposition interior spacing");UpdateSymbolic++;
    }
    free(queue.entry);
    return 0;
unresolved:
    free(queue.entry);
    fprintf(stderr,"ERROR: requested deposition spacing cannot be represented without invalid elements\n");
    return -1;
}

/************************************************************************
 *									*
 *	user_deposit( ) - get the parameters for deposition of a new	*
 *  material.								*
 *									*
 *  Original:	CSR	3/86						*
 *		MEL	4/87   rewrote the work portion			*
 *		LCC     4/90   add 1-D mode				*
 ************************************************************************/
user_deposit( par, param )
char *par;
struct par_str *param;
{
    register int mater, ndiv, i;
    double thick;
    struct line *p;
    struct line *op;
    int nump, numop;
    int square = TRUE;
    int iskel;
    extern string_offset();
    char *sfile;
    int imp, impa, sol, r;
    double conc, space;
    double temp, pressure;
    int oldnn = nn;
    double *construction=NULL;
    int construction_count=0;

    sfile = NULL;
    if( InvalidMeshCheck()) return -1;

    /*increment the process step counter*/
    process_step++;
    GridSave();

#   define chosen(A) (is_specified( param, A) && get_bool( param, A))
    if (chosen ("silicon")) mater = Si;
    if (chosen ("poly")) mater = Poly;
    if (chosen ("oxide")) mater = SiO2;
    if (chosen ("nitride")) mater = SiNi;
    if (chosen ("oxynitride")) mater = OxNi;
    if (chosen ("photoresist")) mater = PhRs;
    if (chosen ("aluminum")) mater = Al;
    if (chosen ("gaas")) mater = GaAs;

    thick = 1e-4 * get_float( param, "thick");
    ndiv = get_int( param, "divisions");

    /* check for a string file for deposit */
    sfile = get_string( param, "file" );

    if (sfile == NULL)  {
	if ( is_specified( param, "space" ) )
	    space = get_float( param, "space" ) * 1.0e-4;
	else
	    space = thick / ndiv;

	/* People have been known to forget vital information */
	if (thick <= 0 || space <= 0) {
	    fprintf( stderr, "Bad deposition parameters: thick %e spacing %e\n",
		thick*1e4, space*1e4);
	    return(-1);
	}

	if ( is_specified( param, "square" ) )
	    square = get_bool( param, "square" );
	else
	    square = TRUE;
    }

    /*work out which impurity to use in the doping*/
    imp = -1;
    if      ( get_bool( param, "arsenic" ))    imp = As;
    else if ( get_bool( param, "phosphorus" )) imp = P;
    else if ( get_bool( param, "antimony" ))   imp = Sb;
    else if ( get_bool( param, "boron" ))      imp = B;
    else if ( get_bool( param, "gallium" ))    imp = Ga;
    else if ( get_bool( param, "beryllium" ))	imp = iBe;
    else if ( get_bool( param, "magnesium" ))	imp = iMg;
    else if ( get_bool( param, "selenium" ))	imp = iSe;
    else if ( get_bool( param, "isilicon" ))	imp = iSi;
    else if ( get_bool( param, "tin" ))		imp = iSn;
    else if ( get_bool( param, "germanium" ))	imp = iGe;
    else if ( get_bool( param, "zinc" ))	imp = iZn;
    else if ( get_bool( param, "carbon" ))	imp = iC;
    else if ( get_bool( param, "generic" ))	imp = iG;

    /*make a new region*/
    r = mk_reg( mater );

    /*get the concentration*/
    conc = get_float( param, "concentrat" );

    /*generate the line of the surface*/
    nump = find_surf( &p );


    /*generate the offset line*/
    if (sfile)  {
        numop = string_offset( sfile, &op );
        thick = p[0].y - op[0].y;
    }
    else
        numop = gen_offset( nump, p, &op, thick, space);

#ifdef SUPREM_DOMAIN_OFFSET
    if(!sfile && mode==TWOD && getenv("SUPREM4_PROCESS_SURFACE") &&
       !strcmp(getenv("SUPREM4_PROCESS_SURFACE"),"1") &&
       getenv("SUPREM4_DEPOSITION_GEOMETRY") &&
       !strcmp(getenv("SUPREM4_DEPOSITION_GEOMETRY"),"domain")) {
        extern int suprem_nominal_offset(int,const double *,double,double,double **);
        double *input=salloc(double,2*nump);
        for(i=0;i<nump;i++){input[2*i]=p[i].x;input[2*i+1]=p[i].y;}
        construction_count=suprem_nominal_offset(nump,input,thick,space,&construction);
        free(input);
        if(construction_count<2)panic("cannot construct nominal deposition surface");
    }
#endif
    process_surface_clear();

    /*build a skeleton region and go at it*/
    iskel = build_skel( nump, p, numop, op, ndiv, r, (double)thick/ndiv, square );

    free(op);

    /*do the triangulation of the skeleton region*/
    grid( iskel, r, FALSE, TRUE );

    /*remove the original skeleton*/
    free_skel(iskel);

    /*clean up mesh*/
    bd_connect_material( "after deposit" );
    if(!sfile && deposition_spacing(r,space)<0)return -1;
    if(construction) {
        process_surface_record(mater,construction_count,construction);
        free(construction);
    }

    /*deposit new nodes with dopant*/
    if ( imp > 0 ) {
	add_impurity( imp, 1.0e5, -1 );

	switch( imp ) {
	case As  : impa = Asa;  break;
	case Sb  : impa = Sba;  break;
	case B   : impa = Ba;   break;
	case P   : impa = Pa;   break;
	case iBe : impa = iBea; break;
	case iMg : impa = iMga; break;
	case iSe : impa = iSea; break;
	case iSi : impa = iSia; break;
	case iSn : impa = iSna; break;
	case iGe : impa = iGea; break;
	case iZn : impa = iZna; break;
	case iC  : impa = iCa;  break;
	case iG  : impa = iGa;  break;
	default  : impa = 0;    break;
	}
	if ( impa && (imptosol[impa] == -1)) {
	    add_impurity(impa, 1.0, -1);
	    SET_FLAGS(impa, ACTIVE);
	}

	sol = imptosol[imp];
	for(i = 0; i < num_nd_reg(r); i++)
	    set_sol_nd( nd_reg(r,i), sol, conc );
    }


    /*make sure things are okey-dokey*/
    mtest1( "after deposit" );
    mtest2( "after deposit" );

    /*if we deposited poly, calculate the initial grain size*/
    if ( mater == Poly ) {
	temp = pressure = 0.0;

	/* get parameters */
	if ( is_specified( param, "temperature" ) )
	    temp = get_float(param, "temperature");
	if ( is_specified( param, "pressure" ) )
	    pressure = get_float(param, "pressure");

	/*add the grain size*/
        add_impurity( GRN, 0.0, -1 );

	/*inititialize the grain size*/
	poly_init( oldnn, temp, pressure, thick, imp,conc,ndiv );
    }
    return(0);
}





/************************************************************************
 *									*
 *	build_skel() - Build a skeleton out of the surface line and	*
 *  the offset outer line.						*
 *									*
 *  Original:	MEL	4/87						*
 *									*
 ************************************************************************/
build_skel( nump, p, numop, op, div, rn, sp, square )
int nump;
struct line *p;
int numop;
struct line *op;
int div;
int rn;
double sp;
int square;
{
    register int i, j;
    int ir;
    int mat = mat_reg(rn);
    int node, n1, gasnd, matnd;
    int e, ln, ip, e1, ef;
    float c[3];
    int sq;

#define min3(A,B,C) (A<B)?((B<C)?B:C):((A<C)?A:C)
#define min2(A,B) ((A<B)?(A):(B))

    /*create the initial skeleton shape*/
    ir = cr_sreg(mat);

    /*go through all the points on the surface and fix the nodes there*/
    ln = -1;
    for( i = 0; i < nump; i++ ) {
	/*find the gas node and maybe the mater node*/
	for( matnd = -1, j = 0; j < num_nd(p[i].p); j++ ) {
	    node = pt[ p[i].p ]->nd[j];
	    if ( mat_nd(node) == GAS )   gasnd = j;
	    if ( mat_nd(node) == mat ) matnd = j;
	}

	/*mark the point as non surface*/
	clr_surf( p[i].p );

	/*if we have no node of mater - change the gas node*/
	if ( matnd == -1 ) {
	    node = nd_pt(p[i].p, gasnd);
	    nd[ node ]->mater = mat;
	    p[i].p = node;
	}
	else {
	    /*destroy the gas node*/
	    n1 = nd_pt(p[i].p, gasnd);
	    fing_nd( n1 );

	    /*save the material node*/
	    node = nd_pt(p[i].p, matnd);

	    /*fix the node list up*/
	    sub_nd_pt(n1, p[i].p);

	    /*store the material node number at p*/
	    p[i].p = node;
	}

	/*build the skeleton by adding edges*/
	if ( ln != -1 ) {
	    /*find the correct edge to add...*/
	    for(ef = -1, e1 = 0; e1 < num_edge_nd(ln) && ef == -1; e1++) {
		e = edge_nd(ln, e1);
		if ( (nd_edg(e,0) == node) || (nd_edg(e,1) == node) ) ef = e;
	    }
	    if ( ef == -1 )
		ef = mk_edg(ln, node);
	    else
		clr( edg[ef], ESURF );

	    ad_edge( ir, ef, sreg[ir]->bnd, MAYBE, BEFORE );
	    set_space(pt_nd(node), len_edg(ef));
	    set_space(pt_nd(ln), min2(len_edg(ef),get_space(pt_nd(ln))));
	}
	else {
	    set_space(pt_nd(node), sp);
	}
	ln = node;
    }

    /*now go up the separator...*/
    for(i = 0; i < div-1; i++) {
	/*figure the coordinates*/
	c[0] = p[nump-1].x + (op[numop-1].x - p[nump-1].x) * (i + 1) / div;
	c[1] = p[nump-1].y + (op[numop-1].y - p[nump-1].y) * (i + 1) / div;

	if ((mode != ONED) && square && (c[0] > p[nump-1].x))
	    c[0] = p[nump-1].x;

	/*make a point and a node*/
	ip = mk_pt( 2, c );
	set_space(ip, sp);
	node = mk_nd( ip, mat );

	/*make the edge and add it*/
	e = mk_edg(ln, node);
	ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
	set_space(pt_nd(node), len_edg(e));
	set_space(pt_nd(ln), min2(len_edg(e),get_space(pt_nd(ln))));

	ln = node;
    }

    /*now walk the offset line*/
    for(i = numop-1; i >= 0; i--) {
	c[0] = op[i].x;
	c[1] = op[i].y;

	if ((mode != ONED) && square && ((c[0] - p[nump-1].x) > 1.0e-8) ) {
	    c[0] = p[nump-1].x;
	    ip = mk_pt( 2, c );
	    sq = FALSE;
	}
	else if ( (mode != ONED) && square && (p[0].x - c[0] > 1.0e-8) ) {
	    c[0] = p[0].x;
	    ip = mk_pt( 2, c );
	    sq = FALSE;
	}
	else {
	    ip = mk_pt( 2, c );
	    set_surf(ip);
	    node = mk_nd( ip, GAS );
	    sq = TRUE;
	}

	node = mk_nd( ip, mat );

	/*make the edge and add it*/
	e = mk_edg(ln, node);
	if ( (mode != ONED) && sq ) set(edg[e], ESURF);
	ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
	set_space(pt_nd(node), len_edg(e));
	set_space(pt_nd(ln), min2(len_edg(e),get_space(pt_nd(ln))));
	ln = node;
    }

    if ( mode == TWOD ) {
	/*finally back down the divider*/
	for(i = div-1; i > 0; i--) {
	    /*figure the coordinates*/
	    c[0] = p[0].x + (op[0].x - p[0].x) * i / div;
	    c[1] = p[0].y + (op[0].y - p[0].y) * i / div;

	    if ((mode != ONED) && square && (c[0] < p[0].x))
		c[0] = p[0].x;

	    /*make a point and a node*/
	    ip = mk_pt( 2, c );
	    set_space(ip, sp);
	    node = mk_nd( ip, mat );

	    /*make the edge and add it*/
	    e = mk_edg(ln, node);
	    ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
	    set_space(pt_nd(node), len_edg(e));
	    set_space(pt_nd(ln), min2(len_edg(e),get_space(pt_nd(ln))));
	    ln = node;
	}

	/*and connect the region*/
	e = mk_edg(ln, p[0].p);
	ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
    }

    /*fix up the region nice*/
    ck_clock(ir, TRUE);

    return(ir);
}



