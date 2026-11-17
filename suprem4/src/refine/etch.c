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
/*   etch.c                Version 5.1     */
/*   Last Modification : 7/3/91 15:40:40 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <global.h>
#include <constant.h>
#include <geom.h>
#include <regrid.h>
#include <material.h>
#include <impurity.h>
#include "refine.h"
#include "skel.h"
#include "process_surface.h"
#include "directional_etch.h"

/* Snapshot the directional cutter before any region is removed. Connectivity
 * of an oxide region to the ambient does not make oxide beneath a poly gate
 * visible. Cast a shadow from the interior of every non-target region. The
 * chord-midpoint construction removes exactly the same TARGET volume as
 * stopping at the first opaque surface, without making the cutter coincident
 * with material interfaces. See etch_ray_opaque_polygon for the invariant. */
static int directional_cutter(float **p,int material,int count,
    struct line *surface,double depth,double xmin,double xmax,
    double ymin,double ymax)
{
    struct etch_ray_front front;
    struct etch_ray_segment *boundary=NULL;
    struct LLedge *first,*bp;
    int r,i,j,a,b,n=0;
    front.n=1;front.s=malloc(sizeof(*front.s));
    if(!front.s)return -1;
    front.s[0].x0=xmin;front.s[0].x1=xmax;
    front.s[0].y0=front.s[0].y1=ymax+depth;
    for(i=1;i<count;i++)
        if(etch_ray_limit(&front,surface[i-1].x,surface[i-1].y+depth,
                               surface[i].x,surface[i].y+depth))goto fail;
    if(material!=-1)for(r=0;r<nreg;r++) {
        if(mat_reg(r)==material || mat_reg(r)==GAS)continue;
        boundary=malloc((size_t)sreg[reg[r]->sreg]->len*sizeof(*boundary));
        if(!boundary)goto fail;
        first=bp=sreg[reg[r]->sreg]->bnd;j=0;
        do {
            if(j>=sreg[reg[r]->sreg]->len)goto fail;
            a=pt_nd(nB(bp));b=pt_nd(nF(bp));
            boundary[j].x0=cordinate(a,0);boundary[j].y0=cordinate(a,1);
            boundary[j].x1=cordinate(b,0);boundary[j++].y1=cordinate(b,1);
            bp=bp->next;
        }while(bp!=first);
        if(etch_ray_opaque_polygon(&front,boundary,j))goto fail;
        free(boundary);boundary=NULL;
    }
    p[0]=salloc(float,2*front.n+2);p[1]=salloc(float,2*front.n+2);
    p[0][n]=xmin-(xmax-xmin);p[1][n++]=ymin-(ymax-ymin);
    for(i=0;i<front.n;i++) {
        p[0][n]=front.s[i].x0;p[1][n++]=front.s[i].y0;
        p[0][n]=front.s[i].x1;p[1][n++]=front.s[i].y1;
    }
    p[0][n]=xmax+(xmax-xmin);p[1][n++]=ymin-(ymax-ymin);
    free(front.s);return n;
fail:
    free(boundary);free(front.s);return -1;
}

static int mesh_diag_enabled()
{
    char *v = (char *)getenv("SUPREM4_MESH_DIAG");
    return(v && *v && *v != '0');
}

/* Optional local trace, in database centimeters.  Include crossing segments,
 * not only endpoints in the requested interval. */
static void etch_trace_ring(stage, region, skeleton)
char *stage;
int region, skeleton;
{
    char *window = getenv("SUPREM4_ETCH_TRACE_X");
    double lo, hi, x0, x1, y0, y1;
    struct LLedge *first, *bp;
    int a, b;
    if (!window || sscanf(window, "%lf %lf", &lo, &hi) != 2) return;
    first = bp = sreg[skeleton]->bnd;
    do {
        a = pt_nd(nB(bp)); b = pt_nd(nF(bp));
        x0 = cordinate(a,0); x1 = cordinate(b,0);
        y0 = cordinate(a,1); y1 = cordinate(b,1);
        if (fmax(x0,x1) >= lo && fmin(x0,x1) <= hi)
            fprintf(stderr, "[etchtrace] %s region=%d skel=%d edge=%d points=%d,%d xy=%.17g,%.17g:%.17g,%.17g\n",
                    stage,region,skeleton,bp->edge,a,b,x0,y0,x1,y1);
        bp = bp->next;
    } while (bp != first);
}

struct etch_interface { float c[2][2]; };

/* Refinement creates live subfaces before neighbor arrays are rebuilt. Find
 * the other owner through shared point identity, not stale parent neighbors.
 * Material interfaces can have distinct nodes at the same physical point. */
static int interface_other_owner(int t,int f,int region)
{
    int p=pt_nd(nd_face_ele(t,0,f)),q=pt_nd(nd_face_ele(t,1,f));
    int a,b,n,u,g,x,y;
    int opposite=pt_nd(vert_tri(t,f));
    int side=orient2d_float_sign(cord_arr(p),cord_arr(q),cord_arr(opposite));
    for(a=0;a<num_nd(p);a++) {
        n=nd_pt(p,a);
        for(b=0;b<num_tri_nd(n);b++) {
            u=tri_nd(n,b);
            if(dead_tri(u) || !leaf(tri[u]) || reg_tri(u)==region)continue;
            for(g=0;g<num_face(u);g++) {
                x=pt_nd(nd_face_ele(u,0,g));y=pt_nd(nd_face_ele(u,1,g));
                if((x==p && y==q)||(x==q && y==p)) {
                    int other=pt_nd(vert_tri(u,g));
                    int other_side=orient2d_float_sign(cord_arr(p),cord_arr(q),cord_arr(other));
                    /* Old and replacement triangles coexist during remeshing.
                     * Their coincident outer edges are not covered interfaces:
                     * a real neighbor lies on the opposite side of the face. */
                    if((side<0 && other_side>0)||(side>0 && other_side<0))return u;
                }
            }
        }
    }
    return -1;
}

#ifdef SUPREM4_TEST_INTERFACE
int test_interface_other_owner(int t,int f,int region)
{ return interface_other_owner(t,f,region); }
#endif

/* Snapshot interfaces with other regions after refinement, before removal.
 * Coordinates survive the database compaction performed by bd_connect. */
static int removed_interfaces(int region, struct etch_interface **interfaces)
{
    int i,t,f,neighbor,j,k,count=0,capacity=0;
    *interfaces=NULL;
    if(mode!=TWOD)return 0;
    for(i=0;i<num_tri_reg(region);i++) {
        t=tri_reg(region,i);
        if(dead_tri(t) || !leaf(tri[t]))continue;
        for(f=0;f<num_face(t);f++) {
            neighbor=interface_other_owner(t,f,region);
            if(neighbor<0)continue;
            if(count==capacity) {
                capacity=capacity ? 2*capacity:32;
                *interfaces=sralloc(struct etch_interface,capacity,*interfaces);
            }
            for(j=0;j<2;j++)for(k=0;k<2;k++)
                (*interfaces)[count].c[j][k]=cordinate(pt_nd(nd_face_ele(t,j,f)),k);
            count++;
        }
    }
    return count;
}

/* A formerly covered face is exposed only if the assembled mesh leaves it
 * without a solid neighbor. Marking every old interface before assembly can
 * mislabel still-covered domain vertices and expose an artificial sidewall. */
static int expose_removed_interfaces(struct etch_interface *interfaces, int count)
{
    int t=0,f,i,p,q,changed=FALSE;
    double x0,y0,x1,y1;
    while(done_tri(t)) {
        for(f=0;f<num_face(t);f++) {
            if(neigh_fc(t,f)!=BC_OFFSET)continue;
            p=pt_nd(nd_face_ele(t,0,f));q=pt_nd(nd_face_ele(t,1,f));
            x0=cordinate(p,0);y0=cordinate(p,1);
            x1=cordinate(q,0);y1=cordinate(q,1);
            for(i=0;i<count;i++) {
                float (*c)[2]=interfaces[i].c;
                if((x0==c[0][0] && y0==c[0][1] && x1==c[1][0] && y1==c[1][1]) ||
                   (x0==c[1][0] && y0==c[1][1] && x1==c[0][0] && y1==c[0][1])) {
                    mk_surf(p);mk_surf(q);set_face_exp(t,f);changed=TRUE;break;
                }
            }
        }
        next_tri(t);
    }
    return changed;
}

static void mesh_diag_poly(name, n, p)
char *name;
int n;
float *p[MAXDIM+1];
{
    int i;
    double dx, dy, len, minlen = MAXFLOAT;

    if (!mesh_diag_enabled()) return;
    for(i = 1; i < n; i++) {
	dx = p[0][i] - p[0][i-1];
	dy = (mode == ONED)?0.0:(p[1][i] - p[1][i-1]);
	len = sqrt(dx * dx + dy * dy);
	if (len < minlen) minlen = len;
    }
    if (n < 2) minlen = 0.0;
    fprintf(stderr, "[meshdiag] %s points=%d min_segment=%e\n", name, n, minlen);
}

/*some useful macros for getting at the parameters*/
#define CHOSEN(A) (is_specified( param, A) && get_bool( param, A))
#define FETCH(N,V,S) V = S*get_float( param, N)
#define FETCHD(N,V,S,D) if(is_specified(param,N)) V = S*get_float( param, N);\
			else V = D;

#define MOVE( A, B, C ) {pt[A]->cordo[0] = pt[A]->cord[0]=B;\
			 pt[A]->cordo[1] = pt[A]->cord[1]=C;\
			 inside[A]=TRUE;}

#define DIS(p1,p2) hypot(pt[p1]->cord[0]-pt[p2]->cord[0], pt[p1]->cord[1]-pt[p2]->cord[1])

/*A macro for point in the etch region*/
#define IN_ETCH( AX, AY, B ) ptnreg( (AX), (AY), (B)->num, (B)->cord )

/*-----------------USER_ETCH--------------------------------------------
 * User interface to etching routine.
 *----------------------------------------------------------------------*/
user_etch( par, param)
    char *par;
    struct par_str *param;
{
    register int i;
    static int mater_save;
    int mater;
    int f;
    float x;
    float xma, xmi, yma, ymi;
    double thick;
    double rate[MAXMAT], time;
    static float *p[MAXDIM+1];	/*list of points on the etch boundary*/
    static int nump;
    static int maxp;
    struct line *surf;
    int nums;
    char *sfile;
    extern string_etch();

    sfile = NULL;
    if( InvalidMeshCheck()) return -1;
    GridSave();

    /* check for a string file for deposit */
    sfile = get_string( param, "file" );

    mater = -1;
    if (CHOSEN ("silicon")) mater = Si;
    if (CHOSEN ("poly")) mater = Poly;
    if (CHOSEN ("oxide")) mater = SiO2;
    if (CHOSEN ("oxynitride")) mater = OxNi;
    if (CHOSEN ("nitride")) mater = SiNi;
    if (CHOSEN ("photoresist")) mater = PhRs;
    if (CHOSEN ("aluminum")) mater = Al;
    if (CHOSEN ("gaas")) mater = GaAs;

    /*get the limits of the given material*/
    xmi = xma = cordinate(pt_nd(0),0);
    ymi = yma = cordinate(pt_nd(0),1);
    for(i = 1; i < nn; i++) {
	if ( xmi > cordinate(pt_nd(i),0) ) xmi = cordinate(pt_nd(i),0);
	if ( xma < cordinate(pt_nd(i),0) ) xma = cordinate(pt_nd(i),0);
	if ( ymi > cordinate(pt_nd(i),1) ) ymi = cordinate(pt_nd(i),1);
	if ( yma < cordinate(pt_nd(i),1) ) yma = cordinate(pt_nd(i),1);
    }

    if ((CHOSEN("continue") || CHOSEN("done")) && !p[0]) {
        fprintf(stderr,"Etch polygon requires a start command.\n");
        return -1;
    }
    if (p[0] && !(CHOSEN("continue") || CHOSEN("done"))) {
        fprintf(stderr,"Finish the current etch polygon before starting another etch.\n");
        return -1;
    }

    /*work out the type of the etch specified*/
    if ( CHOSEN("left") || CHOSEN("right")) {

	/* 1-D dosen't have this option */
        if ( mode != ONED ) {

	    p[0] = salloc(float, 4);
	    p[1] = salloc(float, 4);
	    p[3] = salloc(float, 4);
	    /*get the endpoints of the line*/
	    FETCH ( "p1.x", p[0][0], 1e-4 );
	    FETCHD( "p1.y", p[1][0], 1e-4 , 2*ymi - yma);
	    FETCHD( "p2.x", p[0][1], 1e-4 , p[0][0]);
	    FETCHD( "p2.y", p[1][1], 1e-4 , 2*yma - ymi);
	    FETCHD( "spacing", p[3][0], 1e-4, -1.0 );
	    FETCHD( "spacing", p[3][1], 1e-4, -1.0 );

	    /*figure out which way to go*/
	    if ( get_bool( param, "left" ) )
		x = 2.0 * xmi - xma;
	    else if ( get_bool( param, "right" ) )
		x = 2.0 * xma - xmi;

	    /*set up the final sides*/
	    p[0][2] = x;
	    p[1][2] = p[1][1];
	    p[3][2] = xma - xmi;
	    p[0][3] = x;
	    p[1][3] = p[1][0];
	    p[3][3] = xma - xmi;

	    nump = 4;
	    etch( p, nump, mater );

	    sfree(p[0]);
	    sfree(p[1]);
            sfree(p[3]);
            p[0]=p[1]=p[3]=NULL; nump=maxp=0;
	}
	else {
	    fprintf(stderr, "\n\? NO [left|right] options in 1-D! \n\n");
	}
    }
    else if ( CHOSEN( "dry" ) ) {

        const double *construction=NULL;
        int construction_count=0;
	/*how much do we take off*/
	thick = get_float( param, "thick" ) * 1.0e-4;

	/*find the surface line*/
	nums = find_surf( &surf );
        if(getenv("SUPREM4_PROCESS_SURFACE") && !strcmp(getenv("SUPREM4_PROCESS_SURFACE"),"1"))
            construction_count=process_surface_get(mater,nums,surf,&construction);
        if(construction_count) {
            free(surf);nums=construction_count;surf=salloc(struct line,nums);
            for(i=0;i<nums;i++){surf[i].x=construction[2*i];surf[i].y=construction[2*i+1];}
            fprintf(stderr,"[process-surface] dry etch uses construction vertices=%d\n",nums);
        }
        process_surface_clear();

	p[3] = NULL; /* Dry etches have no per-vertex spacing array. */
	if(mode==TWOD) {
            nump=directional_cutter(p,mater,nums,surf,thick,xmi,xma,ymi,yma);
            if(nump<0){free(surf);fprintf(stderr,"Cannot allocate directional etch front.\n");return -1;}
        }
        else {
	p[0] = salloc( float, nums+10 );
	p[1] = salloc( float, nums+10 );

	/*build the etch region*/
	p[0][0] = xmi - (xma - xmi);
	if ( mode == TWOD ) p[1][0] = ymi - (yma - ymi);

	/*now add the offset surface line*/
	for(i = 0; i < nums; i++) {
	   if ( mode == ONED ) p[0][i+1] = surf[i].x + thick;
	   else {
		p[0][i+1] = surf[i].x;
		p[1][i+1] = surf[i].y + thick;
	   }
	}

	/*add the last to complete the rectangle*/
	if ( mode == ONED ) {
	   nump = nums+1;
	}
	else if ( mode == TWOD ) {
	   p[0][i+1] = xma + (xma - xmi);
	   p[1][i+1] = ymi - (yma - ymi);
	   nump = nums+2;
	}
        }

	process_step++;

	etch( p, nump, mater );

        free(surf);
	sfree(p[0]);
	sfree(p[1]);
        p[0]=p[1]=NULL; nump=maxp=0;

    }
    else if ( CHOSEN( "start" ) ) {

	mater_save = -1;
        if ( mode != ONED ) {
	    /*initialize the etch region*/
	    maxp = 500;
	    p[0] = salloc( float, maxp );
	    p[1] = salloc( float, maxp );
	    p[3] = salloc( float, maxp );
	    nump = 1;
	    FETCH( "x", p[0][0], 1e-4 );
	    FETCH( "y", p[1][0], 1e-4 );
	    FETCHD( "spac", p[3][0], 1e-4, -1.0 );
	}
	else {
	    fprintf(stderr,"\n\?NO [start|continue|done] options in 1-D!\n");
	}
	if ( mater != -1 ) mater_save = mater;
    }
    else if ( CHOSEN( "continue" ) ) {

        if ( mode != ONED ) {

	    if ( nump >= maxp ) {
                maxp += 500;
		p[0] = sralloc( float, maxp, p[0] );
		p[1] = sralloc( float, maxp, p[1] );
		p[3] = sralloc( float, maxp, p[3] );
	    }
	    /*continue defining the etch region*/
	    FETCH( "x", p[0][nump], 1e-4 );
	    FETCH( "y", p[1][nump], 1e-4 );
	    FETCHD( "spac", p[3][nump], 1e-4, -1.0 );
	    nump++;
	}
        else {
	    fprintf(stderr,"\n\?NO [start|continue|done] options in 1-D!\n");
	}
	if ( mater != -1 ) mater_save = mater;
    }
    else if ( CHOSEN( "done" ) ) {

        if ( mode != ONED ) {
	    if ( mater != -1 ) mater_save = mater;
	    mater = mater_save;
	    mater_save = -1;

	    /*finish defining the etch region*/
	    if ( nump >= maxp ) {
		maxp += 500;
		p[0] = sralloc( float, maxp, p[0] );
		p[1] = sralloc( float, maxp, p[1] );
                p[3] = sralloc( float, maxp, p[3] );
	    }
	    FETCH( "x", p[0][nump], 1e-4 );
	    FETCH( "y", p[1][nump], 1e-4 );
	    FETCHD( "spac", p[3][nump], 1e-4, MAXFLOAT );
	    nump++;
	    etch( p, nump, mater );
	    free(p[0]);
	    free(p[1]);
            free(p[3]);
            p[0]=p[1]=p[3]=NULL; nump=maxp=0;
	}
	else {
	    fprintf(stderr,"\n\?NO [start|continue|done] options in 1-D!\n");
	}

    }
    else if ( CHOSEN( "all" ) ) {
	for(i = 0; i < nreg; i++)
	    if ( mat_reg(i) == mater ) rem_reg(i);
	bd_connect("after etch");
    }
    else if ( CHOSEN( "physical" ) ) {
	for( i = 0; i < MAXMAT; i++ ) rate[i] = 0.0;
	FETCHD( "r.silicon", rate[Si], 1.0e-4, 0.0 );
	FETCHD( "r.oxide", rate[SiO2], 1.0e-4, 0.0 );
	FETCHD( "r.oxynitride", rate[OxNi], 1.0e-4, 0.0 );
	FETCHD( "r.nitride", rate[SiNi], 1.0e-4, 0.0 );
	FETCHD( "r.poly", rate[Poly], 1.0e-4, 0.0 );
	FETCHD( "r.photoresist", rate[PhRs], 1.0e-4, 0.0 );
	FETCHD( "r.aluminum", rate[Al], 1.0e-4, 0.0 );
	FETCHD( "r.gaas", rate[GaAs], 1.0e-4, 0.0 );
	FETCH( "time", time, 1.0 );

	phys_etch( rate, time );
    }
    else if (sfile)
	string_etch(sfile, mater);

    return(0);
}



/*-----------------ETCH-------------------------------------------------
 *----------------------------------------------------------------------*/
etch( p, nump, mater )
float *p[MAXDIM+1];
int nump;
int mater;
{
    int i, e;
    int ir, node, ip, ln, fn;
    int sr;
    int r, rn, s[MAXREG];
    int ns;
    float c[MAXDIM];
    struct LLedge *bp;
    int f, n, ninterfaces=0;
    struct etch_interface *interfaces=NULL;

    /* Projection into the mesh's coordinate type can identify consecutive
     * construction vertices. Canonicalize those exact identities before
     * making edges; retain every nonzero edge, however short. */
    if (mode == TWOD) {
        int kept=0,j;
        for(j=0;j<nump;j++) {
            if(kept && p[0][j]==p[0][kept-1] && p[1][j]==p[1][kept-1]) {
                if(p[3] && p[3][j]>0 &&
                   (p[3][kept-1]<=0 || p[3][j]<p[3][kept-1]))
                    p[3][kept-1]=p[3][j];
                continue;
            }
            p[0][kept]=p[0][j];p[1][kept]=p[1][j];
            if(p[3])p[3][kept]=p[3][j];
            kept++;
        }
        if(kept>1 && p[0][kept-1]==p[0][0] && p[1][kept-1]==p[1][0]) {
            if(p[3] && p[3][kept-1]>0 && (p[3][0]<=0 || p[3][kept-1]<p[3][0]))
                p[3][0]=p[3][kept-1];
            kept--;
        }
        nump=kept;
        if(nump<3) {
            fprintf(stderr,"Etch polygon has fewer than three distinct vertices.\n");
            return -1;
        }
    }
    mesh_diag_poly("etch_polygon", nump, p);
    for(r = 0; r < nreg; r++) clr(reg[r], ETCHED);

    /*for every material that we are worried about...*/
    for(r = 0; r < nreg; r++) {
	if ( (mater != -1) && (mat_reg(r) != mater) ) continue;

	/*only etch exposed regions*/
	if ( !ask( reg[r], EXPOS ) ) continue;

	if ( ask(reg[r], ETCHED) ) continue;
	set(reg[r], ETCHED);

	/*build a skeleton region out of the list of points*/
	ir = cr_sreg(-1);

	/*loop through all the points that were passed*/
	ln = -1;
	for(i = 0; i < nump; i++) {
	    c[0] = p[0][i];
	    c[1] = p[1][i];

	    /*make a point and a node*/
	    ip = mk_pt( 2, c );
	    if ( p[3] != NULL && (p[3][i] > 0.0) )
		set_space(ip, p[3][i]);
	    else
		set_space(ip, MAXFLOAT);
	    node = mk_nd( ip, mat_reg(r) );
	    set(pt[ip], SKELP);
	    if (i == 0) fn = node;

	    /*make the edge and add it*/
	    if ( ln != -1 ) {
		e = mk_edg(ln, node);
		ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
	    }
	    ln = node;
	}
	e = mk_edg(ln, fn);
	if (mode != ONED ) ad_edge( ir, e, sreg[ir]->bnd, MAYBE, BEFORE );
#ifdef DEBUG
	reg_pl(ir);
#endif
#define min2(A,B) (A < B)?A:B
	ck_clock(ir, TRUE);
	if (mesh_diag_enabled())
	    fprintf(stderr, "[meshdiag] etch_skeleton region=%d target_region=%d len=%d material=%d\n",
		    ir, r, sreg[ir]->len, mat_reg(r));
	for(f=1, bp=sreg[ir]->bnd; (bp != sreg[ir]->bnd)||f; f=0, bp=bp->next) {
	    n = nB(bp);
	    ip = pt_nd(n);
	    for(i = 0; i < num_edge_nd(n); i++) {
		set_space(ip,min2(get_space(ip),len_edg(edge_nd(n,i))));
	    }
	}

	/*build a skeleton from the material we are worried about*/
	sr = reg[r]->sreg;;
	hinit(sr, MAXFLOAT);
	ck_clock(sr, TRUE);
	etch_trace_ring("subject-before",r,sr);
	etch_trace_ring("clip-before",r,ir);
	ns = sub_skel(sr, ir, s);
	if (mesh_diag_enabled())
	    fprintf(stderr, "[meshdiag] sub_skel target_region=%d target_skel=%d etch_skel=%d outputs=%d\n",
		    r, sr, ir, ns);
	if (ns == -1) {
	    rem_reg(r);
	}
	else if ( ns == 0) {
	}
	else {

	    for(i = 0; i < ns; i++) {
		etch_trace_ring("retained-before-grid",r,s[i]);
#ifdef DEBUG
		reg_pl(s[i]);
#endif
		rn = mk_reg( mat_reg(r) );
		reg[rn]->thick = reg[r]->thick;
		set(reg[rn], ETCHED);
		hinit( s[i], MAXFLOAT );
		grid( s[i], rn, TRUE, TRUE );
		etch_trace_ring("retained-after-grid",rn,s[i]);
		interp( rn, r );
	    }
	    /* ON interface points can become exposed when their cover disappears.
	     * The assembled rebuild clears gas where retained fragments cover them. */
	    ninterfaces=removed_interfaces(r,&interfaces);
	    while(num_tri_reg(r) > 0) rem_ele(tri_reg(r,0), FALSE);
	    fing_reg(r);
	}

	/*clear the skeleton etch region*/
	for(i = 0; i < ne; i++) clr( tri[i], NEIGHB );

	/*mark as dead all the points in the skeleton*/
	for(i = 1, bp = sreg[ir]->bnd;
	    i || bp != sreg[ir]->bnd; i = 0, bp = bp->next) {
	    if (num_tri_nd(nB(bp)) == 0) rem_1nd(nB(bp));
	}
	if ( mode == ONED ) {
	    bp = bp->prev;
	    if (num_tri_nd(nF(bp)) == 0) rem_1nd(nF(bp));
	}

	/*this step renumbers the regions - have to start over...*/
	if ( ns != 0 ) r = -1;
	bd_connect_material("after etch of a single region");
        if(ninterfaces && expose_removed_interfaces(interfaces,ninterfaces)) {
            neigh_dirty=TRUE;
            bd_connect_material("after exposing uncovered interfaces");
        }
        free(interfaces);interfaces=NULL;ninterfaces=0;
    }
    /*clear the temporary skeleton marker*/
    for(i = 0; i < np; i++) clr(pt[i], SKELP);
}



