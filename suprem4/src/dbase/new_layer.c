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
/*   new_layer.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:21:06 */

#include <stdio.h>
#include "suprem_options.h"
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "material.h"
#include "impurity.h"
#include "mesh_trial.h"
#ifdef CONVEX
#define hypot(a,b) sqrt((a)*(a)+(b)*(b))
#endif

#define X 0
#define Y 1
#define Z 2
#define NEUTRAL BC_OFFSET

#define vertex(t,i) tri[t]->nd[i]
#define dot(A,B) ((mode == ONED)?((A)[0]*(B)[0]):((A)[0]*(B)[0]+(A)[1]*(B)[1]))



/*-----------------NativeOxide------------------------------------------
 * Deposits oxide on bare silicon surfaces.
 * Like new_layer but different.
 *
 * Return value: whether geometry changed
 * Diagnostics:  abort on error
 *
 * notes:
 * Uses tri->nb, nd->tri
 * Insinuates oxide under nitride if nitride on bare silicon.
 * Restores all geometry when done.
 *----------------------------------------------------------------------*/
static int native_covered(int node, int silicon)
{
    int i, material;
    for (i=0;i<num_nd(pt_nd(node));i++) {
        material=mat_nd(nd_pt(pt_nd(node),i));
        if (material!=silicon && material!=GAS) return TRUE;
    }
    return FALSE;
}

/* Build a positive-width cap by following the old oxide free boundary
 * to its intersection with the exposed-run offset. Old boundary vertices
 * stay fixed; a new film fills the gas-side polygon. No zero-width taper. */
struct native_cap {int *node,count,exit;};
static void native_cap_append(struct native_cap *cap,int node)
{
    int *next;
    if(cap->count==INT_MAX || (size_t)cap->count+1>((size_t)-1)/sizeof(int))
        panic("native oxide cap allocation overflow");
    next=realloc(cap->node,((size_t)cap->count+1)*sizeof(int));
    if(!next)panic("native oxide cap allocation failed");
    cap->node=next;cap->node[cap->count++]=node;
}
/* The old oxide free boundary is followed from the joined endpoint while its
 * offset-normal projection strictly increases. If it reaches the new film's
 * offset line, that edge is noded at the crossing (bird's-beak style join).
 * If the boundary first turns away from the film (spacer top, flush or
 * recessed oxide) or ends, the film receives a positive-width end face at the
 * last rising vertex M: a new gas-side vertex X=M+(thick-proj(M))*n. The old
 * edge leaving M is recorded in cap->exit; native_cap_fill certifies with
 * exact predicates that it does not enter the cap polygon. */
static int native_cap_end_face(int current,double current_projection,double thick,
                               double nx,double ny,struct native_cap *cap)
{
    float c[MAXDIM];
    int k,point,node,i;
    double rise=thick-current_projection;
    if(!(rise>0))panic("native oxide cap end face has no height");
    for(k=0;k<MAXDIM;k++)c[k]=cordinate(pt_nd(current),k);
    c[X]=(double)c[X]+rise*nx;c[Y]=(double)c[Y]+rise*ny;
    if(c[X]==xcord(pt_nd(current)) && c[Y]==ycord(pt_nd(current)))
        panic("native oxide cap end face is not representable");
    /* A gas-exposed point owns a GAS node, as sp_edge does for ESURF splits. */
    point=mk_pt(3,c);
    set_space(point,get_space(pt_nd(current)));
    set_surf(point);
    (void)mk_nd(point,GAS);
    node=mk_nd(point,SiO2);
    for(i=0;i<MAXIMP;i++)nd[node]->sol[i]=nd[current]->sol[i];
    native_cap_append(cap,node);
    return node;
}
static int native_join_cap(int silicon,struct LLedge *bare,double thick,
                           struct native_cap *cap)
{
    int oxide=node_mat(silicon,SiO2),p=pt_nd(silicon),i,e,q,other,current,previous=-1;
    int chosen,newpoint,t,j,k,a,b,child,walk;
    double dx=xcord(pt_nd(nF(bare)))-xcord(pt_nd(nB(bare)));
    double dy=ycord(pt_nd(nF(bare)))-ycord(pt_nd(nB(bare)));
    double length=hypot(dx,dy),nx,ny,projection,current_projection=0,fraction;
    float c[MAXDIM],v[3][MAXDIM],saved[MAXDIM];
    if(oxide<0 || length<=0)panic("invalid native oxide join");
    nx=dy/length;ny=-dx/length;current=oxide;native_cap_append(cap,current);cap->exit=-1;
    for(walk=0;walk<ned;walk++) {
        chosen=-1;
        for(i=0;i<num_edge_nd(current);i++) {
            e=edge_nd(current,i);if(e==previous || !ask(edg[e],ESURF))continue;
            if(chosen>=0)panic("native oxide cap has ambiguous free boundary");
            chosen=e;
        }
        /* Oxide free boundary ends (meets a covering material): end face. */
        if(chosen<0)return native_cap_end_face(current,current_projection,thick,nx,ny,cap);
        other=nd_edg(chosen,nd_edg(chosen,0)==current?1:0);q=pt_nd(other);
        projection=((double)xcord(q)-xcord(p))*nx+((double)ycord(q)-ycord(p))*ny;
        if(getenv("SUPREM4_MESH_DIAG"))fprintf(stderr,"[native-cap-walk] from=(%.6g,%.6g) to=(%.6g,%.6g) mat=%d projection=%.6g thick=%.6g\n",xcord(pt_nd(current))*1e4,ycord(pt_nd(current))*1e4,xcord(q)*1e4,ycord(q)*1e4,mat_nd(other),projection*1e4,thick*1e4);
        if(projection>=thick)break;
        /* The film joins old oxide only along a strictly rising surface
         * (bird's beak). A flat or falling surface (flush oxide, spacer
         * top) ends the film at its last rising vertex instead of coating
         * the old oxide's top. */
        if(projection<=current_projection) {
            cap->exit=other;
            return native_cap_end_face(current,current_projection,thick,nx,ny,cap);
        }
        if(node_mat(other,mat_nd(silicon))>=0)
            panic("native oxide cap crosses another exposed run; global union required");
        native_cap_append(cap,other);previous=chosen;current=other;current_projection=projection;
    }
    if(walk==ned)panic("native oxide cap boundary does not intersect offset");
    current_projection=((double)xcord(pt_nd(current))-xcord(p))*nx+
                       ((double)ycord(pt_nd(current))-ycord(p))*ny;
    fraction=(thick-current_projection)/(projection-current_projection);
    for(k=0;k<MAXDIM;k++)c[k]=cordinate(pt_nd(current),k)+fraction*
                                      ((double)cordinate(q,k)-cordinate(pt_nd(current),k));
    for(i=0;i<2;i++) {
        for(k=0;k<mode;k++)if(c[k]!=cordinate(pt_edg(chosen,i),k))break;
        if(k==mode) {
            if(pt_edg(chosen,i)==p)panic("native oxide cap cannot be represented");
            native_cap_append(cap,nd_edg(chosen,i));return cap->node[cap->count-1];
        }
    }
    for(i=0;i<num_tri_edg(chosen);i++) {
        t=tri_edg(chosen,i);if(dead_tri(t) || !leaf(tri[t]))continue;
        a=b=-1;
        for(j=0;j<3;j++) {
            if(vert_tri(t,j)==nd_edg(chosen,0))a=j;
            if(vert_tri(t,j)==nd_edg(chosen,1))b=j;
            for(k=0;k<MAXDIM;k++)v[j][k]=cordinate(pt_nd(vert_tri(t,j)),k);
        }
        if(a<0 || b<0)panic("native oxide cap edge is inconsistent");
        for(child=0;child<2;child++) {
            j=child?b:a;
            for(k=0;k<MAXDIM;k++){saved[k]=v[j][k];v[j][k]=c[k];}
            if(orient2d_float_sign(v[0],v[1],v[2])<=0)panic("native oxide cap has an invalid rounded child");
            for(k=0;k<MAXDIM;k++)v[j][k]=saved[k];
        }
    }
    sp_gredg(chosen,c,&newpoint);
    native_cap_append(cap,node_mat(nd_pt(newpoint,0),SiO2));
    return cap->node[cap->count-1];
}
static int native_cap_orient(int a,int b,int c)
{
    return orient2d_float_sign(cord_arr(pt_nd(a)),cord_arr(pt_nd(b)),cord_arr(pt_nd(c)));
}
static int native_cap_on_segment(int p,int a,int b)
{
    return native_cap_orient(a,b,p)==0 &&
        xcord(pt_nd(p))>=fmin(xcord(pt_nd(a)),xcord(pt_nd(b))) &&
        xcord(pt_nd(p))<=fmax(xcord(pt_nd(a)),xcord(pt_nd(b))) &&
        ycord(pt_nd(p))>=fmin(ycord(pt_nd(a)),ycord(pt_nd(b))) &&
        ycord(pt_nd(p))<=fmax(ycord(pt_nd(a)),ycord(pt_nd(b)));
}
static int native_cap_bad_segment(int a,int b,int c,int d)
{
    int ac=native_cap_orient(a,b,c),ad=native_cap_orient(a,b,d);
    int ca=native_cap_orient(c,d,a),cb=native_cap_orient(c,d,b);
    if(ac*ad<0 && ca*cb<0)return 1;
    return native_cap_on_segment(a,c,d) || native_cap_on_segment(b,c,d) ||
           native_cap_on_segment(c,a,b) || native_cap_on_segment(d,a,b);
}
/* Triangulate the complete cap polygon, retaining every old free-boundary
 * subdivision. Prepare all ears before adding any cap element. */
static void native_cap_fill(struct native_cap *cap,int apex,int region)
{
    int n,*poly,*order,*triangles,i,j,k,count,out=0,a,b,c,ear;
    if(cap->count<2 || cap->count==INT_MAX ||
       (size_t)cap->count+1>((size_t)-1)/(3*sizeof(int)))panic("native oxide cap allocation overflow");
    n=cap->count+1;count=n;
    poly=malloc((size_t)n*sizeof(int));order=malloc((size_t)n*sizeof(int));
    triangles=malloc((size_t)3*(n-2)*sizeof(int));
    int bc[3]={BC_OFFSET,BC_OFFSET,BC_OFFSET};long double area=0;
    double x0=xcord(pt_nd(apex)),y0=ycord(pt_nd(apex));
    if(!poly || !order || !triangles)panic("native oxide cap allocation failed");
    poly[0]=apex;for(i=0;i<cap->count;i++)poly[i+1]=cap->node[i];
    for(i=0;i<n;i++) {
        j=(i+1)%n;
        area+=((long double)xcord(pt_nd(poly[i]))-x0)*(ycord(pt_nd(poly[j]))-y0)-
              ((long double)ycord(pt_nd(poly[i]))-y0)*(xcord(pt_nd(poly[j]))-x0);
        for(k=0;k<i;k++)if(pt_nd(poly[k])==pt_nd(poly[i]))panic("native oxide cap repeats a vertex");
    }
    for(i=0;i<n;i++)for(j=i+2;j<n;j++) {
        if(i==0 && j==n-1)continue;
        if(native_cap_bad_segment(poly[i],poly[(i+1)%n],poly[j],poly[(j+1)%n]))
            panic("native oxide cap boundary requires a global union");
    }
    if(area==0)panic("native oxide cap has zero area");
    /* End-face cap: the old free boundary leaves M=poly[n-2] along
     * (M,exit). Certify that this edge lies strictly outside the cap's
     * interior angle at M and crosses no other cap edge, so the new film
     * cannot overlap the old oxide locally. */
    if(cap->exit>=0) {
        int s=area>0?1:-1,m=poly[n-2],pa=poly[n-3],pc=poly[n-1],x=cap->exit,outside;
        int turn=s*native_cap_orient(pa,m,pc),o1=s*native_cap_orient(pa,m,x),o2=s*native_cap_orient(m,pc,x);
        outside=turn>0?(o1<0 || o2<0):(o1<0 && o2<0);
        if(!outside)panic("native oxide cap end face overlaps the old oxide boundary");
        for(i=0;i<n;i++) {
            a=poly[i];b=poly[(i+1)%n];if(a==m || b==m)continue;
            if(native_cap_bad_segment(m,x,a,b))
                panic("native oxide cap end face crosses the old oxide boundary");
        }
    }
    for(i=0;i<n;i++)order[i]=area>0?i:n-1-i;
    while(count>3) {
        ear=-1;
        for(i=0;i<count;i++) {
            a=poly[order[(i+count-1)%count]];b=poly[order[i]];c=poly[order[(i+1)%count]];
            if(native_cap_orient(a,b,c)<=0)continue;
            for(j=0;j<count;j++) {
                k=poly[order[j]];if(k==a || k==b || k==c)continue;
                if(native_cap_orient(a,b,k)>=0 && native_cap_orient(b,c,k)>=0 && native_cap_orient(c,a,k)>=0)break;
            }
            if(j==count){ear=i;break;}
        }
        if(ear<0)panic("native oxide cap polygon needs a global union");
        triangles[3*out]=a;triangles[3*out+1]=b;triangles[3*out+2]=c;out++;
        for(j=ear;j<count-1;j++)order[j]=order[j+1];count--;
    }
    for(i=0;i<3;i++)triangles[3*out+i]=poly[order[i]];
    if(native_cap_orient(triangles[3*out],triangles[3*out+1],triangles[3*out+2])<=0)
        panic("native oxide cap has a degenerate final ear");
    out++;
    /* Every boundary segment has one owner; every internal diagonal has
     * two opposite owners. Together with positive ears this certifies the
     * prepared triangulation's attachment to the unchanged old boundary. */
    for(i=0;i<out;i++)for(j=0;j<3;j++) {
        int forward=0,reverse=0,boundary=0,ti,f,u=triangles[3*i+j],v=triangles[3*i+(j+1)%3];
        for(ti=0;ti<out;ti++)for(f=0;f<3;f++) {
            a=triangles[3*ti+f];b=triangles[3*ti+(f+1)%3];
            forward+=a==u && b==v;reverse+=a==v && b==u;
        }
        for(k=0;k<n;k++) {
            a=poly[k];b=poly[(k+1)%n];if(area<0){c=a;a=b;b=c;}
            if(a==u && b==v)boundary=1;
        }
        if(forward!=1 || reverse!=(boundary?0:1))panic("native oxide cap attachment is invalid");
    }
    for(i=0;i<out;i++)mk_ele(3,&triangles[3*i],3,bc,region,FALSE);
    free(poly);free(order);free(triangles);
}

NativeOxide(thick)
    double thick;		/* That's me */
{
    int in;			/* Silicon node to deposit on */
    int addedOx=0;		/* whether anything happened */
    int ir, sr, r;
    float Normal[MAXDIM];		/* offset for new point */
    int v[MAXVRT];
    struct LLedge *bnd, *beg, *ende, *bp;
    int g;
    int tn;
    int bo, to;
    int join_existing, first_join, last_join;
    int closed_run, run_edges, run_vertex;
    int bsp, tsp;
    const char *cap_option=suprem_option("SUPREM4_NATIVE_JOIN_CAP");
    int cap=cap_option && !strcmp(cap_option,"1"),first_cap=-1,last_cap=-1;
    struct native_cap first_chain={0},last_chain={0};
    int b[3];
    b[1] = b[2] = b[0] = BC_OFFSET;

    switch( mode ) {

    case ONED :
	/*find the exposed node*/
	for( in = 0; in < nn; in++) {
	    if( nd[in]->mater != Si) continue;
	    if( !exposed( in)) continue;

	    /*create an oxide node*/
	    r = mk_reg(SiO2);
	    addedOx++;
	    bo = mk_nd( nd[ in]->pt, SiO2);

	    Normal[X] = - thick;
	    bsp = SplitNode( bo, Normal, Si);
	    v[0] = bo;
	    v[1] = bsp;
	    (void)mk_ele_nd( 2, v, r);
	}
	break;

    case TWOD:
restart_native_scan:
	for(ir = 0; ir < nreg; ir++) {
	    if ( (mat_reg(ir) != Si) && (mat_reg(ir) != Poly) ) continue;
	    sr = reg[ir]->sreg;

	    /*loop through all the edges in this region looking for exposed*/
	    bnd = sreg[sr]->bnd;
	    ende=NULL;
	    for(g = 1, beg = bnd; g || beg != bnd; g=0, beg = beg->next) {
		if ( ! ask(edg[beg->edge], ESURF) ) continue;
		fprintf(stderr, "[DIAG] NativeOxide: ESURF edge at y=%.3f um (region %d)\n",
		        cordinate(pt_nd(nF(beg)), 1) * 1e4, ir);

                /* A skeleton's list head is arbitrary. Start at the actual
                 * beginning of the exposed run, even when it wraps the head.
                 * Otherwise its initial edges are silently left uncoated. */
                bp=beg;
                while(ask(edg[beg->prev->edge],ESURF)) {
                    beg=beg->prev;
                    if(beg==bp)break;
                }
                closed_run=ask(edg[beg->prev->edge],ESURF);
                run_edges=0;ende=beg;
                do {
                    run_edges++;ende=ende->next;
                } while(ende!=beg && ask(edg[ende->edge],ESURF));

		first_join = !closed_run && node_mat(nB(beg),SiO2) >= 0;
		last_join = !closed_run && node_mat(nB(ende),SiO2) >= 0;
		/* A one-edge bare gap with oxide at both ends needs an interior
		 * seed vertex. Tapering both endpoints alone would add no oxide
		 * and leave the same exposed edge to be processed forever. */
		if (!cap && first_join && last_join && beg->next == ende) {
		    float midpoint[MAXDIM];
		    int point, k;
		    for(k=0;k<mode;k++)
		        midpoint[k]=.5*((double)cordinate(pt_nd(nB(beg)),k)+
		                               cordinate(pt_nd(nF(beg)),k));
		    (void)sp_gredg(beg->edge,midpoint,&point);
		    /* Splitting leaves face-neighbor data dirty. Rebuild before
		     * computing the new vertex's interface normal. */
		    bd_connect("native oxide gap subdivision");
		    goto restart_native_scan;
		}
                if(cap) {
                    first_cap=first_join?native_join_cap(nB(beg),beg,thick,&first_chain):-1;
                    last_cap=last_join?native_join_cap(nB(ende),ende->prev,thick,&last_chain):-1;
                }
		r = mk_reg(SiO2);

		/*for every edge...*/
		for(bp = beg->prev,run_vertex=0;
                    run_vertex<run_edges+(closed_run ? 0:1);
                    bp=bp->next,run_vertex++) {
		    tn = nF(bp);
		    
                    /* A capped join retains the inner material vertex and
                     * attaches to the noded old free boundary. The legacy
                     * endpoint taper remains available for comparison. */
                    if(getenv("SUPREM4_MESH_DIAG"))fprintf(stderr,"[native-seed] node=%d point=%d x=%.17g y=%.17g begin=%d end=%d joins=%d,%d oxide=%d\n",tn,pt_nd(tn),(double)cordinate(pt_nd(tn),0),(double)cordinate(pt_nd(tn),1),nB(beg),nB(ende),first_join,last_join,node_mat(tn,SiO2));
		    to = node_mat(tn,SiO2);
		    join_existing = (first_join && tn == nB(beg)) ||
		                    (last_join && tn == nB(ende));
		    if( to < 0)
			to = mk_nd( pt_nd(tn), SiO2);

		    if (join_existing) {
                        tsp=to; /* existing inner Si/oxide node stays fixed */
                        if(cap)to=tn==nB(beg)?first_cap:last_cap;
                    }
		    else {
		        (void)local_normal( tn, -2, Normal);
		        Normal[X] *= thick;
		        Normal[Y] *= thick;
		        tsp = SplitNode( to, Normal, mat_reg(ir));
		    }

                    if(cap && first_join && run_vertex==1)native_cap_fill(&first_chain,to,r);
		    /*make the connection to the prior node*/
		    if ((cap && last_join && run_vertex==run_edges) || (!join_existing && (closed_run || (tn != nB(beg)) || native_covered(nB(bp),mat_reg(ir))))) {
			if ( (bsp = node_mat( nB(bp), SiO2 )) < 0 )
			    bsp = mk_nd( pt_nd(nB(bp)), SiO2);
                        if(cap && last_join && run_vertex==run_edges)native_cap_fill(&last_chain,bsp,r);
                        else {
			v[0] = to; v[1] = tsp; v[2] = bsp;
			b[1] = b[2] = b[0] = BC_OFFSET;
			(void)mk_ele( 3, v, 3, b, r, FALSE);
                        }
		    }

		    if (!join_existing && (closed_run || (tn!=nB(ende)) || native_covered(nF(bp->next),mat_reg(ir)))) {
                        /* Ordinary next vertices move their original oxide
                         * node outward in SplitNode. Joined endpoints keep
                         * that inner node, so use their pre-noded outer cap. */
                        if(cap && last_join && nF(bp->next)==nB(ende))bsp=last_cap;
                        else if((bsp=node_mat(nF(bp->next),SiO2))<0)
                            bsp=mk_nd(pt_nd(nF(bp->next)),SiO2);
			v[0] = to; v[1] = tsp; v[2] = bsp;
			if (is_surf(pt_nd(bsp))) b[1] = EXPOSED;
			(void)mk_ele( 3, v, 3, b, r, FALSE);
		    }
		}
                free(first_chain.node);free(last_chain.node);
                first_chain.node=last_chain.node=NULL;first_chain.count=last_chain.count=0;
		g = 0; beg = bnd->prev;
	    }
	    if ( ende != NULL ) {
		if(cap && mode==TWOD)bd_connect_material("native oxide join union");
            else bd_connect("native oxide deposition");
		ir = -1;
	    }
	}
	break;
    }
    if(cap && mode==TWOD) {
        int bad;
        bd_connect_material("native oxide join union");
        if(!mesh_trial_validate_embedding(0,&bad))panic("invalid native oxide cap embedding");
    } else bd_connect("native oxide deposition");
}



/*-----------------NEW_LAYER--------------------------------------------
 * Inserts a thin new layer of oxide at the Si/SiO2 interface.
 *
 * Return value: whether geometry changed
 * Diagnostics:  abort on error
 *
 *----------------------------------------------------------------------*/
new_layer(Msil, Mox, oxide_grid, dt )
    int Msil, Mox;		/* Layer below, above interface */
    double oxide_grid;
    double dt;
{
    int io;			/* Oxide node being split */
    int i;
    int p;
    int e, ei, noth;
    int eperp;
    double fig;
    double magv;
    double evec[MAXDIM];
    float c[MAXDIM];
    double dis;
    int did_something = FALSE;

    /*
     * Geometry pass:
     * Loop through nodes on the oxide/silicon interface and see if their
     * local spacing is growing too large.
     */
    for( io = 0; io < nn; io++) {
	/*make sure we are at the oxide / silicon interface*/
	if( nd[ io]->mater != Mox) continue;
	if( node_mat( io, Msil) < 0) continue;
	p = pt_nd(io);
	magv = sqrt( dot(vel_arr(p), vel_arr(p)) );
	if ( magv == 0.0 ) continue;

	/*find the edge that is the most perpindicular to the velocity*/
	fig = 0.0;
	eperp = -1;
	for(i = 0; i < num_edge_nd(io); i++) {
	    e = edge_nd(io,i);
	    if ( nd_edg(e,0) == io )
		ei = 1;
	    else
		ei = 0;
	    noth = nd_edg(e,ei);
	    evec[0] = xcord(p) - xcord(pt_nd(noth));
	    evec[1] = ycord(p) - ycord(pt_nd(noth));

	    dis = dot( evec, vel_arr(p) ) / (magv * len_edg(e));

	    if ( dis > fig ) {
		fig = dis;
		eperp = e;
	    }
	}

	/* All adjacent oxide edges can point along the velocity at a mask
	 * corner. Then no edge opposes the motion, fig remains zero, and
	 * there is no layer-spacing candidate to split. */
	if ( (eperp >= 0) && (fig * len_edg(eperp) >= oxide_grid) &&
	     (dt * magv > 20.0e-8) ) {
	    /*split this edge*/
	    if ( nd_edg(eperp,0) == io )
		ei = 1;
	    else
		ei = 0;
	    noth = nd_edg(eperp,ei);
	    evec[0] = xcord(pt_nd(noth)) - xcord(p);
	    evec[1] = ycord(pt_nd(noth)) - ycord(p);
	    c[0] = 5.0e-8 * evec[0] / len_edg(eperp) + xcord(p);
	    c[1] = 5.0e-8 * evec[1] / len_edg(eperp) + ycord(p);
	    sp_edge(eperp, c, &p, FALSE);
	    did_something=TRUE;
	}
    }
    return(did_something);
}
	



/*-----------------SplitNode-------------------------------------------
 * This routine splits a boundary point in two, separated by a given delta.
 * The area between the two points is to be later filled with triangles
 * of material Mox, so a new node Mox is also created.
 * All the old nodes, bar the Msil one, are attached to the new point.
 * The new Mox one is attached to the old point.
 * This is so that the triangles of Mox do not have to be reconnected.
 * The Msil nodes and triangles are unaffected by the proceedings.
 * A side effect is that oxide can sneak in between nitride and silicon.
 * The return value is the new node created.
 *----------------------------------------------------------------------*/
SplitNode( aNode, v, mSil)
    int aNode;			/* Point being split as seen by mOx */
    float *v;			/* Displacement for new point */
    int mSil;			/* The other side of the interface */
{
    int i, OldPt, NewPt, *in, *inn, NewNd;
    float loc[MAXDIM];

    OldPt = nd[ aNode]->pt;
    loc[X] = pt[ OldPt]->cord[X] + v[X];
    loc[Y] = pt[ OldPt]->cord[Y] + v[Y];
    loc[Z] = pt[ OldPt]->cord[Z] + v[Z];

	/*
	 * Make a new point at the specified location
	 */
	NewPt = mk_pt( 3, loc );

	/*
	 * Transfer flags
	 */
	pt[ NewPt]->flags = pt[ OldPt]->flags;
	if( pt[ OldPt]->flags & SURFACE)
	    pt[ OldPt]->flags &= ~SURFACE;

	/*
	 * Attach old nodes to new point, except for the Msil one
	 * which moves down to the zeroth position at the old point.
	 */ 
	for (in = pt[ OldPt]->nd, inn = in + pt[ OldPt]->nn; in < inn; in++)
	    if (nd[ *in]->mater == mSil) pt[ OldPt]->nd[0] = *in;
	    else {
		pt[ NewPt]->nd[ pt[ NewPt]->nn++] = *in;
		nd[ *in]->pt = NewPt;
	    }
	pt[ OldPt]->nn = 1;

	/*
	 * Create the new oxide node and copy solution values
	 */
	NewNd = mk_nd( OldPt, nd[ aNode]->mater);
    
	for (i = 0; i < MAXIMP; i++)
	    nd[ NewNd]->sol[ i] = nd[ aNode ]->sol[i];

    return( NewNd);
}

/*-----------------EXPOSED----------------------------------------------
 * Tell whether a node is on an exposed edge.
 * Algorithm: look at edges by running through n->tri
 * You might think this routine unnecessary, but remember that backside
 * nodes have gas nodes and pt->flags& SURFACE set as well.
 *----------------------------------------------------------------------*/
exposed(in)
    int in;			/* node number */
{
    int ix, ixx, jin, j, nbr;
    
    if( node_mat( in, GAS) < 0) return(0);
    
    for( ixx = 0; ixx < num_tri_nd(in); ixx++) {
	ix = tri_nd(in, ixx);
	jin = nwhich( ix, in);
	for (j = jin; (j = (j+1)%3) != jin; ) {
	    nbr = tri[ ix]->nb[ j];
	    if (nbr == EXPOSED)
		return(1);
	}
    }
    return( 0);
}
	


/*-----------------SPLIT_THICK------------------------------------------
 * A heuristic function to decide when oxide nodes need splitting
 * Returns suitable displacement of node, 0 if none.
 * Diagnostics: return 0 anyway - can't hurt to not add stuff
 *----------------------------------------------------------------------*/
extern float oxide_grid;

split_thick( n, Mox, Msil, vdisp)
    int n;
    int Mox, Msil;
    float vdisp[MAXDIM];
{
    float SilToOx[MAXDIM], OxToSil[MAXDIM], dox, rough_perp();

    /* Compute a local normal */
    if (local_normal( n, Msil, OxToSil)) return(0);
    SilToOx[ X] = -OxToSil[ X];
    SilToOx[ Y] = -OxToSil[ Y];
    vdisp[ X] = 5e-8 * SilToOx[ X];
    vdisp[ Y] = 5e-8 * SilToOx[ Y];

    /* Estimate the oxde grid spacing on both sides */
    dox  = rough_perp( nd[n]->pt, Mox, SilToOx);	if (dox < 0) return(0);

    /* Provide a little hysteresis by adding grid anywhere between 0.25 and 1*/
    /* Danger: don't want slightly less thick oxides triggering just after field oxide */
    if (dox <= 0.25*oxide_grid)
	return(-1);
    else if (dox <= oxide_grid)
	return( 0);
    else
	return (1);
}

/*-----------------ROUGH_PERP-------------------------------------------
 * Get a rough idea of the spacing in some direction at a node.
 * Algorithm: Try intersecting the direction with the node's poly.
 * 	      Use closest node to direction if there is no intersection.
 * Return value: spacing value, scaled by length of direction vector.
 * Diagnostics:  return -1 for bad input or directions way out of the material.
 *----------------------------------------------------------------------*/
float rough_perp( p, matl, dir)
    int p;		/* Point */
    int matl;		/* + material = node */
    float dir[MAXDIM];	/* direction */
{
    int n, it, t, j, k, ot, o1, o2;
    float p_ot[MAXDIM], o1_o2[MAXDIM], dlen, dlen2, olen2, cos2, cosd, xn[MAXDIM];
    float dinter = 0, proj;

    /* Find the matl node at p */
    n = node_mat( pt[ p]->nd[0], matl);
    if (n == -1) return(-1);

    dlen = sqrt( dlen2 = dot( dir, dir));
    if (dlen == 0) return(-1);

    /* For each triangle containing that node */
    for (it = 0; it < num_tri_nd( n); it++) {
	t = tri_nd( n, it);

	/* Get the position of this node */
	j = nwhich( t, n);

	switch( mode ) {
	case ONED :
	    o1 = nd[ tri[ t]->nd[ (j+1)%nvrt]]->pt;
	    o1_o2[ X] = pt[ o1]->cord[ X] - pt[ p]->cord[X];

	    if ( o1_o2[X] * dir[X] > 0.0 ) {
		dinter = fabs(o1_o2[X]);
	    }
	    break;

	case TWOD :

	    o1 = nd[ tri[ t]->nd[ (j+1)%3]]->pt;
	    o2 = nd[ tri[ t]->nd[ (j+2)%3]]->pt;

	    /* The vector along that side */
	    o1_o2[ X] = pt[ o2]->cord[ X] - pt[ o1]->cord[X];
	    o1_o2[ Y] = pt[ o2]->cord[ Y] - pt[ o1]->cord[Y];

	    /* No point in looking for an intersection if the
	       side is parallel to the direction vector. */
	    olen2 = dot( o1_o2, o1_o2);
	    cosd  = dot( dir,   o1_o2);
	    cos2  = cosd*cosd/(dlen2*olen2);
	    if (cos2 > 0.9999) continue;

	    /* Compute the intersection - shouldn't fail.*/
	    assert(!lil (pt[ p]->cord, dir, pt[ o1]->cord, o1_o2, xn));

	    /* If the intersection is in the segment */
	    if (xn[1] > -0.01 && xn[1] < 1.01)
		/* and it's in the positive direction */
		if (xn[0] > 0)
		    /* save it */
		    dinter = xn[0];
	    break;
	}
    }

    /* If we got an intersection, we're outta of here */
    if (dinter > 0)
	return( dinter);

    /* No intersection. The point must be on the boundary and the
       direction vector points out of the material. What does the
       local spacing in that direction mean? Not much unless the
       direction vector is fairly close to the interface direction.
       Then we would say the distance to the point along the boundary
       closest to the direction is some measure of the spacing.
       So try the maximum projection of each nbr point on the direction*/

    /* Each neighbor point */
    for (it = 0; it < num_tri_nd( n); it++) {
	t = tri_nd( n, it);
	j = nwhich( t, n);
	for (k = 1; k <= 2; k++) {
    	    ot = nd[ tri[ t]->nd[ (j+k)%3]]->pt;
	    
	    /* The vector from p to other */
	    p_ot[ X] = pt[ ot]->cord[ X] - pt[ p]->cord[X];
	    p_ot[ Y] = pt[ ot]->cord[ Y] - pt[ p]->cord[Y];

	    /* Project it on the direction if it's reasonably close */
	    olen2 = dot( p_ot, p_ot);
	    cosd = dot( p_ot, dir);
	    if (cosd > 0) {
		cos2 = cosd*cosd / (dlen2*olen2);
		if (cos2 > 0.75) {
		    proj = cosd / dlen2;
		    if (proj > dinter) dinter = proj;
		}
	    }
	}
    }

    /* Hopefully we got a good projection (normalized above) */
    if (dinter > 0) return( dinter);

    /* I quit */
    else {
	fprintf( stderr, "Warning: unknown spacing for point %d[%d] in direction (%e, %e)\n",
		p, matl, dir[0], dir[1]);
	return(-1);
    }
}



/*-----------------LOCAL_NORMAL-----------------------------------------
 * Get some sort of idea of the interface normal at a node.
 * Given matl > 0, uses only edges facing that material
 * Given matl == GAS (or -1 for backward compatibility), only exposed edges
 * Given matl == -2 any edge that isn't reflecting
 * Algorithm: bisects the angle, by averaging the normalized
 * surface vectors.  Not prestissimo but it works.
 * Returns 0 for success, 1 for no normal (probably bad interface)
 * Also has an option to boost corner normals to help grid smoothness.
 *----------------------------------------------------------------------*/
local_normal (n, matl, ln)
    int n;		/* Node we care about */
    int matl;		/* Facing material, can be GAS */
    float ln[MAXDIM];	/* Returned vector */
{
    int it, t, j, nbr, head, tail, is_refl=0, p1, p2, v, nvec=0;
    float arrow[MAXDIM], tmp, lln, len, refl[MAXDIM], cosf, save[MAXDIM];

    switch( mode ) {
    case ONED :
	ln[0] = 1.0;
	ln[1] = 0.0;
	return(0);
	break;

    case TWOD:
	ln[ 0] = ln[ 1] = 0;
	cosf = 1.0;

	/* First figure if the node is on a reflecting boundary */
	/* Look at the edges that involve this node */
	for( it = 0; it < num_tri_nd( n); it++) {
	    t = tri_nd( n, it);
	    v = nwhich( t, n);
	    for( j = v; (j = (j+1)%3) != v; )
		if( tri[ t]->nb[ j] == NEUTRAL) {
		    /*Store the perpendicular to the boundary */
		    /*If there is more than one, well they should be parallel */
		    is_refl++;
		    p1 = nd[ n]->pt;
		    p2 = nd[ vertex( t, 3-j-v)]->pt;
		    refl[X] = pt[ p1]->cord[Y] - pt[ p2]->cord[Y];
		    refl[Y] = pt[ p2]->cord[X] - pt[ p1]->cord[X];
	    }
	}
	    
	/* Go thru triangle sides at the node */
	for( it = 0; it < num_tri_nd( n); it++) {
	    t = tri_nd( n, it);
	    v = nwhich( t, n);
	    for( j = v; (j = (j+1)%3) != v; ) {

		/* Only want interface sides */
		nbr = tri[ t]->nb[ j];
		if (matl > GAS){
		    if (nbr < 0 || mat_reg(reg_tri(nbr)) != matl) continue;}
		else if (matl == GAS || matl == -1) {
		    /*Not just any negative neighbor, would get 45o at corners*/
		    if (nbr != EXPOSED) continue;
		}
		else if (matl == -2) {
		    if (nbr == NEUTRAL ||
			(nbr >= 0 && mat_reg(reg_tri(nbr)) == nd[n]->mater))
			continue;
		}


		/* Get the vector along the interface */
		head = nd[ tri[ t]->nd[ (j+1)%3]]->pt;
		tail = nd[ tri[ t]->nd[ (j+2)%3]]->pt;
		arrow[ 0] = pt[ head]->cord[ 0] - pt[ tail]->cord[ 0];
		arrow[ 1] = pt[ head]->cord[ 1] - pt[ tail]->cord[ 1];

		/* Rotate it 90 degrees, INTO material matl */
		tmp = arrow[ 0];
		arrow[ 0] = -arrow[ 1];
		arrow[ 1] =  tmp;

		/* Normalize and add it into ln */
		len = hypot( arrow[0], arrow[1]);
		arrow[0] /= len; arrow[1] /= len;
		ln[ 0] += arrow[ 0];
		ln[ 1] += arrow[ 1];

		/* Also keep a copy for the "cos factor" */
		if( norm_style != 1) continue;
		++nvec;
		if( nvec == 1) {
		    save[0] = arrow[0];
		    save[1] = arrow[1];
		} else if (nvec == 2) {
		    cosf = fabs( dot( save, arrow)); /* just normalized 'em */
		    assert( cosf <= 1.001);
		} else {
		    fprintf( stderr, "local_normal: more than two edges at a node?\n");
		}
	    }
	}

	/* Project onto the reflecting boundary's tangent directly.  Subtracting
	 * the normal projection from ln loses the zero normal component through
	 * cancellation, and subsequent SplitNode calls turn that residual into
	 * spurious geometry outside an exactly axis-aligned boundary. */
	if( is_refl) {
	    double tx = -(double)refl[Y], ty = (double)refl[X];
	    double scale;
	    cosf = 1;
	    scale = tx*tx + ty*ty;
	    if (scale == 0.0) return(1);
	    scale = ((double)ln[X]*tx + (double)ln[Y]*ty) / scale;
	    ln[X] = scale*tx;
	    ln[Y] = scale*ty;
	}

	/* Normalize */
	lln = (float) hypot( (double) ln[ 0], (double) ln[ 1]);
	if (lln != 0) {
	    ln[ 0] /= lln;  ln[ 1] /= lln;
	    if( norm_style == 1) {
		cosf = sqrt( 0.5* (1 + cosf ));
		ln[ 0] /= cosf; ln[ 1] /= cosf;
	    }
	    return( 0);
	}
	else
	    return(1);
    }
    return(1);
}

/* Enough already */
dlocal_normal (n, matl, dln)
    int n; int matl; double dln[2];
{
    int oops;
    float ln[MAXDIM];
    if( !(oops = local_normal( n, matl, ln))) {
	dln[0] = ln[0];
	dln[1] = ln[1];
    }
    return( oops);
}
	




test_ln(i,m)
    int i, m;
{
    float ln[MAXDIM];
    (void)local_normal(i,m,ln);
    printf("local normal at node %d is [%f, %f]\n", i, ln[0], ln[1]);
}
