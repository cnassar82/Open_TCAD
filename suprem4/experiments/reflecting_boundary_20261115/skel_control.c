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
/*   skel.c                Version 5.1     */
/*   Last Modification : 7/3/91 15:40:51 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <global.h>
#include <constant.h>
#include <geom.h>
#include <regrid.h>
#include <material.h>
#include <impurity.h>
#include "refine.h"
#include "skel.h"

static int mesh_diag_enabled()
{
    char *v = (char *)getenv("SUPREM4_MESH_DIAG");
    return(v && *v && *v != '0');
}

#define endpt(A)	((dir)?(nF(A)):(nB(A)))
#define begpt(A)	((dir)?(nB(A)):(nF(A)))
#define nexedg(A)	((dir)?(A->next):(A->prev))

/*
 *	skel_reg - This routine creates a skeleton region from the
 *  original filled in region.
 */
skel_reg( r )
int r;
{
    register int i, ind, j;
    register int e1, el, e, sr;
    int t, f, n, done, fnd;
    double l_edge();

    /*create the skel*/
    sr = cr_sreg(mat_reg(r));

    if ( mode == ONED ) {
	bnd_reg(r,0,t,f);

	/*get the node on this face*/
	n = nd_face_ele(t,f,0);
	e1 = edge_nd(n,0);
	if ( nd_edg(e1,0) == n ) ind = 1;
	else ind = 0;
    }
    else {
	/*do a linear search to get a start edge*/
	for(fnd = FALSE, i = 0; (i < num_edg_reg(r)) && !fnd; i++) {
	    e1 = edg_reg(r, i);
	    if ( ask(edg[e1], ESURF) || ask(edg[e1], EBACK) )
		fnd = TRUE;
	    else if ( num_tri_edg(e1) == 1 )
		fnd = TRUE;
	    else {
		for(j = 0; j < num_tri_edg(e1); j++)
		    if ( reg_tri(tri_edg(e1,j)) != r ) fnd = TRUE;
	    }
	}
	ind = 1;
    }

    ad_edge( sr, e1, sreg[sr]->bnd, MAYBE, BEFORE );

    /*find the next connected boundary edge...*/
    el = e1;
    done = FALSE;

    while ( !done ) {
	n = nd_edg(el,ind);

	fnd = FALSE;
	for(i = 0; (i < num_edge_nd(n)) && !fnd; i++) {
	    e = edge_nd(n,i);
	    if ( (e != el) && (e != e1) && is_border(e) ) {
		/*make sure it is a border on this region!*/
		for(j = 0; j < num_ele_edg(e) &&
		       reg_tri(tri_edg(e,j)) != r; j++);

		if ( j < num_ele_edg(e) ) {
		    fnd = TRUE;
		    ad_edge(sr, e, sreg[sr]->bnd, MAYBE, BEFORE);
		}
	    }
	}
	if ( nd_edg(e,0) == n ) ind = 1;
	else ind = 0;
	el = e;
	done = !fnd;
    }
    hinit(sr, MAXFLOAT);
    ck_clock(sr, TRUE);
    return(sr);
}


#define IN 1
#define ON 2
#define OUT 0

static struct LLedge *classified_crossing(int, int, float *, int);

struct clip_crossing {
    int subject, clip;
    double ts, tc;
    float c[3];
};

static int same_crossing(float *a, float *b)
{
    double tolerance;
    int k;
    for(k=0;k<2;k++) {
        tolerance=2*FLT_EPSILON*fmax(fabs((double)a[k]),fabs((double)b[k]));
        if(fabs((double)a[k]-b[k])>tolerance)return FALSE;
    }
    return TRUE;
}

static int crossing_subject(const void *va, const void *vb)
{
    const struct clip_crossing *a=va,*b=vb;
    if(a->subject!=b->subject)return a->subject<b->subject ? -1:1;
    return a->ts<b->ts ? -1:a->ts>b->ts;
}

static int crossing_clip(const void *va, const void *vb)
{
    const struct clip_crossing *a=va,*b=vb;
    if(a->clip!=b->clip)return a->clip<b->clip ? -1:1;
    return a->tc<b->tc ? -1:a->tc>b->tc;
}

static int finite_crossing_possible(int a,int b,int c,int d)
{
    int s1=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c));
    int s2=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(d));
    int s3=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(a));
    int s4=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(b));
    return s1*s2<=0 && s3*s4<=0;
}

/* The arc walker needs at most one crossing per source edge. Resolve
 * multiple crossings before classification, rather than silently keeping
 * an arc that passes through the clipping polygon a second time. Split
 * both perimeters at the identical stored coordinate. */
static void node_multiple_crossings(int subject, int clip)
{
    struct clip_crossing *hits=NULL;
    struct LLedge *s,*q,*bnd;
    int count=0,capacity=0,i,j,k,n,original,edge,next,point,selected=0;
    int *pieces;
    double px,py,dx,dy,qx,qy,ex,ey,det,ts,tc;
    float c[3];
    s=sreg[subject]->bnd;
    do {
        px=cordinate(pt_nd(nB(s)),0);py=cordinate(pt_nd(nB(s)),1);
        dx=(double)cordinate(pt_nd(nF(s)),0)-px;
        dy=(double)cordinate(pt_nd(nF(s)),1)-py;
        q=sreg[clip]->bnd;
        do {
            qx=cordinate(pt_nd(nB(q)),0);qy=cordinate(pt_nd(nB(q)),1);
            ex=(double)cordinate(pt_nd(nF(q)),0)-qx;
            ey=(double)cordinate(pt_nd(nF(q)),1)-qy;
            det=dx*ey-dy*ex;
            if(det!=0 && finite_crossing_possible(pt_nd(nB(s)),pt_nd(nF(s)),
                                                  pt_nd(nB(q)),pt_nd(nF(q)))) {
                ts=((qx-px)*ey-(qy-py)*ex)/det;
                tc=((qx-px)*dy-(qy-py)*dx)/det;
                if(ts>0 && ts<1 && tc>=0 && tc<=1) {
                    if(count==capacity) {
                        capacity=capacity ? 2*capacity:16;
                        hits=sralloc(struct clip_crossing,capacity,hits);
                    }
                    hits[count].subject=s->edge;hits[count].clip=q->edge;
                    /* Sorting follows the database edge, not the perimeter. */
                    hits[count].ts=s->iscc ? ts:1-ts;
                    hits[count].tc=q->iscc ? tc:1-tc;
                    hits[count].c[0]=px+ts*dx;
                    hits[count].c[1]=py+ts*dy;hits[count].c[2]=0;
                    if(ex==0)hits[count].c[0]=qx;
                    if(ey==0)hits[count].c[1]=qy;
                    count++;
                }
            }
            q=q->next;
        }while(q!=sreg[clip]->bnd);
        s=s->next;
    }while(s!=sreg[subject]->bnd);
    if(!count){free(hits);return;}
    qsort(hits,count,sizeof(*hits),crossing_subject);
    /* Ignore duplicate reports of the same clipping corner. */
    for(i=0;i<count;i=j) {
        for(j=i;j<count && hits[j].subject==hits[i].subject;j++);
        n=0;
        for(k=i;k<j;k++)
            if(!n || !same_crossing(hits[k].c,hits[i+n-1].c))
                hits[i+n++]=hits[k];
        if(n>1)for(k=i;k<i+n;k++)hits[selected++]=hits[k];
    }
    count=selected;
    if(!count){free(hits);return;}
    for(k=0;k<2;k++) {
        if(k)qsort(hits,count,sizeof(*hits),crossing_clip);
        original=-1;edge=-1;
        for(i=0;i<count;i++) {
            n=k ? hits[i].clip:hits[i].subject;
            if(n!=original){original=n;edge=n;}
            for(j=0;j<2;j++)
                if(cordinate(pt_edg(edge,j),0)==hits[i].c[0] &&
                   cordinate(pt_edg(edge,j),1)==hits[i].c[1])break;
            if(j<2)continue;
            next=sp_gredg(edge,hits[i].c,&point);
            edge=next;
        }
    }
    free(hits);
    /* ON-to-ON pieces need an interior classification sample so the arc
     * walker can distinguish an outside interval from an inside interval. */
    n=sreg[subject]->len;pieces=salloc(int,n);
    bnd=s=sreg[subject]->bnd;i=0;
    do{pieces[i++]=s->edge;s=s->next;}while(s!=bnd);
    for(i=0;i<n;i++) {
        edge=pieces[i];
        if(!on_bound(pt_edg(edge,0),clip) ||
           !on_bound(pt_edg(edge,1),clip))continue;
        for(j=0;j<2;j++)c[j]=.5*((double)cordinate(pt_edg(edge,0),j)+
                                       cordinate(pt_edg(edge,1),j));
        c[2]=0;
        for(j=0;j<2;j++)
            if(cordinate(pt_edg(edge,j),0)==c[0] &&
               cordinate(pt_edg(edge,j),1)==c[1])break;
        if(j==2)sp_gredg(edge,c,&point);
    }
    free(pieces);
}

/*
 *	sub_skel - Subtract skeleton number 2 from number 1.  This
 *  routine can possible make new resulting regions.  The first passed
 *  region is returned modified, a new region number can also be returned.
 */
sub_skel(s1, s2, s)
int s1, s2;
int *s;
{
    int nr, *inside;
    int allin, allout;

    if(mode==TWOD)node_multiple_crossings(s1,s2);

    /*get the inside array*/
    inside = salloc( int, nn );

    /*check for the points in / out */
    check_in( s1, s2, inside, &allin, &allout );
    if (mesh_diag_enabled())
	fprintf(stderr, "[meshdiag] sub_skel classify s1=%d s2=%d s1_allin=%d s1_allout=%d len1=%d len2=%d\n",
		s1, s2, allin, allout, sreg[s1]->len, sreg[s2]->len);

    /*if all the points are in or out, return*/
    if ( allin ) { free(inside); return( -1 ); }
    if ( allout ) { free(inside); return( 0 ); }

    /*repeat the mission for the other region*/
    check_in( s2, s1, inside, &allin, &allout );
    if (mesh_diag_enabled())
	fprintf(stderr, "[meshdiag] sub_skel reverse_classify s1=%d s2=%d s2_allin=%d s2_allout=%d\n",
		s1, s2, allin, allout);

    if ( mode == ONED ) {
	nr = sub_1dskel(s1, s2, s, inside );
    }
    else {
	nr = sub_2dskel(s1, s2, s, &inside );
    }

    free(inside);
    if (mesh_diag_enabled())
	fprintf(stderr, "[meshdiag] sub_skel result s1=%d s2=%d outputs=%d\n", s1, s2, nr);
    return(nr);
}


sub_1dskel( s1, s2, s, inside )
int s1, s2;
int *s;
int *inside;
{
    struct LLedge *bp, *bq, *bnd, *edg_crs();
    int num = -1, e, enew;
    int g, ip;
    float c[MAXDIM];
    int out_node;
    int dir;

    /*compared to the other one, this is easy!*/
    /*walk the original region and find edges of the form in/on to out*/
    bnd = sreg[s1]->bnd;
    for(g = 1, bp = bnd; g || bp != bnd; g=0, bp = bp->next) {

	if ( ( (inside[nB(bp)] != OUT) && (inside[nF(bp)] == OUT) ) ||
	     ( (inside[nF(bp)] != OUT) && (inside[nB(bp)] == OUT) ) ) {

	    if ( (inside[nB(bp)] != OUT) && (inside[nF(bp)] == OUT) )
		dir = 1;
	    else
		dir = 0;

	    /*create a new skeleton...*/
	    num++;
	    s[num] = cr_sreg(sreg[s1]->mat);

	    /*add the fragment of the initial edge*/
	    e = bp->edge;
	    bq = edg_crs( e, s2, c );
	    out_node = endpt(bp);

	    /*shorten the edge*/
	    enew = sp_edge( e, c, &ip, TRUE );
	    mk_surf(ip);

	    /*convolutions to figure out the edge*/
	    if ( enew < 0 ) {
		if ( nd_edg(e,enew+2) != out_node )
		    ad_edge(s[num], e, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    }
	    else if ( (out_node == nd_edg(e,0)) || (out_node == nd_edg(e,1)) )
		ad_edge(s[num], e, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    else
		ad_edge(s[num], enew, sreg[s[num]]->bnd, MAYBE, BEFORE);

	    /*walk around the existing region until we cross again*/
	    for(bq = nexedg(bp);
		(inside[begpt(bq)] == OUT) && (inside[endpt(bq)] != IN);
		bq = nexedg(bq)) {
		e = bq->edge;
		ad_edge(s[num], e, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    }
	}
    }
    return(num+1);
}


#define min2(A,B) ((A)<(B))?(A):(B)

/* Position along the positively oriented clipping perimeter.  Keep the
 * edge identity as well as the position: several material fragments may
 * cross the same clipping edge. */
static double clip_position(sr, edge, node)
int sr, node;
struct LLedge *edge;
{
    struct LLedge *bp = sreg[sr]->bnd;
    double pos = 0, dx, dy, x, y, len, alpha;
    do {
        len = l_edge(bp->edge);
        if (bp == edge) {
            x = cordinate(pt_nd(nB(bp)),0);
            y = cordinate(pt_nd(nB(bp)),1);
            dx = (double)cordinate(pt_nd(nF(bp)),0) - x;
            dy = (double)cordinate(pt_nd(nF(bp)),1) - y;
            if (len == 0) panic("zero length clipping edge");
            alpha = (((double)cordinate(pt_nd(node),0)-x)*dx +
                     ((double)cordinate(pt_nd(node),1)-y)*dy)/(dx*dx+dy*dy);
            if (alpha < 0) alpha = 0;
            if (alpha > 1) alpha = 1;
            return pos + alpha*len;
        }
        pos += len;
        bp = bp->next;
    } while (bp != sreg[sr]->bnd);
    panic("missing clipping edge");
    return 0;
}

static int same_clip_point(int a, int b)
{
    double tolerance, x, y;
    int k;
    /* Intersections and original corners are stored as floats.  At an
     * axis-aligned boundary a computed intersection can differ from its
     * corner by roundoff (e.g. x=-1e-23 instead of zero).  Inserting both
     * creates a reversing edge and changes the polygon's winding number.
     * Compare at storage precision, without a fixed physical length cutoff. */
    for (k=0; k<2; k++) {
        x = fabs((double)cordinate(a,k));
        y = fabs((double)cordinate(b,k));
        tolerance = 2 * FLT_EPSILON * fmax(x,y);
        if (fabs((double)cordinate(a,k)-cordinate(b,k)) > tolerance)
            return FALSE;
    }
    return TRUE;
}

/* Splitting an interface edge can append nodes in several materials.
 * Keep their ON classifications in bounds as the source skeleton grows. */
static void clip_nodes(int **inside, int *capacity, int point)
{
    int j;
    if (nn > *capacity) {
        *inside = sralloc(int, nn, *inside);
        for (j=*capacity; j<nn; j++) (*inside)[j] = ON;
        *capacity = nn;
    }
    for (j=0; j<num_nd(point); j++) (*inside)[nd_pt(point,j)] = ON;
}

sub_2dskel( s1, s2, s, classifications )
int s1, s2;
int *s;
int **classifications;
{
    int i, e, e1, enew, eadd;
    int first_node, last_node, out_node;
    int num = -1;
    int g, ip, in1;
    int *inside = *classifications, capacity = nn;
    int count, base, a, b, k, component, current, next, used[MAXREG];
    int fragments[MAXREG], starts[MAXREG], ends[MAXREG], successor[MAXREG];
    int incoming[MAXREG];
    double startpos[MAXREG], endpos[MAXREG], perimeter, delta, best, walked;
    struct LLedge *endedge[MAXREG], *walk;
    float c[3];
    struct LLedge *bnd, *bp, *bq, *bts, *bte, *edg_crs();
    double hsp,l_edge();

    base = nsreg;
    /* Collect retained arcs first. An arc is not necessarily a component:
     * the inner and outer walls of a U-shaped liner must be joined along
     * the cut, rather than each being closed across the trench void. */
    bnd = sreg[s1]->bnd;
    for(g = 1, bp = bnd; g || bp != bnd; g=0, bp = bp->next) {
	if ( ask(edg[bp->edge], MARKED) ) continue;
	set(edg[bp->edge], MARKED);

	/* Exposure belongs to a finite new cut edge. Merely touching an
	 * existing interface must not add gas nodes beneath a covering layer. */
	/* Removed interior points can expose a different material underneath
	 * (for example oxide beneath a patterned resist opening). */
	if (inside[nB(bp)] == IN) mk_surf(pt_nd(nB(bp)));

	if ( (inside[nB(bp)] != OUT) && (inside[nF(bp)] == OUT) ) {

	    /*create a new skeleton...*/
	    num++;
	    if (num >= MAXREG) panic("too many retained boundary arcs");
	    s[num] = cr_sreg(sreg[s1]->mat);

	    /*add the fragment of the initial edge*/
	    e = bp->edge;
	    out_node = nF(bp);
	    hsp = min2( get_space(pt_edg(e,0)), get_space(pt_edg(e,0)));
	    bte = classified_crossing(e, s2, c,
                                     inside[nB(bp)] == ON ? nB(bp) : -1);

	    /*shorten the edge*/
	    enew = sp_edge_keep( e, c, &ip, TRUE, out_node );
	    clip_nodes(&inside, &capacity, ip);
	    *classifications = inside;
	    set_space( ip, min2(get_space(ip), hsp) );

	    /*save the first node*/
	    first_node = node_mat(nd_pt(ip,0), sreg[s1]->mat);

	    /*convolutions to figure out the edge*/
	    if ( enew < 0 ) {
		if ( nd_edg(e,enew+2) != out_node )
		    ad_edge(s[num], eadd=e, sreg[s[num]]->bnd, MAYBE, BEFORE);
		else
		    eadd = e;
	    }
	    else {
		if ( (out_node == nd_edg(e,0)) || (out_node == nd_edg(e,1)) )
		    eadd = e;
		else
		    eadd = enew;
		ad_edge(s[num], eadd, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    }
	    set(edg[eadd], REGS);
	    set(edg[eadd], MARKED);
	    while( bp->edge != eadd ) {
		set(edg[bp->edge], MARKED);
		bp = bp->next;
	    }

	    /*walk around the existing region until we cross again*/
	    for(bq = bp->next;
		(inside[nF(bq)] == OUT);
		bq = bq->next) {
		e = bq->edge;
		set(edg[e],MARKED);
		ad_edge(s[num], e, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    }
	    bp = bq;

	    /*bq is an edge that cross, add the possible frag*/
	    e = bq->edge;
	    out_node = nB(bq);
	    hsp = min2( get_space(pt_edg(e,0)), get_space(pt_edg(e,0)));
	    bts = classified_crossing(e, s2, c,
                                     inside[nF(bq)] == ON ? nF(bq) : -1);

	    /*shorten the edge*/
	    enew = sp_edge_keep( e, c, &ip, TRUE, out_node );
	    clip_nodes(&inside, &capacity, ip);
	    *classifications = inside;
	    set_space( ip, min2(get_space(ip), hsp) );

	    /*save the last node*/
	    last_node = node_mat(nd_pt(ip,0), sreg[s1]->mat);

	    /*convolutions to figure out the edge*/
	    if ( enew < 0 ) {
		if ( nd_edg(e,enew+2) != out_node )
		    ad_edge(s[num], eadd=e, sreg[s[num]]->bnd, MAYBE, BEFORE);
		else
		    eadd = e;
	    }
	    else {
		if ( (out_node == nd_edg(e,0)) || (out_node == nd_edg(e,1)) )
		    eadd = e;
		else
		    eadd = enew;
		ad_edge(s[num], eadd, sreg[s[num]]->bnd, MAYBE, BEFORE);
	    }
	    set(edg[eadd], MARKED);
	    set(edg[eadd], REGS);

	    fragments[num] = s[num];
	    starts[num] = first_node;
	    ends[num] = last_node;
	    endedge[num] = bts;
	    startpos[num] = clip_position(s2, bte, first_node);
	    endpos[num] = clip_position(s2, bts, last_node);

	}
    }

    count = num+1;
    perimeter = 0;
    walk = sreg[s2]->bnd;
    do { perimeter += l_edge(walk->edge); walk = walk->next; }
    while (walk != sreg[s2]->bnd);
    for (a=0; a<count; a++) used[a] = incoming[a] = 0;
    /* The retained arcs follow the subject boundary positively. The new
     * cut edges follow the clipping boundary negatively, so material stays
     * on the same side throughout the assembled component. */
    for (a=0; a<count; a++) {
        best = MAXFLOAT;
        successor[a] = -1;
        for (b=0; b<count; b++) {
            delta = endpos[a] - startpos[b];
            if (delta < 0) delta += perimeter;
            /* Two retained components can touch at a single ON vertex.
             * A zero-length cut between their arcs is not an interior
             * material connection; close each through its nonzero cut. */
            if (a != b && delta == 0) {
                int incoming_edge, outgoing_edge, before, after, junction;
                double ux,uy,vx,vy,cx,cy,source_turn,clip_turn;
                struct LLedge *cut;
                if(ends[a]!=starts[b])continue;
                junction=ends[a];
                incoming_edge=sreg[fragments[a]]->bnd->prev->edge;
                outgoing_edge=sreg[fragments[b]]->bnd->edge;
                before=nd_edg(incoming_edge,0)==junction ?
                       nd_edg(incoming_edge,1):nd_edg(incoming_edge,0);
                after=nd_edg(outgoing_edge,0)==junction ?
                      nd_edg(outgoing_edge,1):nd_edg(outgoing_edge,0);
                ux=(double)cordinate(pt_nd(junction),0)-cordinate(pt_nd(before),0);
                uy=(double)cordinate(pt_nd(junction),1)-cordinate(pt_nd(before),1);
                vx=(double)cordinate(pt_nd(after),0)-cordinate(pt_nd(junction),0);
                vy=(double)cordinate(pt_nd(after),1)-cordinate(pt_nd(junction),1);
                /* Trace the face with retained material on its left.
                 * At a shared vertex choose the first clockwise outgoing
                 * ray from the reversed incoming ray. This distinguishes
                 * a connected tangent from two components touching at a
                 * reflex vertex, including corners in the clip perimeter. */
                ux=-ux; uy=-uy;
                cut=endedge[a];
                if(pt_nd(nB(cut))==pt_nd(junction) ||
                   same_clip_point(pt_nd(nB(cut)),pt_nd(junction)))cut=cut->prev;
                cx=(double)cordinate(pt_nd(nB(cut)),0)-cordinate(pt_nd(junction),0);
                cy=(double)cordinate(pt_nd(nB(cut)),1)-cordinate(pt_nd(junction),1);
                source_turn=-atan2(ux*vy-uy*vx,ux*vx+uy*vy);
                clip_turn=-atan2(ux*cy-uy*cx,ux*cx+uy*cy);
                if(source_turn<0)source_turn+=2*PI;
                if(clip_turn<0)clip_turn+=2*PI;
                if(source_turn>clip_turn)continue;
            }
            if (delta < best) { best = delta; successor[a] = b; }
        }
        if (successor[a] < 0 || ++incoming[successor[a]] != 1) {
            if(mesh_diag_enabled()) {
                fprintf(stderr,"[meshdiag] ambiguous arc %d -> %d\n",a,successor[a]);
                for(b=0;b<count;b++)
                    fprintf(stderr,"[meshdiag] arc %d nodes %d -> %d positions %.17g -> %.17g\n",
                            b,starts[b],ends[b],startpos[b],endpos[b]);
            }
            panic("ambiguous clipped boundary connectivity");
        }
    }
    component = 0;
    for (a=0; a<count; a++) {
        if (used[a]) continue;
        s[component] = cr_sreg(sreg[s1]->mat);
        current = a;
        do {
            if (used[current]) panic("open clipped boundary");
            used[current] = 1;
            bp = sreg[fragments[current]]->bnd;
            do {
                ad_edge(s[component], bp->edge, sreg[s[component]]->bnd,
                        MAYBE, BEFORE);
                bp = bp->next;
            } while (bp != sreg[fragments[current]]->bnd);
            next = successor[current];
            last_node = ends[current];
            delta = endpos[current] - startpos[next];
            if (delta < 0) delta += perimeter;
            walk = endedge[current];
            walked = endpos[current] -
                     clip_position(s2, walk, nB(walk));
            /* Preserve every intervening clipping corner. Distance-based
             * corner omission changes the material boundary. */
            for (k=0; walked < delta && k<=sreg[s2]->len; k++) {
                ip = pt_nd(nB(walk));
                /* The intersection and a clipping corner may have
                 * different point IDs but identical stored coordinates. */
                if (!same_clip_point(ip,pt_nd(last_node)) &&
                    !same_clip_point(ip,pt_nd(starts[next]))) {
                    mk_surf(ip);
                    in1 = node_mat(nd_pt(ip,0),sreg[s1]->mat);
                    if (in1 == -1) in1 = mk_nd(ip,sreg[s1]->mat);
                    e1 = mk_edg(last_node,in1);
                    mk_surf(pt_nd(last_node));
                    set(edg[e1],REGS); set(edg[e1],ESURF);
                    ad_edge(s[component],e1,sreg[s[component]]->bnd,MAYBE,BEFORE);
                    last_node = in1;
                }
                walk = walk->prev;
                walked += l_edge(walk->edge);
            }
            if (pt_nd(last_node) != pt_nd(starts[next])) {
                mk_surf(pt_nd(last_node));
                mk_surf(pt_nd(starts[next]));
                e1 = mk_edg(last_node,starts[next]);
                set(edg[e1],REGS); set(edg[e1],ESURF); set(edg[e1],MARKED);
                ad_edge(s[component],e1,sreg[s[component]]->bnd,MAYBE,BEFORE);
            }
            current = next;
        } while (current != a);
        component++;
    }
    /* Temporary arc skeletons precede the assembled components. */
    for (i=count-1; i>=0; i--) free_skel(base+i);
    for (i=0; i<component; i++) s[i] -= count;
    for(i=0; i < component; i++) {
	ck_clock(s[i], TRUE);
	if (mesh_diag_enabled())
	    fprintf(stderr, "[meshdiag] sub_2dskel output index=%d sreg=%d len=%d\n",
		    i, s[i], sreg[s[i]]->len);
    }

    return(component);
}



check_in(s1, s2, inside, ai, ao )
int s1, s2;
int *inside;
int *ai, *ao;
{
    struct LLedge *bp, *bnd;
    int allin, allout, cnt;
    int g;

    /*first walk the region points and test in/out*/
    allin = allout = cnt = 0;
    bnd = sreg[s1]->bnd;

    for(g=1, bp = bnd; g || bp != bnd; g = 0,bp = bp->next) {
	if ( on_bound(pt_nd(nB(bp)), s2)) {
	    inside[nB(bp)] = ON;
	    allin++;
	    allout++;
	}
	else {
	    if (pt_in_skel(cord_arr(pt_nd(nB(bp))), s2)) {
		inside[nB(bp)] = IN;
		allin++;
	    }
	    else {
		inside[nB(bp)] = OUT;
		allout++;
	    }
	}
	cnt++;
	if ( mode == ONED ) {
	    if ( on_bound(pt_nd(nF(bp)), s2)) {
		inside[nF(bp)] = ON;
		allin++;
		allout++;
	    }
	    else {
		if (pt_in_skel(cord_arr(pt_nd(nF(bp))), s2)) {
		    inside[nF(bp)] = IN;
		    allin++;
		}
		else {
		    inside[nF(bp)] = OUT;
		    allout++;
		}
	    }
	    cnt++;
	}
    }

    *ai = allin == cnt;
    *ao = allout == cnt;
}



/* The same finite-segment predicate is used for classification and for
 * endpoint intersections.  A fixed one-Angstrom ON band classified nearby
 * parallel edges as crossing even when no intersection could be found. */
static int clip_on_segment(int p, int a, int b)
{
    double px,py,ax,ay,bx,by,dx,dy,len,t,xtol,ytol;
    px=cordinate(p,0); py=cordinate(p,1);
    ax=cordinate(a,0); ay=cordinate(a,1);
    bx=cordinate(b,0); by=cordinate(b,1);
    dx=bx-ax; dy=by-ay; len=hypot(dx,dy);
    /* Propagate coordinate precision through the dot/cross products.
     * A large lateral coordinate must not create an equally large depth
     * tolerance on an almost horizontal interface. */
    xtol=2*FLT_EPSILON*fmax(fabs(px),fmax(fabs(ax),fabs(bx)));
    ytol=2*FLT_EPSILON*fmax(fabs(py),fmax(fabs(ay),fabs(by)));
    if (len==0) return same_clip_point(p,a);
    t=(px-ax)*dx+(py-ay)*dy;
    return t>=-(xtol*fabs(dx)+ytol*fabs(dy)) &&
           t<=len*len+xtol*fabs(dx)+ytol*fabs(dy) &&
           fabs((px-ax)*dy-(py-ay)*dx)<=xtol*fabs(dy)+ytol*fabs(dx);
}

/* An ON endpoint is already the junction chosen by classification. Re-solving
 * its almost parallel intersection after a shared edge has been split can
 * produce a different rounded junction in the adjacent region. Preserve that
 * identity when it belongs to the clipping segment actually intersected. */
static struct LLedge *classified_crossing(int edge, int clip, float *c, int on_node)
{
    struct LLedge *bp, *edg_crs();
    int point;
    bp = edg_crs(edge,clip,c);
    if (on_node >= 0 && bp) {
        point = pt_nd(on_node);
        /* ON can also mean the endpoint lies on a different, collinear
         * clip segment. Do not preserve it when this edge leaves that segment
         * and intersects a different side of the mask farther along. */
        if (clip_on_segment(point,pt_nd(nB(bp)),pt_nd(nF(bp)))) {
            c[0] = cordinate(point,0);
            c[1] = cordinate(point,1);
            c[2] = 0;
        }
    }
    if(mesh_diag_enabled() && bp)fprintf(stderr,
        "[meshdiag] crossing xy=%.17g,%.17g clip=%.17g,%.17g:%.17g,%.17g\n",
        (double)c[0],(double)c[1],(double)cordinate(pt_nd(nB(bp)),0),
        (double)cordinate(pt_nd(nB(bp)),1),(double)cordinate(pt_nd(nF(bp)),0),
        (double)cordinate(pt_nd(nF(bp)),1));
    return bp;
}

int on_bound(p, sr)
int p;
int sr;			/*The coordinates of the points. */
{
    int f;
    float len, alph[MAXDIM];
    struct LLedge *bp, *bnd;
    void d_perp();

    /* If the perpindicular distance to the line is zero, its on the line*/
    bnd = sreg[sr]->bnd;
    for(bp = bnd, f=1; f || bp != bnd; f = 0, bp = bp->next) {
	if ( mode == ONED ) {
	    if ( pdist( p, pt_nd(nB(bp)) ) < 1.0e-8 ||
	         pdist( p, pt_nd(nF(bp)) ) < 1.0e-8 )
		return(TRUE);
	}
	else {
	    if (clip_on_segment(p,pt_nd(nB(bp)),pt_nd(nF(bp))))
		return TRUE;
	}
    }

    /* What a short, strange trip it's been */
    return (FALSE);
}




plot_on(sr,inside)
int sr;
int *inside;
{
    int f;
    struct LLedge *bp, *bnd;

    bnd = sreg[sr]->bnd;
    for(bp = bnd, f=1; f || bp != bnd; f = 0, bp = bp->next) {
	node_pl( nB(bp), inside[nB(bp)]+1 );
    }
}

/*-----------------PTNREG-----------------------------------------------
 * Determine whether the point (x0,y0) is within the region. Algorithm
 * 112 from CACM. This little gem works for clockwise, anticlockwise,
 * or even self-crossing regions.
 *----------------------------------------------------------------------*/
int pt_in_skel(c, sr)
float c[3];
int sr;			/*The coordinates of the points. */
{
    int flag, f;
    double xi, yi, xin, yin;
    struct LLedge *bp, *bnd;

  /*
   * If a line to infinity crosses an odd number of edges, its inside.
   */
    flag = TRUE;
    bnd = sreg[sr]->bnd;

    for(bp = bnd, f=1; f || bp != bnd; f = 0, bp = bp->next) {
	xi = cordinate( pt_edg(bp->edge,0), 0 );
	yi = cordinate( pt_edg(bp->edge,0), 1 );

	xin = cordinate( pt_edg(bp->edge,1), 0 );
	yin = cordinate( pt_edg(bp->edge,1), 1 );

	if ( mode == ONED ) {
	    if ( (c[0] <= xi) == (c[0] > xin) ) flag = !flag;
	}
	else {
	    if ((c[1] <= yi) == (c[1] > yin)) {
		if ((c[0]-xi - (c[1]-yi)*(xin-xi)/(yin-yi) ) < 0) {
		    flag = ! flag;
		}
	    }
	}
    }

    /* What a short, strange trip it's been */
    return (!flag);
}





/*
 * calculate the point at which edge e crosses region s2
 */
struct LLedge *edg_crs( e, s2, c )
int e;
int s2;
float c[MAXDIM];
{
    int e1, g;
    struct LLedge *bp, *bnd;
    float p[MAXDIM], dp[MAXDIM], q[MAXDIM], dq[MAXDIM];
    float alph[MAXDIM], lp, lq;

    /* Search all edges for a true intersection before accepting an ON
     * tolerance match. Float line parameters can select the edge just
     * before a nearby corner, then insert that corner behind the cut and
     * create a small self-crossing loop. Evaluate stored float coordinates
     * in double precision, without enlarging either segment. */
    if (mode == TWOD) {
        double px, py, dx, dy, qx, qy, ex, ey, det, a, b;
        px = cordinate(pt_edg(e,0),0);
        py = cordinate(pt_edg(e,0),1);
        dx = (double)cordinate(pt_edg(e,1),0)-px;
        dy = (double)cordinate(pt_edg(e,1),1)-py;
        bnd = sreg[s2]->bnd;
        bp = bnd;
        do {
            qx = cordinate(pt_edg(bp->edge,0),0);
            qy = cordinate(pt_edg(bp->edge,0),1);
            ex = (double)cordinate(pt_edg(bp->edge,1),0)-qx;
            ey = (double)cordinate(pt_edg(bp->edge,1),1)-qy;
            det = dx*ey-dy*ex;
            if (det != 0 && finite_crossing_possible(pt_edg(e,0),pt_edg(e,1),
                                    pt_edg(bp->edge,0),pt_edg(bp->edge,1))) {
                a = ((qx-px)*ey-(qy-py)*ex)/det;
                b = ((qx-px)*dy-(qy-py)*dx)/det;
                if (a>=0 && a<=1 && b>=0 && b<=1) {
                    c[0] = px+a*dx;
                    c[1] = py+a*dy;
                    if(ex==0)c[0]=qx;
                    if(ey==0)c[1]=qy;
                    /* Preserve an existing ON junction only when this
                     * intersection rounds to that endpoint. An ON endpoint
                     * elsewhere on a collinear overlap is not this crossing. */
                    {
                        int endpoint, point;
                        double xtol=2*FLT_EPSILON*fmax(fabs(px),fabs(px+dx));
                        double ytol=2*FLT_EPSILON*fmax(fabs(py),fabs(py+dy));
                        for(endpoint=0;endpoint<2;endpoint++) {
                            point=pt_edg(e,endpoint);
                            if(fabs((double)c[0]-cordinate(point,0))<=xtol &&
                               fabs((double)c[1]-cordinate(point,1))<=ytol &&
                               clip_on_segment(point,pt_nd(nB(bp)),pt_nd(nF(bp)))) {
                                c[0]=cordinate(point,0); c[1]=cordinate(point,1);
                                break;
                            }
                        }
                    }
                    return bp;
                }
            }
            bp = bp->next;
        } while (bp != bnd);
        /* ON endpoints include collinear overlaps and storage roundoff.
         * Return the existing endpoint, preserving the retained geometry. */
        bp=bnd;
        do {
            int endpoint;
            for (endpoint=0;endpoint<2;endpoint++) {
                int point=pt_edg(e,endpoint);
                if (clip_on_segment(point,pt_nd(nB(bp)),pt_nd(nF(bp)))) {
                    c[0]=cordinate(point,0); c[1]=cordinate(point,1);
                    return bp;
                }
            }
            bp=bp->next;
        } while (bp!=bnd);
        return NULL;
    }

    p[0] = cordinate( pt_edg(e,0), 0);
    p[1] = (mode==ONED)?0.0:cordinate( pt_edg(e,0), 1);
    dp[0] = cordinate( pt_edg(e,1), 0) - p[0];
    dp[1] = (mode==ONED)?0.0:cordinate( pt_edg(e,1), 1) - p[1];
    lp = sqrt(dp[0] * dp[0] + dp[1] * dp[1]);

    bnd = sreg[s2]->bnd;
    for(bp = bnd, g=1; g || bp != bnd; g = 0, bp=bp->next) {
	e1 = bp->edge;
	q[0] = cordinate( pt_edg(e1,0), 0);
	q[1] = (mode==ONED)?0.0:cordinate( pt_edg(e1,0), 1);
	dq[0] = cordinate( pt_edg(e1,1), 0) - q[0];
	dq[1] = (mode==ONED)?0.0:cordinate( pt_edg(e1,1), 1) - q[1];
	lq = sqrt(dq[0] * dq[0] + dq[1] * dq[1]);

	if ( mode == TWOD ) {
	    if ( lil(p,dp,q,dq,alph) == 0 ) {

		/*if they cross, nice and easy*/
		if ((alph[0] >= 0.0) && (alph[0] <= 1.0) &&
		    (alph[1] >= 0.0) && (alph[1] <= 1.0) ) {
		    c[0] = alph[0] * dp[0] + p[0];
		    c[1] = alph[0] * dp[1] + p[1];
		    return( bp );
		}
		/*we have the on condition*/
		else if ( (alph[0] >= 0.0) && (alph[0] <= 1.0) ) {
		    if ( ( fabs(alph[1] * lq) < 1.0e-8 ) ||
		         ( fabs(alph[1] * lq - lq) < 1.0e-8 ) ) {
			c[0] = alph[0] * dp[0] + p[0];
			c[1] = alph[0] * dp[1] + p[1];
			return( bp );
		    }
		}
		else if ( (alph[1] >= 0.0) && (alph[1] <= 1.0) ) {
		    if ( ( fabs(alph[0] * lp) < 1.0e-8 ) ||
		         ( fabs(alph[0] * lp - lp) < 1.0e-8 ) ) {
			c[0] = alph[0] * dp[0] + p[0];
			c[1] = alph[0] * dp[1] + p[1];
			return( bp );
		    }
		}
	    }
	}
	else {
	    if (lil(p,dp,q,dq,alph) == 0) {
		if ((alph[0] >= 0.0) && (alph[0] <= 1.0) &&
		    (alph[1] >= 0.0) && (alph[1] <= 1.0) ) {
		    c[0] = alph[0] * dp[0] + p[0];
		    return( bp );
		}
	    }
	}
    }
    return(NULL);
}
