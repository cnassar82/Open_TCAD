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
/*   kill_node.c                Version 5.1     */
/*   Last Modification : 7/3/91 08:20:57 */

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "constant.h"
#include "global.h"
#include "dbaccess.h"
#include "material.h"
#include "skel.h"

static int kn_diag_on()
{
    char *v = (char *)getenv("SUPREM4_MESH_DIAG");
    return(v && *v && *v != '0');
}

/*-----------------KILL_NODE--------------------------------------------
 * Remove a node from the grid and fix up the triangulation.
 * Algorithm:
 *     if edge node, internalize the node to be removed
 *     call kill_internal_node to do the dirty work
 * Node is passed in as a triangle/vertex pair to avoid O(N) search every
 * time we do this.
 *----------------------------------------------------------------------*/
/* Exact simplicity check for a proposed cavity, including its new closing
 * chord. A boundary fan may be valid while the polygon formed by deleting its
 * center is not. In that case retain the vertex instead of corrupting the mesh. */
static int removal_polygon_valid(int *nodes,int count)
{
    int i,j,k,a,b,c,d,s1,s2,s3,s4,nonflat=FALSE,extreme=0;
    for(i=0;i<count;i++) {
        a=pt_nd(nodes[i]);b=pt_nd(nodes[(i+1)%count]);c=pt_nd(nodes[(i+2)%count]);
        d=pt_nd(nodes[extreme]);
        if(xcord(a)<xcord(d) || (xcord(a)==xcord(d) && ycord(a)<ycord(d)))extreme=i;
        if(xcord(a)==xcord(b) && ycord(a)==ycord(b))return FALSE;
        s1=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c));
        if(s1)nonflat=TRUE;
        else {
            double ux=(double)xcord(b)-xcord(a),uy=(double)ycord(b)-ycord(a);
            double vx=(double)xcord(c)-xcord(b),vy=(double)ycord(c)-ycord(b);
            if(ux*vx+uy*vy<0)return FALSE;
        }
        for(j=i+1;j<count;j++) {
            k=(j+1)%count;
            if(j==(i+1)%count || k==i)continue;
            c=pt_nd(nodes[j]);d=pt_nd(nodes[k]);
            if(fmax(xcord(a),xcord(b))<fmin(xcord(c),xcord(d)) ||
               fmax(xcord(c),xcord(d))<fmin(xcord(a),xcord(b)) ||
               fmax(ycord(a),ycord(b))<fmin(ycord(c),ycord(d)) ||
               fmax(ycord(c),ycord(d))<fmin(ycord(a),ycord(b)))continue;
            s1=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c));
            s2=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(d));
            s3=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(a));
            s4=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(b));
            if(s1*s2<=0 && s3*s4<=0)return FALSE;
        }
    }
    /* A clockwise open-fan polygon lies outside the original star. Letting
     * ck_clock reverse it would fill the neighboring region and overlap
     * existing triangles. Test the convex lexicographic extreme exactly. */
    a=pt_nd(nodes[(extreme+count-1)%count]);b=pt_nd(nodes[extreme]);
    c=pt_nd(nodes[(extreme+1)%count]);
    return nonflat && orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c))>0;
}

#ifdef SUPREM4_TEST_REMOVAL
int test_removal_polygon_valid(int *nodes,int count)
{return removal_polygon_valid(nodes,count);}
#endif

/* The only new boundary edge of an open fan is its closing chord. It must
 * not cross an edge outside the removed point's star or add a third owner
 * to an already existing edge. Check physical points across materials. */
static int removal_chord_valid(int in,int a,int b)
{
    int t,f,j,p=pt_nd(in),c,d,shared,s1,s2,s3,s4,incident;
    a=pt_nd(a);b=pt_nd(b);
    for(t=0;t<ne;t++) {
        if(dead_tri(t) || !leaf(tri[t]))continue;
        incident=FALSE;
        for(j=0;j<3;j++)if(pt_nd(vert_tri(t,j))==p)incident=TRUE;
        if(incident)continue;
        for(f=0;f<3;f++) {
            c=pt_nd(nd_face_ele(t,0,f));d=pt_nd(nd_face_ele(t,1,f));
            if(fmax(xcord(a),xcord(b))<fmin(xcord(c),xcord(d)) ||
               fmax(xcord(c),xcord(d))<fmin(xcord(a),xcord(b)) ||
               fmax(ycord(a),ycord(b))<fmin(ycord(c),ycord(d)) ||
               fmax(ycord(c),ycord(d))<fmin(ycord(a),ycord(b)))continue;
            shared=(a==c || a==d)+(b==c || b==d);
            if(shared==2)return FALSE;
            s1=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c));
            s2=orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(d));
            if(shared==1 && (s1 || s2))continue;
            s3=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(a));
            s4=orient2d_float_sign(cord_arr(c),cord_arr(d),cord_arr(b));
            if(s1*s2<=0 && s3*s4<=0)return FALSE;
        }
    }
    return TRUE;
}

/* Triangulate only the validated cavity interior. Strict exact ear tests
 * preserve collinear boundary vertices and reject diagonals containing another
 * vertex. All output vertices already exist; no snapping or field transfer. */
static int removal_ears(int *nodes,int count,int *ears)
{
    int *active=malloc((size_t)count*sizeof(int));
    int i,j,k,n=count,nt=0,a,b,c,p,found,inside;
    if(!active)panic("out of storage triangulating removal cavity");
    for(i=0;i<count;i++)active[i]=nodes[i];
    while(n>3) {
        found=FALSE;
        for(i=0;i<n;i++) {
            a=pt_nd(active[(i+n-1)%n]);b=pt_nd(active[i]);c=pt_nd(active[(i+1)%n]);
            if(orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c))<=0)continue;
            inside=FALSE;
            for(j=0;j<n;j++) {
                if(j==i || j==(i+n-1)%n || j==(i+1)%n)continue;
                p=pt_nd(active[j]);
                if(orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(p))>=0 &&
                   orient2d_float_sign(cord_arr(b),cord_arr(c),cord_arr(p))>=0 &&
                   orient2d_float_sign(cord_arr(c),cord_arr(a),cord_arr(p))>=0) {inside=TRUE;break;}
            }
            if(inside)continue;
            ears[3*nt]=active[(i+n-1)%n];ears[3*nt+1]=active[i];ears[3*nt+2]=active[(i+1)%n];nt++;
            for(k=i;k<n-1;k++)active[k]=active[k+1];n--;found=TRUE;break;
        }
        if(!found){free(active);return -1;}
    }
    a=pt_nd(active[0]);b=pt_nd(active[1]);c=pt_nd(active[2]);
    if(orient2d_float_sign(cord_arr(a),cord_arr(b),cord_arr(c))<=0){free(active);return -1;}
    for(i=0;i<3;i++)ears[3*nt+i]=active[i];
    free(active);return nt+1;
}

#ifdef SUPREM4_TEST_REMOVAL
int test_removal_ears(int *nodes,int count,int *ears)
{return removal_ears(nodes,count,ears);}
#endif

static int removal_cavity(in,apply)
int in,apply;
{
    int nt=num_tri_nd(in), count=0,i,j,k,t,v,r,sr,e;
    int first,current,last,next,nnodes,closed,remaining,predecessor;
    int *from,*to,*region,*used,*nodes;
    int valid=TRUE,*ears,ntri;
    if(mode==ONED) {
        if(!apply)return TRUE;
        int a=-1,b=-1;
        for(i=0;i<nt;i++) {
            t=tri_nd(in,i);
            if(dead_tri(t))continue;
            j=vert_tri(t,0)==in ? vert_tri(t,1):vert_tri(t,0);
            if(a<0)a=j;else if(j!=a)b=j;
        }
        if(a<0 || b<0)panic("invalid one-dimensional removal cavity");
        t=tri_nd(in,0);r=reg_tri(t);sr=cr_sreg(mat_nd(in));e=mk_edg(a,b);
        ad_edge(sr,e,sreg[sr]->bnd,MAYBE,BEFORE);hinit(sr,MAXFLOAT);
        grid(sr,r,FALSE,FALSE);free_skel(sr);
        while(num_tri_nd(in)>0)rem_ele(tri_nd(in,0),FALSE);
        return TRUE;
    }
    if(mode!=TWOD)panic("unsupported node removal dimension");
    from=malloc((size_t)(nt+1)*sizeof(int));to=malloc((size_t)(nt+1)*sizeof(int));
    region=malloc((size_t)(nt+1)*sizeof(int));used=calloc((size_t)(nt+1),sizeof(int));
    nodes=malloc((size_t)(nt+2)*sizeof(int));
    if(!from || !to || !region || !used || !nodes)panic("out of storage for removal cavity");
    /* The edges opposite the removed vertex are the directed link of its
     * triangle star. Unlike angle sorting, this preserves concave boundary
     * fans and narrow wedges with nearly equal polar angles. Keep region
     * identities rather than assigning the entire cavity to one region. */
    for(i=0;i<nt;i++) {
        t=tri_nd(in,i);
        if(dead_tri(t) || !leaf(tri[t]))continue;
        for(v=0;v<3 && vert_tri(t,v)!=in;v++);
        if(v==3)panic("broken node removal incidence");
        from[count]=vert_tri(t,(v+1)%3);to[count]=vert_tri(t,(v+2)%3);
        region[count++]=reg_tri(t);
    }
    /* A manifold link has at most one incoming and outgoing edge per region.
     * Reject a malformed star before adding replacement triangles. */
    for(i=0;i<count;i++)for(j=i+1;j<count;j++)
        if(region[i]==region[j] && (from[i]==from[j] || to[i]==to[j])) {
            valid=FALSE;
            goto cleanup;
        }
    remaining=count;
    while(remaining) {
        first=-1;
        /* Prefer an open fan's endpoint; otherwise start a closed cycle. */
        for(i=0;i<count;i++)if(!used[i]) {
            if(first<0)first=i;
            predecessor=FALSE;
            for(j=0;j<count;j++)if(!used[j] && region[j]==region[i] && to[j]==from[i])
                predecessor=TRUE;
            if(!predecessor){first=i;break;}
        }
        r=region[first];current=first;nnodes=0;closed=FALSE;
        nodes[nnodes++]=from[first];
        for(k=0;k<=count;k++) {
            if(used[current])panic("repeated removal cavity link");
            used[current]=TRUE;remaining--;last=to[current];
            if(last==nodes[0]){closed=TRUE;break;}
            nodes[nnodes++]=last;next=-1;
            for(j=0;j<count;j++)if(!used[j] && region[j]==r && from[j]==last){next=j;break;}
            if(next<0)break;
            current=next;
        }
        if(kn_diag_on())fprintf(stderr,"[kn] cavity node=%d region=%d vertices=%d closed=%d\n",in,r,nnodes,closed);
        if(nnodes<3)continue;
        if(!removal_polygon_valid(nodes,nnodes)) {valid=FALSE;break;}
        if(!apply && !closed && !removal_chord_valid(in,nodes[0],nodes[nnodes-1])) {valid=FALSE;break;}
        ears=malloc((size_t)3*(nnodes-2)*sizeof(int));
        if(!ears)panic("out of storage for removal ears");
        ntri=removal_ears(nodes,nnodes,ears);
        if(ntri<0){free(ears);valid=FALSE;break;}
        if(apply)for(i=0;i<ntri;i++) {
            t=mk_ele_nd(3,&ears[3*i],r);
            sclk_tri(t);
        }
        free(ears);
    }
    if(apply && valid) {
        while(num_tri_nd(in)>0)rem_ele(tri_nd(in,0),FALSE);
        clkws_dirty=geom_dirty=neigh_dirty=edge_dirty=TRUE;
    }
cleanup:
    free(from);free(to);free(region);free(used);free(nodes);
    return valid;
}

int node_removal_valid(int in) {return removal_cavity(in,FALSE);}
int kill_node_prevalidated(int in) {return removal_cavity(in,TRUE);}

int kill_node(int in)
{
    /* Refuse an optional coarsening operation whose closing chord crosses
     * the cavity. Validate before replacing any material region. */
    if(!node_removal_valid(in))return FALSE;
    return removal_cavity(in,TRUE);
}


/*-----------------TROTATE----------------------------------------------
 * Tricky little routine to walk around a boundary edge looking
 * for the next node.
 *----------------------------------------------------------------------*/
trotate( ie, iv, clock, Pje, Poje, Pj, Pjj)
     int ie;			/* The triangle */
     int iv;			/* The index of the boundary node */
     int clock;			/* 1 for clockwise, 0 for cc */
     int *Pje;			/* The outside neighbor (loop-breaker) */
     int *Poje;			/* The inside neighbor */
     int *Pj;			/* index of rotation node in oje */
     int *Pjj;			/* index of leading node or trailing triangle */
{
    int j, jj, je, oje, delta, ir=tri[ie]->regnum, in=tri[ie]->nd[iv];
    
    /*
     * Rotate around given node (in) until we
     *  a: come back to where we started    -or-
     *  b: hit an interface.
     *  Loop invariants: tri[oje]->nd[j] is the pivot (in)
     *  		 tri[oje]->nb[jj] precedes oje.
     *                   je succeeds oje
     *  The neighbors must be correct for this to work.
     */

    delta = clock? 1 : 2;
    for (oje=ie, j = iv, jj = (j+delta)%3; ; ) {
	 
	 je = tri[ oje]->nb[ 3-j-jj];
	 if ((je < 0) || (je == ie) || (tri[ je]->regnum != ir)) break;

	 jj = twhich (je, oje);
	 j  = nwhich (je, in);
	 oje = je;
     }

    *Pjj = jj;
    *Pj = j;
    *Poje = oje;
    *Pje = je;
}
          
