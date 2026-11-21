/* Reproduce a rejected moving-mesh step: a valid forward edge flip becomes
 * an exterior diagonal at the old coordinates. Exercise the real flipper,
 * geometry routines, list growth, and transactional rollback. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "geom.h"
#include "matrix.h"
#include "impurity.h"
#include "material.h"
#include "mesh_trial.h"

void panic(char *message)
{
    fprintf(stderr, "Unexpected panic: %s\n", message);
    abort();
}

static double area(int t)
{
    float *a = pt[nd[tri[t]->nd[0]]->pt]->cord;
    float *b = pt[nd[tri[t]->nd[1]]->pt]->cord;
    float *c = pt[nd[tri[t]->nd[2]]->pt]->cord;
    return ((double)b[0]-a[0])*(c[1]-a[1]) -
           ((double)b[1]-a[1])*(c[0]-a[0]);
}

static void check_lists(void)
{
    int i, j, t, hits;
    for (i = 0; i < nn; i++) {
        for (j = 0; j < nd[i]->ltri.num; j++) {
            t = nd[i]->ltri.list[j];
            assert(nwhich(t, i) < 3);
        }
        for (j = 0; j < nd[i]->edg.num; j++) {
            t = nd[i]->edg.list[j];
            assert(edg[t]->nd[0] == i || edg[t]->nd[1] == i);
        }
    }
    for (t = 0; t < ne; t++) for (i = 0; i < 3; i++) {
        hits = 0;
        for (j = 0; j < nd[tri[t]->nd[i]]->ltri.num; j++)
            hits += nd[tri[t]->nd[i]]->ltri.list[j] == t;
        assert(hits == 1);
        hits = 0;
        for (j = 0; j < edg[tri[t]->edg[i]]->ele.num; j++)
            hits += edg[tri[t]->edg[i]]->ele.list[j] == t;
        assert(hits == 1);
    }
}

int main(void)
{
    struct pt_str points[6] = {0}, *points_ptr[6];
    struct nd_str nodes[6] = {0}, *nodes_ptr[6];
    struct tri_str triangles[2] = {0}, original[2], *triangles_ptr[2];
    struct edg_str edges[5] = {0}, *edges_ptr[5];
    struct mesh_trial *trial;
    int i, j, iteration;
    const float coordinates[4][2] = {{0,0}, {2,0}, {.5,.5}, {0,2}};
    const int vertices[2][3] = {{0,1,2}, {0,2,3}};
    const int triangle_edges[2][3] = {{2,0,1}, {3,4,0}};
    const int endpoints[5][2] = {{0,2}, {0,1}, {1,2}, {2,3}, {3,0}};
    pt = points_ptr; nd = nodes_ptr; tri = triangles_ptr; edg = edges_ptr;
    np = nn = 4; ne = 2; ned = 5; mode = TWOD; nvrt = nedg = 3;
    for (i = 0; i < 6; i++) {
        pt[i] = &points[i]; nd[i] = &nodes[i]; nd[i]->pt = i;
        if(i<4)memcpy(pt[i]->cord, coordinates[i], sizeof(coordinates[i]));
    }
    for (i = 0; i < 5; i++) {
        edg[i] = &edges[i];
        memcpy(edg[i]->nd, endpoints[i], sizeof(endpoints[i]));
        for (j = 0; j < 2; j++) add_edge_nd(endpoints[i][j], i);
    }
    for (i = 0; i < 2; i++) {
        tri[i] = &triangles[i];
        tri[i]->fath=tri[i]->son=NO_TRIANGLE;
        for (j = 0; j < 3; j++) {
            tri[i]->nd[j] = vertices[i][j];
            tri[i]->edg[j] = triangle_edges[i][j];
            tri[i]->nb[j] = triangle_edges[i][j] == 0 ? 1-i : -1;
            add_tri_nd(vertices[i][j], i);
            add_ele_edg(triangle_edges[i][j], i);
        }
        sclk_tri(i); do_geom(i); original[i] = *tri[i];
    }
    /* PMOS activation: the old slope/intercept predicate reports both
     * endpoints on the same side and leaves this removable surface sliver.
     * Coordinates captured before motion at t=2.61076 s, native64_rjk2jizs.
     * Two triangles share the diagonal between local points 0 and 2. */
    {
        const float quad[4][2] = {
            {0x1.a5f26ap-12f, 0x1.817120p-16f},
            {0x1.a79d72p-12f, 0x1.816c52p-16f},
            {0x1.a2acecp-12f, 0x1.817aacp-16f},
            {0x1.a522aap-12f, 0x1.81737ep-16f}
        };
        double before_min, after_min;
        trial = mesh_trial_save();
        for (i = 0; i < 4; i++) memcpy(pt[i]->cord,quad[i],sizeof quad[i]);
        assert(area(0)>0 && area(1)>0);
        before_min = fmin(area(0),area(1));
        assert(ck_flip(0) == 1);
        tri_flip(0,0,1);
        after_min = fmin(area(0),area(1));
        assert(after_min > before_min*10);
        check_lists();
        mesh_trial_restore(trial);mesh_trial_discard(trial);
        assert(ck_flip(0) == 0); /* Concave quad: exterior new diagonal. */
    }
    for (iteration = 0; iteration < 10; iteration++) {
        trial = mesh_trial_save();
        pt[2]->cord[0] = pt[2]->cord[1] = 2;
        tri_flip(0, 0, 1);
        assert(area(0) > 0 && area(1) > 0);
        check_lists();
        /* The old coordinate-only rollback leaves an inverted triangle. */
        pt[2]->cord[0] = pt[2]->cord[1] = .5;
        assert(area(0) < 0 || area(1) < 0);
        /* Force buffers to move: a snapshot must not restore old pointers. */
        for (i = 0; i < 100; i++) add_tri_nd(0, 123);
        UpdateSymbolic = 0;
        mesh_trial_restore(trial);
        assert(UpdateSymbolic);
        check_lists();
        for (i = 0; i < 2; i++) {
            assert(memcmp(tri[i], &original[i], sizeof(original[i])) == 0);
            assert(area(i) > 0);
            do_geom(i);
        }
        mesh_trial_discard(trial);
    }
    /* Committing a successful trial retains the new valid geometry. */
    trial = mesh_trial_save();
    pt[2]->cord[0] = pt[2]->cord[1] = 2;
    tri_flip(0, 0, 1);
    mesh_trial_discard(trial);
    assert(area(0) > 0 && area(1) > 0);
    assert(edg[0]->nd[0] != 0 && edg[0]->nd[1] != 0);
    check_lists();
    /* A connectivity change must conserve a non-affine chemical field.
     * Exercise the real primitive and its rollback, not just the projection.
     * An affine field must retain its exact nodal values. */
    {
        double before=0,after=0,old[4]={1,7,3,9};
        pt[2]->cord[0]=2;pt[2]->cord[1]=1;
        n_imp=2;soltoimp[0]=B;soltoimp[1]=P;
        for(i=0;i<4;i++) {
            nd[i]->mater=Poly;
            nd[i]->sol[0]=old[i];
            nd[i]->sol[1]=10+pt[i]->cord[0]+2*pt[i]->cord[1];
        }
        for(i=0;i<2;i++)for(j=0;j<3;j++)before+=area(i)/6*nd[tri[i]->nd[j]]->sol[0];
        tri_flip(0,0,1);check_lists();
        for(i=0;i<2;i++)for(j=0;j<3;j++)after+=area(i)/6*nd[tri[i]->nd[j]]->sol[0];
        assert(fabs(after-before)<1e-12*before);
        for(i=0;i<4;i++) {
            assert(nd[i]->sol[0]>=1 && nd[i]->sol[0]<=9);
            assert(nd[i]->sol[1]==10+pt[i]->cord[0]+2*pt[i]->cord[1]);
            old[i]=nd[i]->sol[0];
        }
        /* The implicit integrator owns its history fields and mass terms.
         * A trial flip cannot rewrite the committed input fields. */
        trial=mesh_trial_save();assert(mesh_trial_active());
        assert(tri_flip(0,0,1));
        for(i=0;i<4;i++)assert(nd[i]->sol[0]==old[i]);
        mesh_trial_restore(trial);mesh_trial_discard(trial);
        assert(!mesh_trial_active());check_lists();
        n_imp=0;
    }
    /* Binary32 motion flips a positive captured triangle even though the
     * unrounded motion increases its area. Rejection must be atomic. */
    {
        const float coords[3][2] = {
            {0x1.a5f26ap-12f,0x1.817120p-16f},
            {0x1.a2acecp-12f,0x1.817aacp-16f},
            {0x1.a522aap-12f,0x1.81737ep-16f}
        };
        const float velocity[3]={-0x1.720098p-35f,-0x1.760454p-35f,-0x1.72ff8ap-35f};
        float saved[3][MAXDIM];int bad;
        ne=1;np=3;
        for(i=0;i<3;i++) {
            tri[0]->nd[i]=i;memcpy(pt[i]->cord,coords[i],sizeof coords[i]);
            memcpy(saved[i],pt[i]->cord,sizeof saved[i]);
            pt[i]->vel[0]=0;pt[i]->vel[1]=velocity[i];
        }
        assert(area(0)>0);
        assert(!mesh_trial_move(0.96590769513979402,&bad));assert(bad==0);
        for(i=0;i<3;i++)assert(!memcmp(saved[i],pt[i]->cord,sizeof saved[i]));
        assert(mesh_trial_move(0.48295384756989701,&bad));assert(area(0)>0);
        for(i=0;i<3;i++)assert(pt[i]->cord[1]!=(float)saved[i][1]);
        for(i=0;i<3;i++)memcpy(saved[i],pt[i]->cord,sizeof saved[i]);
        pt[2]->vel[1]=INFINITY;
        assert(!mesh_trial_move(1,&bad));assert(bad==-1);
        for(i=0;i<3;i++)assert(!memcmp(saved[i],pt[i]->cord,sizeof saved[i]));
    }
    /* A translation can preserve every element's area yet collide with a
     * different component. Reject the whole candidate before any writes. */
    {
        const float separated[6][2]={{0,0},{1,0},{0,1},{3,0},{4,0},{3,1}};
        float saved[6][MAXDIM];int bad;
        np=nn=6;ne=2;
        for(i=0;i<6;i++) {
            memcpy(pt[i]->cord,separated[i],sizeof separated[i]);
            memcpy(saved[i],pt[i]->cord,sizeof saved[i]);
            pt[i]->vel[0]=i<3?0:-3;pt[i]->vel[1]=0;
        }
        for(i=0;i<2;i++)for(j=0;j<3;j++) {
            tri[i]->nd[j]=3*i+j;tri[i]->nb[j]=-1022;
        }
        assert(area(0)>0 && area(1)>0);
        assert(!mesh_trial_move(1,&bad));assert(bad>=0);
        for(i=0;i<6;i++)assert(!memcmp(saved[i],pt[i]->cord,sizeof saved[i]));
        assert(mesh_trial_move(.5,&bad));assert(area(0)>0 && area(1)>0);
    }
    for (i = 0; i < nn; i++) { dis_list(&nd[i]->ltri); dis_list(&nd[i]->edg); }
    for (i = 0; i < ned; i++) dis_list(&edg[i]->ele);
    mesh_trial_restore(NULL); mesh_trial_discard(NULL);
    puts("PASS: sliver flip, atomic rounded-motion/collision rejection, smaller valid motion and rollback");
    return 0;
}
