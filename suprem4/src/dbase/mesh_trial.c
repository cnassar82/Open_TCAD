/* Transactional mesh state for rejected moving-boundary diffusion steps.
 * Moving coordinates backward is insufficient: geom() can flip edges in
 * the tentative geometry, and the new diagonal may lie outside the old
 * cavity. Save the original topology and list order as well as positions.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "matrix.h"
#include "mesh_trial.h"

struct mesh_trial {
    int points, nodes, triangles, edges;
    struct pt_str *point;
    struct tri_str *triangle;
    struct edg_str *edge;
    int *lists;
};
static int active_trials;
int mesh_trial_active(void) { return active_trials!=0; }

static void *trial_alloc(size_t count, size_t size)
{
    void *memory;
    if (count > (size_t)-1 / size) panic("Mesh trial allocation overflow");
    memory = malloc((count ? count : 1) * size);
    if (!memory) panic("Out of memory");
    return memory;
}

static void count_list(size_t *count, const struct list_str *list)
{
    size_t extra = (size_t)list->num + 1;
    if (list->num < 0 || *count > (size_t)-1 - extra)
        panic("Mesh trial connectivity overflow");
    *count += extra;
}

static int *save_list(int *out, const struct list_str *list)
{
    *out++ = list->num;
    if (list->num) memcpy(out, list->list, list->num * sizeof(*out));
    return out + list->num;
}

static const int *restore_list(const int *in, struct list_str *list)
{
    int count = *in++, i;
    list->num = 0;
    for (i = 0; i < count; i++) add_list(list, in[i]);
    return in + count;
}

struct mesh_trial *mesh_trial_save(void)
{
    struct mesh_trial *trial = trial_alloc(1, sizeof(*trial));
    size_t count = 0;
    int i, *out;
    trial->points = np; trial->nodes = nn;
    trial->triangles = ne; trial->edges = ned;
    trial->point = trial_alloc(np, sizeof(*trial->point));
    trial->triangle = trial_alloc(ne, sizeof(*trial->triangle));
    trial->edge = trial_alloc(ned, sizeof(*trial->edge));
    for (i = 0; i < np; i++) trial->point[i] = *pt[i];
    for (i = 0; i < ne; i++) trial->triangle[i] = *tri[i];
    for (i = 0; i < nn; i++) {
        count_list(&count, &nd[i]->ltri);
        count_list(&count, &nd[i]->edg);
    }
    for (i = 0; i < ned; i++) {
        trial->edge[i] = *edg[i];
        count_list(&count, &edg[i]->ele);
    }
    out = trial->lists = trial_alloc(count, sizeof(*out));
    for (i = 0; i < nn; i++) {
        out = save_list(out, &nd[i]->ltri);
        out = save_list(out, &nd[i]->edg);
    }
    for (i = 0; i < ned; i++) out = save_list(out, &edg[i]->ele);
    active_trials++;
    return trial;
}

void mesh_trial_restore(const struct mesh_trial *trial)
{
    int i;
    const int *in;
    if (!trial) return;
    if (np != trial->points || nn != trial->nodes ||
        ne != trial->triangles || ned != trial->edges)
        panic("Mesh entity count changed during diffusion trial");
    for (i = 0; i < np; i++) *pt[i] = trial->point[i];
    for (i = 0; i < ne; i++) *tri[i] = trial->triangle[i];
    in = trial->lists;
    for (i = 0; i < nn; i++) {
        in = restore_list(in, &nd[i]->ltri);
        in = restore_list(in, &nd[i]->edg);
    }
    for (i = 0; i < ned; i++) {
        /* List buffers may have grown during flips. Keep their current
         * allocations; never restore shallow copies of those pointers. */
        struct list_str elements = edg[i]->ele, skeleton = edg[i]->skel;
        *edg[i] = trial->edge[i];
        edg[i]->ele = elements;
        edg[i]->skel = skeleton;
        in = restore_list(in, &edg[i]->ele);
    }
    /* The failed trial may have rebuilt the matrix for its flipped mesh. */
    UpdateSymbolic = TRUE;
}

void mesh_trial_discard(struct mesh_trial *trial)
{
    if (!trial) return;
    active_trials--;
    free(trial->point);
    free(trial->triangle);
    free(trial->edge);
    free(trial->lists);
    free(trial);
}

/* Validate the actual stored-coordinate candidate, not just a polynomial in
 * unrounded displacements. Thin triangles can invert solely on float storage.
 * Preparation has no writes or allocation; commit uses the same rounding. */
static float moved_coordinate(int p, int axis, double dt)
{
    return (float)((double)pt[p]->cord[axis] + pt[p]->vel[axis]*dt);
}

int mesh_trial_move(double dt, int *bad_triangle)
{
    int p,k,t,j;
    float candidate[3][2];
    *bad_triangle=-1;
    if(!isfinite(dt))return 0;
    for(p=0;p<np;p++)for(k=0;k<mode;k++)
        if(!isfinite(pt[p]->cord[k]) || !isfinite(pt[p]->vel[k]) ||
           !isfinite(moved_coordinate(p,k,dt)))return 0;
    for(t=0;t<ne;t++) {
        if(dead_tri(t) || !leaf(tri[t]))continue;
        for(j=0;j<num_vert(t);j++) {
            p=pt_nd(vert_tri(t,j));
            for(k=0;k<mode;k++)candidate[j][k]=moved_coordinate(p,k,dt);
        }
        if((mode==TWOD && orient2d_float_sign(candidate[0],candidate[1],candidate[2])<=0) ||
           (mode==ONED && candidate[0][0]<=candidate[1][0])) {
            *bad_triangle=t;return 0;
        }
    }
    if(getenv("SUPREM4_EMBED_TRACE") && dt!=0 && mode==TWOD) {
        double lo,hi;
        FILE *trace;
        if(sscanf(getenv("SUPREM4_EMBED_TRACE"),"%lf,%lf",&lo,&hi)==2 &&
           (trace=fopen("embedding_motion.dat","w"))) {
            fprintf(trace,"dt %.17g\n",dt);
            for(p=0;p<np;p++)if(pt[p]->cord[0]>=lo*1e-4 && pt[p]->cord[0]<=hi*1e-4)
                fprintf(trace,"%d %.17g %.17g %.17g %.17g %.17g %.17g\n",p,
                    (double)pt[p]->cord[0],(double)pt[p]->cord[1],
                    pt[p]->vel[0],pt[p]->vel[1],
                    (double)pt[p]->cord[0]+pt[p]->vel[0]*dt,
                    (double)pt[p]->cord[1]+pt[p]->vel[1]*dt);
            fclose(trace);
        }
    }
    for(p=0;p<np;p++)for(k=0;k<mode;k++)pt[p]->cord[k]=moved_coordinate(p,k,dt);
    return 1;
}
