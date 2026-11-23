/*
 * Geometry core hook (M0 WP6, suprem4/doc/plan_m0_geometry_core.md).
 *
 * Active only when the environment variable SUPREM4_GEOMETRY_CORE is set
 * (and not "0") AND the executable was linked with the core library
 * (build_native64.py --geometry-core); otherwise this is a no-op, so builds
 * and runs without it are unchanged.
 *
 * While active the core owns the topology: oxgrow skips SUPREM's DetectLoop,
 * grid_add and ChooseKillNodes; SUPREM still moves the points inside its
 * moving-mesh time step.
 *
 * At the start of every oxidation step the current structure is written,
 * the core checks it for geometric events (thin / crumb / gap / cavity of
 * size SUPREM4_CORE_DELTA um, default 1 nm) and for broken box-method
 * couplings (worst cot-sum below SUPREM4_CORE_MIN_COUPLING, default -0.5);
 * if needed it repairs the geometry, regenerates the mesh and transfers the
 * fields (dopant inventories conserved), and the result is reloaded exactly
 * like "structure infile=".  The velocity history slots of do_oxide()
 * (sol[MAXIMP-5], sol[MAXIMP-6]) are not in the file; at this point they equal
 * the current velocities, so they are restored from XVEL / YVEL.
 *
 * Test switches: SUPREM4_CORE_RELOAD_ONLY=1 (write and reload the structure
 * unchanged, no core call), SUPREM4_CORE_EVERY=N (force a regeneration every
 * N calls), SUPREM4_CORE_EVENT_DELTA (event size only, um; the motion
 * limit keeps SUPREM4_CORE_DELTA).
 * Other settings: SUPREM4_CORE_TMP (directory for the exchange files,
 * default "."), SUPREM4_CORE_FORCE_FIRST=1 (regenerate on the first call),
 * SUPREM4_CORE_HMAX, SUPREM4_CORE_SCALE (mesh size, read by the core).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "suprem_options.h"
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "dbaccess.h"
#include "material.h"
#include "impurity.h"

extern int suprem_core_oxidation_check(const char *in_path, const char *out_path, double delta_um,
				       double min_coupling, int force, char *msg, int msglen)
    __attribute__((weak));
extern int ig2_read(), ig2_write();

static double env_double(const char *name, double dflt)
{
    const char *v = getenv(name);
    return (v && *v) ? atof(v) : dflt;
}

static int requested = 0;
/* Ask for a regeneration at the next step (elements about to be squeezed). */
void geometry_core_request(void) { requested = 1; }

int geometry_core_enabled(void)
{
    const char *e = getenv("SUPREM4_GEOMETRY_CORE");
    return e && *e && strcmp(e, "0") != 0;
}

/* Returns 1 when the mesh was replaced (the caller must refactor). */
int geometry_core_step(void)
{
    static int calls = 0, warned = 0;
    char in[1024], out[1024], msg[4096];
    const char *dir;
    int rc, i, xv, yv, force;

    if (!geometry_core_enabled()) return 0;
    if (!suprem_core_oxidation_check) {
	if (!warned++) fprintf(stderr, "[geometry-core] SUPREM4_GEOMETRY_CORE is set but this build has no geometry core\n");
	return 0;
    }
    dir = getenv("SUPREM4_CORE_TMP");
    if (!dir || !*dir) dir = ".";
    snprintf(in, sizeof in, "%s/.geometry_core_in_%d.str", dir, (int)getpid());
    snprintf(out, sizeof out, "%s/.geometry_core_out_%d.str", dir, (int)getpid());
    force = (calls == 0 && env_double("SUPREM4_CORE_FORCE_FIRST", 0) != 0);
    {
	int every = (int)env_double("SUPREM4_CORE_EVERY", 0);
	if (every > 0 && calls % every == 0) force = 1;
    }
    if (requested) { force = 1; requested = 0; }
    calls++;
    if (ig2_write(in, 0, 1.0) < 0) {
	fprintf(stderr, "[geometry-core] cannot write %s\n", in);
	return 0;
    }
    if (env_double("SUPREM4_CORE_RELOAD_ONLY", 0) != 0) {
	if (!force) { unlink(in); return 0; }
	if (ig2_read(in, 0, 1.0) < 0) panic("geometry core: cannot reload");
	create_db(nn == 0);
	unlink(in);
	xv = imptosol[XVEL];
	yv = imptosol[YVEL];
	if (xv >= 0 && yv >= 0)
	    for (i = 0; i < nn; i++) {
		nd[i]->sol[MAXIMP - 5] = nd[i]->sol[xv];
		nd[i]->sol[MAXIMP - 6] = nd[i]->sol[yv];
	    }
	fprintf(stderr, "[geometry-core] call %d: reloaded unchanged (test)\n", calls);
	return 1;
    }
    msg[0] = 0;
    rc = suprem_core_oxidation_check(in, out, env_double("SUPREM4_CORE_EVENT_DELTA", env_double("SUPREM4_CORE_DELTA", 1e-3)),
				     env_double("SUPREM4_CORE_MIN_COUPLING", -0.5), force, msg, (int)sizeof msg);
    unlink(in);
    if (rc < 0) {
	fprintf(stderr, "[geometry-core] call %d: %s; continuing with the current mesh\n", calls, msg);
	return 0;
    }
    if (rc == 0) return 0;
    if (ig2_read(out, 0, 1.0) < 0) panic("geometry core: cannot reload the regenerated mesh");
    create_db(nn == 0);
    unlink(out);
    xv = imptosol[XVEL];
    yv = imptosol[YVEL];
    if (xv >= 0 && yv >= 0)
	for (i = 0; i < nn; i++) {
	    nd[i]->sol[MAXIMP - 5] = nd[i]->sol[xv];
	    nd[i]->sol[MAXIMP - 6] = nd[i]->sol[yv];
	}
    fprintf(stderr, "[geometry-core] call %d: mesh regenerated: %s\n", calls, msg);
    {   /* SUPREM's own embedding check on the reloaded mesh */
	extern int mesh_trial_validate_embedding(double, int *);
	int bad = -1, j;
	if (!mesh_trial_validate_embedding(0.0, &bad)) {
	    fprintf(stderr, "[geometry-core] call %d: reloaded mesh FAILS SUPREM embedding check at triangle %d", calls, bad);
	    if (bad >= 0 && bad < ne)
		for (j = 0; j < 3; j++) fprintf(stderr, " (%.7f,%.7f)", 1e4 * xcord(pt_nd(vert_tri(bad, j))), 1e4 * ycord(pt_nd(vert_tri(bad, j))));
	    fprintf(stderr, "\n");
	    if (getenv("SUPREM4_CORE_KEEP_BAD")) ig2_write("core_embedding_fail.str", 0, 1.0);
	}
    }
    return 1;
}

/* Continuous collision limit (core-owned geometry): for every exposed
 * boundary vertex p and exposed edge (a,b) not attached to p, the clearance
 * d0 = dist(p, ab) may shrink by at most half within the step under the
 * current point velocities.  Approaching surfaces therefore never touch: they
 * converge until the clearance is below delta, where the core's gap event
 * merges them.  Pairs already closer than delta/2 are left to the core
 * (a regeneration is requested).  Returns the limited step. */
double geometry_core_collision_dt(double dt)
{
    const char *e = getenv("SUPREM4_CORE_DELTA");
    double delta = 1e-4 * ((e && *e) ? atof(e) : 1e-3);   /* cm */
    int t, j, ne_x = 0, i, k, close = 0;
    int (*edge)[2];
    double lim = dt, vmax = 0, reach;
    if (!geometry_core_enabled() || dt <= 0) return dt;
    edge = malloc(sizeof(*edge) * (size_t)(3 * ne + 1));
    if (!edge) return dt;
    for (t = 0; t < ne; t++) {
	if (dead_tri(t)) continue;
	for (j = 0; j < 3; j++)
	    if (tri[t]->nb[j] == EXPOSED) {
		edge[ne_x][0] = pt_nd(vert_tri(t, (j + 1) % 3));
		edge[ne_x][1] = pt_nd(vert_tri(t, (j + 2) % 3));
		ne_x++;
	    }
    }
    for (i = 0; i < ne_x; i++)
	for (k = 0; k < 2; k++) {
	    int p = edge[i][k];
	    double v = hypot(pt[p]->vel[0], pt[p]->vel[1]);
	    if (v > vmax) vmax = v;
	}
    reach = 2 * vmax * dt + 4 * delta;   /* pairs farther apart cannot matter */
    for (i = 0; i < ne_x; i++) {
	int a = edge[i][0], b = edge[i][1];
	double ax = xcord(a), ay = ycord(a), bx = xcord(b), by = ycord(b);
	double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
	double x0 = fmin(ax, bx) - reach, x1 = fmax(ax, bx) + reach, y0 = fmin(ay, by) - reach, y1 = fmax(ay, by) + reach;
	if (l2 <= 0) continue;
	for (j = 0; j < ne_x; j++)
	    for (k = 0; k < 2; k++) {
		int p = edge[j][k];
		double px, py, u, qx, qy, d0, nx, ny, vqx, vqy, closing, tp;
		if (p == a || p == b) continue;
		px = xcord(p); py = ycord(p);
		if (px < x0 || px > x1 || py < y0 || py > y1) continue;
		u = ((px - ax) * ex + (py - ay) * ey) / l2;
		/* the nearest point is an endpoint joined to p by an edge: a corner, not an approach */
		if (u <= 0 || u >= 1) {
		    int q = (u <= 0) ? a : b, m, joined = 0;
		    for (m = 0; m < ne_x && !joined; m++)
			joined = (edge[m][0] == p && edge[m][1] == q) || (edge[m][0] == q && edge[m][1] == p);
		    if (joined) continue;
		    u = (u <= 0) ? 0 : 1;
		}
		qx = ax + u * ex; qy = ay + u * ey;
		d0 = hypot(px - qx, py - qy);
		if (d0 <= 0 || d0 > reach) continue;
		if (d0 < 0.5 * delta) { close++; continue; }
		nx = (px - qx) / d0; ny = (py - qy) / d0;
		vqx = (1 - u) * pt[a]->vel[0] + u * pt[b]->vel[0];
		vqy = (1 - u) * pt[a]->vel[1] + u * pt[b]->vel[1];
		closing = -((pt[p]->vel[0] - vqx) * nx + (pt[p]->vel[1] - vqy) * ny);
		if (closing <= 0) continue;
		tp = 0.5 * d0 / closing;
		if (tp < lim) lim = tp;
	    }
    }
    free(edge);
    if (lim < 0.2 * dt || close) geometry_core_request();
    if (getenv("SUPREM4_CORE_SPEED_DIAG"))
	fprintf(stderr, "[core-collision] dt %.3e -> %.3e, %d pairs closer than delta/2\n", dt, lim, close);
    return lim;
}
