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
/*   oxgrow.c                Version 5.1     */
/*   Last Modification : 7/3/91 10:52:33 */
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* SUPREM4_OXGROW_TIMES=1: seconds spent in each stage of an oxidation step */
static struct timespec oxt_t0;
#define OXT(what) do { if (getenv("SUPREM4_OXGROW_TIMES")) { struct timespec _t; clock_gettime(CLOCK_MONOTONIC, &_t); \
	fprintf(stderr, "[oxgrow-time] %s %.2f\n", what, oxt_t0.tv_sec ? (_t.tv_sec - oxt_t0.tv_sec) + 1e-9 * (_t.tv_nsec - oxt_t0.tv_nsec) : 0.0); oxt_t0 = _t; } } while (0)
#include <math.h>
#include <stdio.h>
#include <sys/times.h>

#include <global.h>
#include <constant.h>
#include <geom.h>
#include <material.h>
#include <impurity.h>
#include <matrix.h>

/*-----------------OXGROW-----------------------------------------------
 * Controller for oxide growth.
 * Remember: All nodes must be added before
 * computing velocities, otherwise new nodes won't have any.
 *----------------------------------------------------------------------*/
oxgrow (temp, ornt, oxhow, dt )
    float temp;		/* Processing temperature */
    int ornt;		/* Orientation of substrate */
    int oxhow;		/* Dry or wet oxidation */
    double *dt;		/* Time increment */
{
    extern double total;	/* Total diffusion time */
    struct tms before, after;

    times(&before);
    mesh_embedding_trace("oxide entry");

    /* Geometry core (M0 WP6): events, repair, regeneration and field
     * transfer for the whole structure; no-op unless SUPREM4_GEOMETRY_CORE
     * is set in a build linked with the core (oxide/geometry_core.c). */
    { extern int geometry_core_step(void); if (geometry_core_step()) UpdateSymbolic++; }
    OXT("after core");

    /* Compress model: repair the oxide/poly geometry and remesh when the
     * Lagrangian motion has produced a collision, a collapse or a degenerate
     * element (oxide/oxide_repair.c). */
    { extern int oxide_repair(int); if (oxide_repair(0)) UpdateSymbolic++; }

    /* Compute the velocity vector using whatever oxide growth model */
    oxide_vel (temp, ornt, oxhow, *dt);
    OXT("oxide_vel");
    mesh_embedding_trace("oxide velocities");

    /* With the geometry core active it owns the topology: SUPREM's own loop
     * removal, node insertion and node removal are not used. */
    { extern int geometry_core_enabled(void);
    if (!geometry_core_enabled()) {

    /* Check for boundary loops */
    /* Before grid_add to clean up boundary before decorating with new nodes*/
    if(DetectLoop())
	UpdateSymbolic++;
    mesh_embedding_trace("DetectLoop");

    /* Now drop in any new nodes. Velocities, concentrations as old ones */
    if (grid_add(*dt))
	UpdateSymbolic++;
    mesh_embedding_trace("grid_add");
    } }

    /* Make a pointwise velocity list: do this AFTER np has been bumped */
    point_vel( );
    OXT("point_vel");

    /* Geometry core: no triangle may lose more than half its area within
     * the step (predicted from the point velocities); when that limit is far
     * below the step SUPREM wanted, the core regenerates the squeezed
     * elements at the next step instead of the step collapsing. */
    { extern int geometry_core_enabled(void); extern void geometry_core_request(void);
      if (geometry_core_enabled()) {
	  double want = *dt, lim = *dt, x0, y0, x1, y1, x2, y2, u0, v0, u1, v1, u2, v2, a, b, c, A0, s, t, disc;
	  int ie;
	  for (ie = 0; ie < ne; ie++) {
	      if (dead_tri(ie)) continue;
	      x0 = xcord(pt_nd(vert_tri(ie,0))); y0 = ycord(pt_nd(vert_tri(ie,0)));
	      x1 = xcord(pt_nd(vert_tri(ie,1))); y1 = ycord(pt_nd(vert_tri(ie,1)));
	      x2 = xcord(pt_nd(vert_tri(ie,2))); y2 = ycord(pt_nd(vert_tri(ie,2)));
	      u0 = pt[pt_nd(vert_tri(ie,0))]->vel[0]; v0 = pt[pt_nd(vert_tri(ie,0))]->vel[1];
	      u1 = pt[pt_nd(vert_tri(ie,1))]->vel[0]; v1 = pt[pt_nd(vert_tri(ie,1))]->vel[1];
	      u2 = pt[pt_nd(vert_tri(ie,2))]->vel[0]; v2 = pt[pt_nd(vert_tri(ie,2))]->vel[1];
	      /* doubled signed area A(t) = A0 + b t + a t^2 */
	      A0 = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
	      b = (u1-u0)*(y2-y0) + (x1-x0)*(v2-v0) - (v1-v0)*(x2-x0) - (y1-y0)*(u2-u0);
	      a = (u1-u0)*(v2-v0) - (v1-v0)*(u2-u0);
	      if (A0 == 0) continue;
	      s = (A0 > 0) ? 1 : -1;	/* SUPREM's orientation convention */
	      /* smallest t > 0 with s*A(t) = 0.5*s*A0:  a t^2 + b t + 0.5 A0 = 0 */
	      c = 0.5 * A0;
	      t = HUGE_VAL;
	      if (a == 0) { if (b != 0 && -c / b > 0) t = -c / b; }
	      else {
		  disc = b*b - 4*a*c;
		  if (disc >= 0) {
		      double r1 = (-b - sqrt(disc)) / (2*a), r2 = (-b + sqrt(disc)) / (2*a);
		      if (r1 > 0) t = r1;
		      if (r2 > 0 && r2 < t) t = r2;
		  }
	      }
	      (void)s;
	      if (t < lim) lim = t;
	  }
	  if (lim < *dt) *dt = lim;
	  if (lim < 0.2 * want) geometry_core_request();
	  { extern double geometry_core_collision_dt(double); *dt = geometry_core_collision_dt(*dt); }
	  if (getenv("SUPREM4_CORE_SPEED_DIAG"))
	      fprintf(stderr, "[core-area-limit] dt %.3e -> %.3e%s\n", want, *dt, lim < 0.2 * want ? " (regeneration requested)" : "");
      } }

    /* A nonlinear vertical map must sample neighboring constraints
     * consistently before testing volume-element motion. */
    if(oxide_model==O_VERTICAL)vertical_constraint_refine(*dt);

    /* With the repair active, no point moves more than delta/2 per step,
     * so every collision or collapse is seen by the next repair check. */
    { extern int oxide_repair_enabled(void);
      extern double oxide_repair_delta(void), oxide_repair_max_speed(void);
      if (oxide_repair_enabled()) {
	  double v = oxide_repair_max_speed(), lim;
	  if (getenv("SUPREM4_REPAIR_SPEED_DIAG")) {
	      int n, ng = 0, s = imptosol[gas_type]; double cmax = 0, cmin = 1e30;
	      for (n = 0; n < nn; n++) {
		  if (nd[n]->mater == GAS) ng++;
		  if (nd[n]->mater == SiO2 && s >= 0) { if (nd[n]->sol[s] > cmax) cmax = nd[n]->sol[s]; if (nd[n]->sol[s] < cmin) cmin = nd[n]->sol[s]; }
	      }
	      fprintf(stderr, "[repair-speed] vmax %.3e cm/s, gas nodes %d, oxidant %.3e..%.3e\n", v, ng, cmin, cmax);
	  }
	  if (v > 0 && *dt > (lim = 0.5*oxide_repair_delta()/v)) *dt = lim;
      } }

    /* Optional (SUPREM4_CORE_SPEED_LIMIT=1): no point moves more than
     * delta/4 per step.  Off by default since WP6d: the collision prediction
     * above already limits approaching pairs (relative velocity, so also gaps
     * closing from both sides); a global speed cap only throttles the step
     * wherever the oxide grows fastest. */
    { extern int geometry_core_enabled(void); extern double oxide_repair_max_speed(void);
      const char *sl = getenv("SUPREM4_CORE_SPEED_LIMIT");
      if (geometry_core_enabled() && sl && *sl && strcmp(sl, "0") != 0) {
	  const char *e = getenv("SUPREM4_CORE_DELTA");
	  double delta = 1e-4 * ((e && *e) ? atof(e) : 1e-3), v = oxide_repair_max_speed(), lim;
	  if (getenv("SUPREM4_CORE_SPEED_DIAG")) {
	      int p, pm = -1; double vm = 0, s2;
	      for (p = 0; p < np; p++) if (!dead_pt(p)) { s2 = hypot(pt[p]->vel[0], pt[p]->vel[1]); if (s2 > vm) { vm = s2; pm = p; } }
	      if (pm >= 0) fprintf(stderr, "[core-speed] vmax %.3e um/s at (%.4f, %.4f), dt %.3e -> %.3e\n", 1e4 * vm,
				   1e4 * xcord(pm), 1e4 * ycord(pm), *dt, (v > 0 && *dt > 0.25 * delta / v) ? 0.25 * delta / v : *dt);
	  }
	  if (v > 0 && *dt > (lim = 0.25 * delta / v)) *dt = lim;
      } }

    OXT("limits");
    /*cut the time step back */
    back_pedal( dt );
    OXT("back_pedal");

    /* Mark nodes for later removal */
    { extern int geometry_core_enabled(void); if (!geometry_core_enabled()) ChooseKillNodes( *dt); }

    /* some barfola statistics */
    times(&after);
    print_time("Time for oxidation step", &before, &after);

    return;
}

/*-----------------POINT_VEL--------------------------------------------
 * Turns node velocities into point velocities
 *----------------------------------------------------------------------*/

point_vel()
{
    int ip, nx, ix;
    int arr[3];

    arr[0] = imptosol[XVEL];
    arr[1] = imptosol[YVEL];

    /* Move the points */
    for (ip = 0; ip < np; ip++) {
	if (dead_pt(ip) || num_nd(ip) == 0) continue;

	/*
	 * Well, displacement is a function of node and there are
	 * several nodes at a point, so we have to make a decision.
	 * For now we let the point move with its 0th node and
	 * hope all the other nodes at that point wanted to go to the
	 * same place. Except: when there is a silicon node there,
	 * we make sure to go with it. This means the silicon/oxide
	 * boundary moves into silicon, as wanted.
	 */
	if ( (nx = node_mat( nd_pt(ip,0), Si ) ) == -1 )
	    if ( (nx = node_mat( nd_pt(ip,0), Poly ) ) == -1 )
		    nx = nd_pt(ip, 0);

	for(ix = 0; ix < mode; ix++)
	    pt[ip]->vel[ix] = nd[nx]->sol[arr[ix]];
    }
}

/*------------------MOVE_POINT------------------------------------------
 * Turn point velocities into honest-to-god point displacements.
 *----------------------------------------------------------------------*/
move_point( dt )
     double dt;
{
    register int ip, ix;

    for (ip = 0; ip < np; ip++) {
	for( ix = 0; ix < mode; ix++ )
	    pt[ ip]->cord[ ix] += pt[ip]->vel[ix] * dt;
    }
}

/*-----------------TOT_AREAS--------------------------------------------
 * Integrate the total area in each material (volume check)
 *----------------------------------------------------------------------*/
tot_areas( oar, nar )
double *oar, *nar;
{
    double newa[MAXMAT], olda[MAXMAT], delta;
    int i;

    for( i = 0; i < MAXMAT; i++) newa[i] = olda[i] = 0;

    for( i = 0; i < nn; i++) {
	newa[ nd[i]->mater ] += 1e8*nar[i];
	olda[ nd[i]->mater ] += 1e8*oar[i];
    }

    if ( newa[Si] != olda[Si] )
    printf("Silicon consumed %10.5f oxide grown %10.5f (%5.3f)\n", olda[Si] - newa[Si],
	   newa[SiO2] - olda[SiO2], (newa[SiO2] - olda[SiO2])/(olda[Si] - newa[Si]));

    for (i = 0; i < MAXMAT; i++) {
	if (olda[i] != 0){
	    delta = (newa[i] - olda[i])/olda[i];
	    if (delta > 1e-6 || delta < -1e-6)
		printf("change in %s area %8.5f%%\n", MatNames[i]+1, delta*100);
	}
    }
    return;
}


/*-----------------ClockTri---------------------------------------------
 * Searches the grid for clockwise triangles and (optionally) fixes them.
 * Returns total number of clockwise triangles, negative if some are zero.
 * Original: CSR 08/87.
 *----------------------------------------------------------------------*/
#define X 0
#define Y 1
ClockTri( FixEm)
    int FixEm;
{
    int ie, *n, swap, nclock=0, nzero=0;
    float *c0, *c1, *c2, area, area_tri();

    if ( mode == ONED ) return( 0 );

    for (ie = 0; ie < ne; ie++) {
	n = tri[ie]->nd;
	c0 = pt[ nd[ *(n++)]->pt]->cord;
	c1 = pt[ nd[ *(n++)]->pt]->cord;
	c2 = pt[ nd[ *(n++)]->pt]->cord;
	area = area_tri( c0[X], c0[Y], c1[X], c1[Y], c2[X], c2[Y]);
	if( area == 0) nzero++;
	if( area < 0 ) nclock++;
	if (FixEm && area < 0) {
		n = tri[ie]->nd;
		swap = n[1]; n[1] = n[2]; n[2] = swap;
		n = tri[ie]->nb;
		swap = n[1]; n[1] = n[2]; n[2] = swap;
	    }
    }
    return( (nzero+nclock) * (nzero? -1 : 1));
}

