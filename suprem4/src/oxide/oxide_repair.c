/* Geometric repair and remesh for oxidation with the compress model.
 *
 * Lagrangian oxidation of polysilicon changes the topology of the material
 * boundaries: the consumed poly collapses at convex corners and the oxide
 * grown on a gate sidewall meets the oxide on the substrate.  Local mesh
 * surgery cannot handle these events robustly.  Following the approach of
 * commercial simulators (move, repair the geometry, remesh), this module:
 *
 *   1. extracts the boundary loops of every region,
 *   2. repairs the oxide and polysilicon polygons with exact Boolean
 *      operations (refine/oxide_repair.cpp: opening of poly, closing of
 *      oxide, both by delta), leaving every other material untouched,
 *   3. deletes the oxide and polysilicon regions and rebuilds them from the
 *      repaired loops with the skeleton triangulator, reusing the existing
 *      interface points of untouched regions and splitting their boundary
 *      edges wherever a repaired loop ends inside one,
 *   4. interpolates all fields from a snapshot of the old regions and
 *      restores the dopant inventory of each material,
 *   5. validates the new mesh.
 *
 * The repair runs at the start of an oxidation step when it would change the
 * geometry or when an oxide/poly triangle has degenerated; the time step is
 * limited so that no point moves more than delta/2, so every collision or
 * collapse is caught before it can invert an element.
 *
 * Options: SUPREM4_OXIDE_REPAIR (default on; "0" disables),
 *          SUPREM4_REPAIR_DELTA_UM (default 0.002 um),
 *          SUPREM4_REPAIR_QUALITY (default 1e-3: degenerate elements only;
 *            thin oxide layers are legitimately anisotropic),
 *          SUPREM4_MESH_DIAG (report each repair).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "suprem_options.h"
#include "global.h"
#include "constant.h"
#include "geom.h"
#include "dbaccess.h"
#include "material.h"
#include "impurity.h"
#include "skel.h"
#include "matrix.h"

#ifdef SUPREM_DOMAIN_OFFSET
extern int suprem_oxide_repair(int, const int *, const int *, const double *,
    int, int, double, double, double, int *, int **, int **, double **);
#endif
extern int mesh_trial_validate_embedding(double, int *);
#ifdef SUPREM_DOMAIN_OFFSET
extern int suprem_cdt(int, int, const double *, int **);
#endif

#define ACTIVE_MAT(M) ((M) == SiO2 || (M) == Poly)
#define SNAP 3.0e-11			/* cm: Boolean grid is 1e-11 cm */

int oxide_repair_enabled(void)
{
    const char *o = suprem_option("SUPREM4_OXIDE_REPAIR");
#ifdef SUPREM_DOMAIN_OFFSET
    return mode == TWOD && oxide_model == O_ELASTIC && (!getenv("SUPREM4_OXIDE_REPAIR") || o);
#else
    return 0;
#endif
}
double oxide_repair_delta(void)
{
    const char *d = getenv("SUPREM4_REPAIR_DELTA_UM");
    return 1e-4 * ((d && *d) ? atof(d) : 0.002);
}

static double orient(double ax, double ay, double bx, double by, double cx, double cy)
{
    return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
}
static double tri_quality(int t)
{
    double x[3], y[3], a, l = 0;
    int j;
    for (j = 0; j < 3; j++) { x[j] = xcord(pt_nd(vert_tri(t,j))); y[j] = ycord(pt_nd(vert_tri(t,j))); }
    a = orient(x[0], y[0], x[1], y[1], x[2], y[2]);
    for (j = 0; j < 3; j++) { double dx = x[(j+1)%3]-x[j], dy = y[(j+1)%3]-y[j]; l += dx*dx + dy*dy; }
    return l > 0 ? 2*sqrt(3.0)*a/l : 0;		/* 1 equilateral, <= 0 inverted */
}

/*----------------------------------------------------------------------*
 * Region boundary loops (region on the left), as point indices.
 *----------------------------------------------------------------------*/
struct loops { int n, alloc, *mat, *start, *len, *pts, npts, apts; };
static void loops_push_pt(struct loops *L, int p)
{
    if (L->npts == L->apts) { L->apts = L->apts ? 2*L->apts : 1024; L->pts = realloc(L->pts, sizeof(int)*L->apts); }
    L->pts[L->npts++] = p;
}
static void loops_begin(struct loops *L, int mat)
{
    if (L->n == L->alloc) {
	L->alloc = L->alloc ? 2*L->alloc : 64;
	L->mat = realloc(L->mat, sizeof(int)*L->alloc);
	L->start = realloc(L->start, sizeof(int)*L->alloc);
	L->len = realloc(L->len, sizeof(int)*L->alloc);
    }
    L->mat[L->n] = mat; L->start[L->n] = L->npts; L->len[L->n] = 0; L->n++;
}
static void loops_free(struct loops *L) { free(L->mat); free(L->start); free(L->len); free(L->pts); }

static int extract_loops(struct loops *L)
{
    int r, i, t, j, n, a, b, c, ne_b, k, cur, first, *ea, *eb, *used, steps;
    memset(L, 0, sizeof(*L));
    for (r = 0; r < nreg; r++) {
	if (reg[r] == NULL || mat_reg(r) == GAS || num_tri_reg(r) == 0) continue;
	ea = salloc(int, 3*num_tri_reg(r)); eb = salloc(int, 3*num_tri_reg(r)); ne_b = 0;
	for (i = 0; i < num_tri_reg(r); i++) {
	    t = tri_reg(r, i);
	    if (dead_tri(t)) continue;
	    for (j = 0; j < 3; j++) {
		n = neigh_fc(t, j);
		if (n >= 0 && !dead_tri(n) && reg_tri(n) == r) continue;
		a = pt_nd(vert_tri(t,(j+1)%3)); b = pt_nd(vert_tri(t,(j+2)%3)); c = pt_nd(vert_tri(t,j));
		if (orient(xcord(a), ycord(a), xcord(b), ycord(b), xcord(c), ycord(c)) < 0) { k = a; a = b; b = k; }
		ea[ne_b] = a; eb[ne_b] = b; ne_b++;
	    }
	}
	used = scalloc(int, ne_b + 1);
	for (first = 0; first < ne_b; first++) {
	    if (used[first]) continue;
	    loops_begin(L, mat_reg(r));
	    cur = first;
	    for (steps = 0; steps <= ne_b; steps++) {
		used[cur] = 1;
		loops_push_pt(L, ea[cur]); L->len[L->n-1]++;
		for (k = 0; k < ne_b; k++) if (!used[k] && ea[k] == eb[cur]) break;
		if (k == ne_b) break;
		cur = k;
	    }
	    if (eb[cur] != L->pts[L->start[L->n-1]]) {
		fprintf(stderr, "[oxide-repair] open boundary loop in region %d\n", r);
		free(ea); free(eb); free(used); return -1;
	    }
	}
	free(ea); free(eb); free(used);
    }
    return 0;
}

/*----------------------------------------------------------------------*
 * Field snapshot of the active (oxide, poly) triangles.
 *----------------------------------------------------------------------*/
struct snap { int nt, ns, *mat, *slot; double *xy, *val; };
static void snapshot(struct snap *S)
{
    int t, j, k, i;
    S->ns = n_imp + 2;
    S->slot = salloc(int, S->ns);
    for (i = 0; i < n_imp; i++) S->slot[i] = i;
    S->slot[n_imp] = MAXIMP-5; S->slot[n_imp+1] = MAXIMP-6;
    S->nt = 0;
    for (t = 0; t < ne; t++) if (!dead_tri(t) && ACTIVE_MAT(mat_reg(reg_tri(t)))) S->nt++;
    S->mat = salloc(int, S->nt + 1); S->xy = salloc(double, 6*S->nt + 1);
    S->val = salloc(double, 3*S->nt*S->ns + 1);
    for (k = 0, t = 0; t < ne; t++) {
	if (dead_tri(t) || !ACTIVE_MAT(mat_reg(reg_tri(t)))) continue;
	S->mat[k] = mat_reg(reg_tri(t));
	for (j = 0; j < 3; j++) {
	    int n = vert_tri(t,j), p = pt_nd(n);
	    S->xy[6*k+2*j] = xcord(p); S->xy[6*k+2*j+1] = ycord(p);
	    for (i = 0; i < S->ns; i++) S->val[(3*k+j)*S->ns+i] = nd[n]->sol[S->slot[i]];
	}
	k++;
    }
}
static void snap_free(struct snap *S) { free(S->slot); free(S->mat); free(S->xy); free(S->val); }

/* Barycentric interpolation inside a same-material triangle, otherwise the
 * nearest same-material vertex (areas newly filled with oxide). */
static void snap_eval(struct snap *S, int mat, double x, double y, double *out)
{
    int k, j, i, best = -1, bj = 0;
    double bd = HUGE_VAL;
    for (k = 0; k < S->nt; k++) {
	double *c = S->xy + 6*k, a, l0, l1, l2;
	if (S->mat[k] != mat) continue;
	a = orient(c[0], c[1], c[2], c[3], c[4], c[5]);
	if (a <= 0) continue;
	l0 = orient(x, y, c[2], c[3], c[4], c[5]) / a;
	l1 = orient(c[0], c[1], x, y, c[4], c[5]) / a;
	l2 = 1 - l0 - l1;
	if (l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9) {
	    for (i = 0; i < S->ns; i++)
		out[i] = l0*S->val[(3*k)*S->ns+i] + l1*S->val[(3*k+1)*S->ns+i] + l2*S->val[(3*k+2)*S->ns+i];
	    return;
	}
	for (j = 0; j < 3; j++) {
	    double d = hypot(x - c[2*j], y - c[2*j+1]);
	    if (d < bd) { bd = d; best = k; bj = j; }
	}
    }
    for (i = 0; i < S->ns; i++) out[i] = best >= 0 ? S->val[(3*best+bj)*S->ns+i] : 0;
}

static int conserved(int imp)
{
    return imp == As || imp == P || imp == Sb || imp == B || imp == Ga ||
	   (imp >= iBe && imp <= iG && (imp & 1));
}
/* P1 integral of solution slot s over triangles of material m. */
static double inventory_mesh(int m, int s)
{
    int t, j; double sum = 0;
    for (t = 0; t < ne; t++) {
	double a, v = 0;
	if (dead_tri(t) || mat_reg(reg_tri(t)) != m) continue;
	a = 0.5*orient(xcord(pt_nd(vert_tri(t,0))), ycord(pt_nd(vert_tri(t,0))),
		       xcord(pt_nd(vert_tri(t,1))), ycord(pt_nd(vert_tri(t,1))),
		       xcord(pt_nd(vert_tri(t,2))), ycord(pt_nd(vert_tri(t,2))));
	for (j = 0; j < 3; j++) v += nd[vert_tri(t,j)]->sol[s];
	sum += a*v/3;
    }
    return sum;
}

/*----------------------------------------------------------------------*
 * Point lookup by coordinates and T-junction resolution.
 *----------------------------------------------------------------------*/
/* Exact match on the 1e-11 cm Boolean grid (vertices carried through the
 * Booleans keep their integer coordinates; nothing is snapped). */
static int find_point(double x, double y)
{
    int p;
    long long X = llround(x * 1e11), Y = llround(y * 1e11);
    for (p = 0; p < np; p++)
	if (!dead_pt(p) && llround((double)xcord(p) * 1e11) == X && llround((double)ycord(p) * 1e11) == Y) return p;
    return -1;
}
/* Vertex strictly inside segment a-b (within SNAP), parameter in (0,1). */
static int on_segment(double x, double y, double ax, double ay, double bx, double by, double *u)
{
    double dx = bx-ax, dy = by-ay, l2 = dx*dx + dy*dy, d;
    if (l2 <= 0) return 0;
    *u = ((x-ax)*dx + (y-ay)*dy) / l2;
    if (*u <= 1e-9 || *u >= 1 - 1e-9) return 0;
    d = fabs(dx*(y-ay) - dy*(x-ax)) / sqrt(l2);
    return d <= SNAP && hypot(x-ax, y-ay) > SNAP && hypot(x-bx, y-by) > SNAP;
}

/* Split every fixed-region boundary edge that has a repaired-loop vertex in
 * its interior, so untouched regions conform to the rebuilt ones. */
static int split_fixed_edges(int nv, const double *vx, const double *vy)
{
    int ie, i, splits = 0, again = 1, guard;
    float c[MAXDIM];
    double u;
    for (guard = 0; again && guard < 100; guard++) {
	again = 0;
	for (ie = 0; ie < ned && !again; ie++) {
	    int a, b, m;
	    if (dead_edg(ie)) continue;
	    m = mat_nd(nd_edg(ie,0));
	    if (ACTIVE_MAT(m) || m == GAS) continue;
	    a = pt_nd(nd_edg(ie,0)); b = pt_nd(nd_edg(ie,1));
	    if (num_nd(a) < 2 || num_nd(b) < 2) continue;	/* not on a boundary */
	    for (i = 0; i < nv; i++) {
		if (!on_segment(vx[i], vy[i], xcord(a), ycord(a), xcord(b), ycord(b), &u)) continue;
		if (find_point(vx[i], vy[i]) >= 0) continue;
		c[0] = vx[i]; c[1] = vy[i];
		{ int newp; (void)sp_gredg(ie, c, &newp); }
		splits++; again = 1;
		break;
	    }
	}
	if (again) bd_connect("oxide repair: split fixed edge");
    }
    return splits;
}

/*----------------------------------------------------------------------*
 * The repair.
 *----------------------------------------------------------------------*/
/* Node of material m at point p, or -1 (works for points without nodes). */
static int pt_node(int p, int m)
{
    int k;
    for (k = 0; k < num_nd(p); k++) if (mat_nd(nd_pt(p,k)) == m) return nd_pt(p,k);
    return -1;
}
static int edge_between(int n, int m)
{
    int k, ie;
    for (k = 0; k < num_edge_nd(n); k++) {
	ie = edge_nd(n,k);
	if (!dead_edg(ie) && (nd_edg(ie,0) == m || nd_edg(ie,1) == m)) return ie;
    }
    return -1;
}

int oxide_repair(int force)
{
#ifndef SUPREM_DOMAIN_OFFSET
    return 0;
#else
    struct loops L;
    struct snap S;
    int i, j, k, t, nout, changed = 0, *omat = NULL, *ocount = NULL, *flat_mat = NULL, *flat_cnt = NULL;
    double *oxy = NULL, *flat_xy = NULL, delta = oxide_repair_delta(), xmin = HUGE_VAL, xmax = -HUGE_VAL, minq = 1;
    double qtrig = getenv("SUPREM4_REPAIR_QUALITY") ? atof(getenv("SUPREM4_REPAIR_QUALITY")) : 1e-3;
    int s, nv, *lstart, *llen, nl, total, *pidx, ninsert = 0, nint = 0, *imat = NULL;
    double *ix = NULL, *iy = NULL, *gas_template = NULL;
    double *vx, *vy, inv_old[MAXIMP][2];

    if (!oxide_repair_enabled()) return 0;
    for (t = 0; t < ne; t++) if (!dead_tri(t) && ACTIVE_MAT(mat_reg(reg_tri(t)))) {
	double q = tri_quality(t); if (q < minq) minq = q; }
    for (i = 0; i < np; i++) if (!dead_pt(i)) {
	if (xcord(i) < xmin) xmin = xcord(i);
	if (xcord(i) > xmax) xmax = xcord(i); }

    if (extract_loops(&L) < 0) { loops_free(&L); return 0; }
    flat_mat = L.mat; flat_cnt = L.len;
    flat_xy = salloc(double, 2*L.npts + 2);
    for (i = 0; i < L.npts; i++) { flat_xy[2*i] = xcord(L.pts[i]); flat_xy[2*i+1] = ycord(L.pts[i]); }
    nout = suprem_oxide_repair(L.n, flat_mat, flat_cnt, flat_xy, SiO2, Poly, delta, xmin, xmax,
			       &changed, &omat, &ocount, &oxy);
    free(flat_xy);
    if (nout < 0) { loops_free(&L); return 0; }
    if (!force && !changed && minq >= qtrig) {
	free(omat); free(ocount); free(oxy); loops_free(&L); return 0;
    }

    if (getenv("SUPREM4_REPAIR_DUMP")) {
	FILE *f = fopen("repair_loops.txt", "w");
	if (f) {
	    for (i = 0, k = 0; i < nout; k += ocount[i], i++) {
		fprintf(f, "loop %d %d\n", omat[i], ocount[i]);
		for (j = 0; j < ocount[i]; j++) fprintf(f, "%.12g %.12g\n", 1e4*oxy[2*(k+j)], 1e4*oxy[2*(k+j)+1]);
	    }
	    for (i = 0; i < L.n; i++) {
		fprintf(f, "input %d %d\n", L.mat[i], L.len[i]);
		for (j = 0; j < L.len[i]; j++) fprintf(f, "%.12g %.12g\n", 1e4*xcord(L.pts[L.start[i]+j]), 1e4*ycord(L.pts[L.start[i]+j]));
	    }
	    fclose(f);
	}
    }
    /* Output loops: drop holes (sealed voids are filled), collect vertices. */
    lstart = salloc(int, nout + 1); llen = salloc(int, nout + 1);
    for (total = 0, nl = 0, i = 0, k = 0; i < nout; k += ocount[i], i++) {
	double a = 0;
	for (j = 0; j < ocount[i]; j++) {
	    double *p = oxy + 2*(k+j), *q = oxy + 2*(k + (j+1)%ocount[i]);
	    a += p[0]*q[1] - q[0]*p[1];
	}
	if (a <= 0) continue;
	lstart[nl] = k; llen[nl] = ocount[i]; omat[nl] = omat[i]; nl++; total += ocount[i];
    }
    /* Vertex set for T-junctions: repaired loops plus fixed boundaries. */
    vx = NULL; vy = NULL; nv = 0;

    snapshot(&S);
    vx = salloc(double, total + 2*np + 16); vy = salloc(double, total + 2*np + 16);
    for (i = 0; i < nl; i++) for (j = 0; j < llen[i]; j++) {
	vx[nv] = oxy[2*(lstart[i]+j)]; vy[nv] = oxy[2*(lstart[i]+j)+1]; nv++; }
    for (s = 0; s < n_imp; s++) {
	inv_old[s][0] = conserved(soltoimp[s]) ? inventory_mesh(SiO2, s) : 0;
	inv_old[s][1] = conserved(soltoimp[s]) ? inventory_mesh(Poly, s) : 0;
    }

    /* Ambient values carried by gas nodes (uniform); surface points of the
     * deleted regions lose their gas nodes and get new ones. */
    gas_template = salloc(double, MAXIMP);
    for (i = 0; i < MAXIMP; i++) gas_template[i] = 0;
    for (i = 0; i < nn; i++) if (!dead_nd(i) && mat_nd(i) == GAS) {
	for (j = 0; j < MAXIMP; j++) gas_template[j] = nd[i]->sol[j];
	break;
    }

    /* Interior points of the active materials (kept as triangulation
     * vertices where they remain inside the repaired regions). */
    ix = salloc(double, np + 1); iy = salloc(double, np + 1); imat = salloc(int, np + 1); nint = 0;
    for (i = 0; i < np; i++) {
	int m;
	if (dead_pt(i) || num_nd(i) != 1 || is_surf(i)) continue;
	m = mat_nd(nd_pt(i,0));
	if (!ACTIVE_MAT(m)) continue;
	ix[nint] = xcord(i); iy[nint] = ycord(i); imat[nint] = m; nint++;
    }

    /* Delete the active regions. */
    for (i = 0; i < nreg; i++) {
	if (reg[i] == NULL || !ACTIVE_MAT(mat_reg(i))) continue;
	while (num_tri_reg(i) > 0) rem_ele(tri_reg(i,0), FALSE);
	fing_reg(i);
    }
    bd_connect("oxide repair: removed oxide/poly");

    /* Untouched regions conform: split their boundary edges at new vertices,
     * and add their boundary points to the T-junction vertex set. */
    {
	int splits = split_fixed_edges(nv, vx, vy), p;
	if (splits && getenv("SUPREM4_MESH_DIAG"))
	    fprintf(stderr, "[oxide-repair] split %d fixed boundary edges\n", splits);
	/* boundary points of the untouched regions (after the removal every
	 * point that bordered oxide/poly carries a node of one material only
	 * but lies on a region boundary edge) */
	{
	    int *bnd = scalloc(int, np + 1);
	    for (t = 0; t < ne; t++) if (!dead_tri(t))
		for (j = 0; j < 3; j++) if (neigh_fc(t,j) < 0 || reg_tri(neigh_fc(t,j)) != reg_tri(t)) {
		    bnd[pt_nd(vert_tri(t,(j+1)%3))] = 1; bnd[pt_nd(vert_tri(t,(j+2)%3))] = 1; }
	    for (p = 0; p < np; p++) if (!dead_pt(p) && bnd[p]) { vx[nv] = xcord(p); vy[nv] = ycord(p); nv++; }
	    free(bnd);
	}
    }

    /* Rebuild each repaired loop as a region: constrained Delaunay
     * triangulation of the loop (with T-junction vertices of neighbouring
     * loops and untouched regions inserted) plus the old interior points of
     * the same material that are still inside and at least delta from it. */
    pidx = salloc(int, 4*total + 2*np + 16);
    {
    int *extra = NULL, nextra = 0, aextra = 0;
    for (i = 0; i < nl; i++) {
	int m = omat[i], cnt = 0, rn, nin = 0, ntri, *tris = NULL, *ids2;
	double *cxy;
	for (j = 0; j < llen[i]; j++) {
	    double ax = oxy[2*(lstart[i]+j)], ay = oxy[2*(lstart[i]+j)+1];
	    double bx = oxy[2*(lstart[i]+(j+1)%llen[i])], by = oxy[2*(lstart[i]+(j+1)%llen[i])+1];
	    double us[256]; int ids[256], nu = 0, a, b;
	    int p = find_point(ax, ay);
	    if (p < 0) { float c[MAXDIM]; c[0] = ax; c[1] = ay; p = mk_pt(2, c); }
	    if (cnt == 0 || pidx[cnt-1] != p) pidx[cnt++] = p;
	    for (k = 0; k < nv && nu < 256; k++) {
		double u;
		if (on_segment(vx[k], vy[k], ax, ay, bx, by, &u)) { us[nu] = u; ids[nu] = k; nu++; ninsert++; }
	    }
	    for (a = 0; a < nu; a++) for (b = a+1; b < nu; b++) if (us[b] < us[a]) {
		double tu = us[a]; int ti = ids[a]; us[a] = us[b]; ids[a] = ids[b]; us[b] = tu; ids[b] = ti; }
	    for (a = 0; a < nu; a++) {
		int q = find_point(vx[ids[a]], vy[ids[a]]);
		if (q < 0) { float c[MAXDIM]; c[0] = vx[ids[a]]; c[1] = vy[ids[a]]; q = mk_pt(2, c); }
		if (pidx[cnt-1] != q) pidx[cnt++] = q;
	    }
	}
	while (cnt > 1 && pidx[cnt-1] == pidx[0]) cnt--;
	if (cnt < 3) continue;
	triangulate_loop:
	/* A loop that revisits a point (touching at a vertex) is split there
	 * into simple loops; the remainder is queued as another loop. */
	{
	    int a, b2, found = 0;
	    for (a = 0; a < cnt && !found; a++) for (b2 = a+1; b2 < cnt; b2++) if (pidx[a] == pidx[b2]) {
		int len2 = cnt - (b2 - a), q;
		/* queue pidx[b2 .. cnt-1] + pidx[0 .. a-1] as a new loop of points */
		if (nextra + len2 + 1 > aextra) { aextra = 2*(aextra + len2 + 16); extra = realloc(extra, sizeof(int)*aextra); }
		extra[nextra++] = -(len2 + 1) - m*100000;	/* header: -(len+1) - mat*100000 */
		for (q = b2; q < cnt; q++) extra[nextra++] = pidx[q];
		for (q = 0; q < a; q++) extra[nextra++] = pidx[q];
		/* keep pidx[a .. b2-1] */
		for (q = a; q < b2; q++) pidx[q - a] = pidx[q];
		cnt = b2 - a; found = 1; break;
	    }
	    if (found && getenv("SUPREM4_MESH_DIAG")) fprintf(stderr, "[oxide-repair] split a pinched loop at (%.7f,%.7f)\n", 1e4*xcord(pidx[0]), 1e4*ycord(pidx[0]));
	    if (found) {
		if (cnt >= 3) goto triangulate_loop;		/* re-check the kept part */
		goto next_queued;
	    }
	}
	/* interior candidates */
	cxy = salloc(double, 2*(cnt + nint) + 2);
	ids2 = salloc(int, cnt + nint + 1);
	for (j = 0; j < cnt; j++) { cxy[2*j] = xcord(pidx[j]); cxy[2*j+1] = ycord(pidx[j]); }
	for (k = 0; k < nint; k++) {
	    double x = ix[k], y = iy[k], dmin = HUGE_VAL, hloc = HUGE_VAL;
	    int inside = 0, a;
	    if (imat[k] != m) continue;
	    for (a = 0; a < cnt; a++) {
		double x1 = cxy[2*a], y1 = cxy[2*a+1], x2 = cxy[2*((a+1)%cnt)], y2 = cxy[2*((a+1)%cnt)+1];
		double dx = x2-x1, dy = y2-y1, l2 = dx*dx + dy*dy, u, d;
		if ((y1 > y) != (y2 > y) && x < x1 + (y-y1)*(x2-x1)/(y2-y1)) inside = !inside;
		u = l2 > 0 ? ((x-x1)*dx + (y-y1)*dy)/l2 : 0; u = u < 0 ? 0 : (u > 1 ? 1 : u);
		d = hypot(x - x1 - u*dx, y - y1 - u*dy);
		if (d < dmin) { dmin = d; hloc = sqrt(l2); }
	    }
	    if (!inside || dmin < delta || dmin < 0.35*hloc) continue;
	    cxy[2*(cnt+nin)] = x; cxy[2*(cnt+nin)+1] = y; ids2[nin] = k; nin++;
	}
	ntri = suprem_cdt(cnt + nin, cnt, cxy, &tris);
	if (ntri <= 0) {
	    fprintf(stderr, "[oxide-repair] triangulation of a %s loop (%d points) failed\n",
		    m == Poly ? "poly" : "oxide", cnt);
	    (void)ig2_write("oxide_repair_failure.str", 0, 1.0);
	    panic("oxide repair could not triangulate a repaired region");
	}
	{
	    int *node = salloc(int, cnt + nin + 1);
	    for (j = 0; j < cnt; j++) { node[j] = pt_node(pidx[j], m); if (node[j] < 0) node[j] = mk_nd(pidx[j], m); }
	    for (j = 0; j < nin; j++) {
		float c[MAXDIM]; int p;
		c[0] = cxy[2*(cnt+j)]; c[1] = cxy[2*(cnt+j)+1];
		p = mk_pt(2, c); node[cnt+j] = mk_nd(p, m);
	    }
	    rn = mk_reg(m);
	    for (t = 0; t < ntri; t++) {
		int v[3], b[3];
		v[0] = node[tris[3*t]]; v[1] = node[tris[3*t+1]]; v[2] = node[tris[3*t+2]];
		if (orient(xcord(pt_nd(v[0])), ycord(pt_nd(v[0])), xcord(pt_nd(v[1])), ycord(pt_nd(v[1])),
			   xcord(pt_nd(v[2])), ycord(pt_nd(v[2]))) < 0) { int tmp = v[1]; v[1] = v[2]; v[2] = tmp; }
		b[0] = b[1] = b[2] = BC_OFFSET;
		(void)mk_ele(3, v, 3, b, rn, FALSE);
	    }
	    free(node);
	}
	free(tris); free(cxy); free(ids2);
	next_queued:
	/* loops split off pinched loops */
	while (nextra > 0) {
	    int h, len2, q, start = -1;
	    for (q = 0; q < nextra; q++) if (extra[q] < 0) start = q;	/* last header */
	    h = extra[start]; m = (-h) / 100000; len2 = (-h) % 100000 - 1;
	    for (q = 0; q < len2; q++) pidx[q] = extra[start + 1 + q];
	    nextra = start; cnt = len2;
	    if (cnt >= 3) goto triangulate_loop;
	}
    }
    free(extra);
    }
    /* points created for vertices that ended in no triangle */
    for (i = 0; i < np; i++) if (!dead_pt(i) && num_nd(i) == 0) fing_pt(i);
    neigh_dirty = TRUE;
    bd_connect("oxide repair: rebuilt oxide/poly");

    /* Exposure: an active boundary face without a neighbour that is not on a
     * side wall faces the ambient.  Points touch the ambient iff they lie on
     * such a face. */
    {
	int *exp = scalloc(int, np + 1), p, q;
	for (t = 0; t < ne; t++) {
	    if (dead_tri(t) || !ACTIVE_MAT(mat_reg(reg_tri(t)))) continue;
	    for (j = 0; j < 3; j++) {
		int a = pt_nd(vert_tri(t,(j+1)%3)), b = pt_nd(vert_tri(t,(j+2)%3));
		int side = (fabs(xcord(a)-xmin) < SNAP && fabs(xcord(b)-xmin) < SNAP) ||
			   (fabs(xcord(a)-xmax) < SNAP && fabs(xcord(b)-xmax) < SNAP);
		if (neigh_fc(t,j) >= 0 || side) { if (neigh_fc(t,j) < 0) tri[t]->nb[j] = BC_OFFSET; continue; }
		tri[t]->nb[j] = EXPOSED;
		exp[a] = exp[b] = 1;
	    }
	}
	for (p = 0; p < np; p++) {
	    int act = 0;
	    if (dead_pt(p)) continue;
	    for (q = 0; q < num_nd(p); q++) if (ACTIVE_MAT(mat_nd(nd_pt(p,q)))) act = 1;
	    if (!act) continue;
	    if (exp[p]) {
		set_surf(p);
		if (pt_node(p, GAS) < 0) (void)mk_nd(p, GAS);
	    }
	    else {
		int g = pt_node(p, GAS);
		clr_surf(p);
		if (g >= 0) rem_1nd(g);
	    }
	}
	free(exp);
    }
    bd_connect_material("oxide repair: exposure");

    /* Fields: interpolate, then restore each material's dopant inventory. */
    {
	int n2;
	for (n2 = 0; n2 < nn; n2++) if (!dead_nd(n2) && mat_nd(n2) == GAS)
	    for (j = 0; j < MAXIMP; j++) nd[n2]->sol[j] = gas_template[j];
	double *val = salloc(double, S.ns + 1);
	int n;
	for (n = 0; n < nn; n++) {
	    int m;
	    if (dead_nd(n)) continue;
	    m = mat_nd(n);
	    if (!ACTIVE_MAT(m)) continue;
	    snap_eval(&S, m, xcord(pt_nd(n)), ycord(pt_nd(n)), val);
	    for (i = 0; i < S.ns; i++) nd[n]->sol[S.slot[i]] = val[i];
	}
	free(val);
	for (s = 0; s < n_imp; s++) {
	    int mi;
	    if (!conserved(soltoimp[s])) continue;
	    for (mi = 0; mi < 2; mi++) {
		int m = mi ? Poly : SiO2;
		double now = inventory_mesh(m, s), f;
		if (now <= 0 || inv_old[s][mi] <= 0) continue;
		f = inv_old[s][mi] / now;
		for (n = 0; n < nn; n++) if (!dead_nd(n) && mat_nd(n) == m) nd[n]->sol[s] *= f;
	    }
	}
    }

    /* Validate. */
    {
	int bad = -1;
	double q, worst = 1;
	for (t = 0; t < ne; t++) if (!dead_tri(t)) {
	    q = tri_quality(t);
	    if (q <= 0) { bad = t; break; }
	    if (ACTIVE_MAT(mat_reg(reg_tri(t))) && q < worst) worst = q;
	}
	if (bad >= 0 || !mesh_trial_validate_embedding(0.0, &bad)) {
	    fprintf(stderr, "[oxide-repair] invalid rebuilt mesh at triangle %d", bad);
	    if (bad >= 0 && bad < ne) {
		fprintf(stderr, " mat=%d q=%.3g", mat_reg(reg_tri(bad)), tri_quality(bad));
		for (j = 0; j < 3; j++) fprintf(stderr, " (%.7f,%.7f)", 1e4*xcord(pt_nd(vert_tri(bad,j))), 1e4*ycord(pt_nd(vert_tri(bad,j))));
	    }
	    fprintf(stderr, "\n[oxide-repair] %d T-junction insertions, %d candidate vertices\n", ninsert, nv);
	    if (bad >= 0 && bad < ne) for (j = 0; j < 3; j++) {
		int a = pt_nd(vert_tri(bad,j)), b = pt_nd(vert_tri(bad,(j+1)%3)), kk, best = -1; double bd = HUGE_VAL, u;
		for (kk = 0; kk < nv; kk++) {
		    double mx = 0.5*(xcord(a)+xcord(b)), my = 0.5*(ycord(a)+ycord(b)), dd = hypot(vx[kk]-mx, vy[kk]-my);
		    if (dd < bd && hypot(vx[kk]-xcord(a), vy[kk]-ycord(a)) > SNAP && hypot(vx[kk]-xcord(b), vy[kk]-ycord(b)) > SNAP) { bd = dd; best = kk; }
		}
		if (best >= 0) fprintf(stderr, "  edge %d: nearest candidate (%.9f,%.9f) on_segment=%d\n", j, 1e4*vx[best], 1e4*vy[best],
			on_segment(vx[best], vy[best], xcord(a), ycord(a), xcord(b), ycord(b), &u));
	    }
	    (void)ig2_write("oxide_repair_failure.str", 0, 1.0);
	    panic("oxide repair produced an invalid mesh");
	}
	if (getenv("SUPREM4_MESH_DIAG")) {
	    int wt = -1, nexp = 0, nsurf = 0, p2;
	    double wq = 2;
	    for (t = 0; t < ne; t++) if (!dead_tri(t)) {
		for (j = 0; j < 3; j++) if (tri[t]->nb[j] == EXPOSED) nexp++;
		if (ACTIVE_MAT(mat_reg(reg_tri(t))) && (q = tri_quality(t)) < wq) { wq = q; wt = t; }
	    }
	    for (p2 = 0; p2 < np; p2++) if (!dead_pt(p2) && is_surf(p2)) nsurf++;
	    fprintf(stderr, "[oxide-repair] %s: %d loops rebuilt, %d T-junctions, min quality %.3g -> %.3g; %d exposed faces, %d surface points\n",
		    changed ? "geometry repaired" : "remeshed", nl, ninsert, minq, worst, nexp, nsurf);
	    if (wt >= 0) {
		fprintf(stderr, "[oxide-repair] worst mat=%d", mat_reg(reg_tri(wt)));
		for (j = 0; j < 3; j++) { int pp = pt_nd(vert_tri(wt,j));
		    fprintf(stderr, " (%.7f,%.7f nn=%d)", 1e4*xcord(pp), 1e4*ycord(pp), num_nd(pp)); }
		fprintf(stderr, "\n");
	    }
	}
    }
    UpdateSymbolic++;

    free(pidx); free(vx); free(vy); free(lstart); free(llen); free(ix); free(iy); free(imat); free(gas_template);
    free(omat); free(ocount); free(oxy);
    snap_free(&S); loops_free(&L);
    return 1;
#endif
}

/* Largest point speed (cm/s), for the delta/2 per-step motion limit. */
double oxide_repair_max_speed(void)
{
    int p; double v = 0;
    for (p = 0; p < np; p++) if (!dead_pt(p)) {
	double s = hypot(pt[p]->vel[0], pt[p]->vel[1]);
	if (s > v) v = s;
    }
    return v;
}
