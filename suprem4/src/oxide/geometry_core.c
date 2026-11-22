/*
 * Geometry core hook (M0 WP6, suprem4/doc/plan_m0_geometry_core.md).
 *
 * Active only when the environment variable SUPREM4_GEOMETRY_CORE is set
 * (and not "0") AND the executable was linked with the core library
 * (build_native64.py --geometry-core); otherwise this is a no-op, so builds
 * and runs without it are unchanged.
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
 * Other settings: SUPREM4_CORE_TMP (directory for the exchange files,
 * default "."), SUPREM4_CORE_FORCE_FIRST=1 (regenerate on the first call),
 * SUPREM4_CORE_HMAX, SUPREM4_CORE_SCALE (mesh size, read by the core).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
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
    calls++;
    if (ig2_write(in, 0, 1.0) < 0) {
	fprintf(stderr, "[geometry-core] cannot write %s\n", in);
	return 0;
    }
    msg[0] = 0;
    rc = suprem_core_oxidation_check(in, out, env_double("SUPREM4_CORE_DELTA", 1e-3),
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
    return 1;
}
