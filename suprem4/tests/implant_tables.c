/* Exercise actual table lookup and ion/material mapping, not fitted substitutes. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "global.h"
#include "constant.h"
#include "material.h"
#include "impurity.h"
#include "implant.h"

char *IMP_FILE;
int imp_model;
int fgtmom(), setprs();
/* Mapping tests use Gaussian constants; integration is tested by process decks. */
int qpgauss(double dx, struct imp_info *info) { return 0; }
int gtprs() { assert(0); return -1; }
int qpprs() { assert(0); return -1; }

static double range(int ion, int mat, double energy, double *lateral)
{
    struct imp_info info;
    assert(setprs(ion, mat, energy, &info) == 0);
    assert(isfinite(info.vert[GUS_RP]));
    assert(isfinite(info.lat[LATSIG]) && info.lat[LATSIG] > 0);
    *lateral = info.lat[LATSIG];
    return info.vert[GUS_RP];
}

int main(int argc, char **argv)
{
    double rp, lat, sigma, gamma, kurt, expected;
    assert(argc == 2);
    IMP_FILE = argv[1];
    override = 0; imp_model = GAUSS;
    /* This exact first row used to multiply uninitialized lower moments by 0. */
    assert(fgtmom(B, SiO2, 10., &rp, &sigma, &gamma, &kurt, &lat) == 0);
    assert(fabs(rp - .0355) < 1e-12);
    assert(fabs(range(B, SiO2, 60., &lat) - .2070) < 1e-12);
    assert(fgtmom(P, SiO2, 125., &rp, &sigma, &gamma, &kurt, &lat) == 0);
    assert(fabs(range(P, SiO2, 125., &lat) - rp) < 1e-12);
    expected = .07555 + (.08207 - .07555) / 2.;
    assert(fabs(range(BF2, Si, 110., &lat) - expected) < 1e-12);
    /* Seven-column legacy rows use documented isotropic lateral straggle. */
    assert(fabs(lat - (.03702 + .03490) / 2.) < 1e-12);
    expected = .0713 + (.1068 - .0713) * (110. - 89.1) / (134. - 89.1);
    assert(fabs(range(BF2, SiO2, 110., &lat) - expected) < 1e-12);
    assert(fgtmom(999, Si, 125., &rp, &sigma, &gamma, &kurt, &lat) < 0);
    puts("PASS: oxide/P125/BF2 mapping, interpolation, first row, missing table");
    return 0;
}
