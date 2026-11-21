/* Runtime algorithm options and their validated defaults.
 *
 * These algorithms were introduced as opt-in SUPREM4_* environment switches
 * and are now the default: they complete the fresh RIT 0.5 um CMOS NMOS/PMOS
 * flows through metal1 with every checkpoint audit passing (see
 * doc/rit_cmos_full_flow_repair_20261119.md) and pass the geometry/etch
 * regression battery at least as well as the legacy paths.
 *
 * suprem_option(name) returns the effective value of an option:
 *   - the environment value, when set to anything but a disabling word;
 *   - NULL when set to "0", "off", "no", "false", "legacy" or "" (the legacy
 *     algorithm is then used; presence-tested switches are disabled too);
 *   - otherwise the validated default below (NULL if the option has none).
 * Diagnostic and trace switches are not listed and keep plain getenv().
 */
#ifndef SUPREM_OPTIONS_H
#define SUPREM_OPTIONS_H
#include <stdlib.h>
#include <string.h>

static const char *suprem_option(const char *name)
{
    static const char *const defaults[][2]={
        {"SUPREM4_PROCESS_SURFACE","1"},
#ifdef SUPREM_DOMAIN_OFFSET
        /* Requires the Clipper2 deposition backend (--domain-offset). */
        {"SUPREM4_DEPOSITION_GEOMETRY","domain"},
#endif
        {"SUPREM4_MATERIAL_DOMAINS","oxide"},
        {"SUPREM4_MESH_PATCH_REMOVAL","1"},
        {"SUPREM4_FE_OXIDANT_P1","1"},
        {"SUPREM4_FE_BBAR_P1","1"},
        {"SUPREM4_FE_LINEAR_TRACE","1"},
        {"SUPREM4_NATIVE_JOIN_CAP","1"},
        {"SUPREM4_CHARGED_FLUX","sg"},
        /* Minimum element layers across dielectric films before each
         * non-oxidizing diffusion (dbase/thin_film.c). */
        {"SUPREM4_THIN_FILM_LAYERS","8"},
    };
    static const char *const disabling[]={"0","off","no","false","legacy",""};
    const char *value=getenv(name);
    size_t i;
    if(value) {
        for(i=0;i<sizeof(disabling)/sizeof(disabling[0]);i++)
            if(!strcmp(value,disabling[i]))return NULL;
        return value;
    }
    for(i=0;i<sizeof(defaults)/sizeof(defaults[0]);i++)
        if(!strcmp(name,defaults[i][0]))return defaults[i][1];
    return NULL;
}
#endif
