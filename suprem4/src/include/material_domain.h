#ifndef SUPREM_MATERIAL_DOMAIN_H
#define SUPREM_MATERIAL_DOMAIN_H
/* Separate physical equivalence from legacy region/deposition identity.
 * protected_region or differing model keys preserve a physical constraint.
 * This bridge only produces disk regions, as required by the legacy skeleton.
 * Geometry and fields are inputs to neither relabeling nor retirement. */
struct domain_region {int material, model_key, protected_region;};
struct domain_triangle {int vertex[3], region;};
/* 0: validated map; -1: invalid incidence; -2: allocation failure.
 * On failure map is unchanged. On success entries are compact domain IDs. */
int material_domain_plan(int np,int nr,const struct domain_region *regions,
    int nt,const struct domain_triangle *triangles,int *map,int *domains);
/* Opt-in legacy adapter, called only while derived connectivity is rebuilding. */
int material_domain_normalize(void);
#endif
