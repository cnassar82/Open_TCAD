#ifndef SUPREM_MESH_PATCH_H
#define SUPREM_MESH_PATCH_H

/* Pure geometry proposal: input is a complete live triangular star.
 * Coordinates are the stored binary32 values; vertices are local indices.
 * No database objects are created, retired or modified by this API. */
struct mesh_patch_point { float xy[2]; };
struct mesh_patch_triangle { int vertex[3]; int region; };
struct mesh_patch_result {
    int count;
    struct mesh_patch_triangle *triangle;
    double old_min_quality, new_min_quality;
    const char *reason;
};
enum mesh_patch_status {
    MESH_PATCH_OK, MESH_PATCH_PINNED, MESH_PATCH_INVALID, MESH_PATCH_NOMEM
};
int mesh_patch_propose(int center, int point_count,
    const struct mesh_patch_point *point, int triangle_count,
    const struct mesh_patch_triangle *triangle, struct mesh_patch_result *result);
int mesh_patch_remove_point(int point, const char **reason);
void mesh_patch_release(struct mesh_patch_result *result);
#endif
