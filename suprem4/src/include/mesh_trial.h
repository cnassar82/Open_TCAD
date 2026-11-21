#ifndef SUPREM_MESH_TRIAL_H
#define SUPREM_MESH_TRIAL_H

/* A diffusion trial may move points and flip edges, but must not add,
 * remove or renumber mesh entities between save and restore. */
struct mesh_trial;
struct mesh_trial *mesh_trial_save(void);
void mesh_trial_restore(const struct mesh_trial *trial);
void mesh_trial_discard(struct mesh_trial *trial);
/* Implicit diffusion carries inventory in its old/mid/new mass terms.
 * Static rezoning must transfer nd[].sol; an active implicit trial must not. */
int mesh_trial_active(void);

/* Read-only embedding check; dt=0 checks committed coordinates. */
int mesh_trial_validate_embedding(double dt, int *bad_triangle);

/* Returns zero without mutation if rounded motion is invalid. */
int mesh_trial_move(double dt, int *bad_triangle);

#endif
