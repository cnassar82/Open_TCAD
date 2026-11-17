#ifndef SUPREM_PROCESS_SURFACE_H
#define SUPREM_PROCESS_SURFACE_H
#include <stdio.h>
struct line;
/* A construction-precision description of the last deposited free surface.
 * It is valid only while the complete exposed material boundary still matches
 * the recorded mesh embedding. It is not a replacement for buried interfaces. */
void process_surface_clear(void);
void process_surface_record(int material,int count,const double *xy);
int process_surface_get(int material,int actual_count,const struct line *actual,
                        const double **xy);
void process_surface_write(FILE *file,double sx,double sy);
int process_surface_read(const char *line,double sx,double sy);
#endif
