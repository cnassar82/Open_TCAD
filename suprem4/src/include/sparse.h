#ifndef SUPREM_SPARSE_H
#define SUPREM_SPARSE_H

/* Neighbor identifiers use the same width as the mesh database. */
extern int generate(int, int, int, int (*)(int, int *, int *), int **, int *);
extern int get_connect(int, int *, int *);
extern int FEconnect(int, int *, int *);
extern int get_elasconn(int, int *, int *);
extern int clean_list(int *, int *);

#endif
