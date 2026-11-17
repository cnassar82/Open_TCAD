/* Time a production ordering on a captured native-int SUPREM sparse map.
 * A matrix map starts with its row pointers, followed by neighbor IDs. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"

extern int min_ia_fill(int *, int *, int);
void panic(char *message) { fprintf(stderr, "%s\n", message); abort(); }
int node_mat(int node, int material) { abort(); return -1; }

int main(int argc, char **argv)
{
    FILE *input;
    long bytes;
    int *ia, *permutation, *seen, n, i, fill;
    clock_t started;
    if (argc != 2) { fprintf(stderr, "usage: %s captured-map.bin\n", argv[0]); return 2; }
    input = fopen(argv[1], "rb");
    assert(input && fseek(input, 0, SEEK_END) == 0);
    bytes = ftell(input);
    assert(bytes > (long)sizeof(int) && bytes % sizeof(int) == 0);
    rewind(input);
    ia = malloc(bytes);
    assert(ia && fread(ia, 1, bytes, input) == bytes);
    fclose(input);
    n = ia[0] - 1;
    assert(n > 0 && n < bytes / sizeof(int) && ia[n] == bytes / sizeof(int));
    for (i = 0; i < n; i++) assert(ia[i] <= ia[i+1]);
    for (i = ia[0]; i < ia[n]; i++) assert(ia[i] >= 0 && ia[i] < n);
    permutation = salloc(int, n);
    seen = scalloc(int, n);
    assert(permutation && seen);
    mode = TWOD;
    started = clock();
    fill = min_ia_fill(ia, permutation, 0);
    for (i = 0; i < n; i++) {
        assert(permutation[i] >= 0 && permutation[i] < n);
        assert(!seen[permutation[i]]++);
    }
    printf("{\"nodes\":%d,\"adjacency_entries\":%d,\"fill_estimate\":%d,"
           "\"cpu_seconds\":%.6f}\n", n, ia[n]-ia[0], fill,
           (double)(clock()-started) / CLOCKS_PER_SEC);
    free(ia); free(permutation); free(seen);
    return 0;
}
