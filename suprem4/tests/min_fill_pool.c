/* Include the implementation to force exhaustion of its private link pool. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "../src/math/min_fill.c"

void panic(char *message)
{
    fprintf(stderr, "Unexpected panic: %s\n", message);
    abort();
}

int node_mat(int node, int material)
{
    assert(!"1D path is not used in this test");
    return -1;
}

int main(void)
{
    link *head = NULL, *p;
    int i, pass, cursor, row, neighbor;
    int ia[31], reorder[8], seen[8];

    storage_size = 3;
    su = hwm = 0;
    for (i = 0; i < 100; i++) {
        new_link(head);
        head->v = i;
        /* Every growth must preserve every previously handed-out pointer. */
        for (p = head, row = i; p; p = p->n, row--) assert(p->v == row);
    }
    assert(su == 100 && hwm == 100);
    while (head) {
        p = head;
        head = head->n;
        p->n = freelist;
        freelist = p;
        su--;
    }
    assert(su == 0);
    /* Recycling links must also preserve links into older blocks. */
    for (i = 0; i < 100; i++) {
        new_link(head);
        head->v = i;
    }
    for (p = head, i = 99; p; p = p->n, i--) assert(p->v == i);
    free_storage();
    assert(!storage && !freelist && su == 0);

    /* A symmetric chain, including each diagonal: known no-fill ordering. */
    mode = TWOD;
    for (pass = 0; pass < 3; pass++) {
        cursor = 9;
        for (row = 0; row < 8; row++) {
            ia[row] = cursor;
            for (neighbor = max(0, row-1); neighbor <= min(7, row+1); neighbor++)
                ia[cursor++] = neighbor;
            seen[row] = 0;
        }
        ia[8] = cursor;
        assert(cursor == 31);
        assert(min_ia_fill(ia, reorder, 0) == 7);
        for (i = 0; i < 8; i++) {
            assert(reorder[i] >= 0 && reorder[i] < 8);
            assert(!seen[reorder[i]]++);
        }
        assert(!storage && !freelist && su == 0);
    }
    puts("PASS: stable ordering-pool growth, recycling, and repeated no-fill ordering");
    return 0;
}
