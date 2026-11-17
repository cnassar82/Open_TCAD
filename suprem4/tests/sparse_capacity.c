/* Exercise real sparse-map construction with IDs above the signed-short limit
 * and a node whose degree exceeds the initial neighbor scratch capacity. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#define STATIC_ALLOCATION_TIME
#include "global.h"
#include "constant.h"
#include "geom.h"

#define GRAPH_NODES 40001

void panic(char *message)
{
    fprintf(stderr, "Unexpected panic: %s\n", message);
    abort();
}

/* A symmetric star including its diagonal. These supply mesh adjacency;
 * get_connect(), its capacity contract, and generate() are production code. */
int num_neigh_nd(int node)
{
    return node == 0 ? GRAPH_NODES : 2;
}

int neigh_nd(int node, int index)
{
    return node == 0 ? index : (index == 0 ? 0 : node);
}

int main(void)
{
    int sentinel[3] = {123, 456, 789}, available = 1;
    int *map, capacity, row, start, entries, pass;
    assert(get_connect(0, sentinel + 1, &available) == -1);
    assert(available == GRAPH_NODES);
    assert(sentinel[0] == 123 && sentinel[1] == 456 && sentinel[2] == 789);

    for (pass = 0; pass < 3; pass++) {
        /* Deliberately smaller than the row-pointer header. */
        capacity = 1;
        map = salloc(int, capacity);
        assert(generate(GRAPH_NODES, 1, 1, get_connect, &map, &capacity) == 0);
        assert(map[0] == GRAPH_NODES + 1);
        assert(capacity == 3 * GRAPH_NODES - 2);
        start = map[0];
        for (row = 0; row < GRAPH_NODES; row++) assert(map[start + row] == row);
        for (row = 1; row < GRAPH_NODES; row++) {
            start = map[row];
            entries = map[row+1] - start;
            assert(entries == 2 && map[start] == 0 && map[start+1] == row);
        }
        assert(map[map[GRAPH_NODES] - 1] == GRAPH_NODES - 1);
        free(map);
    }
    available = 0;
    clean_list(sentinel, &available);
    assert(available == 0);
    puts("PASS: 40001 node IDs, high-degree capacity, map growth, and repeated cleanup");
    return 0;
}
