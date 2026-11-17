/* Focused memory-safety tests of the production node-removal and neighbor
 * routines. Skeleton operations are mocked: this tests scratch storage and
 * boundary ordering, not the correctness of cavity retriangulation. */
#include <assert.h>
#include <math.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#define STATIC_ALLOCATION_TIME
#include "constant.h"
#include "global.h"
#include "dbaccess.h"
#include "skel.h"

extern int kill_node();
extern int tnabor();
static int valence, edge_count, grid_count, fail_realloc;
static const char *expect_panic;
static int ends[1024][2];
static jmp_buf panic_target;

void panic(message)
char *message;
{
    if (expect_panic) {
        assert(strstr(message, expect_panic) != NULL);
        longjmp(panic_target, 1);
    }
    fprintf(stderr, "Unexpected panic: %s\n", message);
    abort();
}

extern void *__real_realloc(void *, size_t);
void *__wrap_realloc(void *pointer, size_t bytes)
{
    if (fail_realloc) { fail_realloc = 0; return NULL; }
    return __real_realloc(pointer, bytes);
}

int cr_sreg(mat) int mat;
{
    static struct sreg_str skeleton;
    sreg[0] = &skeleton;
    return 0;
}
int mk_edg(a, b) int a, b;
{
    assert(edge_count < 1024);
    ends[edge_count][0] = a;
    ends[edge_count][1] = b;
    return edge_count++;
}
int ad_edge() { return 0; }
int hinit() { return 0; }
int pl_interface() { return 0; }
int free_skel() { return 0; }
int rem_ele() { nd[0]->ltri.num = 0; return 0; }
int grid()
{
    int i;
    double previous = -4.0;
    grid_count++;
    assert(edge_count == valence);
    for (i = 0; i < edge_count; i++) {
        int a = ends[i][0], b = ends[i][1];
        double angle = atan2(pt[a]->cord[1], pt[a]->cord[0]);
        assert(angle >= previous);
        previous = angle;
        assert(b == ends[(i + 1) % edge_count][0]);
    }
    return 0;
}

static void test_node_removal(count)
int count;
{
    int i;
    valence = count;
    edge_count = grid_count = 0;
    nn = np = count + 1;
    ne = count;
    nd = calloc(nn, sizeof(*nd));
    pt = calloc(np, sizeof(*pt));
    tri = calloc(ne, sizeof(*tri));
    for (i = 0; i < nn; i++) {
        double angle = count ? 2.0 * PI * (i - 1) / count : 0;
        nd[i] = calloc(1, sizeof(*nd[i]));
        pt[i] = calloc(1, sizeof(*pt[i]));
        nd[i]->pt = i;
        pt[i]->nn = 1;
        pt[i]->nd[0] = i;
        if (i) {
            pt[i]->cord[0] = cos(angle);
            pt[i]->cord[1] = sin(angle);
        }
    }
    for (i = 0; i < count; i++) {
        tri[i] = calloc(1, sizeof(*tri[i]));
        tri[i]->nd[0] = 0;
        tri[i]->nd[1] = i + 1;
        tri[i]->nd[2] = (i + 1) % count + 1;
        add_list(&nd[0]->ltri, i);
    }
    /* Exercise the real neighbor collector beyond its old 200-entry
     * buffer, then pass that connectivity to the real kill_node. */
    assert(num_neigh_nd(0) == count + 1);
    for (i = 0; i <= count; i++) assert(neigh_nd(0, i) == i);
    assert(num_neigh_pt(0) == count + 1);
    for (i = 0; i <= count; i++) assert(neigh_pt(0, i) == i);
    kill_node(0);
    assert(nd[0]->ltri.num == 0);
    assert(grid_count == (count >= 3));
    dis_list(&nd[0]->ltri);
    for (i = 0; i < nn; i++) { free(nd[i]); free(pt[i]); }
    for (i = 0; i < ne; i++) free(tri[i]);
    free(nd); free(pt); free(tri);
}

static void test_interface_neighbors()
{
    struct nd_str nodes[5];
    struct nd_str *node_ptrs[5];
    struct pt_str point;
    struct pt_str *point_ptrs[] = {&point};
    struct tri_str triangles[2];
    struct tri_str *triangle_ptrs[] = {&triangles[0], &triangles[1]};
    int i, first = 0, second = 1, empty = 0;
    memset(nodes, 0, sizeof(nodes));
    memset(&point, 0, sizeof(point));
    memset(triangles, 0, sizeof(triangles));
    nd = node_ptrs; pt = point_ptrs; tri = triangle_ptrs;
    for (i = 0; i < 5; i++) nd[i] = &nodes[i];
    nodes[2].pt = 1; nodes[3].pt = 2; nodes[4].pt = 3;
    point.nn = 2; point.nd[0] = 0; point.nd[1] = 1;
    nodes[0].ltri.num = nodes[1].ltri.num = 1;
    nodes[0].ltri.list = &first; nodes[1].ltri.list = &second;
    triangles[0].nd[0] = 0; triangles[0].nd[1] = 2; triangles[0].nd[2] = 3;
    triangles[1].nd[0] = 1; triangles[1].nd[1] = 3; triangles[1].nd[2] = 4;
    /* Both material nodes at the interface must contribute neighbors. */
    assert(num_neigh_pt(0) == 4);
    for (i = 0; i < 4; i++) assert(neigh_pt(0, i) == i);
    elim_dups(NULL, &empty);
    assert(empty == 0);
}

static void test_neighbors()
{
    struct nd_str nodes[6];
    struct nd_str *node_ptrs[6];
    struct tri_str triangles[2];
    struct tri_str *triangle_ptrs[2];
    int i, ke, kn;
    memset(nodes, 0, sizeof(nodes));
    memset(triangles, 0, sizeof(triangles));
    nd = node_ptrs;
    tri = triangle_ptrs;
    for (i = 0; i < 6; i++) { nd[i] = &nodes[i]; nd[i]->pt = i; }
    for (i = 0; i < 2; i++) tri[i] = &triangles[i];
    for (i = 0; i < 3; i++) {
        tri[0]->nd[i] = i;
        tri[1]->nd[i] = i + 3;
    }
    assert(!tnabor(0, 1, &ke, &kn));
    nodes[3].pt = 0;
    assert(!tnabor(0, 1, &ke, &kn));
    nodes[4].pt = 1;
    assert(tnabor(0, 1, &ke, &kn));
    assert(ke == 2 && kn == 2);
    /* Three matches used to overflow the two-entry local array. */
    nodes[5].pt = 2;
    expect_panic = "tnabor";
    if (!setjmp(panic_target)) { tnabor(0, 1, &ke, &kn); abort(); }
    /* Degenerate triangles produce all nine pair matches. */
    for (i = 0; i < 6; i++) nodes[i].pt = 0;
    if (!setjmp(panic_target)) { tnabor(0, 1, &ke, &kn); abort(); }
    expect_panic = 0;
    mode = ONED;
    nodes[0].pt = 0; nodes[1].pt = 1;
    nodes[3].pt = 1; nodes[4].pt = 2;
    assert(tnabor(0, 1, &ke, &kn));
    assert(ke == 0 && kn == 1);
}

static void test_large_list()
{
    struct list_str list = {0, 0, NULL};
    int i;
    for (i = 0; i < 100000; i++) add_list(&list, i);
    assert(list.num == 100000 && list.all >= list.num);
    for (i = 0; i < list.num; i++) assert(list.list[i] == i);
    sub_list(&list, 32767);
    assert(list.num == 99999 && list.list[32767] == 32768);
    add_list(&list, 100000);
    assert(list.num == 100000 && list.list[99999] == 100000);
    dis_list(&list);
    assert(list.num == 0 && list.all == 0 && list.list == NULL);
    add_list(&list, 42);
    assert(list.num == 1 && list.list[0] == 42);
    dis_list(&list);
}

static void test_list_failure()
{
    /* Keep this static: its contents must be inspectable after longjmp. */
    static struct list_str list;
    int i, capacity, *original;
    for (i = 0; i < 10; i++) add_list(&list, i);
    original = list.list;
    capacity = list.all;
    fail_realloc = 1;
    expect_panic = "Out of storage";
    if (!setjmp(panic_target)) { add_list(&list, 10); abort(); }
    assert(list.list == original && list.all == capacity && list.num == 10);
    for (i = 0; i < 10; i++) assert(list.list[i] == i);
    expect_panic = NULL;
    add_list(&list, 10);
    assert(list.list[10] == 10);
    dis_list(&list);
    list.num = list.all = INT_MAX;
    expect_panic = "integer capacity";
    if (!setjmp(panic_target)) { add_list(&list, 10); abort(); }
    expect_panic = NULL;
    list.num = list.all = 0;
}

int main()
{
    int i;
    int counts[] = {0, 2, 3, 50, 51, 100, 101, 512};
    mode = TWOD;
    nvrt = 3;
    for (i = 0; i < sizeof(counts) / sizeof(counts[0]); i++)
        test_node_removal(counts[i]);
    test_interface_neighbors();
    test_neighbors();
    test_large_list();
    test_list_failure();
    puts("PASS: large connectivity lists, high-valence nodes, and triangle-neighbor safety");
    return 0;
}
