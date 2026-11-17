/* Exercise byte counts without actually allocating multiple gigabytes. */
#include <assert.h>
#include <limits.h>
#include <setjmp.h>
#include <stdio.h>
#include <string.h>
#include "global.h"

static jmp_buf overflow;
static size_t requested;
static int calls;
static double storage;

void panic(char *message)
{
    assert(strcmp(message, "Allocation size overflow") == 0);
    longjmp(overflow, 1);
}

static void *record_malloc(size_t bytes)
{
    requested = bytes;
    calls++;
    return &storage;
}

static void *record_calloc(size_t count, size_t width)
{
    return record_malloc(count * width);
}

static void *record_realloc(void *pointer, size_t bytes)
{
    assert(pointer == &storage);
    return record_malloc(bytes);
}

#define malloc record_malloc
#define calloc record_calloc
#define realloc record_realloc

int main(void)
{
    size_t count = 3;
    double *pointer;
    pointer = salloc(double, count++);
    assert(pointer == &storage && count == 4 && requested == 3 * sizeof(double));
    pointer = scalloc(double, count++);
    assert(pointer == &storage && count == 5 && requested == 4 * sizeof(double));
    pointer = sralloc(double, count++, pointer);
    assert(pointer == &storage && count == 6 && requested == 5 * sizeof(double));

    if (sizeof(size_t) > sizeof(unsigned int)) {
        count = (size_t)UINT_MAX / sizeof(double) + 2;
        pointer = salloc(double, count);
        assert(pointer == &storage && requested > UINT_MAX);
        pointer = scalloc(double, count);
        assert(pointer == &storage && requested == count * sizeof(double));
        pointer = sralloc(double, count, pointer);
        assert(pointer == &storage && requested == count * sizeof(double));
    }
    calls = 0;
    if (!setjmp(overflow)) {
        pointer = salloc(double, (size_t)-1 / sizeof(double) + 1);
        assert(!"overflow was not rejected");
    }
    assert(calls == 0);
    if (!setjmp(overflow)) {
        pointer = scalloc(double, (size_t)-1 / sizeof(double) + 1);
        assert(!"calloc overflow was not rejected");
    }
    assert(calls == 0);
    if (!setjmp(overflow)) {
        pointer = sralloc(double, (size_t)-1 / sizeof(double) + 1, &storage);
        assert(!"realloc overflow was not rejected");
    }
    assert(calls == 0);
    puts("PASS: allocation widths, single evaluation, and overflow rejection");
    return 0;
}
