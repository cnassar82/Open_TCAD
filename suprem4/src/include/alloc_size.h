#ifndef SUPREM_ALLOC_SIZE_H
#define SUPREM_ALLOC_SIZE_H

#include <stdlib.h>

/* Retain size_t width and reject wrapped byte counts before allocation. */
static __inline__ size_t suprem_array_bytes(size_t count, size_t width)
{
    if (width != 0 && count > (size_t)-1 / width) {
        panic("Allocation size overflow");
        abort();
    }
    return count * width;
}

#endif
