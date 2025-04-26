/*
 * mem_util.h -- DCR 94-10-12
 * ==========================
 * Prototypes for memory allocation and deallocation routines.
 *
 */

#include <copyright.h>

void *Elem(size_t size);
void *Array(size_t n, size_t size);
void **Grid(size_t m, size_t n, size_t size);
void FreeElem(void *ptr);
void FreeArray(void *ptr);
void FreeGrid(size_t m, void **ptr);

/* mem_util.h */
