/*
 * mem_util.c -- DCR 94-10-12
 * ==========================
 * Memory allocation and deallocation routines.
 *
 * Global functions: Elem(), Array(), Grid(), FreeElem(), FreeArray(),
 *   FreeGrid().
 *
 */

/*** Include files ***/

#include <stdio.h>
#include <stdlib.h>
/*#include <malloc.h>*/
#include <mem_util.h>

/*** End of preamble ***/

void *Elem(size_t size)
{
  /* Allocates and returns pointer to one element of memory, size "size" */

  void *ptr;

  if ((ptr = malloc(size)) == NULL) {
    (void) fprintf(stderr, "Elem(): Unable to allocate memory "
		   "(%lu byte(s) requested).\n", (unsigned long) size);
    exit(1);
  }

  return ptr;
}

void *Array(size_t n, size_t size)
{
  /* Allocates array: "n" elements of size "size" */

  void *ptr;

  if ((ptr = malloc(n * size)) == NULL) {
    (void) fprintf(stderr, "Array(): Unable to allocate array "
		   "(%lu element(s) of %lu byte(s) requested).\n",
		   (unsigned long) n, (unsigned long) size);
    exit(1);
  }

  return ptr;
}

void **Grid(size_t m, size_t n, size_t size)
{
  /* Allocates "m" by "n" grid, each element size "size" */

  size_t i;
  void **ptr;

  ptr = (void **) Array(m, sizeof(void *));

  for (i = 0; i < m; i++)
    ptr[i] = Array(n, size);

  return ptr;
}

void FreeElem(void *ptr)
{
  /* Frees memory element pointed to by "ptr" */

  if (ptr == NULL) {
    (void) fprintf(stderr, "FreeElem(): Invalid memory address.\n");
    exit(1);
  }

  free(ptr);
}

void FreeArray(void *ptr)
{
  /* Frees array at "ptr" */

  FreeElem(ptr);
}

void FreeGrid(size_t m, void **ptr)
{
  /* Frees grid ("m" rows) at "ptr" */

  size_t i;

  for (i = 0; i < m; i++)
    FreeArray(ptr[i]);

  FreeArray((void *) ptr);
}

/* mem_util.c */
