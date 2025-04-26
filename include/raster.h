/*
 * raster.h -- DCR 94-08-30
 * ========================
 * Type definitions and function prototypes for raster handling routines.
 *
 */

/* Copyright notice */

#include <copyright.h>

/* Other header files to include */

#include <stdio.h>

#ifndef sparc
#ifdef SunOS
#  include <sys/stdtypes.h> /* For size_t */
#endif
#endif

#include <colors.h>

/* Colormap structure */

typedef struct {
  int type;         /* Colormap type */
  size_t size;      /* Number of colors in colormap */
  COLOR *r, *g, *b; /* Pointers to red, green, and blue colormap entries */
} COLORMAP_T;

/* Image structure */

/*
 * NOTE: w & h should be integers to ensure comparisons with negative
 * numbers do not cause unexpected problems (i.e. promotion of negative
 * values to unsigned...)
 *
 */

typedef struct {
  int w;        /* Width in pixels */
  int h;        /* Height in pixels */
  COLOR **data; /* Pointer to image data */
} IMAGE_T;

/* Function prototypes (cf. raster.c) */

void MakeColormap(COLORMAP_T *colormap);
void AllocImage(IMAGE_T **image, int w, int h);
void DrawPoint(IMAGE_T *image, int x, int y, COLOR c);
void DrawVector(IMAGE_T *image, int x1, int y1, int x2, int y2, COLOR c);
void DumpRaster(COLORMAP_T *colormap, IMAGE_T *image, const char *filename);
void FreeImage(IMAGE_T *image);

/* raster.h */
