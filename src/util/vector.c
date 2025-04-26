/*
 * vector.c -- DCR 94-08-24
 * ========================
 * Some routines for common vector/matrix operations.
 *
 * Global functions: Transpose(), Transform(), MultMat(), GetBasis(),
 *   MakeBasis().
 *
 */

/*** Include files ***/

#include <stdio.h>
#include "vector.h"

/* Local functions */

static void swap(double *, double *);

/* End of preamble */

void Transpose(MATRIX a)
{
  /* Applies transpose operator on matrix "a" */

  swap(&a[X][Y], &a[Y][X]);
  swap(&a[X][Z], &a[Z][X]);
  swap(&a[Y][Z], &a[Z][Y]);
}

void Transform(MATRIX a, VECTOR u, VECTOR v)
{
  /* Applies matrix "a" to vector "u", returning vector "v" */

  v[X] = DOT(a[X], u);
  v[Y] = DOT(a[Y], u);
  v[Z] = DOT(a[Z], u);
}

void MultMat(MATRIX a, MATRIX b, MATRIX c)
{
  /* Multiples matrices "a" and "b", returning the result in matrix "c" */

  c[X][X] = a[X][X] * b[X][X] + a[X][Y] * b[Y][X] + a[X][Z] * b[Z][X];
  c[X][Y] = a[X][X] * b[X][Y] + a[X][Y] * b[Y][Y] + a[X][Z] * b[Z][Y];
  c[X][Z] = a[X][X] * b[X][Z] + a[X][Y] * b[Y][Z] + a[X][Z] * b[Z][Z];
  c[Y][X] = a[Y][X] * b[X][X] + a[Y][Y] * b[Y][X] + a[Y][Z] * b[Z][X];
  c[Y][Y] = a[Y][X] * b[X][Y] + a[Y][Y] * b[Y][Y] + a[Y][Z] * b[Z][Y];
  c[Y][Z] = a[Y][X] * b[X][Z] + a[Y][Y] * b[Y][Z] + a[Y][Z] * b[Z][Z];
  c[Z][X] = a[Z][X] * b[X][X] + a[Z][Y] * b[Y][X] + a[Z][Z] * b[Z][X];
  c[Z][Y] = a[Z][X] * b[X][Y] + a[Z][Y] * b[Y][Y] + a[Z][Z] * b[Z][Y];
  c[Z][Z] = a[Z][X] * b[X][Z] + a[Z][Y] * b[Y][Z] + a[Z][Z] * b[Z][Z];
}

void VecToMat(VECTOR v1, VECTOR v2, MATRIX a)
{
  /* Multiplies column vector "v1" and row vector "v2", forming matrix "a" */

  a[X][X] = v1[X] * v2[X];
  a[X][Y] = v1[X] * v2[Y];
  a[X][Z] = v1[X] * v2[Z];
  a[Y][X] = v1[Y] * v2[X];
  a[Y][Y] = v1[Y] * v2[Y];
  a[Y][Z] = v1[Y] * v2[Z];
  a[Z][X] = v1[Z] * v2[X];
  a[Z][Y] = v1[Z] * v2[Y];
  a[Z][Z] = v1[Z] * v2[Z];
}

void GetBasis(VECTOR v1, VECTOR v2, VECTOR v3)
{
  /* Given vector "v1", this routine returns basis "v1", "v2", "v3" */

  static const MATRIX i = {{1,0,0},  /* Unit matrix */
			   {0,1,0},
			   {0,0,1}};

  /* Get spanning set...first guess: choose i[Y] and i[Z] as 2nd & 3rd vecs */

  COPY_VEC(i[Y], v2);
  COPY_VEC(i[Z], v3);

  /* If v1 is actually null, set 1st vector to i[X] and return */

  if (v1[X] == 0 && v1[Y] == 0 && v1[Z] == 0) {
    COPY_VEC(i[X], v1);
    return;
  }

  /*
   * If v1 does not have an X component, make 2nd vector i[X]. If in
   * addition v1 does not have a Y component, make 3rd vector i[Y].
   * Now v1, v2, and v3 span 3-space.
   *
   */

  if (v1[X] == 0) {
    COPY_VEC(i[X], v2);
    if (v1[Y] == 0)
      COPY_VEC(i[Y], v3);
  }

  /* Construct the orthonormal basis */

  MakeBasis(v1, v2, v3);
}

void MakeBasis(VECTOR v1, VECTOR v2, VECTOR v3)
{
  /*
   * Uses the Gram-Schmidt process to convert the spanning set v1, v2, v3
   * into an orthonormal basis set.
   *
   */

  VECTOR n, v;
  double proj;

  /* Construct first basis vector */

  NORM_VEC(v1, MAG(v1));

  /* Construct second basis vector */

  COPY_VEC(v1, n);
  proj = DOT(v2, n);
  SCALE_VEC(n, proj);
  SUB_VEC(v2, n, v2);
  NORM_VEC(v2, MAG(v2));

  /* Construct third basis vector */

  COPY_VEC(v3, v);
  COPY_VEC(v2, n);
  proj = DOT(v, n);
  SCALE_VEC(n, proj);
  SUB_VEC(v3, n, v3);
  COPY_VEC(v1, n);
  proj = DOT(v, n);
  SCALE_VEC(n, proj);
  SUB_VEC(v3, n, v3);
  NORM_VEC(v3, MAG(v3));
}

static void swap(double *x, double *y)
{
  /* Swaps double values pointed to by "x" and "y" */

  double t = *x;

  *x = *y;
  *y = t;
}

/* vector.c */
