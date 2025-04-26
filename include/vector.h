/*
 * vector.h -- DCR 94-04-18
 * ========================
 * Useful macros for common operations on 3-vectors.
 *
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

/*** Copyright notice ***/

#include <copyright.h>

/*** Other header files to include ***/

#include <math.h>
#include "scalar.h"

/*** Three dimensions... ***/

#ifndef N_DIM
#  define N_DIM 3
#endif

/*** Vector and matrix type definitions ***/

typedef double VECTOR[N_DIM];
typedef double MATRIX[N_DIM][N_DIM];

/*** Vector component definitions ***/

#define X 0
#define Y 1
#define Z 2

/*** Macro definitions ***/

/* Assigns a value (x,y,z) to vector v */

#define SET_VEC(v, x, y, z) {\
  (v)[X] = (x);\
  (v)[Y] = (y);\
  (v)[Z] = (z);\
}

/* Assigns zero to vector v */

#define ZERO_VEC(v) SET_VEC((v), 0, 0, 0)

/* Copies vector v1 to vector v2 */

#define COPY_VEC(v1, v2) {\
  (v2)[X] = (v1)[X];\
  (v2)[Y] = (v1)[Y];\
  (v2)[Z] = (v1)[Z];\
}

/* Adds vectors v1 & v2 and puts the result in vector v */

#define ADD_VEC(v1, v2, v) {\
  (v)[X] = (v1)[X] + (v2)[X];\
  (v)[Y] = (v1)[Y] + (v2)[Y];\
  (v)[Z] = (v1)[Z] + (v2)[Z];\
}

/* Subtracts vector v2 from vector v1 and puts the result in vector v */

#define SUB_VEC(v1, v2, v) {\
  (v)[X] = (v1)[X] - (v2)[X];\
  (v)[Y] = (v1)[Y] - (v2)[Y];\
  (v)[Z] = (v1)[Z] - (v2)[Z];\
}

/* Multiplies vector v by scalar a */

#define SCALE_VEC(v, a) {\
  (v)[X] *= (a);\
  (v)[Y] *= (a);\
  (v)[Z] *= (a);\
}

/* Divides vector v by scalar a */

#define NORM_VEC(v, a) {\
  double _n = 1.0 / (a);\
  SCALE_VEC((v), _n);\
}

/* Returns dot product of vectors v1 & v2 */

#define DOT(v1, v2) ((v1)[X] * (v2)[X] + (v1)[Y] * (v2)[Y] + (v1)[Z] * (v2)[Z])

/* Returns cross product of vectors v1 & v2 in vector v */

#define CROSS(v1, v2, v) {\
  (v)[X] = (v1)[Y] * (v2)[Z] - (v1)[Z] * (v2)[Y];\
  (v)[Y] = (v1)[Z] * (v2)[X] - (v1)[X] * (v2)[Z];\
  (v)[Z] = (v1)[X] * (v2)[Y] - (v1)[Y] * (v2)[X];\
}

/* Returns z-component of cross product of vectors v1 & v2 */

#define CROSS_Z(v1, v2) ((v1)[X] * (v2)[Y] - (v1)[Y] * (v2)[X])

/* Returns square magnitude of vector v */

#define MAG_SQ(v) (DOT((v), (v)))

/* Returns magnitude of vector v */

#define MAG(v) (sqrt(MAG_SQ(v)))

/* Returns square distance between vectors v1 & v2 */

#define DIST_SQ(v1, v2)\
  (SQ((v1)[X] - (v2)[X]) + SQ((v1)[Y] - (v2)[Y]) + SQ((v1)[Z] - (v2)[Z]))

/* Returns distance between vectors v1 & v2 */

#define DIST(v1, v2) (sqrt(DIST_SQ((v1), (v2))))

/* Returns TRUE if square distance between vectors v1 & v2 smaller than
   square of sum of scalars r1 & r2 */

#define OVERLAP(v1, r1, v2, r2) (DIST_SQ((v1), (v2)) < SQ((r1) + (r2)))

/* Zeroes 3x3 matrix a */

#define ZERO_MAT(a) {\
  ZERO_VEC((a)[X]);\
  ZERO_VEC((a)[Y]);\
  ZERO_VEC((a)[Z]);\
}

/* Makes matrix a the unit matrix */

#define UNIT_MAT(a) {\
  ZERO_MAT(a);\
  (a)[X][X] = (a)[Y][Y] = (a)[Z][Z] = 1;\
}

/* Makes a skew-symmetric matrix a out of vector v */

#define SKEW_SYM_MAT(v, a) {\
  ZERO_MAT(a);\
  (a)[X][Y] =   (v)[Z];\
  (a)[X][Z] = - (v)[Y];\
  (a)[Y][X] = - (v)[Z];\
  (a)[Y][Z] =   (v)[X];\
  (a)[Z][X] =   (v)[Y];\
  (a)[Z][Y] = - (v)[X];\
}

/* Copies matrix a to matrix b */

#define COPY_MAT(a, b) {\
  COPY_VEC((a)[X], (b)[X]);\
  COPY_VEC((a)[Y], (b)[Y]);\
  COPY_VEC((a)[Z], (b)[Z]);\
}

/* Multiplies matrix a by scalar b */

#define SCALE_MAT(a, b) {\
  SCALE_VEC((a)[X], (b));\
  SCALE_VEC((a)[Y], (b));\
  SCALE_VEC((a)[Z], (b));\
}

/* Adds matrix a to b and puts the result in c */

#define ADD_MAT(a, b, c) {\
  ADD_VEC((a)[X], (b)[X], (c)[X]);\
  ADD_VEC((a)[Y], (b)[Y], (c)[Y]);\
  ADD_VEC((a)[Z], (b)[Z], (c)[Z]);\
}

/* Subtracts matrix b from a and puts the result in c */

#define SUB_MAT(a, b, c) {\
  SUB_VEC((a)[X], (b)[X], (c)[X]);\
  SUB_VEC((a)[Y], (b)[Y], (c)[Y]);\
  SUB_VEC((a)[Z], (b)[Z], (c)[Z]);\
}

/*** Prototypes ***/

extern void Transpose(MATRIX a);
extern void Transform(MATRIX a, VECTOR u, VECTOR v);
extern void MultMat(MATRIX a, MATRIX b, MATRIX c);
extern void VecToMat(VECTOR v1, VECTOR v2, MATRIX a);
extern void GetBasis(VECTOR v1, VECTOR v2, VECTOR v3);
extern void MakeBasis(VECTOR v1, VECTOR v2, VECTOR v3);

#endif /* _VECTOR_H_ */

/* vector.h */
