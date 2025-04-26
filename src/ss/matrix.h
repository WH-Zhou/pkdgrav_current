#ifndef MATRIX_HINCLUDED
#define MATRIX_HINCLUDED

#ifndef FLOATTYPE_HINCLUDED /* default to double -- see pkdgrav/floattype.h */
typedef double FLOAT;
#endif

#ifndef LINALG_HINCLUDED /* see pkdgrav/linalg.h */
typedef FLOAT Scalar;
typedef Scalar Vector[3];
typedef Vector Matrix[3];
#endif

void matrixInvert(Matrix m);
void matrixDiagonalize(Matrix m, Vector vEvals, Matrix mEvecs, int bSort);

#endif /* !MATRIX_HINCLUDED */
