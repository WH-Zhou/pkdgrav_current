/*
** Linear algebra routines for matrix inversion and diagonalization.
** Includes a standalone test driver.
*/

/*#define STANDALONE*/ /* uncomment to compile standalone test */

#ifdef NUMREC
/* compile with -DNUMREC */
#endif

#ifdef GSL
/* compile with -DGSL -I/path/to/headers -L/path/to/lib -lgsl -lgslcblas */
#endif

#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#ifdef NUMREC
#include <math.h>
#ifndef N_DIM
#define N_DIM 3
typedef Vector VECTOR;
typedef Matrix MATRIX;
#define COPY_VEC(u, v) {\
(v)[0] = (u)[0];\
(v)[1] = (u)[1];\
(v)[2] = (u)[2];\
}
#define COPY_MAT(a, b) {\
COPY_VEC((a)[0], (b)[0]);\
COPY_VEC((a)[1], (b)[1]);\
COPY_VEC((a)[2], (b)[2]);\
}
#endif

/*
** Following routines adapted from Numerical Recipes in C (2nd ed).
*/

static void ludcmp(MATRIX a,int *indx)
{
    /* based on ludcmp(), NRiC(2e) 2.3 */
    
    const int n = N_DIM;
    
    VECTOR vv;
    double big,dum,sum,temp;
    int imax=0,i,j,k;
    
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) {
            fprintf(stderr, "ludcmp(): Singular matrix.\n");
            exit(1);
        }
        vv[i]=1.0/big;
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            sum=a[i][j];
            for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            sum=a[i][j];
            for (k=0;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ((dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=0;k<n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        /*if (a[j][j] == 0.0) a[j][j]=TINY;*/
        if (j != n - 1) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<n;i++) a[i][j] *= dum;
        }
    }
}

static void lubksb(MATRIX a,int *indx,VECTOR b)
{
    /* based on lubksb(), NRiC(2e) 2.3 */
    
    const int n = N_DIM;
    
    double sum;
    int ip,ii=(-1),i,j;
    
    for (i=0;i<n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii >= 0)
            for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum != 0.0) ii=i;
        b[i]=sum;
    }
    for (i=n-1;i>=0;i--) {
        sum=b[i];
        for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}

static void invert(MATRIX a)
{
    /* inverts matrix */

    MATRIX b;
    VECTOR col;
    int indx[N_DIM],i,j;
    
    ludcmp(a,indx);
    
    for (j=0;j<N_DIM;j++) {
        for (i=0;i<N_DIM;i++)
            col[i] = 0;
        col[j] = 1;
        lubksb(a,indx,col);
        for (i=0;i<N_DIM;i++)
            b[i][j] = col[i];
    }
    
    COPY_MAT(b,a);
}

/*
** Following routines adapted from Numerical Recipes (3rd ed).
*/

static void eigsrt(VECTOR d,MATRIX v)
{
    /* based on eigsrt(), NR3 11.1 */
    
    const int n = N_DIM; /* i.e., hard-wired for 3-D */
    
    int i,j,k;
    double p;
    
    for (i=0;i<n-1;i++) {
        p=d[k=i];
        for (j=i;j<n;j++)
            if (d[j] >= p) p=d[k=j];
        if (k != i) {
            d[k]=d[i];
            d[i]=p;
            for (j=0;j<n;j++) {
                p=v[j][i];
                v[j][i]=v[j][k];
                v[j][k]=p;
            }
        }
    }
}

static void rot(MATRIX a,double s,double tau,int i,int j,int k,int l)
{
    /* from struct Jacobi, NR3 11.1 */
    
    double g=a[i][j];
    double h=a[k][l];
    a[i][j]=g-s*(h+g*tau);
    a[k][l]=h+s*(g-h*tau);
}

static void jacobi(MATRIX a,VECTOR d,MATRIX v,int bSort)
{
    /* based on struct jacobi(), NR3 11.1 */
    
    const int n = N_DIM; /* i.e., hard-wired for 3-D */
    const double EPS = 1.1102e-16; /* numeric_limits<Doub>::epsilon() */
    
    int i,j,ip,iq;
    double tresh,theta,tau,t,sm,s,h,g,c;
    VECTOR b,z;
    
    for (ip=0;ip<n;ip++) {
        for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
        v[ip][ip]=1.0;
    }
    for (ip=0;ip<n;ip++) {
        b[ip]=d[ip]=a[ip][ip];
        z[ip]=0.0;
    }
    for (i=1;i<=50;i++) {
        sm=0.0;
        for (ip=0;ip<n-1;ip++) {
            for (iq=ip+1;iq<n;iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0) {
            if (bSort) eigsrt(d,v);
            return;
        }
        if (i < 4)
            tresh=0.2*sm/(n*n);
        else
            tresh=0.0;
        for (ip=0;ip<n-1;ip++) {
            for (iq=ip+1;iq<n;iq++) {
                g=100.0*fabs(a[ip][iq]);
                if (i > 4 && g <= EPS*fabs(d[ip]) && g <= EPS*fabs(d[iq]))
                    a[ip][iq]=0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h=d[iq]-d[ip];
                    if (g <= EPS*fabs(h))
                        t=(a[ip][iq])/h;
                    else {
                        theta=0.5*h/(a[ip][iq]);
                        t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0) t = -t;
                    }
                    c=1.0/sqrt(1+t*t);
                    s=t*c;
                    tau=s/(1.0+c);
                    h=t*a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq]=0.0;
                    for (j=0;j<ip;j++)
                        rot(a,s,tau,j,ip,j,iq);
                    for (j=ip+1;j<iq;j++)
                        rot(a,s,tau,ip,j,j,iq);
                    for (j=iq+1;j<n;j++)
                        rot(a,s,tau,ip,j,iq,j);
                    for (j=0;j<n;j++)
                        rot(v,s,tau,j,ip,j,iq);
                }
            }
        }
        for (ip=0;ip<n;ip++) {
            b[ip] += z[ip];
            d[ip]=b[ip];
            z[ip]=0.0;
        }
    }
    fprintf(stderr, "jacobi(): Too many iterations.\n");
    exit(1);
}

void matrixInvert(Matrix m)
{
    invert(m);
}

void matrixDiagonalize(Matrix m, Vector vEvals, Matrix mEvecs, int bSort)
{
    Matrix mCopy;
    
    COPY_MAT(m, mCopy);
    jacobi(mCopy, vEvals, mEvecs, bSort); /* mCopy is changed by routine */
}

#endif /* NUMREC */

#ifdef GSL

/* hardwired for 3-D vectors and matrices */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

static void gsl_to_vector(const gsl_vector *A, Vector v)
{
    v[0] = gsl_vector_get(A, 0);
    v[1] = gsl_vector_get(A, 1);
    v[2] = gsl_vector_get(A, 2);
}

static void matrix_to_gsl(Matrix m, gsl_matrix *A)
{
    gsl_matrix_set(A, 0, 0, m[0][0]);
    gsl_matrix_set(A, 1, 0, m[1][0]);
    gsl_matrix_set(A, 2, 0, m[2][0]);
    gsl_matrix_set(A, 0, 1, m[0][1]);
    gsl_matrix_set(A, 1, 1, m[1][1]);
    gsl_matrix_set(A, 2, 1, m[2][1]);
    gsl_matrix_set(A, 0, 2, m[0][2]);
    gsl_matrix_set(A, 1, 2, m[1][2]);
    gsl_matrix_set(A, 2, 2, m[2][2]);
}

static void gsl_to_matrix(const gsl_matrix *A, Matrix m)
{
    m[0][0] = gsl_matrix_get(A, 0, 0);
    m[1][0] = gsl_matrix_get(A, 1, 0);
    m[2][0] = gsl_matrix_get(A, 2, 0);
    m[0][1] = gsl_matrix_get(A, 0, 1);
    m[1][1] = gsl_matrix_get(A, 1, 1);
    m[2][1] = gsl_matrix_get(A, 2, 1);
    m[0][2] = gsl_matrix_get(A, 0, 2);
    m[1][2] = gsl_matrix_get(A, 1, 2);
    m[2][2] = gsl_matrix_get(A, 2, 2);
}

void matrixInvert(Matrix m)
{
    /*
     ** Note this comment from the GSL documentation: "It is preferable
     ** to avoid direct use of the inverse whenever possible, as the
     ** linear solver functions can obtain the same result more
     ** efficiently and reliably (consult any introductory textbook on
     ** numerical linear algebra for details)."
     */
    
    gsl_matrix *A = gsl_matrix_alloc(3, 3);
    gsl_matrix *inverse = gsl_matrix_alloc(3, 3);
    gsl_permutation *p = gsl_permutation_alloc(3);
    int i, j, signum;
    
    matrix_to_gsl(m, A);
    if (gsl_linalg_LU_decomp(A, p, &signum) != GSL_SUCCESS) {
        fprintf(stderr, "matrixInvert(): Unable to perform LU decomposition.");
        exit(1);
    }
    if (gsl_linalg_LU_invert(A, p, inverse) != GSL_SUCCESS) {
        fprintf(stderr, "matrixInvert(): Unable to perform matrix inversion.");
        exit(1);
    }
    gsl_matrix_free(A);
    gsl_permutation_free(p);
    gsl_to_matrix(inverse, m);
    gsl_matrix_free(inverse);
}

void matrixDiagonalize(Matrix m, Vector vEvals, Matrix mEvecs, int bSort)
{
    /*
     ** Returns eigenvectors corresponding to eigenvalues vEvals as
     ** columns of mEvecs, sorted in decreasing order of eigenvalue
     ** magnitude if bSort != 0. Assumes m is a real, SYMMETRIC matrix
     ** (e.g., an inertia tensor) -- this means the eigenvectors are
     ** real and orthonormal. Note m is not changed on return.
     */
    
    gsl_matrix *A = gsl_matrix_alloc(3, 3);
    gsl_vector *eval = gsl_vector_alloc(3);
    gsl_matrix *evec = gsl_matrix_alloc(3, 3);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);
    
    matrix_to_gsl(m, A);
    gsl_eigen_symmv(A, eval, evec, w); /* modifies A */
    gsl_matrix_free(A);
    gsl_eigen_symmv_free(w);
    if (bSort != 0)
        gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
    gsl_to_vector(eval, vEvals);
    gsl_to_matrix(evec, mEvecs);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
}

#endif /* GSL */

#ifdef STANDALONE

static void matrix_show(Matrix m)
{
    int i, j;
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            printf("%10.2e", m[i][j]);
        printf("\n");
    }
}

static void matrix_transpose(Matrix m)
{
    double t;
    t = m[1][0];
    m[1][0] = m[0][1];
    m[0][1] = t;
    t = m[2][0];
    m[2][0] = m[0][2];
    m[0][2] = t;
    t = m[2][1];
    m[2][1] = m[1][2];
    m[1][2] = t;
}

static void matrix_apply(const Matrix m, const Vector u, Vector v)
{
    v[0] = m[0][0] * u[0] + m[1][0] * u[1] + m[2][0] * u[2];
    v[1] = m[0][1] * u[0] + m[1][1] * u[1] + m[2][1] * u[2];
    v[2] = m[0][2] * u[0] + m[1][2] * u[1] + m[2][2] * u[2];
}

static double dot(const Vector u, const Vector v)
{
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

static double mag(const Vector v)
{
    return sqrt(dot(v, v));
}

int main(void) {
    Vector vEvals, v;
    Matrix mEvecs, A = {{3.0, 2.0, 1.0}, {2.0, 3.0, 2.0}, {1.0, 2.0, 3.0}};
    int i;
    printf("Input matrix =\n");
    matrix_show(A);
    matrixInvert(A);
    printf("Inverted matrix =\n");
    matrix_show(A);
    matrixInvert(A);
    printf("Re-inverted =\n");
    matrix_show(A);
    matrixDiagonalize(A, vEvals, mEvecs, 1 /* sort */);
    printf("Eigenvalues =\n");
    printf("%10.2e%10.2e%10.2e\n", vEvals[0], vEvals[1], vEvals[2]);
    matrix_transpose(mEvecs); /* columns to rows */
    printf("Eigenvectors (rows: E1, E2, E3) =\n");
    matrix_show(mEvecs);
    for (i = 0; i < 3; i++) {
        matrix_apply(A, mEvecs[i], v);
        printf("Diagonalize check: compare with E%i:\n", i + 1);
        printf("%10.2e%10.2e%10.2e\n", v[0] / vEvals[i], v[1] / vEvals[i],
               v[2] / vEvals[i]);
    }
    printf("Orthonormalization checks:\n");
    for (i = 0; i < 3; i++)
        printf("E%i magnitude = %g\n", i + 1, mag(mEvecs[i]));
    printf("E1 dotted with E2 = %g\n", dot(mEvecs[0], mEvecs[1]));
    printf("E1 dotted with E3 = %g\n", dot(mEvecs[0], mEvecs[2]));
    printf("E2 dotted with E3 = %g\n", dot(mEvecs[1], mEvecs[2]));
    printf("NOTE: eigenvectors only unique to within a sign\n");
    return 0;
}

#endif /* STANDALONE */
