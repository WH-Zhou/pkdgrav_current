/*
** Routines for finding roots of polynomials with real coefficients.
** Includes a high-precision version of the quadratic formula and
** limited versions of cubic and quartic solvers (for real roots
** only), along with a full solver using complex arithmetic, and a
** standalone test driver.
*/

/*#define STANDALONE*/ /* uncomment to compile standalone test */

#ifdef NUMREC
/* compile with -DNUMREC */
#endif

#ifdef GSL
/* compile with -DGSL -I/path/to/headers -L/path/to/lib -lgsl */
#endif

#define EPS (2.0e-12) /* for deciding if root is real */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "polyroots.h" /* not really needed at present */

#ifdef GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#endif

#ifdef NUMREC

/*
** Following routines adapted from Numerical Recipes in C (2nd ed).
*/

typedef struct {
  double r,i;
  } dComplex;

static dComplex Complex(double re,double im)
{
    dComplex c;
    
    c.r = re;
    c.i = im;
    
    return c;
    }

static dComplex ComplexAdd(dComplex a,dComplex b)
{
    dComplex c;
    
    c.r = a.r + b.r;
    c.i = a.i + b.i;
    
    return c;
    }

static dComplex ComplexSub(dComplex a,dComplex b)
{
    dComplex c;
    
    c.r = a.r - b.r;
    c.i = a.i - b.i;
    
    return c;
    }

static dComplex ComplexMul(dComplex a,dComplex b)
{
    dComplex c;
    
    c.r = a.r*b.r - a.i*b.i;
    c.i = a.i*b.r + a.r*b.i;

    return c;
    }

static dComplex ComplexDiv(dComplex a,dComplex b)
{
    dComplex c;
    double r,den;
    
    if (fabs(b.r) >= fabs(b.i)) {
	r = b.i/b.r;
	den = b.r + r*b.i;
	c.r = (a.r + r*a.i)/den;
	c.i = (a.i - r*a.r)/den;
	}
    else {
	r = b.r/b.i;
	den = b.i + r*b.r;
	c.r = (a.r*r + a.i)/den;
	c.i = (a.i*r - a.r)/den;
	}
    
    return c;
    }

static double ComplexAbs(dComplex z)
{
    double x,y,ans,temp;
    
    x = fabs(z.r);
    y = fabs(z.i);
    
    if (x == 0.0)
	ans = y;
    else if (y == 0.0)
	ans = x;
    else if (x > y) {
	temp = y/x;
	ans = x*sqrt(1.0 + temp*temp);
	}
    else {
	temp = x/y;
	ans = y*sqrt(1.0 + temp*temp);
	}
    
    return ans;
    }

static dComplex ComplexSqrt(dComplex z)
{
    dComplex c;
    double x,y,w,r;
    
    if (z.r == 0.0 && z.i == 0.0) {
	c.r = c.i = 0.0;
	return c;
	}
    
    x = fabs(z.r);
    y = fabs(z.i);
    if (x >= y) {
	r = y/x;
	w = sqrt(x)*sqrt(0.5*(1.0 + sqrt(1.0 + r*r)));
	}
    else {
	r = x/y;
	w = sqrt(y)*sqrt(0.5*(r + sqrt(1.0 + r*r)));
	}
    if (z.r >= 0.0) {
	c.r = w;
	c.i = z.i/(2.0*w);
	}
    else {
	c.i = ((z.i >= 0) ? w : -w);
	c.r = z.i/(2.0*c.i);
	}
    
    return c;
    }

static dComplex RealComplexMul(double x,dComplex a)
{
    dComplex c;
    
    c.r = x*a.r;
    c.i = x*a.i;
    
    return c;
    }

static double dMax(double a,double b)
{
    return (a > b ? a : b);
    }

#define EPSS 1.0e-15 /* estimated fractional roundoff error */
#define MR 8 /* try this many different fractional values... */
#define MT 10 /* ...once every this number of steps */
#define MAXIT (MT*MR) /* max allowed iterations */

static void laguer(dComplex a[],int m,dComplex *x,int *its)
/*
** Based on laguer(). NRiC(2e) 9.5.
**
** Given the degree m and the m+1 complex coefficients a[0..m] of the
** polynomial (a[0] is the constant term), and given a complex value
** x, this routine improves x by Laguerre's method until it converges,
** within the achievable roundoff limit, to a root of the given
** polynomial.  The number of iterations taken is returned as "its".
*/
{
    int iter,j;
    double abx,abp,abm,err;
    dComplex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
    static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0}; /* fractions used to break a limit cycle */
    
    for (iter=1;iter<=MAXIT;iter++) { /* loop over iterations up to allowed maximum */
	*its = iter;
	b = a[m];
	err = ComplexAbs(b);
	d = f = Complex(0.0,0.0);
	abx = ComplexAbs(*x);
	for (j=m-1;j>=0;j--) { /* efficient computation of the polynomial and its first two derivatives */
	    f = ComplexAdd(ComplexMul(*x,f),d);
	    d = ComplexAdd(ComplexMul(*x,d),b);
	    b = ComplexAdd(ComplexMul(*x,b),a[j]);
	    err = ComplexAbs(b) + abx*err;
	    }
	err *= EPSS; /* estimate of roundoff error in evaluating polynomial */
	if (ComplexAbs(b) <= err)
	    return; /* we are on the root */
	g = ComplexDiv(d,b); /* the generic case: use Laguerre's formula */
	g2 = ComplexMul(g,g);
	h = ComplexSub(g2,RealComplexMul(2.0,ComplexDiv(f,b)));
	sq = ComplexSqrt(RealComplexMul((double) (m - 1),ComplexSub(RealComplexMul((double) m,h),g2)));
	gp = ComplexAdd(g,sq);
	gm = ComplexSub(g,sq);
	abp = ComplexAbs(gp);
	abm = ComplexAbs(gm);
	if (abp < abm)
	    gp = gm;
	dx = ((dMax(abp,abm) > 0.0 ? ComplexDiv(Complex((double) m,0.0),gp) :
	       RealComplexMul(exp(log(1.0 + abx)),Complex(cos((double) iter),sin((double) iter)))));
	x1 = ComplexSub(*x,dx);
	if (x->r == x1.r && x->i == x1.i)
	    return; /* converged */
	/*
	** Every so often we take a fractional step, to break any
	** limit cycle (itself a rare occurrence)...
	*/
	if (iter % MT)
	    *x = x1;
	else
	    *x = ComplexSub(*x,RealComplexMul(frac[iter/MT],dx));
	}
    /*
    ** Very unusual to get here: can occur only for complex roots.
    ** Try a different starting guess for the root.
    */
    fprintf(stderr,"laguer(): too many iterations.\n");
    exit(1);
    }

#undef MAXIT
#undef MT
#undef MR
#undef EPSS

#endif /* NUMREC */

static double dsgn(double x)
{
	return (x < 0.0 ? -1.0 : 1.0); /* zero taken as positive */
	}

#ifdef GSL
int polyQuadSolve(double a, double b, double c, double *x1, double *x2)
{
    /*
     ** Finds real roots of quadratic equation with real coefficients.
     ** Returns 1 (and does not alter x1 and x2) if no real roots.
     ** Note *x1 <= *x2 on return.
     */
    
    int nRoots;
    nRoots = gsl_poly_solve_quadratic(a, b, c, x1, x2);
    if (nRoots == 0)
        return 1; /* no real roots */
    if (nRoots == 1)
        *x2 = *x1;
    return 0;
    }
#else
int polyQuadSolve(double a,double b,double c,double *x1,double *x2)
{
    /*
    ** Optimized quadratic formula for real coefficients and roots.
    ** Returns 1 (and does not alter x1 and x2) if no real roots.
    ** Note *x1 <= *x2 on return, by construction;
    */
    
    double D,q;
    
    assert(a != 0.0);
    
    if (b == 0.0 && c == 0.0) { /* special case x^2 = 0 */
	*x1 = *x2 = 0.0;
	return 0;
	}
    
    D = b*b - 4.0*a*c;
    
    if (D < 0.0)
	return 1; /* no real roots */
    
    q = -0.5*(b + dsgn(b)*sqrt(D));
    assert(q != 0.0);
    *x1 = q/a;
    *x2 = c/q;
    
    if (*x1 > *x2) { /* sort in ascending order */
	double dTmp = *x2;
	*x2 = *x1;
	*x1 = dTmp;
	}
    
    return 0;
    }
#endif /* !GSL */

#ifdef GSL
int polyCubicSolveLimited(double a1, double a2, double a3, double *x1, double *x2, double *x3)
{
    /* real coefficients and roots only (x^3 + a1x^2 + a2x + a3 = 0) */
    
    int nRoots;
    nRoots = gsl_poly_solve_cubic(a1, a2, a3, x1, x2, x3);
    if (nRoots == 1)
        return 1; /* 1 real, 2 complex roots (x2, x3 unmodified) */
    assert(nRoots == 3);
    return 0;
    }
#else
int polyCubicSolveLimited(double a1,double a2,double a3,double *x1,double *x2,double *x3)
{
    /* real coefficients and roots only (x^3 + a1x^2 + a2x + a3c = 0) -- not recommended */
    
    double Q,R,dTheta;
    
    Q = (a1*a1 - 3*a2)/9;
    R = ((2*a1*a1 - 9*a2)*a1 + 27*a3)/54;
    
    if (Q == 0.0 && R == 0.0) { /* special case: perfect cube */
	*x1 = *x2 = *x3 = -a1/3;
	return 0;
	}
    
    if (R*R >= Q*Q*Q)
	return 1; /* 1 real, 2 complex roots */
    
    assert(Q > 0.0);
    
    dTheta = acos(R/sqrt(Q*Q*Q));
    
    *x1 = -2*sqrt(Q)*cos(dTheta/3) - a1/3;
    *x2 = -2*sqrt(Q)*cos(dTheta/3 + 2*M_PI/3) - a1/3;
    *x3 = -2*sqrt(Q)*cos(dTheta/3 - 2*M_PI/3) - a1/3;
    
    return 0;
    }
#endif /* !GSL */

int polyQuarticSolveLimited(double a1,double a2,double a3,double a4,double *x1,double *x2,double *x3,double *x4)
{
    /* real coefficients and roots, for special cases only (x^4 + a1x^3 + a2x^2 + a3x + a4 = 0) */
    
    double cr1,cr2,cr3,cr;
    
    if (a1 == 0.0 && a2 == 0.0 && a3 == 0.0 && a4 == 0.0) {
        /* trivial case */
        *x1 = *x2 = *x3 = *x4 = 0.0;
        return 0;
	}
    
    if (a4 == 0.0) {
        /* 0 is a root if the last coefficient is zero */
        if (polyCubicSolveLimited(a1,a2,a3,x2,x3,x4))
            return 1;
        *x1 = 0.0;
        return 0;
	}
    
    if (polyCubicSolveLimited(-a2,a1*a3 - 4*a4,4*a2*a4 - a3*a3 - a1*a1*a4,&cr1,&cr2,&cr3))
	return 1;
    
    if ((a1*a1 - 4*a2 + 4*cr1 >= 0.0 && cr1*cr1 - 4*a4 > 0.0))
	cr = cr1;
    else if (a1*a1 - 4*a2 + 4*cr2 >= 0.0 && cr2*cr2 - 4*a4 > 0.0)
	cr = cr2;
    else if (a1*a1 - 4*a2 + 4*cr3 >= 0.0 && cr3*cr3 - 4*a4 > 0.0)
	cr = cr3;
    else
	return 1;
    
    if (polyQuadSolve(1.0,0.5*(a1 + sqrt(a1*a1 - 4*a2 + 4*cr)),0.5*(cr - sqrt(cr*cr - 4*a4)),x1,x2))
	return 1;
    if (polyQuadSolve(1.0,0.5*(a1 - sqrt(a1*a1 - 4*a2 + 4*cr)),0.5*(cr + sqrt(cr*cr - 4*a4)),x3,x4))
	return 1;
    
    return 0;
    }

#ifdef NUMREC

void polyFindRealRoots(int nDegree,double dCoefs[],double dRoots[],int *nRealRoots)
{
    /*
    ** Based on zroots(), NRiC(2e) 9.5.
    **
    ** Uses Laguerre's method (via a modified version of zroots()) to
    ** get all real roots of the nDegree polynomial equation with real
    ** coefficients dCoefs (dCoefs[0] is the constant term).  The
    ** nRealRoots found are stored in ascending order in dRoots (which
    ** should have space for nDegree roots).
    */
    
    dComplex *a,*roots,x,b,c,*ad;
    double r;
    int i,j,jj,its;
    
    assert(dCoefs[nDegree] != 0.0); /* degenerate case */
    
    /* store coefficients in complex array */

    a = (dComplex *) malloc((nDegree + 1)*sizeof(dComplex));
    assert(a != NULL);
    
    for (i=0;i<=nDegree;i++)
	a[i] = Complex(dCoefs[i],0.0);
    
    /* allocate space for the roots */
    
    roots = (dComplex *) malloc(nDegree*sizeof(dComplex));
    assert(roots != NULL);
    
    *nRealRoots = 0; /* initialize */
    
    /* following code modeled after zroots()... */
    
    ad = (dComplex *) malloc((nDegree + 1)*sizeof(dComplex));
    assert(ad != NULL);
    
    for (j=0;j<=nDegree;j++)
	ad[j] = a[j]; /* copy of coefficients for successive deflation */
    for (j=nDegree;j>=1;j--) { /* loop over each root to be found */
	x = Complex(0.0,0.0); /* start at zero to favor convergence to smallest remaining root */
	laguer(ad,j,&x,&its); /* find the root */
	if (fabs(x.i) <= 2.0*EPS*fabs(x.r))
	    x.i = 0.0;
	roots[j-1] = x;
	/* forward deflation */
	b = ad[j];
	for (jj=j-1;jj>=0;jj--) {
	    c = ad[jj];
	    ad[jj] = b;
	    b = ComplexAdd(ComplexMul(x,b),c);
	    }
	}
    free((void *) ad);
    for (j=0;j<nDegree;j++) /* polish the (real) roots using the undeflated coefficients */
	if (roots[j].i == 0.0) {
	    laguer(a,nDegree,&roots[j],&its);
	    dRoots[(*nRealRoots)++] = roots[j].r;
	    }
    free((void *) roots);
    free((void *) a);
    for (j=1;j<*nRealRoots;j++) { /* sort (real) roots by straight insertion */
	r = dRoots[j];
	for (i=j-1;i>=0;i--) {
	    if (dRoots[i] <= r)
		break;
	    dRoots[i+1] = dRoots[i];
	    }
	dRoots[i+1] = r;
	}
    }

#endif /* NUMREC */

#ifdef GSL

void polyFindRealRoots(int nDegree, double dCoefs[], double dRoots[], int *nRealRoots)
{
    /* general polynomial root solver -- only real roots returned */
    
    double *z; /* space for complex roots */
    int nCoefs, iRet, i;
    assert(nDegree > 0);
    z = malloc(2 * nDegree * sizeof(double));
    if (z == NULL) {
        fprintf(stderr, "polyFindRealRoot(): Unable to allocate space for roots.");
        exit(1);
    }
    nCoefs = nDegree + 1;
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(nCoefs);
    if (w == NULL) {
        fprintf(stderr, "polyFindRealRoots(): Unable to allocate workspace.");
        exit(1);
    }
    iRet = gsl_poly_complex_solve(dCoefs, nCoefs, w, z);
    gsl_poly_complex_workspace_free(w);
    if (iRet != GSL_SUCCESS) {
        fprintf(stderr, "polyFindRealRoots(): Unable to converge to solution.");
        exit(1);
    }
    *nRealRoots = 0;
    for (i = 0; i < nDegree; i++)
        if (fabs(z[2 * i + 1]) <= EPS * fabs(z[2 * i])) {
            dRoots[i] = z[2 * i];
            ++(*nRealRoots);
        }
    free((void *) z);
}

#endif /* GSL */

#ifdef STANDALONE

#define NUM_ITER 100000 /* if > 1, invokes timing */

#if (NUM_ITER > 1)
#include <time.h>
static double getCPU(void)
{
    return (double) clock()/CLOCKS_PER_SEC; /* convert to seconds */
}
#endif

/* remaining functions for testing purposes */

static int quad_solve_std(double a,double b,double c,double *x1,double *x2)
{
    /* standard quadratic formula -- not recommended */
    
    double D;
    
    assert(a != 0.0);
    
    D = b*b - 4.0*a*c;
    
    if (D < 0.0)
	return 1; /* no real roots */
    
    *x1 = (-b + sqrt(D))/(2.0*a);
    *x2 = (-b - sqrt(D))/(2.0*a);
    
    return 0;
    }

static void usage(const char *progname)
{
    fprintf(stderr, "Usage: %s degree highest-power-coef [ next-power-coef ... constant-coef ]\n", progname);
    exit(1);
}

#ifdef NUMREC

#define MAXM 100 /* maximum anticipated value of m */
static void zroots(dComplex a[],int m,dComplex roots[],int polish)
/*
** Based on zroots(), NRiC(2e) 9.5.
** 
** Given the degree m and the m+1 complex coefficients a[0..m] of the
** polynomial (a[0] is the constant term), this routine successively
** calls laguer() and finds all m complex roots in roots[0..m-1].  The
** boolean variable "polish" should be input as true (1) if polishing
** (also by Laguerre's method) is desired, false (0) if the roots will
** be subsequently polished by other means.
*/
{
    int i,its,j,jj;
    dComplex x,b,c,ad[MAXM];
    
    for (j=0;j<=m;j++)
	ad[j] = a[j]; /* copy of coefficients for successive deflation */
    for (j=m;j>=1;j--) { /* loop over each root to be found */
	x = Complex(0.0,0.0); /* start at zero to favor convergence to smallest remaining root */
	laguer(ad,j,&x,&its); /* find the root */
	if (fabs(x.i) <= 2.0*EPS*fabs(x.r))
	    x.i = 0.0;
	roots[j-1] = x;
	/* forward deflation */
	b = ad[j];
	for (jj=j-1;jj>=0;jj--) {
	    c = ad[jj];
	    ad[jj] = b;
	    b = ComplexAdd(ComplexMul(x,b),c);
	    }
	}
    if (polish)
	for (j=0;j<m;j++) /* polish the roots using the undeflated coefficients */
	    laguer(a,m,&roots[j],&its);
    for (j=1;j<m;j++) { /* sort roots by their real parts by straight insertion */
	x = roots[j];
	for (i=j-1;i>=0;i--) {
	    if (roots[i].r <= x.r)
		break;
	    roots[i+1] = roots[i];
	    }
	roots[i+1] = x;
	}
    }

#undef MAXM

static void laguerre_solve(double dCoefs[],int nCoefs,dComplex cRoots[],int bPolish)
{
    dComplex *cCoefs;
    int i;
    
    cCoefs = (dComplex *) malloc(nCoefs*sizeof(dComplex));
    assert(cCoefs != NULL);
    
    for (i=0;i<nCoefs;i++)
	cCoefs[i] = Complex(dCoefs[i],0.0);
    
    zroots(cCoefs,nCoefs-1,cRoots,bPolish);
    
    free((void *) cCoefs);
    }

static double eval_poly(double c[],int n,double x)
{
    double p;
    int i;
    
    p = c[i=n];
    while (i > 0)
	p = p*x + c[--i];
    
    return p;
    }

static void show_roots(double dCoefs[],int nDegree,dComplex cRoots[])
{
    int i,bComplex=0;
    
    printf("Real roots:\n");
    for (i=0;i<nDegree;i++)
	if (cRoots[i].i == 0.0)
	    printf("   %.6e (evaluates to %.6e)\n",cRoots[i].r,eval_poly(dCoefs,nDegree,cRoots[i].r));
    
    printf("Complex roots:\n");
    for (i=0;i<nDegree;i++)
	if (cRoots[i].i != 0.0) {
	    printf("   %.6e + i(%.6e)\n",cRoots[i].r,cRoots[i].i);
	    bComplex = 1;
	    }
    
    if (!bComplex)
	printf("   none\n");
    }

int main(int argc,char *argv[])
{
    dComplex *cRoots;
    double *dCoefs;
    int i,nDegree,nCoefs,iFlag;
    
#if (NUM_ITER > 1)
    double t0;
#endif
    
    setbuf(stdout,(char *) NULL);
    
    if (argc < 2)
        usage(argv[0]);
    
    nDegree = atoi(argv[1]);
    
    if (nDegree < 1) {
	fprintf(stderr,"Polynomial degree must be positive integer\n");
	usage(argv[0]);
	}
    
    nCoefs = nDegree + 1;
    
    if (argc != nCoefs + 2) {
	fprintf(stderr,"Must supply %i coefficients for %i-degree polynomial\n",nCoefs,nDegree);
        usage(argv[0]);
	}
    
    dCoefs = (double *) malloc(nCoefs*sizeof(double));
    assert(dCoefs != NULL);
    
    for (i=0;i<nCoefs;i++)
	dCoefs[i] = atof(argv[argc - i - 1]); /* reverse order */
    
    if (dCoefs[nCoefs - 1] == 0.0) {
	fprintf(stderr,"Leading coefficient cannot be zero (e.g., ax^2 + bx + c = 0, a != 0)\n");
	free((void *) dCoefs);
        usage(argv[0]);
	}
    
    cRoots = (dComplex *) malloc(nDegree*sizeof(dComplex));
    assert(cRoots != NULL);
    
#if (NUM_ITER > 1)
    t0 = getCPU();
#endif
    
    for (i=0;i<NUM_ITER;i++) {
	
	switch (nDegree) {
	case 1:
	    /* trivial case */
	    cRoots[0] = Complex(-dCoefs[0]/dCoefs[1],0.0);
	    if (i == 0)
		show_roots(dCoefs,nDegree,cRoots);
	    break;
	case 2:
	    cRoots[0].i = cRoots[1].i = 0.0; /* initialize for formula solvers */
	    if (i == 0)
		printf("Standard quadratic formula\n");
	    iFlag = quad_solve_std(dCoefs[2],dCoefs[1],dCoefs[0],&cRoots[0].r,&cRoots[1].r);
	    if (i == 0) {
		if (iFlag)
		    printf("   No real roots\n");
		else
		    show_roots(dCoefs,nDegree,cRoots);
		printf("Optimized quadratic formula\n");
		}
	    iFlag = polyQuadSolve(dCoefs[2],dCoefs[1],dCoefs[0],&cRoots[0].r,&cRoots[1].r);
	    if (i == 0) {
		if (iFlag)
		    printf("   No real roots\n");
		else
		    show_roots(dCoefs,nDegree,cRoots);
		}
	    break;
	case 3:
            cRoots[0].i = cRoots[1].i = cRoots[2].i = 0.0; /* initialize for formula solvers */
	    if (i == 0)
		printf("Cubic formula\n");
	    iFlag = polyCubicSolveLimited(dCoefs[2]/dCoefs[3],dCoefs[1]/dCoefs[3],dCoefs[0]/dCoefs[3],&cRoots[0].r,&cRoots[1].r,&cRoots[2].r);
	    if (i == 0) {
		if (iFlag)
		    printf("   1 real root, 2 complex roots -- skipped\n");
		else
		    show_roots(dCoefs,nDegree,cRoots);
		}
	    break;
	case 4:
            cRoots[0].i = cRoots[1].i = cRoots[2].i = cRoots[3].i = 0.0; /* initialize for formula solvers */
	    if (i == 0)
		printf("Quartic formula\n");
	    iFlag = polyQuarticSolveLimited(dCoefs[3]/dCoefs[4],dCoefs[2]/dCoefs[4],dCoefs[1]/dCoefs[4],dCoefs[0]/dCoefs[4],&cRoots[0].r,&cRoots[1].r,&cRoots[2].r,&cRoots[3].r);
	    if (i == 0) {
		if (iFlag)
		    printf("   Unable to solve -- skipped\n");
		else
		    show_roots(dCoefs,nDegree,cRoots);
		}
	    break;
	    }
	
	}
    
#if (NUM_ITER > 1)
    printf("EXECUTION TIME = %g sec\n",getCPU() - t0);
    t0 = getCPU();
#endif
    
    for (i=0;i<NUM_ITER;i++) {
	if (i == 0)
	    printf("Laguerre's method (no polish)\n");
	laguerre_solve(dCoefs,nCoefs,cRoots,0);
	if (i == 0)
	    show_roots(dCoefs,nDegree,cRoots);
	}
    
#if (NUM_ITER > 1)
    printf("EXECUTION TIME = %g sec\n",getCPU() - t0);
    t0 = getCPU();
#endif
    
    for (i=0;i<NUM_ITER;i++) {
	if (i == 0)
	    printf("Laguerre's method (with polish)\n");
	laguerre_solve(dCoefs,nCoefs,cRoots,1);
	if (i == 0)
	    show_roots(dCoefs,nDegree,cRoots);
	}
    
#if (NUM_ITER > 1)
    printf("EXECUTION TIME = %g sec\n",getCPU() - t0);
#endif
    
    free((void *) cRoots);
    free((void *) dCoefs);
    
    return 0;
    }

#endif /* NUMREC */

#ifdef GSL

static void show_roots(const double dCoefs[], int nDegree, const double dRoots[])
{
    double dReal, dImag;
    int i, bFound;
    
    printf("Real roots:\n");
    bFound = 0;
    for (i = 0; i < nDegree; i++) {
        dReal = dRoots[2 * i];
        dImag = dRoots[2 * i + 1];
        if (fabs(dImag) <= EPS * fabs(dReal)) {
            printf("   %.6e (evaluates to %.6e)\n", dReal,
                   gsl_poly_eval(dCoefs, nDegree + 1, dReal));
            bFound = 1;
        }
    }
    if (bFound == 0)
        printf("   none\n");
    
    printf("Complex roots:\n");
    bFound = 0;
    for (i = 0; i < nDegree; i++) {
        dReal = dRoots[2 * i];
        dImag = dRoots[2 * i + 1];
        if (fabs(dImag) > EPS * fabs(dReal)) {
            printf("   %.6e + i(%.6e)\n", dReal, dImag);
            bFound = 1;
        }
    }
    if (bFound == 0)
        printf("   none\n");
}

int main(int argc, char *argv[])
{
    double *dCoefs, *dRoots;
    int i, nDegree, nCoefs, iFlag;
    
#if (NUM_ITER > 1)
    double t0;
#endif
    
    setbuf(stdout, (char *) NULL);
    
    if (argc < 2)
        usage(argv[0]);
    
    nDegree = atoi(argv[1]);
    
    if (nDegree < 1) {
        fprintf(stderr, "Polynomial degree must be positive integer\n");
        usage(argv[0]);
    }
    
    nCoefs = nDegree + 1;
    
    if (argc != nCoefs + 2) {
        fprintf(stderr, "Must supply %i coefficients for %i-degree polynomial\n", nCoefs, nDegree);
        usage(argv[0]);
    }
    
    dCoefs = (double *) malloc(nCoefs * sizeof(double));
    assert(dCoefs != NULL);
    
    for (i = 0; i < nCoefs; i++)
        dCoefs[i] = atof(argv[argc - i - 1]); /* reverse order */
    
    if (dCoefs[nCoefs - 1] == 0.0) {
        fprintf(stderr, "Leading coefficient cannot be zero (e.g., ax^2 + bx + c = 0, a != 0)\n");
        free((void *) dCoefs);
        usage(argv[0]);
    }
    
    dRoots = (double *) malloc(2 * nDegree * sizeof(double));
    assert(dRoots != NULL);
    
#if (NUM_ITER > 1)
    t0 = getCPU();
#endif
    
    for (i = 0; i < NUM_ITER; i++) {
        
        switch (nDegree) {
            case 1:
                /* trivial case */
                dRoots[0] = - dCoefs[0] / dCoefs[1];
                dRoots[1] = 0.0;
                if (i == 0)
                    show_roots(dCoefs, nDegree, dRoots);
                break;
            case 2:
                dRoots[1] = dRoots[3] = 0.0; /* zero imaginary parts for output */
                if (i == 0)
                    printf("Standard quadratic formula\n");
                iFlag = quad_solve_std(dCoefs[2], dCoefs[1], dCoefs[0], &dRoots[0], &dRoots[2]);
                if (i == 0) {
                    if (iFlag != 0)
                        printf("   No real roots\n");
                    else
                        show_roots(dCoefs, nDegree, dRoots);
                    printf("Quadratic solver\n");
                }
                iFlag = polyQuadSolve(dCoefs[2], dCoefs[1], dCoefs[0], &dRoots[0], &dRoots[2]);
                if (i == 0) {
                    if (iFlag != 0)
                        printf("   No real roots\n");
                    else
                        show_roots(dCoefs, nDegree, dRoots);
                }
                break;
            case 3:
                dRoots[1] = dRoots[3] = dRoots[5] = 0.0;
                if (i == 0)
                    printf("Cubic solver\n");
                 iFlag = polyCubicSolveLimited(dCoefs[2] / dCoefs[3], dCoefs[1] / dCoefs[3], dCoefs[0] /dCoefs[3], &dRoots[0], &dRoots[2], &dRoots[4]);
                if (i == 0) {
                    if (iFlag != 0)
                        printf("   1 real root, 2 complex roots -- skipped\n");
                    else
                        show_roots(dCoefs, nDegree, dRoots);
                }
                break;
            case 4:
                dRoots[1] = dRoots[3] = dRoots[5] = dRoots[7] = 0.0;
                if (i == 0)
                    printf("Quartic solver\n");
                iFlag = polyQuarticSolveLimited(dCoefs[3] / dCoefs[4], dCoefs[2] / dCoefs[4], dCoefs[1] / dCoefs[4], dCoefs[0] / dCoefs[4], &dRoots[0], &dRoots[2], &dRoots[4], &dRoots[6]);
                if (i == 0) {
                    if (iFlag != 0)
                        printf("   Unable to solve -- skipped\n");
                    else
                        show_roots(dCoefs, nDegree, dRoots);
                }
                break;
        }
        
    }
    
#if (NUM_ITER > 1)
    printf("EXECUTION TIME = %g sec\n", getCPU() - t0);
    t0 = getCPU();
#endif
    
    for (i = 0; i < NUM_ITER; i++) {
        if (i == 0)
            printf("General method\n");
        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(nCoefs);
        assert(w != NULL);
        gsl_poly_complex_solve(dCoefs, nCoefs, w, dRoots);
        gsl_poly_complex_workspace_free(w);
        if (i == 0)
            show_roots(dCoefs, nDegree, dRoots);
    }
    
#if (NUM_ITER > 1)
    printf("EXECUTION TIME = %g sec\n",getCPU() - t0);
#endif
    
    free((void *) dRoots);
    free((void *) dCoefs);
    
    return 0;
}

#endif /* GSL */

#endif /* STANDALONE */
