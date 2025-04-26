#ifndef POLYROOTS_HINCLUDED
#define POLYROOTS_HINCLUDED

int polyQuadSolve(double a, double b, double c, double *x1, double *x2);
int polyCubicSolveLimited(double a1, double a2, double a3, double *x1,
			  double *x2, double *x3);
int polyQuarticSolveLimited(double a1, double a2, double a3, double a4,
			    double *x1, double *x2, double *x3, double *x4);
void polyFindRealRoots(int nDegree, double dCoefs[], double dRoots[],
		       int *nRealRoots);

#endif /* !POLYROOTS_HINCLUDED */
