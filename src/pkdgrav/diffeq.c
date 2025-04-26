/*
** Routines for integrating ordinary differential equations.
** Includes a standalone test driver.
*/

/*#define STANDALONE*/ /* uncomment to compile standalone test */

#ifdef NUMREC
/* compile with -DNUMREC */
#endif

#ifdef GSL
/* compile with -DGSL -I/path/to/headers -L/path/to/lib -lgsl */
#endif

#include <stdio.h>
#include <stdlib.h>
#include "diffeq.h"

#ifdef NUMREC

#include <math.h>
#include <assert.h>

/*
** Following routines adapted from Numerical Recipes in C (2nd ed).
** Note these are revised versions -- see online book.
*/

static
void rkck(FLOAT y[], FLOAT dydx[], int n, FLOAT x, FLOAT h, FLOAT yout[],
	  FLOAT yerr[], diffeqDerivsT derivs, void *user_data)
/*
** Based on rkck(), NRiC(2e) 16.2.
**
** Given values for n variables y[1..n] and their derivatives
** dydx[1..n] known at x, use the fifth-order Cash-Karp Runge-Kutta
** method to advance the solution over an interval h and return the
** incremented variables as yout[1..n]. Also return an estimate of the
** local truncation error in yout using the embedded fourth-order
** method. The user supplies the routine derivs, which returns
** derivatives dydx at x (user_data is passed to derivs).
*/
{
    int i;
    static FLOAT a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
      b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
      b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
      b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
      b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
      c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
      dc5 = -277.00/14336.0;
    FLOAT dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
      dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    FLOAT *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
    ak2=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ak2 != NULL);
    ak3=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ak3 != NULL);
    ak4=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ak4 != NULL);
    ak5=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ak5 != NULL);
    ak6=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ak6 != NULL);
    ytemp=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ytemp != NULL);
    for (i=0;i<n;i++) /* First step. */
	ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2,user_data); /* Second step. */
    for (i=0;i<n;i++)
	ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3,user_data); /* Third step. */
    for (i=0;i<n;i++)
	ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4,user_data); /* Fourth step. */
    for (i=0;i<n;i++)
	ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5,user_data); /* Fifth step. */
    for (i=0;i<n;i++)
	ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6,user_data); /* Sixth step. */
    for (i=0;i<n;i++) /* Accumulate increments with proper weights. */
	yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=0;i<n;i++) /* Estimate error as difference between fourth and fifth order methods. */
	yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    free((void *) ytemp);
    free((void *) ak6);
    free((void *) ak5);
    free((void *) ak4);
    free((void *) ak3);
    free((void *) ak2);
    }

#define SAFETY 0.9 /* SAFETY, PGROW, and PSHRNK relate to Eq. (16.2.10) */
#define PGROW  (-0.2)
#define PSHRNK (-0.25)
#define ERRCON 1.89e-4 /* ERRCON is (5/SAFETY)^(1/PGROW); see use below */

FLOAT FMAX(FLOAT,FLOAT);
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
		   (maxarg1) : (maxarg2))

FLOAT FMIN(FLOAT,FLOAT);
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
		   (minarg1) : (minarg2))

void rkqs(FLOAT y[], FLOAT dydx[], int n, FLOAT *x, FLOAT htry, FLOAT eps,
	  FLOAT yscal[], FLOAT *hdid, FLOAT *hnext, diffeqDerivsT derivs,
	  void *user_data)
/*
** Based on rkqs(), NRiC(2e) 16.2.
**
** Fifth-order Runge-Kutta step with monitoring of local truncation
** error to ensure accuracy and adjust stepsize. Input are the
** dependent variable vector y[1..n] and its derivative dydx[1..n] at
** the starting value of the independent variable x. Also input are
** the stepsize to be attempted htry, the required accuracy eps, and
** the vector yscal[1..n] against which the error is scaled. On
** output, y and x are replaced by their new values, hdid is the
** stepsize that was actually accomplished, and hnext is the estimated
** next stepsize. derivs is the user-supplied routine that computes
** the right-hand side derivatives. user_data is passed through to
** derivs.
*/
{
    int i;
    FLOAT errmax,h,htemp,xnew,*yerr,*ytemp;
    FLOAT maxarg1,maxarg2,minarg1,minarg2; /* for FMAX() and FMIN() */
    yerr=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(yerr != NULL);
    ytemp=(FLOAT *) malloc(n*sizeof(FLOAT));
    assert(ytemp != NULL);
    h=htry; /* Set stepsize to the initial trial value. */
    for (;;) {
	rkck(y,dydx,n,*x,h,ytemp,yerr,derivs,user_data); /* Take a step. */
	errmax=0.0; /* Evaluate accuracy. */
	for (i=0;i<n;i++){
		errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));}
	errmax /= eps; /* Scale relative to required tolerance. */
	if (errmax <= 1.0) break; /* Step succeeded. Compute size of next step. */
	htemp=SAFETY*h*pow(errmax,PSHRNK);
	/* Truncation error too large, reduce stepsize. */
	h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h)); /* No more than a factor of 10. */
	xnew=(*x)+h;
	if (xnew == *x) {
	    fprintf(stderr,"rkqs(): stepsize underflow\n");
	    exit(1);
	    }
	}
    if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
    else *hnext=5.0*h; /* No more than a factor of 5 increase. */
    *x += (*hdid=h);
    for (i=0;i<n;i++) y[i]=ytemp[i];
    free((void *) ytemp);
    free((void *) yerr);
    }

#undef FMIN
#undef FMAX
#undef ERRCON
#undef PSHRNK
#undef PGROW
#undef SAFETY

#define MAXSTP 10000
#define TINY 1.0e-15 /* was 1.0e-30 (else get step underflow with rigid aggs) */

FLOAT SIGN(FLOAT,FLOAT);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* code for storage of intermediate results omitted */

void odeint(FLOAT ystart[], int nvar, FLOAT x1, FLOAT x2, FLOAT eps, FLOAT h1,
            FLOAT hmin, int *nok, int *nbad, diffeqDerivsT derivs, void *user_data,
            void (*stepper)(FLOAT [], FLOAT [], int, FLOAT *, FLOAT, FLOAT, FLOAT [],
                            FLOAT *, FLOAT *, diffeqDerivsT, void *))
/*
** Based on odeint(), NRiC(2e) 16.2.
**
** Runge-Kutta driver with adaptive stepsize control. Integrate
** starting values ystart[1..nvar] from x1 to x2 with accuracy eps. h1
** should be set as a guessed first stepsize, hmin as the minimum
** desired stepsize (can be zero). On output nok and nbad are the
** number of good and bad (but retried and fixed) steps taken, and
** ystart is replaced by values at the end of the integration
** interval. derivs is the user-supplied routine for calculating the
** right-hand side derivative, while stepper is the name of the
** stepper routine to be used. user_data is passed 2through to derivs.
*/
{
    int nstp,i;
    FLOAT x,hnext,hdid,h;
    FLOAT *yscal,*y,*dydx;
    yscal=(FLOAT *) malloc(nvar*sizeof(FLOAT));
    assert(yscal != NULL);
    y=(FLOAT *) malloc(nvar*sizeof(FLOAT));
    assert(y != NULL);
    dydx=(FLOAT *) malloc(nvar*sizeof(FLOAT));
    assert(dydx != NULL);
    x=x1;
    h=SIGN(h1,x2-x1);
    *nok = *nbad = 0;
    for (i=0;i<nvar;i++) y[i]=ystart[i];
    for (nstp=0;nstp<MAXSTP;nstp++) { /* Take at most MAXSTP steps. */
	(*derivs)(x,y,dydx,user_data);
	for (i=0;i<nvar;i++)
	    /* Scaling used to monitor accuracy. This general-purpose
	       choice can be modified if need be. */
	    yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
	if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x; /* If stepsize can overshoot, decrease. */
	(*stepper)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,user_data);
	if (hdid == h) ++(*nok); else ++(*nbad);
	if ((x-x2)*(x2-x1) >= 0.0) { /* Are we done? */
	    for (i=0;i<nvar;i++) ystart[i]=y[i];
	    free((void *) dydx);
	    free((void *) y);
	    free((void *) yscal);
	    return; /* Normal exit. */
	    }
	if (fabs(hnext) <= hmin) {
	    fprintf(stderr,"WARNING: odeint(): small step size (%g)\n",fabs(hnext));
	    /*exit(1);*/
	    }
	h=hnext;
	assert(fabs(h) > 0.0);
	}
    fprintf(stderr,"odeint(): too many steps\n");
    exit(1);
    }

#undef SIGN
#undef TINY
#undef MAXSTP

void diffeqIntegrate(FLOAT fVars[], int nVars, FLOAT fStart, FLOAT fStop,
                     FLOAT fAccuracy, FLOAT fStepTry, FLOAT fStepMin,
                     diffeqDerivsT funcDerivs, void *pUserData) {
    int nok, nbad; /* ignored */
    odeint(fVars, nVars, fStart, fStop, fAccuracy, fStepTry, fStepMin,
           &nok, &nbad, funcDerivs, pUserData, rkqs);
}

#endif /* NUMREC */

#ifdef GSL

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#define MAX_GSL_ODEIV2_STEPS 1000

void diffeqIntegrate(FLOAT fVars[], int nVars, FLOAT fStart, FLOAT fStop,
                     FLOAT fAccuracy, FLOAT fStepTry, FLOAT fStepMin,
                     diffeqDerivsT funcDerivs, void *pUserData) {
    int iRet;
    gsl_odeiv2_system sys = {funcDerivs, NULL, nVars, pUserData};
    gsl_odeiv2_driver *d =
        gsl_odeiv2_driver_alloc_standard_new(&sys, gsl_odeiv2_step_rkck,
                                             fStepTry, fAccuracy, fAccuracy,
					     1.0, 1.0);
    gsl_odeiv2_driver_set_hmin(d, fStepMin); /* use 0 for no minimum */
    gsl_odeiv2_driver_set_hmax(d, fStepTry); /* default is GSL_DBL_MAX */
    gsl_odeiv2_driver_set_nmax(d, MAX_GSL_ODEIV2_STEPS); /* 0 for no limit */
    iRet = gsl_odeiv2_driver_apply(d, &fStart, fStop, fVars);
    if (iRet != GSL_SUCCESS) {
        fprintf(stderr, "diffeqIntegrate(): driver error code %i\n", iRet);
        exit(1);
    }
    gsl_odeiv2_driver_free(d);
}

#undef MAX_GSL_ODEIV2_STEPS

#endif /* GSL */

#ifdef STANDALONE

#include <math.h>

int test_derivs(FLOAT fDum, const FLOAT fVars[], FLOAT fDerivs[], void *vDum)
{
    /* simple harmonic oscillator, amplitude 1, period 2 pi, phase 0 */

    fDerivs[0] = fVars[1];
    fDerivs[1] = - fVars[0];
    return 0;
    }

int main(void)
{
    FLOAT fVars[2], t, dt;
    int i, nSteps = 100;

    fVars[0] = 0.0;
    fVars[1] = 1.0;

    dt = 10.0 / nSteps;
    for (i = 0; i < nSteps; i++) {
        t = i * dt;
        if (i == 0)
            printf("Initial conditions (t = %g): x = %g, v = %g\n",
                   t, fVars[0], fVars[1]);
        diffeqIntegrate(fVars, 2, t, t + dt, 1.0e-9, dt, 0.0, test_derivs, NULL);
        printf("Step %2i (t = %4.1f): x, v = %6.3f, %6.3f (expect %6.3f, %6.3f)\n",
               i, t + dt, fVars[0], fVars[1], sin(t + dt), cos(t + dt));
	}

    return 0;
    }

#endif /* STANDALONE */
