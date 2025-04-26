/* suppress pedantic "empty translation unit" warning if !DEM */
typedef int make_iso_compilers_happy;

#ifdef DEM

#include <math.h> /* for log10() */
#include "collision.h" /* includes pkd.h and so dem.h */
#include "linalg.h"
#include "dem.h"

#define MIN_LOG_OVERLAP (-6.0)
#define MAX_LOG_OVERLAP ( 0.0)

#define MIN_COS_A (-1.0)
#define MAX_COS_A ( 1.0)

#define MIN_LOG_S (-5.0)
#define MAX_LOG_S ( 1.0)

#ifndef SQ
double SQ(double);
#define SQ(x) ((x)*(x))
#endif

void demWipeOverlapElement(DEM_ELEMENT *e) {
	e->iOrder = -1;
	vectorZero(e->vShear);
	vectorZero(e->vnOld);
#ifdef DEM_ROTATION_DASHPOT
        vectorZero(e->vRoll);
        e->dTwist = 0.0;
#endif /* DEM_ROTATION_DASHPOT */
#ifdef DEM_DIAG
                vectorZero(e->vFn);
                vectorZero(e->vFt);
#endif /* DEM_DIAG */
	e->liOverlapCounter = 0;
	}

void demAssignParticleInertia(int p1stuck,int p2stuck,double mass1,double mass2,double rad1,double rad2,double rad1sq,double rad2sq,
							  double *mass_inv1,double *mass_inv2,double *mom1,double *mom2,double *mom_inv1,double *mom_inv2) {

	/* inertia multipliers for application of forces and torques based on linear and angular accelerations */

	if (!p1stuck) {
		*mass_inv1 = 1./mass1;
		*mom1 = 0.4*mass1*rad1sq;
		*mom_inv1 = 1./(*mom1);
		}

	if (!p2stuck) {
		if (mass1 == mass2 && !p1stuck) { /* save time if masses are equal */
			*mass_inv2 = (*mass_inv1);
			if (rad1 == rad2) { /* save time if both masses and radii are equal */
				*mom2 = (*mom1);
				*mom_inv2 = (*mom_inv1);
				}
			else {
				*mom2 = 0.4*mass2*rad2sq;
				*mom_inv2 = 1./(*mom2);
				}
			}
		else {
			*mass_inv2 = 1./mass2;
			*mom2 = 0.4*mass2*rad2sq;
			*mom_inv2 = 1./(*mom2);
			}
		}
	}

int demCheckOverlapPoint(const Vector vRelPos,double dRadSq) {
	return vectorMagSq(vRelPos) <= dRadSq;
	}

#ifdef DEM_TRACK_ORIENT

#define EPS (1.0e-6) /* for orthogonality check */

void demUpdateOrient(PARTICLE *p, double dt)
{
	/*
	** Updates principal axis orientations of particle p over time
	** interval dt. This function is called during the pkdgrav
	** drift step, so dt = dDelta.  Because the particles are
	** spheres, the equations are decoupled and can be solved
	** during the regular pkdgrav leapfrog (kicks to particle
	** spins are unchanged compared to the procedure without
	** orientation tracking).  See Richardson (1995) for the
	** equations, which simplify here because the inertia tensor
	** is diagonal with equal elements.
	*/

	Vector v;
	Matrix R, O;
	double w;

	w = vectorMag(p->w);

	if (w == 0.0)
		return;

	vectorCopy(p->w, v);
	vectorNorm(v); /* normalized rotation axis */

	matrixQuatRot(v, w * dt, R); /* get rotation matrix */

	matrixCopy(p->orient, O);
	matrixMultiply(O, R, p->orient); /* rotate the principal axes */

#ifdef OBSOLETE
	double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z;
	Vector vSpin;
	Matrix mSpaceToBody;

	/* principal axes are stored as the columns of the orient matrix */

	matrixTranspose(p->orient, mSpaceToBody); /* columns to rows */
	matrixTransform(mSpaceToBody, p->w, vSpin); /* get spin in body frame */

	p1x = p->orient[0][0] + (vSpin[2]*p->orient[0][1] - vSpin[1]*p->orient[0][2])*dt;
	p1y = p->orient[1][0] + (vSpin[2]*p->orient[1][1] - vSpin[1]*p->orient[1][2])*dt;
	p1z = p->orient[2][0] + (vSpin[2]*p->orient[2][1] - vSpin[1]*p->orient[2][2])*dt;
	p2x = p->orient[0][1] + (vSpin[0]*p->orient[0][2] - vSpin[2]*p->orient[0][0])*dt;
	p2y = p->orient[1][1] + (vSpin[0]*p->orient[1][2] - vSpin[2]*p->orient[1][0])*dt;
	p2z = p->orient[2][1] + (vSpin[0]*p->orient[2][2] - vSpin[2]*p->orient[2][0])*dt;
	p3x = p->orient[0][2] + (vSpin[1]*p->orient[0][0] - vSpin[0]*p->orient[0][1])*dt;
	p3y = p->orient[1][2] + (vSpin[1]*p->orient[1][0] - vSpin[0]*p->orient[1][1])*dt;
	p3z = p->orient[2][2] + (vSpin[1]*p->orient[2][0] - vSpin[0]*p->orient[2][1])*dt;

	vectorSet(p->orient[0], p1x, p2x, p3x);
	vectorSet(p->orient[1], p1y, p2y, p3y);
	vectorSet(p->orient[2], p1z, p2z, p3z);
	/*
	** "orient" should be an orthonormal matrix, but the ODE
	** solution above will not necessarily guarantee this.  To
	** mitigate, we use a Gram-Schmidt orthonormalization process.
	*/

	matrixOrthonormalize(p->orient, 0);

	/*
	** renormalize the columns of the matrix then
	** assert that they are orthogonal within a tolerance.  (Also
	** see aggsRungeAdvance - previously in msrAggsAdvance().)
	*/

	/*

	vectorNorm(p->orient[0]);
	vectorNorm(p->orient[1]);
	vectorNorm(p->orient[2]);

	assert(fabs(vectorDot(p->orient[0], p->orient[1])) < EPS);
	assert(fabs(vectorDot(p->orient[0], p->orient[2])) < EPS);
	assert(fabs(vectorDot(p->orient[1], p->orient[2])) < EPS);
	*/
#endif /* OBSOLETE */

	}

#undef EPS

#endif /* DEM_TRACK_ORIENT */

#ifdef DEM_WALLS_REACT

void pkdDEMWallsReact(PKD pkd,double dTotalXForceFromParticles[],double dTotalYForceFromParticles[],double dTotalZForceFromParticles[], double dTotalXTorqueFromParticles[], double dTotalYTorqueFromParticles[], double dTotalZTorqueFromParticles[]) {
	PARTICLE *p;
	int i,j,nLocal = pkdLocal(pkd);
	for (i=0;i<MAX_NUM_WALL_ASSEMBLIES;i++) {
		dTotalXForceFromParticles[i] = 0.0;
		dTotalYForceFromParticles[i] = 0.0;
		dTotalZForceFromParticles[i] = 0.0; //structure definition in pst.h

		dTotalXTorqueFromParticles[i] = 0.0;
		dTotalYTorqueFromParticles[i] = 0.0;
		dTotalZTorqueFromParticles[i] = 0.0;
		}

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		for (j=0;j<MAX_NUM_WALL_ASSEMBLIES;j++) {
			dTotalXForceFromParticles[j] += p->dXForceOnWalls[j]; //dZForceOnWalls calculated in smoothfcn.c
			dTotalYForceFromParticles[j] += p->dYForceOnWalls[j];
			dTotalZForceFromParticles[j] += p->dZForceOnWalls[j];

			dTotalXTorqueFromParticles[j] += p->dXTorqueOnWalls[j]; //dXTorqueOnWalls calculated in smoothfcn.c
			dTotalYTorqueFromParticles[j] += p->dYTorqueOnWalls[j];
			dTotalZTorqueFromParticles[j] += p->dZTorqueOnWalls[j];
			}
		}
	}

#endif /* DEM_WALLS_REACT */

void pkdDEMZeroSmallMotions(PKD pkd,double dAccCritSq,double dDeltaSq)
{
	/*
	** Zeroes accelerations and velocities (and torques per unit mass
	** and spins) for particles with these quantities below a
	** user-specified threshold.  The strategy is that the
	** acceleration quantities must be oppositely directed to the
	** velocity quantities, and both must have magnitudes below the
	** threshold (a critical acceleration, e.g., a small fraction of
	** the ambient gravity) to be zeroed.  Note: at this point p->v is
	** the half-step velocity, but that's ok because if we zero both
	** p->a and p->v, then the full-step velocity will also be zero.
	** The same is true for p->w.  This routine should be called after
	** DEM forces have been fully calculated and before the closing
	** velocity kick.
	*/

	PARTICLE *p;
	double r2,dp,a2,v2;
	int i,n;

	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		r2 = RADIUS(p)*RADIUS(p);
		assert(r2 > 0.0);
		dp = p->a[0]*p->v[0] + p->a[1]*p->v[1] + p->a[2]*p->v[2]; /* dv/dt dot v */
		if (dp < 0.0) {
			a2 = p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2];
			if (a2 > 0.0 && a2 < dAccCritSq) {
				v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
				if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq) {
					p->v[0] = p->v[1] = p->v[2] = 0.0;
					p->a[0] = p->a[1] = p->a[2] = 0.0;
					}
				}
			}
		dp = p->wDot[0]*p->w[0] + p->wDot[1]*p->w[1] + p->wDot[2]*p->w[2]; /* dw/dt dot w */
		if (dp < 0.0) {
			a2 = p->wDot[0]*p->wDot[0] + p->wDot[1]*p->wDot[1] + p->wDot[2]*p->wDot[2];
			if (a2 > 0.0 && a2 < dAccCritSq*r2) {
				v2 = p->w[0]*p->w[0] + p->w[1]*p->w[1] + p->w[2]*p->w[2];
				if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq/r2) {
					p->w[0] = p->w[1] = p->w[2] = 0.0;
					p->wDot[0] = p->wDot[1] = p->wDot[2] = 0.0;
					}
				}
			}
		}
	}

void pkdDEMStats(PKD pkd,DEM_STATS *ds) {
	PARTICLE *p;
	double dWidth,d;
	int i,nLocal = pkdLocal(pkd);

	/* find minimum particle separation & maximum relative overlap speed on this processor */

	ds->fDistMin = FLOAT_MAXVAL;
	ds->fSpeedMax = 0.0;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fDistMin < ds->fDistMin)
			ds->fDistMin = p->DEMStats.fDistMin;
		if (p->DEMStats.fSpeedMax > ds->fSpeedMax)
			ds->fSpeedMax = p->DEMStats.fSpeedMax;
		}

	assert(ds->fDistMin >= 0.0); /* sanity check */
	assert(ds->fSpeedMax >= 0.0);

	/* initialize histogram bins */

	assert(DEM_NUM_OVERLAP_BINS > 0);

	for (i=0;i<DEM_NUM_OVERLAP_BINS;i++)
		ds->pOverlapHist[i] = 0;

	dWidth = (MAX_LOG_OVERLAP - MIN_LOG_OVERLAP)/DEM_NUM_OVERLAP_BINS;

	/* fill overlap histogram with local particle data */

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fOverlap > 0.0) {
			d = log10(p->DEMStats.fOverlap);
			if (d < MIN_LOG_OVERLAP)
				++ds->pOverlapHist[0];
			else if (d >= MAX_LOG_OVERLAP) /* > MAX only possible if small particle fully inside large particle */
				++ds->pOverlapHist[DEM_NUM_OVERLAP_BINS - 1];
			else
				++ds->pOverlapHist[(int) ((d - MIN_LOG_OVERLAP)/dWidth)];
			}
		}

	/* repeat for cos(alpha) histogram */

	assert(DEM_NUM_COS_A_BINS > 0);

	for (i=0;i<DEM_NUM_COS_A_BINS;i++)
		ds->pCosAHist[i] = 0;

	dWidth = (MAX_COS_A - MIN_COS_A)/DEM_NUM_COS_A_BINS;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		d = p->DEMStats.fCosA;
		if (d >= MIN_COS_A && d <= MAX_COS_A) {
			if (d == MAX_COS_A)
				++ds->pCosAHist[DEM_NUM_COS_A_BINS - 1];
			else
				++ds->pCosAHist[(int) ((d - MIN_COS_A)/dWidth)];
			}
		}

	/* ditto for S vectors */

	assert(DEM_NUM_S_BINS > 0);

	for (i=0;i<DEM_NUM_S_BINS;i++)
		ds->pSHist[i] = 0;

	dWidth = (MAX_LOG_S - MIN_LOG_S)/DEM_NUM_S_BINS;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->DEMStats.fS2 > 0.0) {
			d = 0.5*log10(p->DEMStats.fS2/SQ(RADIUS(p)));
			if (d < MIN_LOG_S)
				++ds->pSHist[0];
			else if (d >= MAX_LOG_S)
				++ds->pSHist[DEM_NUM_S_BINS - 1];
			else
				++ds->pSHist[(int) ((d - MIN_LOG_S)/dWidth)];
			}
		}
	}

#ifdef DEM_FIXED_BALL

void pkdDEMSetBall(PKD pkd, double dFixedRadius, int iFixedBallOption)
{
	static double epsilon = 0.1;

	PARTICLE *p;
	double x,r;
	int i,k,nLocal;

	nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		r = RADIUS(p); /* special case */
		switch (iFixedBallOption) {
		case 0:
			x = (2.0 + epsilon)*dFixedRadius;
			break;
		case 1:
			x = r + (1.0 + epsilon)*dFixedRadius;
			break;
		case 2:
			x = r + (1.0 + epsilon)*p->dNextLargestRad;
			break;
		default:
			assert(0); /* shouldn't get here */
			}
		/* make sure ball isn't bigger than periodic cell */
		for (k=0;k<3;k++)
			assert(x < 0.5*pkd->fPeriod[k]);
		p->fBallMax = x; /* for expanding tree node size */
		p->fBall2 = x*x; /* for searching around particle */
		p->cpStart = 0; /* force smooth search to start at tree root */
		}
	}

void pkdDEMGetBallRadius(PKD pkd, int iFixedBallOption, double *dFixedRadius)
{
	PARTICLE *p;
	double r;
	int i,nLocal;

	*dFixedRadius = 0.0;
	nLocal = pkdLocal(pkd);
	
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		r = RADIUS(p);
		if (iFixedBallOption == 0) {
			if (r > *dFixedRadius)
				*dFixedRadius = r;
			}
		else if (iFixedBallOption == 1)
			*dFixedRadius += r;
		else
			assert(0);
		}
	}

void pkdDEMFixedBallNextLargest(PKD pkd)
{
	PARTICLE *p;
	PARTICLE *p_next = NULL; /* initialized to suppress compiler warning */
	double r = 0.0; /* ditto */
	int i,nLocal;
	nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal-1;i++) {
		p = &pkd->pStore[i];
		p_next = &pkd->pStore[i+1];
		r = RADIUS(p_next);
		p->dNextLargestRad = r;
		}
	p_next->dNextLargestRad = r;
}

#endif /* DEM_FIXED_BALL */

#endif /* DEM */
