/* suppress pedantic "empty translation unit" warning if !WALLS */
typedef int make_iso_compilers_happy;

#ifdef WALLS

#include <math.h>
#include "pkd.h"
#include "collision.h" /* includes walls.h */
#include "polyroots.h"

/*#define DEBUG_WALLS*/ /*DEBUG to assist with debugging -- comment out for production runs! */

#ifdef DEBUG_WALLS
#define TRACE_PART (-1) /* set to -1 to disable */
int DEBUG_TRACE(int iOrder);
#define DEBUG_TRACE(iOrder) ((iOrder) == TRACE_PART)
#endif

/* time magnitudes less than this times dDelta are taken to be ZERO */

#ifdef TINY
#undef TINY
#endif

#define TINY (1.0e-6) /* what should this be?... */

/* following so we can use TRUE/FALSE with confidence... */

#ifdef TRUE
#undef TRUE
#endif

#define TRUE 1

#ifdef FALSE
#undef FALSE
#endif

#define FALSE 0

const int nWallFaces[] = { /* cf. walls.h */
  2/*WallPlane*/,
  8/*WallTriangle*/,
  10/*WallRectangle*/,
  4/*WallDisk*/,
  2/*WallCylinderInfinite*/,
  4/*WallCylinderFinite*/,
  3/*WallShell*/
};

void wallsRotate(const WALL_PARAMS *WP,int iWallID,Vector vPos,Vector vVel,Vector vOmegaV,double dt)
{
	const WALL_DATA *w;
	Vector vRelPos,vTmp,vRelPosOld;
	double dRdotZ;

	assert(iWallID >= 0 && iWallID < WP->nWalls);
	assert(iWallID == WP->pWalls[iWallID].iWallID);
	w = &WP->pWalls[iWallID].wd;

	if (w->dAngSpeed == 0.0)
		return;

	assert(w->iType != WallTriangle && w->iType != WallRectangle);

	vectorSub(vPos,w->vOrigin,vRelPos); /* vector from wall origin to particle center */
	dRdotZ = vectorDot(vRelPos,w->vOrient); /* component along rotation axis */
	vectorScale(w->vOrient,dRdotZ,vTmp); /* vector along axis */
	vectorSub(vRelPos,vTmp,vRelPos); /* vector perpendicular to axis pointing at particle */
	vectorCopy(vRelPos,vRelPosOld); /* save a pre-rotation copy */
	vectorRotate(vRelPos,w->vOrient,w->dAngSpeed*dt); /* rotate vector by requisite amount */
	vectorSub(vRelPos,vRelPosOld,vTmp); /* vector change in particle position */
	vectorAdd(vPos,vTmp,vPos); /* update particle position */
	vectorScale(w->vOrient,w->dAngSpeed,vTmp); /* angular velocity vector */
	vectorCross(vTmp,vRelPos,vVel); /* update particle velocity (due to rotation) */
	vectorCross(vTmp,vVel,vOmegaV); /* centripetal (2nd-order) term */
	}

void pkdWallsUpdateStuckParticles(PKD pkd,const WALL_PARAMS *WP,int bUpdatePos,double dDelta)
{
	/* note dDelta is only used for rotation */

	PARTICLE *p;
	const WALL *w;
	int i,iWallID,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (PARTICLE_STUCK(p)) {
			iWallID = PARTICLE_WALL_ID(p);
			assert(iWallID >= 0 && iWallID < WP->nWalls);
			w = &WP->pWalls[iWallID];
			assert(iWallID == w->iWallID);
			if (bUpdatePos)
				vectorAdd(p->r,w->vTravel,p->r);
			vectorCopy(w->vTotVel,p->v);
			if (PARTICLE_STUCK_ON_ROTATING_WALL(p,WP)) {
				vectorScale(w->wd.vOrient,w->wd.dAngSpeed,p->w);
				if (dDelta != 0.0)
					wallsRotate(WP,iWallID,p->r,p->v,p->omega_v,dDelta);
				}
			else
				vectorZero(p->w);
#ifdef NEED_VPRED
			vectorCopy(p->v,p->vPred);
			vectorCopy(p->w,p->wPred);
#endif
			}
		}
	}

static double dsgn(double x)
{
	return (x < 0.0 ? -1.0 : 1.0); /* zero taken as positive */
	}

static void wallGetRelPosAndVel(const WALL_DATA *w,const Vector vOscVel,const Vector r,const Vector v,double dt,Vector vRelPos,Vector vRelVel)
{
	/*
	** Returns position and velocity of particle relative to wall
	** origin, accounting for linear and oscillatory motion during
	** step.
	*/

	Vector vTmp;

	vectorSub(r,w->vOrigin,vRelPos);
	vectorScale(w->vVel,dt,vTmp); /* displacement vector due to linear motion */
	vectorSub(vRelPos,vTmp,vRelPos);
	vectorSub(v,w->vVel,vRelVel);
	if (w->dOscAmp != 0.0) { /* oscillatory component -- for simplicity we assume constant velocity during step */
		vectorScale(vOscVel,dt,vTmp);
		vectorSub(vRelPos,vTmp,vRelPos);
		vectorSub(vRelVel,vOscVel,vRelVel);
		}
	}

static void wallGetRectangleGeometry(const WALL_DATA *w,double *dWidth,double *dLength,Matrix m)
{
	*dWidth = vectorMag(w->vVertex1);
	assert(*dWidth > 0.0);
	*dLength = vectorMag(w->vVertex2);
	assert(*dLength > 0.0);

	/* m transforms space frame to body frame for rectangle */

	/* body frame: vertex 1 is on +ve x-axis; 2 on +ve y-axis */

	vectorCopy(w->vVertex1,m[0]); /* 1st row */
	vectorScale(m[0],1.0/(*dWidth),m[0]); /* normalize */
	vectorCopy(w->vVertex2,m[1]); /* 2nd row */
	vectorScale(m[1],1.0/(*dLength),m[1]); /* normalize */
	vectorCopy(w->vOrient,m[2]); /* 3rd row, already normalized (cross product of first 2 rows) */
	}

static int wallParticleIntersectsRectangle(const Vector vPos,double dRadius,const WALL_DATA *w)
{
	/* assumes particle already intersecting/touching plane containing rectangle (i.e. dR2 >= 0.0 using definition below) -- save times and avoids roundoff problems! */

	/* does NOT consider touching sides */

	/* vPos is relative to rectangle origin (lower-left corner) */

	Matrix m;
	Vector v;
	double dWidth,dLength,dx,dy,x2,y2,dx2,dy2,dR2;

	wallGetRectangleGeometry(w,&dWidth,&dLength,m);

	matrixTransform(m,vPos,v); /* get position in body frame */

	/* handy abbreviations */

	dx = v[0] - dWidth;
	dy = v[1] - dLength;

	/* quick test: if projected position inside rectangle, have intersection (or face touch) */

	if (v[0] > 0.0 && dx < 0.0 && v[1] > 0.0 && dy < 0.0)
		return 1;

	/* to proceed further, particle must be penetrating plane */

	if (fabs(v[2]) >= dRadius)
		return 0;

	/* if projected position far enough outside rectangle, no intersection */

	if (v[0] <= -dRadius || dx >= dRadius || v[1] <= -dRadius || dy >= dRadius)
		return 0;

	/* more handy abbreviations */

	x2 = v[0]*v[0];
	y2 = v[1]*v[1];
	dx2 = dx*dx;
	dy2 = dy*dy;

	dR2 = dRadius*dRadius - v[2]*v[2]; /* i.e. R^2 - z^2 */

	/* check sides */

	if ((v[0] > 0.0 && dx < 0.0 && (y2 < dR2 || dy2 < dR2)) ||
		(v[1] > 0.0 && dy < 0.0 && (x2 < dR2 || dx2 < dR2)))
		return 1;

	/* check corners */

	if (x2 + y2 < dR2 || x2 + dy2 < dR2 || dx2 + y2 < dR2 || dx2 + dy2 < dR2)
		return 1;

	/* no intersection */

	return 0;
	}

static int wallParticleIntersectsDisk(Vector vPos,double dRadPart,const WALL_DATA *w,double rdotn)
{
	/* vPos is relative to disk origin */

	Vector v,vProj;
	double rProj;

	vectorScale(w->vOrient,rdotn,v); /* vector from origin to particle projected along symmetry axis */
	vectorSub(vPos,v,vProj); /* vector from origin to particle projected on disk plane */
	rProj = vectorMag(vProj); /* distance from origin (projected on plane) */

	/*
	** Quick test: if projected position is inside disk (accounting
	** for any inner hole), have intersection (or face touch).
	*/
	if (rProj < w->dRadius && (w->dHoleRadius == 0.0 || rProj > w->dHoleRadius))
		return 1;

	/* to proceed further, particle must be penetrating plane */

	if (fabs(rdotn) >= dRadPart)
		return 0;

	/* if projected position is far enough outside disk, no intersection */

	if (rProj >= w->dRadius + dRadPart)
		return 0;

	/* if projected position is far enough inside any disk hole, no intersection */

	if (w->dHoleRadius > 0.0 && rProj <= w->dHoleRadius - dRadPart)
		return 0;

	/* have intersection by now if projected distance is zero */

	if (rProj == 0.0)
		return 1;

	/* otherwise check intersection with perimeter(s) */

	vectorScale(vProj,w->dRadius/rProj,v);
	vectorSub(vPos,v,v); /* vector from particle center to nearest point on outer edge of disk */

	if (vectorMagSq(v) < dRadPart*dRadPart)
		return 1;

	if (w->dHoleRadius > 0.0) {
		vectorScale(vProj,w->dHoleRadius/rProj,v);
		vectorSub(vPos,v,v); /* vector from particle center to nearest point on inner edge of disk (if applicable) */

		if (vectorMagSq(v) < dRadPart*dRadPart)
			return 1;
		}

	return 0;
	}

/*
** Following functions determine possible intersect times for various
** wall types.  General strategy: if r and v are position and velocity
** of particle relative to wall origin, seek dt such that r-prime = r
** + v dt is in contact with the wall.  The contact condition varies
** for each geometry.  Wall and particle data are input; face minimum
** impact times and a flag whether the particle is overlapping the
** wall are output.  Wall ID is passed for debugging purposes only.
*/

static void wallIntersectFlatFace(int iWallID,const WALL_DATA *w,const Vector vOscVel,const PARTICLE *p,double dStart,double dTmax,double *dtUpper,double *dtLower,int *bOverlap)
{
	/*
	** Returns times to impact either face of flat surface (any
	** perimeter handled as separate intersection test by calling
	** function).  Strategy: set |r-prime dot n-hat| = RADIUS(p).
	** Note for this and subsequent intersect functions: dStart is
	** needed in order to determine the exact location of any moving
	** walls, and dTmax is used when deciding whether a collision
	** interval is short enough to be considered to be zero (cf. TINY
	** macro).
	*/

	/* following definitions work for all flat geometries */

	const int nFaces = 2; /* upper & lower */
	enum {Upper=0,Lower};

	Vector vRelPos,vRelVel;
	double rdotn,vdotn,*dt;
	int i;

	assert(w->iType == WallPlane || w->iType == WallRectangle ||
		   w->iType == WallTriangle || w->iType == WallDisk);

	*dtUpper = *dtLower = DBL_MAX; /* DBL_MAX indicates no collision */
	*bOverlap = FALSE;

	/* get position and velocity relative to wall origin, accounting for any wall motion */

	wallGetRelPosAndVel(w,vOscVel,p->r,p->v,dStart,vRelPos,vRelVel);

	rdotn = vectorDot(vRelPos,w->vOrient); /* height above or below disk */
	vdotn = vectorDot(vRelVel,w->vOrient); /* speed toward or away from disk */

	/* check easy return conditions BEFORE checking overlap */

	if (w->iType == WallPlane && vdotn == 0.0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (type %i): no motion relative to (infinite) plane\n",p->iOrder,iWallID,w->iType);
#endif
		return; /* particle moving parallel to wall (or stationary) -- no collision */
		}

	if (rdotn*vdotn > 0.0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (type %i): moving apart\n",p->iOrder,iWallID,w->iType);
#endif
		return; /* particle moving away from wall -- no collision */
		}

	/* check for particle overlap with wall */

	if (fabs(rdotn) <= RADIUS(p)) {
		/* for finite geometries, check for actual intersection */
		switch (w->iType) {
		case WallPlane:
			*bOverlap = TRUE;
			break;
		case WallTriangle:
			assert(0); /* not ready yet */
		case WallRectangle:
			*bOverlap = wallParticleIntersectsRectangle(vRelPos,RADIUS(p),w);
			break;
		case WallDisk:
			*bOverlap = wallParticleIntersectsDisk(vRelPos,RADIUS(p),w,rdotn);
		    break;
		default:
			assert(0);
			}
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder) && *bOverlap)
			fprintf(stderr,"%i & wall %i (type %i) overlapping/touching: rdotn %g R %g vdotn %g\n",p->iOrder,iWallID,w->iType,rdotn,RADIUS(p),vdotn);
#endif
		}

	/* abort now if impact not possible */

	if (vdotn == 0.0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (type %i): no face impact possible (vdotn = 0)\n",p->iOrder,iWallID,w->iType);
#endif
		return;
		}

	/* consider each face in turn (doing it this way avoids having to solve a quadratic) */

	for (i=0;i<nFaces;i++) {
		switch (i) {
		case Upper:
			*(dt = dtUpper) = (RADIUS(p) - rdotn)/vdotn; /* time to impact */
			break;
		case Lower:
			*(dt = dtLower) = - (RADIUS(p) + rdotn)/vdotn;
			break;
		default:
			assert(0);
			}
		if (*dt < TINY*dTmax) /* handle possible roundoff error */
			*dt = 0.0;
		/* for finite geometries, check for actual contact (should just be touching plane) */
		if (w->iType == WallTriangle || w->iType == WallRectangle || w->iType == WallDisk) {
			Vector vNewPos;
			vectorScale(vRelVel,*dt,vNewPos); /* relative displacement */
			vectorAdd(vNewPos,vRelPos,vNewPos); /* new relative position */
			switch (w->iType) {
			case WallTriangle:
				assert(0);
			case WallRectangle:
				if (!wallParticleIntersectsRectangle(vNewPos,RADIUS(p),w))
					*dt = DBL_MAX;
				break;
			case WallDisk:
				if (!wallParticleIntersectsDisk(vNewPos,RADIUS(p),w,vectorDot(vNewPos,w->vOrient)))
					*dt = DBL_MAX;
				break;
			default:
				assert(0);
				}
			}
		}

	/* check again for wall face overlap (roundoff-induced...) */

	if ((*dtUpper < 0.0 && *dtLower >= 0.0) || (*dtUpper >= 0.0 && *dtLower < 0.0))
		*bOverlap = TRUE;

#ifdef DEBUG_WALLS
	if (DEBUG_TRACE(p->iOrder))
		fprintf(stderr,"%i & wall %i (type %i) final contact times: upper %g lower %g; dStart %g bOverlap %i\n",p->iOrder,iWallID,w->iType,*dtUpper,*dtLower,dStart,*bOverlap);
#endif
	}

static void wallIntersectPlane(const WALL *w,const PARTICLE *p,double dStart,double dTmax,double dtCol[],int *bOverlap)
{
	wallIntersectFlatFace(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[PlaneUpper],&dtCol[PlaneLower],bOverlap);
	}

static void wallIntersectRing(int iWallID,const WALL_DATA *w,const Vector vOscVel,const PARTICLE *p,double dStart,double dTmax,double *dtRing,int *bOverlap)
{
	/*
	** Returns time to impact ring (e.g. perimeter of finite disk, or
	** ends of finite cylinder).  Usually requires solving a quartic!
	*/

	Matrix mOrient;
	Vector vRelPos,vRelVel,vDprime,vPprime,n;
	double r2,dAlpha,R2,dBeta,dGamma,dz,pz,dCoefs[5],dRoots[4],dProjDist;
	int i,nRealRoots,bTouching;

	/* don't assert on types here, since ring is a derived type */

	*dtRing = DBL_MAX;
	*bOverlap = FALSE;

	/* form relative position and velocity vectors, accounting for motion of the ring */

	wallGetRelPosAndVel(w,vOscVel,p->r,p->v,dStart,vRelPos,vRelVel);

	r2 = RADIUS(p)*RADIUS(p); /* particle radius squared */

	if (w->dRadius == 0.0) { /* special case: single point */
		dCoefs[2] = vectorMagSq(vRelVel);
		if (dCoefs[2] == 0.0) {
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (point): no relative motion\n",p->iOrder,iWallID);
#endif
			return; /* no relative motion -- note overlap not checked for */
			}
		dCoefs[1] = 2.0*vectorDot(vRelPos,vRelVel);
		if (dCoefs[1] >= 0.0) {
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (point): moving apart\n",p->iOrder,iWallID);
#endif
			return; /* not approaching */
			}
		dCoefs[0] = vectorMagSq(vRelPos) - r2;
		if (polyQuadSolve(dCoefs[2],dCoefs[1],dCoefs[0],&dRoots[0],&dRoots[1])) {
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (point): no real roots\n",p->iOrder,iWallID);
#endif
			return; /* no real solutions ==> no collision */
			}
		if (fabs(dRoots[0]) < TINY*dTmax)
			dRoots[0] = 0.0;
		if (fabs(dRoots[1]) < TINY*dTmax)
			dRoots[1] = 0.0;
		if (dCoefs[0] > 0.0) { /* not overlapping -- choose smallest positive root */
			if (dRoots[0] > 0.0) /* roots guaranteed to have dRoots[0] <= dRoots[1] */
				*dtRing = dRoots[0];
			else if (dRoots[1] > 0.0)
				*dtRing = dRoots[1];
			else {
				if (dRoots[0] == 0.0 || dRoots[1] == 0.0)
					*bOverlap = TRUE; /* roundoff-induced overlap */
				else
					assert(0); /* all roots negative ==> moving away, but already checked for */
				}
			}
		else { /* overlap! -- choose smallest-magnitude non-positive root */
			*bOverlap = TRUE;
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (point) overlapping/touching: dCoefs[0] = %g\n",p->iOrder,iWallID,dCoefs[0]);
#endif
			if (dRoots[1] <= 0.0) /* recall roots guaranteed to have dRoots[0] <= dRoots[1] */
				*dtRing = dRoots[1];
			else if (dRoots[0] <= 0.0)
				*dtRing = dRoots[0];
			else
				assert(0); /* means dt > 0 but overlapping -- shouldn't happen! */
			}
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (point) final contact time: %g; bOverlap %i\n",p->iOrder,iWallID,*dtRing,*bOverlap);
#endif
		return;
		}

	/*
	** General case.  Strategy: treat this as a ray/torus intersection
	** following the development of Max Aubrey Wagner (2004)...
	*/

	*dtRing = DBL_MAX; /* assume no solution by default */

	/*
	** Construct matrix to transform ray into frame in which torus
	** lies in the xy-plane, with z-hat pointing up.  So mOrient below
	** transforms a real-space vector into a torus-space vector,
	** i.e. mOrient times vOrient = (0,0,1).  Note that rotation
	** around the torus axis is irrelevant, since we're only
	** interested in the time to intersect, so we can use
	** vectorGetBasis() here.
	*/

	vectorCopy(w->vOrient,mOrient[2]); /* loads the last row */
	vectorGetBasis(mOrient[2],mOrient[1],mOrient[0]);

	/*
	** Need to know point on ring closest to particle: if particle and
	** ring are actually touching, this is contact point.  Otherwise,
	** this helps us decide if particle already overlapping ring.
	*/

	vectorScale(w->vOrient,vectorDot(vRelPos,w->vOrient),n);
	vectorSub(vRelPos,n,n); /* projected vector from ring center to particle center */
	dProjDist = vectorMag(n);
	if (dProjDist == 0.0) {
		/* particle lies on ring axis -- nearest contact point arbitrary */
		vectorScale(mOrient[0],w->dRadius,n); /* mOrient[1] would work just as well */
		}
	else {
		vectorScale(n,w->dRadius/dProjDist,n); /* vector from ring center to point on ring nearest particle */
		}
	vectorSub(vRelPos,n,n); /* vector from nearest point on ring to particle center */

	if (dProjDist == 0.0 && fabs(vectorDot(vRelVel,w->vOrient)) == vectorMag(vRelVel)) {
		/* special case -- particle on and moving parallel to ring axis -- quartic solver barfs */
		if (vectorDot(vRelPos,vRelVel) >= 0.0)
			return; /* moving away */
		if (w->dRadius > RADIUS(p)) {
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (ring): particle passing through ring on axis\n",p->iOrder,iWallID);
#endif
			return; /* passing through without touching ring */
			}
		/* development very similar to point case above... */
		dCoefs[2] = vectorMagSq(vRelVel);
		assert(dCoefs[2] > 0.0);
		dCoefs[1] = 2.0*vectorDot(n,vRelVel);
		dCoefs[0] = vectorMagSq(n) - r2;
		if (polyQuadSolve(dCoefs[2],dCoefs[1],dCoefs[0],&dRoots[0],&dRoots[1]))
			assert(0); /* shouldn't happen! */
		if (fabs(dRoots[0]) < TINY*dTmax)
			dRoots[0] = 0.0;
		if (fabs(dRoots[1]) < TINY*dTmax)
			dRoots[1] = 0.0;
		if (dCoefs[0] > 0.0) { /* not overlapping -- choose smallest positive root */
			if (dRoots[0] > 0.0) /* roots guaranteed to have dRoots[0] <= dRoots[1] */
				*dtRing = dRoots[0];
			else if (dRoots[1] > 0.0)
				*dtRing = dRoots[1];
			else {
				if (dRoots[0] == 0.0 || dRoots[1] == 0.0)
					*bOverlap = TRUE; /* roundoff-induced overlap */
				else
					assert(0); /* all roots negative ==> moving away, but already checked for */
				}
			}
		else { /* overlap! -- choose smallest-magnitude non-positive root */
			*bOverlap = TRUE;
			if (dRoots[1] <= 0.0) /* recall roots guaranteed to have dRoots[0] <= dRoots[1] */
				*dtRing = dRoots[1];
			else if (dRoots[0] <= 0.0)
				*dtRing = dRoots[0];
			else
				assert(0); /* means dt > 0 but overlapping -- shouldn't happen! */
			}
		return;
		}

	matrixTransform(mOrient,vRelPos,vPprime); /* ray origin in torus frame */
	matrixTransform(mOrient,vRelVel,vDprime); /* ray direction in torus frame */

	/* useful factors, following Wagner's notation */

	dAlpha = vectorMagSq(vDprime);

	if (dAlpha == 0.0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (ring): no relative motion\n",p->iOrder,iWallID);
#endif
		return; /* no relative motion -- note overlap not checked for */
		}

	R2 = w->dRadius*w->dRadius; /* ring/torus radius squared */

	dBeta = 2.0*vectorDot(vPprime,vDprime); /* too complex to check for moving-away condition here */
	dGamma = vectorMagSq(vPprime) - r2 - R2;

	dz = vDprime[2];
	pz = vPprime[2];

	dCoefs[4] = dAlpha*dAlpha; /* t^4 coefficient */
	dCoefs[3] = 2.0*dAlpha*dBeta; /* t^3 coefficient */
	dCoefs[2] = dBeta*dBeta + 2.0*dAlpha*dGamma + 4.0*R2*dz*dz; /* t^2 coefficient */
	dCoefs[1] = 2.0*dBeta*dGamma + 8.0*R2*pz*dz; /* t coefficient */
	dCoefs[0] = dGamma*dGamma + 4.0*R2*(pz*pz - r2); /* constant */

	polyFindRealRoots(4,dCoefs,dRoots,&nRealRoots);

	if (nRealRoots == 0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (ring): no real roots\n",p->iOrder,iWallID);
#endif
		return; /* no real roots ==> no collision */
		}

	bTouching = FALSE;

	for (i=0;i<nRealRoots;i++)
		if (fabs(dRoots[i]) < TINY*dTmax) { /* allow for round-off error */
			dRoots[i] = 0.0;
			bTouching = TRUE;
			}

	if (bTouching) {
		/*
		** Need to handle special case of particle just touching ring
		** because of possible legitimate second interior bounce
		** (cf. cylinder case, which is simpler).  If particle moving
		** away from point of contact, but projected particle center
		** at or beyond ring radius, no collision occurs (presumably
		** this case arose from earlier collision during the step;
		** unfortunately, because of possible round-off error, we
		** perform this test *after* computing roots).  Otherwise, if
		** moving away and interior, there may be a second bounce, so
		** choose smallest positive root (if any).  If particle moving
		** toward point of contact, this is an overlap condition.
		*/
		if (dProjDist == 0.0) {
			/* very special case: particle touching entire ring! (note motion NOT along axis) */
			if (vectorDot(vRelPos,vRelVel) < 0.0) {
				/* attempting to penetrate ring: overlap condition */
				*bOverlap = TRUE;
				*dtRing = 0.0;
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (ring): touching but colliding\n",p->iOrder,iWallID);
#endif
				}
			else {
				/* do nothing */
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (ring): touching but moving apart\n",p->iOrder,iWallID);
#endif
				}
			}
		else {
			if (vectorDot(n,vRelVel) >= 0.0) { /* (n is vector from point of contact to particle center) */
				/* particle moving away from point of contact */
				if (dProjDist >= w->dRadius) {
					/* particle moving away from ring */
					for (i=0;i<nRealRoots;i++)
						assert(dRoots[i] <= 0.0); /* no positive roots should be found */
#ifdef DEBUG_WALLS
					if (DEBUG_TRACE(p->iOrder))
						fprintf(stderr,"%i & wall %i (ring): moving apart\n",p->iOrder,iWallID);
#endif
					}
				else {
					/* particle may hit ring a second time (interior bounce) */
					for (i=0;i<nRealRoots;i++)
						if (dRoots[i] > 0.0) /* array comes pre-sorted */
							*dtRing = dRoots[i];
					}
				}
			else {
				/* particle moving toward ring, but already touching ==> overlap condition */
				*bOverlap = TRUE;
				*dtRing = 0.0;
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (ring): touching but colliding\n",p->iOrder,iWallID);
#endif
				}
			}
		}
	else {
		/* not touching, but might be overlapping (interpenetrating)!... */
		if (vectorMag(n) > RADIUS(p)) { /* not overlapping -- choose smallest positive root */
			for (i=0;i<nRealRoots;i++)
				if (dRoots[i] > 0.0) { /* array comes pre-sorted */
					*dtRing = dRoots[i];
					break;
					}
			}
		else { /* overlap! -- choose smallest-magnitude non-positive root */
			*bOverlap = TRUE;
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (ring) overlapping: |n| %g > R %g\n",p->iOrder,iWallID,vectorMag(n),RADIUS(p));
#endif
			for (i=nRealRoots-1;i>=0;i--)
				if (dRoots[i] <= 0.0) { /* array comes pre-sorted */
					*dtRing = dRoots[i];
					break;
					}
			assert(*dtRing != DBL_MAX); /* must find a solution */
			}
		assert(*dtRing != 0.0); /* handled earlier (touching case) */
		}

#ifdef DEBUG_WALLS
	if (DEBUG_TRACE(p->iOrder))
		fprintf(stderr,"%i & wall %i (ring) final contact time: %g; bOverlap %i\n",p->iOrder,iWallID,*dtRing,*bOverlap);
#endif
	}

static void wallIntersectDisk(const WALL *w,const PARTICLE *p,double dStart,double dTmax,int bWallsEdgeDetect,double dtCol[],int *bOverlap)
{
	dtCol[DiskUpper] = dtCol[DiskLower] = dtCol[DiskRing] = dtCol[DiskHole] = DBL_MAX;

	*bOverlap = FALSE;
	
	if (w->wd.dRadius > 0.0)
		wallIntersectFlatFace(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[DiskUpper],&dtCol[DiskLower],bOverlap);

	/* check for perimeter intersection, if desired (REQUIRED for point detection, i.e. degenerate disk) */

	if (bWallsEdgeDetect) {
		int bOverlapEdge;

		wallIntersectRing(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[DiskRing],&bOverlapEdge);

		if (bOverlapEdge)
			*bOverlap = TRUE;

		if (w->wd.dHoleRadius > 0.0) {
			WALL_DATA wd;
			/* simplest approach: temporarily replace the "wall" radius with the hole radius... */
			wd = w->wd; /* struct copy */
			wd.dRadius = w->wd.dHoleRadius;
			wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[DiskHole],&bOverlapEdge);
			if (bOverlapEdge)
				*bOverlap = TRUE;
			}
		}
	}

static void wallIntersectCylinder(int iWallID,const WALL_DATA *w,const Vector vOscVel,const PARTICLE *p,double dStart,double dTmax,double *dtInner,double *dtOuter,int *bOverlap)
{
	/*
	** Returns times to impact inner and outer surface of cylinder
	** (ends of finite cylinder handled as ring intersection test by
	** calling function).  Strategy: find dt such that |r-prime -
	** (r-prime dot z-hat) z-hat| = w->dRadius +/- RADIUS(p).  NOTE:
	** taper (for finite cylinders) makes the test much more complex,
	** so that case is handled separately.
	*/

	const int nFaces = 2; /* inner & outer */
	enum {Inner=0,Outer};

	Vector vRelPos,vRelVel,vTmp,vRproj,vVproj;
	double rdotz,vdotz,r2,rdotv,v2,dProjDist,dContactInner,dContactOuter,*dt,c,dt1,dt2;
	int i,bAboveOrBelow;

	assert(w->iType == WallCylinderInfinite || w->iType == WallCylinderFinite);

	assert(dtInner != NULL && dtOuter != NULL && bOverlap != NULL);

	*dtInner = *dtOuter = DBL_MAX;
	*bOverlap = FALSE;

	/* get position and velocity relative to wall origin, accounting for any wall motion */

	wallGetRelPosAndVel(w,vOscVel,p->r,p->v,dStart,vRelPos,vRelVel);

	/* construct relative position and velocity vectors projected on plane perpendicular to cylinder axis */

	rdotz = vectorDot(vRelPos,w->vOrient);
	vectorScale(w->vOrient,rdotz,vTmp);
	vectorSub(vRelPos,vTmp,vRproj);

	vdotz = vectorDot(vRelVel,w->vOrient);
	vectorScale(w->vOrient,vdotz,vTmp);
	vectorSub(vRelVel,vTmp,vVproj);
		
	r2 = vectorMagSq(vRproj);
	v2 = vectorMagSq(vVproj);

	/* check for return conditions */

	if (v2 == 0.0) {
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (cylinder): no motion relative to (infinite) axis\n",p->iOrder,iWallID);
#endif
		return; /* moving parallel to cylinder axis (or stationary) -- no collision */
		}

	rdotv = vectorDot(vRproj,vVproj);

	dProjDist = sqrt(r2);
	dContactInner = w->dRadius - RADIUS(p);
	dContactOuter = w->dRadius + RADIUS(p);

	if (rdotv > 0.0 && dProjDist >= w->dRadius) { /* dRadius instead of dContactOuter to allow for round-off error */
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (cylinder): moving apart\n",p->iOrder,iWallID);
#endif
		return; /* moving away from cylinder -- no collision */
		}

	/* check for overlap */

	bAboveOrBelow = (w->iType == WallCylinderFinite && fabs(rdotz) >= 0.5*w->dLength); /* handy for finite case (>=, not >, because we don't check for end intersection here) */

	if (dProjDist >= dContactInner && dProjDist <= dContactOuter && !bAboveOrBelow) {
		*bOverlap = TRUE;
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i overlapping/touching wall %i (cylinder): r-CI %g r-CO %g\n",p->iOrder,iWallID,dProjDist - dContactInner,dProjDist - dContactOuter);
#endif
		}

	/*
	** Need to check both inside and outside cylinder, because even
	** though a particle inside an infinite cylinder can't collide
	** with the outside face, and vice versa (whereas for a finite
	** cylinder, an exterior particle can still hit an inside wall),
	** there is still the possibility of an overlap and both faces
	** will need to be considered in order to get the minimal required
	** correction, if allowed.
	*/

	for (i=0;i<nFaces;i++) {
		switch (i) {
		case Inner:
			if (dContactInner <= 0.0) {
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (cylinder): no inside face (%i) contact possible\n",p->iOrder,iWallID,i);
#endif
				continue; /* special case of particle diameter larger than or equal to cylinder diameter */
				}
			c = r2 - dContactInner*dContactInner; /* constant term in quadratic equation */
			dt = dtInner;
			break;
		case Outer:
			if (dProjDist < dContactInner) {
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (cylinder): no outside face (%i) contact possible\n",p->iOrder,iWallID,i);
#endif
				continue; /* can't hit outside face from inside, even with finite cylinder */
				}
			c = r2 - dContactOuter*dContactOuter;
			dt = dtOuter;
			break;
		default:
			assert(0);
			}
		if (polyQuadSolve(v2,2.0*rdotv,c,&dt1,&dt2)) {
			assert(dProjDist >= dContactInner);
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i (cylinder) face %i: no real roots\n",p->iOrder,iWallID,i);
#endif
			continue; /* no real roots ==> no collision (if outside or overlapping/touching only) */
			}
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"%i & wall %i (cylinder) contact times (face %i): dt1 %g dt2 %g; dStart %g\n",p->iOrder,iWallID,i,dt1,dt2,dStart);
#endif
		if (fabs(dt1) < TINY*dTmax)
			dt1 = 0.0;
		if (fabs(dt2) < TINY*dTmax)
			dt2 = 0.0;
		/* treat cases individually */
		if (dProjDist < dContactInner && !bAboveOrBelow) {
			/* completely inside */
			assert(dt1 <= 0.0 && dt2 >= 0.0); /* recall roots are sorted numerically; "=" because of round-off correction above */
			*dt = dt2; /* true for either face */
			}
		else if ((dt1 == 0.0 || dt2 == 0.0) && i == Inner && !bAboveOrBelow) {
			/* just touching inside -- there could be legit solution (2nd bounce on inner face) */
			if (dt2 > 0.0) {
				*dt = dt2;
				*bOverlap = FALSE; /* turn off flag to allow this collision to be processed */
				}
			else
				*dt = 0.0;
			}
		else if (dProjDist > dContactOuter && !bAboveOrBelow) {
			/* completely outside */
			assert(dt1 >= 0.0 && dt2 >= 0.0);
			*dt = dt1; /* true for either face */
			}
		else if ((dt1 == 0.0 || dt2 == 0.0) && i == Outer && !bAboveOrBelow) {
			/* just touching outside -- must be a missed collision (because already checked that particle is moving inward) */
#ifdef DEBUG_WALLS
			if (!((dt1 == 0.0 && dt2 > 0.0) || (dt1 > 0.0 && dt2 == 0.0)))
				fprintf(stderr,"%i & wall %i ASSERT: dt1 %g dt2 %g\n",p->iOrder,iWallID,dt1,dt2);
#endif
			assert((dt1 == 0.0 && dt2 > 0.0) || (dt1 > 0.0 && dt2 == 0.0));
			*dt = 0.0;
			}
		else if (bAboveOrBelow) {
			/* outside cylinder, but could hit either face -- take smallest positive time to contact, if any */
			if (dt1 > 0.0)
				*dt = dt1;
			else if (dt2 > 0.0)
				*dt = dt2;
			else {
				*dt = DBL_MAX; /* not approaching cylinder */
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (cylinder) face %i: above or below cylinder and not approaching\n",p->iOrder,iWallID,i);
#endif
				}
			}
		else if (*bOverlap) { /* flag may have been turned off for special case (inner face) above */
			/* overlapping (or touching face) */
			if (dt2 <= 0.0)
				*dt = dt2;
			else
				*dt = dt1;
#ifdef DEBUG_WALLS
			if (*dt > 0.0)
				fprintf(stderr,"%i & wall %i ASSERT: iFace %i dt1 %g dt2 %g r-CI %.16e r-CO %.16e\n",p->iOrder,iWallID,i,dt1,dt2,dProjDist - dContactInner,dProjDist - dContactOuter);
#endif
			assert(*dt <= 0.0);
			}
		if (w->iType == WallCylinderFinite) { /* check for finite cylinder intersection */
			vectorScale(vRelVel,*dt,vTmp); /* change in position to point of contact */
			vectorAdd(vTmp,vRelPos,vTmp); /* position relative to cylinder origin at point of contact */
			if (fabs(vectorDot(vTmp,w->vOrient)) >= 0.5*w->dLength) {
				*dt = DBL_MAX; /* missed cylinder face -- no collision */
#ifdef DEBUG_WALLS
				if (DEBUG_TRACE(p->iOrder))
					fprintf(stderr,"%i & wall %i (cylinder) face %i: missed (finite) cylinder face\n",p->iOrder,iWallID,i);
#endif
				}
			}
		if (dt1 == 0.0 && dt2 > 0.0 && i == Inner && !bAboveOrBelow)
			++i; /* skip outer face for this special case */
		}

	/* check again for cylinder overlap (roundoff-induced...) */

	if ((*dtInner < 0.0 && *dtOuter >= 0.0) || (*dtInner >= 0.0 && *dtOuter < 0.0) || (dContactInner <= 0.0 && *dtOuter == 0.0 /* (touching line) */))
		*bOverlap = TRUE;

	if (*bOverlap) {
#ifdef DEBUG_WALLS
		if ((*dtInner == DBL_MAX && dContactInner > 0.0) || *dtOuter == DBL_MAX)
			fprintf(stderr,"%i & wall %i ASSERT: no overlap backstep found\n",p->iOrder,iWallID);
#endif
/*DEBUG!! don't understand why this trips so often		assert((*dtInner < DBL_MAX || dContactInner <= 0.0) && *dtOuter < DBL_MAX);*/ /* otherwise, unphysical outcome likely (large corrective step) */
		}

#ifdef DEBUG_WALLS
	if (DEBUG_TRACE(p->iOrder))
		fprintf(stderr,"%i & wall %i (cylinder) final contact times: inner %g outer %g; dStart %g bOverlap %i\n",p->iOrder,iWallID,*dtInner,*dtOuter,dStart,*bOverlap);
#endif
	}

static void wallIntersectRectangle(const WALL *w,const PARTICLE *p,double dStart,double dTmax,int bWallsEdgeDetect,double dtCol[],int *bOverlap)
{
	int i;

	for (i=0;i<nWallFaces[WallRectangle];i++)
		dtCol[i] = DBL_MAX;

	*bOverlap = FALSE;

	wallIntersectFlatFace(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[RectangleUpper],&dtCol[RectangleLower],bOverlap);

	/* check for perimeter intersection */

	if (bWallsEdgeDetect) {
		/*
		** Strategy: treat as 4 point intersection tests (for corners)
		** and 4 line intersection tests (for sides).
		*/
		WALL_DATA wd;
		double dWidth,dLength,dDum;
		int bOverlapEdge;
		wd = w->wd; /* struct copy */
		wd.iType = WallCylinderFinite; /* needed for sides */
		wd.dRadius =  0.0; /* degenerate ring (point)/cylinder (line) */
		dWidth = vectorMag(wd.vVertex1); /* sides 1 & 3 */
		dLength = vectorMag(wd.vVertex2); /* sides 2 & 4 */
		assert(dWidth > 0.0 && dLength > 0.0);
		/* first corner: origin */
		wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[RectangleCorner1],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* first side: origin to vertex 1 */
		vectorScale(w->wd.vVertex1,0.5,wd.vOrient); /* half vector for origin calc */
		vectorAdd(w->wd.vOrigin,wd.vOrient,wd.vOrigin); /* (origin at center of cylinder) */
		vectorScale(w->wd.vVertex1,1.0/dWidth,wd.vOrient); /* unit vector */
		wd.dLength = dWidth;
		wallIntersectCylinder(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dDum,&dtCol[RectangleSide1],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* second corner: vertex 1 */
		vectorAdd(w->wd.vOrigin,w->wd.vVertex1,wd.vOrigin);
		wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[RectangleCorner2],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* second side: origin to vertex 2 */
		vectorScale(w->wd.vVertex2,0.5,wd.vOrient);
		vectorAdd(w->wd.vOrigin,wd.vOrient,wd.vOrigin);
		vectorScale(w->wd.vVertex2,1.0/dLength,wd.vOrient);
		wd.dLength = dLength;
		wallIntersectCylinder(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dDum,&dtCol[RectangleSide2],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* third corner: vertex 2 */
		vectorAdd(w->wd.vOrigin,w->wd.vVertex2,wd.vOrigin);
		wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[RectangleCorner3],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* fourth corner: opposite origin */
		vectorAdd(wd.vOrigin,w->wd.vVertex1,wd.vOrigin);
		wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[RectangleCorner4],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* third side: vertex 1 to corner opposite origin */
		vectorScale(w->wd.vVertex2,0.5,wd.vOrient);
		vectorSub(wd.vOrigin,wd.vOrient,wd.vOrigin);
		vectorScale(w->wd.vVertex2,1.0/dLength,wd.vOrient);
		wallIntersectCylinder(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dDum,&dtCol[RectangleSide3],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		/* fourth side: vertex 2 to corner opposite origin */
		vectorAdd(w->wd.vOrigin,w->wd.vVertex2,wd.vOrigin);
		vectorScale(w->wd.vVertex1,0.5,wd.vOrient);
		vectorAdd(wd.vOrigin,wd.vOrient,wd.vOrigin);
		vectorScale(w->wd.vVertex1,1.0/dWidth,wd.vOrient);
		wd.dLength = dWidth;
		wallIntersectCylinder(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dDum,&dtCol[RectangleSide4],&bOverlapEdge);
		if (bOverlapEdge)
			*bOverlap = TRUE;
		}
	}

static void wallIntersectCylinderInfinite(const WALL *w,const PARTICLE *p,double dStart,double dTmax,double dtCol[],int *bOverlap)
{
	wallIntersectCylinder(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[CylinderInfiniteInner],&dtCol[CylinderInfiniteOuter],bOverlap);
	}

#ifdef BROKEN /*DEBUG! omitting this for now (funnel)*/
static void wallIntersectFunnel(int WallID,const Vector vOrigin,const Vector vZhat,double R/*radius*/,double L/*length*/,double T/*taper*/,const PARTICLE *p,double dStart/*not used*/,double dTmax,double *dtInner,double *dtOuter,int *bOverlap)
{
	/* need comment */

	const int nFaces = 2; /* inner & outer */
	enum {Inner=0,Outer};

	Vector vRelPos,vRelVel;
	double *dt,r2,v2,rdotv,rdotz,rz2,vdotz,vz2,rzvz,dTanTheta,dTheta,dSinTheta,dCosTheta,dGamma,dDelta,d,e,a,b,c,dt1,dt2,dOffset;
	int i;

	*dtInner = *dtOuter = DBL_MAX;
	*bOverlap = FALSE;

	/* preliminaries */

	vectorSub(p->r,vOrigin,vRelPos);
	r2 = vectorMagSq(vRelPos);
	vectorCopy(p->v,vRelVel); /* currently funnels are stationary */
	v2 = vectorMagSq(vRelVel);
	if (v2 == 0.0)
		return; /* note: overlap not checked for */
	rdotv = vectorDot(vRelPos,vRelVel);
	rdotz = vectorDot(vRelPos,vZhat);
	rz2 = rdotz*rdotz;
	vdotz = vectorDot(vRelVel,vZhat);
	vz2 = vdotz*vdotz;
	rzvz = rdotz*vdotz;
	assert(L > 0.0 && L < DBL_MAX);
	dTanTheta = R*T/L;
	dTheta = atan(dTanTheta); /* slope angle of funnel, in radians */
	dSinTheta = sin(dTheta);
	dCosTheta = cos(dTheta);

	dGamma = -dTanTheta/dCosTheta;

	d = 1.0 - dGamma*(dSinTheta - dGamma);
	e = dSinTheta - 2.0*dGamma;

	a = v2 - d*vz2;
	b = 2.0*(rdotv - d*rzvz); /* plus a face-dependent term */
	c = r2 - d*rz2; /* plus two face-dependent terms */

	/* solve quadratic for both faces */

	for (i=0;i<nFaces;i++) {
		switch (i) {
		case Inner:
			dt = dtInner;
			dDelta = R*(1 - 0.5*T - RADIUS(p)*T*dSinTheta/L)/dCosTheta;
			dOffset = dDelta - RADIUS(p);
			break;
		case Outer:
			dt = dtOuter;
			dDelta = R*(1 - 0.5*T + RADIUS(p)*T*dSinTheta/L)/dCosTheta;
			dOffset = dDelta + RADIUS(p);
			break;
		default:
			assert(0);
			}

		if (polyQuadSolve(a,b + e*dOffset*vdotz,c + e*dOffset*rdotz - dOffset*dOffset,&dt1,&dt2))
			dt1 = dt2 = DBL_MAX; /* no solutions ==> no collision */

		if (fabs(dt1) < TINY*dTmax)
			dt1 = 0.0;

		if (fabs(dt2) < TINY*dTmax)
			dt2 = 0.0;

		*dt = (dt1 < dt2 && dt1 >= 0.0 ? dt1 : dt2 >= 0.0 ? dt2 : DBL_MAX);

		{
			Vector vNewPos;
			double f,dNewrdotz;
			vectorScale(p->v,*dt,vNewPos);
			vectorAdd(vNewPos,p->r,vNewPos);
			vectorSub(vNewPos,vOrigin,vNewPos); /* assumes stationary cylinder */
			dNewrdotz = vectorDot(vNewPos,vZhat);
			if (i == Inner)
				f = (dNewrdotz + 0.5*L + RADIUS(p)*dSinTheta)/L;
			else
				f = (dNewrdotz + 0.5*L - RADIUS(p)*dSinTheta)/L;
			if (f <= 0.0 || f >= 1.0) *dt = DBL_MAX;
			}
		}

	/* check for funnel overlap (roundoff error happens...) */

	if ((*dtInner < 0.0 && *dtOuter >= 0.0) || (*dtInner >= 0.0 && *dtOuter < 0.0))
		*bOverlap = TRUE;
	}
#endif /*DEBUG! BROKEN (funnel)*/

static void wallIntersectCylinderFinite(const WALL *w,const PARTICLE *p,double dStart,double dTmax,int bWallsEdgeDetect,double dtCol[],int *bOverlap)
{
	dtCol[CylinderFiniteInner] = dtCol[CylinderFiniteOuter] = dtCol[CylinderFiniteRingUpper] = dtCol[CylinderFiniteRingLower] = DBL_MAX;

	*bOverlap = FALSE;

	if (w->wd.dLength > 0.0) {
		if (w->wd.dRadius > 0.0 && w->wd.dTaper > 0.0)
			assert(0);
#ifdef BROKEN /*DEBUG! (funnel)*/
			wallIntersectFunnel(w->iWallID,wcf->vOrigin,wcf->vZhat,wcf->dRadius,wcf->dLength,wcf->dTaper,p,dStart,dTmax,&dtCol[CylinderFiniteInner],&dtCol[CylinderFiniteOuter],bOverlap);
#endif /*DEBUG! (funnel)*/
		else
			wallIntersectCylinder(w->iWallID,&w->wd,w->vOscVel,p,dStart,dTmax,&dtCol[CylinderFiniteInner],&dtCol[CylinderFiniteOuter],bOverlap);
		}

	/* now check intersections with ends of cylinders, if desired (REQUIRED for ring/point detection, i.e. degenerate cylinder) */

	if (bWallsEdgeDetect) {
		WALL_DATA wd;
		Vector vOffset;
		int bOverlapRing;

		wd = w->wd; /* struct copy */
		vectorScale(wd.vOrient,-0.5*wd.dLength,vOffset); /* lower (base) ring first */
		vectorAdd(w->wd.vOrigin,vOffset,wd.vOrigin);
		wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[CylinderFiniteRingLower],&bOverlapRing);
		if (bOverlapRing)
			*bOverlap = TRUE;

		if (wd.dLength > 0.0) {
#ifdef DEBUG_WALLS
			if (DEBUG_TRACE(p->iOrder))
				fprintf(stderr,"%i & wall %i: checking second (upper) ring...\n",p->iOrder,w->iWallID);
#endif
			vectorSub(w->wd.vOrigin,vOffset,wd.vOrigin);
			wd.dRadius *= (1.0 - wd.dTaper);
			wallIntersectRing(w->iWallID,&wd,w->vOscVel,p,dStart,dTmax,&dtCol[CylinderFiniteRingUpper],&bOverlapRing);
			if (bOverlapRing)
				*bOverlap = TRUE;
			}
		}
	}

void wallsGetTimesToIntersect(const WALL_PARAMS *WP,const WALL *w,PARTICLE *p,double dStart,double dTmax,double dtCol[],int *bOverlap)
{
	/*
	** Determines possible impact times of particle p with wall w and
	** sets a flag if the particle is overlapping the wall.  The
	** maximum timestep (dTmax) is used to deal with roundoff error
	** (see TINY macro).
	*/

	Vector vOrgPos;
	int i,iType = w->wd.iType;

	/*
	** Because a particle, when backdrifted, may be intersecting a
	** wall, such that it doesn't check for a subsequent collision
	** with the same wall (which is possible with, e.g. cylindrical
	** walls), we *forward drift* the particles now to the start time
	** of the current collision search (supplied by
	** smoothfcn.c:CheckForCollisions()), and later restore the
	** start-of-step positions, accounting for the extra time interval
	** in the returned collision time predictions.  We could restrict
	** this to just walls that could be affected by this issue
	** (e.g. cylinders, not planes), but it's quick to do and provides
	** an extra error check (because after forward drifting, the
	** particle should NOT be intersecting anything), so we do it for
	** all cases here.
	**
	** NOTE: because we advance forward, at machine precision the
	** particle will usually be in perfect contact with a wall after
	** its last collision with that wall -- this is another reason why
	** we use bOverlap (in this case we get a dt of zero).
	**
	** Also: walls can move -- this needs to be taken into account to
	** get the prediction right if dStart > 0, so we pass dStart to
	** all functions.
	*/

	vectorCopy(p->r,vOrgPos);
	for (i=0;i<3;i++)
		p->r[i] += p->v[i]*dStart;

	switch (iType) {
	case WallPlane:
		wallIntersectPlane(w,p,dStart,dTmax,dtCol,bOverlap);
		break;
	case WallTriangle:
	    /*wallTriangleIntersect(w,p,dStart,dTmax,WP->bWallsEdgeDetect,dtCol,bOverlap);*/
	    fprintf(stderr, "Triangle wall not available for hard sphere.\n");
	    exit(1);
	    break;
	case WallRectangle:
		wallIntersectRectangle(w,p,dStart,dTmax,WP->bWallsEdgeDetect,dtCol,bOverlap);
		break;
	case WallDisk:
		wallIntersectDisk(w,p,dStart,dTmax,WP->bWallsEdgeDetect,dtCol,bOverlap);
		break;
	case WallCylinderInfinite:
		wallIntersectCylinderInfinite(w,p,dStart,dTmax,dtCol,bOverlap);
		break;
	case WallCylinderFinite:
		wallIntersectCylinderFinite(w,p,dStart,dTmax,WP->bWallsEdgeDetect,dtCol,bOverlap);
		break;
	case WallShell:
	    /*wallIntersectShell(w,p,dStart,dTmax,WP->bWallsEdgeDetect,dtCol,bOverlap);*/
	    fprintf(stderr, "Shell wall not available for hard sphere.\n");
	    exit(1);
	    break;
	default:
		assert(0);
		}

	if (!*bOverlap) { /* sanity check */
		int bForwardStepFound = FALSE;
		for (i=0;i<nWallFaces[iType];i++)
			if (dtCol[i] > 0.0) {
				bForwardStepFound = TRUE;
				break;
				}
#ifdef SPRINGS
		(void) bForwardStepFound; /* suppress set-but-not-used warning */
#else
		/*DEBUG! at the moment, this assert needs to be commented out to work w/springs */
		assert(bForwardStepFound);
#endif
		}

	for (i=0;i<nWallFaces[iType];i++)
		if (dtCol[i] < DBL_MAX)
			dtCol[i] += dStart;

#ifdef DEBUG_WALLS
	if (DEBUG_TRACE(p->iOrder)) {
		fprintf(stderr,"%i & wall %i reported contact times: ",p->iOrder,w->iWallID);
		for (i=0;i<nWallFaces[iType];i++)
			fprintf(stderr,"%g (%g) ",dtCol[i],dtCol[i]-dStart);
		fprintf(stderr,"\n");
		}
#endif

	vectorCopy(vOrgPos,p->r);
	}

/*
** Following functions resolve wall-particle collisions for various
** wall types.  General strategy: use equations (15) and (17) from
** Richardson (1994) with particle 1 representing the wall (with
** infinite mass, radius, etc.).  The trick is to determine n-hat for
** each case, which is a function of the wall geometry (and the impact
** point).  Note these routines assume the particle has been advanced
** to the point of contact.  Output are the new velocity and spin of
** the particle as a result of the collision.  Particle iOrder is
** passed for debugging purposes.
*/

static void wallBounce(int iWallID,const WALL_DATA *w,int iOrder,const COLLISION_PARAMS *CP,int bTinyStep,const Vector n,const Vector v,const Vector s1,double dRadius,Vector vVel,Vector vSpin)
{
	/*
	** Once the normal (n) has been determined, the equations are
	** solved the same way for all geometries via this function.  Also
	** input are the relative velocity, the velocity contribution from
	** any wall spin, radius of the particle, and its velocity and
	** spin (both overwritten).  The bTinyStep flag indicates whether
	** inelastic collapse prevention is in effect.
	*/

	Vector R,s2,s,u,un,ut,Rxu;
	double vdotn,dEpsN,dEpsT,udotn;

	vdotn = vectorDot(v,n);
#ifdef DEBUG_WALLS
	if (vdotn >= 0.0)
		fprintf(stderr,"%i & wall %i ASSERT: vdotn %g v %g %g %g n %g %g %g\n",iOrder,iWallID,vdotn,v[0],v[1],v[2],n[0],n[1],n[2]);
#endif
	assert(vdotn < 0.0); /* otherwise not approaching */
	if (bTinyStep) {
		dEpsN = CP->dCollapseEpsN;
		dEpsT = CP->dCollapseEpsT;
		}
	else if (CP->iSlideOption == MaxTrv && fabs(vdotn) < CP->dSlideLimit) {
		dEpsN = CP->dSlideEpsN;
		dEpsT = CP->dSlideEpsT;
		}
	else {
		dEpsN = w->dEpsN;
		dEpsT = w->dEpsT;
		}
	/*
	** For special case of particle hitting wall, restitution
	** equations for particle reduce to:
	**
	**    v' = v - (1 + e_n) u_n - (2/7) (1 - e_t) u_t
	**    w' = w - (5/7) 1/R^2 (1 - e_t) (R X u)
	**
	** where R is vector in direction (- n-hat) with magnitude
	** equal to particle radius.
	*/
	vectorScale(n,-dRadius,R);
	vectorCross(vSpin,R,s2);
	vectorSub(s2,s1,s);
	vectorAdd(v,s,u); /* total relative velocity at impact point */
	udotn = vectorDot(u,n);
	vectorScale(n,udotn,un); /* normal component */
	vectorSub(u,un,ut); /* tangential component */
	vectorScale(un,1.0 + dEpsN,un);
	vectorSub(vVel,un,vVel);
	vectorScale(ut,(2.0/7.0)*(1.0 - dEpsT),ut);
	vectorSub(vVel,ut,vVel); /* new velocity */
	vectorCross(R,u,Rxu);
	vectorScale(Rxu,(5.0/7.0)*(1.0 - dEpsT)/(dRadius*dRadius),Rxu);
	vectorSub(vSpin,Rxu,vSpin); /* new spin */
	}

static void wallGetNormalToRectanglePerimeter(const Vector vRelPos,const WALL_DATA *w,Vector vNormal)
{
	/* only overwrites vNormal if actually in contact with perimeter */

	Matrix m;
	Vector v;
	double dWidth,dLength;

	wallGetRectangleGeometry(w,&dWidth,&dLength,m);

	matrixTransform(m,vRelPos,v); /* rel pos in rectangle body frame */

	if (v[0] <= 0.0 || v[0] >= dWidth || v[1] <= 0.0 || v[1] >= dLength) {
		/* in contact with rectangle perimeter */
		Matrix mTrans;
		Vector vTmp;
		/* check corners first */
		if (v[0] <= 0.0 && v[1] <= 0.0)
			vectorCopy(v,vNormal);
		else if (v[0] >= dWidth && v[1] <= 0.0)
			vectorSet(vNormal,v[0] - dWidth,v[1],v[2]);
		else if (v[0] <= 0.0 && v[1] >= dLength)
			vectorSet(vNormal,v[0],v[1] - dLength,v[2]);
		else if (v[0] >= dWidth && v[1] >= dLength)
			vectorSet(vNormal,v[0] - dWidth,v[1] - dLength,v[2]);
		/* now sides */
		else if (v[0] <= 0.0)
			vectorSet(vNormal,v[0],0.0,v[2]);
		else if (v[0] >= dWidth)
			vectorSet(vNormal,v[0] - dWidth,0.0,v[2]);
		else if (v[1] <= 0.0)
			vectorSet(vNormal,0.0,v[1],v[2]);
		else if (v[1] >= dLength)
			vectorSet(vNormal,0.0,v[1] - dLength,v[2]);
		else
			assert(0);

		/* transform back to space frame and normalize */

		matrixTranspose(m,mTrans);
		matrixTransform(mTrans,vNormal,vTmp);
		vectorNorm(vTmp);
		vectorCopy(vTmp,vNormal);
		}
	}

static void wallGetNormalToDiskPerimeter(const Vector vRelPos,const Vector vProj,const WALL_DATA *w,Vector vNormal)
{
	/* only overwrites vNormal if actually in contact with perimeter */

	double dRadDist;

	dRadDist = vectorMag(vProj); /* projected particle center distance from disk origin */
	if (dRadDist >= w->dRadius) {
		/* in contact with disk perimeter */
		if (dRadDist == 0.0)
			vectorCopy(vRelPos,vNormal); /* special case: "disk" is a singularity! */
		else {
			Vector vTmp;
			vectorScale(vProj,w->dRadius/dRadDist,vTmp);
			vectorSub(vRelPos,vTmp,vNormal); /* vector from point of contact to particle center */
			}
		vectorNorm(vNormal);
		}
	else if (dRadDist > 0.0 && w->dHoleRadius > 0.0 && dRadDist <= w->dHoleRadius) {
		/* in contact with hole edge */
		Vector vTmp;
		vectorScale(vProj,w->dHoleRadius/dRadDist,vTmp);
		vectorSub(vRelPos,vTmp,vNormal);
		vectorNorm(vNormal);
		}
	}

static void wallBounceFlatFace(int iWallID,const WALL_DATA *w,const Vector vOscVel,const COLLISION_PARAMS *CP,COLLIDER *c1,double dt)
{
	/*
	** For this case, normal at impact point is just +/- vOrient
	** (depending on whether the particle hits the upper or lower face
	** of the wall), except that for finite geometries, if the
	** particle center lies at or beyond the perimeter, contact is
	** with the perimeter, in which case the normal points from the
	** perimeter contact point to the particle center.
	*/

	Vector vRelPos,vRelVel,n,vProj,vWallSpin;
	double rdotn;

	assert(w->iType == WallPlane || w->iType == WallRectangle ||
		   w->iType == WallTriangle || w->iType == WallDisk);

	wallGetRelPosAndVel(w,vOscVel,c1->r,c1->v,dt,vRelPos,vRelVel);

	/*
	** The assert below applies only to the infinite plane case: if
	** rdotn is zero, the particle center is exactly in the plane, so
	** the normal to use is ambiguous.  For finite geometries, rdotn
	** of zero implies the particle is in contact with the perimeter,
	** so the correct normal (which will be in the plane) will be
	** calculated later on.
	*/

	rdotn = vectorDot(vRelPos,w->vOrient);
	assert(rdotn != 0.0 || w->iType != WallPlane);
	vectorScale(w->vOrient,dsgn(rdotn),n); /* may be zero vector */

	vectorScale(w->vOrient,rdotn,vProj);
	vectorSub(vRelPos,vProj,vProj); /* relative position vector projected on wall, needed in some cases below */

	if (w->dAngSpeed != 0.0) { /* get any spin contribution to wall velocity */
		Vector vTmp;
		vectorScale(w->vOrient,w->dAngSpeed,vTmp);
		vectorCross(vTmp,vProj,vWallSpin);
		}
	else
		vectorZero(vWallSpin);

	/* for finite geometries, check for perimeter bounce */

	switch (w->iType) {
	case WallPlane:
		break;
	case WallTriangle:
		assert(0);
	case WallRectangle:
		wallGetNormalToRectanglePerimeter(vRelPos,w,n);
		break;
	case WallDisk: /* zero radius ok */
		wallGetNormalToDiskPerimeter(vRelPos,vProj,w,n);
		break;
	default:
		assert(0);
		}

	assert(vectorMag(n) > 0.0); /* must have found a legitimate normal by now */
	wallBounce(iWallID,w,c1->id.iOrder,CP,c1->bTinyStep,n,vRelVel,vWallSpin,c1->fRadius,c1->v,c1->w);
	}

static void wallBounceCylinder(int iWallID,const WALL_DATA *w,const Vector vOscVel,const COLLISION_PARAMS *CP,COLLIDER *c1,double dt)
{
	/*
	** Here the normal points away from the cylinder wall (either
	** inside or outside), but the wall may be at an angle (if the
	** cylinder is tapered) and there is the possibility the contact
	** point is with either end of the cylinder.
	*/

	Vector vRelPos,vRelVel,vPerp,n,wxR;
	double rdotz,dRadDist,dRadLim,dRotRad;

	assert(w->iType == WallCylinderInfinite || w->iType == WallCylinderFinite);

	wallGetRelPosAndVel(w,vOscVel,c1->r,c1->v,dt,vRelPos,vRelVel);

	rdotz = vectorDot(vRelPos,w->vOrient);
	vectorScale(w->vOrient,rdotz,vPerp);
	vectorSub(vRelPos,vPerp,vPerp); /* perpendicular vector from cylinder axis to particle center */
	dRadDist = vectorMag(vPerp);

	if (w->iType == WallCylinderInfinite || w->dTaper == 0.0 || w->dLength == 0.0)
		dRadLim = w->dRadius;
	else
		dRadLim = w->dRadius*(1.0 - (rdotz + 0.5*w->dLength)*w->dTaper/w->dLength);

	dRotRad = dRadLim;

	/* note: not checking for most cylinder overlap conditions here... */

	vectorZero(n); /* to keep compiler happy */

	/* compute n-hat, which points from point of contact to particle center */

	if (w->iType == WallCylinderInfinite || fabs(rdotz) < 0.5*w->dLength) {
		assert(w->dLength < DBL_MAX || dRadDist > 0.0); /* cannot be on cylinder axis, inside (infinite) cylinder, and hitting cylinder all at the same time */
		if (dRadDist < dRadLim)
			vectorScale(vPerp,-1.0/dRadDist,n); /* inside cylinder radius */
		else
			vectorScale(vPerp,1.0/dRadDist,n); /* outside cylinder radius */
		}

	/* extra complications for finite cylinders */

	if (w->iType == WallCylinderFinite) {
		double dLenLim;
		if (w->dTaper == 0.0 || w->dLength == 0.0)
			dLenLim = 0.5*w->dLength;
		else
			dLenLim = 0.5*w->dLength*sqrt(1.0 + w->dRadius*w->dRadius*w->dTaper*w->dTaper/(w->dLength*w->dLength));
		if (fabs(rdotz) >= dLenLim) { /* in contact with cylinder end */
			if (dRadDist == 0.0)
				vectorCopy(vRelPos,n); /* bouncing on cylinder axis */
			else {
				Vector vTmp1,vTmp2;
				/* need to reconstuct vectors at exact end of cylinder */
				vectorScale(w->vOrient,0.5*w->dLength*dsgn(rdotz),vTmp1); /* vector from cylinder origin to end of cylinder, along axis */
				vectorScale(vPerp,dRadLim/dRadDist,vTmp2); /* perpendicular vector from ring axis to point of contact */
				vectorAdd(vTmp1,vTmp2,vTmp1); /* vector from cylinder origin to point of contact on ring */
				vectorSub(vRelPos,vTmp1,n); /* vector from point of contact to particle center */
				}
			vectorNorm(n); /* or divide by c1->fRadius */
			}
		else if (w->dTaper > 0.0 && w->dLength > 0.0 && w->dRadius > 0.0) { /* cone (frustum really) -- must rotate normal */
			Vector zXr;
			double dTheta;
			assert(w->dLength < DBL_MAX); /* must have finite length */
			vectorCross(w->vOrient,vRelPos,zXr);
			vectorNorm(zXr); /* we will be rotating normal about this axis */
			dTheta = -atan(w->dRadius*w->dTaper/w->dLength); /* slope angle of funnel; this is minus the rotation angle */
			vectorRotate(n,zXr,dTheta);
			dRotRad = vectorMag(vPerp) - c1->fRadius*cos(dTheta);
			}
		}

	if (dRadDist == 0.0)
		vectorZero(wxR);
	else {
		Vector vOmegaR;
		vectorScale(w->vOrient,w->dAngSpeed*dRotRad/dRadDist,vOmegaR);
		vectorCross(vOmegaR,vPerp,wxR); /* cylinder spin: Omega X (cylinder radius in the Vperp direction) */
		}

	wallBounce(iWallID,w,c1->id.iOrder,CP,c1->bTinyStep,n,vRelVel,wxR,c1->fRadius,c1->v,c1->w);
	}

void pkdWallsDoCollision(PKD pkd,const COLLISION_PARAMS *CP,COLLIDER *c1,const COLLIDER *c2,double dt,int *piOutcome)
{
	/*
	** Handles collision between particle (stored in c1) and wall (in
	** c2).  The particle (and wall, indirectly) are advanced to the
	** impact time dt.  The outcome type is returned in piOutcome
	** (currently death, merge, or bounce).  NOTE: just to be
	** confusing, the restitution equations in wallBounce() take the
	** wall to have index 1 and the particle to have index 2, which is
	** the opposite of the c1 and c2 convention used here...
	*/

	const WALL_PARAMS *WP = &CP->WP;
	const WALL *w;
	const WALL_DATA *wd;
	int i,iWallID;

	iWallID = -1 - c2->id.iOrder;
	assert(iWallID >= 0 && iWallID < WP->nWalls);
	w = &WP->pWalls[iWallID];
	assert(iWallID == w->iWallID);
	wd = &w->wd;
	for (i=0;i<3;i++)
		c1->r[i] += c1->v[i]*dt; /* move particle to point of contact */
#ifdef DEBUG_WALLS
	if (DEBUG_TRACE(c1->id.iOrder))
		fprintf(stderr,"%i & wall %i bounce before: r %g %g %g v %g %g %g dt %g\n",c1->id.iOrder,iWallID,c1->r[0],c1->r[1],c1->r[2],c1->v[0],c1->v[1],c1->v[2],dt);
#endif
	if (wd->dEpsN < 0.0) /* "death" wall! (particle removed) */
		*piOutcome = DEATH;
	else if (wd->dEpsN == 0.0) { /* sticky wall! */
		for (i=0;i<3;i++)
			c1->v[i] = c1->w[i] = 0.0; /* v overwritten in pkdWallsUpdateStuckParticles(); w may be overwritten below */
		c1->iColor = -1 - iWallID; /* or c2->id.iOrder */
		pkd->pStore[c1->id.iIndex].iColor = c1->iColor; /* needed because pkdPutColliderInfo() doesn't store color *//*DEBUG maybe it should??*/
		if (COLLIDER_STUCK_ON_ROTATING_WALL(c1,WP)) {
			/* transfer wall spin to particle spin */
			assert(wd->iType != WallTriangle && wd->iType != WallRectangle);
			vectorScale(wd->vOrient,wd->dAngSpeed,c1->w);
			}
		*piOutcome = MERGE;
		}
	else { /* bounce */
		switch (wd->iType) {
		case WallPlane:
		case WallTriangle:
		case WallRectangle:
		case WallDisk:
			wallBounceFlatFace(w->iWallID,wd,w->vOscVel,CP,c1,dt);
			break;
		case WallCylinderInfinite:
		case WallCylinderFinite:
			wallBounceCylinder(w->iWallID,wd,w->vOscVel,CP,c1,dt);
			break;
		case WallShell:
			/*wallBounceShell(w->iWallID,wd,w->vOscVel,CP,c1,dt);*/
			assert(0);/*DEBUG! (shell bounce not implemented yet)*/
			break;
		default:
			assert(0);
			}
#ifdef DEBUG_WALLS
		if (DEBUG_TRACE(c1->id.iOrder))
			fprintf(stderr,"%i & wall %i bounce after:  r %g %g %g v %g %g %g\n",c1->id.iOrder,iWallID,c1->r[0],c1->r[1],c1->r[2],c1->v[0],c1->v[1],c1->v[2]);
#endif
		for (i=0;i<3;i++)
			c1->r[i] -= c1->v[i]*dt; /* move particle back to start of step */
		*piOutcome = BOUNCE;
		}
	}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wallsio.h"
#include "aggs.h"
#include "collision.h"

#ifdef DEM_WALLS_REACT

/*** derivatives for aggregate spin and orientation updates ***/

int wallsEulerDerivs(FLOAT t,const FLOAT vars[],FLOAT derivs[],void *wd_as_void)
{
	/* t unused */
	FLOAT *torque = ((WALL_DATA *) wd_as_void)->torque;
	FLOAT *moments = ((WALL_DATA *) wd_as_void)->moments;

	/* omega[0, 1, 2] */
	derivs[0] = (torque[0] + vars[1]*vars[2]*(moments[1] - moments[2]))/moments[0];
	derivs[1] = (torque[1] + vars[2]*vars[0]*(moments[2] - moments[0]))/moments[1];
	derivs[2] = (torque[2] + vars[0]*vars[1]*(moments[0] - moments[1]))/moments[2];
 
	/* q1[0, 1, 2] */
	derivs[3] = vars[2]*vars[6] - vars[1]*vars[9];
	derivs[4] = vars[2]*vars[7] - vars[1]*vars[10];
	derivs[5] = vars[2]*vars[8] - vars[1]*vars[11];
 
	/* q2 */
	derivs[6] = vars[0]*vars[9]  - vars[2]*vars[3];
	derivs[7] = vars[0]*vars[10] - vars[2]*vars[4];
	derivs[8] = vars[0]*vars[11] - vars[2]*vars[5];
 
	/* q3 */
	derivs[9]  = vars[1]*vars[3] - vars[0]*vars[6];
	derivs[10] = vars[1]*vars[4] - vars[0]*vars[7];
	derivs[11] = vars[1]*vars[5] - vars[0]*vars[8];
    
    return 0;
	}

#endif /* DEM_WALLS_REACT */

#endif /* WALLS */
