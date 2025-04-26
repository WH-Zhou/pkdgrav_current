#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#ifdef COLLISIONS

#include "pkd.h" /* for PARTICLE struct */
#include "ssio.h" /* for SSDATA struct */

/*#define SPECIAL_PARTICLES*/ /* should be added as a compile option */

#ifdef SPECIAL_PARTICLES
#include "special.h"
#endif

#ifdef AGGS
#include "aggs.h" /* for Aggregate struct */
#endif

#ifdef RUBBLE_ZML
#include "rubble.h"
#elif defined(COLLMOD_ZML)
#include "collmod.h"
#endif

#ifdef CHARGE
typedef struct {
  double dQ; /* total charge */
  double dZ; /* first moment of charge distribution (z-component) */
  double dSigZ; /* second moment of charge distribution (z-component) */
  double dArea; /* area of the cylinder end caps -- this is a hack */
  } CHARGE_PARAMS;
void pkdChargeZGetMoments(PKD pkd,CHARGE_PARAMS *CP);
void pkdChargeZApplyMoments(PKD pkd,const CHARGE_PARAMS *CP);
#endif /* CHARGE */

/* define PARTICLE_ID and COLLIDER before including walls... */

typedef struct {
	int iPid;
	int iOrder;
	int iIndex;
	int iOrgIdx;
	} PARTICLE_ID;

typedef struct {
  PARTICLE_ID id;
  FLOAT fMass;
  FLOAT fRadius;
  FLOAT r[3];
  FLOAT v[3];
  FLOAT w[3];
  FLOAT a[3];
  FLOAT dt;
  int iColor;
  int iRung;
  int bTinyStep;
#ifdef AGGS
  Aggregate agg;
#endif
#ifdef WALLS
  FLOAT omega_v[3]; /* needed for particles stuck on rotating cylinders */
#endif
#ifdef SLIDING_PATCH
  double dPy;
#endif
#ifdef ORIGIN_HISTOGRAM
	FLOAT origin_bins[NUM_ORIGIN_BINS];
#endif
#ifdef COLLMOD_ZML
  int bFrag;
#endif
	} COLLIDER;

#ifdef WALLS
#include "walls.h"
#endif

FLOAT RADIUS(PARTICLE *p);
#define RADIUS(p) ((p)->fSoft) /* currently PARTICLE has no radius field;
								  to ensure all interactions are strictly
								  Newtonian, set the softening length no
								  larger than the physical radius. */

int REJECT(PARTICLE *p);
#define REJECT(p) ((p)->dtCol < 0.0)

int COLLISION(double t);
#define COLLISION(t) ((t) < DBL_MAX)

unsigned int BIT(unsigned int n);
#define BIT(n) (1 << (n))

#define COLL_LOG_NONE		0
#define COLL_LOG_VERBOSE	1
#define COLL_LOG_TERSE		2

#define MISS	      0
#define MERGE	      BIT(0)
#define BOUNCE	      BIT(1)
#define FRAG	      BIT(2)

#ifdef RORY_EXTENSION
#define BINARY_MERGE  BIT(3)
#endif

#ifdef COLLMOD_ZML
#define EXPLODE BIT(4)
#endif

#ifdef WALLS
#define DEATH (-1)
#endif

#define MAX_NUM_FRAG 100	/* PJC: 14/07/14 Must be >4 if you will generate fragment chains with more than 4 members. */

enum {ConstEps,PowerLaw,Compacted,Borderies}; /* EpsN options */

enum {EscVel,MaxTrv}; /* slide options */

enum {OverlapIgnore=-1,OverlapIsError,OverlapBackstep,OverlapAdjPos,OverlapRepel,OverlapMerge}; /* overlap strategies */

typedef struct {
  int iOutcomes;
  double dMergeLimit;
  double dDensity;
  int iDensityAltCol;
  double dDensityAltVal;
  int iEpsNOption;
  double dEpsN;
  double dEpsT;
  double dEpsNCoef;
  double dEpsNExp;
  double dEpsNVStar;
  double dEpsNMin;
  int iSlideOption;
  double dSlideLimit;
  double dSlideLimit2;
  double dSlideEpsN;
  double dSlideEpsT;
  double dCollapseLimit;
  double dCollapseEpsN;
  double dCollapseEpsT;
  double dCrushLimit;
  double dCrushEpsN;
  double dCrushEpsT;
  int iOverlapOption;
  int bStrictOverlap;
  double dBackstepLimit;
  double dAdjPosLimit;
  double dRepelFac;
  double dFragLimit;
#ifdef RUBBLE_ZML
  int iRubForcedOutcome;
  int iRubColor;
  double dRubbleMinFracMass;
  int bDoRubbleKDKRestart;
  DUST_BINS_PARAMS DB;
#endif
#ifdef COLLMOD_ZML
	int bDoCollModKDKRestart;
	DUST_BINS_PARAMS DB;
#endif
#ifdef AGGS
  int bAggsSolveQuartic;
  STRENGTH_PARAMS SP;
#endif
#ifdef WALLS
  WALL_PARAMS WP;
#endif
#ifdef CHARGE
  CHARGE_PARAMS CP;
#endif
} COLLISION_PARAMS;

/* handy macros */

FLOAT SOFT(COLLIDER *c);
#define SOFT(c) ((c)->fRadius)

FLOAT SOFT_FROM_SSDATA(SSDATA *d);
#define SOFT_FROM_SSDATA(d) ((d)->radius) /* used in pkd.c */

#ifdef AGGS /* do these here to get around self-ref loop (COLLIDER needs Aggregate) */
int COLLIDER_IS_AGG(COLLIDER *c);
#define COLLIDER_IS_AGG(c) ((c)->id.iOrgIdx < 0)

int COLLIDER_AGG_IDX(COLLIDER *c);
#define COLLIDER_AGG_IDX(c) (-1 - (c)->id.iOrgIdx)
#endif /* AGGS */

int ALLOW_COLOR_CHANGE(PARTICLE *p);
/*
** This macro does not handle all cases.  For example, if particle
** merging is enabled, and iDensityAltCol is non-zero, you probably
** want added restrictions on color changes.
*/
#if defined(AGGS) && defined(WALLS)
#define ALLOW_COLOR_CHANGE(p) (!AGGS_NO_COLOR_CHANGE && !PARTICLE_STUCK(p))
#elif defined(AGGS)
#define ALLOW_COLOR_CHANGE(p) (!AGGS_NO_COLOR_CHANGE) /* p not used */
#elif defined(WALLS)
#define ALLOW_COLOR_CHANGE(p) (!PARTICLE_STUCK(p))
#elif defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
#define ALLOW_COLOR_CHANGE(p) (FALSE)
#elif defined(RORY_EXTENSION)
#define ALLOW_COLOR_CHANGE(p) (FALSE)
#elif defined(DEM)
#define ALLOW_COLOR_CHANGE(p) (FALSE)
#else
#define ALLOW_COLOR_CHANGE(p) (TRUE)
#endif

/* function prototypes */

void pkdNextCollision(PKD pkd, double *dtCol, int *iOrder1, int *iOrder2);
FLOAT pkdApplyBCs(PKD pkd,double dTime,const COLLIDER *c1,COLLIDER *c2,FLOAT fOffset[]);
void pkdGetColliderInfo(PKD pkd, int iOrder, COLLIDER *c);
void pkdPutColliderInfo(PKD pkd,const COLLIDER *c,int iOrder2,PARTICLE *p,double dt);
void pkdDoCollision(
					PKD pkd, double dt, const COLLIDER *c1, const COLLIDER *c2,
					int bPeriodic, int iStartStep, double dDelta, double dTime,
					const COLLISION_PARAMS *CP, int *piOutcome,double *dT,
					COLLIDER *cOut, int *pnOut
#ifdef COLLMOD_ZML
					,const DUST_BINS_PARAMS *DB,int *iDustBin,DustBins *DustBin
#endif
					);
#ifdef AGGS
#ifdef SLIDING_PATCH
void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *c1,const COLLIDER *c2,
						int bPeriodic,const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,COLLIDER *cOut,int *pnOut,
						FLOAT *dx,FLOAT *dy,FLOAT *dvy);
#else
void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *c1,const COLLIDER *c2,
						int bPeriodic,const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,COLLIDER *cOut,int *pnOut);
#endif
#endif /* AGGS */
#ifdef WALLS
void pkdWallsDoCollision(PKD pkd,const COLLISION_PARAMS *CP,COLLIDER *c1,const COLLIDER *c2,double dt,int *piOutcome);
#endif
void pkdResetColliders(PKD pkd, int iOrder1, int iOrder2);

#ifdef COLLMOD_ZML
void pkdPutColliderInfoSC(PKD pkd,const COLLIDER *c,int iOrder2,int iOrder3,PARTICLE *p,double dt);
#endif

#ifdef RORY_EXTENSION
double LastKickTime(int iRung, double dBaseStep, double dTimeNow);
void SetMergerRung(const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
				   double dBaseStep,double dTimeNow,int iTime0);
void pkdFindTightestBinary(PKD pkd,double *dBindEn,int *iOrder1,int *iOrder2,
						   int *n);
void pkdMergeBinary(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
					int bPeriodic,double dBaseStep,double dTimeNow,int iTime0,
					double dDensity,int *bool);
void pkdSetBall(PKD pkd,double dDelta,double fac);
#ifdef SLIDING_PATCH
#define MAXLARGEMASS 25	 /* Maximum number of particles to randomize */
#define MAXNEIGHBORS 250 /* Maximum neighbors surrounding a large particle */
void pkdFindLargeMasses(PKD,double,double,double,double,PARTICLE *p,double *,int *);
void pkdGetNeighborParticles(PKD,double *,double,int,double,PARTICLE *p,double *,int *);
#endif /* SLIDING_PATCH */
#endif /* RORY_EXTENSION */

#ifdef ORIGIN_HISTOGRAM
void MergeHistograms(FLOAT m1, FLOAT* bins1, FLOAT m2, FLOAT const* bins2);
void pkdInitializeOriginBins(PKD pkd,const DUST_BINS_PARAMS *DB);
#endif /* ORIGIN_HISTOGRAM */

#endif /* COLLISIONS */

#endif /* !COLLISION_HINCLUDED */
