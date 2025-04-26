#ifndef DEM_HINCLUDED
#define DEM_HINCLUDED

#ifdef DEM

#include "linalg.h"

/*
** For parallel jobs, there can be problems with the cache size due to
** the relatively large particle struct sizes.  The demand on the
** cache size will be affected by MAX_NUM_OVERLAPS_PER_PARTICLE (and
** the walls-specific value, if applicable), defined below.  If a
** value needs to be increased, if more parameters are added to the
** DEM_PARAMS struct, or if the number of particles per processor is
** high, consider increasing MDL_CACHE_SIZE in mpi/mdl.h (and maybe
** pthread/mdl.h).
*/

#define MAX_NUM_OVERLAPS_PER_PARTICLE 12

#ifdef WALLS
#define MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS 3
#endif

#ifdef DEM_WALLS_REACT
#define MAX_NUM_WALL_ASSEMBLIES 2 /* if changed, also change in wallsio.h */
#endif

typedef struct {
  int iOrder;
  double vShear[3];
  double vnOld[3];
#ifdef DEM_ROTATION_DASHPOT
  double vRoll[3];   /* spring-dashpot-slider rolling model */
  double dTwist;     /* spring-dashpot-slider twisting model */
#endif /* DEM_ROTATION_DASHPOT */
#if defined(__APPLE__) && defined(__LP64__) /* to prevent compiler warning on MacOS X */
  int liOverlapCounter;
#else
  long int liOverlapCounter;
#endif
#ifdef DEM_DIAG
  double vFn[3]; /* normal force due to contact */
  double vFt[3]; /* tangential force due to contact */
#endif /* DEM_DIAG */
  } DEM_ELEMENT;

#ifdef DEM_COHESION
typedef struct {
  int nCohesiveCoeff;
  int nConnetCoeff;
  int *iColor;
  double *dCohesiveCoeff;
  } DEM_COHESION_VAR_PARAMS;
#endif /* DEM_COHESION */

typedef struct {
  /* dEpsN, dEpsT, and dDelta are same as in parameter file (for now) */
  double dEpsN;
  double dEpsT;
  double dDelta;
  /* user-supplied parameters */
  double dKn;
  double dKt;
  double dMuS;
  double dMuR;
  double dMuT;
  double dAccCrit;
  double dMinorFrac,dMajorFrac,dErrorFrac;
  int bReadDEMData;
  int bReadOrient;
  int nDEMOutputs;
  int iDEMStatsInterval;
#ifdef DEM_FIXED_BALL
  int iFixedBallOption;
#endif /* DEM_FIXED_BALL */
#ifdef DEM_ROTATION_DASHPOT
  double dBetaS; /* shape parameter for spring-dashpot-slider rotation model */
#endif /* DEM_ROTATION_DASHPOT */
  int bUseLegacyRoll; /* use nonuniform torques for rolling friction? */
  int bUseContactRadius; /* use contact radius instead of reduced radius? */
#ifdef DEM_COHESION
  int iCohesionModel; /* 0=microscopic contact model, 1=granular bridge model */
  double dCohesiveCoeff; /* the cohesive coefficient */
  DEM_COHESION_VAR_PARAMS sCohesionVarParams; /* particles have various cohesion */
#endif /* DEM_COHESION */
#ifdef DEM_TWOLAYERS
  double dKnOuter;
  double dKtOuter;
  double dKnInner;
  double dKtInner;
  double dInnerOverlapBndry;
#endif /* DEM_TWOLAYERS */
  double dTangentialSpringDrag; /*DEBUG: might be timestep dependent parameter - consider revising */
  /*DEBUG
  ** could make dTangentialSpringDrag a rate of decay (e.g., see note
  ** in master.c) or binary (i.e., no decay, instant decay)
  */
  } DEM_PARAMS;

/*
** Header elements for .dem files:
**  time (double)
**  n_data (int)
**  MAX_NUM_OVERLAPS_PER_PARTICLE (int)
**  MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS (int), -1 if WALLS not defined
**  flag (int), 1 if DEM_ROTATION_DASHPOT defined, else 0
*/

#define DEMHEAD_SIZE (sizeof(double) + 4*sizeof(int))

/*
** Data elements for .dem files:
**  vPred[3] (doubles)
**  wPred[3] (doubles)
**  and MAX_NUM_OVERLAPS_PER_PARTICLE repetitions of:
**   iOrder (int)
**   vShear[3] (doubles)
**   vnOld[3] (doubles)
**   [DEM_ROTATION_DASHPOT only] vRoll[3] (doubles)
**   [DEM_ROTATION_DASHPOT only] dTwist (double)
**   liOverlapCounter (long int but xdr_long() only does 4 bytes)
**  [walls only] and MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS repetitions of:
**   iOrder (int)
**   vShear[3] (doubles)
**   vnOld[3] (doubles)
**   [DEM_ROTATION_DASHPOT only] vRoll[3] (doubles)
**   [DEM_ROTATION_DASHPOT only] dTwist (double)
**   liOverlapCounter (long int but xdr_long() only does 4 bytes)
*/

#ifdef WALLS

#ifdef DEM_ROTATION_DASHPOT
#define DEMDATA_SIZE (6*sizeof(double) + MAX_NUM_OVERLAPS_PER_PARTICLE*(sizeof(int) + 10*sizeof(double) + 4) + MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS*(sizeof(int) + 10*sizeof(double) + 4))
#else
#define DEMDATA_SIZE (6*sizeof(double) + MAX_NUM_OVERLAPS_PER_PARTICLE*(sizeof(int) + 6*sizeof(double) + 4) + MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS*(sizeof(int) + 6*sizeof(double) + 4))
#endif /* DEM_ROTATION_DASHPOT */

#else

#ifdef DEM_ROTATION_DASHPOT
#define DEMDATA_SIZE (6*sizeof(double) + MAX_NUM_OVERLAPS_PER_PARTICLE*(sizeof(int) + 10*sizeof(double) + 4))
#else
#define DEMDATA_SIZE (6*sizeof(double) + MAX_NUM_OVERLAPS_PER_PARTICLE*(sizeof(int) + 6*sizeof(double) + 4))
#endif /* DEM_ROTATION_DASHPOT */

#endif /* WALLS */

#define DEM_NUM_OVERLAP_BINS 100
#define DEM_NUM_COS_A_BINS   100
#define DEM_NUM_S_BINS       100

typedef struct {
  FLOAT fDistMin,fSpeedMax;
  int pOverlapHist[DEM_NUM_OVERLAP_BINS];
  int pCosAHist[DEM_NUM_COS_A_BINS];
  int pSHist[DEM_NUM_S_BINS];
} DEM_STATS;

typedef struct {
  FLOAT fDistMin,fSpeedMax,fOverlap,fCosA,fS2;
} DEM_STATS_PARTICLE;

#ifdef DEM_TIDAL_SPACE
/* specify marker iOrder numbers here; particles should not be coplanar */
#define IORDER_MARKER1 0
#define IORDER_MARKER2 62
#define IORDER_MARKER3 101
#define PLANET_COLOR 4 /* blue */
#endif

#ifdef DEM_TIDAL_LOCAL
#define PLANET_MASS 3.0034689055351672e-06 /* 1 M_Earth in M_Sun */
/* specify local frame origin & gravity vector here (use get_local.tcsh!) */
#define LOCAL_ORIG_X (-1.1224161397443008e-08)
#define LOCAL_ORIG_Y (-1.2960545208483053e-09)
#define LOCAL_ORIG_Z (1.8328978809586364e-09)
#define LOCAL_GRAV_X (0.092058741745686)
#define LOCAL_GRAV_Y (0.046172926339131)
#define LOCAL_GRAV_Z (-0.060803261269753)
#define LOCAL_SURF_X (-0.5)
#define LOCAL_SURF_Y (-0.4)
#define LOCAL_SURF_Z (-0.3)
#endif

#if defined(DEM_TIDAL_SPACE) || defined(DEM_TIDAL_LOCAL)
#define DEM_TIDAL_FILENAME "dem_tidal_acc.dat"
typedef struct {
  double dTime;
  double vPlanetPos[3];
  double vPlanetVel[3]; /* not used currently */
  double vAggPos[3]; /* starts at zero, but will drift */
  double vAggVel[3]; /* starts at zero, but will drift; not used */
  double vAggSpin[3]; /* in space frame */
  double vAggAcc[3]; 
  double vAggSpinDot[3]; /* in space frame; this is torque per unit mass */
  double vMarker1Pos[3];
  double vMarker1Vel[3];
  double vMarker2Pos[3];
  double vMarker2Vel[3];
  double vMarker3Pos[3];
  double vMarker3Vel[3];
  } DEM_TIDAL;
#endif

/* function prototypes */

int demCheckOverlapPoint(const Vector vRelPos,double dRadSq);
void demWipeOverlapElement(DEM_ELEMENT *e);
void demAssignParticleInertia(int p1stuck,int p2stuck,double mass1,double mass2,double rad1,double rad2,double rad1sq,double rad2sq,double *mass_inv1,double *mass_inv2,double *mom1,double *mom2,double *mom_inv1,double *mom_inv2);

#endif /* DEM */

#endif /* !DEM_HINCLUDED */
