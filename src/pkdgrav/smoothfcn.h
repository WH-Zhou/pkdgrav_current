#ifndef SMOOTHFCN_HINCLUDED
#define SMOOTHFCN_HINCLUDED

#include "pkd.h"
#include "floattype.h"

#ifdef STARFORM
#include "supernova.h"
#endif

#ifdef WALLS
#include "walls.h"
#endif

#include "foverride.h"

/* Benz Method (Default) */
#if !defined(PRES_MONAGHAN) && !defined(PRES_HK)
#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)
#endif

/* Monaghan Method */
#ifdef PRES_MONAGHAN
#define PRES_PDV(a,b) ((a+b)*0.5)
#define PRES_ACC(a,b) (a+b)
#endif

/* HK */
#ifdef PRES_HK
#define PRES_PDV(a,b) sqrt(a*b)
#define PRES_ACC(a,b) (sqrt(a*b)*2)
#endif

#define SMF_SMOOTHAGAIN    1

typedef struct smfParameters {
  double H;
  double a;
  double dDeltaAccelFac;
  double dSinkRadius;
  double dSinkBoundOrbitRadius;
  double dBHSinkEddFactor;
  double dBHSinkAlphaFactor;
  double dBHSinkFeedbackFactor;
  double dSinkCurrentDelta;
#ifdef GASOLINE
  double alpha;
  double beta;
  double gamma;
  double algam;
  int bGeometric;
  int bCannonical;
  int bGrowSmoothList;
#endif
  int bSinkThermal;
  int iSmoothFlags; /* Read/Write locally.  Master sets initial value. */
#if defined(STARFORM) || defined(CHECKSOFT)
  double dTime;
#endif
#ifdef STARFORM
  double dMinMassFrac;
  double dRadPreFactor;
  double dTimePreFactor;
  int bShortCoolShutoff;
  int bSNTurnOffCooling;
  int bSmallSNSmooth;
  double dSecUnit;
  double dGmUnit;
  struct snContext sn;
#endif    
#ifdef COLLISIONS
  double dTime;
  double dCentMass; /* for FindRejects() only */
  double dStart; /* collision search time interval */
  double dEnd;
  int bAllowSimulColl;
  double dCollapseLimit;
  int iOverlapOption;
  int bStrictOverlap;
  double dBackstepLimit;
  double dAdjPosLimit;
  int bOverlap;
#endif
#ifdef AGGS
  int bAggsSolveQuartic;
#endif
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
#ifdef WALLS
  WALL_PARAMS WP;
#endif
#ifdef RORY_EXTENSION
  double dDelta;
  double dMaxBinaryEcc;
#endif
#if defined(SPRINGS) || defined(DEM)
  FOVERRIDE_PARAMS FO;
#endif
  PKD pkd; /* useful for diagnostics, etc. */
  } SMF;

typedef struct nNeighbor {
	int iPid;
	int iIndex;
	PARTICLE *pPart;
	FLOAT fDist2;
	FLOAT dx;
	FLOAT dy;
	FLOAT dz;
	} NN;


enum smx_smoothtype {
  SMX_NULL,
  SMX_DENSITY,
  SMX_MARKDENSITY,
  SMX_MARKIIDENSITY,
  SMX_MARK,
  SMX_MEANVEL,
  SMX_DELTAACCEL,
  SMX_SINKACCRETE,
  SMX_BHSINKACCRETE,
#ifdef GASOLINE
  SMX_SPHPRESSURETERMS,
  SMX_DIVVORT,
  SMX_SHOCKTRACK,
  SMX_HKPRESSURETERMS,
  SMX_SPHPRESSURE,
  SMX_SPHVISCOSITY,
  SMX_HKVISCOSITY,
#ifdef STARFORM
  SMX_DIST_DELETED_GAS,
  SMX_DELETE_GAS,
  SMX_DIST_SN_ENERGY,
#endif
#ifdef SIMPLESF
  SMX_SIMPLESF_FEEDBACK,
#endif
#endif
#ifdef COLLISIONS
  SMX_REJECTS,
  SMX_COLLISION,
#ifdef SPRINGS
  SMX_ASSIGN_SPRINGS,
  SMX_COLOR_SPRINGS,
  SMX_MOVESPECIALS_SPRINGS,
  SMX_DO_SPRINGS,
#endif /* SPRINGS */
#ifdef DEM
  SMX_ASSIGN_DEM,
  SMX_DO_DEM,
  SMX_DEM_STATS,
#endif /* DEM */
  SMX_FINDBINARY,
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
  SMX_FIND_OVERLAPS,
#endif /* SLIDING_PATCH */
  };


/*  SMX_NULL */
void NullSmooth(PARTICLE *,int,NN *,SMF *);

/* SMX_DENSITY */
void initDensity(void *);
void combDensity(void *,void *);
void Density(PARTICLE *,int,NN *,SMF *);
void DensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARKDENSITY */
void initParticleMarkDensity(void *);
void initMarkDensity(void *);
void combMarkDensity(void *,void *);
void MarkDensity(PARTICLE *,int,NN *,SMF *);
void MarkDensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARKIIDENSITY */
void initParticleMarkIIDensity(void *);
void initMarkIIDensity(void *);
void combMarkIIDensity(void *,void *);
void MarkIIDensity(PARTICLE *,int,NN *,SMF *);
void MarkIIDensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARK */
void initMark(void *);
void combMark(void *,void *);

/* SMX_MEANVEL */
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

/* SMX_DELTAACCEL */
void DeltaAccel(PARTICLE *,int,NN *,SMF *);
void initDeltaAccel(void *);
void combDeltaAccel(void *,void *);

/* SMX_SINKACCRETE */
void SinkAccrete(PARTICLE *,int,NN *,SMF *);
void initSinkAccrete(void *);
void combSinkAccrete(void *,void *);

/* SMX_BHSINKACCRETE */
void BHSinkAccrete(PARTICLE *,int,NN *,SMF *);
void initBHSinkAccrete(void *);
void combBHSinkAccrete(void *,void *);

#ifdef GASOLINE

/* SMX_SPHPRESSURETERMS */
void initSphPressureTermsParticle(void *);
void initSphPressureTerms(void *);
void combSphPressureTerms(void *,void *);
void SphPressureTerms(PARTICLE *,int,NN *,SMF *);
void SphPressureTermsSym(PARTICLE *,int,NN *,SMF *);

/* SMX_DIVVORT */
void initDivVort(void *);
void combDivVort(void *,void *);
void DivVort(PARTICLE *,int,NN *,SMF *);
void DivVortSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SHOCKTRACK */
void initShockTrack(void *);
void combShockTrack(void *,void *);
void ShockTrack(PARTICLE *,int,NN *,SMF *);
void ShockTrackSym(PARTICLE *,int,NN *,SMF *);

/* SMX_HKPRESSURETERMS */
void initHKPressureTermsParticle(void *);
void initHKPressureTerms(void *);
void combHKPressureTerms(void *,void *);
void HKPressureTerms(PARTICLE *,int,NN *,SMF *);
void HKPressureTermsSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SPHPRESSURE */
void initSphPressureParticle(void *);
void initSphPressure(void *);
void combSphPressure(void *,void *);
void postSphPressure(PARTICLE *,SMF *);
void SphPressure(PARTICLE *,int,NN *,SMF *);
void SphPressureSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SPHVISCOSITY */
void initSphViscosityParticle(void *);
void initSphViscosity(void *);
void combSphViscosity(void *,void *);
void SphViscosity(PARTICLE *,int,NN *,SMF *);
void SphViscositySym(PARTICLE *,int,NN *,SMF *);

/* SMX_HKVISCOSITY */
void initHKViscosityParticle(void *);
void initHKViscosity(void *);
void combHKViscosity(void *,void *);
void HKViscosity(PARTICLE *,int,NN *,SMF *);
void HKViscositySym(PARTICLE *,int,NN *,SMF *);

#ifdef STARFORM

/* SMX_DIST_DELETED_GAS */
void initDistDeletedGas(void *p1);
void combDistDeletedGas(void *p1,void *p2);
void DistDeletedGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DELETE_GAS */
void DeleteGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DIST_SN_ENERGY */
void initTreeParticleDistSNEnergy(void *p1);
void initDistSNEnergy(void *p1);
void combDistSNEnergy(void *p1,void *p2);
void DistSNEnergy(PARTICLE *p, int, NN *, SMF *);
void postDistSNEnergy(PARTICLE *p1, SMF *smf);

#endif

#ifdef SIMPLESF
/* SMX_SIMPLESF_FEEDBACK */
void initSimpleSF_Feedback(void *p1);
void combSimpleSF_Feedback(void *p1,void *p2);
void SimpleSF_Feedback(PARTICLE *, int, NN *, SMF *);

#endif

#endif

#ifdef COLLISIONS

/* SMX_REJECTS */
void initFindRejects(void *p);
void combFindRejects(void *p1,void *p2);
void FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

/* SMX_COLLISION */
void CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

#ifdef SPRINGS

/* SMX_SPRINGS */
void AssignSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void ColorSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void MoveSpecialSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void initDoSpringsParticle(void *p);
void initDoSprings(void *p);
void combDoSprings(void *p1,void *p2);
void postDoSprings(PARTICLE *p,SMF *smf);
void DoSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
/* from springs.c (not SMX, should move??) */
void springsAssignColor(PARTICLE *p,int nSmooth,NN *nnList,int ALLOW_COLOR_CHANGE);
void springsMoveSpecials(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void springsForceSpecialsTensile(PARTICLE *p);
void springsForceSpecialsShear(PARTICLE *p);
void springsForceSpecialsTwist(PARTICLE *p);
void springsForceSpecialsRadial(PARTICLE *p);

#endif /* SPRINGS */

#ifdef DEM

/* SMX_DEM */
void AssignDEM(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void initDoDEMParticle(void *p);
void initDoDEM(void *p);
void combDoDEM(void *p1,void *p2);
void postDoDEM(PARTICLE *p,SMF *smf);
void DoDEM(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void DEMStats(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
#endif /* DEM */

void FindBinary(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

/* SMX_FIND_OVERLAPS */
void initFindOverlaps(void *p);
void combFindOverlaps(void *p1,void *p2);
void FindOverlaps(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

#endif /* SLIDING_PATCH */

#endif
