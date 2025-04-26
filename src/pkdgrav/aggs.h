#ifndef AGGS_HINCLUDED
#define AGGS_HINCLUDED

#ifdef AGGS

/*
** aggs.h
**
** Author: Kenneth W. Flynn
**         flynnk@astro.umd.edu
** Mods:   Derek C. Richardson
**         dcr@astro.umd.edu
**
** Modified: 01/28/01; DCR: 07/10/02, 05/29/03, lots since.
**
** This is the main header file for the aggregate expansion to
** pkdgrav.  This file and other associated changes extend the pkdgrav
** code to handle one or more aggregate structures.  Each aggregate is
** a collection of particles that do not mutually interact but that
** can interact with outside particles via collisions and/or gravity
** torques.  Rotation of the aggregate is computed using Euler's
** equations.
**
** To enable this functionality, define the AGGS flag.
**
** Aggregates themselves require few changes to the original pkdgrav.
** During the initial setup, the aggregrate structures defined below
** need to be filled out.  Calculation of gravity is unaffected, we
** merely make aggregrate particles "inactive" during a kick, making
** use of pre-existing functionality.  During drift, we update
** particle positions in the aggregate "by hand" under the influence
** of their rotation, computed via the Euler equations.  Collision
** handling is complicated by aggregate rotation, and some
** approximations need to be made.
**
** Combining AGGS with periodic boundary conditions however requires
** significant changes to the original pkdgrav, mostly because of the
** possibly of aggregates straddling boundaries.
*/
 
#include "linalg.h"
#include "pkd.h" /* for PARTICLE struct (etc.) */

#define AGGS_NO_COLOR_CHANGE FALSE /* set this to TRUE to prevent particles from changing color in aggs */

#define AGGS_INIT_BUFF_SIZE 16 /* allocate this many aggs minimum */

/* Constant in I = 2/5 mR^2 (for spheres).  Change for other bodies. */
#define AGGS_PARTICLE_INERTIA_PREFACTOR 0.4

/* aggregate macro prototypes */
int IS_AGG(const PARTICLE *p);
int AGG_IDX(const PARTICLE *p);
int AGG_SET_IDX(PARTICLE *p,int IDX);

/* particle is an aggregate if its original index is negative */
#define IS_AGG(p) ((p)->iOrgIdx < 0)

/* aggregate index = -1 - (particle original index) */
#define AGG_IDX(p) (-1 - (p)->iOrgIdx)

#define AGG_SET_IDX(p,IDX) {assert((IDX) >= 0); (p)->iOrgIdx = -1 - (IDX);}

/***************************** Structures *****************************/

/* Main structure that defines the aggregate. */
typedef struct
{
  /* Flag showing whether any particles are assigned to this agg. */
  int bAssigned;

  /* Mass of the aggregate. */
  FLOAT mass;

  /* Maximum "radius" of the aggregate (used for strength calculations). */
  FLOAT rad_max;

  /* Effective radius of the aggregate (ditto). */
  FLOAT rad_eff;

  /* COM position of the aggregate. */
  Vector r_com;

  /* COM velocity of the aggregate. */
  Vector v_com;

  /* COM acceleration of the aggregate. */
  Vector a_com;

  /* Torques acting on the aggregate in the body frame. */
  Vector torque;
 
  /* Angular velocity of the aggregate in the body frame. */
  Vector omega;
 
  /*
  ** Moments of inertia for the aggregate.  Note this is not a spatial
  ** vector, but a collection of three values.
  */
  Vector moments;
 
  /*
  ** Transformation matrix from body to space frame, called lambda.
  ** This is simply a matrix whose columns are the principal axes of
  ** inertia. A call to matrixTranspose() will yield the inverse
  ** conversion.
  */
  Matrix lambda;

  /* Time of last update (within step interval) for "drifting" */

  double dLastUpdate;

#ifdef AGGS_IN_PATCH
  double dPy_com; /* for the symplectic patch integrator */
#endif

  } Aggregate;

/* strength parameters (see collision.h) */

typedef struct {
  double dTensileCoef;
  double dTensileExp;
  double dShearCoef;
  double dShearExp;
  } STRENGTH_PARAMS;

/***************************** Functions *****************************/

void pkdAggsFind(PKD pkd,int *iMaxIdx,int *iMaxOrg);

void pkdAggsCountPart(PKD pkd,int iAggIdx,int *nPart);

#ifdef AGGS_IN_PATCH
void pkdAggsDelete(PKD pkd,int iAggIdx,double dStepTime,double dEventTime,int *bFound);
#else
void pkdAggsDelete(PKD pkd,int iAggIdx,int *bFound);
#endif

void pkdAggsMerge(PKD pkd,int iOldIdx,int iNewIdx);

void pkdAggsBackDrift(PKD pkd,int iAggIdx,double dt);
#ifdef FAST_AGGS
void pkdAggsGetRef(PKD pkd,int iAggIdx,Vector r_ref);

void pkdAggsGetPeriodicCOM(PKD pkd,int iAggIdx,const Vector r_ref,Scalar *m,Vector mr,Vector mv);

void pkdAggsSetUnwrapped(PKD pkd,int iAggIdx,const Vector r_com);
#endif /* FAST_AGGS */
#ifdef AGGS_IN_PATCH
void pkdAggsGetCOM(PKD pkd,int iAggIdx,const PATCH_PARAMS *PP,Scalar *m,Vector mr,Vector mv);
#else
void pkdAggsGetCOM(PKD pkd,int iAggIdx,Scalar *m,Vector mr,Vector mv);
#endif

void pkdAggsGetAxesAndSpin(PKD pkd,int iAggIdx,const Vector r_com,
						   const Vector v_com,Matrix I,Vector L,
						   Scalar *rad_max_sq,Scalar *volume);

void pkdAggsSetBodyPos(PKD pkd,int iAggIdx,Matrix spaceToBody);

#ifdef AGGS_IN_PATCH
void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,const Vector r_com,Matrix lambda,double dTime,const PATCH_PARAMS *PP);
#else
void pkdAggsSetSpacePos(PKD pkd,int iAggIdx,const Vector r_com,Matrix lambda);
#endif

void pkdAggsSetSpaceVel(PKD pkd,int iAggIdx,const Vector v_com,const Vector a_com,const Vector omega,const Vector omegadot,Matrix lambda,double dDelta);

void pkdAggsSetSpaceSpins(PKD pkd,int iAggIdx,const Vector omega);

#ifndef DEM
void pkdAggsGetAccel(PKD pkd,int iAggIdx,Scalar *m,Vector ma);

void pkdAggsCheckStress(PKD pkd,int iAggIdx,
						Scalar mass,Scalar rad_max,Scalar rad_eff,
						const Vector r_com,const Vector a_com,
						const Vector omega,STRENGTH_PARAMS *SP,
						int *nLost,int *nLeft);

void pkdAggsGetTorque(PKD pkd,int iAggIdx,const Vector r_com,const Vector a_com,
					  Vector torque);
#else /* DEM */
void pkdAggsGetAccelAndTorqueDEM(PKD pkd,int iAggIdx,const Vector r_com,Scalar *m,
					  Vector ma, Vector torque);
#endif /* DEM */

#ifdef AGGS_IN_PATCH
void pkdAggsInPatchGetRef(PKD pkd,int iAggIdx,Scalar *m_max,Vector r_max);
void pkdAggsInPatchGetUnwrapped(PKD pkd,int iAggIdx,const Vector r_ref,double dTime,const PATCH_PARAMS *PP);
void pkdAggsInPatchOffset(PKD pkd,int iAggIdx,FLOAT dx,FLOAT dy,FLOAT dvy,int bDoUnwrapped);
void aggsRotateForPatch(Vector w,double dDirection);
#endif /* AGGS_IN_PATCH */

/** Function that returns the values for the derivatives of the Euler 
 *   equations of motion.  For more information, refer to the accompanying
 *   document.
 *
 *  Parameters (in):
 *   t - Time to evaluate at (ignored)
 *   vars - Current value of variables being integrated.  In detail:
 *            0 = omega[0]
 *            1 = omega[1]
 *            2 = omega[2]
 *            3 = q1[0]
 *            4 = q1[1]
 *            5 = q1[2]
 *            6 = q2[0]
 *            7 = q2[1]
 *            8 = q2[2]
 *            9 = q3[0]
 *           10 = q3[1]
 *           11 = q3[2]
 *   agg_as_void - The aggregate structure, used to obtain the principal 
 *                  moments of inertia
 *
 *  Parameters (out):
 *   derivs - Values for the derivatives of the aforementioned variables
 */
int aggsEulerDerivs(FLOAT t,const FLOAT vars[],FLOAT derivs[],void* agg_as_void);

/*
** Integrate aggregate spin and orientation over interval dt.
**
** Parameters (in):
**  agg - Pointer to aggregate
**  dt - Integration interval (time step)
**
** Parameters (out):
**  agg - Spin and orientation updated
*/
void aggsRungeAdvance(Aggregate *agg,double dDelta,double dt);

#endif /* AGGS */

#endif /* !AGGS_HINCLUDED */
