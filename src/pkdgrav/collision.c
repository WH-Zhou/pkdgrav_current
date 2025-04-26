/* suppress pedantic "empty translation unit" warning if !COLLISIONS */
typedef int make_iso_compilers_happy;

#ifdef COLLISIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "collision.h"
#include "ssdefs.h" /* for PLANETESIMAL */

#ifdef RUBBLE_ZML
#include "rubble.h"
#elif defined(COLLMOD_ZML)
#include "collmod.h"
#endif

#define BOUNCE_OK 0
#define NEAR_MISS 1

#ifdef RUBBLE_ZML	
/* here instead of rubble.h because it uses COLLIDER... */
extern void rubRubbleCollide(const COLLIDER *col_a,const COLLIDER *col_b,
							 const COLLISION_PARAMS *cp,COLLIDER **col_c,
							 double *dMassInDust,int *nOut);
#endif

#ifdef COLLMOD_ZML
/* here instead of collmod.h because it uses COLLIDER... */
extern void collCollModCollide(const COLLIDER *col_a,const COLLIDER *col_b,
							   const COLLISION_PARAMS *cp,COLLIDER **col_c,
							   double *dMassInDust,int *nOut, PKD pkd);
#endif

void pkdNextCollision(PKD pkd,double *dt,int *iOrder1,int *iOrder2)
{
	/*
	 ** Returns time and iOrder of particles for earliest predicted
	 ** collision on this processor.  Argument initialization must
	 ** occur in calling routine.
	 */

	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue; /* skip over inactive particles */
		if (p->iOrder < 0) continue; /* skip over deleted particles */
		if (p->dtCol < *dt) {
			*dt = p->dtCol;
			*iOrder1 = p->iOrder;
			*iOrder2 = p->iOrderCol;
			}
		}
	}

#ifdef RUBBLE_ZML
void pkdMergerReorder(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,const COLLIDER *c,COLLIDER *cOut,double dt,int bDiagInfo,const COLLISION_PARAMS *CP,double dMassInDust)
#else
void pkdMergerReorder(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,const COLLIDER *c,COLLIDER *cOut,double dt,int bDiagInfo)
#endif
{
	/*
	** NOTE:
	**   c1 = (input) pointer to 1st particle's pre-collision info
	**   c2 = (input) pointer to 2nd particle's pre-collision info
	**    c = (input) pointer to array of post-collision particle
	**        info, adjusted to start of step
	** cOut = (output) pointer to array of post-collision particle
	**        info, intended for collision log
	*/

	const PARTICLE_ID *pMrg=NULL,*pDel=NULL,*pOth=NULL;

	assert(c1 != NULL);
	assert(c2 != NULL);
	assert(c != NULL);
	assert(cOut != NULL);
#ifdef RUBBLE_ZML
	assert(CP != NULL);
#endif

	if (c1->id.iPid == pkd->idSelf) { /* local particle */
		/*
		** Keep this particle if it has larger mass, or, if both
		** particles are the same mass, keep this particle if it has a
		** lower original index, or, if both particles have same
		** original index, keep this particle if it has a lower
		** processor number, or, if both particles are on same
		** processor, keep particle with lower iOrder.
		*/
		if (c1->fMass > c2->fMass ||
			(c1->fMass == c2->fMass &&
			 (c1->id.iOrgIdx < c2->id.iOrgIdx ||
			  (c1->id.iOrgIdx == c2->id.iOrgIdx &&
			   (c1->id.iPid < c2->id.iPid ||
				(c1->id.iPid == c2->id.iPid &&
				 c1->id.iOrder < c2->id.iOrder))))))
			pMrg = &c1->id; /* new centre-of-mass particle */
		else
			pDel = &c1->id; /* this particle is removed */
		pOth = &c2->id;
		}
	if (c2->id.iPid == pkd->idSelf) { /* same thing for other particle */
		if (c2->fMass > c1->fMass ||
			(c2->fMass == c1->fMass &&
			 (c2->id.iOrgIdx < c1->id.iOrgIdx ||
			  (c2->id.iOrgIdx == c1->id.iOrgIdx &&
			   (c2->id.iPid < c1->id.iPid ||
				(c2->id.iPid == c1->id.iPid &&
				 c2->id.iOrder < c1->id.iOrder))))))
			pMrg = &c2->id;
		else
			pDel = &c2->id;
		pOth = &c1->id;
		}
	/* 
	** Store or delete particle as appropriate, and make sure master
	** knows ID of surviving particle for tracking purposes.  INT_MAX
	** is used to indicate that the other collider is now gone.
	*/
	if (pMrg != NULL) {
		pkdPutColliderInfo(pkd,c,INT_MAX,&pkd->pStore[pMrg->iIndex],dt);
		if (bDiagInfo)
			cOut->id = *pMrg; /* struct copy */
#ifdef RUBBLE_ZML
		{
			/*DEBUG -- RUBBLE_ZML: CHECK FOR ABNORMAL DENSITY*/
			double rho;
			rho = c->fMass/(4.0/3.0*M_PI*CUBE(c->fRadius));
			/*	  if (rho > 4.0e6 && pkd->pStore[pMrg->iIndex].iColor == PLANETESIMAL)*/ /* f=1 */
			if (rho > 2.0e4 && pkd->pStore[pMrg->iIndex].iColor == PLANETESIMAL) /* f=6 */
				(void) printf("BAD MERGED DENSITY: iOrder=%i iOrgIdx=%i rho=%g\n",pMrg->iOrder,pMrg->iOrgIdx,rho);
			}
		if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {
			/*
			** Store any dust created by an interpolated collision
			** into particle dust parameter.
			*/
			pkd->pStore[pMrg->iIndex].dDustMass = dMassInDust; /* dealt with immediately after collision */
			}
#endif /* RUBBLE_ZML */
		}

	if (pDel != NULL) {
		pkdDeleteParticle(pkd,&pkd->pStore[pDel->iIndex]);
		if (bDiagInfo && pMrg == NULL)
			cOut->id = *pOth; /* may need merger info */
		}
	}

#ifdef RORY_EXTENSION

double LastKickTime(int iRung,double dBaseTime,double dTimeNow) 
{
	/*
	 ** Determines the last time a particle received a force update.
	 */

	double dRungTime;
	int nRungSteps;

	dRungTime = dBaseTime/(1<<iRung);
	nRungSteps = (int)(dTimeNow/dRungTime);

	return dTimeNow - (dRungTime*nRungSteps);
}

void SetMergerRung(const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
				   double dBaseStep,double dTimeNow,int iTime0)
{
  double lastkick1,lastkick2;
  FLOAT m1,m2,m;
  int k;

  m1 = c1->fMass;
  m2 = c2->fMass;
  m = m1 + m2;

  if (c1->iRung!=c2->iRung) {
    lastkick1=LastKickTime(c1->iRung,dBaseStep,(dTimeNow-iTime0));
    lastkick2=LastKickTime(c2->iRung,dBaseStep,(dTimeNow-iTime0));
    if (m1 > m2) {
      c->iRung=c1->iRung;
      if (lastkick1>lastkick2) {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m1/m)*c1->a[k]*(lastkick1-lastkick2);
	}
      } else {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m2/m)*c2->a[k]*(lastkick2-lastkick1);
	}
      }
    } else { /* c2->fMass > c1->fMass */
      c->iRung=c2->iRung;
      if (lastkick1>lastkick2) {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m1/m)*c1->a[k]*(lastkick1-lastkick2);
	}
      } else {
	for (k=0;k<3;k++) {
	  c->v[k]+=(m2/m)*c2->a[k]*(lastkick2-lastkick1);
	}
      }
    }
  } else { /* c1->iRung == c2->iRung */
    c->iRung=c1->iRung;
  }
  /* This was the old way of fixing the resultant particle's rung
   * c->iRung = (c2->fMass > c1->fMass ? c2->iRung : c1->iRung); */
}

void pkdFindTightestBinary(PKD pkd,double *dBindEn,int *iOrder1,int *iOrder2,int *n)
{
  /* Find the local binary with the highest binding energy. */

  PARTICLE *p;
  int i,nLocal = pkdLocal(pkd);

  for (i=0;i<nLocal;i++) {
    p = &pkd->pStore[i];
    if (!TYPEQueryACTIVE(p)) continue;
    if (p->iOrder<0) continue;
    if (p->dtCol < *dBindEn) { /* This particle is more tightly bound */
      *n=2;
      *dBindEn = p->dtCol;
      *iOrder1 = p->iOrder;
      *iOrder2 = p->iOrderCol;
    }
  }
}

void pkdMergeBinary(PKD pkd,const COLLIDER *pc1,const COLLIDER *pc2,COLLIDER *cOut,
					int bPeriodic,double dBaseStep,double dTimeNow,int iTime0,
					double dDensity,int *bool)
{
  /*DEBUG ignores iDensityAltCol and dDensityAltVal*/
  COLLIDER c1,c2,c;
  FLOAT rc1[3],rc2[3],vc1[3],vc2[3],ac1[3],ac2[3];
  FLOAT m1,m2,m;
  int k,bDiagInfo;
  FLOAT rsq1=0,rsq2=0,vsq1=0,vsq2=0,r2,r1,ang_mom,ke;
  const double dDenFac = 4.0/3*M_PI;  
  FLOAT fOffset[3];

  c1 = *pc1;
  c2 = *pc2;

  /* First find center of mass position and velocity */

  m1=c1.fMass;
  m2=c2.fMass;
  m=m1+m2;
  r1=c1.fRadius;
  r2=c2.fRadius;

  bDiagInfo = (c1.id.iPid == pkd->idSelf); 
 
  if (bPeriodic) 
	  (void) pkdApplyBCs(pkd,dTimeNow,&c1,&c2,fOffset);/*DEBUG fOffset not used--DCR*/

  for (k=0;k<3;k++) {
    c.r[k] = (m1*c1.r[k] + m2*c2.r[k])/m;
    c.v[k] = (m1*c1.v[k] + m2*c2.v[k])/m;
    c.a[k] = (m1*c1.a[k] + m2*c2.a[k])/m;
    c.w[k] = (m1*c1.w[k] + m2*c2.w[k])/m;/*DEBUG can't do this with spins, can you?--DCR*/
    rc1[k] = c1.r[k] - c.r[k];
    rc2[k] = c2.r[k] - c.r[k];
    vc1[k] = c1.v[k] - c.v[k];
    vc2[k] = c2.v[k] - c.v[k];
    ac1[k] = c1.a[k] - c.a[k];
    ac2[k] = c2.a[k] - c.a[k];
    rsq1+=rc1[k]*rc1[k];
    vsq1+=vc1[k]*vc1[k];
    rsq2+=rc2[k]*rc2[k];
    vsq2+=vc2[k]*vc2[k];
  }
  c.fMass=m;
  if (dDensity)
    c.fRadius = pow(m/(dDenFac*dDensity),1.0/3);
  else
    c.fRadius = (m2 > m1 ? r2*pow(m/m2,1.0/3) : r1*pow(m/m1,1.0/3));

  /* This much angular momentum and kinetic energy are lost. *//*DEBUG why is ang mom lost?--DCR*/
  ang_mom = m1*sqrt(rsq1)*sqrt(vsq1) + m2*sqrt(rsq2)*sqrt(vsq1);
  ke=0.5*(m1*vsq1 + m2*vsq2);

  /* For now both colliders are the same type */
  c.iColor = c1.iColor;

  SetMergerRung(&c1,&c2,&c,dBaseStep,dTimeNow,iTime0);
#ifdef RUBBLE_ZML
  assert(0); /* really shouldn't combine these options! */
#else
  pkdMergerReorder(pkd,&c1,&c2,(const COLLIDER *) &c,&c,-1.0,bDiagInfo);
#endif

  *cOut = c;

  /* For parallel purposes, let's now set the boolean to 1 */
  *bool=1;
}

void pkdSetBall(PKD pkd,double dDelta,double fac)
{
  PARTICLE *p;
  double x;
  int i,k,nLocal;

  nLocal = pkdLocal(pkd);

  for (i=0;i<nLocal;i++) {
    p = &pkd->pStore[i];
    x = 0.0; /* dummy variable: square of absolute speed */
    for (k=0;k<3;k++) 
      x += p->v[k]*p->v[k];
	/* search radius = factor * dt * v + R */
    x = fac*dDelta*sqrt(x) + RADIUS(p); /* dummy: search ball radius */
	/* make sure ball isn't bigger than periodic cell */
	for (k=0;k<3;k++)
	  assert(x < 0.5*pkd->fPeriod[k]);
	p->fBallMax = x; /* for expanding tree node size */
	p->fBall2 = x*x; /* for searching around particle */
	p->cpStart = 0; /* force smooth search to start at tree root */
  }
}

#ifdef SLIDING_PATCH

void pkdFindLargeMasses(PKD pkd,double dMass,double dCentralMass,double dHelio,double hill,PARTICLE *p,double *r,int *n)
{
    int i,nLocal = pkdLocal(pkd);

    *n=0;
    for (i=0;i<nLocal;i++) {
	if (pkd->pStore[i].fMass > dMass) {
	    p[*n]=pkd->pStore[i];
	    r[*n]=hill*(pow((pkd->pStore[i].fMass/(3*dCentralMass)),(1.0/3))*(dHelio+pkd->pStore[i].r[0])); /* This assumes z height is negligible */ 
	    *n+= 1;
	    assert(*n < MAXLARGEMASS);
	    }
	}	
    }

void pkdGetNeighborParticles(PKD pkd,double *r,double dRadius2,int id,double dTime,PARTICLE *p,double *sep2,int *n) 
{
    double dDist2,dist2;
    double dx,dx1,dy,dy1,dz,dz1,sx,sy,sz;
    int i,j,nLocal = pkdLocal(pkd);

    /* This subroutine assumes that it is part of the sliding
       patch. Although there may be future uses for this that would be
       more general. Because of the unusual geometry of the sliding
       patch (namely that the ghost cells themselves are shearing),
       this algorithm is surprisingly tricky, and is modelled after
       Derek's INTERSECT macro in smooth.h. --Rory 07/12/04 
    */

    *n=0;
    i=0;
    while (i < nLocal) {
	/* First make sure we're not dealing with the large mass
	   particle */
	if (pkd->pStore[i].iOrder == id) {
	    i++;
	    continue;
	    }
	
	dDist2=0;	
	dx = pkd->pStore[i].r[0] - r[0];
	dx1 = r[0] - pkd->pStore[i].r[0];
	sy = r[1];
	if (dx > 0.0 ) {
	    dx1 += pkd->fPeriod[0];
	    if (dx1 < dx) {
		dDist2 = dx1*dx1;
		sx = r[0] + pkd->fPeriod[0];
		sy += SHEAR(1,dTime,pkd->PP);
		if (sy >= 0.5*pkd->fPeriod[1]) sy -= pkd->fPeriod[1];
		else if (sy < -0.5*pkd->fPeriod[1]) sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 = dx*dx;
		sx = r[0];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	
	else if (dx1 > 0.0) {
	    dx += pkd->fPeriod[0];
	    if (dx < dx1) {
		dDist2 = dx*dx;
		sx =  r[0] - pkd->fPeriod[0];
		sy += SHEAR(-1,dTime,pkd->PP);
		if (sy >= 0.5*pkd->fPeriod[1]) sy -= pkd->fPeriod[1];
		else if (sy < -0.5*pkd->fPeriod[1]) sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 = dx1*dx1;
		sx = r[0];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	
	    }
	else {
	    dDist2 = 0.0;
	    sx = r[0];
	    }
	
	dy = pkd->pStore[i].r[1] - sy;
	dy1 = sy - pkd->pStore[i].r[1];
	if (dy > 0.0) {
	    dy1 += pkd->fPeriod[1];
	    if (dy1 < dy) {
		dDist2 += dy1*dy1;
		sy += pkd->fPeriod[1];
		}
	    else {
		dDist2 += dy*dy;
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	else if (dy1 > 0) {
	    dy += pkd->fPeriod[1];
	    if (dy < dy1) {
		dDist2 += dy*dy;
		sy -= pkd->fPeriod[1];
		}
	    else {
		dDist2 += dy1*dy1;
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}
	    }
	else {
	    }
	
	/* Is this necessary in sliding_patch? */
	dz = pkd->pStore[i].r[2] - r[2];
	dz1 = r[2] - pkd->pStore[i].r[2];
	if (dz > 0.0) {
	    dz1 += pkd->fPeriod[2];
	    if (dz1 < dz) {
		dDist2 += dz1*dz1;
		sz = r[2]+pkd->fPeriod[2];
		}
	    else {
		dDist2 += dz*dz;
		sz = r[2];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}
	    }
	else if (dz1 > 0.0) {
	    dz += pkd->fPeriod[2];
	    if (dz < dz1) {
		dDist2 += dz*dz;
		sz = r[2]-pkd->fPeriod[2];
		}
	    else {
		dDist2 += dz1*dz1;
		sz = r[2];
		}
	    if (dDist2 > dRadius2) {
		i++;
		continue;
		}	    
	    }
	else {
	    sz = r[2];
	    }
    
	/* The particle is within the sphere */
	assert(dDist2 < dRadius2);
    
	/* Just for safe keeping here is the version in the
	   non-sliding patch, non periodic case */
	dist2=0;	
	for (j=0;j<3;j++) 
	    dist2+= (pkd->pStore[i].r[j]-r[j])*(pkd->pStore[i].r[j]-r[j]);
	   /* if ((dist2 < radius*radius) && (dist2 != 0)) { */
	
	p[*n] = pkd->pStore[i];
	sep2[*n] = dDist2;	
	*n+= 1;
	assert(*n < MAXNEIGHBORS);
	i++;
	}    
    }

#endif /* SLIDING_PATCH */

#endif /* RORY_EXTENSION */

FLOAT pkdApplyBCs(PKD pkd,double dTime,const COLLIDER *c1,COLLIDER *c2,FLOAT fOffset[])
{
	/*
	 ** Determine whether second particle is supposed to be an image,
	 ** and adjust accordingly.  This is assumed to be the case if the
	 ** particle positions are separated by more than half a box
	 ** length.  This may fail for small box sizes or fast particles!
	 */

	int i;
	FLOAT fShear = 0.0;

	fOffset[0] = fOffset[1] = fOffset[2] = 0.0;

	for (i=0;i<3;i++) {
		if (pkd->fPeriod[i] == FLOAT_MAXVAL)
			continue; /* skip if no BCs in this direction */
		if (c2->r[i] - c1->r[i] >= 0.5*pkd->fPeriod[i]) {
			fOffset[i] -= pkd->fPeriod[i];
			c2->r[i] -= pkd->fPeriod[i];
			}
		else if (c2->r[i] - c1->r[i] < - 0.5*pkd->fPeriod[i]) {
			fOffset[i] += pkd->fPeriod[i];
			c2->r[i] += pkd->fPeriod[i];
			}
#ifdef SLIDING_PATCH
		if (i == 0 && pkd->PP->bPatch) {
			if (fOffset[0] < 0.0) {
				fShear = 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth;
				fOffset[1] = SHEAR(-1,dTime,pkd->PP);
				}
			else if (fOffset[0] > 0.0) {
				fShear = - 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth;
				fOffset[1] = SHEAR(1,dTime,pkd->PP);
				}
			c2->r[1] += fOffset[1];
			c2->v[1] += fShear;
			c2->dPy -= fShear/3.0;
			}
#endif
		}

	/* return the azimuthal shear in case it's needed */

	return fShear;
	}

void pkdGetColliderInfo(PKD pkd,int iOrder,COLLIDER *c)
{
	/*
	 ** Returns collider info for particle with matching iOrder.
	 */

	PARTICLE *p;
	int i,j,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {
			assert(TYPEQueryACTIVE(p));
			c->id.iPid = pkd->idSelf;
			c->id.iOrder = iOrder;
			c->id.iIndex = i;
			c->id.iOrgIdx = p->iOrgIdx;
			c->fMass = p->fMass;
			c->fRadius = RADIUS(p);
			for (j=0;j<3;j++) {
				c->r[j] = p->r[j];
				c->v[j] = p->v[j];
				c->w[j] = p->w[j];
				c->a[j] = p->a[j];
				}
			c->iColor = p->iColor;
			c->dt = p->dt;
			c->iRung = p->iRung;
			c->bTinyStep = p->bTinyStep;
#ifdef WALLS
			vectorCopy(p->omega_v,c->omega_v); /* unnecessary step (gets overwritten if needed) */
#endif
#ifdef SLIDING_PATCH
			c->dPy = p->dPy;
#endif
#ifdef ORIGIN_HISTOGRAM
			for(j = 0; j < NUM_ORIGIN_BINS; ++j)
				c->origin_bins[j] = p->origin_bins[j];
#endif
			return;
			}
		}
	}

void pkdPutColliderInfo(PKD pkd,const COLLIDER *c,int iOrder2,PARTICLE *p,double dt)
{
	/*
	 ** Stores collider info in particle structure (except id &
	 ** color).  Also, dt is stored in dtPrevCol for inelastic
	 ** collapse checks.
	 **
	 ** dt is set to -1 in a binary merge.
	 **
	 ** NOTE: Because colliding particles have their positions traced
	 ** back to the start of the step using their NEW velocities, it
	 ** is possible for a neighbor to no longer lie inside a previous
	 ** collider's search ball.  To fix this, the ball radius is
	 ** expanded to compensate.  This is needed so we can use resmooth
	 ** instead of smooth to save time after a collision occurs during
	 ** the search interval.
	 **
	 ** ALSO: for periodic BCs, would need to subtract period length
	 ** if there was a wrap!  But, this would require passing all the
	 ** BC info.  Instead, limit fBall2 growth to 50%.
	 **
	 ** IMPORTANT: The contents of "c" are NOT checked here -- it is
	 ** up to the caller to ensure that anything in c that is copied
	 ** to p actually makes sense!
	 */

	int i;
	double dx,dy,dz,r,fBall;
#ifdef COLLMOD_ZML
	if (c->bFrag == 1) { /*DEBUG this really needs a comment...*/
		/*printf("p->iOrgIdx %d p->iColor %d\n", p->iOrgIdx, p->iColor);*/
		p->iColor = 12;
		p->ctimer = 50;
	}
#endif

	p->fMass = c->fMass;
	p->fSoft = SOFT(c);
	dx = c->r[0] - p->r[0];
	dy = c->r[1] - p->r[1];
	dz = c->r[2] - p->r[2];
	for (i=0;i<3;i++) {
		p->r[i] = c->r[i];
		p->v[i] = c->v[i];
		p->w[i] = c->w[i];
		p->a[i] = c->a[i];
#ifdef NEED_VPRED
		p->vPred[i] = p->v[i] - dt*p->a[i];
#endif
		}
	r = sqrt(dx*dx + dy*dy + dz*dz);
	assert(p->fBall2 >= 0.0);
	fBall = sqrt(p->fBall2);
	if (r > 0.5*fBall)
		r = 0.5*fBall; /* only allow a 50% increase */
	p->fBall2 = (fBall + r)*(fBall + r);
	p->dtPrevCol = dt; /* N.B. this is set to -1 in a binary merge */
	p->iPrevCol = iOrder2; /* stored to avoid false collisions */
	p->iRung = c->iRung;
#ifdef WALLS
	vectorCopy(c->omega_v,p->omega_v); /* needed only for particles stuck on cylinders */
#endif
#ifdef SLIDING_PATCH
	p->dPy = c->dPy;
#endif
#ifdef ORIGIN_HISTOGRAM
	for(i = 0; i < NUM_ORIGIN_BINS; ++i)
		p->origin_bins[i] = c->origin_bins[i];
#endif
}

void pkdMerge(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
			  const COLLISION_PARAMS *cp,COLLIDER **cOut,int *pnOut,
			  int iStartStep,double dDelta,double dTime)
{
#ifdef SLIDING_PATCH
	assert(0); /* dPy handling needed (specifically, c->dPy needs to be assigned the correct value) */
#endif

	/*
	** Merges colliders into new centre-of-mass particle (cOut[0])
	** with radius determined by the supplied fixed bulk density.  If
	** cp->dDensity is zero, the radius of the merged particle is
	** instead set by the requirement of conserving total volume.
	** NEW: if cp->iDensityAltCol is also non-zero, particles of that
	** color will be assigned cp->dDensityAltVal, or a value weighted
	** by mass if one particle has color PLANETESIMAL.
	*/

	const double dDenFac = 4.0/3*M_PI;

	COLLIDER *c;
	FLOAT com_pos[3],com_vel[3],rc1[3],rc2[3],vc1[3],vc2[3],ang_mom[3];
	FLOAT m1,m2,m,r1,r2,r,i1,i2,i;
	int k;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	/* handle the bulk density options */
	if (cp->dDensity > 0.0) {
		double dDensity = cp->dDensity;
		if (cp->iDensityAltCol > 0) {
			if (c1->iColor == cp->iDensityAltCol && c2->iColor == cp->iDensityAltCol)
				dDensity = cp->dDensityAltVal;
			else if (c1->iColor == PLANETESIMAL && c2->iColor == cp->iDensityAltCol)
				dDensity = (m1*cp->dDensity + m2*cp->dDensityAltVal)/m;
			else if (c1->iColor == cp->iDensityAltCol && c2->iColor == PLANETESIMAL)
				dDensity = (m1*cp->dDensityAltVal + m2*cp->dDensity)/m;
			}
		r = pow(m/(dDenFac*dDensity),1.0/3);
		}
	else
		r = pow(r1*r1*r1 + r2*r2*r2,1.0/3); /* conserves volume */
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	i = 0.4*m*r*r;

	for (k=0;k<3;k++) {
		com_pos[k] = (m1*c1->r[k] + m2*c2->r[k])/m;
		rc1[k] = c1->r[k] - com_pos[k];
		rc2[k] = c2->r[k] - com_pos[k];
		com_vel[k] = (m1*c1->v[k] + m2*c2->v[k])/m;
		vc1[k] = c1->v[k] - com_vel[k];
		vc2[k] = c2->v[k] - com_vel[k];
	}

	ang_mom[0] = m1*(rc1[1]*vc1[2] - rc1[2]*vc1[1]) + i1*c1->w[0] +
		         m2*(rc2[1]*vc2[2] - rc2[2]*vc2[1]) + i2*c2->w[0];
	ang_mom[1] = m1*(rc1[2]*vc1[0] - rc1[0]*vc1[2]) + i1*c1->w[1] +
				 m2*(rc2[2]*vc2[0] - rc2[0]*vc2[2]) + i2*c2->w[1];
	ang_mom[2] = m1*(rc1[0]*vc1[1] - rc1[1]*vc1[0]) + i1*c1->w[2] +
				 m2*(rc2[0]*vc2[1] - rc2[1]*vc2[0]) + i2*c2->w[2];

	*pnOut = 1;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	c = &(*cOut)[0];

	/* NOTE: id info set in pkdDoCollision(), used only for log */
#ifdef WALLS
	assert(0); /*DEBUG! actually need to set at least the color here, otherwise the macro PARTICLE_STUCK_TO_ROTATING_WALL may seg fault!  of course, we probably don't want the merge option with walls anyway...*/
#endif

	c->fMass = m;
	c->fRadius = r;

	for (k=0;k<3;k++) {
		c->r[k] = com_pos[k];
		c->v[k] = com_vel[k];
		c->w[k] = ang_mom[k]/i;
		c->a[k] = (m1*c1->a[k] + m2*c2->a[k])/m;
	}

#ifdef ORIGIN_HISTOGRAM
	for (k = 0; k < NUM_ORIGIN_BINS; ++k)
		c->origin_bins[k] = c1->origin_bins[k];
        MergeHistograms(m1, c->origin_bins, m2, c2->origin_bins);
#endif /* ORIGIN_HISTOGRAM */
	/* Set merger's timestep to iRung of largest mass. */
#ifdef RORY_EXTENSION
	/* We have to do some messy logic to correct for the possibility
	   of the two colliders being on different rungs, and hence having
	   received their last kick at different times. */
	SetMergerRung(c1,c2,c,dDelta,dTime,iStartStep);
#else
	c->iRung = (c2->fMass > c1->fMass ? c2->iRung : c1->iRung);
#endif
}

int pkdBounce(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
			  double dEpsN,double dEpsT,int iOverlapOption,
			  COLLIDER **cOut,int *pnOut)
{
	/* bounces colliders, preserving particle order */

	COLLIDER *co1,*co2;
	FLOAT n[3],s1[3],s2[3],v[3],s[3],u[3],un[3],ut[3],p[3],q[3];
	FLOAT m1,m2,m,r1,r2,i1,i2,mu,alpha,beta;
	FLOAT a,b,c,d;
	int i;

	*pnOut = 2;
	*cOut = (COLLIDER *) malloc(*pnOut*sizeof(COLLIDER));
	assert(*cOut != NULL);

	(*cOut)[0] = *c1;
	(*cOut)[1] = *c2;

	m1 = c1->fMass;
	m2 = c2->fMass;
	m = m1 + m2;
	r1 = c1->fRadius;
	r2 = c2->fRadius;
	i1 = 0.4*m1*r1*r1;
	i2 = 0.4*m2*r2*r2;
	mu = m1*m2/m;
	alpha = 2.5*(1.0/m1 + 1.0/m2);
	beta = 1.0/(1.0 + alpha*mu); /* always 2/7 for uniform spheres */

	a = 0;
	for (i=0;i<3;i++) {
		n[i] = (c2->r[i] - c1->r[i]);
		a += n[i]*n[i];
		}
	a = 1.0/sqrt(a);
	for (i=0;i<3;i++)
		n[i] *= a;

	s1[0] = r1*(c1->w[1]*n[2] - c1->w[2]*n[1]);
	s1[1] = r1*(c1->w[2]*n[0] - c1->w[0]*n[2]);
	s1[2] = r1*(c1->w[0]*n[1] - c1->w[1]*n[0]);

	s2[0] = r2*(c2->w[2]*n[1] - c2->w[1]*n[2]);
	s2[1] = r2*(c2->w[0]*n[2] - c2->w[2]*n[0]);
	s2[2] = r2*(c2->w[1]*n[0] - c2->w[0]*n[1]);

	for (i=0;i<3;i++) {
		v[i] = c2->v[i] - c1->v[i];
		s[i] = s2[i] - s1[i];
		}

	for (i=0;i<3;i++)
		u[i] = v[i] + s[i];

	a = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];

	if (a >= 0) {
#if (INTERNAL_WARNINGS != 0)
		static int nWarn = 1;
		if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
			(void) fprintf(stderr,"NEAR MISS WARNING #%i (pid=%i): %i & %i (a = %g)\n",
						   nWarn,pkd->idSelf,c1->id.iOrder,c2->id.iOrder,a);
		++nWarn;
#endif /* INTERNAL WARNINGS */
		if (iOverlapOption == OverlapIsError)
			assert(0); /* near miss not allowed */
		else
			return NEAR_MISS; /* particles remain unchanged */
		}
	for (i=0;i<3;i++) {
		un[i] = a*n[i];
		ut[i] = u[i] - un[i];
		}

#ifdef WALLS
	/*
	** For the special case of a particle hitting another particle
	** stuck on a wall, the collision equations reduce to simpler form
	** (essentially the stuck particle is treated as having infinite
	** mass).  NOTE: it is assumed that the velocity and spin of the
	** stuck particle are sensible!
	*/
	if (COLLIDER_STUCK(c1) || COLLIDER_STUCK(c2)) {
		assert(!(COLLIDER_STUCK(c1) && COLLIDER_STUCK(c2))); /* they can't BOTH be stuck! */
		beta = 2.0/7.0; /* redundant */
		if (COLLIDER_STUCK(c1))
			mu = m2;
		else
			mu = m1;
		}
#endif /* WALLS */

	a = (1.0 + dEpsN);
	b = beta*(1.0 - dEpsT);
	for (i=0;i<3;i++)
		p[i] = a*un[i] + b*ut[i];

	a = mu*b;
	q[0] = a*(n[1]*u[2] - n[2]*u[1]);
	q[1] = a*(n[2]*u[0] - n[0]*u[2]);
	q[2] = a*(n[0]*u[1] - n[1]*u[0]);

	a =   m2/m;
	b = - m1/m;
	c = r1/i1;
	d = r2/i2;

#ifdef WALLS
	if (COLLIDER_STUCK(c1)) {
		a = 0.0;
		b = -1.0;
		c = 0.0;
		}
	else if (COLLIDER_STUCK(c2)) {
		a = 1.0;
		b = 0.0;
		d = 0.0;
		}
#endif /* WALLS */

	co1 = &(*cOut)[0];
	co2 = &(*cOut)[1];

	for (i=0;i<3;i++) {
		co1->v[i] += a*p[i];
		co2->v[i] += b*p[i];
		co1->w[i] += c*q[i];
		co2->w[i] += d*q[i];
		}

#ifndef WALLS
	/* dT check...
	{
		double dT;
		dT = -mu*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2]) +
		  0.5*mu*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]) +
		  (r1*c1->w[0] + r2*c2->w[0])*q[0] +
		  (r1*c1->w[1] + r2*c2->w[1])*q[1] +
		  (r1*c1->w[2] + r2*c2->w[2])*q[2] +
		  0.5*alpha*(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
		(void) printf("COLLISION %i & %i e_n=%f e_t=%f dT %e\n",
					  c1->id.iOrder,c2->id.iOrder,dEpsN,dEpsT,dT);
		}
	*/
#endif /* !WALLS */

	return BOUNCE_OK;
	}

void pkdFrag(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,
			 COLLIDER **cOut,int *pnOut)
{
	/*
	 ** Fragments colliders into at most MAX_NUM_FRAG pieces.
	 ** Not implemented yet.  Development notes:
	 ** - be sure new particles are ACTIVE
	 ** - may need to assert(*pnOut <= MAX_NUM_FRAG)
	 ** - remember to set id info for logging purposes
	 ** - and decide on starting iRung values
	 ** - and assign colors
	 */

	assert(0);
	*cOut = NULL;
	*pnOut = 0;
#ifdef ORIGIN_HISTOGRAM
	/* if implemented, copy origin_bins to each new particle */
#endif
	}

void pkdDoCollision(
					PKD pkd,double dt,const COLLIDER *pc1,const COLLIDER *pc2,
					int bPeriodic,int iStartStep,double dDelta,double dTime,
					const COLLISION_PARAMS *CP,int *piOutcome,double *dT,
					COLLIDER *cOut,int *pnOut
#ifdef COLLMOD_ZML
					,const DUST_BINS_PARAMS *DB,int *iDustBin,DustBins *DustBin
#endif
					)
{
	/*
	 ** Processes collision by advancing particle coordinates to
	 ** impact time, determining the collision outcome, and moving the
	 ** resulting particles back to their (new) start-of-step
	 ** positions.  Diagnostic info is returned if the first collider
	 ** resides on the processor (the second collider may not be
	 ** local, in which case this routine is called twice in parallel,
	 ** once for each particle).  Local copies of the collider info
	 ** are used because the particle coordinates need to be advanced
	 ** to the moment of collision.  With periodic boundary
	 ** conditions, the second collider may be offset by one period as
	 ** well.  The offset is removed before calling
	 ** pkdPutColliderInfo().
	 */

	COLLIDER c1,c2,*c = NULL;
	FLOAT fOffset[3],fShear=0.0/*to keep compiler happy*/;
	double v2,ve2,dMergeLimit2,dFragLimit2;
	int bDiagInfo,iOverlap,iOutcome,i,j,n;
#ifdef SPINS_IN_SPACE_FRAME /*DEBUG there may be a better approach*/
#ifdef SLIDING_PATCH
	FLOAT ct,st,wx,wy;
#endif
#endif
#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
	double dMassInDust;
#endif

    /* prevent possible error */
	n = -1;

#ifdef COLLMOD_ZML
	*iDustBin = -1;
	DustBin->dMass = 0.0;
#endif

	/* verbose collision output...
	(void) printf("COLLISION %i [proc %i] & %i [proc %i] (dTime=%g dt=%g)\n",
				  pc1->id.iOrder,pc1->id.iPid,pc2->id.iOrder,pc2->id.iPid,dTime,dt);
	*/

	/* get local copies of collider info */

	c1 = *pc1; /* struct copy */
	c2 = *pc2;

	bDiagInfo = (c1.id.iPid == pkd->idSelf);

	if (bDiagInfo && dT != NULL) *dT = 0.0;

#ifdef WALLS
	if (c2.id.iOrder < 0) { /* wall collision */
		assert(bDiagInfo); /* wall is not a particle, so only 1 CPU involved */
		assert(!COLLIDER_STUCK(&c1));
		pkdWallsDoCollision(pkd,CP,&c1,&c2,dt,piOutcome);
		if (*piOutcome == DEATH) {
			pkdDeleteParticle(pkd,&pkd->pStore[c1.id.iIndex]);
			cOut[0] = c2;
			*pnOut = 1;
			}
		else {
			pkdPutColliderInfo(pkd,&c1,c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
			cOut[0] = c1;
			cOut[1] = c2;
			*pnOut = 2;
			}
		return;
		}
#endif /* WALLS */

	if (bPeriodic)
		fShear = pkdApplyBCs(pkd,dTime,&c1,&c2,fOffset); /* note: NOT dTime + dt, because c1 & c2 (and the collision circumstances) are defined at/from dTime (start of step), not dTime + dt */

	/* handle overlap cases if applicable */

	iOverlap = OverlapIgnore;

	if (dt <= 0.0) {
		iOverlap = CP->iOverlapOption;
		switch (iOverlap) {
		case OverlapIgnore:
		case OverlapBackstep:
			break; /* do nothing */
		case OverlapIsError:
			assert(0);
		case OverlapAdjPos:
			dt = 0.0; /* effect is instantaneous */
			break;
		case OverlapRepel:
			assert(0); /* shouldn't be here */
		case OverlapMerge:
			dt = 0.0;
			break;
		default:
			assert(0);
			}
		}

	/* advance coordinates to impact time */

	for (i=0;i<3;i++) {
		c1.r[i] += c1.v[i]*dt;
		c2.r[i] += c2.v[i]*dt;
	}

#ifdef WALLS
	if (COLLIDER_STUCK_ON_ROTATING_WALL(&c1,&CP->WP)) {
		vectorCopy(pc1->r,c1.r); /* override previous update */
		wallsRotate(&CP->WP,COLLIDER_WALL_ID(&c1),c1.r,c1.v,c1.omega_v,dt);
		}
	else if (COLLIDER_STUCK_ON_ROTATING_WALL(&c2,&CP->WP)) {
		vectorCopy(pc2->r,c2.r); /* override previous update */
		wallsRotate(&CP->WP,COLLIDER_WALL_ID(&c2),c2.r,c2.v,c2.omega_v,dt);
		}
#endif

#ifdef SPINS_IN_SPACE_FRAME /*DEBUG 2/13/08: now taking spins to be in rotating frame, so this conversion no longer needed*/
	                    /* RP-DEBUG 6/2/09: Really?  Spins are in the space frame, think I (though we aren't yet consistent) */
#ifdef SLIDING_PATCH
	/* for sliding patch, need to transform spins to patch frame as well */

	ct = cos(pkd->PP->dOrbFreq*(dTime + dt)); /* could also use pkd->dTime */
	st = sin(pkd->PP->dOrbFreq*(dTime + dt));
	wx =  ct*c1.w[0] + st*c1.w[1];
	wy = -st*c1.w[0] + ct*c1.w[1];
	c1.w[0] = wx;
	c1.w[1] = wy;
	wx =  ct*c2.w[0] + st*c2.w[1];
	wy = -st*c2.w[0] + ct*c2.w[1];
	c2.w[0] = wx;
	c2.w[1] = wy;
#endif
#endif

	/* determine collision outcome */

	switch (iOverlap) {	/* handle overlap condition, if needed */
	case OverlapIgnore:
	case OverlapBackstep:
		break;
	case OverlapAdjPos:
	  {
		  /* move particles apart along line of centers until just touching */
		  double r_hat[3],r_mag,dr,m_tot;
		  int k;
		  r_mag = 0.0;
		  for (k=0;k<3;k++) {
			  r_hat[k] = c1.r[k] - c2.r[k]; /* points from c2 to c1 */
			  r_mag += r_hat[k]*r_hat[k];
			  }
		  r_mag = sqrt(r_mag); /* separation */
		  assert(r_mag > 0.0); /* can't handle exactly colocated bodies */
		  dr = c1.fRadius + c2.fRadius;
		  assert(r_mag <= dr); /* this is the whole point: they're overlapping */
		  dr -= r_mag; /* expansion factor */
		  assert(dr >= 0.0);
		  m_tot = c1.fMass + c2.fMass;
		  assert(m_tot > 0.0);
		  for (k=0;k<3;k++) {
			  r_hat[k] /= r_mag;
			  c1.r[k] += (c2.fMass/m_tot)*dr*r_hat[k];
			  c2.r[k] -= (c1.fMass/m_tot)*dr*r_hat[k];
			  }
		  /* handle as normal collision (or miss) from here on */
		  break;
		  }
	case OverlapMerge:
		/* force merge */
		assert(CP->iOutcomes & MERGE);
		iOutcome = MERGE;
		pkdMerge(pkd,&c1,&c2,CP,&c,&n,iStartStep,dDelta,dTime);
		assert(c != NULL);
		assert(n == 1);
#ifdef COLLMOD_ZML
		dMassInDust = 0.0;
#endif
		goto finish;
	default:
		assert(0);
		}

#ifdef RUBBLE_ZML	

	/* should this be merged with stuff below? (i.e. after "CP->iOutcomes & FRAG" line) */
	if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {
		printf("c1.id.iOrder = %i c2.id.iOrder = %i\n",
			   c1.id.iOrder, c2.id.iOrder);
		rubRubbleCollide(&c1,&c2,CP,&c,&dMassInDust,&n); /*dMinMass to be part CP*/

		if (n == 1)
			iOutcome = MERGE; /* could have produced dust */
		else if (n == 2) {
			assert(0); /* this should never happen */
			iOutcome = BOUNCE; /* ditto */
			/* would need to call pkdBounce() here */
			}
		else if (n > 2)
			iOutcome = FRAG; /* colliders turned into rubble piles */
		else
			assert(0); /* not possible */
		
		/*printf("iOutcome = %s\n",iOutcome);*/

		goto finish; /* need to call add or delete particle routine */
		}

#endif /* RUBBLE_ZML */

#ifdef COLLMOD_ZML
	printf("c1.id.iOrder, iPid, iIndex, iOrgIdx = %i %i %i %i c2.id.iOrder, iPid, iIndex, iOrgIdx = %i %i %i %i\n",
	       c1.id.iOrder, c1.id.iPid, c1.id.iIndex, c1.id.iOrgIdx, c2.id.iOrder, c2.id.iPid, c2.id.iIndex, c2.id.iOrgIdx);
	collCollModCollide(&c1,&c2,CP,&c,&dMassInDust,&n,pkd); /*dMinMass to be part CP*/

	/* in this collision model all collision outcomes could produce dust */

#ifdef ORIGIN_HISTOGRAM
	/*
	** Make dust origin histogram now.  For now this is just the
	** mass-weighted origin histograms of both colliders.  This should
	** work for most cases, but we might want to do something
	** different for a disruptive hit-and-run.  NOTE: we do this even
	** if no mass goes to dust, and even if we're not on collider 1's
	** processor, because in the disruptive case of more than 2
	** remnants, we want to set the origin histograms of *all*
	** fragments to be these mass-weighted values.
	*/
	for (i=0;i<NUM_ORIGIN_BINS;i++)
		DustBin->origin_bins[i] = c1.origin_bins[i];
	MergeHistograms(c1.fMass,DustBin->origin_bins,c2.fMass,c2.origin_bins);
#endif /* ORIGIN_HISTOGRAM */

	if (dMassInDust > 0.0) {
		collLocateBin(&pkd->pStore[c1.id.iIndex],DB,iDustBin,NULL/*ignore r*/);
		DustBin->dMass = dMassInDust;
		}

	if (n == 0)
		iOutcome = EXPLODE; /* collision resulted in only unresolved dust (use better label?) */
	else if (n == 1)
		iOutcome = MERGE; /* this could still be an erosive collision but the remnants are below resolution */
	else if ((n == 2) && (dMassInDust == 0)) {
		double dEpsN = 0.5;
		double dEpsT = 1.0;
		iOutcome = BOUNCE; /* hit-and-run */
		/* call pkdBounce here - this is a true hit and run */
		assert(c == NULL);
		pkdBounce(pkd,&c1,&c2,dEpsN,dEpsT,CP->iOverlapOption,&c,&n);
		}
	else if (n == 2)
		iOutcome = BOUNCE; /* not a true hit and run - at least some damage */
	else if (n > 2)
		iOutcome = FRAG; /* disruptive collision with multiple remnants above resolution limit */
	else
		assert(0); /* no other possibilities! */

	goto finish; /* need to call add or delete particle routine */
		
#endif /* COLLMOD_ZML */
	
	v2 = 0.0;
  
	for (i=0;i<3;i++)
		v2 += (c2.v[i] - c1.v[i])*(c2.v[i] - c1.v[i]);

	ve2 = 2.0*(c1.fMass + c2.fMass)/(c1.fRadius + c2.fRadius);

	dMergeLimit2 = CP->dMergeLimit*CP->dMergeLimit;
	if (CP->dMergeLimit > 0.0)
		dMergeLimit2 *= ve2; /* apply scaling, if any */

	dFragLimit2 = CP->dFragLimit*CP->dFragLimit;
	if (CP->dFragLimit > 0.0)
		dFragLimit2 *= ve2;

	/* note: dImpactEnergy = 0.5*c1.fMass*c2.fMass/(c1.fMass + c2.fMass)*v2; */

	iOutcome = MISS;

#ifdef RUBBLE_ZML  
	if (CP->iRubForcedOutcome == RUB_FORCED_MERGE)
#else
	if (CP->iOutcomes == MERGE ||
		((CP->iOutcomes & MERGE) && v2 <= dMergeLimit2))
#endif
	  {
		iOutcome = MERGE;
		pkdMerge(pkd,&c1,&c2,CP,&c,&n,iStartStep,dDelta,dTime);
		assert(c != NULL);
		assert(n == 1);
#ifndef RUBBLE_ZML /* if NOT rubble!... */
		if (CP->iOutcomes & (BOUNCE | FRAG)) { /* check if spinning too fast */
			double w2max,w2=0.0;
			w2max = c->fMass/(c->fRadius*c->fRadius*c->fRadius);
			for (i=0;i<3;i++)
				w2 += c->w[i]*c->w[i];
			if (w2 > w2max) {
				int rv;
				free((void *) c);
				iOutcome = BOUNCE;
				rv = pkdBounce(pkd,&c1,&c2,CP->dEpsN,CP->dEpsT,CP->iOverlapOption,&c,&n);
				assert(rv == BOUNCE_OK);
				assert(c != NULL);
				assert(n == 2);
				}
			}
#endif
		}
#ifdef RUBBLE_ZML
	else if (CP->iRubForcedOutcome == RUB_FORCED_BOUNCE)
#else
    else if (CP->iOutcomes & BOUNCE)
#endif
	  {
		/*DEBUG be sure any changes here also reflected in pkdAggsDoCollision()*/
		double dEpsN,dEpsT;
		iOutcome = BOUNCE;
#ifdef SPRINGS
	bounce:
#endif
		dEpsN = dEpsT = 1.0;
		if (c1.bTinyStep || c2.bTinyStep) {
			dEpsN = CP->dCollapseEpsN;
			dEpsT = CP->dCollapseEpsT;
			}
		else if ((CP->iSlideOption == EscVel && v2 < CP->dSlideLimit2*ve2) ||
				 (CP->iSlideOption == MaxTrv && v2 < CP->dSlideLimit2)) {
			dEpsN = CP->dSlideEpsN;
			dEpsT = CP->dSlideEpsT;
			}
		else if (v2 > CP->dCrushLimit) { /* add frag check */
			dEpsN = CP->dCrushEpsN;
			dEpsT = CP->dCrushEpsT;
			}
		else if (CP->iEpsNOption == ConstEps) {
			dEpsN = CP->dEpsN;
			dEpsT = CP->dEpsT;
			}
		else {
			const double ls = 1.49597892e13;
			const double ts = 5.0226355648e6;
			double rn[3],vn=0,dTmp;
			for (i=0;i<3;i++) {
				rn[i] = c2.r[i] - c1.r[i];
				vn += (c2.v[i] - c1.v[i])*rn[i];
				}
			/* note vn >= 0 ==> near miss -- pkdBounce() will deal with it */
			if (vn < 0) {
				vn = - vn/sqrt(rn[0]*rn[0] + rn[1]*rn[1] + rn[2]*rn[2]);
				vn *= ls/ts; /* conversion to cm/s */
				switch (CP->iEpsNOption) {
				case PowerLaw:
					dEpsN = CP->dEpsNCoef*pow(vn,CP->dEpsNExp);
					break;
				case Compacted: /* Hatzes et al. 1988 */
					dTmp = c1.fRadius*c2.fRadius/(c1.fRadius + c2.fRadius)*ls;
					dEpsN = -1.90*exp(-(-0.01*dTmp + 0.41)*vn);
					break;
				case Borderies: /* Borderies et al. 1984 */
					if (vn >= CP->dEpsNVStar) {
						dTmp = CP->dEpsNVStar/vn;
						dTmp *= dTmp;
						dEpsN = sqrt(-2.0/3.0*dTmp + sqrt(10.0/3.0*dTmp - 5.0/9.0*dTmp*dTmp));
						}
					else
						dEpsN = 1.0;
					break;
				default:
					assert(0);
					}
				if (dEpsN < CP->dEpsNMin)
					dEpsN = CP->dEpsNMin;
				else if (dEpsN > 1.0)
					dEpsN = 1.0;
				dEpsT = CP->dEpsT;
				}
			}
		if (pkdBounce(pkd,&c1,&c2,dEpsN,dEpsT,CP->iOverlapOption,&c,&n) == NEAR_MISS)
			iOutcome = MISS;
		assert(c != NULL);
		assert(n == 2);
		}
	else if ((CP->iOutcomes & FRAG) && v2 > dFragLimit2) {
		iOutcome = FRAG;
#ifdef SPRINGS
		/* treat outcome as bounce, but remove all springs attached to affected particles (handled in master) */
		goto bounce;
#endif
		pkdFrag(pkd,&c1,&c2,&c,&n);
		assert(c != NULL);
		assert(n <= MAX_NUM_FRAG);
		}

	finish:

#ifdef COLLMOD_ZML
	assert(n >= 0);
#else
	assert(n > 0);
#endif

	if (CP->iOverlapOption == OverlapIsError)
		assert(iOutcome != MISS); /* SOMETHING has to happen... */

	if (bDiagInfo) {
		*piOutcome = iOutcome;
		if (dT != NULL) {
			double moi;
			*dT = 0.0; /* redundant */
			for (j=0;j<n;j++) {
				moi = 0.4*c[j].fMass*c[j].fRadius*c[j].fRadius;
				for (i=0;i<3;i++)
					*dT += c[j].fMass*(c[j].v[i]*c[j].v[i]) +
						moi*(c[j].w[i]*c[j].w[i]);
				}
			moi = 0.4*c1.fMass*c1.fRadius*c1.fRadius;
			for (i=0;i<3;i++)
				*dT -= c1.fMass*(c1.v[i]*c1.v[i]) + moi*(c1.w[i]*c1.w[i]);
			moi = 0.4*c2.fMass*c2.fRadius*c2.fRadius;
			for (i=0;i<3;i++)
				*dT -= c2.fMass*(c2.v[i]*c2.v[i]) + moi*(c2.w[i]*c2.w[i]);
			*dT *= 0.5;
			}
		/*printf("Compare: dT = %e\n",*dT);*/
		for (i=0;i<(n < MAX_NUM_FRAG ? n : MAX_NUM_FRAG);i++) /* only first MAX_NUM_FRAG reported */
			cOut[i] = c[i];
		*pnOut = n;
		}

	/* trace particles back to start of step */

	for (i=0;i<n;i++) /* only actually affects local particles */
		for (j=0;j<3;j++)
			c[i].r[j] -= c[i].v[j]*dt;

#ifdef WALLS
	for (i=0;i<n;i++)
		if (COLLIDER_STUCK_ON_ROTATING_WALL(&c[i],&CP->WP)) {
			for (j=0;j<3;j++)
				c[i].r[j] += c[i].v[j]*dt; /* undo traceback */
			wallsRotate(&CP->WP,COLLIDER_WALL_ID(&c[i]),c[i].r,c[i].v,c[i].omega_v,-dt);
			}
#endif

#ifdef SPINS_IN_SPACE_FRAME
#ifdef SLIDING_PATCH
	/* and transform spin(s) back to space frame */
	
	for (i=0;i<n;i++) {
		wx = ct*c[i].w[0] - st*c[i].w[1];
		wy = st*c[i].w[0] + ct*c[i].w[1];
		c[i].w[0] = wx;
		c[i].w[1] = wy;
		}
#endif
#endif

	/* handle output cases */

	if (n == 1) { /* merge */
#ifdef SLIDING_PATCH
		/* set y-momentum to be center-of-mass value */
		c[0].dPy = (c1.fMass*c1.dPy + c2.fMass*c2.dPy)/(c1.fMass + c2.fMass);
#endif
#ifdef RUBBLE_ZML
		pkdMergerReorder(pkd,&c1,&c2,c,cOut,dt,bDiagInfo,CP,dMassInDust);
#else
		pkdMergerReorder(pkd,&c1,&c2,c,cOut,dt,bDiagInfo);
#endif
		}
	else if (n == 2) { /* bounce (or miss) */
#ifdef SLIDING_PATCH
		/* adjust y-momenta according to velocity and position (backdrift) change */
		c[0].dPy += (c[0].v[1] - c1.v[1]) + 2.0*pkd->PP->dOrbFreq*(c[0].r[0] - (c1.r[0] - c1.v[0]*dt));
		c[1].dPy += (c[1].v[1] - c2.v[1]) + 2.0*pkd->PP->dOrbFreq*(c[1].r[0] - (c2.r[0] - c2.v[0]*dt));
#endif
		if (c1.id.iPid == pkd->idSelf)
			pkdPutColliderInfo(pkd,&c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
		if (c2.id.iPid == pkd->idSelf) {
			if (bPeriodic) {
				for (i=0;i<3;i++)
					c[1].r[i] -= fOffset[i];
				c[1].v[1] -= fShear; /* if applicable */
#ifdef SLIDING_PATCH
				c[1].dPy += fShear/3.0;
#endif
				}
			pkdPutColliderInfo(pkd,&c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
			}
		}
	else { /* fragmentation */
#ifdef SLIDING_PATCH
		assert(0); /* need to worry about dPy changes */
#endif
#ifdef RUBBLE_ZML
		if (CP->iRubForcedOutcome == RUB_FORCED_NONE) {
			PARTICLE p;
			if (c1.id.iPid == pkd->idSelf) {
				pkdPutColliderInfo(pkd,&c[0],c2.id.iOrder,&pkd->pStore[c1.id.iIndex],dt);
				pkd->pStore[c1.id.iIndex].iColor = CP->iRubColor; /* override color */
				}
			if (c2.id.iPid == pkd->idSelf) {
				pkdPutColliderInfo(pkd,&c[1],c1.id.iOrder,&pkd->pStore[c2.id.iIndex],dt);
				pkd->pStore[c2.id.iIndex].iColor = CP->iRubColor;
				}
			/*
			 ** For now, create new particles on 1st processor.
			 ** Hopefully domain decomposition will rebalance quickly.
			 ** Note that rubble pile particles created by 2nd
			 ** processor are ignored.  Why do we bother making them
			 ** at all in that case?  Good question.
			 */
			if (c1.id.iPid == pkd->idSelf) {
				double a[3];
				double temp_mass;

				p = pkd->pStore[c1.id.iIndex]; /* quick and dirty initialization */
				p.iColor = CP->iRubColor;
				assert(p.iColor != PLANETESIMAL);
				p.dtPrevCol = 0;
				p.iPrevCol = INT_MAX;

				temp_mass = 1/(c[1].fMass+c[2].fMass);
				a[0] = temp_mass*(c[1].fMass*c[1].a[0] + c[2].fMass*c[2].a[0]);
				a[1] = temp_mass*(c[1].fMass*c[1].a[1] + c[2].fMass*c[2].a[1]);
				a[2] = temp_mass*(c[1].fMass*c[1].a[2] + c[2].fMass*c[2].a[2]);

				for (i=2;i<n;i++) {
					/* rubble pieces do not have acceleration in RUBBLE_PILE struct */
					c[i].a[0] = a[0];
					c[i].a[1] = a[1];
					c[i].a[2] = a[2];
					pkdPutColliderInfo(pkd,&c[i],-1,&p,dt);
					/* fix up some things pkdPutColliderInfo() doesn't take care of */
					p.iOrgIdx = c[i].id.iOrgIdx;
					pkdNewParticle(pkd,p);
					}
				}
			}
#elif defined(COLLMOD_ZML)
		{

			PARTICLE *p;

			if (n == 0) {

				if (c1.id.iPid == pkd->idSelf) {
					printf("No resolved remnants -- all mass going to dust (%e)\n", dMassInDust);
					p = &pkd->pStore[c1.id.iIndex];
					pkdDeleteParticle(pkd,p);
					}

				if (c2.id.iPid == pkd->idSelf) {
					p = &pkd->pStore[c2.id.iIndex];
					pkdDeleteParticle(pkd,p);
					}

				}

			else if (n > 2) {

				/*
				** When more than 2 remnants are resolved, need to add
				** particles and put mass in dust.
				*/

				if (c1.id.iPid == pkd->idSelf) {
					/*puts("TAG 1");*/
					p = &pkd->pStore[c1.id.iIndex];
					pkdPutColliderInfo(pkd,&c[0],c2.id.iOrder,p,dt);
#ifdef ORIGIN_HISTOGRAM
					/* all new fragments inherit mass-weighted origin bins of colliders */
					for (i=0;i<NUM_ORIGIN_BINS;i++)
						p->origin_bins[i] = DustBin->origin_bins[i];
#endif /* ORIGIN_HISTOGRAM */
					}
				if (c2.id.iPid == pkd->idSelf) {
					/*puts("TAG 2");*/
					p = &pkd->pStore[c2.id.iIndex];
					pkdPutColliderInfo(pkd,&c[1],c1.id.iOrder,p,dt);
#ifdef ORIGIN_HISTOGRAM
					for (i=0;i<NUM_ORIGIN_BINS;i++)
						p->origin_bins[i] = DustBin->origin_bins[i];
#endif /* ORIGIN_HISTOGRAM */
					}

				/*
				** For now, create new particles on 1st processor.
				** Hopefully domain decomposition will rebalance
				** quickly...
				*/

				if (c1.id.iPid == pkd->idSelf) {
					/*puts("TAG 3");*/
					PARTICLE pNew;

					pNew = pkd->pStore[c1.id.iIndex]; /* quick and dirty initialization (struct copy) */
					/*printf("c1.id.iIndex %d\n", c1.id.iIndex);*/
					pNew.iPrevCol = INT_MAX;
					/*puts("---");*/
					for (i=2;i<n;i++) {
						printf("NEW %d\n", i);
						pkdPutColliderInfo(pkd,&c[i],c1.id.iOrder,&pNew,dt); /*pNew is re-used*/
						pNew.iOrgIdx = pkd->pStore[c1.id.iIndex].iOrgIdx; /*Always take c[i] index//c[i].id.iOrgIdx;*/
						printf("pNew.iOrder %d Mass %E\n",pNew.iOrder, pNew.fMass);
						pkdNewParticle(pkd,pNew);
						/*puts("---");*/
						}
					}
				}
			}

#else /* (COLLMOD_ZML) */
	assert(0); /* not implemented yet */
	/* note in this case new iOrders start at pkd->nMaxOrderDark */
#endif /* !RUBBLE_ZML && !COLLMOD_ZML */

		} /* fragmentation */

#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
	if (bDiagInfo && n > 0) {
		printf("Mass going to dust = %e\n", dMassInDust);
		if (n > 1 && c[0].fMass < c[1].fMass){
			printf("Largest remnant mass = %e\n", c[1].fMass);
			}
		else{
			printf("Largest remnant mass = %e\n", c[0].fMass);
			}
	}
#endif /* RUBBLE_ZML || COLLMOD_ZML */

	/* Free resources */

	free((void *) c); /* (if c is NULL, op is ignored) */
}

void pkdResetColliders(PKD pkd,int iOrder1,int iOrder2)
{
	/*
	 ** Set SMOOTHACTIVE flag for those particles that were involved
	 ** in the collision or that were about to collide with one of the
	 ** particles involved in the collision. Reset otherwise. This
	 ** means new collision circumstances will only be determined for
	 ** these particles, while the remaining (active) particles will
	 ** rely on previously detected potential collisions.
	 */

	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
#ifdef WALLS
		/*
		** Since walls don't react to collisions, fewer resets needed.
		** The logic is slightly more complicated, but worth it for
		** potential large savings.
		*/
		if (p->iOrder == iOrder1 || p->iOrder == iOrder2 ||
			(p->iOrderCol == iOrder1 && iOrder1 >= 0) ||
			(p->iOrderCol == iOrder2 && iOrder2 >= 0))
			TYPESet(p,TYPE_SMOOTHACTIVE);
		else
			TYPEReset(p,TYPE_SMOOTHACTIVE);
#else
		if (p->iOrder == iOrder1 || p->iOrder == iOrder2 ||
			p->iOrderCol == iOrder1 || p->iOrderCol == iOrder2)
			TYPESet(p,TYPE_SMOOTHACTIVE);
		else
			TYPEReset(p,TYPE_SMOOTHACTIVE);
#endif /*!WALLS*/
		}
	}


#ifdef ORIGIN_HISTOGRAM /*This has been copied from rubble.c - still need origin histograms with COLLMOD. */
/*
 ** Implement histograms of "region of material origin" per particle.
 ** 28 April 2006, gwl.  The fixed-length array of floating-point
 ** values, 'origin_bins', gives the fraction of material that
 ** originated in a specific radial bin.  That is, the sum of the
 ** elements of the array is always 1.  Initially, all the material
 ** came from the location of particle.  Subsequent interactions
 ** (collisions, fragmentation, mergers, dust) modify the composition
 ** histogram. Modified slightly from those in rubble.c by ZML 26.04.12.
 */

void MergeHistograms(FLOAT m1,FLOAT bins1[],FLOAT m2,const FLOAT bins2[])
{
	/*
	 ** Merge two histograms together, weighted by masses.  The first
	 ** histogram is modified.  NOTE: if either m1 or m2 (not both!) is
	 ** zero, this function effectively copies the non-zero data.
	 */
	
	int i = 0;
	
	assert(m1 > 0.0 || m2 > 0.0);
	
	for (i=0;i<NUM_ORIGIN_BINS;i++) {
		bins1[i] = (m1*bins1[i] + m2*bins2[i])/(m1 + m2);
	}
}

void pkdInitializeOriginBins(PKD pkd,const DUST_BINS_PARAMS* DB)
{
	/*DEBUG this function should perhaps be made consistent with the
	  dust-bin location functions rubLocateBin() and collLocateBin(),
	  to ensure the same formula is being used...*/
	
	/*
	** Locate "home" bin of particles and set that origin bin to 1.0.
	** This should be called exactly once at the beginning of the
	** simulation.
	*/
	
	PARTICLE *p;
	FLOAT r_xy;
	int i,j,iBin,nLocal = pkdLocal(pkd);
	
	for (i=0;i<nLocal;i++) {
		
		p = &pkd->pStore[i];
		
		for(j=0;j<NUM_ORIGIN_BINS;j++)
			p->origin_bins[j] = 0.0;
		
		/*
		** Use distance in the projected x-y plane to find origin bin the
		** particle belongs to.  Note that the number of origin bins is
		** fixed at compile time, while the number of dust bins is a
		** run-time parameter.  The origin bins use the same inner and
		** outer bounds, however.
		*/
		
		r_xy = sqrt(p->r[0]*p->r[0] + p->r[1]*p->r[1]);
		
		iBin = NUM_ORIGIN_BINS*(r_xy - DB->dDustBinsInner)/(DB->dDustBinsOuter - DB->dDustBinsInner);
		
		if (iBin < 0 || iBin >= NUM_ORIGIN_BINS) {
			fprintf(stderr,"WARNING: Particle outside bin range: iOrder = %i r_xy = %g inner = %g outer = %g\n",
					p->iOrder,r_xy,DB->dDustBinsInner,DB->dDustBinsOuter);
			/* for now, assign origin bin to closest bin */
			if (iBin < 0)
				p->origin_bins[0] = 1.0;
			else
				p->origin_bins[NUM_ORIGIN_BINS - 1] = 1.0;
			}
		else
			p->origin_bins[iBin] = 1.0; /* entire initial particle came from this bin */
		
	}
}

#endif /*ORIGIN_HISTOGRAM*/

#endif /* COLLISIONS */
