#include <math.h>
#include <assert.h>
#include "smoothfcn.h"
#include "matrix.h"
/*DEBUG following functions should be collected somewhere else...*/
int sign(double A);
#define sign(A) ((A)<0.?-1:(A)>0.?1:0)
double min(double A,double B);
#define min(A,B) ((A) > (B) ? (B) : (A))
double max(double A,double B);
#define max(A,B) ((A) > (B) ? (A) : (B))

#ifdef COLLISIONS
#include "ssdefs.h"
#include "collision.h"
#include "polyroots.h"
#endif

#ifdef AGGS
#include "aggs.h"
#endif

#ifdef WALLS
#include "walls.h"
#endif

#ifdef SPRINGS
#include "random.h"
#endif

/*
 Change the way the Balsara Switch is applied:
*/
/*
#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define SWITCHCOMBINE(a,b) (a->BalsaraSwitch > b->BalsaraSwitch ? a->BalsaraSwitch : b->BalsaraSwitch)
#define SWITCHCOMBINE(a,b) (a->BalsaraSwitch*b->BalsaraSwitch)
#define SWITCHCOMBINE(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch > 0.5 || \
           (a->BalsaraSwitch > 0.5 && (dx*a->aPres[0]+dy*a->aPres[1]+dz*a->aPres[2]) > 0) || \
           (b->BalsaraSwitch > 0.5 && (dx*b->aPres[0]+dy*b->aPres[1]+dz*b->aPres[2]) < 0)) ? 1 : 0)
#define SWITCHCOMBINEA(a,b) (a->BalsaraSwitch >= 1 || b->BalsaraSwitch >= 1 ? 1 : 0)
#define SWITCHCOMBINEA(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch)*(a->ShockTracker > b->ShockTracker ? a->ShockTracker : b->ShockTracker))
*/

#define ACCEL(p,j) (((PARTICLE *)(p))->a[j])
#define KPCCM 3.085678e21

#ifdef SHOCKTRACK
/* Shock Tracking on: p->ShockTracker and p->aPres are defined */

#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define SWITCHCOMBINEA(a,b) ((a->BalsaraSwitch*b->BalsaraSwitch)*a->ShockTracker*b->ShockTracker)
#define SWITCHCOMBINEB(a,b) (a->BalsaraSwitch*b->BalsaraSwitch)

#define ACCEL_PRES(p,j) (((PARTICLE *)(p))->aPres[j])
#define ACCEL_COMB_PRES(p,j) ((((PARTICLE *)(p))->a[j])+=(((PARTICLE *)(p))->aPres[j]))

#else

#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define SWITCHCOMBINEA(a,b) SWITCHCOMBINE(a,b)
#define SWITCHCOMBINEB(a,b) SWITCHCOMBINE(a,b)

#define ACCEL_PRES(p,j) (((PARTICLE *)(p))->a[j])
#define ACCEL_COMB_PRES(p,j) 

#endif

#ifdef PEAKEDKERNEL
/* Standard M_4 Kernel */
#define BALL2(a) ((a)->fBall2)
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else ak = 0.25*ak*ak*ak; \
        }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
            if (adk < 2./3.) { \
               if (adk > 0) adk = -1/adk; \
			   } \
            else { \
               adk = -3 + 2.25*adk; \
			   } \
			} \
		else { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		}
#else
#ifdef M43D
/* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
#define BALL2(a) ((a)->fBall2)
#define KERNEL(ak,ar2) { \
		ak = sqrt(ar2); \
		if (ar2 < 1.0) ak = 6.*0.25/350./3. *(1360+ar2*(-2880 \
			 +ar2*(3528+ak*(-1890+ak*(-240+ak*(270-6*ar2)))))); \
		else ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
	                 +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
	                 +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
		else adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
	                 +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
                }

#else
#ifdef HSHRINK
/* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
#define dSHRINKFACTOR        0.805995977
#define BALL2(a) ((a)->fBall2*(dSHRINKFACTOR*dSHRINKFACTOR))
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
		else ak = 0.0; \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else if (ar2 < 4.0) { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
		else adk = 0.0; \
                }

#else
/* Standard M_4 Kernel */
#define BALL2(a) ((a)->fBall2)
#define KERNEL(ak,ar2) { \
		ak = 2.0 - sqrt(ar2); \
		if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
		else ak = 0.25*ak*ak*ak; \
                }
#define DKERNEL(adk,ar2) { \
		adk = sqrt(ar2); \
		if (ar2 < 1.0) { \
			adk = -3 + 2.25*adk; \
			} \
		else { \
			adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
			} \
                }
#endif
#endif
#endif

void NullSmooth(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf) {
}

void initDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	}

void combDensity(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	}

void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	FLOAT ih2,r2,rs,fDensity;
	int i;

	ih2 = 4.0/BALL2(p);
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		fDensity += rs*nnList[i].pPart->fMass;
		}
	p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}

void DensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;

	ih2 = 4.0/(BALL2(p));
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		}
	}

/* Mark functions:
   These functions calculate density and mark the particles in some way
   at the same time.
   Only SmoothActive particles do this gather operation (as always).

   MarkDensity:
	 Combination: if porig DensZeroed combine porig+pcopy into porig
                  if pcopy DensZeroed set porig = pcopy 
     Smooth:
       Active particles get density, nbrofactive
	   All gather neighbours are labelled as DensZeroed, get density
       --> effectively all particles and gather neighbours get density and are labelled DensZeroed
       --> These densities will potentially be lacking scatter neighbours so only correct
           if all particles involved in this operation OR scatter later added
       Gather/Scatter Neighbours of Active Particles get nbrofactive

   MarkSmooth:
     Go through full tree looking for particles than touch a smooth active particle
     and mark them with specified label: eg. TYPE_Scatter

   MarkIIDensity:
     Init:        Densactive particles not dens zeroed get dens zeroed
	 Combination: If pcopy is active make porig active
	              Densactive particles only:
                        if porig DensZeroed combine density porig+pcopy into porig
                        if pcopy DensZeroed set density  porig = pcopy 
     Smooth:
       Densactive: get density
                   Densactive gather neighbours get density, DensZeroed (Nbrofactive if reqd)
	   Not Densactive, but Active:
                   get nbrofactive
			       Densactive gather neighbours get density, DensZeroed	(Nbrofactive)
       Not Densactive, Not Active 
                   get nbrofactive if gather neighbour active
			       Densactive gather neighbours get density, DensZeroed	
*/

void initParticleMarkDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	TYPESet((PARTICLE *) p,TYPE_DensZeroed);
	}

void initMarkDensity(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	}

void combMarkDensity(void *p1,void *p2)
{
	if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) 
		((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
		((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
		}
	((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
	}

void MarkDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	assert(0);
}

void MarkDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i;
	unsigned int qiActive;

	assert(TYPETest(p,TYPE_GAS));
	ih2 = 4.0/(BALL2(p));
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	if (TYPETest(p,TYPE_ACTIVE)) {
		TYPESet(p,TYPE_NbrOfACTIVE);
		for (i=0;i<nSmooth;++i) {
			r2 = nnList[i].fDist2*ih2;
			KERNEL(rs,r2);
			rs *= fNorm;
			q = nnList[i].pPart;
			assert(TYPETest(q,TYPE_GAS));
			p->fDensity += rs*q->fMass;
			if (TYPETest(q,TYPE_DensZeroed)) 
				q->fDensity += rs*p->fMass;
			else {
				q->fDensity = rs*p->fMass;
				TYPESet(q, TYPE_DensZeroed);
				}
			TYPESet(q,TYPE_NbrOfACTIVE);
			}
		} 
	else {
		qiActive = 0;
		for (i=0;i<nSmooth;++i) {
			r2 = nnList[i].fDist2*ih2;
			KERNEL(rs,r2);
			rs *= fNorm;
			q = nnList[i].pPart;
			assert(TYPETest(q,TYPE_GAS));
			if (TYPETest(p,TYPE_DensZeroed)) 
				p->fDensity += rs*q->fMass;
			else {
				p->fDensity = rs*q->fMass;
				TYPESet(p,TYPE_DensZeroed);
				}
			if (TYPETest(q,TYPE_DensZeroed)) 
				q->fDensity += rs*p->fMass;
			else {
				q->fDensity = rs*p->fMass;
				TYPESet(q, TYPE_DensZeroed);
				}
			qiActive |= q->iActive;
			}
		if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
		}
	}

void initParticleMarkIIDensity(void *p)
{
	if (TYPEFilter((PARTICLE *) p,TYPE_DensACTIVE|TYPE_DensZeroed,
				   TYPE_DensACTIVE)) {
		((PARTICLE *)p)->fDensity = 0.0;
		TYPESet((PARTICLE *)p,TYPE_DensZeroed);
/*		if (((PARTICLE *)p)->iOrder == CHECKSOFT) fprintf(stderr,"Init Zero Particle 3031A: %g \n",((PARTICLE *) p)->fDensity);*/
		}
	}
/* copies of remote particles */
void initMarkIIDensity(void *p)
    {
    ((PARTICLE *) p)->fDensity = 0.0;
/*    if (((PARTICLE *)p)->iOrder == CHECKSOFT) fprintf(stderr,"Init Cache Zero Particle 3031A: %g \n",((PARTICLE *) p)->fDensity);*/
    }

void combMarkIIDensity(void *p1,void *p2)
    {
    if (TYPETest((PARTICLE *) p1,TYPE_DensACTIVE)) {
	if (TYPETest((PARTICLE *) p1,TYPE_DensZeroed)) 
	    ((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	else if (TYPETest((PARTICLE *) p2,TYPE_DensZeroed)) {
	    ((PARTICLE *)p1)->fDensity = ((PARTICLE *)p2)->fDensity;
	    TYPESet((PARTICLE *)p1,TYPE_DensZeroed);
	    }
	}
    ((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
    }

void MarkIIDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    assert(0);
    }

void MarkIIDensitySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
    {
    PARTICLE *q;
    FLOAT fNorm,ih2,r2,rs;
    int i;
    unsigned int qiActive;

    ih2 = 4.0/(BALL2(p));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    if (TYPETest(p,TYPE_DensACTIVE)) {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    p->fDensity += rs*q->fMass;
/*	    if (p->iOrder == CHECKSOFT) fprintf(stderr,"DensActive Particle %iA: %g %i  %g\n",p->iOrder,p->fDensity,q->iOrder,q->fMass);*/
	    if (TYPETest(q,TYPE_DensACTIVE)) {
		if (TYPETest(q,TYPE_DensZeroed)) {
		    q->fDensity += rs*p->fMass;
/*		    if (q->iOrder == CHECKSOFT) fprintf(stderr,"qDensActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		    }
		else {
		    q->fDensity = rs*p->fMass;
		    TYPESet(q,TYPE_DensZeroed);
/*		    if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qDensActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		    }
		}
	    if (TYPETest(p,TYPE_ACTIVE)) TYPESet(q,TYPE_NbrOfACTIVE);
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    else if (TYPETest(p,TYPE_ACTIVE)) {
	TYPESet( p,TYPE_NbrOfACTIVE);
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    TYPESet(q,TYPE_NbrOfACTIVE);
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) {
		q->fDensity += rs*p->fMass;
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"qActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q,TYPE_DensZeroed);
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qActive Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    }
	}
    else {
	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    qiActive |= q->iActive;
	    if (!TYPETest(q,TYPE_DensACTIVE)) continue;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs *= fNorm;
	    if (TYPETest(q,TYPE_DensZeroed)) {
		q->fDensity += rs*p->fMass;
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"qOther Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    else {
		q->fDensity = rs*p->fMass;
		TYPESet(q,TYPE_DensZeroed);
/*		if (q->iOrder == CHECKSOFT) fprintf(stderr,"zero qOther Particle %iA: %g %i \n",q->iOrder,q->fDensity,p->iOrder);*/
		}
	    }
	if (qiActive & TYPE_ACTIVE) TYPESet(p,TYPE_NbrOfACTIVE);
	}
    }


void initMark(void *p)
{
        }

void combMark(void *p1,void *p2)
{
	((PARTICLE *)p1)->iActive |= ((PARTICLE *)p2)->iActive;
	}


void initDeltaAccel(void *p)
{
	}

void combDeltaAccel(void *p1,void *p2)
{
    if (TYPEQueryACTIVE((PARTICLE *) p1) && ((PARTICLE *)p2)->dt < ((PARTICLE *)p1)->dt) 
	    ((PARTICLE *)p1)->dt = ((PARTICLE *)p2)->dt; 
    }

void DeltaAccel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	int i;
    FLOAT dax,da2,r2,dt;
	PARTICLE *q;

#ifdef DELTACCELCAP
	FLOAT pSoft2,qSoft2;
	pSoft2 = p->fSoft*p->fSoft;
#endif
	/*	assert(TYPEQueryACTIVE((PARTICLE *) p)); */

	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2;
	    if (r2 > 0) {
		  q = nnList[i].pPart;
#ifdef DELTAACCELACTIVE
		  if (!TYPEQueryACTIVE((PARTICLE *) q)) continue;
#endif
		  dax = p->a[0]-q->a[0];
		  da2 = dax*dax;
		  dax = p->a[1]-q->a[1];
		  da2 += dax*dax;
		  dax = p->a[2]-q->a[2];
		  da2 += dax*dax;
		  if (da2 > 0) {
#ifdef DELTACCELCAP
			if (r2 < pSoft2) r2 = pSoft2;
			qSoft2 = q->fSoft*q->fSoft;
			if (r2 < qSoft2) r2 = qSoft2;
#endif

   		    dt = smf->dDeltaAccelFac*sqrt(sqrt(r2/da2));  /* Timestep dt = Eta sqrt(deltar/deltaa) */
		    if (dt < p->dt) p->dt = dt;
		    if (
#ifndef DELTAACCELACTIVE
				TYPEQueryACTIVE((PARTICLE *) q) && 
#endif
				(dt < q->dt)) q->dt = dt;
		    }
		  }
	    }
    }

void initSinkAccrete(void *p)
{
	}

void combSinkAccrete(void *p1,void *p2)
{
    if (!(TYPETest( ((PARTICLE *) p1), TYPE_DELETED )) &&
        TYPETest( ((PARTICLE *) p2), TYPE_DELETED ) ) {
		((PARTICLE *) p1)-> fMass = ((PARTICLE *) p2)-> fMass;
	    pkdDeleteParticle( NULL, p1 );
	    }
    }

void SinkAccrete(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	int i;
	double dSinkRadius2 = smf->dSinkRadius*smf->dSinkRadius, 
	       EBO,Eq,r2,dvx,dv2,ifMass;
	PARTICLE *q;

	/* G = 1 
	 p is sink particle
	 q is gas particle */
	EBO = -0.5*p->fMass/smf->dSinkBoundOrbitRadius;

	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2;
        if (r2 > 0 && r2 <= dSinkRadius2) {
		  q = nnList[i].pPart;
		  if (TYPETest( q, TYPE_GAS )) {
			dvx = p->v[0]-q->v[0];
			dv2 = dvx*dvx;
			dvx = p->v[1]-q->v[1];
			dv2 += dvx*dvx;
			dvx = p->v[2]-q->v[2];
			dv2 += dvx*dvx;
			Eq = -p->fMass/sqrt(r2) + 0.5*dv2;
#ifdef GASOLINE
			if (smf->bSinkThermal) Eq+= q->u;
#endif
			if (Eq < EBO) {
			   ifMass = 1./(p->fMass + q->fMass);
			   p->r[0] = ifMass*(p->fMass*p->r[0]+q->fMass*q->r[0]);
			   p->r[1] = ifMass*(p->fMass*p->r[1]+q->fMass*q->r[1]);
			   p->r[2] = ifMass*(p->fMass*p->r[2]+q->fMass*q->r[2]);
			   p->v[0] = ifMass*(p->fMass*p->v[0]+q->fMass*q->v[0]);
			   p->v[1] = ifMass*(p->fMass*p->v[1]+q->fMass*q->v[1]);
			   p->v[2] = ifMass*(p->fMass*p->v[2]+q->fMass*q->v[2]);
			   p->a[0] = ifMass*(p->fMass*p->a[0]+q->fMass*q->a[0]);
			   p->a[1] = ifMass*(p->fMass*p->a[1]+q->fMass*q->a[1]);
			   p->a[2] = ifMass*(p->fMass*p->a[2]+q->fMass*q->a[2]);
			   p->fMass += q->fMass;
			   assert(q->fMass != 0);
			   q->fMass = 0;
			   pkdDeleteParticle(smf->pkd, q);
			   }
		    }
          }   
	    }
    }

/* Cached Tree Active particles */
void initBHSinkAccrete(void *p)
{
#ifdef GASOLINE
    if (TYPEQueryTREEACTIVE((PARTICLE *) p))
	((PARTICLE *)p)->u = 0.0;
#endif
    }

void combBHSinkAccrete(void *p1,void *p2)
{
#ifdef GASOLINE
    if (!(TYPETest( ((PARTICLE *) p1), TYPE_DELETED )) &&
        TYPETest( ((PARTICLE *) p2), TYPE_DELETED ) ) {
	((PARTICLE *) p1)-> fMass = ((PARTICLE *) p2)-> fMass;
	pkdDeleteParticle( NULL, p1 );
	    }
    else if (TYPEQueryTREEACTIVE((PARTICLE *) p1)) {
	((PARTICLE *)p1)->u += ((PARTICLE *)p2)->u;
	}
#endif
}

void BHSinkAccrete(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef GASOLINE
	PARTICLE *q = NULL;

	FLOAT ih2,r2,rs,fDensity;
	FLOAT v[3],cs,fW,dv2,dv;
	FLOAT mdot, mdotEdd, dm, dmq, dE, ifMass;
	int i;

	ih2 = 4.0/BALL2(p);
	fDensity = 0.0; cs = 0;
	v[0] = 0; v[1] = 0; v[2] = 0;
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    q = nnList[i].pPart;
	    /* Sink should NOT be part of this list */
	    assert(TYPETest(q,TYPE_GAS));
	    fW = rs*q->fMass;
	    fDensity += fW;
	    v[0] += fW*q->v[0];
	    v[1] += fW*q->v[1];
	    v[2] += fW*q->v[2];
	    cs += fW*q->c;
	    }
        dv2 = 0;
	for (i=0;i<3;i++) {
	    dv = v[i]/fDensity-p->v[i];
	    dv2 += dv*dv;
	    }
	cs = cs/fDensity;

        /* Scale density after using it to normalize averages */
        fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	printf("BHSink:  Density: %f C_s: %f dv: %f\n",fDensity,cs,sqrt(dv2));

        /* Bondi-Hoyle rate: cf. di Matteo et al 2005 (G=1) */
	mdot = smf->dBHSinkAlphaFactor*p->fMass*p->fMass*fDensity*pow(cs*cs+dv2,-1.5);
	/* Eddington Limit Rate */
	mdotEdd = smf->dBHSinkEddFactor*p->fMass;
	printf("BHSink:  mdot (BH): %f mdot (Edd): %f\n",mdot,mdotEdd);

	if (mdot > mdotEdd) mdot = mdotEdd;
	dm = mdot*smf->dSinkCurrentDelta;
	dE = smf->dBHSinkFeedbackFactor*dm;
	printf("BHSink:  Delta: %f dm: %f dE %f\n",smf->dSinkCurrentDelta,dm,dE);

	for (;;) {
	    FLOAT r2min = -1;
	    for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2;
		if (r2 < r2min && nnList[i].pPart->fMass > 0) {
		    r2min = r2;
		    q = nnList[i].pPart;
		    }
		}
	    assert( r2min >= 0 && q != p && q->fMass > 0);
	    /* We have our victim */
	    dmq = (dm > q->fMass ? q->fMass : dm);
	    ifMass = 1./(p->fMass + dmq);
	    /* Adjust sink properties (conserving momentum etc...) */
	    p->r[0] = ifMass*(p->fMass*p->r[0]+dmq*q->r[0]);
	    p->r[1] = ifMass*(p->fMass*p->r[1]+dmq*q->r[1]);
	    p->r[2] = ifMass*(p->fMass*p->r[2]+dmq*q->r[2]);
	    p->v[0] = ifMass*(p->fMass*p->v[0]+dmq*q->v[0]);
	    p->v[1] = ifMass*(p->fMass*p->v[1]+dmq*q->v[1]);
	    p->v[2] = ifMass*(p->fMass*p->v[2]+dmq*q->v[2]);
	    p->a[0] = ifMass*(p->fMass*p->a[0]+dmq*q->a[0]);
	    p->a[1] = ifMass*(p->fMass*p->a[1]+dmq*q->a[1]);
	    p->a[2] = ifMass*(p->fMass*p->a[2]+dmq*q->a[2]);
	    p->fMass += dmq;
	    q->fMass -= dmq;
	    dm -= dmq;
	    if (q->fMass < 1e-3*dmq) {
		q->fMass = 0;
		pkdDeleteParticle(smf->pkd, q);
		}
	    if (dm < 1e-3*dmq) break;
	    }   

	/* Recalculate Normalization */
	fDensity = 0.0; 
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    fDensity += rs*nnList[i].pPart->fMass;
	    }

	/* Low order: just adding energy directly to u */
	/* We should do sink calcs often enough so that du << u */
	fW = dE/fDensity;
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].pPart;
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    if (q->fMass > 0) {
		q->u += fW*rs*q->fMass;
		}
	    }
#endif
}


#ifdef SUPERCOOL
void initMeanVel(void *p)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p)->vMean[j] = 0.0;
		}
	}

void combMeanVel(void *p1,void *p2)
{
	int j;

	for (j=0;j<3;++j) {
		((PARTICLE *)p1)->vMean[j] += ((PARTICLE *)p2)->vMean[j];
		}
	}

void MeanVel(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/BALL2(p);
	fNorm = M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			}
		}
	}

void MeanVelSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs;
	int i,j;

	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		KERNEL(rs,r2);
		rs *= fNorm;
		q = nnList[i].pPart;
		for (j=0;j<3;++j) {
			p->vMean[j] += rs*q->fMass/q->fDensity*q->v[j];
			q->vMean[j] += rs*p->fMass/p->fDensity*p->v[j];
			}
		}
	}
#endif /* SUPER_COOL */


#ifdef GASOLINE
/* Original Particle */
void initSphPressureTermsParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef PDVDEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		}
	}

/* Cached copies of particle */
void initSphPressureTerms(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef PDVDEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combSphPressureTerms(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->PdV += ((PARTICLE *)p2)->PdV;
#ifdef PDVDEBUG
		((PARTICLE *)p1)->PdVvisc += ((PARTICLE *)p2)->PdVvisc;
		((PARTICLE *)p1)->PdVpres += ((PARTICLE *)p2)->PdVpres;
#endif
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only version -- untested */
void SphPressureTerms(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pPdV,pa[3],pmumax;
	FLOAT ph,pc,pDensity,visc,hav,vFac,absmu;
	FLOAT fNorm,fNorm1,fNorm2;
	int i;

	if (!TYPEQueryACTIVE(p)) return;

	pc = p->c;
	pDensity = p->fDensity;
	pPoverRho2 = p->PoverRho2;
	pmumax = p->mumax;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
	fNorm = M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	
	fNorm2 = fNorm1*(smf->a);    /* Comoving accelerations */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	pPdV=0.0;
	pa[0]=0.0;
	pa[1]=0.0;
	pa[2]=0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;
		if (dvdotdr>0.0) {
			pPdV += rs1 * PRES_PDV(pPoverRho2,q->PoverRho2) * dvdotdr;
			rs1 *= PRES_ACC(pPoverRho2,q->PoverRho2);
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
		else {
			hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu 
			   multiply by a to be consistent with physical c */
			absmu = -hav*dvdotdr*smf->a 
			    / (nnList[i].fDist2+0.01*hav*hav);
			/* mu terms for gas time step */
			if (absmu>pmumax) pmumax=absmu;
			/* viscosity term */

			visc = SWITCHCOMBINE(p,q)*
				(smf->alpha*(pc + q->c) + smf->beta*2*absmu) 
					*absmu/(pDensity + q->fDensity);
		        pPdV += rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
			rs1 *= PRES_ACC(pPoverRho2,q->PoverRho2) + visc;
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
 		}
	p->PdV += fNorm1*pPdV;
	p->mumax = pmumax;
	ACCEL(p,0) += fNorm2*pa[0];
	ACCEL(p,0) += fNorm2*pa[1];
	ACCEL(p,0) += fNorm2*pa[2];
	}

void SphPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pPdV,pa[3],pMass,pmumax;
	FLOAT ph,pc,pDensity,visc,hav,absmu;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

#ifdef PDVCHECK
	char ach[456];
#endif

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pPdV=0.0;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass;
				if (dvdotdr>0.0) {
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || !finite(rq * pPoverRho2 * dvdotdr) || !finite(rp * q->PoverRho2 * dvdotdr) || fabs(rq * pPoverRho2 * dvdotdr * 1e-5) > p->u || fabs(rp * q->PoverRho2 * dvdotdr *1e-5) > q->u) {
						sprintf(ach,"PDV-ERR-1 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g \n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
					pPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
					q->PdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#ifdef PDVDEBUG
					p->PdVpres += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
					q->PdVpres += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#endif
					rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
					rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
					rp *= aFac; /* convert to comoving acceleration */
					rq *= aFac;
					pa[0] -= rq * dx;
					pa[1] -= rq * dy;
					pa[2] -= rq * dz;
					ACCEL(q,0) += rp * dx;
					ACCEL(q,1) += rp * dy;
					ACCEL(q,2) += rp * dz;
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || fabs(rq) > 1e50 || fabs(rp) > 1e50 || fabs(ACCEL(p,0))+fabs(ACCEL(p,1))+fabs(ACCEL(p,2))+fabs(pa[0])+fabs(pa[1])+fabs(pa[2]) > 1e50 || fabs(ACCEL(q,0))+fabs(ACCEL(q,1))+fabs(ACCEL(q,2)) > 1e50) {
						sprintf(ach,"PDV-ACC-1 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
              		}
				else {
             		/* h mean - using just hp probably ok */
					hav=0.5*(ph+sqrt(0.25*BALL2(q)));
					/* mu multiply by a to be consistent with physical c */
					absmu = -hav*dvdotdr*smf->a 
						/(nnList[i].fDist2+0.01*hav*hav);
					/* mu terms for gas time step */
					if (absmu>pmumax) pmumax=absmu;
					if (absmu>q->mumax) q->mumax=absmu;
					/* viscosity term */

					visc = SWITCHCOMBINE(p,q)*
					  (smf->alpha*(pc + q->c) + smf->beta*2*absmu) 
					  *absmu/(pDensity + q->fDensity);
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || !finite(rq * (pPoverRho2 + 0.5*visc) * dvdotdr) || !finite(rp * (q->PoverRho2 + 0.5*visc) * dvdotdr) || fabs(rq * (pPoverRho2 + 0.5*visc) * dvdotdr * 1e-5) > p->u || fabs(rp * (q->PoverRho2 + 0.5*visc) * dvdotdr *1e-5) > q->u) {
						sprintf(ach,"PDV-ERR-2 %d - %d: Den %g - %g u %g %g PdV+ %g %g %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * (pPoverRho2 + 0.5*visc) * dvdotdr,rq * (pPoverRho2) * dvdotdr,rq * (0.5*visc) * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax);
						mdlDiag(smf->pkd->mdl,ach);
						sprintf(ach,"PDV-ERR-2 PdV %g %g %g Parts %g %g %g %g %g %g %g %g %g uPred %g %g %g %g %d %d \n",rq * (pPoverRho2 + 0.5*visc) * dvdotdr,rq * (pPoverRho2) * dvdotdr,rq * (0.5*visc) * dvdotdr, visc, dvdotdr, pc, q->c, absmu, hav, smf->a,p->BalsaraSwitch,q->BalsaraSwitch,p->uPred,q->uPred,p->uDot,q->uDot,p->iRung,q->iRung);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
					pPdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
					q->PdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2) + 0.5*visc)*dvdotdr;
#ifdef PDVDEBUG			
					p->PdVpres += rq*(PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
					q->PdVpres += rp*(PRES_PDV(q->PoverRho2,pPoverRho2))*dvdotdr;
					p->PdVvisc += rq*(0.5*visc)*dvdotdr;
					q->PdVvisc += rp*(0.5*visc)*dvdotdr;
#endif
					rq *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
					rp *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
					rp *= aFac; /* convert to comoving acceleration */
					rq *= aFac;

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz;
					ACCEL(q,0) += rp*dx;
					ACCEL(q,1) += rp*dy;
					ACCEL(q,2) += rp*dz;
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || fabs(rq) > 1e50 || fabs(rp) > 1e50 || fabs(ACCEL(p,0))+fabs(ACCEL(p,1))+fabs(ACCEL(p,2))+fabs(pa[0])+fabs(pa[1])+fabs(pa[2]) > 1e50 || fabs(ACCEL(q,0))+fabs(ACCEL(q,1))+fabs(ACCEL(q,2)) > 1e50) {
						sprintf(ach,"PDV-ACC-2 %d - %d: Den %g - %g u %g - %g PdV+ %g a %g %g %g %g %g %g vPred %g %g %g a %g %g %g vPred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,ACCEL(p,0),ACCEL(p,1),ACCEL(p,2),p->vPred[0],p->vPred[1],p->vPred[2],ACCEL(q,0),ACCEL(q,1),ACCEL(q,2),q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
              		}
				}
			else {
				/* q not active */
				if (dvdotdr>0.0) {
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || !finite(rq * pPoverRho2 * dvdotdr) || !finite(rp * q->PoverRho2 * dvdotdr) || fabs(rq * pPoverRho2 * dvdotdr * 1e-5) > p->u || fabs(rp * q->PoverRho2 * dvdotdr *1e-5) > q->u) {
						sprintf(ach,"PDV-ERR-3 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g \n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif

					pPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
#ifdef PDVDEBUG
  					p->PdVpres += rq*(PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
#endif
					rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
					rq *= aFac; /* convert to comoving acceleration */

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz;
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || fabs(rq) > 1e50 || fabs(ACCEL(p,0))+fabs(ACCEL(p,1))+fabs(ACCEL(p,2))+fabs(pa[0])+fabs(pa[1])+fabs(pa[2]) > 1e50) {
						sprintf(ach,"PDV-ACC-3 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
              		}
				else {
             		/* h mean */
					hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
					/* mu multiply by a to be consistent with physical c */
					absmu = -hav*dvdotdr*smf->a 
						/(nnList[i].fDist2+0.01*hav*hav);
					/* mu terms for gas time step */
					if (absmu>pmumax) pmumax=absmu;
					/* viscosity term */

					visc = SWITCHCOMBINE(p,q)*
					  (smf->alpha*(pc + q->c) + smf->beta*2*absmu) 
					  *absmu/(pDensity + q->fDensity);
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || !finite(rq * (pPoverRho2 + 0.5*visc) * dvdotdr) || !finite(rp * (q->PoverRho2 + 0.5*visc) * dvdotdr) || fabs(rq * (pPoverRho2 + 0.5*visc) * dvdotdr * 1e-5) > p->u || fabs(rp * (q->PoverRho2 + 0.5*visc) * dvdotdr *1e-5) > q->u) {
						sprintf(ach,"PDV-ERR-4 %d - %d: Den %g - %g u %g - %g PdV+ %g %g %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g \n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * (pPoverRho2 + 0.5*visc) * dvdotdr,rq * (pPoverRho2) * dvdotdr,rq * (0.5*visc) * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif

					pPdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
#ifdef PDVDEBUG			
  					p->PdVpres += rq*(PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
					p->PdVvisc += rq*(0.5*visc)*dvdotdr;
#endif
					rq *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
					rq *= aFac; /* convert to comoving acceleration */

					pa[0] -= rq*dx;
					pa[1] -= rq*dy;
					pa[2] -= rq*dz; 
#ifdef PDVCHECK
					if (p->iOrder==880556 || q->iOrder==880556 || fabs(rq) > 1e50 || fabs(ACCEL(p,0))+fabs(ACCEL(p,1))+fabs(ACCEL(p,2))+fabs(pa[0])+fabs(pa[1])+fabs(pa[2]) > 1e50 ) {
						sprintf(ach,"PDV-ACC-4 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
						mdlDiag(smf->pkd->mdl,ach);
						}
#endif
             		}
				}
	        }
		p->PdV += pPdV;
		p->mumax = pmumax;
		ACCEL(p,0) += pa[0];
		ACCEL(p,1) += pa[1];
		ACCEL(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */

	        r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;
			if (dvdotdr>0.0) {
#ifdef PDVCHECK
#endif
				q->PdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#ifdef PDVDEBUG			
				q->PdVpres += rp*(PRES_PDV(q->PoverRho2,pPoverRho2))*dvdotdr;
#endif

				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= aFac; /* convert to comoving acceleration */

		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
#ifdef PDVCHECK
				if (p->iOrder==880556 || q->iOrder==880556 || fabs(rp) > 1e50 || fabs(ACCEL(q,0))+fabs(ACCEL(q,1))+fabs(ACCEL(q,2)) > 1e50) {
					sprintf(ach,"PDV-ACC-5 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
					mdlDiag(smf->pkd->mdl,ach);
					}
#endif
				}
			else {
				/* h mean */
		        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu multiply by a to be consistent with physical c */
		        absmu = -hav*dvdotdr*smf->a 
					/(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity */

				visc = SWITCHCOMBINE(p,q)*
				  (smf->alpha*(pc + q->c) + smf->beta*2*absmu) 
				  *absmu/(pDensity + q->fDensity);
#ifdef PDVCHECK
#endif
				q->PdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2) + 0.5*visc)*dvdotdr;
#ifdef PDVDEBUG			
				q->PdVpres += rp*(PRES_PDV(q->PoverRho2,pPoverRho2))*dvdotdr;
				q->PdVvisc += rp*(0.5*visc)*dvdotdr;
#endif
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
				rp *= aFac; /* convert to comoving acceleration */

		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
#ifdef PDVCHECK
				if (p->iOrder==880556 || q->iOrder==880556 || fabs(rp) > 1e50 || fabs(ACCEL(q,0))+fabs(ACCEL(q,1))+fabs(ACCEL(q,2)) > 1e50) {
					sprintf(ach,"PDV-ACC-6 %d - %d: Den %g - %g u %g - %g PdV+ %g v %g %g %g Pred %g %g %g v %g %g %g Pred %g %g %g fB %g %g - %g %g a %g - %g\n",p->iOrder,q->iOrder,p->fDensity,q->fDensity,p->u,q->u,rq * pPoverRho2 * dvdotdr,p->v[0],p->v[1],p->v[2],p->vPred[0],p->vPred[1],p->vPred[2],q->v[0],q->v[1],q->v[2],q->vPred[0],q->vPred[1],q->vPred[2],sqrt(p->fBall2),p->fBallMax,sqrt(q->fBall2),q->fBallMax,rp,rq);
					mdlDiag(smf->pkd->mdl,ach);
					}
#endif
				}
	        }
		} 
	}

/* NB: ACCEL_PRES used here -- 
   with shock tracking disabled: #define NOSHOCKTRACK
   it is: a->aPres
   otherwise it is identical to p->a 
   The postSphPressure function combines p->a and p->aPres
*/

/* Original Particle */
void initSphPressureParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef PDVDEBUG
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		}
	}

/* Cached copies of particle */
void initSphPressure(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef PDVDEBUG
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		ACCEL_PRES(p,0) = 0.0;
		ACCEL_PRES(p,1) = 0.0;
		ACCEL_PRES(p,2) = 0.0;
		}
	}

void combSphPressure(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->PdV += ((PARTICLE *)p2)->PdV;
#ifdef PDVDEBUG
		((PARTICLE *)p1)->PdVpres += ((PARTICLE *)p2)->PdVpres;
#endif
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL_PRES(p1,0) += ACCEL_PRES(p2,0);
		ACCEL_PRES(p1,1) += ACCEL_PRES(p2,1);
		ACCEL_PRES(p1,2) += ACCEL_PRES(p2,2);
		}
	}

void postSphPressure(PARTICLE *p, SMF *smf)
{
        if ( TYPEQuerySMOOTHACTIVE((PARTICLE *)p) ) {
            if ( TYPEQueryACTIVE((PARTICLE *)p) ) {
                    ACCEL_COMB_PRES(p,0);
                    ACCEL_COMB_PRES(p,1);
                    ACCEL_COMB_PRES(p,2);
                    }
            }
	}

/* Gather only version -- untested */
void SphPressure(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pPdV,pa[3],pmumax;
	FLOAT ph,pc,pDensity,visc,hav,vFac,absmu;
	FLOAT fNorm,fNorm1,fNorm2;
	int i;

	assert(0);

	if (!TYPEQueryACTIVE(p)) return;

	pc = p->c;
	pDensity = p->fDensity;
	pPoverRho2 = p->PoverRho2;
	pmumax = p->mumax;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
	fNorm = M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	
	fNorm2 = fNorm1*(smf->a);    /* Comoving accelerations */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	pPdV=0.0;
	pa[0]=0.0;
	pa[1]=0.0;
	pa[2]=0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;
		if (dvdotdr>0.0) {
			pPdV += rs1 * PRES_PDV(pPoverRho2,q->PoverRho2) * dvdotdr;
#ifdef PDVDEBUG
			p->PdVpres += fNorm1*rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2))*dvdotdr;
#endif
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
		else {
			hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu 
			   multiply by a to be consistent with physical c */
			absmu = -hav*dvdotdr*smf->a 
			    / (nnList[i].fDist2+0.01*hav*hav);
			/* mu terms for gas time step */
			if (absmu>pmumax) pmumax=absmu;
			/* viscosity term */

			visc = SWITCHCOMBINE(p,q)*
			  (smf->alpha*(pc + q->c) + smf->beta*2*absmu) 
			  *absmu/(pDensity + q->fDensity);
		        pPdV += rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
#ifdef PDVDEBUG
			p->PdVpres += fNorm1*rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
#endif
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
			pa[0] -= rs1 * dx;
			pa[1] -= rs1 * dy;
			pa[2] -= rs1 * dz;
			}
 		}
	p->PdV += fNorm1*pPdV;
	p->mumax = pmumax;
	ACCEL_PRES(p,0) += fNorm2*pa[0];
	ACCEL_PRES(p,0) += fNorm2*pa[1];
	ACCEL_PRES(p,0) += fNorm2*pa[2];
	}

void SphPressureSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pPdV,pa[3],pMass,pmumax;
	FLOAT ph,pc,pDensity;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pPdV=0.0;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
			        rp = rs1 * pMass;
				pPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
				q->PdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#ifdef PDVDEBUG
				p->PdVpres += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
				q->PdVpres += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#endif
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= aFac; /* convert to comoving acceleration */
				rq *= aFac;
				pa[0] -= rq * dx;
				pa[1] -= rq * dy;
				pa[2] -= rq * dz;
				ACCEL_PRES(q,0) += rp * dx;
				ACCEL_PRES(q,1) += rp * dy;
				ACCEL_PRES(q,2) += rp * dz;
              		        }
			else {
				/* q not active */
			        pPdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
#ifdef PDVDEBUG
			        p->PdVpres += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
#endif
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rq *= aFac; /* convert to comoving acceleration */

				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz;
              		        }
             		}
		p->PdV += pPdV;
		p->mumax = pmumax;
		ACCEL_PRES(p,0) += pa[0];
		ACCEL_PRES(p,1) += pa[1];
		ACCEL_PRES(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	                q = nnList[i].pPart;
                        if (!TYPEQueryACTIVE(q)) continue; /* neither active */

                        r2 = nnList[i].fDist2*ih2;
                        DKERNEL(rs1,r2);
                        rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;
			q->PdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#ifdef PDVDEBUG
			q->PdVpres += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
#endif
			rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			rp *= aFac; /* convert to comoving acceleration */

		        ACCEL_PRES(q,0) += rp*dx;
		        ACCEL_PRES(q,1) += rp*dy;
		        ACCEL_PRES(q,2) += rp*dz;
	                }
		} 
	}

/* Original Particle */
void initSphViscosityParticle(void *p)
{
#ifdef PDVDEBUG
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		 ((PARTICLE *)p)->PdVvisc = 0.0;
	         }
#endif
}

/* Cached copies of particle */
void initSphViscosity(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef PDVDEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
#endif
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combSphViscosity(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->PdV += ((PARTICLE *)p2)->PdV;
#ifdef PDVDEBUG
		((PARTICLE *)p1)->PdVvisc += ((PARTICLE *)p2)->PdVvisc;
#endif
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only */
void SphViscosity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        assert(0);
	}

/* Symmetric Gather/Scatter version */
void SphViscositySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPdV,pa[3],pMass,pmumax;
	FLOAT ph,pc,pDensity,visc,hav,absmu;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	ph = sqrt(0.25*BALL2(p));
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);    /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE(p)) {
		/* p active */
		pmumax = p->mumax;
		pPdV=0.0;
		pa[0]=0.0;
		pa[1]=0.0;
		pa[2]=0.0;
		for (i=0;i<nSmooth;++i) {
			q = nnList[i].pPart;
			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rq = rs1 * q->fMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;
			if (dvdotdr > 0.0) continue;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass;
				/* h mean - using just hp probably ok */
				hav=0.5*(ph+sqrt(0.25*BALL2(q)));
				/* mu multiply by a to be consistent with physical c */
				absmu = -hav*dvdotdr*smf->a 
				  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>pmumax) pmumax=absmu;
				if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity term */

				visc = 
				  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c) 
				   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu) 
				  *absmu/(pDensity + q->fDensity);

				pPdV += rq*0.5*visc*dvdotdr;
				q->PdV += rp*0.5*visc*dvdotdr;
#ifdef PDVDEBUG
				p->PdVvisc += rq*0.5*visc*dvdotdr;
				q->PdVvisc += rp*0.5*visc*dvdotdr;
#endif
				rq *= visc;
				rp *= visc;
				rp *= aFac; /* convert to comoving acceleration */
				rq *= aFac;
				
				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz;
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			else {
				/* q not active */
			        /* h mean */
			        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
				/* mu multiply by a to be consistent with physical c */
				absmu = -hav*dvdotdr*smf->a 
				  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
				if (absmu>pmumax) pmumax=absmu;
				/* viscosity term */

				visc = 
				  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c) 
				   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu) 
				  *absmu/(pDensity + q->fDensity);
				
				pPdV += rq*0.5*visc*dvdotdr;
#ifdef PDVDEBUG
				p->PdVvisc += rq*0.5*visc*dvdotdr;
#endif
				rq *= visc;
				rq *= aFac; /* convert to comoving acceleration */
				
				pa[0] -= rq*dx;
				pa[1] -= rq*dy;
				pa[2] -= rq*dz; 
				}
	                }
		p->PdV += pPdV;
		p->mumax = pmumax;
		ACCEL(p,0) += pa[0];
		ACCEL(p,1) += pa[1];
		ACCEL(p,2) += pa[2];
		}
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
           	        q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */

	                r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm1;
			rp = rs1 * pMass;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;
			if (dvdotdr > 0.0) continue;

			/* h mean */
		        hav = 0.5*(ph+sqrt(0.25*BALL2(q)));
			/* mu multiply by a to be consistent with physical c */
		        absmu = -hav*dvdotdr*smf->a 
			  /(nnList[i].fDist2+0.01*hav*hav);
				/* mu terms for gas time step */
			if (absmu>q->mumax) q->mumax=absmu;
				/* viscosity */

			visc = 
			  (SWITCHCOMBINEA(p,q)*smf->alpha*(pc + q->c) 
			   + SWITCHCOMBINEB(p,q)*smf->beta*2*absmu) 
			  *absmu/(pDensity + q->fDensity);

			q->PdV += rp*0.5*visc*dvdotdr;
#ifdef PDVDEBUG
			q->PdVvisc += rp*0.5*visc*dvdotdr;
#endif
			rp *= visc;
			rp *= aFac; /* convert to comoving acceleration */

			ACCEL(q,0) += rp*dx;
			ACCEL(q,1) += rp*dy;
			ACCEL(q,2) += rp*dz;
	                }
		} 
	}



void initDivVort(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p )) {
		((PARTICLE *)p)->divv = 0.0;
		((PARTICLE *)p)->curlv[0] = 0.0;
		((PARTICLE *)p)->curlv[1] = 0.0;
		((PARTICLE *)p)->curlv[2] = 0.0;
#ifdef SHOCKTRACK
		((PARTICLE *)p)->divrhov = 0.0;
		((PARTICLE *)p)->gradrho[0] = 0.0;
		((PARTICLE *)p)->gradrho[1] = 0.0;
		((PARTICLE *)p)->gradrho[2] = 0.0;
#endif
		}
	}

void combDivVort(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1 )) {
		((PARTICLE *)p1)->divv += ((PARTICLE *)p2)->divv;
		((PARTICLE *)p1)->curlv[0] += ((PARTICLE *)p2)->curlv[0];
		((PARTICLE *)p1)->curlv[1] += ((PARTICLE *)p2)->curlv[1];
		((PARTICLE *)p1)->curlv[2] += ((PARTICLE *)p2)->curlv[2];
#ifdef SHOCKTRACK
		((PARTICLE *)p1)->divrhov += ((PARTICLE *)p2)->divrhov;
		((PARTICLE *)p1)->gradrho[0] += ((PARTICLE *)p2)->gradrho[0];
		((PARTICLE *)p1)->gradrho[1] += ((PARTICLE *)p2)->gradrho[1];
		((PARTICLE *)p1)->gradrho[2] += ((PARTICLE *)p2)->gradrho[2];
#endif
		}
	}

/* Gather only version -- untested */
void DivVort(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pcurlv[3],pdivv;
	FLOAT pDensity;
	FLOAT fNorm,vFac,a2;
	int i;

	if (!TYPEQueryACTIVE(p)) return;

	pDensity = p->fDensity;
	ih2 = 4.0/BALL2(p);
	a2 = (smf->a*smf->a);
	fNorm = M_1_PI*ih2*ih2; 
	vFac = (smf->bCannonical ? 1./a2 : 1.0); /* converts v to xdot */

	pdivv=0.0;
	pcurlv[0]=0.0;
	pcurlv[1]=0.0;
	pcurlv[2]=0.0;
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= q->fMass/pDensity;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;
		pdivv += rs1*dvdotdr;
		pcurlv[0] += rs1*(dvz*dy - dvy*dz);
		pcurlv[1] += rs1*(dvx*dz - dvz*dx);
		pcurlv[2] += rs1*(dvy*dx - dvx*dy);
 		}
	p->divv -=  fNorm*pdivv;  /* physical */
	p->curlv[0] += fNorm*vFac*pcurlv[0];
	p->curlv[1] += fNorm*vFac*pcurlv[1];
	p->curlv[2] += fNorm*vFac*pcurlv[2];
	}

/* Output is physical divv and curlv -- thus a*h_co*divv is physical */
void DivVortSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pMass,pDensity;
	FLOAT fNorm,dv,vFac,a2;
	int i;
 
	mdlassert(smf->pkd->mdl, TYPETest(p,TYPE_GAS));
	
	pDensity = p->fDensity;
	pMass = p->fMass;
	ih2 = 4.0/BALL2(p);
	a2 = (smf->a*smf->a);
	fNorm = 0.5*M_1_PI*ih2*ih2*sqrt(ih2); 
	vFac = (smf->bCannonical ? 1./a2 : 1.0); /* converts v to xdot */

	if (TYPEQueryACTIVE( p )) {
		/* p active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
		mdlassert(smf->pkd->mdl, TYPETest(q,TYPE_GAS));
	        r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *= fNorm;
			rq = rs1 * q->fMass/pDensity;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) +
				nnList[i].fDist2*smf->H;

			if (TYPEQueryACTIVE(q)) {
				/* q active */
				rp = rs1 * pMass/q->fDensity;
				p->divv -= rq*dvdotdr;
				q->divv -= rp*dvdotdr;
				dv=vFac*(dvz*dy - dvy*dz);
				p->curlv[0] += rq*dv;
				q->curlv[0] += rp*dv;
				dv=vFac*(dvx*dz - dvz*dx);
				p->curlv[1] += rq*dv;
				q->curlv[1] += rp*dv;
				dv=vFac*(dvy*dx - dvx*dy);
				p->curlv[2] += rq*dv;
				q->curlv[2] += rp*dv;

#ifdef SHOCKTRACK
				p->divrhov -= rs1*dvdotdr*q->fMass;
				q->divrhov -= rs1*dvdotdr*pMass;
				p->gradrho[0] += rs1*q->fMass*dx;
				q->gradrho[0] -= rs1*pMass*dx;
				p->gradrho[1] += rs1*q->fMass*dy;
				q->gradrho[1] -= rs1*pMass*dy;
				p->gradrho[2] += rs1*q->fMass*dz;
				q->gradrho[2] -= rs1*pMass*dz;
#endif
		        }
			else {
		        /* q inactive */
				p->divv -= rq*dvdotdr;
				dv=vFac*(dvz*dy - dvy*dz);
				p->curlv[0] += rq*dv;
				dv=vFac*(dvx*dz - dvz*dx);
				p->curlv[1] += rq*dv;
				dv=vFac*(dvy*dx - dvx*dy);
				p->curlv[2] += rq*dv;

#ifdef SHOCKTRACK
				p->divrhov -= rs1*dvdotdr*q->fMass;
				p->gradrho[0] += rs1*q->fMass*dx;
				p->gradrho[1] += rs1*q->fMass*dy;
				p->gradrho[2] += rs1*q->fMass*dz;
#endif
		        }
	        }
		} 
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
	        q = nnList[i].pPart;
		mdlassert(smf->pkd->mdl, TYPETest(q,TYPE_GAS));
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */

			r2 = nnList[i].fDist2*ih2;
			DKERNEL(rs1,r2);
			rs1 *=fNorm;
			rp = rs1*pMass/q->fDensity;

			dx = nnList[i].dx;
			dy = nnList[i].dy;
			dz = nnList[i].dz;
			dvx = p->vPred[0] - q->vPred[0];
			dvy = p->vPred[1] - q->vPred[1];
			dvz = p->vPred[2] - q->vPred[2];
			dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
				+ nnList[i].fDist2*smf->H;
			/* q active */
			q->divv -= rp*dvdotdr;
			dv=vFac*(dvz*dy - dvy*dz);
			q->curlv[0] += rp*dv;
			dv=vFac*(dvx*dz - dvz*dx);
			q->curlv[1] += rp*dv;
			dv=vFac*(dvy*dx - dvx*dy);
			q->curlv[2] += rp*dv;

#ifdef SHOCKTRACK
			q->divrhov -= rs1*dvdotdr*pMass;
			q->gradrho[0] -= rs1*pMass*dx;
			q->gradrho[1] -= rs1*pMass*dy;
			q->gradrho[2] -= rs1*pMass*dz;
#endif
	        }
		} 
	}

void initShockTrack(void *p)
{
#ifdef SHOCKTRACK
	if (TYPEQueryACTIVE((PARTICLE *) p )) {
		((PARTICLE *)p)->divrhov = 0.0;
		}
#endif
	}

void combShockTrack(void *p1,void *p2)
{
#ifdef SHOCKTRACK
	if (TYPEQueryACTIVE((PARTICLE *) p1 )) {
	        if (((PARTICLE *)p2)->divrhov > ((PARTICLE *)p2)->divrhov) 
		  ((PARTICLE *)p1)->divrhov = ((PARTICLE *)p2)->divrhov;
		}
#endif
	}

/* Gather only version -- untested */
void ShockTrack(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	assert(0);
	}

/* Output is physical divv and curlv -- thus a*h_co*divv is physical */
void ShockTrackSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
#ifdef SHOCKTRACK
	PARTICLE *q;
	double Mach,dv;
	int i,j;

	if (TYPEQueryACTIVE( p )) {
		/* p active */
		for (i=0;i<nSmooth;++i) {
	                q = nnList[i].pPart;
			Mach = 0;
			for (j=0;j<3;j++) {
			  dv = p->vPred[j] - q->gradrho[j];
			  Mach += dv*dv;
			}
			Mach = sqrt(Mach)/p->c;
			if (Mach > p->divrhov) p->divrhov = Mach;
			if (TYPEQueryACTIVE(q)) {
			  Mach = 0;
			  for (j=0;j<3;j++) {
			    dv = q->vPred[j] - p->gradrho[j];
			    Mach += dv*dv;
			  }
			  Mach = sqrt(Mach)/q->c;
			  if (Mach > q->divrhov) q->divrhov = Mach;
		        }
	        }
		} 
	else {
		/* p not active */
		for (i=0;i<nSmooth;++i) {
 	                q = nnList[i].pPart;
			if (!TYPEQueryACTIVE(q)) continue; /* neither active */
			Mach = 0;
			for (j=0;j<3;j++) {
			  dv = q->vPred[j] - p->gradrho[j];
			  Mach += dv*dv;
			}
			Mach = sqrt(Mach)/q->c;
			if (Mach > q->divrhov) q->divrhov = Mach;
	        }
		} 
#endif
	}

/* Original Particle */
void initHKPressureTermsParticle(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef DEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		}
	}

/* Cached copies of particle */
void initHKPressureTerms(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
#ifdef DEBUG
		((PARTICLE *)p)->PdVvisc = 0.0;
		((PARTICLE *)p)->PdVpres = 0.0;
#endif
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combHKPressureTerms(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->PdV += ((PARTICLE *)p2)->PdV;
#ifdef DEBUG
		((PARTICLE *)p1)->PdVvisc += ((PARTICLE *)p2)->PdVvisc;
		((PARTICLE *)p1)->PdVpres += ((PARTICLE *)p2)->PdVpres;
#endif
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only version -- (untested)  */
void HKPressureTerms(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,pMass,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	if (!TYPEQueryACTIVE(p)) return;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity);
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1 * q->fMass;;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr>0.0) {
			p->PdV += rs1*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2));
			rs1 *= aFac;
			ACCEL(p,0) -= rs1*dx;
			ACCEL(p,1) -= rs1*dy;
			ACCEL(p,2) -= rs1*dz;
			}
		else {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 : 
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);
			visc = pQonRho2 + qQonRho2;
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh+ph);
			absmu = -hav*dvdotdr*smf->a
				/(nnList[i].fDist2+0.01*hav*hav);
			if (absmu>p->mumax) p->mumax=absmu;
			p->PdV += rs1 * (PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc) * dvdotdr;
			rs1 *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
			rs1 *= aFac; /* convert to comoving acceleration */
			ACCEL(p,0) -= rs1*dx;
			ACCEL(p,1) -= rs1*dy;
			ACCEL(p,2) -= rs1*dz;
			}
		}
	}

/* Bulk viscosity and standard pressure forces */
void HKPressureTermsSym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pPoverRho2,pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,pMass,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	pPoverRho2 = p->PoverRho2;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity );
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1;
		rq = rs1*q->fMass;
		rp = rs1*pMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr>0.0) {
			if (TYPEQueryACTIVE(p)) {
		                p->PdV += rq*PRES_PDV(pPoverRho2,q->PoverRho2)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rq *= aFac;
				ACCEL(p,0) -= rq*dx;
				ACCEL(p,1) -= rq*dy;
				ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				q->PdV += rp*PRES_PDV(q->PoverRho2,pPoverRho2)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2));
				rp *= aFac; /* convert to comoving acceleration */
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			}
		else {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 : 
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);
			visc = pQonRho2 + qQonRho2;
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh + ph);
			absmu = -hav*dvdotdr*smf->a/(nnList[i].fDist2+0.01*hav*hav);
			if (TYPEQueryACTIVE(p)) {
				if (absmu>p->mumax) p->mumax=absmu;
				p->PdV += rq*(PRES_PDV(pPoverRho2,q->PoverRho2) + 0.5*visc)*dvdotdr;
				rq *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
				rq *= aFac; /* convert to comoving acceleration */
				ACCEL(p,0) -= rq*dx;
				ACCEL(p,1) -= rq*dy;
				ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				if (absmu>q->mumax) q->mumax=absmu;
				q->PdV += rp*(PRES_PDV(q->PoverRho2,pPoverRho2) + 0.5*visc)*dvdotdr;
				rp *= (PRES_ACC(pPoverRho2,q->PoverRho2) + visc);
				rp *= aFac; /* convert to comoving acceleration */
				ACCEL(q,0) += rp*dx;
				ACCEL(q,1) += rp*dy;
				ACCEL(q,2) += rp*dz;
				}
			}
		}
	}

/* Original Particle */
void initHKViscosityParticle(void *p)
{
	ACCEL(p,0) += ACCEL_PRES(p,0);
	ACCEL(p,1) += ACCEL_PRES(p,1);
	ACCEL(p,2) += ACCEL_PRES(p,2);
	}

/* Cached copies of particle */
void initHKViscosity(void *p)
{
	if (TYPEQueryACTIVE((PARTICLE *) p)) {
		((PARTICLE *)p)->mumax = 0.0;
		((PARTICLE *)p)->PdV = 0.0;
		ACCEL(p,0) = 0.0;
		ACCEL(p,1) = 0.0;
		ACCEL(p,2) = 0.0;
		}
	}

void combHKViscosity(void *p1,void *p2)
{
	if (TYPEQueryACTIVE((PARTICLE *) p1)) {
		((PARTICLE *)p1)->PdV += ((PARTICLE *)p2)->PdV;
		if (((PARTICLE *)p2)->mumax > ((PARTICLE *)p1)->mumax)
			((PARTICLE *)p1)->mumax = ((PARTICLE *)p2)->mumax;
		ACCEL(p1,0) += ACCEL(p2,0);
		ACCEL(p1,1) += ACCEL(p2,1);
		ACCEL(p1,2) += ACCEL(p2,2);
		}
	}

/* Gather only version */
void HKViscosity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        assert(0);
	}

/* Bulk viscosity */
void HKViscositySym(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT ih2,r2,rs1,rq,rp;
	FLOAT dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	FLOAT pQonRho2,qQonRho2,qhdivv;
	FLOAT ph,pc,pDensity,visc,absmu,qh,pMass,hav;
	FLOAT fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c;
	pDensity = p->fDensity;
	pMass = p->fMass;
	ph = sqrt(0.25*BALL2(p));
	/* QonRho2 given same scaling with a as PonRho2 */
	pQonRho2 = (p->divv>0.0 ? 0.0 : fabs(p->divv)*ph*smf->a
				*(smf->alpha*pc + smf->beta*fabs(p->divv)*ph*smf->a)/pDensity );
	ih2 = 4.0/BALL2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = (smf->a);        /* comoving acceleration factor */
	vFac = (smf->bCannonical ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		r2 = nnList[i].fDist2*ih2;
		DKERNEL(rs1,r2);
		rs1 *= fNorm1;
		rq = rs1*q->fMass;
		rp = rs1*pMass;

		dx = nnList[i].dx;
		dy = nnList[i].dy;
		dz = nnList[i].dz;
		dvx = p->vPred[0] - q->vPred[0];
		dvy = p->vPred[1] - q->vPred[1];
		dvz = p->vPred[2] - q->vPred[2];
		dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + nnList[i].fDist2*smf->H;

		if (dvdotdr<0.0) {
			qh=sqrt(0.25*BALL2(q));
			qhdivv = qh*fabs(q->divv)*smf->a; /* units of physical velocity */
			qQonRho2 = (qhdivv>0.0 ? 0.0 : 
						qhdivv*(smf->alpha*q->c + smf->beta*qhdivv)/q->fDensity);

			visc = SWITCHCOMBINE(p,q)*(pQonRho2 + qQonRho2);
			/* mu -- same timestep criteria as standard sph above (for now) */
			hav=0.5*(qh + ph);
			absmu = -hav*dvdotdr*smf->a/(nnList[i].fDist2+0.01*hav*hav);
			if (TYPEQueryACTIVE(p)) {
				if (absmu>p->mumax) p->mumax=absmu;
		        p->PdV += rq*0.5*visc*dvdotdr;
				rq *= visc;
				rq *= aFac; /* convert to comoving acceleration */
		        ACCEL(p,0) -= rq*dx;
		        ACCEL(p,1) -= rq*dy;
		        ACCEL(p,2) -= rq*dz;
		        }
			if (TYPEQueryACTIVE(q)) {
				if (absmu>q->mumax) q->mumax=absmu;
				q->PdV += rp*0.5*visc*dvdotdr;
				rp *= visc;
				rp *= aFac; /* convert to comoving acceleration */
		        ACCEL(q,0) += rp*dx;
		        ACCEL(q,1) += rp*dy;
		        ACCEL(q,2) += rp*dz;
				}
			}
		}
	}

#ifdef STARFORM
void initDistDeletedGas(void *p1)
{
	if(!TYPETest(((PARTICLE *)p1), TYPE_DELETED)) {
    /*
     * Zero out accumulated quantities.
     */
		((PARTICLE *)p1)->fMass = 0;
		((PARTICLE *)p1)->v[0] = 0;
		((PARTICLE *)p1)->v[1] = 0;
		((PARTICLE *)p1)->v[2] = 0;
		((PARTICLE *)p1)->u = 0;
		((PARTICLE *)p1)->uDot = 0.0;
		((PARTICLE *)p1)->fMetals = 0.0;
		((PARTICLE *)p1)->fMFracOxygen = 0.0;
		((PARTICLE *)p1)->fMFracIron = 0.0;
		}
    }

void combDistDeletedGas(void *vp1,void *vp2)
{
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.  
     */
	PARTICLE *p1 = vp1;
	PARTICLE *p2 = vp2;

	if(!TYPETest((p1), TYPE_DELETED)) {
		FLOAT delta_m = p2->fMass;
		FLOAT m_new,f1,f2;
		FLOAT fTCool; /* time to cool to zero */
		
		m_new = p1->fMass + delta_m;
		if (m_new > 0) {
			f1 = p1->fMass /m_new;
			f2 = delta_m  /m_new;
			if(p1->uDot < 0.0)
				fTCool = p1->uPred/p1->uDot; 
			
			p1->fMass = m_new;
			p1->u = f1*p1->u + f2*p2->u;
			p1->uPred = f1*p1->uPred + f2*p2->uPred;
#ifdef COOLDEBUG
			assert(p1->u >= 0.0);
#endif
			p1->v[0] = f1*p1->v[0] + f2*p2->v[0];            
			p1->v[1] = f1*p1->v[1] + f2*p2->v[1];            
			p1->v[2] = f1*p1->v[2] + f2*p2->v[2];            
			p1->fMetals = f1*p1->fMetals + f2*p2->fMetals;
			p1->fMFracOxygen = f1*p1->fMFracOxygen + f2*p2->fMFracOxygen;
			p1->fMFracIron = f1*p1->fMFracIron + f2*p2->fMFracIron;
			if(p1->uDot < 0.0)
				p1->uDot = p1->uPred/fTCool;
			}
		}
    }

void DistDeletedGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs,rstot,delta_m,m_new,f1,f2;
	FLOAT fTCool; /* time to cool to zero */
	int i;

	assert(TYPETest(p, TYPE_GAS));
	ih2 = 4.0/BALL2(p);
        rstot = 0;        
	for (i=0;i<nSmooth;++i) {
            q = nnList[i].pPart;
	    if(TYPETest(q, TYPE_DELETED)) continue;
	    assert(TYPETest(q, TYPE_GAS));
            r2 = nnList[i].fDist2*ih2;            
            KERNEL(rs,r2);
            rstot += rs;
        }
	if(rstot <= 0.0) {
	    if(p->fMass == 0.0) /* the particle to be deleted has NOTHING */
		return;
	    /* we have a particle to delete and nowhere to put its mass
	     * => we will keep it around */
	    pkdUnDeleteParticle(smf->pkd, p);
	    return;
	    }
	assert(rstot > 0.0);
	fNorm = 1./rstot;
	assert(p->fMass >= 0.0);
	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
		if(TYPETest(q, TYPE_DELETED)) continue;
		
		r2 = nnList[i].fDist2*ih2;            
		KERNEL(rs,r2);
	    /*
	     * All these quantities are per unit mass.
	     * Exact if only one gas particle being distributed or in serial
	     * Approximate in parallel (small error).
	     */
		delta_m = rs*fNorm*p->fMass;
		m_new = q->fMass + delta_m;
		/* Cached copies can have zero mass: skip them */
		if (m_new == 0) continue;
		f1 = q->fMass /m_new;
		f2 = delta_m  /m_new;
		q->fMass = m_new;
		if(q->uDot < 0.0)
			fTCool = q->uPred/q->uDot; 
		
                /* Only distribute the properties
                 * to the other particles on the "home" machine.
                 * u, v, and fMetals will be distributed to particles
                 * that come through the cache in the comb function.
                 */
		q->u = f1*q->u+f2*p->u;
		q->uPred = f1*q->uPred+f2*p->uPred;
#ifdef COOLDEBUG
		assert(q->u >= 0.0);
#endif
		q->v[0] = f1*q->v[0]+f2*p->v[0];            
		q->v[1] = f1*q->v[1]+f2*p->v[1];            
		q->v[2] = f1*q->v[2]+f2*p->v[2];            
		q->fMetals = f1*q->fMetals + f2*p->fMetals;
                q->fMFracOxygen = f1*q->fMFracOxygen + f2*p->fMFracOxygen;
                q->fMFracIron = f1*q->fMFracIron + f2*p->fMFracIron;
		if(q->uDot < 0.0) /* make sure we don't shorten cooling time */
			q->uDot = q->uPred/fTCool;
        }
}

void DeleteGas(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
/* flag low mass gas particles for deletion */
	PARTICLE *q;
	FLOAT fMasstot,fMassavg;
	int i;

	assert(TYPETest(p, TYPE_GAS));
	fMasstot = 0;
#ifdef COOLDEBUG
	assert(p->fMass >= 0.0);
#endif

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
	    assert(TYPETest(q, TYPE_GAS));
		fMasstot += q->fMass;
        }

	fMassavg = fMasstot/(FLOAT) nSmooth;

	if (p->fMass < smf->dMinMassFrac*fMassavg) {
		pkdDeleteParticle(smf->pkd, p);
        }
	else {
		assert (p->fMass > 0.0);
		}
        
}

void initTreeParticleDistSNEnergy(void *p1)
{
    /* Convert energy and metals to non-specific quantities (not per mass)
     * to make it easier to divvy up SN energy and metals.  
     */
    
    if(TYPETest((PARTICLE *)p1, TYPE_GAS)){
        ((PARTICLE *)p1)->fESNrate *= ((PARTICLE *)p1)->fMass;
        ((PARTICLE *)p1)->fMetals *= ((PARTICLE *)p1)->fMass;    
        ((PARTICLE *)p1)->fMFracOxygen *= ((PARTICLE *)p1)->fMass;    
        ((PARTICLE *)p1)->fMFracIron *= ((PARTICLE *)p1)->fMass;    
        }
    
    }

void initDistSNEnergy(void *p1)
{
    /*
     * Warning: kludgery.  We need to accumulate mass in the cached
     * particle, but we also need to keep the original mass around.
     * Let's use another (hopefully unused) field to hold the original
     * mass.
     */
    ((PARTICLE *)p1)->PdV = ((PARTICLE *)p1)->fMass;

    /*
     * Zero out accumulated quantities.
     */
    ((PARTICLE *)p1)->fESNrate = 0.0;
    ((PARTICLE *)p1)->fMetals = 0.0;
    ((PARTICLE *)p1)->fMFracOxygen = 0.0;
    ((PARTICLE *)p1)->fMFracIron = 0.0;
    }

void combDistSNEnergy(void *p1,void *p2)
{
    /*
     * See kludgery notice above.
     */
    FLOAT fAddedMass = ((PARTICLE *)p2)->fMass - ((PARTICLE *)p2)->PdV;
    
    ((PARTICLE *)p1)->fMass += fAddedMass;
    ((PARTICLE *)p1)->fESNrate += ((PARTICLE *)p2)->fESNrate;
    ((PARTICLE *)p1)->fMetals += ((PARTICLE *)p2)->fMetals;
    ((PARTICLE *)p1)->fMFracOxygen += ((PARTICLE *)p2)->fMFracOxygen;
    ((PARTICLE *)p1)->fMFracIron += ((PARTICLE *)p2)->fMFracIron;
    ((PARTICLE *)p1)->fTimeCoolIsOffUntil = max( ((PARTICLE *)p1)->fTimeCoolIsOffUntil,
                ((PARTICLE *)p2)->fTimeCoolIsOffUntil );
    }

void DistSNEnergy(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_Pres,fAveDens,fNewBall,f2h2;
        FLOAT delta_m,m_new,f1,f2,fBlastRadius,fShutoffTime,fmind;
	int i,counter,imind;

	if ( (p->fMSN == 0.0) && (p->fESNrate == 0.0)){return;}
	assert(TYPETest(p, TYPE_STAR));
	ih2 = 4.0/BALL2(p);
        f2h2=BALL2(p);
        rstot = 0.0;  
        fNorm_u = 0.0;
        fNorm_Pres = 0.0;
	fAveDens = 0.0;
	
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
            r2 = nnList[i].fDist2*ih2;            
            KERNEL(rs,r2);
            q = nnList[i].pPart;
            fNorm_u += q->fMass*rs;
            rs *= fNorm;
            fAveDens += q->fMass*rs;
            fNorm_Pres += q->fMass*q->uPred*rs;
	    assert(TYPETest(q, TYPE_GAS));
            }
        fNorm_Pres *= (smf->gamma-1.0);
        /* from McKee and Ostriker (1977) ApJ 218 148 */
        fBlastRadius = smf->dRadPreFactor*pow(p->fNSN,0.32)*pow(fAveDens,-0.16)*pow(fNorm_Pres,-0.2);
	if (smf->bShortCoolShutoff){
        /* End of snowplow phase */
            fShutoffTime = smf->dTimePreFactor*pow(p->fNSN,0.31)*pow(fAveDens,0.27)*pow(fNorm_Pres,-0.64);
	} else{        /* McKee + Ostriker 1977 t_{max} */
            fShutoffTime = smf->dTimePreFactor*pow(p->fNSN,0.32)*pow(fAveDens,0.34)*pow(fNorm_Pres,-0.70);
	}
	
   if (smf->bSmallSNSmooth) {
        fmind = BALL2(p);
        imind = 0;
        if ( p->fESNrate > 0.0 ) {
            /* Change smoothing radius to blast radius 
             * so that we only distribute mass, metals, and energy
             * over that range. 
             */
            f2h2 = fBlastRadius*fBlastRadius;
            ih2 = 4.0/f2h2;
            rstot = 0.0;  
            fNorm_u = 0.0;
            for (i=0;i<nSmooth;++i) {
                if ( nnList[i].fDist2 < fmind ){imind = i; fmind = nnList[i].fDist2;}
                if ( nnList[i].fDist2 < f2h2 ) {
                    r2 = nnList[i].fDist2*ih2;            
                    KERNEL(rs,r2);
                    q = nnList[i].pPart;
                    fNorm_u += q->fMass*rs;
                    assert(TYPETest(q, TYPE_GAS));
                    }
                }
            }
        
        /* If there's no gas particle within blast radius,
           give mass and energy to nearest gas particle. */
        if (fNorm_u ==0.0){
            r2 = nnList[imind].fDist2*ih2;            
            KERNEL(rs,r2);
	    /*
	     * N.B. This will be NEGATIVE, but that's OK since it will
	     * cancel out down below.
	     */
            fNorm_u = nnList[imind].pPart->fMass*rs;
            }
	 }
	assert(fNorm_u != 0.0);
        fNorm_u = 1./fNorm_u;
counter=0;
	for (i=0;i<nSmooth;++i) {
            q = nnList[i].pPart;
            if (smf->bSmallSNSmooth) {
                if ( (nnList[i].fDist2 <= f2h2) || (i == imind) ) {
                    if( p->fNSN > 0.0 && smf->bSNTurnOffCooling && 
                       (fBlastRadius*fBlastRadius >= nnList[i].fDist2)) {
                        q->fTimeCoolIsOffUntil = max(q->fTimeCoolIsOffUntil,
                            smf->dTime + fShutoffTime);}
counter++;  
                    r2 = nnList[i].fDist2*ih2;  
                    KERNEL(rs,r2);
                    /* Remember: We are dealing with total energy rate and total metal
                     * mass, not energy/gram or metals per gram.  
                     * q->fMass is in product to make units work for fNorm_u.
                     */
                    q->fESNrate += rs*fNorm_u*q->fMass*p->fESNrate;
                    q->fMetals += rs*fNorm_u*q->fMass*p->fSNMetals;
                    q->fMFracOxygen += rs*fNorm_u*q->fMass*p->fMOxygenOut;
                    q->fMFracIron += rs*fNorm_u*q->fMass*p->fMIronOut;
                    q->fMass += rs*fNorm_u*q->fMass*p->fMSN;
                    }
            } else {
                r2 = nnList[i].fDist2*ih2;  
                KERNEL(rs,r2);
                /* Remember: We are dealing with total energy rate and total metal
                 * mass, not energy/gram or metals per gram.  
                 * q->fMass is in product to make units work for fNorm_u.
                 */
                q->fESNrate += rs*fNorm_u*q->fMass*p->fESNrate;
                q->fMetals += rs*fNorm_u*q->fMass*p->fSNMetals;
                q->fMFracOxygen += rs*fNorm_u*q->fMass*p->fMOxygenOut;
                q->fMFracIron += rs*fNorm_u*q->fMass*p->fMIronOut;
                
                if ( p->fESNrate > 0.0 && smf->bSNTurnOffCooling && 
                     (fBlastRadius*fBlastRadius >= nnList[i].fDist2)){
                    q->fTimeCoolIsOffUntil = max(q->fTimeCoolIsOffUntil,
                        smf->dTime + fShutoffTime);
    counter++;
                    }
                /*	update mass after everything else so that distribution
                    is based entirely upon initial mass of gas particle */
                q->fMass += rs*fNorm_u*q->fMass*p->fMSN;
                } 
            }
/*if(counter>0) printf("%i ",counter);
if (p->fNSN >0) printf("%i E51:  %g  Dens:  %g  P:  %g  R:  %g shutoff time: %g   \n",counter,p->fNSN,fAveDens,fNorm_Pres,fBlastRadius,fShutoffTime);
if(p->fNSN!= 0.0)printf("E51:  %g  Dens:  %g  P:  %g  R:  %g shutoff time: %g  \n",p->fNSN,fAveDens,fNorm_Pres,fBlastRadius,fShutoffTime);*/
}

void postDistSNEnergy(PARTICLE *p1, SMF *smf)
{
    /* Convert energy and metals back to specific quantities (per mass)
       because we are done with our conservative calculations */
    
    if(TYPETest(p1, TYPE_GAS)){
        p1->fESNrate /= p1->fMass;
        p1->fMetals /= p1->fMass;    
        p1->fMFracIron /= p1->fMass;    
        p1->fMFracOxygen /= p1->fMass;    
        }
    
    }

#endif /* STARFORM */

#ifdef SIMPLESF
void initSimpleSF_Feedback(void *p1)
{
    /*
     * Zero out accumulated quantities.
     */
    ((PARTICLE *)p1)->u = 0.0;
    ((PARTICLE *)p1)->fMetals = 0.0;
    ((PARTICLE *)p1)->fTimeForm = 0.0;
    }

void combSimpleSF_Feedback(void *p1,void *p2)
{
    ((PARTICLE *)p1)->u += ((PARTICLE *)p2)->u;
    ((PARTICLE *)p1)->fMetals += ((PARTICLE *)p2)->fMetals;
    ((PARTICLE *)p1)->fTimeForm += ((PARTICLE *)p2)->fTimeForm;
    }

void SimpleSF_Feedback(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	PARTICLE *q;
	FLOAT fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_t;
	int i;

	assert(TYPETest(p, TYPE_STAR));
	ih2 = 4.0/BALL2(p);
	rstot = 0;        
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;            
		KERNEL(rs,r2);
		rstot += rs;
        }
	
	fNorm = 1./rstot;
	fNorm_u = fNorm*p->fMass*p->fESN;
	assert(fNorm_u > 0.0);

	fNorm_t = fNorm*p->PdV; /* p->PdV store the cooling delay dtCoolingShutoff */
	assert(fNorm_t > 0.0);

	for (i=0;i<nSmooth;++i) {
		q = nnList[i].pPart;
	    assert(TYPETest(q, TYPE_GAS));
		r2 = nnList[i].fDist2*ih2;            
		KERNEL(rs,r2);
		q->u += rs*fNorm_u/q->fMass;
		q->fMetals += rs*fNorm;
		q->fTimeForm += rs*fNorm_t;
        }
	}

#endif /* SIMPLESF */

#endif /* GASOLINE */

#ifdef COLLISIONS

void initFindRejects(void *p)
{
	((PARTICLE *) p)->dtCol = 0.0;
	}

void combFindRejects(void *p1,void *p2)
{
	if (REJECT((PARTICLE *) p2))
		((PARTICLE *) p1)->dtCol = -1.0; /* reject if dtCol < 0.0 */
	}

void FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Checks "nSmooth" neighbours of "p" for overlap (physical, Hill
	** radius, or a combination).  When the particles in question are
	** of different "size", the one with the smallest size is rejected
	** first, otherwise the one with the higher iOrder is rejected.
	** Note that only planetesimal neighbors not already rejected are
	** considered.  This procedure uses a combiner cache so that the
	** neighbours of "p" can be flagged with maximum efficiency.  Note
	** the number of rejects flagged may be different when running in
	** parallel (see msrFindRejects()).
	*/
	
	PARTICLE *pn;
	double r,rn,sr;
	int i;
	
#ifndef SLIDING_PATCH
	double a=0.0,r2,v2,an,rh;
#endif
	
	if (REJECT(p))
		return; /* already rejected */
	
	r = RADIUS(p);
	
#ifndef SLIDING_PATCH
	if (smf->dCentMass > 0.0) {
		r2 = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2];
		v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
		assert(r2 > 0.0); /* particle must not be at origin */
		a = 2.0/sqrt(r2) - v2/(smf->dCentMass + p->fMass);
		assert(a != 0.0); /* can't handle parabolic orbits */
		a = 1.0/a;
		}
#endif
	
	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		if (pn->iOrder == p->iOrder || REJECT(pn))
			continue;
		rn = RADIUS(pn);
#ifndef SLIDING_PATCH
		if (smf->dCentMass > 0.0) {
			r2 = pn->r[0]*pn->r[0] + pn->r[1]*pn->r[1] + pn->r[2]*pn->r[2];
			v2 = pn->v[0]*pn->v[0] + pn->v[1]*pn->v[1] + pn->v[2]*pn->v[2];
			assert(r2 > 0.0);
			an = 2.0/sqrt(r2) - v2/(smf->dCentMass + pn->fMass);
			assert(an != 0.0);
			an = 1.0/an;
			rh = 0.5*(a + an)*pow((p->fMass + pn->fMass)/(3.0*smf->dCentMass),1.0/3.0);
			if (rh > r)
				r = rh;
			if (rh > rn)
				rn = rh;
			}
#endif
		if (rn > r || (rn == r && pn->iOrder < p->iOrder))
			continue;
		sr = r + rn;
		if (nnList[i].fDist2 <= sr*sr)
			pn->dtCol = -1.0; /* see REJECT() macro */
		}
	}

void _CheckForCollapse(PARTICLE *p,double dt,double rdotv,double r2,SMF *smf)
{
	/*
	** Sets bTinyStep flag of particle "p" to 1 if distance travelled
	** since last collision small (factor "smf->dCollapseLimit")
	** compared to particle radius.  Here, "rdotv" is the dot product
	** of the relative position and relative velocity (note rdotv must
	** be negative), and "r2" is the square of the distance. This
	** routine should only be called by CheckForCollision().
	*/
  
	double dRatio;

	assert(rdotv <= 0.0 && p->dtPrevCol <= dt);
	dRatio = rdotv*(p->dtPrevCol - dt)/(sqrt(r2)*RADIUS(p));
	if (dRatio < smf->dCollapseLimit) {
#if (INTERNAL_WARNINGS != 0)
		static int nWarn = 1;
		if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
			(void) fprintf(stderr,"TINY STEP WARNING #%i (pid=%i) [T=%g]: %i & %i (dt=%g, dRatio=%g)\n",
						   nWarn,smf->pkd->idSelf,smf->dTime,p->iOrder,p->iOrderCol,dt,dRatio);
		++nWarn;
#endif /* INTERNAL_WARNINGS */
		p->bTinyStep = 1;
	}
}

#ifdef RORY_EXTENSION

void FindBinary(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/* This subroutine looks for binaries. Search only those particle
	 * on rungs higher than iMinBinaryRung. The particles have been 
	 * predetermined to be on high rungs which may result from being
	 * bound. The particles are considered bound if they have a total 
	 * energy < 0.
	 */

	PARTICLE *pn;
	int i,j;
	FLOAT *x,*v,v_circ2,ke=0,pe,r=0,vsq=0;

	x=malloc(3*sizeof(FLOAT));
	v=malloc(3*sizeof(FLOAT));		

	for (i=0;i<nSmooth;i++) {
	  pn=nnList[i].pPart;
	  if (p == pn) continue;
	  /* Now transform to p's rest frame */
	  for (j=0;j<3;j++) {		
	    v[j]=p->v[j] - pn->v[j];
	    x[j]=p->r[j] - pn->r[j];
	    vsq+=v[j]*v[j];
	    ke+=0.5*vsq;
	    r+=x[j]*x[j];
	  }
	  r=sqrt(r);

#ifdef SLIDING_PATCH
	  if (p->r[0] > pn->r[0] && nnList[i].dx < 0)
	    v[1] += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
	  else if (p->r[0] < pn->r[0] && nnList[i].dx > 0)
	    v[1] -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
#endif
	
	  /* First cut: Is the particle bound? */
	  pe=-p->fMass/r;
	  /* WARNING! Here I am overloading the dtCol field. In this 
	   ** situation I am going to fill it with the binding energy of
	   ** the binary. This is used later in pst/pkdFindTightestBinary.
	   */
	  p->dtCol = FLOAT_MAXVAL;
	  if ((ke+pe) < 0 ) {
	    v_circ2=ke-pe;
	    /* This is quick, but not optimal. We need to be sure we don't 
	       automatically merge any bound particles, as those on highly
	       eccentric orbits may be distrupted at apocenter. Therefore
	       we assume that these particles will only reach this point of
	       the code near pericenter (due to iMinBinaryRung), and we can 
	       safely exclude particles that do not meet the following 
	       criterion. Some fiddling with iMinBinaryRung and dMaxBinaryEcc
	       may be necessary to acheive an appropriate balance of merging.
	    */
	    if (vsq < sqrt((1+smf->dMaxBinaryEcc)/(1-smf->dMaxBinaryEcc))*v_circ2) {
	      p->dtCol = ke+pe;
	      p->iOrderCol = pn->iOrder;
	    }
	  }
	}
	free((void *) x);
	free((void *) v);
}

#endif /* RORY_EXTENSION */

/* ugly macros to assist with aggregate and/or particles-stuck-on-wall collision detection (hard-sphere only) */

#if defined(DEM)
#define NEED_2NDORD(p) (FALSE)
#define USE_QUARTIC(p,pn) (FALSE)
#elif defined(AGGS)
int NEED_2NDORD(const PARTICLE *p);
int USE_QUARTIC(const PARTICLE *p,const PARTICLE *pn);
#ifdef WALLS
#define NEED_2NDORD(p) (IS_AGG(p) || PARTICLE_STUCK_ON_ROTATING_WALL(p,&smf->WP))
#define USE_QUARTIC(p,pn) ((smf->bAggsSolveQuartic && (IS_AGG(p) || IS_AGG(pn))) || (smf->WP.bWallsSolveQuartic && (PARTICLE_STUCK_ON_ROTATING_WALL(p,&smf->WP) || PARTICLE_STUCK_ON_ROTATING_WALL(pn,&smf->WP))))
#else /* WALLS */
#define NEED_2NDORD(p) (IS_AGG(p))
#define USE_QUARTIC(p,pn) (smf->bAggsSolveQuartic && (IS_AGG(p) || IS_AGG(pn)))
#endif /* !WALLS */
#elif defined(WALLS)
#define NEED_2NDORD(p) (PARTICLE_STUCK_ON_ROTATING_WALL(p,&smf->WP))
#define USE_QUARTIC(p,pn) (smf->WP.bWallsSolveQuartic && (PARTICLE_STUCK_ON_ROTATING_WALL(p,&smf->WP) || PARTICLE_STUCK_ON_ROTATING_WALL(pn,&smf->WP)))
#endif

void CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Checks whether particle "p" will collide with any of its "nSmooth"
	 ** nearest neighbors in "nnList" during drift interval smf->dStart to
	 ** smf->dEnd, relative to current time.  If any collisions are found,
	 ** the relative time to the one that will occur first is noted in
	 ** p->dtCol along with the iOrder of the collider in p->iOrderCol.
	 ** Note that it can happen that only one particle of a colliding pair
	 ** will "know" about the collision, since the other may be inactive,
	 ** or, for some reason, does not have the other particle on its
	 ** neighbour list.
	 */

	/*
	** NOTE: ideally, when a future collision is detected, we would like to
	** advance (drift) *all* particles up to the collision time, so that we
	** are always moving forward in time, and there is no need to worry
	** about "ghost" particles (post-collision particles drifted back to
	** the start of the step).  BUT, this would require doing a full smooth
	** over all particles after each collision, instead of just recomputing
	** the circumstances of the particles involved in the collision.  This
	** would be far too expensive when collisions are frequent.
	*/

	PARTICLE *pn;
	FLOAT vx,vy,vz,rdotv,v2,sr,dr2,dt;
	int i,bCheckForOverlap,bOverlap;

#if defined(WALLS) || defined(AGGS)
	FLOAT qx,qy,qz;
#endif

	assert(TYPEQueryACTIVE(p)); /* just to be sure */
	assert(smf->dStart >= 0.0 || smf->iOverlapOption == OverlapBackstep);
	assert(smf->dStart < smf->dEnd || (smf->dStart == 0.0 && smf->dEnd == 0.0)); /* should only be equal when doing IC check */

	/* following flag used to simplify overlap handling logic */
	bCheckForOverlap = (smf->bOverlap || (smf->dStart == 0.0 && smf->iOverlapOption != OverlapIgnore && smf->iOverlapOption != OverlapRepel));

	p->dtCol = DBL_MAX; /* initialize */
	p->iOrderCol = -1;
	p->bTinyStep = 0;
	p->fOverlap = 0.0;
	/*
	** The bOverlap flag is set in msrDoCollisions() when an overlap
	** condition was previously detected requiring adjust position or
	** merge.  This flag allows us to keep looking for these types of
	** overlap conditions without advancing the collision clock.
	*/
	if (smf->dStart == 0.0 && !smf->bOverlap) {
		p->dtPrevCol = 0.0; /* these are set in PutColliderInfo() */
		p->iPrevCol = INT_MAX;
		}

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		/*
		 ** Consider all valid particles, even if inactive.  We can
		 ** skip this neighbor if it's the last particle we collided
		 ** with in this search interval, eliminating the potential
		 ** problem of "ghost" particles colliding due to roundoff
		 ** error (but note for aggs, or particles stuck to walls,
		 ** there could be legitimate multiple collisions between the
		 ** particles in the interval, without intervening collisions
		 ** involving *both* of those particles---in the case of aggs,
		 ** this possibility is ignored [see pkdAggsDoCollision()];
		 ** for stuck particles, the free particle must hit something
		 ** else before hitting the stuck particle again, so the free
		 ** particle will recognize that the repeated collision is
		 ** legitimate, even if the stuck particle does not, and the
		 ** collision will be carried out).  We can also ignore
		 ** neighbors scheduled for deletion (negative iOrder).
		 */
		if (pn == p || pn->iOrder < 0 || pn->iOrder == p->iPrevCol)
			continue;
		assert(pn->iOrder != p->iOrder); /* paranoia check */
#ifdef AGGS
		/*
		 ** Do not consider collisions between particles in the same
		 ** aggregate.
		 */
		if (IS_AGG(p) && IS_AGG(pn) && AGG_IDX(p) == AGG_IDX(pn))
			continue;
#endif
#ifdef WALLS
		/* two stuck particles cannot collide (even if on different walls...) */
		if (PARTICLE_STUCK(p) && PARTICLE_STUCK(pn))
			continue;
#endif
		vx = p->v[0] - pn->v[0];
		vy = p->v[1] - pn->v[1];
		vz = p->v[2] - pn->v[2];
#if defined(WALLS) || defined(AGGS)
		/* centripetal terms (vx,vy,vz above contain omega x r term) */
		qx = qy = qz = 0.0;
		if (NEED_2NDORD(p)) {
			qx = p->omega_v[0];
			qy = p->omega_v[1];
			qz = p->omega_v[2];
			}
		if (NEED_2NDORD(pn)) {
			qx -= pn->omega_v[0];
			qy -= pn->omega_v[1];
			qz -= pn->omega_v[2];
			}
#endif /* WALLS || AGGS */
#ifdef SLIDING_PATCH
		if (smf->PP.bPatch) {
			if (p->r[0] > pn->r[0] && nnList[i].dx < 0.0)
				vy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
			else if (p->r[0] < pn->r[0] && nnList[i].dx > 0.0)
				vy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
			}
#endif
		/* for aggs or stuck particles on rotating walls, "rdotv" is really r dot u */
		rdotv = nnList[i].dx*vx + nnList[i].dy*vy + nnList[i].dz*vz;
#if defined(WALLS) || defined(AGGS)
		if (!USE_QUARTIC(p,pn)) /* skip rdotv test if using 4th-order solver */
#endif /* WALLS || AGGS */
		if ((!bCheckForOverlap || !smf->bStrictOverlap) && rdotv >= 0.0)
			continue; /* skip if particles not approaching */
		v2 = vx*vx + vy*vy + vz*vz;
#if defined(WALLS) || defined(AGGS)
		v2 += nnList[i].dx*qx + nnList[i].dy*qy + nnList[i].dz*qz;
#endif /* WALLS || AGGS */
		sr = RADIUS(p) + RADIUS(pn);
		//JD: Changed dr2 to not depend upon sr, the inflation was breaking it
		dr2 = nnList[i].fDist2 - sr*sr;//(RADIUS(p) + RADIUS(pn))*(RADIUS(p) + RADIUS(pn));
		bOverlap = (dr2 <= 0.0); /* zero means just touching */
		dt = DBL_MAX;
#if defined(WALLS) || defined(AGGS)
		if (USE_QUARTIC(p,pn)) {
			double dCoefs[5],dRoots[4]; /* for 4th-order polynomial */
			int j,nRealRoots;

			dCoefs[4] = 0.25*(qx*qx + qy*qy + qz*qz); /* t^4 coefficient */
			dCoefs[3] = vx*qx + vy*qy + vz*qz; /* t^3 coefficient */
			dCoefs[2] = v2; /* t^2 coefficient */
			dCoefs[1] = 2.0*rdotv; /* t coefficient */
			dCoefs[0] = dr2; /* constant */
			if (dCoefs[4] == 0.0) { /* can happen if neither object has rotation (for example) */
				assert(dCoefs[3] == 0.0); /* i.e. reduced to quadratic equation */
				if (dCoefs[2] == 0.0 || polyQuadSolve(dCoefs[2],dCoefs[1],dCoefs[0],&dRoots[0],&dRoots[1]))
					nRealRoots = 0; /* no relative motion, or no solution to quadratic */
				else
					nRealRoots = 2;
				}
			else
				polyFindRealRoots(4,dCoefs,dRoots,&nRealRoots);
			if (!bOverlap && nRealRoots == 0)
				continue; /* no real roots ==> no collision */
			if (nRealRoots > 0) {
				if (!bOverlap) { /* not overlapping -- choose smallest positive root */
					for (j=0;j<nRealRoots;j++)
						if (dRoots[j] > 0.0) { /* array comes pre-sorted */
							dt = dRoots[j];
							break;
							}
					if (dt == DBL_MAX)
						continue; /* surfaces moving away */
					}
				else { /* overlap! -- choose smallest-magnitude non-positive root */
					for (j=nRealRoots-1;j>=0;j--)
						if (dRoots[j] <= 0.0) { /* array comes pre-sorted */
							dt = dRoots[j];
							break;
							}
					assert(dt != DBL_MAX); /* must find a solution *//*DEBUG! not necessarily: the intersection equation we're solving is only an approximation of the true motion---leave this as an assert for now, but we may need to allow for this possibility, like we do below for the quadratic case */
					}
				}
			}
		else
#endif /* WALLS || AGGS */
		  {
			  double dt1,dt2;

			  if (!bOverlap && rdotv >= 0.0)
				  continue; /* particles not approaching, touching, or overlapping */
			  if (v2 != 0.0) { /* v2 can be negative for quadratic case as it contains an extra signed term */
				  if (polyQuadSolve(v2,2.0*rdotv,dr2,&dt1,&dt2))
					  continue; /* no real solutions ==> no collision */
				  if (!bOverlap) { /* not overlapping -- choose smallest positive root */
					  if (dt1 > 0.0) /* roots guaranteed to have dt1 <= dt2 */
						  dt = dt1;
					  else if (dt2 > 0.0)
						  dt = dt2;
					  else
#if defined(WALLS) || defined(AGGS)
						  continue; /* may be dealing with approximate motion, so we just give up */
#else
						  assert(0); /* means dt <= 0 but rdotv < 0 -- shouldn't happen! */
#endif
					  }
				  else { /* overlap! -- choose smallest-magnitude non-positive root */
					  if (dt2 <= 0.0) /* roots guaranteed to have dt1 <= dt2 */
						  dt = dt2;
					  else if (dt1 <= 0.0)
						  dt = dt1;
					  else
#if defined(WALLS) || defined(AGGS)
						  continue; /* may be dealing with approximate motion, so we just give up */
#else
						  assert(0); /* means dt > 0 but overlapping -- shouldn't happen! */
#endif
					  }
				  assert(dt != DBL_MAX);
				  }
			  }
		/*
		 ** Normally there should be no touching or overlapping particles
		 ** at the start of the step.  But inelastic collapse and other
		 ** numerical problems may make it necessary to relax this...
		 */
		if (bCheckForOverlap && bOverlap) {
			FLOAT fOverlap;
			if (dt == DBL_MAX) { /* no relative motion */
				assert(v2 == 0.0 && rdotv == 0.0);
				dt = 0.0;
				}
			assert(dt <= 0.0);
			assert(nnList[i].fDist2 >= 0.0);
			assert(sr > 0.0);
			fOverlap = 1.0 - sqrt(nnList[i].fDist2)/sr; //sr replaced with (RADIUS(p) + RADIUS(pn))
			if (smf->iOverlapOption == OverlapIsError &&
				(!smf->bAllowSimulColl || dt < 0.0)) {
				fprintf(stderr,"OVERLAP/CONTACT [t=%e]:\n"
						"%i (r=%g,%g,%g,iRung=%i) &\n"
						"%i (rn=%g,%g,%g,iRung=%i)\n"
						"fDist=%g v=%g,%g,%g v_mag=%g\n"
						"rv=%g sr=%g overlap=%g%% dt=%g\n",
						smf->dTime,
						p->iOrder,p->r[0],p->r[1],p->r[2],p->iRung,
						pn->iOrder,pn->r[0],pn->r[1],pn->r[2],pn->iRung,
						sqrt(nnList[i].fDist2),vx,vy,vz,sqrt(vx*vx+vy*vy+vz*vz),
						rdotv,sr,100.0*fOverlap,dt);
				assert(0); /* particle not allowed to touch or overlap initially */
				}
			else {
#if (INTERNAL_WARNINGS != 0)
				if (p->iOrder < pn->iOrder &&
					(smf->iOverlapOption != OverlapAdjPos || fOverlap >= smf->dAdjPosLimit)) {
					static int nWarn = 1;
					if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
						fprintf(stderr,"OVERLAP/CONTACT WARNING #%i (pid=%i) [t=%e]: %i & %i (%g%%, dt=%g)\n",
								nWarn,smf->pkd->idSelf,smf->dTime,p->iOrder,pn->iOrder,100.0*fOverlap,dt);
					++nWarn;
					}
#endif /* INTERNAL_WARNINGS */
				/* in some cases we just ignore overlap */
				switch (smf->iOverlapOption) {
				case OverlapIsError:
					continue; /* ignore; means dt = 0 but simultaneous collisions allowed */
				case OverlapBackstep:
					if (fabs(dt) > smf->dBackstepLimit)
						continue; /* ignore if backstep too large */
					if (dt == 0.0) {
#if (INTERNAL_WARNINGS != 0)
						static int nWarn = 1;
						if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
							fprintf(stderr,"ZERO BACKSTEP WARNING #%i (pid=%i) [t=%e]: %i & %i\n",
									nWarn,smf->pkd->idSelf,smf->dTime,p->iOrder,pn->iOrder);
						++nWarn;
#endif /* INTERNAL_WARNINGS */
						continue;
						}
					break;
				case OverlapRepel:
					assert(0); /* shouldn't be here --- overlap dealt with in gravity calculation */
				case OverlapAdjPos:
					if (fOverlap < smf->dAdjPosLimit)
						continue; /* ignore if adjustment too small */
				case OverlapMerge:
					break; /* dealt with in pkdDoCollision() or pkdAggsDoCollision() */
				default:
					assert(0);
					}
				/*
				** Take largest overlap (for this particle).  But for
				** backstep, take largest negative step, which may not
				** be the same collision!
				*/
				if ((smf->iOverlapOption == OverlapBackstep && dt < p->dtCol) ||
					(smf->iOverlapOption != OverlapBackstep && fOverlap > p->fOverlap)) {
					if (smf->iOverlapOption == OverlapBackstep)
						p->dtCol = dt;
					else
						/* NOTE: fOverlap is a relative measure ---
						   might be better to choose largest absolute
						   distance correction instead of largest
						   fOverlap (e.g. small particles inside big
						   moon)... */
						p->dtCol = -fOverlap; /* ensures largest overlap handled first */
					p->iOrderCol = pn->iOrder;
					p->fOverlap = fOverlap;
					}
				continue;
				}
			} /* if overlap */
		/* finally, decide if this collision should be stored */
		if (dt > smf->dStart && dt <= smf->dEnd) {
			if (dt >= p->dtCol) /* skip if not an earlier collision */
				continue; /* note: simultaneous collisions handled in msrDoCollisions() */
			p->dtCol = dt;
			p->iOrderCol = pn->iOrder;
			/* rdotv slightly different for aggs in following---ok?... */
			if (smf->dCollapseLimit > 0.0)
				_CheckForCollapse(p,dt,rdotv,nnList[i].fDist2,smf);
#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
			/*
			 ** At the start of the top step we need to know if any
			 ** *planetesimals* (i.e., not rubble pieces) are predicted
			 ** to collide during the top step so that we can demote
			 ** them to the lowest timestep rungs (since the rubble
			 ** pieces need to be treated much more carefully, and the
			 ** collision will generally occur sometime in the middle
			 ** of the step).  Note that we only check at the beginning
			 ** of the step (i.e. when dStart is exactly zero), since
			 ** there are no resmoothing circumstances that can lead
			 ** two planetesimals to collide that have *not* already
			 ** undergone a collision themselves during the step.
			 */
			if (smf->dStart == 0.0 && p->iColor == PLANETESIMAL && pn->iColor == PLANETESIMAL){
				p->bMayCollide = 1;
				//puts("HERE1");
				}
			/* note: flag is reset with call to pkdRubbleResetColFlag() */
#endif /* RUBBLE_ZML || COLLMOD_ZML */
			}
		}

#ifdef WALLS
	/* now check for any wall collisions, if applicable */
	{
		WALL_PARAMS *WP = &smf->WP;
		WALL *w;
		double dtCol[WALL_MAX_NUM_FACES],dtColPrev=DBL_MAX;
		int j,iColPrev=-1;

		if (PARTICLE_STUCK(p))
			return; /* already stuck on wall */
		for (i=0;i<WP->nWalls;i++) {
			w = &WP->pWalls[i];
			wallsGetTimesToIntersect(WP,w,p,smf->dStart,smf->dEnd,dtCol,&bOverlap);
			for (j=0;j<nWallFaces[w->wd.iType];j++) {
				dt = dtCol[j];
				/*
				** Accept this collision if it occurs during the
				** search interval (including endpoints) and: 1) it
				** occurs at or before any previously determined
				** collision time for this particle during this
				** search; or 2) it's a backstep and the previously
				** determined collision time wasn't; or 3) it's a
				** backstep, as was the previously determined
				** collision time, but it's a smaller backstep.
				*/
				if ((dt >= smf->dStart || (bOverlap && fabs(dt) < smf->dBackstepLimit)) && dt <= smf->dEnd &&
					((dt >= smf->dStart && p->dtCol >= smf->dStart && dt <= p->dtCol) ||
					 (dt <= smf->dStart && p->dtCol >= smf->dStart) ||
					 (dt <= smf->dStart && p->dtCol <= smf->dStart && fabs(dt - smf->dStart) <= fabs(dt - p->dtCol)))) {
					dtColPrev = p->dtCol; /* store previous values in case of simultaneous collision (see below) */
					iColPrev = p->iOrderCol;
					p->dtCol = dt;
					p->iOrderCol = -1 - i; /* wall index encoding */
					if (smf->dCollapseLimit > 0.0 && !bOverlap && dt > 0.0 /*DEBUG! bOverlap should have been set!*/) {
						rdotv = -dt*vectorMagSq(p->v); /* minus sign because objects are approaching */
						dr2 = -dt*rdotv; /* minus sign to get magnitude */
						_CheckForCollapse(p,dt,rdotv,dr2,smf);
					}
				}
			if (bOverlap) {
				/* moved overlap check after dtCol loop to allow for simultaneous collisions */
				assert(smf->iOverlapOption == OverlapIsError || smf->iOverlapOption == OverlapBackstep); /* only these 2 overlap options currently supported */
				if (smf->iOverlapOption == OverlapIsError &&
					(!smf->bAllowSimulColl || p->dtCol < smf->dStart)) {
					fprintf(stderr,"WALL OVERLAP/CONTACT [t=%e]: %i (r=%g,%g,%g,iRung=%i) & wall %i, dt=%g\n",
							smf->dTime,p->iOrder,p->r[0],p->r[1],p->r[2],p->iRung,i,p->dtCol - smf->dStart);
					assert(0);
					}
				else {
#if (INTERNAL_WARNINGS != 0)
					static int nWarn = 1;
					if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
						fprintf(stderr,"WALL OVERLAP/CONTACT WARNING #%i (pid=%i) [t=%e]: %i & wall %i, dt=%g\n",
								nWarn,smf->pkd->idSelf,smf->dTime,p->iOrder,i,p->dtCol - smf->dStart);
					++nWarn;
#endif /* INTERNAL_WARNINGS */
					if (smf->iOverlapOption == OverlapIsError) {
						/* revert to previous soonest collision, if any */
						p->dtCol = dtColPrev;
						p->iOrderCol = iColPrev;
						assert(p->dtCol > smf->dStart); /* otherwise, too many simultaneous collisions for this particle */
						}
					}
				} /* if overlap */
			}
		}
	}
#endif /* WALLS */
}

#undef NEED_2NDORD

#ifdef SPRINGS

void AssignSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** This routine should be called once (for each processor) at the
	** very beginning of the simulation in order to determine which
	** particles are within a linkage length
	** (smf->FO.SP.dLinkageLength) of each other initially and will
	** therefore be subject to spring forces.  This quantity is
	** measured in effective particle radii (see NOTE below).  If
	** smf->FO.SP.dZeroStrainLength is zero, the initial strain
	** between all linked pairs is set to zero, regardless of their
	** separations (meaning, their zero strain lengths are set to
	** their component particle separations in units of the effective
	** radius).  Otherwise, the initial strain is equal to the excess
	** separation relative to the zero strain length (again measured
	** in effective particle radii).  If a different initial strain
	** distribution is desired, run pkdgrav for only 1 step, convert
	** the resulting springs file to human-readable text, and edit the
	** file by hand.  Conversions can be done using the spr2txt and
	** txt2spr utilities in the ss_core package.
	**
	** NOTE: for cases where there is a size distribution (particles
	** with different radii), we choose to use an effective radius
	** equal to the average: Re = (R1 + R2)/2.  It's not clear what
	** the correct choice should be.  Other possibilities include the
	** geometric mean Re = sqrt(R1*R2) or the reduced radius Re =
	** (R1*R2)/(R1 + R2), but these latter choices have possibly
	** undesirable limiting behaviors when R2/R1 -> 0 (specifically,
	** Re --> 0), whereas the first choice does not (Re --> R1/2).
	**
	** IMPORTANT: nSmooth must be at least as large as the maximum
	** number of particles that can pack within the strain length
	** times the effective radius, otherwise springs may not be
	** assigned symmetrically, resulting in loss of momentum
	** conservation.  This will be especially tricky to enforce if
	** particles have a range of sizes...
	*/

	const PARTICLE *pn;
	const SPRING_PARAMS *SP = &smf->FO.SP;
	SPRING *ps;
	double dEffRad,d,dStrainLimit;
	int i;

	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);
	mdlassert(smf->pkd->mdl,nSmooth <= MAX_NUM_SPRINGS_PER_PARTICLE);

	for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
#ifdef YOUNGS_INCREASE
		if (1)
#elif (REFORM_SPRINGS == 0 || REFORM_SPRINGS == 1 || REFORM_SPRINGS == 3) /* YOUNGS_INCREASE */
	    if (i >= nSmooth)
#elif (REFORM_SPRINGS == 2) /* YOUNGS_INCREASE */
		if (i >= nSmooth || p->iColor < 0) /* can safely skip wall particles */
#endif /* YOUNGS_INCREASE */
		{
			ps = &p->springs[i]; /* array order doesn't matter */
			ps->iOrder = -1; /* unused neighbor slot */
			ps->fZeroStrainLength = 0.;
			ps->fYoungsModulus = 0.;
			ps->fStressLimit = 0.;
			continue;
			}
		ps = &p->springs[i]; /* array order doesn't matter */
		pn = nnList[i].pPart;
		assert(pn->iOrder >= 0);

#if (REFORM_SPRINGS == 0)
		if (p->iOrder >= pn->iOrder)
#elif (REFORM_SPRINGS == 1) /* REFORM_SPRINGS */
		if (p->iOrder >= pn->iOrder || (p->iColor < 0 && pn->iColor < 0)) /* save disk space by not forming springs between two wall particles */
#elif (REFORM_SPRINGS == 2) /* REFORM_SPRINGS */
		if (p->iOrder >= pn->iOrder || pn->iColor < 0) /* can safely skip neighbors that are wall particles */
#elif (REFORM_SPRINGS == 3)  /* REFORM_SPRINGS */
		if (p->iOrder >= pn->iOrder || (p->iColor >= 0 && pn->iColor >= 0) || (p->iColor < 0 && pn->iColor < 0)) /* allow spring to form only between one wall particle and one non-wall particle */
#endif /* REFORM_SPRINGS */
		{
			ps->iOrder = -1; /* unused neighbor slot */
			ps->fZeroStrainLength = 0.;
			ps->fYoungsModulus = 0.;
			ps->fStressLimit = 0.;
			continue;
			}

		/* PARTICLE PAIR QUALIFIES: FORM SPRING BETWEEN PARTICLE AND THIS NEIGHBOR IF INITIAL SEPARATION WARRANTS */
		ps->iOrder = pn->iOrder;
		ps->fZeroStrainLength = 0.; /* initialize: zero means no springs */
		ps->fYoungsModulus = 0.;
		ps->fStressLimit = 0.;
		dStrainLimit = 0.;
		assert(nnList[i].fDist2 > 0.);
		d = sqrt(nnList[i].fDist2); /* could optimize by working with square instead */
		dEffRad = 0.5*(RADIUS(p) + RADIUS(pn)); /* effective radius = average radius */
#if (REFORM_SPRINGS == 0)
		/* check to see if neighbor is close enough such that spring should be attached */
		if (d > SP->dLinkageLength*dEffRad || SP->dMeanYoungsModulus == 0.) {
			ps->iOrder = -1;
			continue;
			}

		/* assign equilibrium separation */
		if (SP->dZeroStrainLength == 0.)
			ps->fZeroStrainLength = d/dEffRad;
		else
			ps->fZeroStrainLength = SP->dZeroStrainLength;

		/*
		** Assign Young's modulus using specified standard deviation
		** around the mean value and throw out values zero or less and
		** 2x MeanY and greater.
		*/
		while (ps->fYoungsModulus <= 0. || ps->fYoungsModulus >= (SP->dMeanYoungsModulus + SP->dMeanYoungsModulus)) /* avoid negatives & maintain mean by making 0 < Y_range < 2*MeanY */
			ps->fYoungsModulus = SP->dMeanYoungsModulus + SP->dYoungsStdDev*(double)randGaussian();

		/* Glass Beads Hack */
#ifdef GLASS_BEADS/*DEBUG!!!*/
		static const double dConvert = (1.49597870e11/1.9891e30*(365.25*24*3600)*(365.25*24*3600)/(2.0*M_PI)/(2.0*M_PI)); /* Pascals --> pkdgrav units (multiply) [1.897217e-06] */
		double dCrossSection = M_PI*dEffRad*dEffRad*smf->FO.SP.dPorosityFactor;
		double fZ = ps->fZeroStrainLength;
		ps->fZeroStrainLength /= 1. + (smf->FO.DP.dKn*(dEffRad + dEffRad - d)/(dCrossSection*ps->fYoungsModulus*dConvert));
		dStrainLimit = 2./ps->fZeroStrainLength - 1.;
		ps->fStressLimit = dStrainLimit*ps->fYoungsModulus;
		printf("%d %d d=%f%% R_eff=%e old_FZ=%f fZ=%f Y=%e StrainLim=%e fStressLim=%e --> %f N\n",
			   p->iOrder,ps->iOrder,100. - 50.*d/dEffRad,dEffRad,fZ,ps->fZeroStrainLength,ps->fYoungsModulus,dStrainLimit,ps->fStressLimit,ps->fStressLimit*dCrossSection*1.49597870e11*1.49597870e11);
		continue;
#endif

		/*
		** Assign strain limit using specified standard deviation
		** around the mean value and throw out values zero or less,
		** 2x(MeanY) and greater, greater than max allowed limit, and
		** less than values 2x mean value the distance from the max
		** allowed strain limit.
		*/
		while (
			   /* avoid negatives & maintain mean by making 0 < Strain_range < 2*MeanStrainLimit: */
			   (SP->dMeanStrainLimit <= 0.5*SP->dMaxStrainLimit && (dStrainLimit <= 0. || dStrainLimit >= (SP->dMeanStrainLimit + SP->dMeanStrainLimit)))
			   /* avoid negatives & maintain mean by making MaxStrainLimit-2*(MaxStrainLimit-MeanStrainLimit) < Strain_range < 1: */
			   || (SP->dMeanStrainLimit > 0.5*SP->dMaxStrainLimit && (dStrainLimit > SP->dMaxStrainLimit || dStrainLimit < (SP->dMeanStrainLimit + SP->dMeanStrainLimit - SP->dMaxStrainLimit))))
			dStrainLimit = SP->dMeanStrainLimit + SP->dStrainStdDev*(double)randGaussian();

		/* assign stress limit */
		ps->fStressLimit = dStrainLimit*ps->fYoungsModulus;

#else /* REFORM_SPRINGS */

		/* assign equilibrium separation */
		ps->fZeroStrainLength = ZERO_STRAIN_LENGTH;
		if (d > (ps->fZeroStrainLength*(1. + SP->dMeanStrainLimit)*dEffRad)) {
			ps->iOrder = -1;
			continue;
			}

		/* assign Young's modulus */
		while (ps->fYoungsModulus <= 0. || ps->fYoungsModulus >= (SP->dMeanYoungsModulus + SP->dMeanYoungsModulus))
			ps->fYoungsModulus = SP->dMeanYoungsModulus + SP->dYoungsStdDev*(double)randGaussian();

		/* assign stress limit */
		ps->fStressLimit = SP->dMeanStrainLimit*ps->fYoungsModulus; /* DEBUG: for now, no variation is allowed in strain limit using REFORM_SPRINGS...may want to change this (1/2) */
#endif /* REFORM_SPRINGS */
		/*
		** if (p->iColor == 4 || pn->iColor == 4) ps->fStressLimit = 10.*dStrainLimit*ps->fYoungsModulus; HACK! (for stress tests : less likely to fission near wall)
		*/
		}
	}

void ColorSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Colors particles based upon number of attached springs:
	**
	**   yellow: 2 springs attached
	**   orange: 1 spring attached
	**   red: unlinked particle
	*/

	springsAssignColor(p,nSmooth,nnList,ALLOW_COLOR_CHANGE(p));
}

void MoveSpecialSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Uses particle color to identify 'special particles' that can be
	** either repositioned or subject to external forces.
	**
	** Types:
	**
	** I.  Apply forces to special particles: (keep track of force applied to special particles, calculate area --> pressure)
	**     A.  pullplane - force outermost x-y(y-z,z-x) planes in x(y,z)-dimension away from x(y,z)=0 plane (tensile stress)
	**     B.  shearplane - force outermost x-y planes into z-dimension (opposite directions)
	**     C.  twistplane - tangential force applied to outermost x-y planes around x-axis
	**     D.  pullradial - applied force directed out radially from origin
	** II. Move special particles: (measure restoring forces exerted on special particles, calculate area --> pressure)
	**     A.  moveplane - move outermost x-y(y-z,z-x) planes in x(y,z)-dimension away from x(y,z)=0 plane (tensile stress)
        */

#if (SPRINGS_STRESS_TEST)
	/* POTENTIALLY USEFUL HACK: the following "if" bracket --> halts integration if any particle is far away from origin */
	if ((sign(p->r[0])*p->r[0] > 5e-7) || (sign(p->r[1])*p->r[1] > 5e-7) || (sign(p->r[2])*p->r[2] > 5e-7)) {
		printf("particle %d at (%e,%e,%e) is far from origin... integration halting here.\n",p->iOrder,p->r[0],p->r[1],p->r[2]);
		exit(1);
	}
#endif /* SPRINGS_STRESS_TEST */

	/* the following check should be applied if gravity is off (i.e., this check expects that no kicks have been performed on the particle this step):
	if (p->a[0] != 0. || p->a[1] != 0. || p->a[2] != 0.) {
	        printf("\"moveplane\" assumes that the particle's acceration is zero before adding the spring force.  The particle's acceleration is nonzero... exiting.");
		exit(1); took out for gravity
	}
	*/

	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);
	mdlassert(smf->pkd->mdl,nSmooth <= MAX_NUM_SPRINGS_PER_PARTICLE);

#ifdef SPRINGS_STRESS_MOVE
	springsMoveSpecials(p,nSmooth,nnList,smf);
#elif (SPRINGS_STRESS_TEST == SPRINGS_STRESS_FORCE_TENSILE)
	springsForceSpecialsTensile(p);
#elif (SPRINGS_STRESS_TEST == SPRINGS_STRESS_FORCE_SHEAR)
	springsForceSpecialsShear(p);
#elif (SPRINGS_STRESS_TEST == SPRINGS_STRESS_FORCE_TWIST)
	springsForceSpecialsTwist(p);
#elif (SPRINGS_STRESS_TEST == SPRINGS_STRESS_FORCE_RADIAL)
	springsForceSpecialsRadial(p);
#endif /* SPRINGS_STRESS_TEST */
}

void initDoSpringsParticle(void *p)
{

	/* initialize each particle */

	int k;

	for (k=0;k<3;k++)
		((PARTICLE *)p)->dDeltaAccel[k] = 0.;
	}

void initDoSprings(void *p)
{

	/* initialize each cached particle (may not be necessary?...) */

	int k;

	for (k=0;k<3;k++)
		((PARTICLE *)p)->dDeltaAccel[k] = 0.;
	}

void combDoSprings(void *p1,void *p2)
{

	/* update each particle with accumulated information contained in the cache */

	int k;

	for (k=0;k<3;k++)
		((PARTICLE *)p1)->dDeltaAccel[k] += ((PARTICLE *)p2)->dDeltaAccel[k];
	}


void postDoSprings(PARTICLE *p,SMF *smf)
{

	/* save new accelerations for each particle */

	/* (note: we use dDeltaAccel because we may have multiple processors updating the same particles) */

	int k;

	for(k=0;k<3;k++)
		p->a[k] += p->dDeltaAccel[k];
	}

#if (REFORM_SPRINGS == 0)
void DoSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Checks neighbors to see whether a spring exists between them
	** (see AssignSprings()), and if so, whether a restoring force
	** should be applied, or if the link should be broken.  If a
	** particle no longer has any springs attached to it, its color is
	** changed to 2 (red) for visualization purposes.
	**
	** In the current formulation, a stress limit is equivalent to a
	** maximum strain, and the strain depends on the particle
	** separation.  In particular, the maximum particle separation
	** (beyond which the restoring force is zero) is D(1 + L/Y), where
	** D is the separation at zero strain, L is the stress limit, and
	** Y is the Young's modulus.  E.g., if L = Y/2, then the maximum
	** separation is 3D.  In our model, D is measured in "effective
	** particle radii" (see NOTE in AssignSprings() above).
	**
	** IMPORTANT: it is assumed that particles in springs always see
	** each other (i.e. they are symmetric) -- see note about symmetry
	** in AssignSprings() above.
	*/

	static const double dConvert = (1.49597870e11/1.9891e30*(365.25*24*3600)*(365.25*24*3600)/(2.0*M_PI)/(2.0*M_PI)); /* Pascals --> pkdgrav units (multiply) [1.897217e-06] */
	static long int iCounter = 0; ++iCounter;

	/* DEBUG!!! HACK!!! ANGMOM HACK (START) */
	/*
	static double Lx = 0., Lx0 = 0.; Lx += p->v[0];
	static double Ly = 0., Ly0 = 0.; Ly += p->v[1];
	static double Lz = 0., Lz0 = 0.; Lz += p->v[2];
	static double Lx_pred = 0., Lx_pred0 = 0.; Lx_pred += p->vPred[0];
	static double Ly_pred = 0., Ly_pred0 = 0.; Ly_pred += p->vPred[1];
	static double Lz_pred = 0., Lz_pred0 = 0.; Lz_pred += p->vPred[2];
	if (iCounter == 686) {Lx0 = Lx ; Ly0 = Ly ; Lz0 = Lz ; Lx_pred0 = Lx_pred ; Ly_pred0 = Ly_pred ; Lz_pred0 = Lz_pred;}
	if (iCounter % 686 == 0) printf("%d %e %e %e %e %e %e\n",p->iOrder,(Lx-Lx0)/Lx0,(Ly-Ly0)/Ly0,(Lz-Lz0)/Lz0,(Lx_pred-Lx_pred0)/Lx_pred0,(Ly_pred-Ly_pred0)/Ly_pred0,(Lz_pred-Lz_pred0/Lz_pred0));
	if (iCounter % 686 == 0) Lx = Ly = Lz = Lx_pred = Ly_pred = Lz_pred = 0.;
	*/
	/* DEBUG!!! HACK!!! ANGMOM HACK (END) */

#if (SPRINGS_STRESS_TEST)
	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	static const int iParticles = 1176; /* # of particles */
	static const int iTugFreq = iParticles * 2500; /* increase acceleration (or distance: moveplane) of special particle by dTugStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
	        start_of_step = 1;
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

#endif /* SPRINGS_STRESS_TEST */

	PARTICLE *pn;
	SPRING *ps;
	double d,dEffRad,dStrain,dStress,dCrossSection,vx,vy,vz;
	double dReducedMass,dConvert_CrossSection,dRestoringQuantity,dDampingQuantity;
	double dRestoringForce_Particle,dRestoringForce_Neighbor,dDampingForce_Particle,dDampingForce_Neighbor,dDampingCrit;
	int i,j,bFoundNbr;

#if (SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE)
	static double dPressure_right = 0.;
	static double dPressure_left = 0.;
	if (start_of_step)
		dPressure_left = dPressure_right = 0.;
#endif
	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);
	/*
	printf("%i %li %i\n",p->iOrder,smf->pkd->mdl,nSmooth);
	mdlassert(smf->pkd->mdl,nSmooth <= MAX_NUM_SPRINGS_PER_PARTICLE);
	*/

#ifdef WALLS
	/* following time-saver will not work for liberated particles that
	   are stuck to a wall, because stuck particles have negative
	   colors -- so we'll just enter the spring search loop and find
	   nothing: inefficient, but ok */
	/* note: may also skip if p->iColor < 0 [srs] */
#endif
	if (p->iColor == 2) /* also had needed to add: && (iCounter % iParticles) in this if-statement for stress-test purposes - not sure why - hoping it's superfluous in current state.
			    ** planning to leave this comment in until confident...
			    */
		return; /* save time by ignoring previously liberated particles - allow last particle in step to proceed */

	/* loop over this particle's springs */
	for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
	        ps = &p->springs[i];
		if (ps->iOrder < 0 || ps->fZeroStrainLength == 0.0)
		        continue;
		bFoundNbr = 0;
		/* see if other spring particle is one of the current neighbors */
		for (j=0;j<nSmooth;j++) {
		        pn = nnList[j].pPart;
			if (ps->iOrder != pn->iOrder)
			        continue;
			if (pn->iColor == 2)
			        continue;
			bFoundNbr = 1;
			/* perform on neighbor's call if neighbor has lower iOrder number or skip if neighbor is particle */
			if (p->iOrder >= pn->iOrder)
			        continue;
			/* compare stress to stress limit */
			assert(nnList[j].fDist2 > 0.0);
			d = sqrt(nnList[j].fDist2); /* could optimize by working with square instead */
			dEffRad = 0.5*(RADIUS(p) + RADIUS(pn)); /* effective radius = average radius */
			dStrain = (d/dEffRad - ps->fZeroStrainLength)/ps->fZeroStrainLength;
			/* no repulsion */
			if (dStrain < 0.0)
				dStrain = 0.0;
			dStress = ps->fYoungsModulus*dStrain;
			if (dStress > ps->fStressLimit) {
				ps->fZeroStrainLength = 0.0; /* this bond is broken */
				continue;
			}
			vx = p->vPred[0] - pn->vPred[0]; /* DEBUG!!! HACK!!! vPred --> v !! */
			vy = p->vPred[1] - pn->vPred[1]; /* DEBUG!!! HACK!!! vPred --> v !! */
			vz = p->vPred[2] - pn->vPred[2]; /* DEBUG!!! HACK!!! vPred --> v !! */
					
			/* apply restoring force in direction toward neighbor */
			dCrossSection = M_PI*dEffRad*dEffRad*smf->FO.SP.dPorosityFactor;
			dReducedMass = (p->fMass*pn->fMass)/(p->fMass+pn->fMass);

			/* define intermediate terms to get units of acceleration per length in dRestoringForce
			   and units of acceleration per length per velocity in dDamping */
			dConvert_CrossSection = dConvert*dCrossSection;
			dRestoringQuantity = dConvert_CrossSection*dStress;

			/* define critical damping value */
			dDampingCrit = sqrt(4.*dReducedMass*dConvert_CrossSection*ps->fYoungsModulus/(ps->fZeroStrainLength*dEffRad));
			dDampingQuantity = smf->FO.SP.dDamp*dDampingCrit*(vx*nnList[j].dx+vy*nnList[j].dy+vz*nnList[j].dz);

			/*** this particle ***/
			dRestoringForce_Particle = dRestoringQuantity/(p->fMass*d);
			dDampingForce_Particle = dDampingQuantity/(nnList[j].fDist2*p->fMass);

			p->dDeltaAccel[0] -= (dRestoringForce_Particle+dDampingForce_Particle)*nnList[j].dx; /* sign: (dx,dy,dz) = r - rn */
			p->dDeltaAccel[1] -= (dRestoringForce_Particle+dDampingForce_Particle)*nnList[j].dy;
			p->dDeltaAccel[2] -= (dRestoringForce_Particle+dDampingForce_Particle)*nnList[j].dz;

			/*** neighbor ***/
			dRestoringForce_Neighbor = dRestoringQuantity/(pn->fMass*d);
			dDampingForce_Neighbor = dDampingQuantity/(nnList[j].fDist2*pn->fMass);

			pn->dDeltaAccel[0] += (dRestoringForce_Neighbor+dDampingForce_Neighbor)*nnList[j].dx;
			pn->dDeltaAccel[1] += (dRestoringForce_Neighbor+dDampingForce_Neighbor)*nnList[j].dy;
			pn->dDeltaAccel[2] += (dRestoringForce_Neighbor+dDampingForce_Neighbor)*nnList[j].dz;

#if (SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE) /* save values for pressure approximation */
			if (p->iColor == 4 && pn->iColor != 4) { /* this particle is special, neighbor is not */
				if (nnList[j].dx < 0.) dPressure_left += (dRestoringQuantity/d + dDampingQuantity/nnList[j].fDist2) * nnList[j].dx;
				if (nnList[j].dx > 0.) dPressure_right += (dRestoringQuantity/d + dDampingQuantity/nnList[j].fDist2) * nnList[j].dx;
				}
			if (p->iColor != 4 && pn->iColor == 4) { /* neighbor is special, this particle is not */
				if (nnList[j].dx > 0.) dPressure_left -= (dRestoringQuantity/d + dDampingQuantity/nnList[j].fDist2) * nnList[j].dx;
				if (nnList[j].dx < 0.) dPressure_right -= (dRestoringQuantity/d + dDampingQuantity/nnList[j].fDist2) * nnList[j].dx;
				}
#endif /* SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE */
			break;
			}
		if (!bFoundNbr) {
			ps->fZeroStrainLength = 0.0;
#if (INTERNAL_WARNINGS != 0)
			static int nWarn = 1;
			if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
				(void) fprintf(stderr,"SPRINGS WARNING #%i (pid=%i) [T=%g]: %i & %i not neighbors, spring removed\n(momentum may not be conserved!!)\n",
							   nWarn,smf->pkd->idSelf,smf->dTime,p->iOrder,ps->iOrder);
			++nWarn;
#endif /* INTERNAL_WARNINGS */
			}
        }

#if (SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE)

	/* springsCalcPressure(p,iParticles,fZeroStrainLength,dCubeStrain,start_of_step,dConvert); */

#endif /* SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE */

	/*
	** take out statement below once confident extra condition isn't needed above
	
	if (p->iColor == 2 && !(iCounter % iParticles))
	return; */ /* skip next line for last particle of step, which "should have" been skipped above before looping over springs */
}

#elif (REFORM_SPRINGS > 0) /* REFORM_SPRINGS */
void DoSprings(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Checks neighbors to see whether a spring exists between them
	** (see AssignSprings()), and if so, whether a restoring force
	** should be applied, or if the link should be broken.  If a
	** particle no longer has any springs attached to it, its color is
	** changed to 2 (red) for visualization purposes.
	**
	** In the current formulation, a stress limit is equivalent to a
	** maximum strain, and the strain depends on the particle
	** separation.  In particular, the maximum particle separation
	** (beyond which the restoring force is zero) is D(1 + L/Y), where
	** D is the separation at zero strain, L is the stress limit, and
	** Y is the Young's modulus.  E.g., if L = Y/2, then the maximum
	** separation is 3D.  In our model, D is measured in "effective
	** particle radii" (see NOTE in AssignSprings() above).
	**
	** IMPORTANT: it is assumed that particles in springs always see
	** each other (i.e. they are symmetric) -- see note about symmetry
	** in AssignSprings() above.
	*/

	static const double dConvert = (1.49597870e11/1.9891e30*(365.25*24*3600)*(365.25*24*3600)/(2.0*M_PI)/(2.0*M_PI)); /* Pascals --> pkdgrav units (multiply) [1.897217e-06] */

	const SPRING_PARAMS SP = &smf->FO.SP;

	PARTICLE *pn;
	SPRING *ps;
	double d,dEffRad,dStrain,dStress,dCrossSection,vx,vy,vz;
	double dReducedMass,dConvert_CrossSection,dRestoringQuantity,dDampingQuantity;
	double dRestoringAccel_Particle,dRestoringAccel_Neighbor,dDampingAccel_Particle,dDampingAccel_Neighbor,dDampingCrit,dZeroStrainLength;
	int i,j,bFoundNbr,bApplyForce,iSpring;
	static unsigned long int liCounter; liCounter++;

	/* DEBUG!!! HACK!!! ANGMOM HACK (START) */
	/*
	static double Lx = 0., Lx0 = 0.; Lx += p->v[0];
	static double Ly = 0., Ly0 = 0.; Ly += p->v[1];
	static double Lz = 0., Lz0 = 0.; Lz += p->v[2];
	static double Lx_pred = 0., Lx_pred0 = 0.; Lx_pred += p->vPred[0];
	static double Ly_pred = 0., Ly_pred0 = 0.; Ly_pred += p->vPred[1];
	static double Lz_pred = 0., Lz_pred0 = 0.; Lz_pred += p->vPred[2];
	if (liCounter == 686) {Lx0 = Lx ; Ly0 = Ly ; Lz0 = Lz ; Lx_pred0 = Lx_pred ; Ly_pred0 = Ly_pred ; Lz_pred0 = Lz_pred;}
	if (liCounter % 68600 == 0) printf("%d %e %e %e %e %e %e\n",p->iOrder,(Lx-Lx0)/Lx0,(Ly-Ly0)/Ly0,(Lz-Lz0)/Lz0,(Lx_pred-Lx_pred0)/Lx_pred0,(Ly_pred-Ly_pred0)/Ly_pred0,(Lz_pred-Lz_pred0/Lz_pred0));
	if (liCounter % 686 == 0) Lx = Ly = Lz = Lx_pred = Ly_pred = Lz_pred = 0.;
	*/
	/* DEBUG!!! HACK!!! ANGMOM HACK (END) */

	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);

#ifdef YOUNGS_INCREASE
	if (liCounter == 1)
	        (void) printf("SPRINGS: YOUNGS_INCREASE defined, YOUNGS_INCREASE_PER_STEP = %e\n",YOUNGS_INCREASE_PER_STEP);
#endif /* YOUNGS_INCREASE */

#if (REFORM_SPRINGS == 2)
	if (p->iColor < 0) /* can safely skip wall particles */
	        return;
#endif /* REFORM_SPRINGS == 2 */

	/* loop over this particle's neighbors - check separation to see if cohesive forces should be applied */
	for (j=0;j<nSmooth;j++) {
	        pn = nnList[j].pPart;
#if (REFORM_SPRINGS == 1)
		if (p->iOrder >= pn->iOrder || (p->iColor < 0 && pn->iColor < 0)) /* perform on neighbor's call if neighbor has lower iOrder number, ignore if neighbor is particle itself, save disk space by skipping wall-particle to wall-particle spring assignment */
#elif (REFORM_SPRINGS == 2) /* REFORM_SPRINGS == 1 */
		if (p->iOrder >= pn->iOrder || pn->iColor < 0) /* can safely skip neighbors that are wall particles */
#elif (REFORM_SPRINGS == 3) /* REFORM_SPRINGS == 1 */
		if (p->iOrder >= pn->iOrder || (p->iColor >= 0 && pn->iColor >= 0) || (p->iColor < 0 && pn->iColor < 0)) /* allow spring formation only between one non-wall particle and one wall particle */
#endif /* REFORM_SPRINGS == 1 */
		        continue;

		/* compare stress to stress limit */
		assert(nnList[j].fDist2 > 0.0);
		d = sqrt(nnList[j].fDist2); /* could optimize by working with square instead */
		dEffRad = 0.5*(RADIUS(p) + RADIUS(pn)); /* effective radius = average radius */
		dStrain = (d/dEffRad - ZERO_STRAIN_LENGTH)/ZERO_STRAIN_LENGTH;
		/* if (dStrain < 0.) dStrain = 0.; */ /* uncomment to turn repulsion OFF */
		if (dStrain > SP->dMeanStrainLimit)
		        bApplyForce = 0;
		else
		        bApplyForce = 1;
		bFoundNbr = 0;
		/* check to see if this neighbor is listed as one of particle's springs */
		for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
		        ps = &p->springs[i];
			if (ps->iOrder == pn->iOrder) {
			        assert(!bFoundNbr); /* ensure no duplicate springs */
				bFoundNbr = 1;
				iSpring = i;
			}
		}
		if (bFoundNbr && bApplyForce) {
		        ps = &p->springs[iSpring];
#ifdef YOUNGS_INCREASE
			/* HACK to make Y increase with time - (incompatable with dYoungsStdDev != 0) */
			if (ps->fYoungsModulus < SP->dMeanYoungsModulus) {
			       ps->fYoungsModulus += YOUNGS_INCREASE_PER_STEP*SP->dMeanYoungsModulus;
			       if (ps->fYoungsModulus > SP->dMeanYoungsModulus)
				       ps->fYoungsModulus = SP->dMeanYoungsModulus;
			}
#endif /* YOUNGS_INCREASE */
		}
		else if (bFoundNbr && !bApplyForce) {
			ps = &p->springs[iSpring];
			ps->iOrder = -1;
			ps->fZeroStrainLength = 0.;
			continue;
		}
		else if (!bFoundNbr && bApplyForce) {
			for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
				ps = &p->springs[i];
				if (ps->iOrder == -1) {
					iSpring = i;
					ps->iOrder = pn->iOrder;
					ps->fZeroStrainLength = ZERO_STRAIN_LENGTH;
#ifdef YOUNGS_INCREASE
					ps->fYoungsModulus = YOUNGS_INCREASE_PER_STEP*SP->dMeanYoungsModulus; /* HACK to allow Y to increase in time */
#else /* YOUNGS_INCREASE */
					ps->fYoungsModulus = 0.;
					if (SP->dMeanYoungsModulus == 0.)
					        break; /* skip if Young's modulus is zero (doing this here so that damping forces might still be applied in such a case) */
					while (ps->fYoungsModulus <= 0. || ps->fYoungsModulus >= (SP->dMeanYoungsModulus + SP->dMeanYoungsModulus)) /* avoid negatives, maintain mean by requiring: 0 < particle's Y < 2*Mean_Y */
					        ps->fYoungsModulus = SP->dMeanYoungsModulus + SP->dYoungsStdDev*(double)randGaussian();
#endif /* YOUNGS_INCREASE */
					ps->fStressLimit = SP->dMeanStrainLimit*ps->fYoungsModulus; /* DEBUG: for now, no variation is allowed in strain limit using REFORM_SPRINGS...may want to change this (2/2) */
					break;
				}
			}
		}
		else if (!bFoundNbr && !bApplyForce)
			continue;

		ps = &p->springs[iSpring]; /* ps now points to spring connecting this particle to its neighbor */
		dStress = ps->fYoungsModulus*dStrain;

		vx = p->vPred[0] - pn->vPred[0];
		vy = p->vPred[1] - pn->vPred[1];
		vz = p->vPred[2] - pn->vPred[2];
					
		/* define intermediate terms to get units of acceleration per length in dRestoringAccel */
		dCrossSection = M_PI*dEffRad*dEffRad*SP->dPorosityFactor;
		dConvert_CrossSection = dConvert*dCrossSection;
		dRestoringQuantity = dConvert_CrossSection*dStress;
		dRestoringAccel_Particle = dRestoringQuantity/(p->fMass*d); /* this particle */
		dRestoringAccel_Neighbor = dRestoringQuantity/(pn->fMass*d); /* neighbor */

		/* calculate critical damping and damping forces if dDamp is nonzero, otherwise skip to save time */
		if (SP->dDamp) {
			dReducedMass = (p->fMass*pn->fMass)/(p->fMass+pn->fMass);
			dDampingCrit = sqrt(4.*dReducedMass*dConvert_CrossSection*ps->fYoungsModulus/(ps->fZeroStrainLength*dEffRad)); /* critical damping value */
			dDampingQuantity = SP->dDamp*dDampingCrit*(vx*nnList[j].dx+vy*nnList[j].dy+vz*nnList[j].dz);
			dDampingAccel_Particle = dDampingQuantity/(nnList[j].fDist2*p->fMass); /* this particle */
			dDampingAccel_Neighbor = dDampingQuantity/(nnList[j].fDist2*pn->fMass); /* neighbor */
		}
		else
			dDampingAccel_Particle = dDampingAccel_Neighbor = 0.;

		/*** this particle ***/
		p->dDeltaAccel[0] -= (dRestoringAccel_Particle+dDampingAccel_Particle)*nnList[j].dx; /* sign: (dx,dy,dz) = r - rn */
		p->dDeltaAccel[1] -= (dRestoringAccel_Particle+dDampingAccel_Particle)*nnList[j].dy;
		p->dDeltaAccel[2] -= (dRestoringAccel_Particle+dDampingAccel_Particle)*nnList[j].dz;

		/*** neighbor ***/
		pn->dDeltaAccel[0] += (dRestoringAccel_Neighbor+dDampingAccel_Neighbor)*nnList[j].dx;
		pn->dDeltaAccel[1] += (dRestoringAccel_Neighbor+dDampingAccel_Neighbor)*nnList[j].dy;
		pn->dDeltaAccel[2] += (dRestoringAccel_Neighbor+dDampingAccel_Neighbor)*nnList[j].dz;
	}
}
#endif /* REFORM_SPRINGS */

#endif /* SPRINGS */

#ifdef DEM

//#define DEBUG_DEM /*DEBUG!*/ /* define for debugging */

int DEBUG_TRACE(int iOrder);
#ifdef DEBUG_DEM
//#define DEBUG_TRACE(iOrder) (iOrder % 100000 == 0)
#define DEBUG_TRACE(iOrder) ((iOrder) == 6429)
#else
#define DEBUG_TRACE(iOrder) (0)
#endif

void AssignDEM(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*DEBUG why is this a smooth routine? it doesn't loop over neighbors*/

	int i;

	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);

	/* particle-particle DEM element */

	for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++)
		demWipeOverlapElement(&((PARTICLE *)p)->overlaps[i]); /* array order doesn't matter */

#ifdef WALLS
	/* particle-wall DEM element */

	for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;i++)
		demWipeOverlapElement(&((PARTICLE *)p)->walloverlaps[i]); /* array order doesn't matter */
#endif /* WALLS */
	}

void initDoDEMParticle(void *p)
{

	if (DEBUG_TRACE(((PARTICLE *)p)->iOrder)) fprintf(stderr,"Entering initDoDEMParticle(): iOrder = %i\n",((PARTICLE *)p)->iOrder);

	/* initialize each particle */

	int k;

	((PARTICLE *)p)->dPressure = 0.; /* reset for next step */

	for (k=0;k<3;k++) {
		((PARTICLE *)p)->dDeltaAccel[k] = 0.;
		((PARTICLE *)p)->wDot[k] = 0.;
#ifdef AGGS
        ((PARTICLE *)p)->DEM_aggaccels[k] = 0.;
        ((PARTICLE *)p)->DEM_aggtorques[k] = 0.;
#endif
		}
	}

void initDoDEM(void *p)
{

	/* initialize each cached particle (may not be necessary?...) */

	if (DEBUG_TRACE(((PARTICLE *)p)->iOrder)) fprintf(stderr,"Entering initDoDEM(): iOrder = %i\n",((PARTICLE *)p)->iOrder);

	int k;

	((PARTICLE *)p)->dPressure = 0.;
	for (k=0;k<3;k++) {
		((PARTICLE *)p)->dDeltaAccel[k] = 0.;
		((PARTICLE *)p)->wDot[k] = 0.;
#ifdef AGGS 
        ((PARTICLE *)p)->DEM_aggaccels[k] = 0.;
        ((PARTICLE *)p)->DEM_aggtorques[k] = 0.;
#endif /* AGGS */
		}

	/*
	for (k=0;k<MAX_NUM_OVERLAPS_PER_PARTICLE;k++)
		((PARTICLE *)p)->overlaps[k].iOrder = ((PARTICLE *)p)->overlaps[k].iOrder;
	*/
	}

void combDoDEM(void *p1,void *p2)
{

	/* update each particle with accumulated information contained in the cache */

	if (DEBUG_TRACE(((PARTICLE *)p1)->iOrder)) fprintf(stderr,"Entering combDoDEM(): iOrder = %i\n",((PARTICLE *)p1)->iOrder);

	int k;

	((PARTICLE *)p1)->dPressure += ((PARTICLE *)p2)->dPressure;
	((PARTICLE *)p1)->iColor = ((PARTICLE *)p2)->iColor; /* **Hack to give the ability to particles to delete their neighbors*/

	for (k=0;k<3;k++) {
		((PARTICLE *)p1)->dDeltaAccel[k] += ((PARTICLE *)p2)->dDeltaAccel[k];
		((PARTICLE *)p1)->wDot[k] += ((PARTICLE *)p2)->wDot[k];
#ifdef AGGS 
        ((PARTICLE *)p1)->DEM_aggaccels[k] += ((PARTICLE *)p2)->DEM_aggaccels[k];
        ((PARTICLE *)p1)->DEM_aggtorques[k] += ((PARTICLE *)p2)->DEM_aggtorques[k];
#endif /* AGGS */
		}

	/* for deleted particles, may want to do this w/o combiner cache - only doing this for iOrder to save time, other variables get reset upon new contact */
	/*
	for (k=0;k<MAX_NUM_OVERLAPS_PER_PARTICLE;k++)
		((PARTICLE *)p1)->overlaps[k].iOrder = ((PARTICLE *)p2)->overlaps[k].iOrder;
	*/
	}

/* #define DEM_PRESSURE */
void postDoDEM(PARTICLE *p,SMF *smf)
{
	/* save new accelerations for each particle */

	/* (note: we use dDeltaAccel because we may have multiple processors updating the same particles) */

	if (DEBUG_TRACE(((PARTICLE *)p)->iOrder)) fprintf(stderr,"Entering postDoDEM(): iOrder = %i(%i)\n",((PARTICLE *)p)->iOrder,smf->pkd->idSelf);

	if (p->iColor == 666) pkdDeleteParticle(smf->pkd,p); /* **Hack to give the ability to particles to delete their neighbors*/

	int k;

	for (k=0;k<3;k++)
		p->a[k] += p->dDeltaAccel[k];

#ifdef AGGS
	if (IS_AGG(p))
		for (k=0;k<3;k++)
			p->wDot[k] = 0.0; /* just in case */
#endif /* AGGS */

#ifdef DEM_PRESSURE 
	/* recolor particles based upon the sum of the normal forces that act on them --> if we keep this, put into dem.c */
	if (p->dPressure == 0.) {
		p->iColor = 2;
		}
	else {
		p->dPressure /= (RADIUS(p)*RADIUS(p));
		p->iColor = floor(30. + 150.*p->dPressure);
		if (p->iColor > 255) p->iColor = 255;
		}
#endif
	}

#ifdef CHARGE

/* following function disabled! problem: how to remove contribution of
   neighbors from mean-field contribution to particle's acceleration? */

/*
void _do_charge_force(PARTICLE *p,const PARTICLE *pn,const CHARGE_PARAMS *CP)
{
	}
*/

void _do_charge_xfer(PARTICLE *p,double dRadius,double dRadiusNbr)
{
	const double eps = 3.0e-5; /* should really depend on # timesteps per collision */
	const double sig = 2.238e36; /* approx transport coefficient in pkdgrav units */

	/* applies charge transfer algorithm for overlapping (colliding) particles */

	p->dCharge += eps*(0.5*sig*(dRadius*dRadius - dRadiusNbr*dRadiusNbr) - p->dCharge);
	}

#endif /* CHARGE */

/*#define DEM_OVERLAP_OUTPUT*/ /* uncomment for diagnostic overlap output */
#ifdef DEM_OVERLAP_OUTPUT
#define nParticles 1 /*DEBUG! must hardcode the number of particles for correct angmom/energy/overlap data files output */
#endif /* DEM_OVERLAP_OUTPUT */

void DoDEM(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** some parameters:
	**
	** p:       particle
	** pn:      neighbor
	** r1:      particle's radius
	** r2:      neighbor's radius
	** l1:      scalar dist from particle's center to point of intersection (lever-arm)
	** l2:      scalar dist from neighbor's center to point of intersection
	** d:       scalar dist between centers
	** x:       penetration (scalar)
	** kn:      (tensile) spring constant (defines max penetration, might use [typical impact energy / maximum desired penetration^2] --> [E_k / x_max^2], where x_max ~ 1% of Radii)
        ** kt:      (shear) spring constant
	*/
	/*
	** spring-dashpot-slider rolling model:
	** kr:         rolling spring coefficient
	** Cr:         rolling viscous damping coefficient
	** MrCrit:     critical rolling torque
	** dReducedRadius: effective radius
	** beta_s:     the shape parameter used in Jiang et al. (2015)
	** vRelW[3]:   relative angular velocity
	** vRelRol[3]: relative rolling rate
	**
	** spring-dashpot-slider twisting model:
	** kq:         twisting spring coefficient
	** Cq:         twisting viscous damping coefficient
	** MtCrit:     critical twisting torque
	** vRelTwi[3]: relative twisting rate
	*/

	enum {DEM_NOERR,DEM_MINOR,DEM_MAJOR,DEM_ERROR}; /* for overlap warning/error reporting */

	PARTICLE *pn=NULL; /* initialized to suppress compiler warning */
	DEM_ELEMENT *pe;
	int i,j,bFoundNbr,bApplyForce,iOverlap,bUseContactRadius,bUseLegacyRoll;
	double dReducedMass,dReducedRadius=0.0,dContactRadius=0.0,dRollLever=0.0,dRollLeverNeighbor=0.0;
#ifdef DEM_ROTATION_DASHPOT
	double kr,Cr,MrCrit,kq,Cq,MtCrit,beta_s,dRelTwi;
	double dtemp1,dtemp2;
	Vector vRelW,vRelRol,vRelTwi;
	Vector vtemp1;
#else
	double dRollingFriction,dTwistRate,dTwistingFriction=0.,Rt;
	Vector r,rt;
#endif /* DEM_ROTATION_DASHPOT */
	Vector vRft,vTfn;
	Vector n,t,v,s1,s2,s,u,un,ut,Fn,Ft,Ftn,vDum,vDum2;
#ifdef AGGS 
    Vector vBodyPos,vAggLever,vAggForces,vAggTorques,vAggTemp,vAggTemp2,vRollForces,vTwistForces;
#endif
	double m1,m2,r1,r2,l1,l2,r1sq,r2sq,kn,kt,i1,i2,d,x,Un,Ut,mu_s,mu_r,mu_t,Cn,Ct;
	double dDelta,dEpsN,dEpsT,dTangentialSpringDrag,dSqrtReducedMass,F_outershell=0.;
	double a=0.0,a1,a2,b,b1,b2,N=HUGE_VAL,T,F,F2,dLnEpsN=0.0,dLnEpsT,dLnEpsNsq=-1.0,dLnEpsTsq;
	static double CnPreFac,CtPreFac;

#ifdef WALLS
	int bOutcome;
#endif

#ifdef DEM_TWOLAYERS
	double kn_outer,kt_outer,kn_inner,kt_inner,dInBnd,dInBndFrac,dEffRad,x_inner;
	static double CnOuterPreFac,CtOuterPreFac,CnInnerPreFac,CtInnerPreFac;
#endif
#ifdef DEM_COHESION
	int cohe_model;
    int ncohe_coeff; /* number of cohesive coefficient */
	double cohe_coeff; /* cohesive coefficient */
	double dtemp_cohe;
#endif

#ifdef DEM_ROTATION_DASHPOT
	/*DEBUG! (DCR 6/7/18) this is to suppress "set but not used"
	  compilation warnings.  The issue is that these variables are
	  currently only used if DEM_ROTATION_DASHPOT is *not* defined,
	  but they are set at the same time as dContactRadius, which *can*
	  be used.  See docs/README_SSDEM for more information.
	  Obviously, this needs to be cleaned up.*/
	(void) dRollLever;
	(void) dRollLeverNeighbor;
#endif
	
	if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Entering DoDEM(): iOrder = %i(%i)\n",p->iOrder,smf->pkd->idSelf);

#ifdef DEM_OVERLAP_OUTPUT /* prints to file duration of each overlap */
	static unsigned long int liCounter = 0; ++liCounter;
	static unsigned long int liStep = 0; if (liCounter % nParticles == 1) ++liStep;
	static FILE *overlapOUT;
	if (liCounter == 1) overlapOUT = fopen( "overlap.dat", "w" );
#endif /* DEM_OVERLAP_OUTPUT */

	dEpsN = smf->FO.DP.dEpsN;
	dEpsT = smf->FO.DP.dEpsT;
	dDelta = smf->FO.DP.dDelta;
	kn = smf->FO.DP.dKn;
	kt = smf->FO.DP.dKt;
	mu_s = smf->FO.DP.dMuS; /* coefficient of static friction */
	mu_r = smf->FO.DP.dMuR; /* coefficient of rolling friction */
	mu_t = smf->FO.DP.dMuT; /* coefficient of twisting friction */
#ifdef AGGS	/* Currently, AGGS does not incorporate rolling/twisting */
	mu_r = 0.0;
	mu_t = 0.0;
#endif
#ifdef DEM_ROTATION_DASHPOT
	beta_s = smf->FO.DP.dBetaS; /* the shape parameter used in spring-dashpot-slider rotation model */
#endif /* DEM_ROTATION_DASHPOT */
	bUseLegacyRoll = smf->FO.DP.bUseLegacyRoll; /* use original rolling friction (nonuniform torques) */
	bUseContactRadius = smf->FO.DP.bUseContactRadius; /* use contact radius instead of reduced radius in rolling friction calculation */
#ifdef DEM_COHESION
	cohe_model = smf->FO.DP.iCohesionModel;
    ncohe_coeff = smf->FO.DP.sCohesionVarParams.nCohesiveCoeff;
#endif
#ifdef DEM_TWOLAYERS
	dInBndFrac = smf->FO.DP.dInnerOverlapBndry;
	kn_outer = smf->FO.DP.dKnOuter;
	kt_outer = smf->FO.DP.dKtOuter;
	kn_inner = smf->FO.DP.dKnInner;
	kt_inner = smf->FO.DP.dKtInner;
#endif
	dTangentialSpringDrag = 0.; /*DEBUG! overrides any read-in value*/
	r1 = RADIUS(p);
	r1sq = r1*r1;

	mdlassert(smf->pkd->mdl,smf->FO.iForceOverrideOption == FO_STRENGTH);

	const double pi_sq = M_PI*M_PI;
	static int bCnCtPreFacCalculated = 0; /* since dEpsN and dEpsT are identical for all particles in sim */
	if (!bCnCtPreFacCalculated) {
		/* damping term: normal */
		mdlassert(smf->pkd->mdl,dEpsN > 0.0);
		dLnEpsN = log(dEpsN);
		dLnEpsNsq = dLnEpsN*dLnEpsN;
		a = pi_sq + dLnEpsNsq;
		CnPreFac = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn/a);
		CnPreFac += CnPreFac;
#ifdef DEM_TWOLAYERS
		CnOuterPreFac = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn_outer/a);
		CnOuterPreFac += CnOuterPreFac;
		CnInnerPreFac = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn_inner/a);
		CnInnerPreFac += CnInnerPreFac;
#endif

		/* damping term: tangential */
		mdlassert(smf->pkd->mdl,dEpsT > 0.0);
		dLnEpsT = log(dEpsT); /* could this be log(fabs(dEpsT))? */
		dLnEpsTsq = dLnEpsT*dLnEpsT;
		a = pi_sq + dLnEpsTsq;
		CtPreFac = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt/a);
		CtPreFac += CtPreFac;
#ifdef DEM_TWOLAYERS
		CtOuterPreFac = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt_outer/a);
		CtOuterPreFac += CtOuterPreFac;
		CtInnerPreFac = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt_inner/a);
		CtInnerPreFac += CtInnerPreFac;
#endif
		(void) printf("CnPreFac and CtPreFac calculated on pid %i\n",smf->pkd->idSelf); /*DEBUG*/
		bCnCtPreFacCalculated = 1;
		}

	/*DEBUG*/
	/*
	  static int count;
	  ++count;
	*/

	/*DEBUG!!!
	  if (p->iOrder == 59921 || p->iOrder == 62359) {
	  int ii;
	  printf("%i [%g]: x=%23.16e y=%23.16e vy=%23.16e\n",p->iOrder,smf->dTime,p->r[0],p->r[1],p->v[1]);
	  printf("%i [%g] nbrs:",p->iOrder,smf->dTime);
	  for (ii=0;ii<nSmooth;ii++)
	  printf(" %i",nnList[ii].pPart->iOrder);
	  printf("\n%i [%g] ovlp:",p->iOrder,smf->dTime);
	  for (ii=0;ii<MAX_NUM_OVERLAPS_PER_PARTICLE;ii++)
	  printf(" %i",p->overlaps[ii].iOrder);
	  printf("\n");
	  }
	*/

/* #define PHOBOS */
#ifdef PHOBOS   /* here's a hack to put in tidal forces accnt`ing for Mars and to delete particles that hit it */
	double dW_orb,dMars_absposSq,dMars_relposSq,dMars_absposSqINV,dMars_relposSqINV,dMartianGravFac;
	Vector dMars_abspos,dMars_relpos,dMars_absposNORM,dMars_relposNORM;
	dW_orb = 1.1452821663139312e+03; /* freq in sec / 2pi --> 2.280329865e-4*5022440.7616851525 */
	dMars_abspos[0] = 6.2674634687629510e-05*cos(dW_orb*smf->dTime); /* vectors point to Mars */
	dMars_abspos[1] = 6.2674634687629510e-05*sin(dW_orb*smf->dTime);
	dMars_abspos[2] = 0.;
	vectorSub(dMars_abspos,p->r,dMars_relpos);
	dMars_absposSq = 3.9281098332278120e-09; /* vectorMagSq(dMars_abspos) */
	dMars_relposSq = vectorMagSq(dMars_relpos);
	if (dMars_relposSq < (2.2657388467760265e-05 + r1)*(2.2657388467760265e-05 + r1)) {
		p->iColor = 666;  /* death */
		return;
		}
	dMars_absposSqINV = 2.5457536638640234e+08; /* 1./vectorMagSq(dMars_abspos) */
	dMars_relposSqINV = 1./dMars_relposSq;
	vectorScale(dMars_abspos,sqrt(dMars_absposSqINV),dMars_absposNORM);
	vectorScale(dMars_relpos,sqrt(dMars_relposSqINV),dMars_relposNORM);

	/* apply martian tidal forces */
	dMartianGravFac = 3.2264216234999999e-07; /* GM = 1. * 6.4171e23*5.02785e-31 */
	for (i=0;i<3;i++) {
		p->dDeltaAccel[i] += dMartianGravFac*(dMars_relposNORM[i]*dMars_relposSqINV - dMars_absposNORM[i]*dMars_absposSqINV);
		/* printf("%d %.16e %.16e %.16e\n",p->iOrder,p->r[i],dMars_relpos[i],dMars_abspos[i]); */
		}
#endif /* PHOBOS */

	//printf("%d %.8g\n",p->iOrder,p->fBall2);
	//int bGhost=0;
	/*
	  FILE *f_dbr_26nbr, *f_dbr_917nbr, *f_dbr_All1, *f_dbr_All;
	  if ((p->iOrder == 26 || p->iOrder == 917) && smf->dTime == 0*dDelta) {
	  f_dbr_26nbr = fopen("List26nbr.out","a");
	  f_dbr_917nbr = fopen("List917nbr.out","a");
	  
	  for (ii=0;ii<nSmooth;ii++){
	  for (j=0; j < smf->pkd->nLocal; j++){
      pn = &smf->pkd->pStore[j];
      if (pn->iOrder == nnList[ii].pPart->iOrder){
	  if (p->iOrder == 26)
	  fprintf(f_dbr_26nbr," %.8g %i %23.16e %23.16e %23.16e %23.16e\n",smf->dTime, pn->iOrder, pn->r[0]*1.49597871e13, pn->r[1]*1.49597871e13, pn->r[2]*1.49597871e13, sqrt(nnList[ii].fDist2)*1.49597871e13);
	  if (p->iOrder == 917)
	  fprintf(f_dbr_917nbr," %.8g %i %23.16e %23.16e %23.16e %23.16e\n",smf->dTime, pn->iOrder, pn->r[0]*1.49597871e13, pn->r[1]*1.49597871e13, pn->r[2]*1.49597871e13, sqrt(nnList[ii].fDist2)*1.49597871e13);
	  break;
	  }
      }
	  }
	  fclose(f_dbr_26nbr);
	  fclose(f_dbr_917nbr);
	  }
	  
	  if (p->iOrder == 26 && smf->dTime == 0*dDelta){
	  f_dbr_All1     = fopen("ListAll1.out","a");
	  for (j=0; j < smf->pkd->nLocal; j++){
	  pn = &smf->pkd->pStore[j];
	  fprintf(f_dbr_All1," %.8g %i %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",smf->dTime, pn->iOrder, pn->r[0], pn->r[1], pn->r[2], pn->v[0], pn->v[1], pn->v[2], pn->w[0], pn->w[1], pn->w[2]);
	  }
	  fclose(f_dbr_All1);
	  }
	  
	  if (p->iOrder == 26 && smf->dTime == 1*dDelta){
	  f_dbr_All     = fopen("ListAll.out","a");
	  for (j=0; j < smf->pkd->nLocal; j++){
	  pn = &smf->pkd->pStore[j];
	  fprintf(f_dbr_All," %.8g %i %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",smf->dTime, pn->iOrder, pn->r[0], pn->r[1], pn->r[2], pn->v[0], pn->v[1], pn->v[2], pn->w[0], pn->w[1], pn->w[2]);
	  }
	  fclose(f_dbr_All);
	  }
	*/

	/* check to see if all overlapping particles are in neighbor list */
	for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++) {
		pe = &p->overlaps[i];
		if (pe->iOrder != -1) {
			for (j=0;j<nSmooth;j++) {
				pn = nnList[j].pPart;
				if (pn->iOrder == pe->iOrder)
					break;
				}
			if (j == nSmooth) {
				const double au = 1.495978707e11; /* 1 au in m, from ss.h */
				double dPnRadius=RADIUS(pn);
				int ii;
				fprintf(stderr,"DEM ERROR: particle %i in overlap on previous step no longer in neighbor list of particle %i (pid=%i)\n",pe->iOrder,p->iOrder,smf->pkd->idSelf);
				printf("%i [time %.8g]: radius=%23.16e x=%23.16e y=%23.16e z=%23.16e fBall=%23.16e m\n",p->iOrder,smf->dTime,dPnRadius*au,p->r[0]*au,p->r[1]*au,p->r[2]*au,sqrt(p->fBall2)*au);
				printf("%i [time %.8g] nbrs (d in m):",p->iOrder,smf->dTime);
				for (ii=0;ii<nSmooth;ii++)
					printf(" (%i, %g)",nnList[ii].pPart->iOrder,sqrt(nnList[ii].fDist2)*au);
				printf("\n%i [time %.8g] ovlp:",p->iOrder,smf->dTime);
				for (ii=0;ii<MAX_NUM_OVERLAPS_PER_PARTICLE;ii++)
					printf(" %i",p->overlaps[ii].iOrder);
				printf("\n");
				for (j=0;j<smf->pkd->nLocal;j++){
					pn = &smf->pkd->pStore[j];
					dPnRadius = RADIUS(pn);
					if (pn->iOrder == pe->iOrder)
						printf("Funny Particle: %i [time %.8g]: radius=%23.16e x=%23.16e y=%23.16e z=%23.16e fBall=%23.16e m\n",pn->iOrder,smf->dTime,dPnRadius*au,pn->r[0]*au,pn->r[1]*au,pn->r[2]*au,sqrt(p->fBall2)*au);
					}
#ifdef PHOBOS
				demWipeOverlapElement(pe); /* simply erase overlap history of lost neighbor */
#else
				mdlassert(smf->pkd->mdl,j < nSmooth);
#endif /* PHOBOS */
				}
			}
		}

	/* loop over neighbors */

	if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Starting loop over neighbors: iOrder = %i(%i)\n",p->iOrder,smf->pkd->idSelf);

	/*DEBUG
	** It's too complicated to deal with predicted velocities and
	** spins of particles that are either stuck to walls or part of
	** rigid aggregates.  So for now we just use the unpredicted
	** values.  This is not as accurate, but doing better is probably
	** not worth the effort.  This is a brute-force hack!
	*/
#ifdef WALLS
	if (PARTICLE_STUCK(p)) {
		vectorCopy(p->v,p->vPred);
		vectorCopy(p->w,p->wPred);
		}
#endif /* WALLS */

	for (j=0;j<nSmooth;j++) {
		pn = nnList[j].pPart;

#ifdef CHARGE
		/*_do_charge_force(p,pn,&smf->CP.CP <-- need to add to smf!);*/
		/*also, only call if pn != p*/
#endif

#ifdef WALLS
		if (PARTICLE_STUCK(p) && PARTICLE_STUCK(pn))
			continue;
		if (PARTICLE_STUCK(pn)) { /*DEBUG see note above about predicted velocities and spins for stuck particles */
			vectorCopy(pn->v,pn->vPred);
			vectorCopy(pn->w,pn->wPred);
			}
#endif /* WALLS */

#ifdef AGGS // JCM DEBUG
        if (IS_AGG(p) && IS_AGG(pn) && AGG_IDX(p) == AGG_IDX(pn))
			continue;
#endif /* AGGS */

		bApplyForce = bFoundNbr = 0;
		iOverlap = -1; /* for safety */

#ifdef DEM_OVERLAP_OUTPUT /* this will only work if nParticles is small, otherwise it should seg fault b/c it can't allocate memory */
		static double x_max[nParticles][nParticles];
		if (liStep == 1)
			x_max[p->iOrder][j] = 0;
#endif /* DEM_OVERLAP_OUTPUT */

		/* perform on neighbor's call if neighbor has lower iOrder number or skip if neighbor is particle */
		if (p->iOrder >= pn->iOrder) /*DEBUG could check if pn->iOrder == p->iOrder right at the start...*/
			continue;

		/* check to see if particle is overlapping with neighbor */
        /* Exactly Overlapping? Redundant? -RB 4/20/2015*/
		mdlassert(smf->pkd->mdl,nnList[j].fDist2 > 0.0);

		r2 = RADIUS(pn);
#ifdef GLASS_BEADS
		SPRING *ps;
		double dConjoinedSize;
		for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) { /* might save time by only doing this loop once to find the spring (another springs loop below) */
			ps = &p->springs[i];
			if (ps->iOrder == pn->iOrder && ps->fZeroStrainLength == 0. && ps->fStressLimit < 0) {
				dConjoinedSize = ps->fStressLimit; /* separation between centers at time of tangential fracture */
				r1 = -r1*dConjoinedSize/(r1 + r2); /*DEBUG:check this for KE gain*/
				r2 = -r1 - dConjoinedSize;
				r1sq = r1*r1;
				}
			break;
			}
#endif /* GLASS_BEADS */

		/* if there's an overlap, a force should be applied */
		if ((r1 + r2)*(r1 + r2) > nnList[j].fDist2)
			bApplyForce = 1;

#ifdef CHARGE
		if (bApplyForce)
			_do_charge_xfer(p,r1,r2);
#endif /* CHARGE */

		for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++) {
			pe = &p->overlaps[i];
			if (pe->iOrder == pn->iOrder) {
				mdlassert(smf->pkd->mdl,!bFoundNbr); /* ensure no duplicate overlaps */
				bFoundNbr = 1;
				iOverlap = i;
				}
			}

		/**Hack to give particles the ability to delete their neighbors next two lines (Red Mars)*/
		/*
		if (bApplyForce && p->iColor == 2 && pn->iColor > 0) pn->iColor = 666;
			continue;
		*/

		/* four cases: */
		if (bFoundNbr && bApplyForce) {
			pe = &p->overlaps[iOverlap];
#ifndef DEM_TWOLAYERS
			++pe->liOverlapCounter; /* if two layers, take care of the counter once determined which layer we're in */
#endif
			}

		else if (bFoundNbr && !bApplyForce) {
			pe = &p->overlaps[iOverlap];
#ifdef DEM_OVERLAP_OUTPUT
			fprintf(overlapOUT,"%g %d %d %i %e\n",smf->dTime,p->iOrder,pe->iOrder,pe->liOverlapCounter,x_max[p->iOrder][j]);
			x_max[p->iOrder][j] = 0.;
#endif /* DEM_OVERLAP_OUTPUT */
			demWipeOverlapElement(pe);
			continue;
			}

		else if (!bFoundNbr && bApplyForce) {
			for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++) {
				pe = &p->overlaps[i];
				if (pe->iOrder == -1) {
					iOverlap = i;
					pe->iOrder = pn->iOrder;
					vectorZero(pe->vShear);
					vectorZero(pe->vnOld);
#ifdef DEM_ROTATION_DASHPOT
					vectorZero(pe->vRoll);
                                        pe->dTwist = 0.0;
#endif /* DEM_ROTATION_DASHPOT */
#ifdef DEM_DIAG
					vectorZero(pe->vFn);
					vectorZero(pe->vFt);
#endif /* DEM_DIAG */
#ifndef DEM_TWOLAYERS
					pe->liOverlapCounter = 1; /* if two layers, take care of the counter once determined which layer we're in */
#endif
					break;
					}
				}
#ifdef PHOBOS
			if (i >= MAX_NUM_OVERLAPS_PER_PARTICLE) (void) fprintf(stderr,"DEM WARNING (pid=%i;nid=%i) [T=%g]: overlaps exceed nSmooth! no free spots, skipping contact\n",smf->pkd->idSelf,pn->iOrder,smf->dTime);
			continue;
#else
			mdlassert(smf->pkd->mdl,i < MAX_NUM_OVERLAPS_PER_PARTICLE); /* ensure free spot was indeed found and used */
#endif /* PHOBOS */
			}
		else if (!bFoundNbr && !bApplyForce)
			continue;

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Starting force calculation on neighbor: iOrder = %i(%i), neighbor = %i, neighbor in DEM struct = %i\n",p->iOrder,smf->pkd->idSelf,pn->iOrder,pe->iOrder);

		mdlassert(smf->pkd->mdl,iOverlap >= 0 && iOverlap < MAX_NUM_OVERLAPS_PER_PARTICLE);

		pe = &p->overlaps[iOverlap]; /* pe now points to the overlap of this particle and its neighbor */

		m1 = p->fMass;
		m2 = pn->fMass;
                r2sq = r2*r2;

#ifdef WALLS
		if (PARTICLE_STUCK(p))
			dReducedMass = m2;
		else if (PARTICLE_STUCK(pn))
			dReducedMass = m1;
		else
#endif /* WALLS */

		if (m1 == m2) /* save time if masses are equal */
			dReducedMass = 0.5*m1;
		else
			dReducedMass = m1*m2/(m1 + m2);
		dSqrtReducedMass = sqrt(dReducedMass);
		d = sqrt(nnList[j].fDist2);
		if (r1 == r2) /* save time if radii are equal */
			l1 = l2 = 0.5*d;
		else {
			l1 = (r1sq - r2sq + nnList[j].fDist2)/(d + d);
			l2 = d - l1;
#ifdef PHOBOS
			if (l1 < 0.0 || l2 < 0.0) {
				(void) fprintf(stderr,"DEM WARNING (pid=%i;nid=%i) [T=%g]: l1 = %.16e && l2 = %.16e\n",smf->pkd->idSelf,pn->iOrder,smf->dTime,l1,l2);
				if (l1 <= 0.0) {
					l1 = 0;
					l2 = d;
					}
				else if (l2 <= 0.0) {
					l1 = d;
					l2 = 0;
					}
				else
					mdlassert(smf->pkd->mdl,0);
				}
#else
			/* mdlassert(smf->pkd->mdl,l1 >= 0.0 && l2 >= 0.0); -- move after overlap check */
#endif /* PHOBOS */
			}
		x = r1 + r2 - d;
		if (x < 0.0) {
			(void) fprintf(stderr,"DEM ERROR (pid=%i;nid=%i) [T=%g]: x = %.16e, d = %.16e r1sq =%.16e, r2sq = %.16e, r1 = %.16e, r2 = %.16e, (r1 + r2)*(r1 + r2) = %.16e, nnList[j].fDist2 = %.16e bApplyForce = %d\n",
					smf->pkd->idSelf,pn->iOrder,smf->dTime,x,d,r1sq,r2sq,r1,r2,(r1 + r2)*(r1 + r2),nnList[j].fDist2,bApplyForce);
			mdlassert(smf->pkd->mdl,x >= 0.0);
			}
#ifdef DEM_TWOLAYERS /* if an inner boundary is defined where stiffness changes */
		if (r1 == r2) /* save time if radii are equal */
			dEffRad = r1;
		else
			dEffRad = 0.5*(r1 + r2);
		dInBnd = dEffRad*dInBndFrac;
		if (x > dInBnd) {
			kn = kn_inner;
			kt = kt_inner;
			x_inner = x - dInBnd;
			mdlassert(smf->pkd->mdl,x_inner >= 0.0);
			F_outershell = dInBnd*kn_outer;
			if (pe->liOverlapCounter <= 0)
				--pe->liOverlapCounter;
			else { /* scale static spring if boundary layer was just crossed (conserve energy in tang spring: 1/2 kt*(vShear)^2 */
#ifdef DEM_OVERLAP_OUTPUT
				fprintf(overlapOUT,"inner boundary crossed - %g %d %d %li %e\n",smf->dTime,p->iOrder,pe->iOrder,pe->liOverlapCounter,x_max[p->iOrder][j]);
#endif
				vectorScale(pe->vShear,sqrt(kt_outer/kt_inner),pe->vShear); /* vnOld remains the same */
				pe->liOverlapCounter = -1;
				}
			}
		else {
			kn = kn_outer;
			kt = kt_outer;
			x_inner = -1.;
			F_outershell = 0.;
			if (pe->liOverlapCounter >= 0)
				++pe->liOverlapCounter;
			else { /* scale static spring... */
#ifdef DEM_OVERLAP_OUTPUT
				fprintf(overlapOUT,"inner boundary crossed - %g %d %d %li %e\n",smf->dTime,p->iOrder,pe->iOrder,pe->liOverlapCounter,x_max[p->iOrder][j]);
#endif
				vectorScale(pe->vShear,sqrt(kt_inner/kt_outer),pe->vShear);
				pe->liOverlapCounter = 1;
				}
			}
#endif /* DEM_TWOLAYERS */

#ifdef DEM_OVERLAP_OUTPUT
		if (x > x_max[p->iOrder][j])
			x_max[p->iOrder][j] = x;
#endif /* DEM_OVERLAP_OUTPUT */

		{
		double r_min = (r1 <= r2 ? r1 : r2);
		int iWarningType = DEM_NOERR;
		if (x > smf->FO.DP.dErrorFrac*r_min)
			iWarningType = DEM_ERROR;
		else if (x > smf->FO.DP.dMajorFrac*r_min) {
			iWarningType = DEM_MAJOR;
#ifdef PHOBOS
			if (p->iColor == -1 || pn->iColor == -1) p->iColor = pn->iColor = -1;
#endif
			}
		else if (x > smf->FO.DP.dMinorFrac*r_min)
			iWarningType = DEM_MINOR;

		if (iWarningType != DEM_NOERR) {
#if (INTERNAL_WARNINGS != 0)
			static int nWarn = 1; /*DEBUG may not work in MPI?*/
#ifndef GLASS_BEADS
			if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0 || iWarningType != DEM_MINOR) {
				double v_mag = sqrt((p->v[0] - pn->v[0])*(p->v[0] - pn->v[0]) +
									(p->v[1] - pn->v[1])*(p->v[1] - pn->v[1]) +
									(p->v[2] - pn->v[2])*(p->v[2] - pn->v[2]));
				fprintf(stderr,"DEM %s #%i (pid=%i,T=%g) %s overlap: %i & %i, %f%% of%s radius (relative v = %g)\n",
						iWarningType == DEM_ERROR ? "ERROR" : "WARNING",nWarn,smf->pkd->idSelf,smf->dTime,
						iWarningType == DEM_MINOR ? "minor" : iWarningType == DEM_MAJOR ? "major" : "fatal",
						p->iOrder,pe->iOrder,100.0*x/r_min,r1 != r2 ? " smaller" : "",v_mag); /* max overlap = 200% */
				}
#endif
			++nWarn;
#endif /* INTERNAL_WARNINGS */
				mdlassert(smf->pkd->mdl,iWarningType != DEM_ERROR);
			}
		}
		a = 1./d;
		vectorSet(n,-nnList[j].dx*a,-nnList[j].dy*a,-nnList[j].dz*a); /* sign: n[i] = unit vector of [rn - r], points from particle's center toward neighbor's center */

		vectorCross(p->wPred,n,s1);
		vectorCross(pn->wPred,n,s2);
		mdlassert(smf->pkd->mdl,l1 >= 0.0 && l2 >= 0.0); /* should be caught earlier by overlap check */
		vectorScale(s1,l1,s1);
		vectorScale(s2,-l2,s2);

                Un = Ut = 0.;
                for (i=0;i<3;i++)
                        v[i] = pn->vPred[i] - p->vPred[i];

#ifdef SLIDING_PATCH
		/*DEBUG! is "v[i] = pn->vPred[i] - p->vPred[i]" above correct for SLIDING_PATCH?*/
                if (smf->PP.bPatch) {
                        /* adjust neighbor y-velocity if it's a ghost! */
                        /* note: DEM uses pn->v - p->v while collisions does the opposite, so the sign of the adjustment to v[1] is swapped below... */
                        if (p->r[0] > pn->r[0] && nnList[j].dx < 0.0)
                                v[1] -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
                        else if (p->r[0] < pn->r[0] && nnList[j].dx > 0.0)
                                v[1] += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
                        }
#endif /* SLIDING PATCH */

#ifdef DEM_ROTATION_DASHPOT
                for (i=0;i<3;i++) {
                        s[i] = s2[i] - s1[i];
                        u[i] = v[i] + s[i];
                        Un += u[i]*n[i];
                        }
                for (i=0;i<3;i++) {
                        un[i] = Un*n[i];
                        ut[i] = u[i] - un[i];
                        Ut += ut[i]*ut[i];
                        }
#else /* Directional constant torque models */
                Rt = 0.;
                if (mu_r != 0.) {
                        for (i=0;i<3;i++) {
                                s[i] = s2[i] - s1[i];
                                r[i] = s2[i] + s1[i]; /* rolling */
                                u[i] = v[i] + s[i];
                                Un += u[i]*n[i];
                                }

                        for (i=0;i<3;i++) {
                                un[i] = Un*n[i];
                                ut[i] = u[i] - un[i];
                                Ut += ut[i]*ut[i];
                                Rt += r[i]*r[i]; /* rolling */
                                }
                        Rt = sqrt(Rt);
                        }
                else { /* let's make 2 functions in dem.c that do this dependent on mu_r == 0? to clean this up [srs] */
                        for (i=0;i<3;i++) {
                                s[i] = s2[i] - s1[i];
                                u[i] = v[i] + s[i];
                                Un += u[i]*n[i];
                                }
                        for (i=0;i<3;i++) {
                                un[i] = Un*n[i];
                                ut[i] = u[i] - un[i];
                                Ut += ut[i]*ut[i];
                                }
/*DEBUG! if (p->iOrder == 89773) fprintf(stderr,"*** %i v after = %23.16e,%23.16e,%23.16e\n",p->iOrder,v[0],v[1],v[2]);*/
/*DEBUG! if (p->iOrder == 89773) fprintf(stderr,"*** %i un = %23.16e,%23.16e,%23.16e\n",p->iOrder,un[0],un[1],un[2]);*/
                        }
#endif /* (!)DEM_ROTATION_DASHPOT */
                Ut = sqrt(Ut);

		/* unit tangential velocity vector (ensure that this is needed, also want a un and unPred... - [srs]) */
		if (Ut != 0.0)
			vectorScale(ut,1./Ut,t);
		else
			vectorZero(t); /* avoid dividing by zero */

		/*
		** relative to neighbor - i.e. in frame where velocity at contact point of particle
		** is zero, this is the tangential velocity at the contact point of neighbor
		*/

#ifndef DEM_ROTATION_DASHPOT
		/* unit rolling vector */
		if (Rt != 0.0)
			vectorScale(r,1./Rt,rt); /* rt[i] = unit vector of velocity of spinning component at contact point */
		else /* save time and avoid dividing by zero */
			vectorZero(rt);
#endif  /* Directional constant torque models */

		/* damping terms */

#ifndef DEM_TWOLAYERS
		Cn = CnPreFac*dSqrtReducedMass;
		Ct = CtPreFac*dSqrtReducedMass; /* ...relating dEpsN,dEpsT to Cn,Ct comes from Cleary etal. '98 */
#else  /* if an inner boundary is defined where stiffness changes */
		if (x > dInBnd) {
			Cn = CnInnerPreFac*dSqrtReducedMass;
		   	Ct = CtInnerPreFac*dSqrtReducedMass;
			}
		else {
			Cn = CnOuterPreFac*dSqrtReducedMass;
		   	Ct = CtOuterPreFac*dSqrtReducedMass;
			}
#endif /* DEM_TWOLAYERS */

#ifdef WALLS
		demAssignParticleInertia(PARTICLE_STUCK(p),PARTICLE_STUCK(pn),m1,m2,r1,r2,r1sq,r2sq,&a1,&a2,&i1,&i2,&b1,&b2);
#else
		demAssignParticleInertia(0,0,m1,m2,r1,r2,r1sq,r2sq,&a1,&a2,&i1,&i2,&b1,&b2);
#endif

/* Initialization */
#ifndef DEM_ROTATION_DASHPOT
		if (Rt != 0. || mu_t != 0.) {
#endif
			if (bUseLegacyRoll) {
				dRollLever = l1;
				dRollLeverNeighbor = l2;
				}
			else if (bUseContactRadius) {
				a = r1sq - l1*l1;
				mdlassert(smf->pkd->mdl,a >= 0.0);
				dRollLever = dRollLeverNeighbor = dContactRadius = sqrt(a);
				}
			else {
				if (r1 == r2)
					dRollLever = dRollLeverNeighbor = dReducedRadius = r1;
				else
					dRollLever = dRollLeverNeighbor = dReducedRadius = r1*r2/(r1+r2);
				}
#ifndef DEM_ROTATION_DASHPOT
			}
#else
		if (!bUseContactRadius) {
			dReducedRadius= r1*r2/(r1+r2);
			dtemp1 = beta_s * dReducedRadius;
			}
		else
			dtemp1 = beta_s * dContactRadius;
		dtemp1 *= dtemp1;
		kr = kn * dtemp1; /* spring term */
		Cr = Cn * dtemp1; /* damping term */
		dtemp2 = 2.0 * dtemp1;
		kq = kt * dtemp2;
		Cq = Ct * dtemp2;
		for (i=0;i<3;i++)
			vRelW[i] = pn->wPred[i] - p->wPred[i];
		dRelTwi = vectorDot(n,vRelW);
		dtemp2 = 0.5*(l2-l1)/(l1*l2);
		vectorCross(n,u,vtemp1);
		for (i=0;i<3;i++) {
			vRelTwi[i] = dRelTwi * n[i];
			vRelRol[i] = vRelW[i] - vRelTwi[i] + dtemp2*vtemp1[i];
			}
#endif /* DEM_ROTATION_DASHPOT */

		/* Begin Hack to save time for mu_s == 0.
		if (mu_s == 0.) {
			for (i=0;i<3;i++) Fn[i] = -kn*x*n[i] + Cn*un[i];
			vectorZero(Ft);
			vectorZero(Ftn);
			}
			else { */
		/* End Hack to save time for mu_s == 0. */

		/* Tangential displacement vector */
		
		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Starting calculation on tangential displacement vector: iOrder = %i(%i), neighbor = %i\n",p->iOrder,smf->pkd->idSelf,pe->iOrder);

		if (vectorMagSq(pe->vnOld) != 0.0) {
#ifndef DEM_TWOLAYERS /* if boundary is crossed, vnOld will be defined and |liOverlapCounter| == 1 */
			//mdlassert(smf->pkd->mdl,pe->liOverlapCounter > 1); /* should not be first step in overlap *//*DEBUG assert can trip on restart from DEM file*/
#endif

			/*
			** What we really would want here is to track the motion of a
			** point between the current contact point and the equilibrium
			** contact point, then use the changes in the motion and
			** orientation of that point to determine the direction and
			** magnitude of the vector that points from the current
			** contact point to the equilibrium contact point (or vice
			** versa). For example, if one particle twists around the
			** current normal, but the other particle does not exhibit
			** this twisting relative to the motion at the midpoint, then
			** this other particle's equilibrium contact point should
			** really not change, though it would feel the compulsion to
			** rotate around its equilibrium contact point (where its
			** neighbor is).
			**
			** We'll take an approximation, both that the tangential
			** displacements stay small, and that contacts don't last long
			** compared to the motions of the particles (and where they
			** do, as in the case of grains grinding against one another,
			** they might exhibit some creep).
			*/

			/*
			** Either store this stuff in memory the previous step (current implemenation),
			** or at the very least ensure only computing vPrOld once per particle
			** If vnOld is not stored in the DEM element struct, we need to compute it somehow,
			** this would be an example of how to start:

			for (i=0;i<3;i++) {
				vPrOld[i] = p->r[i] - p->v[i]*dDelta;
				vPnrOld[i] = pn->r[i] - pn->v[i]*dDelta;
				vnOld[i] = vPnrOld[i] - vPrOld[i];
				}

			vectorNorm(vnOld); // sign: vnOld = unit vector of [vPrnOld - vPrOld], points from particle's center toward neighbor's center
			*/

			/*
			** We first consider rotation around the normal - for this,
			** one approach is to simply rotate the tangential spring
			** around the normal (we'll use the previous normal, vnOld) by
			** the average amount that the two particles rotated around
			** vnOld in the previous step. However, intuitively, there are
			** drawbacks, especially in the case of particles interacting
			** with stationary wall faces.
			*/

			/*
			** vDum below will be the average of the spins of the two particles from 1/2 step ago, midway through the previous full drift
			** step, because we want to know the amount of rotation taken around the normal over the last step.  Thus we use the spins
			** from the midpoint.
			**
			** Then we perform rotation of vShear around the normal at the midpoint between the previous step and this step by first taking
			** an average of n-hat and the old n-hat (vnOld).  vDum2 below is this normalized vector.
			*/
			  
			vectorAdd(p->w,pn->w,vDum);
			vectorAdd(pe->vnOld,n,vDum2);
			vectorNorm(vDum2);

			/* could consider checking that this is nonzero before expensive next op, but only in special cases will it be nonzero... */
			vectorRotate(pe->vShear,vDum2,0.5*dDelta*vectorDot(vDum,vDum2));
#ifdef DEM_ROTATION_DASHPOT
			vectorRotate(pe->vRoll,vDum2,0.5*dDelta*vectorDot(vDum,vDum2));
#endif /* DEM_ROTATION_DASHPOT */
			/*
			** Next, we change the orientation of the displacement vector
			** in the same way that the orientation of the normal has
			** changed over the previous step, and assume that the normal
			** did not change by more than 90 degrees.
                        */

			/*
			vectorScale(n,vectorDot(n,pe->vShear),vDum);
			vectorSub(pe->vShear,vDum,vDum);
			if ((a = vectorMagSq(vDum)) != 0.0)
				vectorScale(vDum,sqrt(vectorMagSq(pe->vShear)/a),pe->vShear);
			*/

			vectorCross(pe->vnOld,n,vDum); // axis of rotation (not normalized)
			if ((a = vectorMagSq(vDum)) != 0.0) {
				a = sqrt(a);
				vectorScale(vDum,1./a,vDum);
				vectorRotate2(pe->vShear,vDum,a,vectorDot(n,pe->vnOld)); // converting to an angle and using vectorRotate() failed when recomputing the sin and cos
#ifdef DEM_ROTATION_DASHPOT
				vectorRotate2(pe->vRoll,vDum,a,vectorDot(n,pe->vnOld));
#endif /* DEM_ROTATION_DASHPOT */
				}
			}
 
		/* set vnOld to current n for use in next step in the case that the overlap persists to the next step */
		vectorCopy(n,pe->vnOld);

		N = T = 0.0; /* normal, tangential forces (scalar) */
 
#ifdef DEM_TWOLAYERS
		if (x_inner >= 0.0)
			x = x_inner; /* if an inner boundary is defined where stiffness changes and is penetrated */
#endif

		for (i=0;i<3;i++) { /*DEBUG: optimize this function for cases when !dMuS)*/
			pe->vShear[i] += ut[i]*dDelta; /* this step's contribution to the tangential spring */
			Fn[i] = -(kn*x + F_outershell)*n[i] + Cn*un[i]; /* check to ensure that all other forces are computed first (e.g. gravity) */
			Ft[i] = kt*pe->vShear[i] + Ct*ut[i]; /*      ""      */
#ifdef DEM_ROTATION_DASHPOT
			pe->vRoll[i] += vRelRol[i]*dDelta;
#endif /* DEM_ROTATION_DASHPOT */
			N += Fn[i]*Fn[i]; /* square of normal force */
			T += Ft[i]*Ft[i]; /* square of tangential force */
                        }
#ifdef DEM_ROTATION_DASHPOT
               pe->dTwist += dRelTwi*dDelta;
#endif /* DEM_ROTATION_DASHPOT */
			
		/* Glass Beads Hack */
#ifdef _GLASS_BEADS/*DEBUG: old version - for now, skip what is just below and simply break spring in the case that static friction is exceeded (further down) */
		/*
		** hack for DEM with SPRINGS: if there is a spring connecting
		** the particles and the DEM tangential stress exceeds (1/3 ?)
		** the (normal) spring stress limit, then break the spring.
		*/

		static const double dConvert = (1.49597870e11/1.9891e30*(365.25*24*3600)*(365.25*24*3600)/(2.0*M_PI)/(2.0*M_PI)); /* Pascals --> pkdgrav units (multiply) [1.897217e-06] */
		SPRING *ps;
		double dBending_sq = 1./mu_s/mu_s; /* how many times stronger tensile strength is over shear strength squared */
		for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
			ps = &p->springs[i];
			if (ps->iOrder == pn->iOrder && ps->fZeroStrainLength > 0.) {
				a = ps->fStressLimit*dConvert*M_PI*0.5*(r1 + r2)*0.5*(r1 + r2)*smf->FO.SP.dPorosityFactor;
				if (dBending_sq*T > a*a) {
					ps->fZeroStrainLength = 0.; /* Break this link */
					printf("SPRING BROKEN VIA DEM TANGENTIAL EXTENSION (%i & %i)\n",p->iOrder,ps->iOrder);
					}
				}
			}
#endif /* _GLASS_BEADS */

		F2 = N*mu_s*mu_s; /* square of max allowed static friction */
		N = sqrt(N); /* N is now defined as the normal force (not squared) */

#ifdef DEM_ROTATION_DASHPOT
		if (!bUseContactRadius) {
			MrCrit = mu_r * beta_s * dReducedRadius * N;
			MtCrit = mu_t * beta_s * dReducedRadius * mu_s * N;
			}
		else {
			MrCrit = mu_r * beta_s * dContactRadius * N;
			MtCrit = mu_t * beta_s * dContactRadius * mu_s * N;
			}
#endif /* DEM_ROTATION_DASHPOT */

#ifdef DEM_PRESSURE
		p->dPressure += N; /*DEBUG: hack to sum total pressure on particle*/
		pn->dPressure += N; /*DEBUG: hack to sum total pressure on neighbor*/
#endif

		/*DEBUG*/
		/*
		printf("DoDEM() call = %d neighbor = %d pe->liOverlapCounter = %li x/r = %e N = %e |S|/r = %e S/|S| dot n = %e\n",count,pn->iOrder,pe->liOverlapCounter,x/r1,N,vectorMag(pe->vShear)/r2,vectorDot(pe->vShear,n)/vectorMag(pe->vShear));
		printf("S[0]/r = %e t[0] = %e  ;  F = %e  ;  Ct*Ut = %e  ;  Ft[0] = %e p->w[0] = %e\n",pe->vShear[0]/r1,t[0],sqrt(F2),Ct*Ut,Ft[0],p->w[0]);
		printf("S[1]/r = %e t[1] = %e  ;  b1 = %e  ;  Ct*Ut = %e  ;  Ft[1] = %e p->w[1] = %e\n",pe->vShear[1]/r1,t[1],b1,Ct*Ut,Ft[1],p->w[1]);
		printf("S[2]/r = %e t[2] = %e  ;  l1 = %e  ;  Ct*Ut = %e  ;  Ft[2] = %e p->w[2] = %e\n",pe->vShear[2]/r1,t[2],l1,Ct*Ut,Ft[2],p->w[2]);
		*/

		if (F2 < T) {
			if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Static friction limit exceeded: iOrder = %i(%i), neighbor = %i\n, tangential force = %.16e, limit = %.16e\n",p->iOrder,smf->pkd->idSelf,pe->iOrder,sqrt(T),sqrt(F2));

			/* Glass Beads Hack */
#ifdef GLASS_BEADS/*DEBUG*/
			/*
			** hack for DEM with SPRINGS: if there is a spring connecting the
			** particles and the static limit is exceed, then break the spring.
			*/

			SPRING *ps;
			int bSpringJustBroken;
			bSpringJustBroken = 0;
			for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
				ps = &p->springs[i];
				if (ps->iOrder == pn->iOrder && ps->fZeroStrainLength != 0.) {
					ps->fZeroStrainLength = 0.; /* Break this link */
					/*DEBUG:
					** the following line is a hack to save the distance when the spring broke. done to mitigate that fact that after a spring
					** breaks due to tangential forces, particles can find themselves in overlap and subject to significant DEM repulsive forces,
					** which is not desirable. also note that this will not get saved in *.spr files (but _probably_ in *.chk files).
					*/
					ps->fStressLimit = -d;
					bSpringJustBroken = 1;
					printf("SPRING BROKEN VIA DEM TANGENTIAL EXTENSION (%i & %i)\n",p->iOrder,ps->iOrder);
					}
				}
			if (bSpringJustBroken) continue; /* skip the DEM force calculation between this particle and its neighbor */
#endif /* GLASS_BEADS */
			/*
			** Once static friction is exceeded:
			** scale vShear to make Ft = dTangentialSpringDrag*F*Ft-hat with the new vShear-hat = Ft-hat,
			** where Ft-hat will be some linear combination of the prior vShear-hat and the current t.
			**
			** for dTangentialSpringDrag = 1, once slippage occurs, strain is reset to a value such that:
			**   tangential restoring force
			** + tangential damping force
			** + other tangential forces (e.g. grav)
			**   --------------------------------------------
			** = maximum frictional force (mu_s*Normal Force)
			**
			** Let's clarify:
			** In the case of F < Ct*Ut (tangential kinetic damping exceeds max allowed static friciton):
			**    Ft[i] = F*t[i].
			**
			** In the case of F > Ct*Ut (tangential kinetic damping does not exceed max allowed static friciton, but F < Ct*Ut + kt*vShear):
			**    Ft[i] = Ct*ut[i] + dTangentialSpringDrag*(F - Ct*Ut)*vShear-hat[i].
			**
			** Then we reset vShear to be equal to the second argument on the right side of the above equation
			** divided by kt and re-oriented in the direction of Ft, as determined by the above equation.
			** (Or simply zeroed in the former case of F < Ct*Ut.)
			**
			** Note that |Ft| will be F either if Ct*Ut > N or if dTangentialSpringDrag = 1.
			** In all other cases, some degree of slip-stick behavior will be exhibited.
			*/

			F = N*mu_s; /* max allowed static friction (not squared) */
			if (dTangentialSpringDrag != 0.0) {

				/* compute the new magnitude of vShear, but if tangential kintetic friction alone exceeds static friction, vShear will be zeroed */
				a = F - Ct*Ut;
				b = a <= 0. ? 0. : vectorMag(pe->vShear);
				a = b <= 0. ? 0. : dTangentialSpringDrag*a/(kt*b);

				/* set vShear */
				vectorScale(pe->vShear,a,pe->vShear); /* vShear for this step */
				vectorScale(t,min(F,Ct*Ut),vDum); /* the force due to tangential kinetic friction */
				vectorScale(pe->vShear,kt,vDum2);
				vectorAdd(vDum,vDum2,vDum); /* total tangential force due to both static and kinetic friction */

				mdlassert(smf->pkd->mdl,T > 0.);
				if (a)
					vectorScale(Ft,a*b/sqrt(T),pe->vShear);
				else
					vectorZero(pe->vShear);					
				vectorCopy(vDum,Ft);
				}

			else {
				vectorZero(pe->vShear);
				vectorScale(t,min(F,Ct*Ut),Ft);
				}
			}
		vectorCross(Ft,n,Ftn); /* tangential DEM induced torque per unit length */
		/* Begin Hack to save time for mu_s == 0. */
		/*	}*/
		/* End Hack to save time for mu_s == 0. */

#ifdef DEM_ROTATION_DASHPOT
		for (i=0;i<3;i++) /* rolling */
			vRft[i] = kr*pe->vRoll[i];
		dtemp1 = vectorMag(vRft);
		if (dtemp1 > MrCrit) {
			vectorScale(vRft, 1./dtemp1, vtemp1);
			for (i=0;i<3;i++) {
				vRft[i] = MrCrit * vtemp1[i];
				pe->vRoll[i] -= vRelRol[i] * dDelta;
				}
			}
		else
			for (i=0;i<3;i++)
				vRft[i] += Cr * vRelRol[i];
		dtemp1 = kq * pe->dTwist; /* twisting */
		if ( fabs(dtemp1) > MtCrit ) {
			if ( dtemp1 > 0.0 )
				for (i=0;i<3;i++)
					vTfn[i] = MtCrit * n[i];
			else
				for (i=0;i<3;i++)
					vTfn[i] = - MtCrit * n[i];
			pe->dTwist -= dRelTwi * dDelta;
			}
		else
			for (i=0;i<3;i++)
				vTfn[i] = dtemp1 * n[i] + Cq * vRelTwi[i];
#else /* Directional constant torque models */
		if (Rt != 0.) /* Rt always zero if mu_r is zero (above), so !Rt condition is better here */
			dRollingFriction = mu_r*N;
		else
			dRollingFriction = 0.;

		/* Below, we use "twisting friction" to damp out spin along the normal axis connecting the particles' centers */
		if (mu_t != 0.) {

			/* make [dTwistRate * n-hat] be the relative differential spin around normal axis: */
			vectorSub(pn->wPred,p->wPred,vDum);
			dTwistRate = vectorDot(vDum,n);

			dTwistingFriction = mu_t*N*sign(dTwistRate);

			/* alternate damping method, [Cn]*[contact area]*[differential spin]*n-hat, may want "angular" spring for this implementation */

			/* vectorScale(n,Cn*a*bb,vTfn); */
			/* above, a is contact area of interaction (this is the area inside a circle drawn by the overlap at particle surfaces, with contact point at center) */
			}

		/* moments (tangential DEM torque per unit length computed above) */
		if (Rt != 0.) {
			vectorCross(rt,n,vDum); /* unit vector giving angular direction of rolling friction torque */
			vectorScale(vDum,dRollingFriction,vRft); /* rolling friction torque per unit length */
			/* below we check to see if rolling friction impulse is greater than rolling momentum. if so, make equal */
			a = vectorMagSq(vRft);
			b = dDelta*(b1*l1*l1 + b2*l2*l2);
			b *= a*b;
			//			printf("%e %e %e %e %e %e %e\n",Rt,dDelta*b1*vRft[0],sqrt(a),dDelta*b1*dDelta*b1*a,b,vectorMag(p->w)*d,vectorMag(p->wPred)*d);
			if (Rt*Rt < b) {
				vectorScale(vRft,(vectorMag(p->w)*l1 + vectorMag(pn->w)*l2)/sqrt(b),vRft);
				//				printf("geraldine ferraro %e\n",vRft[0]);
				}
			}
		else
			vectorZero(vRft);
		if (mu_t != 0.)
			vectorScale(n,dContactRadius*dTwistingFriction,vTfn); /* twisting friction torque */
		else
			vectorZero(vTfn);
#endif /* (!)DEM_ROTATION_DASHPOT */

		/* apply forces, torques */

		if (DEBUG_TRACE(p->iOrder))
			fprintf(stderr,"applying forces and moments: iOrder = %i(%i), neighbor = %i\n",p->iOrder,smf->pkd->idSelf,pe->iOrder);

#ifdef DEM_COHESION
        /* find cohesive coefficient: cohe_coeff */
        if ( ncohe_coeff == 0 ) { /* constant cohesion */
            cohe_coeff = smf->FO.DP.dCohesiveCoeff;
            }
        else { /* particle various cohesion (denoted by iColor) */
            double dpCohe_coeff; /* cohesive coefficient of current particle */
            dpCohe_coeff = smf->FO.DP.dCohesiveCoeff; /* if no specify */
            for (i=0;i<ncohe_coeff;i++) { /* hack: may affect efficiency if ncohe_coeff is too large */
                if (p->iColor == smf->FO.DP.sCohesionVarParams.iColor[i]) {
                    dpCohe_coeff = smf->FO.DP.sCohesionVarParams.dCohesiveCoeff[i];
                    break;
                    }
                }
            if (p->iColor == pn->iColor)
                cohe_coeff = dpCohe_coeff;
            else { /* particle and its neighbour have different colors */
                double dpnCohe_coeff; /* cohesive coefficient of neighbor particle */
                dpnCohe_coeff = smf->FO.DP.dCohesiveCoeff; /* if no specify */
                for (i=0;i<ncohe_coeff;i++) { /* hack: may affect efficiency if ncohe_coeff is too large */
                    if (pn->iColor == smf->FO.DP.sCohesionVarParams.iColor[i]) {
                        dpnCohe_coeff = smf->FO.DP.sCohesionVarParams.dCohesiveCoeff[i];
                        break;
                        }
                    }
                switch (smf->FO.DP.sCohesionVarParams.nConnetCoeff) {
                    case 0: /* no cohesion between p and pn */
                        cohe_coeff = 0.0;
                        break;
                    case 1: /* use averaged cohesion of p and pn */
                    cohe_coeff = 0.5 * (dpCohe_coeff + dpnCohe_coeff);
                        break;
                    case 2: /* use maximum cohesion */
                    cohe_coeff = max(dpCohe_coeff,dpnCohe_coeff);
                        break;
                    case 3: /* use minimum cohesion */
                    cohe_coeff = min(dpCohe_coeff,dpnCohe_coeff);
                        break;
                    }
                }
            }
		if (cohe_model == 0)
			if (!bUseContactRadius)
				dtemp_cohe = cohe_coeff * 4.0 * beta_s * beta_s * dReducedRadius * dReducedRadius;
			else
				dtemp_cohe = cohe_coeff * M_PI * beta_s * beta_s * dContactRadius*dContactRadius;
		else
			dtemp_cohe = cohe_coeff * 2.0*M_PI * r1sq * r2sq / (r1sq + r2sq); /* Sanchez's method */
		Fn[0] +=  dtemp_cohe*n[0];
		Fn[1] +=  dtemp_cohe*n[1];
		Fn[2] +=  dtemp_cohe*n[2];
#endif /* DEM_COHESION */

#ifdef DEM_DIAG
		pe->vFn[0] = Fn[0];
		pe->vFn[1] = Fn[1];
		pe->vFn[2] = Fn[2];
		pe->vFt[0] = Ft[0];
		pe->vFt[1] = Ft[1];
		pe->vFt[2] = Ft[2];
#endif /* DEM_DIAG */

#ifndef AGGS // JCM 10.10.20
#ifdef WALLS
		for (i=0;i<3;i++) {
			if (!PARTICLE_STUCK(p)) {
				p->dDeltaAccel[i] += a1*(Fn[i] + Ft[i]);
#ifdef DEM_ROTATION_DASHPOT
                p->wDot[i] -= b1*(l1*Ftn[i] - vRft[i] - vTfn[i]);
#else
				p->wDot[i] -= b1*(l1*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); /* might try a factor of 10xdMu_r to make Contact Radius commensurate with legacy RF */
#endif /* !DEM_ROTATION_DASHPOT */
				}
			if (!PARTICLE_STUCK(pn)) {
				pn->dDeltaAccel[i] -= a2*(Fn[i] + Ft[i]);
#ifdef DEM_ROTATION_DASHPOT
				pn->wDot[i] -= b2*(l2*Ftn[i] + vRft[i] + vTfn[i]);
#else
				pn->wDot[i] -= b2*(l2*Ftn[i] + dRollLeverNeighbor*vRft[i] + vTfn[i]);
#endif /* !DEM_ROTATION_DASHPOT */
				}
			}
		/*DEBUG!{
			double FnMag,FtMag;
			FnMag = a1*sqrt(Fn[0]*Fn[0] + Fn[1]*Fn[1] + Fn[2]*Fn[2]);
			FtMag = a1*sqrt(Ft[0]*Ft[0] + Ft[1]*Ft[1] + Ft[2]*Ft[2]);
			if (FnMag > 2000 || FtMag > 2000)
				printf("SOUND %e %i %i %g %g\n",
					   smf->dTime,p->iOrder,pn->iOrder,FnMag,FtMag);
					   }*/
#else /* WALLS */
		for (i=0;i<3;i++) {
			p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
			pn->dDeltaAccel[i] -= (Fn[i] + Ft[i])*a2;
#ifdef DEM_ROTATION_DASHPOT
			p->wDot[i] -= b1*(l1*Ftn[i] - vRft[i] - vTfn[i]);
			pn->wDot[i] -= b2*(l2*Ftn[i] + vRft[i] + vTfn[i]);
#else
			p->wDot[i] -= b1*(l1*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); /* might try a factor of 10xdMu_r to make Contact Radius commensurate with legacy RF */
			pn->wDot[i] -= b2*(l2*Ftn[i] + dRollLeverNeighbor*vRft[i] + vTfn[i]);
#endif /* !DEM_ROTATION_DASHPOT */
			}
#endif /* !WALLS */


#else /* AGGS */
#ifdef WALLS
        if (!PARTICLE_STUCK(p)) {
            if (p->iOrgIdx < 0) {
                vectorScale(n,l1,vAggTemp); /* Vector from particle center to contact point. */
                vectorSub(p->r_unwrapped,p->agg_COM,vBodyPos); /* Location of particle center relative to agg COM. JVD */
                vectorAdd(vBodyPos,vAggTemp,vAggLever); /* Location of contact point in body frame. JVD - r_body*/
                vectorAdd(Fn,Ft,vAggForces); /* Collect tangential and normal forces. */
                vectorScale(vAggForces,a1,vAggTemp2); /* Get DEM acceleration */
                vectorCross(vAggLever,vAggForces,vAggTorques); /* Calculate torques from norm. and tan. forces. */
                vectorAdd(vAggTemp2,p->DEM_aggaccels,p->DEM_aggaccels); /* Add new accels from DoDEM. */
                vectorAdd(vAggTorques,p->DEM_aggtorques,p->DEM_aggtorques); /* Ditto for torques. */
#ifdef DEM_ROTATION_DASHPOT
                /* ROLLING */
                vectorCross(vRft,vAggTemp,vAggTemp2); /* Vector parallel to rolling friction force. */
                if (vectorMag(vRft) == 0 || vectorMag(vAggTemp2) == 0) {
#if (INTERNAL_WARNINGS != 0)
                    static int nWarn = 1;
                    if (mu_r != 0 && (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)) {
                        fprintf(stderr,"Particle: %d: WARNING: Rolling torque too small to apply to aggregate. Setting rolling torque to zero.\n", p->iOrder);
                    }
                    ++nWarn;
#endif /* INTERNAL_WARNINGS */
                    vectorZero(vRollForces); /* In fact, we just skip the rolling torque calculation in this case. We zero out vRollForces as a marker for pn (see below). */
                } else {
                    vectorNorm(vAggTemp2); /* Unit vector in direction of rolling force. */
                    vectorScale(vAggTemp2,vectorMag(vRft)/vectorMag(vAggTemp),vRollForces); /* The net rolling force. */
                    vectorCross(vAggLever,vRollForces,vAggTorques); /* Net rolling torque on the agg. */
                    vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
                }
                /* TWISTING */
                vectorCopy(vTfn,vAggTorques);
                vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
#else

#endif /* !DEM_ROTATION_DASHPOT */
            } else {
                for (i=0;i<3;i++) { 
                    p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
#ifdef DEM_ROTATION_DASHPOT
                    p->wDot[i] -= b1*(l1*Ftn[i] - vRft[i] - vTfn[i]);
#else
                    p->wDot[i] -= b1*(l1*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); 
#endif /* !DEM_ROTATION_DASHPOT */
                }
            }
        }

        if (!PARTICLE_STUCK(pn)) {
            if (pn->iOrgIdx < 0) {
                vectorScale(n,-l2,vAggTemp); /* Repeat above calculation for neighbor particle, if in an agg. */
                vectorSub(pn->r_unwrapped,pn->agg_COM,vBodyPos); //JVD - New vector math
                vectorAdd(vBodyPos,vAggTemp,vAggLever);
                vectorScale(vAggForces,-a2,vAggTemp2);
                vectorScale(vAggForces,-1.0,vAggForces); /* Force should be equal and *opposite* */
                vectorCross(vAggLever,vAggForces,vAggTorques);
                vectorAdd(vAggTemp2,pn->DEM_aggaccels,pn->DEM_aggaccels);
                vectorAdd(vAggTorques,pn->DEM_aggtorques,pn->DEM_aggtorques);
#ifdef DEM_ROTATION_DASHPOT
                /* ROLLING */
                if (vectorMag(vRollForces) == 0) {
#if (INTERNAL_WARNINGS != 0)
                    static int nWarn = 1;
                    if (mu_r != 0 && (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)) {
                        fprintf(stderr,"Particle: %d: WARNING: Rolling torque too small to apply to aggregate. Setting rolling torque to zero.\n", pn->iOrder);
                    }
                    ++nWarn;
#endif /* INTERNAL_WARNINGS */
                } else {
                vectorCross(vAggLever,vRollForces,vAggTorques); /* Net rolling torque on the neighbor agg */
                vectorAdd(pn->DEM_aggtorques,vAggTorques,pn->DEM_aggtorques);
                }
                /* TWISTING */
                vectorCopy(vTfn,vAggTorques);
                vectorAdd(pn->DEM_aggtorques,vAggTorques,pn->DEM_aggtorques);
#else

#endif /* !DEM_ROTATION_DASHPOT */
            } else {
                for (i=0;i<3;i++) {
                    pn->dDeltaAccel[i] -= (Fn[i] + Ft[i])*a2;
#ifdef DEM_ROTATION_DASHPOT
                    pn->wDot[i] -= b2*(l2*Ftn[i] + vRft[i] + vTfn[i]);
#else
                    pn->wDot[i] -= b2*(l2*Ftn[i] + dRollLeverNeighbor*vRft[i] + vTfn[i]);
#endif /* !DEM_ROTATION_DASHPOT */
                }
            }
        }

#else /* !WALLS */
        if (p->iOrgIdx < 0) {
            vectorScale(n,l1,vAggTemp); /* Vector from particle center to contact point. JVD */
            vectorSub(p->r_unwrapped,p->agg_COM,vBodyPos); //JVD - New vector math
            vectorAdd(vBodyPos,vAggTemp,vAggLever);
            vectorAdd(Fn,Ft,vAggForces); /* Collect tangential and normal forces. */
            vectorScale(vAggForces,a1,vAggTemp2); /* Get DEM acceleration */
            vectorCross(vAggLever,vAggForces,vAggTorques); /* Calculate torques from norm. and tan. forces. */
            vectorAdd(vAggTemp2,p->DEM_aggaccels,p->DEM_aggaccels); /* Add new accels from DoDEM. */
            vectorAdd(vAggTorques,p->DEM_aggtorques,p->DEM_aggtorques); /* Ditto for torques. */
#ifdef DEM_ROTATION_DASHPOT
            /* ROLLING */
            vectorCross(vRft,vAggTemp,vAggTemp2); /* Vector parallel to rolling friction force. */
            if (vectorMag(vRft) == 0 || vectorMag(vAggTemp2) == 0) {
#if (INTERNAL_WARNINGS != 0)
                static int nWarn = 1;
                if (mu_r != 0 && (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)) {
                    fprintf(stderr,"Particle: %d: WARNING: Rolling torque too small to apply to aggregate. Setting rolling torque to zero.\n", p->iOrder);
                }
                ++nWarn;
#endif /* INTERNAL_WARNINGS */
                vectorZero(vRollForces); /* In fact, we just skip the rolling torque calculation in this case. We zero out vRollForces as a marker for pn (see below). */
            } else {
                vectorNorm(vAggTemp2); /* Unit vector in direction of rolling force. */
                vectorScale(vAggTemp2,vectorMag(vRft)/vectorMag(vAggTemp),vRollForces); /* The net rolling force. */
                vectorCross(vAggLever,vRollForces,vAggTorques); /* Net rolling torque on the agg. */
                vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
            }
            /* TWISTING */
            vectorCopy(vTfn,vAggTorques);
            vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
#else

#endif /* !DEM_ROTATION_DASHPOT */
        } else {
            for (i=0;i<3;i++) {
                p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
#ifdef DEM_ROTATION_DASHPOT
                p->wDot[i] -= b1*(l1*Ftn[i] - vRft[i] - vTfn[i]);
#else
				p->wDot[i] -= b1*(l1*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); 
#endif /* !DEM_ROTATION_DASHPOT */
            }
        }

        if (pn->iOrgIdx < 0) {
            vectorScale(n,-l2,vAggTemp); /* Repeat above calculation for neighbor particle, if in an agg. */
            vectorSub(pn->r_unwrapped,pn->agg_COM,vBodyPos); //JVD - New vector math
            vectorAdd(vBodyPos,vAggTemp,vAggLever);
            vectorScale(vAggForces,-a2,vAggTemp2);
            vectorScale(vAggForces,-1.0,vAggForces); /* Force should be equal and *opposite* */
            vectorCross(vAggLever,vAggForces,vAggTorques);
            vectorAdd(vAggTemp2,pn->DEM_aggaccels,pn->DEM_aggaccels);
            vectorAdd(vAggTorques,pn->DEM_aggtorques,pn->DEM_aggtorques);
#ifdef DEM_ROTATION_DASHPOT
            /* ROLLING */
            if (vectorMag(vRollForces) == 0) {
#if (INTERNAL_WARNINGS != 0)
                static int nWarn = 1;
                if (mu_r != 0 && (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)) {
                    fprintf(stderr,"Particle: %d: WARNING: Rolling torque too small to apply to aggregate. Setting rolling torque to zero.\n", pn->iOrder);
                }
                ++nWarn;
#endif /* INTERNAL_WARNINGS */
            } else {
                vectorCross(vAggLever,vRollForces,vAggTorques); /* Net rolling torque on the neighbor agg */
                vectorAdd(pn->DEM_aggtorques,vAggTorques,pn->DEM_aggtorques);
            }
            /* TWISTING */
            vectorCopy(vTfn,vAggTorques);
            vectorAdd(pn->DEM_aggtorques,vAggTorques,pn->DEM_aggtorques);
#else

#endif /* !DEM_ROTATION_DASHPOT */
        } else {
            for (i=0;i<3;i++) {
                pn->dDeltaAccel[i] -= (Fn[i] + Ft[i])*a2;
#ifdef DEM_ROTATION_DASHPOT
				pn->wDot[i] -= b2*(l2*Ftn[i] + vRft[i] + vTfn[i]);
#else
				pn->wDot[i] -= b2*(l2*Ftn[i] + dRollLeverNeighbor*vRft[i] + vTfn[i]);
#endif /* !DEM_ROTATION_DASHPOT */
            }
        }
#endif /* !WALLS */
#endif /* AGGS */
#ifdef GLASS_BEADS /* for safety */
		r1 = RADIUS(p);
		r1sq = r1*r1;
#endif
		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"particle-particle forces and moments applied: iOrder = %i(%i), neighbor = %i, current tally of particle force: %.16e, moment: %.16e, neighbor force: %.16e, moment: %.16e\n",p->iOrder,smf->pkd->idSelf,pe->iOrder,vectorMag(p->dDeltaAccel),vectorMag(p->wDot),vectorMag(pn->dDeltaAccel),vectorMag(pn->wDot));
		}
/******************* Impactor Acceleration profile *************/
/*Note:
* Tracking of acceleration didn't work because they aren't updated till later. Only p->dDeltaAccel is updated here*/
 /*
  double F_total_track[3];
  FILE *f_Impactor;
  if (p->iOrder == 0) {
    f_Impactor = fopen("ImpactorData.out","a");
    for (i=0;i<3;i++)
      F_total_track[i] = p->dDeltaAccel[i];
    fprintf(f_Impactor," %.8g %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e\n",smf->dTime, p->r[0], p->r[1], p->r[2], p->v[0], p->v[1], p->v[2], F_total_track[0], F_total_track[1], F_total_track[2]);
    fclose(f_Impactor);
    }
 */
/******************* Impactor Acceleration profile *************/


#ifdef WALLS
	if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"beginning particle-walls calculations: iOrder = %i(%i)\n",p->iOrder,smf->pkd->idSelf);
	if (PARTICLE_STUCK(p))
		return; /* skip walls portion of DoDEM() if particle is stuck to any wall */

#ifdef DEM_WALLS_REACT
  /* reset tally of force from particle onto non-infinitely massive wall(s assemblege) */
  for (i=0;i<MAX_NUM_WALL_ASSEMBLIES;i++) {
    p->dXForceOnWalls[i] = 0.0;
    p->dYForceOnWalls[i] = 0.0;
    p->dZForceOnWalls[i] = 0.0;
    p->dXTorqueOnWalls[i] = 0.0;
    p->dYTorqueOnWalls[i] = 0.0;
    p->dZTorqueOnWalls[i] = 0.0;
    }
  Vector rArm,rFnT,rFtT,rWrist;
  //rArm: vector that will be distance between particle and TAG Wrist
  //rFnT: normal component of torque'ing rxF
  //rFtT: tangential component of torque'ing rxF
  //rWrist: position of TAG Wrist: Joing between TAG head and TAG arm
  vectorZero(rWrist);
#endif /* DEM_WALLS_REACT */

  /* q: position vector of wall intersection with line connecting particle to wall origin */

	DEM_ELEMENT *we;
	WALL_PARAMS *WP = &smf->WP;
	WALL *w;
	WALL_DATA *wd;
	double a3,dRadius,dRadNarrow,dLength,dTaper,dOpenAngle,dAngle;
    double dTmp = 0;
	Vector O,q,vDum3,v1,v2,vVinyl;
	Matrix A;
	static double CnPreFacWalls[MAX_NUM_WALLS],CtPreFacWalls[MAX_NUM_WALLS];

#ifdef DEM_TWOLAYERS
	static double CnOuterPreFacWalls[MAX_NUM_WALLS],CtOuterPreFacWalls[MAX_NUM_WALLS],CnInnerPreFacWalls[MAX_NUM_WALLS],CtInnerPreFacWalls[MAX_NUM_WALLS];
#endif

	/* check to see if particle is overlapping with wall boundaries */

	/* loop over walls */
	for (j=0;j<WP->nWalls;j++) {
		w = &WP->pWalls[j]; /* wall */
		wd = &w->wd; /* wall's data struct */

		/*
		if (p->iOrder == 0) printf("nWalls: %d\ni: %d\niWallID: %d\niType: %d\nvOrigin: %e, %e, %e\nvOrient: %e, %e, %e\nvVertex1: %e, %e, %e\nvVertex2: %e, %e, %e\nvVel: %e, %e, %e\ndStep: %e\ndOscAmp: %e\ndOscFreq: %e\ndOscPhase: %e\nvOscVec: %e, %e, %e\ndRadius: %e\ndHoleRadius: %e\ndLength: %e\ndTaper: %e\ndOpenangle: %e\ndAngSpeed: %e\ndEpsN: %e\ndEpsT: %e\niColor: %d\ndTrans: %e\nvOscVel: %e, %e, %e\nvTravel: %e, %e, %e\nvTotVel: %e, %e, %e\n",WP->nWalls,j,w->iWallID,w->wd.iType,w->wd.vOrigin[0],w->wd.vOrigin[1],w->wd.vOrigin[2],w->wd.vOrient[0],w->wd.vOrient[1],w->wd.vOrient[2],w->wd.vVertex1[0],w->wd.vVertex1[1],w->wd.vVertex1[2],w->wd.vVertex2[0],w->wd.vVertex2[1],w->wd.vVertex2[2],w->wd.vVel[0],w->wd.vVel[1],w->wd.vVel[2],w->wd.dStep,w->wd.dOscAmp,w->wd.dOscFreq,w->wd.dOscPhase,w->wd.vOscVec[0],w->wd.vOscVec[1],w->wd.vOscVec[2],w->wd.dRadius,w->wd.dHoleRadius,w->wd.dLength,w->wd.dTaper,w->wd.dOpenAngle,w->wd.dAngSpeed,w->wd.dEpsN,w->wd.dEpsT,w->wd.iColor,w->wd.dTrans,w->vOscVel[0],w->vOscVel[1],w->vOscVel[2],w->vTravel[0],w->vTravel[1],w->vTravel[2],w->vTotVel[0],w->vTotVel[1],w->vTotVel[2]);
		wd->vVel[2] += 1.;
		*/

		if (wd->iType != WallPlane && wd->iType != WallDisk && wd->iType != WallRectangle &&
		    wd->iType != WallTriangle && wd->iType != WallCylinderInfinite &&
		    wd->iType != WallCylinderFinite && wd->iType != WallShell)
			mdlassert(smf->pkd->mdl,0);

		if (wd->dEpsN == PARTICLE_OVERRIDE_EPSN)
			dEpsN = smf->FO.DP.dEpsN;
		else
			dEpsN = wd->dEpsN;
		if (wd->dEpsT == PARTICLE_OVERRIDE_EPST)
			dEpsT = smf->FO.DP.dEpsT;
		else
			dEpsT = wd->dEpsT;
		if (wd->dKn == PARTICLE_OVERRIDE_KN)
			kn = smf->FO.DP.dKn;
		else
			kn = wd->dKn;
		if (wd->dKt == PARTICLE_OVERRIDE_KT)
			kt = smf->FO.DP.dKt;
		else
			kt = wd->dKt;
		if (wd->dMuS == PARTICLE_OVERRIDE_MUS)
			mu_s = smf->FO.DP.dMuS;
		else
			mu_s = wd->dMuS;
		if (wd->dMuR == PARTICLE_OVERRIDE_MUR)
			mu_r = smf->FO.DP.dMuR;
		else
			mu_r = wd->dMuR;
		if (wd->dMuT == PARTICLE_OVERRIDE_MUT)
			mu_t = smf->FO.DP.dMuT;
		else
			mu_t = wd->dMuT;
#ifdef DEM_ROTATION_DASHPOT
		if (wd->dBetaS == PARTICLE_OVERRIDE_BETAS)
			beta_s = smf->FO.DP.dBetaS;
		else
			beta_s = wd->dBetaS;
#endif
#ifdef DEM_COHESION
		if (wd->iCohesionModel == PARTICLE_OVERRIDE_COMO)
			cohe_model = smf->FO.DP.iCohesionModel;
		else
			cohe_model = wd->iCohesionModel;
		if (wd->dCohesiveCoeff == PARTICLE_OVERRIDE_COCO)
			cohe_coeff = smf->FO.DP.dCohesiveCoeff;
		else
			cohe_coeff = wd->dCohesiveCoeff;
#endif
#ifdef DEM_TWOLAYERS
		if (wd->dInnerOverlapBndry == PARTICLE_OVERRIDE_IOB)
			dInBndFrac = smf->FO.DP.dInnerOverlapBndry;
		else
			dInBndFrac = wd->dInnerOverlapBndry;
		if (wd->dKnOuter == PARTICLE_OVERRIDE_KNOUTER)
			kn_outer = smf->FO.DP.dKnOuter;
		else
			kn_outer = wd->dKnOuter;
		if (wd->dKtOuter == PARTICLE_OVERRIDE_KTOUTER)
			kt_outer = smf->FO.DP.dKtOuter;
		else
			kt_outer = wd->dKtOuter;
		if (wd->dKnInner == PARTICLE_OVERRIDE_KNINNER)
			kn_inner = smf->FO.DP.dKnInner;
		else
			kn_inner = wd->dKnInner;
		if (wd->dKtInner == PARTICLE_OVERRIDE_KTINNER)
			kt_inner = smf->FO.DP.dKtInner;
		else
			kt_inner = wd->dKtInner;
#endif /* DEM_TWOLAYERS */
		/* dTangentialSpringDrag = 0; *//* set to zero above (disabled) */

#ifdef DEM_OVERLAP_OUTPUT /* this will only work if nParticles is small, otherwise it should seg fault b/c it can't allocate memory */
		static double x_maxWall[nParticles][MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS];
		if (liStep == 1)
			x_maxWall[p->iOrder][j] = 0;
#endif /* DEM_OVERLAP_OUTPUT */

		iOverlap = -1; /* for safety */
		F_outershell = 0.;

		vectorCopy(wd->vOrigin,O);
		vectorCopy(wd->vOrient,n);
		bApplyForce = bFoundNbr = 0;

		/* determine bApplyForce */
		if (wd->iType == WallPlane || wd->iType == WallDisk || wd->iType == WallRectangle ||
		    wd->iType == WallTriangle || wd->iType == WallCylinderFinite) {

			/*
			** find distance from particle center to plane.
			**
			** strategy:
			** find eqn of plane in form: (A*x + B*y + C*z = D).
			** eqn of line through particle center and closest point to plane is r + d*n,
			** where r is the particle center, n is the normal, and d*n is the distance.
			**
			** Combine and solve for d.
			*/

			d = vectorDot(p->r,n) - vectorDot(O,n); /* particle center lies a distance d off plane */ /* edit here to save time by subtracting vectors first then perform single dot product */
			if (fabs(d) >= r1 && wd->iType != WallCylinderFinite)
				bApplyForce = 0;
			else {
				vectorScale(n,d,vDum);
				vectorSub(p->r,vDum,q); /* q is now projection of particle center onto plane */
				switch (wd->iType) {
				case WallCylinderFinite: /* move this entire case outside this switch?? */

					if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with finite cylinder: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

					/*
					** note: For particle radius larger than cylinder radius that lies inside cylinder,
					** only the force from the point on the cylinder closest to center of particle
					** will be felt. If needed, a correction might be added by scaling down this force
					** to account for the counterbalancing net force from the opposite direction. For
					** particles centered inside the cylinder and directly on the cylinder axis, no
					** no force is felt from the sides of the cylinder. On the axis inside a tapered
					** cylinder, a net force is applied "down" the axis toward the wide end via the
					** creation of a phantom point (described below). On the axis inside a uniform
					** no force is applied whatsoever.
					*/

					dRadius = wd->dRadius;
					dLength = wd->dLength;
					dTaper = wd->dTaper;

					dRadNarrow = dRadius*(1. - dTaper); /* cylinder's narrow-end ring (for taper) */
					vectorSub(q,O,vDum); /* points from cylinder axis to particle center parallel to orthogonal plane */
					a2 = vectorMagSq(vDum);
					a1 = dRadius + r1;
					a3 = max(dRadNarrow - r1,0.0);

					/*
					** start with a quick overinclusive proximity check:
					** does any part of the particle lie between
					** cylinder of radius dRadius*[1-dTaper] and
					** cylinder of radius dRadius (both of length dLength)?
					** (if not --> !bApplyForce)
					*/

					if (fabs(d + d) <= dLength + r1 + r1 && a2 <= a1*a1 && a2 >= a3*a3) {
						if (a2 != 0.0) {
							vectorScale(vDum,1./sqrt(a2),t); /* unit vector from cyl axis to p perp to cyl axis */
							vectorScale(t,dRadNarrow,vDum);
							vectorScale(n,0.5*dLength,vDum3);
							vectorAdd(vDum,vDum3,vDum); /* vector from O to point on narrow rim closest to p */
							vectorScale(t,dRadius,vDum2);
							vectorSub(vDum2,vDum3,vDum2); /* vector from O to point on wide rim closest to p */
							vectorSub(p->r,O,vDum3); /* vector pointing from O to p */
							vectorGetClosestPointOnSegment(vDum,vDum2,vDum3,vDum);
							vectorAdd(O,vDum,q); /* point on cylinder closest to p */
							vectorSub(p->r,q,vDum);
							bApplyForce = demCheckOverlapPoint(vDum,r1sq);
							}
						else { /* special case: p lies on cylinder axis */

							/*
							** (for the time being at least) some cleanliness is sacrificed for speed,
							** since the following is mostly 2D:
							**
							**
							** If there is an overlap in this case, the overlap occurs over a ring that
							** is symmetric around the cylinder axis, and therefore the overlap is
							** radially symmetric. The net force should point along the axis, so the
							** strategy used here is to create a "phantom" overlap point on the axis that
							** penetrates the particle by the same amount as it is penetrated by the ring.
							*/

							/* a vector pointing from a point on the narrow rim to p (cyl coords): */
							vectorSet(vDum2,0.,-dRadNarrow,d - 0.5*dLength);

							/* a vector pointing from a point on the wide rim to p (cyl coords): */
							vectorSet(vDum,0.,-dRadius,d + 0.5*dLength);
							
							if (d + d > dLength) {
								if (r1 >= dRadNarrow) { /* skips the next calculation if not needed */
									b = vDum2[1]*vDum2[1] + vDum2[2]*vDum2[2]; /* vectorMagSq(vDum2) */
									if ((bApplyForce = (b <= r1sq))) { /*DEBUG: modularize this and next few lines...*/
										vectorScale(n,d - sqrt(b),vDum);
										vectorAdd(O,vDum,q);
										}
									}
								}

							/*
							** only *tapered* cylinders will push out a particle whose center is
							** 'inside' the cylinder and centered directly on the axis, so we'll
							** ignore particles inside a non-tapered cylinder. Strictly speaking,
							** tapered and non-tapered cylinders should provide frictional
							** stability to particles with radii larger than the cylinders
							** that lie inside the cylinder and directly on the axis. If this
							** proves to be an issue, it can be dealt with as a special case.
							*/
							else if (dTaper || (d + d < -dLength)) {
								vectorZero(vDum3);
								vectorGetClosestPointOnSegment(vDum,vDum2,vDum3,vDum);
								/*
								** vDum is in fact now a *vector* that points from the closest
								** point on the cylinder to the particle's center, in cylinder
								** coordinates.
								*/

								b = vectorMagSq(vDum);
								if ((bApplyForce = (b <= r1sq))) {
									vectorScale(n,d + sqrt(b),vDum);
									vectorAdd(O,vDum,q);
									}
								}
							}
						}
					break;

				case WallPlane:

					if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"overlap with infinite plane: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

					bApplyForce = 1;
					break;

				case WallDisk:

					if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with disk: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

					/*
					** find distance to disk.
					**
					** strategy:
					** compare | q - wall_origin | to disk radius, if less q is on disk.
					** if more, find the point on rim closest to q, then
					** compare | p->r - this point | to particle radius, if less particle is touching rim.
					**
					** to find such a point on the rim closest to q, add to origin:
					** (disk radius)*(normalized vector pointing from origin to q).
					*/

					vectorSub(q,O,vDum); /* radial component of vector from origin to particle */
					b = vectorMagSq(vDum);
					a1 = wd->dRadius*wd->dRadius;
					a2 = wd->dHoleRadius*wd->dHoleRadius;
					bOutcome = (b >= a2) + (b >= a2) + (b <= a1);

					/* hole radius bigger than disk radius */
					if (!bOutcome) {
						mdlassert(smf->pkd->mdl,0);
						}

					/* contact is flush on plane of disk */
					else if (bOutcome == 3) {
						bApplyForce = 1;
						}

					/* check for contact with outer rim of disk */
					else if (bOutcome == 2) {
						vectorNorm(vDum); /* could save time by scaling once, not thrice */
						vectorScale(vDum,wd->dRadius,vDum);
						vectorAdd(O,vDum,q); /* q is now point on disk rim closest to particle */
						vectorSub(p->r,q,vDum);
						d = vectorMagSq(vDum);
						if (d <= r1sq) {
							bApplyForce = 1;
							}
						}

					/* check for contact with inner rim of disk */
					else if (bOutcome == 1) {
						if (b != 0.) {
							vectorNorm(vDum); /* could save time by scaling once, not thrice */
							vectorScale(vDum,wd->dHoleRadius,vDum);
							vectorAdd(O,vDum,q); /* q is now point on disk rim closest to particle */
							vectorSub(p->r,q,vDum);
							d = vectorMagSq(vDum);
							if (d <= r1sq) {
								bApplyForce = 1;
								}
							}
						else { /* particle center on normal axis of disk with hole */
							if (r1 >= wd->dHoleRadius) { /* skips the next calculation if not needed */
								b = wd->dHoleRadius*wd->dHoleRadius + d*d; /* d is still distance off plane perp to plane */
								if ((bApplyForce = (b <= r1sq))) {
									vectorScale(n,d - sign(d)*sqrt(b),vDum);
									vectorAdd(O,vDum,q);
									}
								}
							}
						}
					break;

				case WallRectangle:

					if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with rectangle: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

					/*
					** find distance to rectangle.
					**
					** strategy:
					** first see if q falls inside rectangle by solving for a & b.
					** if both have values between 0 and 1 (inclusive), q is on rectangle.
					** otherwise, set to closest point on rectangle, then check proximity.
					**
					** on which edge or corner this point lies depends upon the values of a & b.
					*/

					vectorCopy(wd->vVertex1,v1);
					vectorCopy(wd->vVertex2,v2);
					vectorSub(q,O,vDum);

					for (i=0;i<3;i++)
						vectorSet(A[i],v1[i],v2[i],n[i]);

					matrixInvert(A);
					a = vectorDot(A[0],vDum);
					b = vectorDot(A[1],vDum); /* note: vectorDot(A[2],vDum) should be zero */

					if (a >= 0. && a <= 1. && b >= 0. && b <= 1.)
						bApplyForce = 1;
					else { /* eight cases: four edges, four corners */

						/* segment joining points O and (O + v2) */
						if (a < 0. && b >= 0. && b <= 1.) {
							vectorScale(v2,b,vDum);
							vectorAdd(O,vDum,q);
							}

						/* segment joining points (O + v1) and (O + v1 + v2) */
						else if (a > 1. && b >= 0. && b <= 1.) {
							vectorScale(v2,b,vDum);
							vectorAdd(O,vDum,vDum);
							vectorAdd(v1,vDum,q);
							}

						/* segment joining points O and (O + v1) */
						else if (b < 0. && a >= 0. && a <= 1.) {
							vectorScale(v1,a,vDum);
							vectorAdd(O,vDum,q);
							}

						/* segment joining points (O + v2) and (O + v1 + v2) */
						else if (b > 1. && a >= 0. && a <= 1.) {
							vectorScale(v1,a,vDum);
							vectorAdd(O,vDum,vDum);
							vectorAdd(v2,vDum,q);
							}

						/* corner at point of origin */
						else if (a < 0. && b < 0.)
							vectorCopy(O,q);

						/* corner at point (O + v1) */
						else if (a > 1. && b < 0.)
							vectorAdd(O,v1,q);

						/* corner at point (O + v2) */
						else if (a < 0. && b > 1.) {
							vectorAdd(O,v2,q);
							/* printf("a=%e b=%e\n",a,b); */
							}

						/* corner at point (O + v1 + v2) */
						else if (a > 1. && b > 1.) {
							vectorAdd(O,v1,vDum);
							vectorAdd(vDum,v2,q);
							}

						else
							mdlassert(smf->pkd->mdl,0); /* keep things honest */

						vectorSub(p->r,q,vDum);
						bApplyForce = demCheckOverlapPoint(vDum,r1sq);
						}
					break;

				case WallTriangle:

					if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with triangle: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

					/*
					** find distance to triangle.
					**
					** strategy:
					** first see if q falls inside triangle by solving for a & b.
					** if both are non-negative and they sum to a value no greater than unity, q is on triangle.
					** otherwise, set q to closest point on triangle, then check proximity.
					**
					** on which edge or corner this point lies depends upon the values of a & b.
					*/

					vectorCopy(wd->vVertex1,v1);
					vectorCopy(wd->vVertex2,v2);
					vectorSub(q,O,vDum);

					for (i=0;i<3;i++)
						vectorSet(A[i],v1[i],v2[i],n[i]);

					matrixInvert(A);
					a = vectorDot(A[0],vDum);
					b = vectorDot(A[1],vDum); /* note: vectorDot(A[2],vDum) should be zero */

					if (a >= 0. && b >= 0. && (a + b) <= 1.)
						bApplyForce = 1;
					else { /* six cases: three edges, three corners (but we combine two corners and an edge) */

						/* segment joining points O and (O + v2) */
						if (a < 0. && b >= 0. && b <= 1.) {
							vectorScale(v2,b,vDum);
							vectorAdd(O,vDum,q);
							}

						/* segment joining points O and (O + v1) */
						else if (b < 0. && a >= 0. && a <= 1.) {
							vectorScale(v1,a,vDum);
							vectorAdd(O,vDum,q);
							}

						/* corner at point of origin */
						else if (a < 0. && b < 0.)
							vectorCopy(O,q);

						/* segment joining points (O + v1) and (O + v2), including points v1 and v2 */
						/* DEBUG: separate these three cases for increased speed */
						else {
							vectorSub(q,O,vDum);
							vectorGetClosestPointOnSegment(v1,v2,vDum,vDum);
							vectorAdd(O,vDum,q);
							}

						vectorSub(p->r,q,vDum);
						bApplyForce = demCheckOverlapPoint(vDum,r1sq);
						}
					break;

				default:
					mdlassert(smf->pkd->mdl,0);
					}
				}
			}

		else if (wd->iType == WallCylinderInfinite) {

			if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with infinite cylinder: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

			vectorSub(O,p->r,vDum);
			vectorCross(n,vDum,vDum2);
			d = vectorMagSq(vDum2); /* sqare of distance from particle center to cylinder axis */
			a1 = wd->dRadius + r1;
			a2 = wd->dRadius - r1;
			a3 = max(a2,0.0);
			if (d != 0.0 && d <= a1*a1 && d >= a3*a3) {
				a = -vectorDot(n,vDum); /* scalar projection of segment (p->r - O) onto cylinder axis */
				vectorScale(n,a,vDum); /* vector pointing from O to point on cylinder axis closest to particle center */
				vectorSub(p->r,vDum,vDum2);
				vectorSub(vDum2,O,vDum2); /* vector pointing from cylinder axis to particle center */
				d = sqrt(d); /* equivalent to magnitude of current vDum2 vector */
				vectorScale(vDum2,wd->dRadius/d,vDum2);
				vectorAdd(O,vDum,vDum);
				vectorAdd(vDum,vDum2,q); /* point on cylinder closest to center of particle */
				bApplyForce = 1;
				}
			/*
			** note: For particle radius larger than cylinder radius
			** that lies inside cylinder, only the force from the
			** point on the cylinder closest to center of particle
			** will be felt. If this case for some reason has a need
			** to be explored, consider scaling down this force to
			** account for counterbalancing net force in the opposite
			** direction. Also, relatedly, for such a particle
			** centered directly on the cylinder axis, rather than the
			** particle being "pinned" to its spot inside the
			** cylinder, no force is felt. Both of these conditions
			** also apply to finite cylinders of uniform width. In the
			** case of the tapered finite cylinder, similar scenarios
			** are treated with more care (discussed in that section
			** of the code).
			**
			** again, if such a case proves salient somehow, revamp
			** accordingly.
			*/
			}

		else if (wd->iType == WallShell) {

			if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"checking for overlap with shell: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,w->iWallID);

			/*
			** note: the closer the particle's cross-sectional radius
			** (along the plane of contact perpendicular to the
			** shell's axis of symmetry) is to the shell's
			** cross-sectional radius on this same plane, the less
			** accurate the treatment stricly is. For example,
			** consider the difference between an impact on the inside
			** vs. an impact on the outside of the shell - a
			** penetration of 'x' from the inside implies a greater
			** degree of overlap than a penetration of 'x' on from the
			** outside, but the treatment is identical (similar in
			** the case of a cylinder).
			*/

			dOpenAngle = wd->dOpenAngle*M_PI/180.; /* use radians */
			dRadius = wd->dRadius;
			vectorSub(p->r,O,vDum);
			a = vectorMagSq(vDum);
			a1 = dRadius + r1;
			a3 = max(dRadius - r1,0.0);
			if (a != 0.0 && a >= a3*a3 && a <= a1*a1) {
				a = sqrt(a);
				d = vectorDot(n,vDum);
				dAngle = acos(d/a);
				if (dAngle >= dOpenAngle) {
					vectorScale(vDum,dRadius/a,vDum2);
					vectorAdd(O,vDum2,q);
					bApplyForce = 1;
					}
				else {
					/* find point on rim closest to p */
					vectorScale(n,d,vDum2);
					vectorSub(vDum,vDum2,vDum2); /* (p - O) projected onto plane perp to n containing (0,0) */
					b = vectorMagSq(vDum2); /* |(p - O)|sin(dAngle) squared */
					if (b == 0.) { /* case where particle center lies directly on normal axis */
						b2 = dRadius*sin(dOpenAngle);
						if (b2 <= r1) {
							b1 = dRadius*cos(dOpenAngle);
							a2 = b2*b2 + (d - b1)*(d - b1); /* dist from center of particle to rim */
							if ((bApplyForce = (a2 <= r1sq))) {
								if (d > b1)
									vectorScale(n,d - sqrt(a2),vDum);
								else if (d < b1)
									vectorScale(n,d + sqrt(a2),vDum);
								else
									mdlassert(smf->pkd->mdl,0);
								vectorAdd(O,vDum,q); /* q is the "phantom point" (see finite cylinder for description) */
								}
							}
						}
					else {
						b1 = cos(dOpenAngle);
						vectorScale(vDum2,dRadius*sqrt((1. - b1*b1)/b),vDum2); /* norm(vDum2)*dRadius*sin(dOpenAngle) */
						vectorScale(n,dRadius*b1,vDum3);
						
						/*
						** vDum2 is now projection of desired point on rim onto plane perp to n containing (0,0)
						** vDum3 is now projection of desired point on rim onto n
						*/

						vectorAdd(vDum3,vDum2,vDum3);
						vectorSub(vDum,vDum3,vDum2);
						if ((bApplyForce = demCheckOverlapPoint(vDum2,r1sq)))
							vectorAdd(O,vDum3,q);
						}
					}
				}
			}

		else
			mdlassert(smf->pkd->mdl,0); /* if this is tripped, ensure wall type is supported */		

		if (dEpsN < 0.0) { /* death wall */
			if (bApplyForce) {

				/* clear particle-particle DEM elements */
				/*DEBUG why bother? the particle will just be deleted...*/
				for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++)
					demWipeOverlapElement(&p->overlaps[i]);

				/* clear particle-wall DEM elements */
				for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;i++)
					demWipeOverlapElement(&p->walloverlaps[i]);

				/* clear neighbors' contact lists of deceased particle -- alternatively, this could be done in a new function call in master just before msrAddDelParticle() */
				/*DEBUG check to see if this works*/
				DEM_ELEMENT *pne;
				int k;
				/* printf("here %i \n",p->iOrder); *//*DEBUG*/
				for (k=0;k<nSmooth;k++) {
					pn = nnList[k].pPart;
					if (pn->iOrder < p->iOrder)
						/* printf("now here %i %i\n",p->iOrder,pn->iOrder); *//*DEBUG*/
						for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++) {
							/*DEBUG! will this even work in parallel? need combiner cache...?*/
							pne = &pn->overlaps[i];
							/* printf("and here %i %i %i\n",p->iOrder,pne->iOrder,pn->iOrder); *//*DEBUG*/
							if (pne->iOrder == p->iOrder) {
#ifdef DEM_OVERLAP_OUTPUT
								fprintf(overlapOUT,"%g %d %d %i %e death\n",smf->dTime,pn->iOrder,pne->iOrder,pne->liOverlapCounter,x_maxWall[pn->iOrder][i]);
								x_maxWall[pn->iOrder][i] = 0.;
#endif /* DEM_OVERLAP_OUTPUT */
								/* printf("%i %i %i\n",p->iOrder,pne->iOrder,pn->iOrder); */
								demWipeOverlapElement(pne);
								if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"iOrder %i(%i) cleared from overlap list of neighbor %i\n",p->iOrder,smf->pkd->idSelf,pn->iOrder);
								break;
								}
							}
					}
				pkdDeleteParticle(smf->pkd,p);
				}
			continue;
			}

		/* determine bFoundNbr */

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"bApplyForce = %i, determining bFoundNbr: iOrder = %i(%i), wallID = %i\n",bApplyForce,p->iOrder,smf->pkd->idSelf,w->iWallID);

		for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;i++) {
			we = &p->walloverlaps[i];
			if (we->iOrder == w->iWallID) {
				mdlassert(smf->pkd->mdl,!bFoundNbr); /* ensure no duplicate overlaps */
				bFoundNbr = 1;
				iOverlap = i;
				}
			}

		/* two binaries, four cases: */
		if (bFoundNbr && bApplyForce) {
			we = &p->walloverlaps[iOverlap];
#ifndef DEM_TWOLAYERS
			++we->liOverlapCounter;
#endif
			}

		else if (bFoundNbr && !bApplyForce) {
			we = &p->walloverlaps[iOverlap];
#ifdef DEM_OVERLAP_OUTPUT
			fprintf(overlapOUT,"%g W_%d %d %i %e\n",smf->dTime,we->iOrder,p->iOrder,we->liOverlapCounter,x_maxWall[p->iOrder][j]);
			x_maxWall[p->iOrder][j] = 0.;
#endif /* DEM_OVERLAP_OUTPUT */
			demWipeOverlapElement(we);
			continue;
			}

		else if (!bFoundNbr && bApplyForce) {
			for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;i++) {
				we = &p->walloverlaps[i];
				if (we->iOrder == -1) {
					iOverlap = i;
					we->iOrder = w->iWallID;
					vectorZero(we->vShear);
					vectorZero(we->vnOld);
#ifdef DEM_ROTATION_DASHPOT
					vectorZero(we->vRoll);
                                        we->dTwist = 0.0;
#endif /* DEM_ROTATION_DASHPOT */
#ifndef DEM_TWOLAYERS
					we->liOverlapCounter = 1;
#endif
#ifdef DEM_DIAG
					vectorZero(we->vFn);
					vectorZero(we->vFt);
#endif /* DEM_DIAG */
					break;
					}
				}
			mdlassert(smf->pkd->mdl,i < MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS); /* ensure free spot was indeed found and used */
			}

		else if (!bFoundNbr && !bApplyForce)
			continue;

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"bApplyForce = %i, bFoundNbr = %i: iOrder = %i(%i), wallID = %i\n",bApplyForce,bFoundNbr,p->iOrder,smf->pkd->idSelf,we->iOrder);

		mdlassert(smf->pkd->mdl,iOverlap >= 0 && we->iOrder == w->iWallID && iOverlap < MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS);

		we = &p->walloverlaps[iOverlap]; /* we now points to the overlap of this particle and its neighbor */

		m1 = p->fMass;
		dReducedMass = m1; /*DEBUG assumes non-reactive (i.e., infinite-mass) wall*/
		dSqrtReducedMass = sqrt(dReducedMass);
		vectorSub(q,p->r,n);
		d = vectorMag(n); /* scalar distance to wall */
		x = r1 - d; /* probably want to be checking x/d for avg and max, warning if max is large... */
		i1 = 0.4*m1*r1*r1;

#ifdef DEM_TWOLAYERS

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"iOrder = %i(%i), wallID = %i, vShear before boundary cross check = %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,vectorMag(we->vShear));

		dInBnd = r1*dInBndFrac;
		if (x > dInBnd) {
			kn = kn_inner;
			kt = kt_inner;
			x_inner = x - dInBnd;
			mdlassert(smf->pkd->mdl,x_inner >= 0.0);
			F_outershell = dInBnd*kn_outer;
			if (we->liOverlapCounter <= 0)
				--we->liOverlapCounter;
			else { /* scale static spring if boundary layer was just crossed (conserve energy in tang spring: 1/2 kt*(vShear)^2 */
#ifdef DEM_OVERLAP_OUTPUT
				fprintf(overlapOUT,"inner boundary crossed - %g W_%d %d %li %e\n",smf->dTime,we->iOrder,p->iOrder,we->liOverlapCounter,x_maxWall[p->iOrder][j]);
#endif
				vectorScale(we->vShear,sqrt(kt_outer/kt_inner),we->vShear); /* vnOld remains the same */
				we->liOverlapCounter = -1;
				}
			}
		else {
			kn = kn_outer;
			kt = kt_outer;
			x_inner = -1.;
			F_outershell = 0.;
			if (we->liOverlapCounter >= 0)
				++we->liOverlapCounter;
			else { /* scale static spring... */
#ifdef DEM_OVERLAP_OUTPUT
				fprintf(overlapOUT,"inner boundary crossed - %g W_%d %d %li %e\n",smf->dTime,we->iOrder,p->iOrder,we->liOverlapCounter,x_maxWall[p->iOrder][j]);
#endif
				vectorScale(we->vShear,sqrt(kt_inner/kt_outer),we->vShear);
				we->liOverlapCounter = 1;
				}
			}
		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"iOrder = %i(%i), wallID = %i, vShear after boundary cross check = %.16e, overlap = %.16e, inner boundary = %.16e, overlap with inner layer = %.16e, kt_inner = %.16e, kt_outer = %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,vectorMag(we->vShear),x,dInBnd,x_inner,kt_inner,kt_outer);

#endif /* DEM_TWOLAYERS */
		
#ifdef DEM_OVERLAP_OUTPUT
		if (x > x_maxWall[p->iOrder][j])
			x_maxWall[p->iOrder][j] = x;
#endif /* DEM_OVERLAP_OUTPUT */

		{
		int iWarningType = DEM_NOERR;
		assert(!PARTICLE_STUCK(p));
		if (x > smf->FO.DP.dErrorFrac*r1)
			iWarningType = DEM_ERROR;
		else if (x > smf->FO.DP.dMajorFrac*r1)
			iWarningType = DEM_MAJOR;
		else if (x > smf->FO.DP.dMinorFrac*r1)
			iWarningType = DEM_MINOR;
		if (iWarningType != DEM_NOERR) {
#if (INTERNAL_WARNINGS != 0)
			static int nWarnWall = 1;
			if (nWarnWall == 1 || nWarnWall%INTERNAL_WARNINGS == 0 || iWarningType != DEM_MINOR) {
				double v_mag = sqrt((p->v[0] - w->vTotVel[0])*(p->v[0] - w->vTotVel[0]) +
									(p->v[1] - w->vTotVel[1])*(p->v[1] - w->vTotVel[1]) +
									(p->v[2] - w->vTotVel[2])*(p->v[2] - w->vTotVel[2]));
				fprintf(stderr,"DEM WALL %s #%i (pid=%i,T=%g) %s overlap: %i & Wall %i, %f%% (relative v = %g)\n",
						iWarningType == DEM_ERROR ? "ERROR" : "WARNING",nWarnWall,smf->pkd->idSelf,smf->dTime,
							   iWarningType == DEM_MINOR ? "minor" : iWarningType == DEM_MAJOR ? "major" : "fatal",
							   p->iOrder,w->iWallID,100.0*x/r1,v_mag); /* maximum overlap is 100% (because walls have no thickness) */
				}
			++nWarnWall;
#endif /* INTERNAL_WARNINGS */
			mdlassert(smf->pkd->mdl,iWarningType != DEM_ERROR);
			}
		}
		
		static int bCnCtPreFacCalculatedWalls[MAX_NUM_WALLS]; /* since dEpsN and dEpsT are identical for all particles in sim (but not for each wall) */

		static int bCnCtPreFacCalculatedWallsInitialized = 0;
		if (!bCnCtPreFacCalculatedWallsInitialized) {
			for (i=0;i<MAX_NUM_WALLS;i++)
				bCnCtPreFacCalculatedWalls[i] = 0;
			bCnCtPreFacCalculatedWallsInitialized = 1;
			}

		if (!bCnCtPreFacCalculatedWalls[w->iWallID]) {

			/* damping term: normal */
			mdlassert(smf->pkd->mdl,dEpsN >= 0.0); /* zero = sticky wall */
			if (dEpsN > 0.0) { 
				dLnEpsN = log(dEpsN);
				dLnEpsNsq = dLnEpsN*dLnEpsN;
				a = pi_sq + dLnEpsNsq;
				CnPreFacWalls[w->iWallID] = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn/a);
				CnPreFacWalls[w->iWallID] += CnPreFacWalls[w->iWallID];
#ifdef DEM_TWOLAYERS
				CnInnerPreFacWalls[w->iWallID] = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn_inner/a);
				CnInnerPreFacWalls[w->iWallID] += CnInnerPreFacWalls[w->iWallID];
				CnOuterPreFacWalls[w->iWallID] = -sign(dLnEpsN)*sqrt(dLnEpsNsq*kn_outer/a);
				CnOuterPreFacWalls[w->iWallID] += CnOuterPreFacWalls[w->iWallID];
#endif
				}
			else {
				CnPreFacWalls[w->iWallID] = 0.0;
#ifdef DEM_TWOLAYERS
				CnInnerPreFacWalls[w->iWallID] = CnOuterPreFacWalls[w->iWallID] = 0.0-sign(dLnEpsN)*sqrt(dLnEpsNsq*kn_outer/a);
#endif
				}

			/* damping term: tangential */
			mdlassert(smf->pkd->mdl,dEpsT >= 0.0);
			if (dEpsT > 0.0) {
				dLnEpsT = log(dEpsT); /* could this be log(fabs(dEpsT))? */
				dLnEpsTsq = dLnEpsT*dLnEpsT;
				a = pi_sq + dLnEpsTsq;
				CtPreFacWalls[w->iWallID] = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt/a);
				CtPreFacWalls[w->iWallID] += CtPreFacWalls[w->iWallID];
#ifdef DEM_TWOLAYERS
				CtInnerPreFacWalls[w->iWallID] = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt_inner/a);
				CtInnerPreFacWalls[w->iWallID] += CtInnerPreFacWalls[w->iWallID];
				CtOuterPreFacWalls[w->iWallID] = -sign(dLnEpsT)*sqrt(dLnEpsTsq*kt_outer/a);
				CtOuterPreFacWalls[w->iWallID] += CtOuterPreFacWalls[w->iWallID];
#endif
				}
			else {
				CtPreFacWalls[w->iWallID] = 0.0;
#ifdef DEM_TWOLAYERS
				CtInnerPreFacWalls[w->iWallID] = CtOuterPreFacWalls[w->iWallID] = 0.0;
#endif
				}
			printf("CnPreFacWalls[%d] and CtPreFacWalls[%d] calculated on pid %i\n",w->iWallID,w->iWallID,smf->pkd->idSelf); /*DEBUG*/
			bCnCtPreFacCalculatedWalls[w->iWallID] = 1;
			}

		/* components of a unit normal vector from particle center to wall */
		mdlassert(smf->pkd->mdl,d != 0.0);
		vectorScale(n,1./d,n);

		/* velocity of particle at contact point w.r.t. its COM */
		vectorCross(p->wPred,n,vDum);
		vectorScale(vDum,d,s1);

		/* velocity of wall at contact point */
#ifdef DEM_WALLS_REACT
		if (wd->dMass > 0.0) {
		    if (wd->iType == WallTriangle || wd->iType == WallRectangle) {
			vectorCross(wd->vVertex1, wd->vVertex2, vVinyl); /* vector normal to wall */
			vectorNorm(vVinyl); /* unit vector */
			vectorScale(vVinyl, vectorDot(wd->vRot, vVinyl), vVinyl); /* component of rotation around normal */
			}
		    else if (wd->iType == WallCylinderFinite || wd->iType == WallDisk) {
			vectorScale(wd->vOrient, vectorDot(wd->vRot, wd->vOrient), vVinyl);
			}
		    else if (wd->iType == WallShell) {
			vectorCopy(wd->vRot, vVinyl);
			}
		    else
			assert(0); /* infinite walls cannot be reactive walls */
		    vectorSub(q, wd->vCoG, vDum2); /* lever-arm (vector) */
		    vectorCross(wd->vRot, vDum2, v); /* relative velocity at contact point due to wall rotation */
		    vectorAdd(v, w->vTotVel, v); /*DEBUG! w->vTotVel defined at half step---see msrWallsMove()*/
		    }
		else {	/* normal wall */
			if (wd->dAngSpeed != 0.) { /* save time if wall is not spinning */
				vectorScale(wd->vOrient,wd->dAngSpeed,vVinyl); /* spin of wall (vector) */
				vectorSub(q,O,vDum2); /* lever-arm (vector) */
				vectorCross(vVinyl,vDum2,v); /*DEBUG: not fully tested! at this point, v should be velocity at CP due to rotation of wall*/
				vectorAdd(v,w->vTotVel,v); /*DEBUG: and at this point, v should be total velocity of wall at the contact point*/
				}
			else {
				vectorCopy(w->vTotVel,v);
				vectorZero(vVinyl); /* used later when computing torques */
				}
		}
#else
		if (wd->dAngSpeed != 0.) { /* save time if wall is not spinning */
			vectorScale(wd->vOrient,wd->dAngSpeed,vVinyl); /* spin of wall (vector) */
			vectorSub(q,O,vDum2); /* lever-arm (vector) */
			vectorCross(vVinyl,vDum2,v); /*DEBUG: not fully tested! at this point, v should be velocity at CP due to rotation of wall*/
			vectorAdd(v,w->vTotVel,v); /*DEBUG: and at this point, v should be total velocity of wall at the contact point*/
			}
		else {
			vectorCopy(w->vTotVel,v);
			vectorZero(vVinyl); /* used later when computing torques */
			}
#endif /* !DEM_WALLS_REACT */
		/*DEBUG: test*/
		if (dEpsN == 0.0) { /* sticky wall */
#ifdef AGGS
			mdlassert(smf->pkd->mdl,!IS_AGG(p)); /* way too complicated to let aggs stick to walls... *//*DEBUG if we ever allow particles to merge to form aggs with DEM, be aware of the possibility of one of the particles already being stuck to a wall...*/
#endif /* AGGS */
			if (p->iColor >= 0) p->iColor = -1 - w->iWallID;
			vectorZero(p->a); /* not essential: stuck particles are not kicked */
			vectorCopy(w->vTotVel,p->v); /* transfer wall motion to particle */
			vectorZero(p->w);
			if (COLLIDER_STUCK_ON_ROTATING_WALL(p,WP)) {
				/* transfer wall spin to particle */
				mdlassert(smf->pkd->mdl,wd->iType != WallTriangle && wd->iType != WallRectangle); /* not supported at the moment */
				vectorScale(wd->vOrient,wd->dAngSpeed,p->w);
				}
			return;
			}

                Un = Ut = 0.;
#ifdef DEM_ROTATION_DASHPOT
                for (i=0;i<3;i++) {
                        v[i] -= p->vPred[i]; /* particle velocity */ /*DEBUG! make sure no sliding patch with walls*/
                        s[i] = -s1[i]; /* particle spin */
                        u[i] = v[i] + s[i];
                        Un += u[i]*n[i];
                        }
                for (i=0;i<3;i++) {
                        un[i] = Un*n[i];
                        ut[i] = u[i] - un[i];
                        Ut += ut[i]*ut[i];
                        }
#else /* Directional constant torque models */
                Rt = 0.;
		if (mu_r != 0.) {
			for (i=0;i<3;i++) {
				v[i] -= p->vPred[i]; /* particle velocity */ /*DEBUG! make sure no sliding patch with walls*/
				s[i] = -s1[i]; /* particle spin */
				r[i] = s1[i]; /* rolling */ /*DEBUG: how to treat wall motion/rotation affect on spin? (no effect for now)*/
				u[i] = v[i] + s[i];
				Un += u[i]*n[i];
				}
			for (i=0;i<3;i++) {
				un[i] = Un*n[i];
				ut[i] = u[i] - un[i];
				Ut += ut[i]*ut[i];
				Rt += r[i]*r[i];
				}
			Rt = sqrt(Rt);
			}
		else {
			for (i=0;i<3;i++) {
				v[i] -= p->vPred[i]; /* particle velocity */ /*DEBUG! make sure no sliding patch with walls*/
				s[i] = -s1[i]; /* particle spin */
				u[i] = v[i] + s[i];
				Un += u[i]*n[i];
				}
			for (i=0;i<3;i++) {
				un[i] = Un*n[i];
				ut[i] = u[i] - un[i];
				Ut += ut[i]*ut[i];
				}
			}
#endif
                Ut = sqrt(Ut);
		/* unit tangential velocity vector */
		if (Ut != 0.0)
			vectorScale(ut,1./Ut,t);
		else
			vectorZero(t); /* avoid dividing by zero */

#ifndef DEM_ROTATION_DASHPOT
		/* unit rolling vector */
		if (Rt != 0.0) /* Rt always zero if mu_r is zero (above), so !Rt condition is better here */
			vectorScale(r,1./Rt,rt); /* rt[i] = unit vector of spinning component at contact point */
		else
			vectorZero(rt); /* rt[i] = unit vector of spinning component at contact point */
#endif  /* Directional constant torque models */

		/* damping terms */
#ifndef DEM_TWOLAYERS
		Cn = CnPreFacWalls[w->iWallID]*dSqrtReducedMass;
		Ct = CtPreFacWalls[w->iWallID]*dSqrtReducedMass; /* ...relating dEpsN,dEpsT to Cn,Ct comes from Cleary etal. '98 */
#else  /* if an inner boundary is defined where stiffness changes */
		if (x > dInBnd) {
			Cn = CnInnerPreFacWalls[w->iWallID]*dSqrtReducedMass;
		   	Ct = CtInnerPreFacWalls[w->iWallID]*dSqrtReducedMass;
			}
		else {
			Cn = CnOuterPreFacWalls[w->iWallID]*dSqrtReducedMass;
		   	Ct = CtOuterPreFacWalls[w->iWallID]*dSqrtReducedMass;
			}
#endif /* DEM_TWOLAYERS */

		a1 = 1./m1;
		b1 = 1./i1;

#ifndef DEM_ROTATION_DASHPOT
		if (Rt != 0. || mu_t != 0.) {
#endif
			if (bUseLegacyRoll)
				dRollLever = d;
			else if (bUseContactRadius) {
				a = r1sq - d*d;
				mdlassert(smf->pkd->mdl,a >= 0.0);
				dRollLever = dContactRadius = sqrt(a);
				}
			else
			    dRollLever = dReducedRadius = r1; /*DEBUG assumes non-reactive (i.e., infinite-mass/radius) wall*//*also, dReducedRadius not used below...*/
#ifndef DEM_ROTATION_DASHPOT
			}
#else
		if (!bUseContactRadius)
			dtemp1 = beta_s * r1;
		else
			dtemp1 = beta_s * dContactRadius;
		dtemp1 *= dtemp1;
		kr = kn * dtemp1; /* spring term */
		Cr = Cn * dtemp1; /* damping term */
		dtemp2 = 2.0 * dtemp1;
		kq = kt * dtemp2; /* spring term */
		Cq = Ct * dtemp2; /* damping term */
		for (i=0;i<3;i++)
			vRelW[i] = vVinyl[i] - p->wPred[i];  /*DEBUG needs testing, e.g., use orientation tracking to see if sphere comes to equilibrium with rotating floor*/
		dRelTwi = vectorDot(n,vRelW);
		dtemp2 = 0.5/d;
		vectorCross(n,u,vtemp1);
		for (i=0;i<3;i++) {
			vRelTwi[i] = dRelTwi * n[i];
			vRelRol[i] = vRelW[i] - vRelTwi[i] + dtemp2*vtemp1[i];
			}
#endif /* DEM_ROTATION_DASHPOT */

		/* Tangential displacement vector */

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Starting calculation on tangential displacement vector: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,we->iOrder);

		if (vectorMagSq(we->vnOld) != 0.0) {
#ifndef DEM_TWOLAYERS /* if boundary is crossed, vnOld will be defined and |liOverlapCounter| == 1 */
			//mdlassert(smf->pkd->mdl,we->liOverlapCounter > 1); /* should not be first step in overlap *//*DEBUG assert can trip on restart from dem file*/
#endif

			/* see comments in particle-particle section above */

			/*
			** If vnOld is not stored in the DEM element struct, we need to compute it somehow.
			** This is straightforward only in the case where vnOld can be easily predicted, as
			** in the case where a particle will interact only with the face of a wall, where
			** n is constant in time (n = vnOld). Otherwise *very* complicated!
			*/

			/*
			** We first consider rotation around the normal - for this,
			** we simply rotate the tangential spring around the average
			** of the previous normal and the current normal.  So we use
			** the normal and the spins of 1/2 step ago (midpoint).
			*/

			vectorAdd(p->w,vVinyl,vDum);
			vectorAdd(we->vnOld,n,vDum2);
			vectorNorm(vDum2);

			/* could consider checking that this is nonzero before expensive next op, but only in special cases will it be nonzero... */
			vectorRotate(we->vShear,vDum2,0.5*dDelta*vectorDot(vDum,vDum2));
#ifdef DEM_ROTATION_DASHPOT
			vectorRotate(we->vRoll,vDum2,0.5*dDelta*vectorDot(vDum,vDum2));
#endif /* DEM_ROTATION_DASHPOT */
			/*
			** Next, we change the orientation of the displacement vector
			** in the same way that the orientation of the normal has
			** changed over the previous step, and assume that the normal
			** did not change by more than 90 degrees
			*/

			/*
			vectorScale(n,vectorDot(n,we->vShear),vDum);
			vectorSub(we->vShear,vDum,vDum);
			if ((a = vectorMagSq(vDum)) != 0.0)
				vectorScale(vDum,sqrt(vectorMagSq(we->vShear)/a),we->vShear);
			*/

			vectorCross(n,we->vnOld,vDum); // axis of rotation (not normalized)
			if ((a = vectorMagSq(vDum)) != 0.0) {
				a = sqrt(a);
				vectorScale(vDum,1./a,vDum);
				vectorRotate2(we->vShear,vDum,a,vectorDot(n,we->vnOld)); // converting to an angle and using vectorRotate() failed when recomputing the sin and cos
#ifdef DEM_ROTATION_DASHPOT
                                vectorRotate2(we->vRoll,vDum,a,vectorDot(n,we->vnOld));
#endif /* DEM_ROTATION_DASHPOT */
                                }
			}

		/* set vnOld to current n for use in next step in the case that the overlap persists to the next step */
		vectorCopy(n,we->vnOld);

		N = T = 0.; /* normal force */

#ifdef DEM_TWOLAYERS
		if (x_inner >= 0.0)
			x = x_inner; /* if an inner boundary is defined where stiffness changes and is penetrated */
#endif

		/* forces */

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Starting calculation on static friction effects: iOrder = %i(%i), wallID = %i\n",p->iOrder,smf->pkd->idSelf,we->iOrder);

		for (i=0;i<3;i++) {
			we->vShear[i] += ut[i]*dDelta; /* this step's contribution to the tangential spring */
			Fn[i] = -(kn*x + F_outershell)*n[i] + Cn*un[i]; /* check to ensure that all other forces are computed first (e.g. gravity) */
			Ft[i] = kt*we->vShear[i] + Ct*ut[i]; /*      ""      */
#ifdef DEM_ROTATION_DASHPOT
			we->vRoll[i] += vRelRol[i]*dDelta;
#endif /* DEM_ROTATION_DASHPOT */
			N += Fn[i]*Fn[i]; /* square of normal force */
			T += Ft[i]*Ft[i]; /* square of tangential force */

			if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"applying forces and moments: iOrder = %i(%i), wallID = %i, Ft[%i] = %.16e, vShear[%i] = %.16e, un[%i] = %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,i,Ft[i],i,we->vShear[i],i,un[i]);
			}
#ifdef DEM_ROTATION_DASHPOT
		we->dTwist += dRelTwi*dDelta;
#endif /* DEM_ROTATION_DASHPOT */
		
		F2 = N*mu_s*mu_s; /* square of max allowed static friction */
		N = sqrt(N); /* scalar magnitude of Normal force (not squared) */

#ifdef DEM_ROTATION_DASHPOT
		if (!bUseContactRadius) {
			MrCrit = mu_r * beta_s * r1 * N;
			MtCrit = mu_t * beta_s * r1 * mu_s * N;
			}
		else {
			MrCrit = mu_r * beta_s * dContactRadius * N;
			MtCrit = mu_t * beta_s * dContactRadius * mu_s * N;
			}

#endif /* DEM_ROTATION_DASHPOT */

#ifdef DEM_PRESSURE
		p->dPressure += N; /*DEBUG: hack to sum total pressure on particle*/
#endif

		/*DEBUG*/
		/*
		printf("DoDEM() call = %d wall = %d pe->liOverlapCounter = %li x/r = %e N = %e |S|/r = %e S/|S| dot n = %e\n",count,w->iWallID,we->liOverlapCounter,x/r1,N,vectorMag(we->vShear)/r1,vectorDot(we->vShear,n)/vectorMag(we->vShear));
		printf("S[0]/r = %e t[0] = %e  ;  F = %e  ;  Ct*Ut = %e  ;  Ft[0] = %e p->w[0] = %e\n",we->vShear[0]/r1,t[0],sqrt(F2),Ct*Ut,Ft[0],p->w[0]);
		printf("S[1]/r = %e t[1] = %e  ;  b1 = %e  ;  Ut = %e  ;  Ft[1] = %e p->w[1] = %e\n",we->vShear[1]/r1,t[1],b1,Ut,Ft[1],p->w[1]);
		printf("S[2]/r = %e t[2] = %e  ;  l1 = %e  ;  Ct*Ut = %e  ;  Ft[2] = %e Fn[2] = %e\n",we->vShear[2]/r1,t[2],d,Ct*Ut,Ft[2],Fn[2]);
		printf("n[0] = %e n[1] = %e n[2] = %e p=%e %e %e w=%e %e %e\n",n[0],n[1],n[2],p->r[0],p->r[1],p->r[2],wd->vOrigin[0],wd->vOrigin[1],wd->vOrigin[2]);
		*/

		if (F2 < T) {
			if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"Static friction limit exceeded: iOrder = %i(%i), wallID = %i\n, tangential force = %.16e, limit = %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,sqrt(T),sqrt(F2));


			F = N*mu_s; /* mag of max allowed force du friction statique (not squared) */
			if (dTangentialSpringDrag != 0.0) {

				/* compute the new magnitude of vShear, but if tangential kintetic friction alone exceeds static friction, vShear will be zeroed */
				a = F - Ct*Ut;
				b = a <= 0. ? 0. : vectorMag(we->vShear);
				a = b <= 0. ? 0. : dTangentialSpringDrag*a/(kt*b);

				/* set vShear */
				vectorScale(we->vShear,a,we->vShear); /* vShear for this step */
				vectorScale(t,min(F,Ct*Ut),vDum); /* the force due to tangential kinetic friction */
				vectorScale(we->vShear,kt,vDum2);
				vectorAdd(vDum,vDum2,vDum); /* total tangential force due to both static and kinetic friction */

				mdlassert(smf->pkd->mdl,T > 0.);
				if (a != 0.0)
					vectorScale(Ft,a*b/sqrt(T),we->vShear);
				else
					vectorZero(we->vShear);					
				vectorCopy(vDum,Ft);
				}

			else {
				vectorZero(we->vShear);
				vectorScale(t,min(F,Ct*Ut),Ft);
				}
                        }
                vectorCross(Ft,n,Ftn); /* tangential DEM induced torque per unit length */
#ifdef DEM_ROTATION_DASHPOT
                for (i=0;i<3;i++)  /* rolling */
                        vRft[i] = kr*we->vRoll[i];
                dtemp1 = vectorMag(vRft);
                if (dtemp1 > MrCrit) {
                        vectorScale(vRft, 1./dtemp1, vtemp1);
                        for (i=0;i<3;i++) {
                                vRft[i] = MrCrit * vtemp1[i];
                                we->vRoll[i] -= vRelRol[i] * dDelta;
                                }
                        }
                else
                        for (i=0;i<3;i++)
                                vRft[i] += Cr * vRelRol[i];
                dtemp1 = kq*we->dTwist;  /* twisting */
                if ( fabs(dtemp1) > MtCrit ) {
                        if ( dtemp1 > 0.0 )
                                for (i=0;i<3;i++)
                                        vTfn[i] = MtCrit * n[i];
                        else
                                for (i=0;i<3;i++)
                                        vTfn[i] = - MtCrit * n[i];
                        we->dTwist -= dRelTwi * dDelta;
                        }
                else
                        for (i=0;i<3;i++)
                                vTfn[i] = dtemp1 * n[i] + Cq * vRelTwi[i];
#else  /* Directional constant torque models */
                if (Rt != 0.)
                        dRollingFriction = mu_r*N;
                else
                        dRollingFriction = 0.;

                /* Below, we use "twisting friction" to damp out spin along the normal axis connecting the particles' centers */
                if (mu_t != 0.) {

                        /* make [dTwistRate * n-hat] be the relative differential spin around normal axis: */
                        vectorSub(vVinyl,p->wPred,vDum);
                        dTwistRate = vectorDot(vDum,n);

                        dTwistingFriction = mu_t*N*sign(dTwistRate); /* this dotted with n-hat is the twisting friction force */
                        /* printf("WALL dContactRadius=%e\nd*d=%e\nr1sq=%e\n",dContactRadius,d*d,r1sq); */ /*DEBUG:make this debug_trace*/
                        }

                /* moments */
                if (Rt != 0.) {
                        vectorCross(rt,n,vDum); /* unit vector giving angular direction of rolling friction torque */
                        vectorScale(vDum,dRollingFriction,vRft); /* rolling friction torque per unit length */
                        /* below we check to see if rolling friction impulse is greater than rolling momentum. if so, make equal */
                        a = vectorMagSq(vRft);
                        b = dDelta*b1*d*d;
                        b *= a*b;
                        //			printf("%e %e %e %e %e %e %e\n",Rt,dDelta*b1*vRft[0],sqrt(a),dDelta*b1*dDelta*b1*a,b,vectorMag(p->w)*d,vectorMag(p->wPred)*d);
                        if (Rt*Rt < b) {
                                vectorScale(vRft,vectorMag(p->w)*d/sqrt(b),vRft);
                                //				printf("geraldine ferraro %d %e\n",p->iOrder,vRft[0]);
                                }
                        }
                else
                        vectorZero(vRft);
                if (mu_t != 0.)
                        vectorScale(n,dContactRadius*dTwistingFriction,vTfn); /* twisting friction torque */
                else
                        vectorZero(vTfn);
#endif /* DEM_ROTATION_DASHPOT */

		/* apply forces, torques */

		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"applying forces and moments: iOrder = %i(%i), wallID = %i, forces on particle COM: normal = %.16e, tangential = %.16e, moments: tangential = %.16e, rolling friction = %.16e, twisting friction = %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,vectorMag(Fn),a1*vectorMag(Ft),(-1.)*b1*d*vectorMag(Ftn),b1*d*vectorMag(vRft),b1*d*vectorMag(vTfn));

		/*DEBUG*/
		/*
		printf("b DoDEM() call = %d wall = %d pe->liOverlapCounter = %li x/r = %e N = %e |S|/r = %e S/|S| dot n = %e\n",count,w->iWallID,we->liOverlapCounter,x/r1,N,vectorMag(we->vShear)/r1,vectorDot(we->vShear,n)/vectorMag(we->vShear));
		printf("b S[0]/r = %e t[0] = %e  ;  F = %e  ;  Ct*Ut = %e  ;  Ft[0] = %e p->w[0] = %e\n",we->vShear[0]/r1,t[0],sqrt(F2),Ct*Ut,Ft[0],p->w[0]);
		printf("b S[1]/r = %e t[1] = %e  ;  b1 = %e  ;  Ct*Ut = %e  ;  Ft[1] = %e p->w[1] = %e\n",we->vShear[1]/r1,t[1],b1,Ct*Ut,Ft[1],p->w[1]);
		printf("b S[2]/r = %e t[2] = %e  ;  l1 = %e  ;  Ct*Ut = %e  ;  Ft[2] = %e p->w[2] = %e\n",we->vShear[2]/r1,t[2],d,Ct*Ut,Ft[2],p->w[2]);
		*/
#ifdef DEM_COHESION
		if (cohe_model == 0)
			if (!bUseContactRadius)
				dtemp_cohe = cohe_coeff * 4.0 * beta_s * beta_s * r1sq;
			else
				dtemp_cohe = cohe_coeff * M_PI * beta_s * beta_s * dContactRadius*dContactRadius;
		else
			dtemp_cohe = cohe_coeff * 2.0*M_PI * r1sq; /* Sanchez's method */
		Fn[0] +=  dtemp_cohe*n[0];
		Fn[1] +=  dtemp_cohe*n[1];
		Fn[2] +=  dtemp_cohe*n[2];
#endif /* DEM_COHESION */
#ifndef AGGS // JCM 10.13.20
		for (i=0;i<3;i++) {
			p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
#ifdef DEM_ROTATION_DASHPOT
			p->wDot[i] -= b1*(d*Ftn[i] - vRft[i] - vTfn[i]);
#else
			p->wDot[i] -= b1*(d*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); /* might try a factor of 10xdMu_r to make Contact Radius commensurate with legacy RF */
#endif /* !DEM_ROTATION_DASHPOT */
			//printf("Fn[%d]=%e Ft[%d]=%e\n",i,Fn[i],i,Ft[i]);
#ifdef DEM_DIAG
			we->vFn[i] = Fn[i];
			we->vFt[i] = Ft[i];
#endif /* DEM_DIAG */
			}
		/*DEBUG!{
			double FnMag,FtMag;
			FnMag = a1*sqrt(Fn[0]*Fn[0] + Fn[1]*Fn[1] + Fn[2]*Fn[2]);
			FtMag = a1*sqrt(Ft[0]*Ft[0] + Ft[1]*Ft[1] + Ft[2]*Ft[2]);
			if (FnMag > 2000 || FtMag > 2000)
				printf("SOUND %e %i %i %g %g\n",
					   smf->dTime,p->iOrder,-1 - w->iWallID,FnMag,FtMag);
					   }*/
#else /* AGGS */
            if (p->iOrgIdx < 0) {
                vectorScale(n,d,vAggTemp); /* Vector from particle center to contact point. */
                vectorSub(p->r_unwrapped,p->agg_COM,vBodyPos); /* Location of contact point in space. */
                vectorAdd(vBodyPos,vAggTemp,vAggLever); /* Vector from agg COM to contact. */
                vectorAdd(Fn,Ft,vAggForces); /* Collect tangential and normal forces. */
                vectorScale(vAggForces,a1,vAggTemp2); /* Get acceleration from force */
                vectorCross(vAggLever,vAggForces,vAggTorques); /* Calculate torques from norm. and tan. forces. */
                vectorAdd(p->DEM_aggaccels,vAggTemp2,p->DEM_aggaccels); /* Add in acceleration from walls to total tally. */
                vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques); /* Ditto for torques. */
#ifdef DEM_ROTATION_DASHPOT
                /* ROLLING */
                vectorCross(vRft,vAggTemp,vAggTemp2); /* Vector parallel to rolling friction force. */
                if (vectorMag(vRft) == 0 || vectorMag(vAggTemp2) == 0) {
#if (INTERNAL_WARNINGS != 0)
                    static int nWarn = 1;
                    if (mu_r != 0 && (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)) {
                        fprintf(stderr,"Particle: %d: WARNING: Rolling torque too small to apply to aggregate. Setting rolling torque to zero.\n", p->iOrder);
                    }
                    ++nWarn;
#endif /* INTERNAL_WARNINGS */
                } else {
                    vectorNorm(vAggTemp2); /* Unit vector in direction of rolling friction force. */
                    vectorScale(vAggTemp2,vectorMag(vRft)/vectorMag(vAggTemp),vRollForces); /* The net rolling force. */
                    vectorCross(vAggLever,vRollForces,vAggTorques); /* Net rolling torque on the agg. */
                    vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
                }
                /* TWISTING */
                vectorCopy(vTfn,vAggTorques);
                vectorAdd(p->DEM_aggtorques,vAggTorques,p->DEM_aggtorques);
#else

#endif /* !DEM_ROTATION_DASHPOT */
            } else {
                for (i=0;i<3;i++) {
                    p->dDeltaAccel[i] += (Fn[i] + Ft[i])*a1;
#ifdef DEM_ROTATION_DASHPOT
                    p->wDot[i] -= b1*(d*Ftn[i] - vRft[i] - vTfn[i]);
#else
                    p->wDot[i] -= b1*(d*Ftn[i] - dRollLever*vRft[i] - vTfn[i]); 
#endif /* !DEM_ROTATION_DASHPOT */
                }
            }
#endif /* AGGS */
		if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"particle-wall forces and moments applied: iOrder = %i(%i), wallID = %i, current tally of particle force: %.16e, moment: %.16e\n",p->iOrder,smf->pkd->idSelf,we->iOrder,vectorMag(p->dDeltaAccel),vectorMag(p->wDot));

#ifdef DEM_WALLS_REACT
		//This ifdef is inside of the wall loop (loops through each wall)
		vectorCopy(wd->vCoG,rWrist);

		vectorZero(rFnT); //Initialize Vector components to Zero;
		vectorZero(rFtT);
		vectorZero(rArm); //Moment arm
		if (wd->dMass != 0.0){ //only do this for walls that have mass -rb ... wd=Wall Data.. defined in master
		    p->dXForceOnWalls[wd->iAssembly] -= (Fn[0] + Ft[0]); //negative reaction force -rb
		    p->dYForceOnWalls[wd->iAssembly] -= (Fn[1] + Ft[1]);
		    p->dZForceOnWalls[wd->iAssembly] -= (Fn[2] + Ft[2]);
		    vectorSub(p->r,rWrist,rArm);
		    vectorCross(rArm,Fn,rFnT); //Fxr -- negative force since we want the reaction force
		    vectorCross(rArm,Ft,rFtT);

#ifdef DEM_ROTATION_DASHPOT
		    p->dXTorqueOnWalls[wd->iAssembly] += (d*Ftn[0] - vRft[0] - vTfn[0]);
		    p->dYTorqueOnWalls[wd->iAssembly] += (d*Ftn[1] - vRft[1] - vTfn[1]);
		    p->dZTorqueOnWalls[wd->iAssembly] += (d*Ftn[2] - vRft[2] - vTfn[2]);
#else
		    p->dXTorqueOnWalls[wd->iAssembly] += (d*Ftn[0] - dRollLever*vRft[0] - vTfn[0]);
		    p->dYTorqueOnWalls[wd->iAssembly] += (d*Ftn[1] - dRollLever*vRft[1] - vTfn[1]);
		    p->dZTorqueOnWalls[wd->iAssembly] += (d*Ftn[2] - dRollLever*vRft[2] - vTfn[2]);
#endif /* DEM_ROTATION_DASHPOT */
		    p->dXTorqueOnWalls[wd->iAssembly] -= (rFnT[0] + rFtT[0]);
		    p->dYTorqueOnWalls[wd->iAssembly] -= (rFnT[1] + rFtT[1]);
		    p->dZTorqueOnWalls[wd->iAssembly] -= (rFnT[2] + rFtT[2]);
		  }
#endif /* DEM_WALLS_REACT */
    }

#endif /* WALLS */

	if (DEBUG_TRACE(p->iOrder)) fprintf(stderr,"function DoDEM() complete: iOrder = %i(%i)\n",p->iOrder,smf->pkd->idSelf);

/*#define CLEAR_BELOW_HOPPER*/ /*DEBUG: relevant for hopper only*/
#ifdef CLEAR_BELOW_HOPPER
	if (p->r[2] < -1.804836963060900027e-11/*-270*/) pkdDeleteParticle(smf->pkd,p);
#endif

	/*DEBUG!!!:Hackity-Hack*/
	/*
	double W;
	W = vectorMagSq(p->w);
	if (r1*W > 1.e-7) {
		vectorScale(p->w,(1./W),p->w);
		printf("particle spins limited to ~3mm sec DEBUG!!!\n");
		}
	*/
	/*vectorZero(p->w);*/
	}

#ifndef SQ
double SQ(double);
#define SQ(x) ((x)*(x))
#endif

void DEMStats(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	** Intelligent comment goes here.
	**
	** NOTE: these are only particle-particle stats; no particle-wall.
	*/

	PARTICLE *pn = NULL;
	double r,dDist,dVx,dVy,dVz,dSpeed;
	int i,iClosest=0;
	
	mdlassert(smf->pkd->mdl,nSmooth > 0);
	
	r = RADIUS(p);
	mdlassert(smf->pkd->mdl,r > 0.0); /* paranoia check */
	p->DEMStats.fDistMin = FLOAT_MAXVAL;
	p->DEMStats.fSpeedMax = 0.0;
	p->DEMStats.fOverlap = 0.0;
	p->DEMStats.fS2 = 0.0;

#ifdef WALLS
	if (PARTICLE_STUCK(p))
		return; /* don't accumulate statistics for stuck particles */  /* ...stuck particles hopefully don't get here anyway [srs] */
#endif

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		if (pn->iOrder == p->iOrder)
			continue;
		dDist = sqrt(nnList[i].fDist2);
		/* find closest neighbor and distance */
		if (dDist < p->DEMStats.fDistMin) {
			iClosest = i;
			p->DEMStats.fDistMin = dDist;
			}
		{
			/* overlap statistics for histogram */
			double rn,sr,dOverlap=0.0;
			rn = RADIUS(pn);
			mdlassert(smf->pkd->mdl,rn > 0.0);
			/* avoiding the double count */
			if (rn > r)
				continue; /* skip if neighbor is bigger */
			if (rn == r && pn->iOrder < p->iOrder)
				continue; /* skip if neighbor is equal sized and has smaller iOrder */
			sr = r + rn;
			if (dDist <= sr) {
				/* overlap! */
				dOverlap = (r + rn - dDist)/(2.0*rn);
				if (dOverlap > p->DEMStats.fOverlap) {
					p->DEMStats.fOverlap = dOverlap; /* store if this is indeed the maximum overlap */
					/*if (dOverlap > 0.001) p->iColor = 4;*/ /*DEBUG:HACK TO RECOLOR LARGE OVERLAPS*/
				}
				/* calculate relative speed */
				dVx = p->v[0] - pn->v[0];
				dVy = p->v[1] - pn->v[1];
				dVz = p->v[2] - pn->v[2];
#ifdef SLIDING_PATCH
				if (smf->PP.bPatch) {
					/* adjust neighbor y-velocity if it's a ghost! */
					if (p->r[0] > pn->r[0] && nnList[i].dx < 0.0)
						dVy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
					else if (p->r[0] < pn->r[0] && nnList[i].dx > 0.0)
						dVy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
					}
#endif /* SLIDING PATCH */
				dSpeed = sqrt(dVx*dVx + dVy*dVy + dVz*dVz);
				/* find highest relative overlap speed */
				if (dSpeed > p->DEMStats.fSpeedMax)
					p->DEMStats.fSpeedMax = dSpeed;
				}
			}
		}
	/* calculate angle between the position vector and the relative velocity vector between two closest particles */
	pn = nnList[iClosest].pPart;
	/* find cos(alpha) */
	{
		dDist = p->DEMStats.fDistMin;
		if (dDist == 0.0)
			p->DEMStats.fCosA = -DBL_MAX;
		else {
			dVx = p->v[0] - pn->v[0];
			dVy = p->v[1] - pn->v[1];
			dVz = p->v[2] - pn->v[2];
#ifdef SLIDING_PATCH
			if (smf->PP.bPatch) {
				if (p->r[0] > pn->r[0] && nnList[i].dx < 0.0)
					dVy += 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
				else if (p->r[0] < pn->r[0] && nnList[i].dx > 0.0)
					dVy -= 1.5*smf->PP.dOrbFreq*smf->PP.dWidth;
				}
#endif /* SLIDING PATCH */
			dSpeed = sqrt(dVx*dVx + dVy*dVy + dVz*dVz);
			if (dSpeed == 0.0)
				p->DEMStats.fCosA = DBL_MAX;
			else {
				double dRdotV = nnList[iClosest].dx*dVx + nnList[iClosest].dy*dVy + nnList[iClosest].dz*dVz;
				p->DEMStats.fCosA = dRdotV/(dDist*dSpeed);
				}
			}
		}
	/* maximum S vector... */
	{
		DEM_ELEMENT *pe;
		double dS2;
		for (i=0;i<MAX_NUM_OVERLAPS_PER_PARTICLE;i++) {
			pe = &p->overlaps[i];
			if (pe->iOrder != -1) {
				dS2 = vectorMagSq(pe->vShear);
				if (dS2 > p->DEMStats.fS2)
					p->DEMStats.fS2 = dS2;
				}
			}
		}
	}

#endif /* DEM */

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

void initFindOverlaps(void *p)
{
	((PARTICLE *)p)->dtCol = 0.0;
	}

void combFindOverlaps(void *p1,void *p2)
{
	if (((PARTICLE *)p2)->dtCol < 0.0) ((PARTICLE *)p1)->dtCol = -1.0;
	}

void FindOverlaps(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	/*
	 ** Streamlined version of FindRejects() designed specifically
	 ** for the sliding patch when we want to randomize particle
	 ** data following an azimuthal boundary wrap.  As part of the
	 ** randomization procedure, we need to make sure we don't
	 ** overlap any particles.  That's what's checked for here.
	 */

	PARTICLE *pn = NULL;
	double r,rn,sr;
	int i;

	mdlassert(smf->pkd->mdl,nSmooth > 1); /* for now */

	mdlassert(smf->pkd->mdl,p->dtCol >= 0.0); /* can't already be rejected */

	mdlassert(smf->pkd->mdl,p->bAzWrap == 1); /* must have wrapped */

	r = RADIUS(p);

	for (i=0;i<nSmooth;i++) {
		pn = nnList[i].pPart;
		if (pn->iOrder == p->iOrder) continue;
		rn = RADIUS(pn);
		sr = r + rn;
		if (nnList[i].fDist2 <= sr*sr) p->dtCol = -1.0; /* cf REJECT() macro */
	}

	if (p->dtCol >= 0.0) p->bAzWrap = 0; /* if not rejected, do not need to regenerate */
	if (p->bAzWrap) printf("FindOverlaps(): particle %i overlaps particle %i\n",p->iOrder,pn->iOrder);
}
#endif /* SLIDING_PATCH */
