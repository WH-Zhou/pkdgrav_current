#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <stdio.h> /* needed for MacOS X to get xdrstdio_create() prototype */
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#include "pkd.h"
#include "ewald.h"
#include "grav.h"
#include "walk.h"
#include "opentype.h"
#include "mdl.h"
#include "tipsydefs.h"

#ifdef COLLISIONS
#include "ssdefs.h" /* in turn includes ssio.h */
#include "collision.h"
#endif

#ifdef SLIDING_PATCH
#include "random.h"
#endif

#ifdef AGGS
#include "aggs.h"
#endif

#ifdef JOHNNY
#include "linalg.h" /* for Vector type and operations */
#include "random.h"
#endif

#ifdef _LARGE_FILES
#define fseek fseeko
#endif

/*
 BalsaraSwitch hackery
*/
/*
#define VISCSWITCHCUT 0.6
*/

double pkdGetTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].sec);
	}

double pkdGetSystemTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].system_sec);
	}

double pkdGetWallClockTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].wallclock_sec);
	}


void pkdClearTimer(PKD pkd,int iTimer)
{
	int i;

	if (iTimer >= 0) {
		pkd->ti[iTimer].sec = 0.0;
		pkd->ti[iTimer].system_sec = 0.0;
		pkd->ti[iTimer].wallclock_sec = 0.0;
		pkd->ti[iTimer].iActive = 0;
		}
	else {
		for (i=0;i<MAX_TIMERS;++i) {
			pkd->ti[i].sec = 0.0;
			pkd->ti[i].system_sec = 0.0;
			pkd->ti[i].wallclock_sec = 0.0;
			pkd->ti[i].iActive = 0;
			}
		}
	}


void pkdStartTimer(PKD pkd,int iTimer)
{
	struct timeval tv;
	pkd->ti[iTimer].iActive++;

	if (pkd->ti[iTimer].iActive == 1) {
		pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
		gettimeofday(&tv,NULL);
		pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
#ifndef _CRAYMPP
		{
	    struct rusage ru;
	    
	    getrusage(0,&ru);
	    pkd->ti[iTimer].system_stamp = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
		}
#endif
		}
	}


void pkdStopTimer(PKD pkd,int iTimer)
{
	double sec;
	struct timeval tv;

	sec = -pkd->ti[iTimer].stamp;
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
	sec += pkd->ti[iTimer].stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].sec += sec;

	sec = -pkd->ti[iTimer].wallclock_stamp;
	gettimeofday( &tv, NULL );
	pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double)tv.tv_usec;
	sec += pkd->ti[iTimer].wallclock_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].wallclock_sec += sec;

#ifndef _CRAYMPP
	{
	struct rusage ru;

	sec = -pkd->ti[iTimer].system_stamp;
	getrusage(0,&ru);
	pkd->ti[iTimer].system_stamp = ((double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec);
	sec += pkd->ti[iTimer].system_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].system_sec += sec;
	}
#endif
	pkd->ti[iTimer].iActive--;
	}


void pkdInitialize(PKD *ppkd,MDL mdl,int iOrder,int nStore,int nLvl,
				   FLOAT *fPeriod,int nDark,int nGas,int nStar)
{
	PKD pkd;
	int j;
	
	pkd = (PKD)malloc(sizeof(struct pkdContext));
	mdlassert(mdl,pkd != NULL);
	pkd->mdl = mdl;
	pkd->iOrder = iOrder;
	pkd->idSelf = mdlSelf(mdl);
	pkd->nThreads = mdlThreads(mdl);
	pkd->nStore = nStore;
	pkd->nLocal = 0;
	pkd->nDark = nDark;
	pkd->nGas = nGas;
	pkd->nStar = nStar;
	pkd->nMaxOrderGas = nGas - 1;
	pkd->nMaxOrderDark = nGas + nDark - 1;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
	/*
	 ** Allocate the main particle store.
	 ** Need to use mdlMalloc() since the particles will need to be
	 ** visible to all other processors thru mdlAquire() later on.
	 **
	 ** We need one EXTRA storage location at the very end to use for 
	 ** calculating acceleration on arbitrary positions in space, for example
	 ** determining the force on the sun. The easiest way to do this is to
	 ** allocate one hidden particle, which won't interfere with the rest of
	 ** the code (hopefully). pkd->pStore[pkd->nStore] is this particle.
	 */
	pkd->pStore = mdlMalloc(pkd->mdl,(nStore+1)*sizeof(PARTICLE));
	mdlassert(pkd->mdl,pkd->pStore != NULL);
	pkd->kdNodes = NULL;
	pkd->piLeaf = NULL;
	pkd->kdTop = NULL;
	/*
	 ** Allocate initial interaction lists
	 */
	pkd->nMaxPart = 500;
	pkd->nMaxCellSoft = 500;
	pkd->nMaxCellNewt = 500;
	pkd->nSqrtTmp = 500;
	pkd->ilp = malloc(pkd->nMaxPart*sizeof(ILP));
	mdlassert(pkd->mdl,pkd->ilp != NULL);
	pkd->ilcs = malloc(pkd->nMaxCellSoft*sizeof(ILCS));
	mdlassert(pkd->mdl,pkd->ilcs != NULL);
	pkd->ilcn = malloc(pkd->nMaxCellNewt*sizeof(ILCN));
	mdlassert(pkd->mdl,pkd->ilcn != NULL);
	pkd->sqrttmp = malloc(pkd->nSqrtTmp*sizeof(double));
	mdlassert(pkd->mdl,pkd->sqrttmp != NULL);
	pkd->d2a = malloc(pkd->nSqrtTmp*sizeof(double));
	mdlassert(pkd->mdl,pkd->d2a != NULL);
	/*
	 ** Ewald stuff!
	 */
	pkd->nMaxEwhLoop = 100;
	pkd->ewt = malloc(pkd->nMaxEwhLoop*sizeof(EWT));
	mdlassert(pkd->mdl,pkd->ewt != NULL);
	/*
	 * Cooling
	 */
#ifdef GASOLINE
#ifndef NOCOOLING
	pkd->Cool = CoolInit();
#endif	
#endif	

	*ppkd = pkd;


	}


void pkdFinish(PKD pkd)
{
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	if (pkd->kdTop) free(pkd->kdTop);
	if (pkd->piLeaf) free(pkd->piLeaf);
	free(pkd->ilp);
	free(pkd->ilcs);
	free(pkd->ilcn);
	free(pkd->sqrttmp);
	free(pkd->d2a);
	free(pkd->ewt);
#if defined(GASOLINE) && !defined(NOCOOLING)
	CoolFinalize(pkd->Cool);
#endif
	mdlFree(pkd->mdl,pkd->pStore);
	free(pkd);
	}


void pkdSeek(PKD pkd,FILE *fp,int nStart,int bStandard)
{
	long lStart;

	/*
	 ** Seek according to true XDR size structures when bStandard is true.
	 ** This may be a bit dicey, but it should work as long
	 ** as no one changes the tipsy binary format!
	 */
	if (bStandard) lStart = 32;
	else lStart = sizeof(struct dump);
	if (nStart > pkd->nGas) {
		if (bStandard) lStart += pkd->nGas*48L;
		else lStart += pkd->nGas*sizeof(struct gas_particle);
		nStart -= pkd->nGas;
		if (nStart > pkd->nDark) {
			if (bStandard) lStart += pkd->nDark*36L;
			else lStart += pkd->nDark*sizeof(struct dark_particle);
			nStart -= pkd->nDark;
			if (bStandard) lStart += nStart*44L;
			else lStart += nStart*sizeof(struct star_particle);
			}
		else {
			if (bStandard) lStart += nStart*36L;
			else lStart += nStart*sizeof(struct dark_particle);
			}
		}
	else {
		if (bStandard) lStart += nStart*48L;
		else lStart += nStart*sizeof(struct gas_particle);
		}
	assert(fseek(fp,lStart,0) == 0);
	}


void pkdGenericSeek(PKD pkd,FILE *fp,int nStart,int iHeader,int iElement)
{
	long lStart;

	lStart = iHeader + nStart*iElement;
	fseek(fp,lStart,0);
	}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal,
				  int bStandard,int iReadIOrder,double dvFac,double dTuFac)
{
	FILE *fp,*fp2 = NULL;
#ifdef SIMPLESF
	FILE *fpmStar = NULL, *fptCoolAgain = NULL;
#endif
	int i,j, iSetMask;
	PARTICLE *p;
	struct dark_particle dp;
	struct gas_particle gp;
	struct star_particle sp;
	float fTmp;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization.
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
#ifdef GASOLINE
		p->u = 0.0;
		p->uPred = 0.0;
#ifdef STARFORM
		p->fESNrate = 0.0;
		p->fNSN = 0.0;
		p->fNSNtot = 0.0;
		p->fMOxygenOut = 0.0;
		p->fMIronOut = 0.0;
		p->fMFracOxygen = 0.0;
		p->fMFracIron = 0.0;
                p->fTimeCoolIsOffUntil = 0.0;
#endif
#ifdef SIMPLESF
		p->fMassStar = 0;
#endif
#ifdef SHOCKTRACK
		p->ShockTracker = 0.0;
#endif
#ifndef NOCOOLING		
		/* Place holders -- later fixed in pkdInitEnergy */
        CoolDefaultParticleData( &p->CoolParticle );
#endif
		p->c = 0.0;
		p->fMetals = 0.0;
		p->fTimeForm = 0.0;
#endif
#ifdef NEED_VPRED
		for (j=0;j<3;++j) {
			p->vPred[j] = 0.0;
			}
#endif
		}
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	mdlassert(pkd->mdl,fp != NULL);
	/*
	 ** Seek to right place in file.
	 */
	pkdSeek(pkd,fp,nStart,bStandard);

	/* Open iOrder file if requested */
	if (iReadIOrder) {
		char atmp[160];
		sprintf(atmp,"%s.iord",pszFileName);
		fp2 = fopen(atmp,"r");
		mdlassert(pkd->mdl,fp2 != NULL);
		/*
		 ** Seek to right place in file
		 */
		switch(iReadIOrder) {
		case 1:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(int));
			break;
		case 2:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(long long));
			break;
		case 3:
			pkdGenericSeek(pkd,fp2,nStart,sizeof(int),sizeof(pkd->pStore[0].iOrder));
			break;
		default:
			fprintf(stderr,"Don't understand iOrder format: %d\n",iReadIOrder);
			mdlassert(pkd->mdl,0);
			}
		}

#ifdef SIMPLESF
	{
		char atmp[160];
		sprintf(atmp,"%s.tCoolAgain",pszFileName);
		fptCoolAgain = fopen(atmp,"r");
		if (fptCoolAgain!=NULL) {
			/*
			 ** Seek to right place in file
			 */
			pkdGenericSeek(pkd,fptCoolAgain,nStart,sizeof(int),sizeof(float));
			}
		else {
			fprintf(stderr, "Could not open %s,  skipped.\n",atmp);
			}

		sprintf(atmp,"%s.mStar",pszFileName);
		fpmStar = fopen(atmp,"r");
		if (fpmStar!=NULL) {
			/*
			 ** Seek to right place in file
			 */
			pkdGenericSeek(pkd,fpmStar,nStart,sizeof(int),sizeof(float));
			}
		else {
			fprintf(stderr, "Could not open %s,  skipped.\n",atmp);
			}
		}
#endif

	/*
	 ** Read Stuff!
	 */
	if (bStandard) {
		FLOAT vTemp;
		long long LongTmp;
		int IntTmp;
		XDR xdrs;
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i; /* temporary */
			if (pkdIsDarkByOrder(pkd,p)) {
				iSetMask = TYPE_DARK;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass >= 0.0);
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
					}
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT				
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else if (pkdIsGasByOrder(pkd,p)) {
				iSetMask = TYPE_GAS;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				assert(p->fMass > 0.0);
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
#ifdef NEED_VPRED
					p->vPred[j] = dvFac*vTemp;
#endif
					}
#ifdef GASOLINE
				xdr_float(&xdrs,&fTmp);
				p->fDensity = fTmp;
				/*
				 ** Convert Temperature to Thermal energy.
				 */
				xdr_float(&xdrs,&fTmp);
				vTemp = fTmp;
				p->u = dTuFac*vTemp;
				p->uPred = dTuFac*vTemp;
#ifdef COOLDEBUG
				if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
				assert(p->u >= 0.0);
				assert(p->uPred >= 0.0);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
#else
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else if (pkdIsStarByOrder(pkd,p)) {
				iSetMask = TYPE_STAR;
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
#ifdef STARFORM
				p->fMassForm = fTmp;
#endif
				assert(p->fMass >= 0.0);
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
					}
#ifdef GASOLINE
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
				xdr_float(&xdrs,&fTmp);
				p->fTimeForm = fTmp;
#else
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
#ifdef CHANGESOFT
				p->fSoft0 = fTmp;
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else mdlassert(pkd->mdl,0);

			TYPESet(p,iSetMask);
			switch (iReadIOrder) {
			case 0:
				break;
			case 1:
				fread(&IntTmp,sizeof(IntTmp),1,fp2);
				p->iOrder = IntTmp;
				break;
			case 2:
				fread(&LongTmp,sizeof(LongTmp),1,fp2);
				p->iOrder = LongTmp;
				break;
			case 3:
				fread(&p->iOrder,sizeof(p->iOrder),1,fp2);
				break;
				}
#ifdef SIMPLESF
			if (fptCoolAgain!=NULL) {
				fread(&fTmp,sizeof(float),1,fptCoolAgain);
				if (pkdIsGasByOrder(pkd,p)) p->fTimeForm = fTmp;
				}
			if (fpmStar!=NULL) {
				fread(&fTmp,sizeof(float),1,fpmStar);
				if (pkdIsGasByOrder(pkd,p)) p->fMassStar = fTmp;
				}
#endif
			}
		xdr_destroy(&xdrs);
		}
	else {
		long long LongTmp;
		int IntTmp;
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i;

			if (pkdIsDarkByOrder(pkd,p)) {
				iSetMask = TYPE_DARK;
				fread(&dp,sizeof(struct dark_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = dp.pos[j];
					p->v[j] = dvFac*dp.vel[j];
					}
				p->fMass = dp.mass;
				assert(p->fMass >= 0.0);
				p->fSoft = dp.eps;
#ifdef CHANGESOFT
				p->fSoft0 = dp.eps;
#endif
				p->fPot = dp.phi;
				}
			else if (pkdIsGasByOrder(pkd,p)) {
				iSetMask = TYPE_GAS;
				fread(&gp,sizeof(struct gas_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = gp.pos[j];
					p->v[j] = dvFac*gp.vel[j];
#ifdef NEED_VPRED
					p->vPred[j] = dvFac*gp.vel[j];
#endif
					}
				p->fMass = gp.mass;
				assert(p->fMass >= 0.0);
				p->fSoft = gp.hsmooth;
#ifdef CHANGESOFT
				p->fSoft0 = gp.hsmooth;
#endif
				p->fPot = gp.phi;
#ifdef GASOLINE
				p->fDensity = gp.rho;
				p->u = dTuFac*gp.temp;
				p->uPred = dTuFac*gp.temp;
#ifdef COOLDEBUG
				if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
				assert(p->u >= 0.0);
				assert(p->uPred >= 0.0);
#endif
				p->fMetals = gp.metals;
#endif
				}
			else if (pkdIsStarByOrder(pkd,p)) {
				iSetMask = TYPE_STAR;
				fread(&sp,sizeof(struct star_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = sp.pos[j];
					p->v[j] = dvFac*sp.vel[j];
					}
				p->fMass = sp.mass;
#ifdef STARFORM
				p->fMassForm = sp.mass;
#endif
				assert(p->fMass >= 0.0);
				p->fSoft = sp.eps;
#ifdef CHANGESOFT
				p->fSoft0 = sp.eps;
#endif
				p->fPot = sp.phi;
#ifdef GASOLINE
				p->fMetals = sp.metals;
				p->fTimeForm = sp.tform;		
#endif
				}
			else mdlassert(pkd->mdl,0);

			TYPESet(p,iSetMask); /* needed to get max order info */
			switch (iReadIOrder) {
			case 0:
				p->iOrder = nStart + i;
				break;
			case 1:
				fread(&IntTmp,sizeof(IntTmp),1,fp2);
				p->iOrder = IntTmp;
				break;
			case 2:
				fread(&LongTmp,sizeof(LongTmp),1,fp2);
				p->iOrder = LongTmp;
				break;
			case 3:
				fread(&p->iOrder,sizeof(p->iOrder),1,fp2);
				break;
				}

#ifdef SIMPLESF
			if (fptCoolAgain!=NULL) {
				fread(&fTmp,sizeof(float),1,fptCoolAgain);
				if (pkdIsGasByOrder(pkd,p)) p->fTimeForm = fTmp;
				}
			if (fpmStar!=NULL) {
				fread(&fTmp,sizeof(float),1,fpmStar);
				if (pkdIsGasByOrder(pkd,p)) p->fMassStar = fTmp;
				}
#endif
			}
		}
	if (fp2!=NULL) fclose(fp2);
	fclose(fp);
	}


void pkdCalcBound(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive, BND *pbndBall)
{
        /* Faster by assuming active order */

	int i,j;
	FLOAT fBall,r;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		pbndBall->fMin[j] = FLOAT_MAXVAL;
		pbndBall->fMax[j] = -FLOAT_MAXVAL;
		}
	/*
	 ** Calculate Active Bounds assume TreeActive is the same as Active (Gnah!)
	 */
	for (i=0;i<pkd->nActive;++i) {
	        fBall = pkd->pStore[i].fBallMax;
	        for (j=0;j<3;++j) {
		        r = pkd->pStore[i].r[j];
			if (r < pbnd->fMin[j]) pbnd->fMin[j] = r;
			if (r > pbnd->fMax[j]) pbnd->fMax[j] = r;
	                if (r < pbndActive->fMin[j]) pbndActive->fMin[j] = r;
			if (r > pbndActive->fMax[j]) pbndActive->fMax[j] = r;
			if (r-fBall < pbndBall->fMin[j]) pbndBall->fMin[j] = r-fBall;
			if (r+fBall > pbndBall->fMax[j]) pbndBall->fMax[j] = r+fBall;
		        }
		}
	/*
	 ** Calculate Local Bounds.
	 */
	for (i=pkd->nActive;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
		        r = pkd->pStore[i].r[j];
			if (r < pbnd->fMin[j]) pbnd->fMin[j] = r;
			if (r > pbnd->fMax[j]) pbnd->fMax[j] = r;
			}
		}

	for (j=0;j<3;++j) {
		pbndTreeActive->fMin[j] = pbndActive->fMin[j];
		pbndTreeActive->fMax[j] = pbndActive->fMax[j];
		}
	}


void pkdCalcBound_old(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive, BND *pbndBall)
{
	int i,j;
	FLOAT fBall;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		pbndTreeActive->fMin[j] = FLOAT_MAXVAL;
		pbndTreeActive->fMax[j] = -FLOAT_MAXVAL;
		pbndBall->fMin[j] = FLOAT_MAXVAL;
		pbndBall->fMax[j] = -FLOAT_MAXVAL;
		}
	/*
	 ** Calculate Local Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
			if (pkd->pStore[i].r[j] < pbnd->fMin[j]) 
				pbnd->fMin[j] = pkd->pStore[i].r[j];
			if (pkd->pStore[i].r[j] > pbnd->fMax[j])
				pbnd->fMax[j] = pkd->pStore[i].r[j];
			}
		}
	/*
	 ** Calculate Active Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndActive->fMin[j]) 
					pbndActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndActive->fMax[j])
					pbndActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	/*
	 ** Calculate TreeActive Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryTREEACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndTreeActive->fMin[j]) 
					pbndTreeActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndTreeActive->fMax[j])
					pbndTreeActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	/*
	 ** Calculate fBall Bounds on TreeActive Particles
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryTREEACTIVE(&(pkd->pStore[i]))) {
		        fBall = pkd->pStore[i].fBallMax;
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j]-fBall < pbndBall->fMin[j]) 
					pbndBall->fMin[j] = pkd->pStore[i].r[j]-fBall;
				if (pkd->pStore[i].r[j]+fBall > pbndBall->fMax[j])
					pbndBall->fMax[j] = pkd->pStore[i].r[j]+fBall;
				}
			}
		}
	}


void pkdGasWeight(PKD pkd)
{
    PARTICLE *p;
    int i;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (TYPETest( p, TYPE_GAS )) p->fWeight=1.0;
		else p->fWeight=0.0;
		}
    }

void pkdRungDDWeight(PKD pkd, int iMaxRung, double dWeight)
{
    PARTICLE *p;
    int i;
    float fRungWeight[50],sum;

    mdlassert(pkd->mdl,iMaxRung<50);
    fRungWeight[0]=1.0;
    sum=1.0;
    for (i=1;i<=iMaxRung;i++) {
                sum*=2.0;
                fRungWeight[i] = dWeight* sum + (1-dWeight);
    }
  
    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		p->fWeight *= fRungWeight[p->iRung];
		}
    }

/*
 ** Partition particles between iFrom and iTo into those < fSplit and
 ** those >= to fSplit.  Find number and weight in each partition.
 */
int pkdWeight(PKD pkd,int d,FLOAT fSplit,int iSplitSide,int iFrom,int iTo,
			 int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
{
	int i,iPart;
	FLOAT fLower,fUpper;

	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	/*
	 ** Calculate the lower weight and upper weight BETWEEN the particles
	 ** iFrom to iTo!
	 */
	fLower = 0.0;
	for (i=iFrom;i<iPart;++i) {
		fLower += pkd->pStore[i].fWeight;
		}
	fUpper = 0.0;
	for (i=iPart;i<=iTo;++i) {
		fUpper += pkd->pStore[i].fWeight;
		}
	if (iSplitSide) {
		*pfLow = fUpper;
		*pfHigh = fLower;
		}
	else {
		*pfLow = fLower;
		*pfHigh = fUpper;
		}
	return(iPart);
	}


/*
 ** Partition particles between iFrom and iTo into those < fSplit and
 ** those >= to fSplit.  Find number and weight in each partition.
 */
int pkdWeightWrap(PKD pkd,int d,FLOAT fSplit,FLOAT fSplit2, int iSplitSide,int iFrom,int iTo,
			  int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
{
	int i,iPart;
	FLOAT fLower,fUpper;

	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (!iSplitSide) {
		iPart = pkdLowerPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	else {
		iPart = pkdUpperPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
		*pnHigh = iPart;
		*pnLow = pkdLocal(pkd)-iPart;
		}
	/*
	 ** Calculate the lower weight and upper weight BETWEEN the particles
	 ** iFrom to iTo!
	 */
	/* Not needed */
	fLower = 0.0;
	for (i=iFrom;i<iPart;++i) {
		fLower += pkd->pStore[i].fWeight;
		}
	fUpper = 0.0;
	for (i=iPart;i<=iTo;++i) {
		fUpper += pkd->pStore[i].fWeight;
		}
	if (!iSplitSide) {
		*pfLow = fLower;
		*pfHigh = fUpper;
		}
	else {
		*pfLow = fUpper;
		*pfHigh = fLower;
		}

	return(iPart);
	}


int pkdOrdWeight(PKD pkd,int iOrdSplit,int iSplitSide,int iFrom,int iTo,
				 int *pnLow,int *pnHigh)
{
	int iPart;

	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	return(iPart);
	}

int pkdLowerPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;
	/*
	int i0=i,j0=j,iold;
	*/
	mdlassert(pkd->mdl,d < 3);

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && pkd->pStore[i].r[d] >= fSplit)
		;
        while (i<--j && pkd->pStore[j].r[d] < fSplit)
		;
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (pkd->pStore[++i].r[d] >= fSplit);
	     while (pkd->pStore[--j].r[d] < fSplit);
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*     iold=i;
     i=i0;
     j=j0;
	
	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] >= fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] < fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
	
	
 done:
    mdlassert(pkd->mdl,i==iold);
	*/
    return(i);
	}


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;
	/*	
	int i0=i,j0=j,iold;
	*/

	mdlassert(pkd->mdl,d < 3);

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && pkd->pStore[i].r[d] < fSplit);
	while (i<--j && pkd->pStore[j].r[d] >= fSplit);
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (pkd->pStore[++i].r[d] < fSplit);
	     while (pkd->pStore[--j].r[d] >= fSplit);
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*
     iold=i;
     i=i0;
     j=j0;
	
	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] < fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] >= fSplit)
            if (i > --j) goto done;

		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    mdlassert(pkd->mdl,i==iold);
	*/
    return(i);
	}


int pkdLowerPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int i,int j)
{
	PARTICLE pTemp;

	mdlassert(pkd->mdl,d < 3);

	if (fSplit1 > fSplit2) {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;

	     while (1) {
	       while (pkd->pStore[++i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1);
	       while (pkd->pStore[--j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	else {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;

	     while (1) {
	       while (pkd->pStore[++i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1);
	       while (pkd->pStore[--j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	

 done1:
	return(i);
        }

int pkdUpperPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int i,int j)
{
	PARTICLE pTemp;

	mdlassert(pkd->mdl,d < 3);

	if (fSplit1 > fSplit2) {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1);
	       while (pkd->pStore[--j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	else {
	     if (i > j) goto done1;
	     i--;
	     j++;
	     while (++i<j && (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1));
	     while (i<--j && (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1));
	     if (i>=j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
	     
	     while (1) {
	       while (pkd->pStore[++i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1);
	       while (pkd->pStore[--j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1);
	       if (i > j) goto done1;
	       pTemp = pkd->pStore[i];
	       pkd->pStore[i] = pkd->pStore[j];
	       pkd->pStore[j] = pTemp;
               }
	     }
	

 done1:
	return(i);
        }



int pkdLowerOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder >= nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder < nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdUpperOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder < nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder >= nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}

int pkdActiveTypeOrder(PKD pkd, unsigned int iTestMask)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;
	/*
	int i0=i,j0=j,iold;
	*/

	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && TYPETest(&(pkd->pStore[i]), iTestMask ))
		;
	while (i<--j && !TYPETest(&(pkd->pStore[j]), iTestMask ))
		;
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (TYPETest(&(pkd->pStore[++i]), iTestMask ));
	     while (!TYPETest(&(pkd->pStore[--j]), iTestMask ));
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
 done1:
	/*	iold=i;
	i=i0;
	j=j0;
	
        if (i > j) goto done;
    while (1) {
        while (TYPETest(&(pkd->pStore[i]), iTestMask ))
            if (++i > j) goto done;
        while (!TYPETest(&(pkd->pStore[j]), iTestMask ))
            if (i > --j) goto done;
                pTemp = pkd->pStore[i];
                pkd->pStore[i] = pkd->pStore[j];
                pkd->pStore[j] = pTemp;
        }
	
 done:
	assert (i==iold);
	*/
	if ( iTestMask & TYPE_ACTIVE )       pkd->nActive = i;
	if ( iTestMask & TYPE_TREEACTIVE )   pkd->nTreeActive = i;
	if ( iTestMask & TYPE_SMOOTHACTIVE ) pkd->nSmoothActive = i;

	return (i);
	}


int pkdActiveOrder(PKD pkd)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;
	/*
	int i0=i,j0=j,iold;
	*/
	
	if (i > j) goto done1;
	i--;
	j++;
	while (++i<j && TYPEQueryACTIVE(&(pkd->pStore[i])))
		;
	while (i<--j && !TYPEQueryACTIVE(&(pkd->pStore[j])))
		;
	if (i>=j) goto done1;
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;

	while (1) {
	     while (TYPEQueryACTIVE(&(pkd->pStore[++i])))
			 ;
	     while (!TYPEQueryACTIVE(&(pkd->pStore[--j])))
			 ;
	     if (i > j) goto done1;
	     pTemp = pkd->pStore[i];
	     pkd->pStore[i] = pkd->pStore[j];
	     pkd->pStore[j] = pTemp;
             }
	
 done1:
	/*
	iold=i;
	i=i0;
	j=j0;

        if (i > j) goto done;
    while (1) {
        while (TYPEQueryACTIVE(&(pkd->pStore[i])))
            if (++i > j) goto done;
        while (!TYPEQueryACTIVE(&(pkd->pStore[j])))
            if (i > --j) goto done;
                pTemp = pkd->pStore[i];
                pkd->pStore[i] = pkd->pStore[j];
                pkd->pStore[j] = pTemp;
        }
	
 done:
        assert (i==iold);
	*/
	return (pkd->nActive = i);
	}


int pkdColRejects_Active_Inactive(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide)
{
	PARTICLE pTemp;
	int nSplit,nSplitInactive,iRejects,i,j;

	mdlassert(pkd->mdl,pkd->nRejects == 0);
	if (iSplitSide) {
		nSplit = pkdLowerPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	else {
		nSplit = pkdUpperPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	if (iSplitSide) {
		nSplitInactive = pkdLowerPart(pkd,d,fSplitInactive,
									  pkdActive(pkd),pkdLocal(pkd)-1);
		}
	else {
		nSplitInactive = pkdUpperPart(pkd,d,fSplitInactive,
									  pkdActive(pkd),pkdLocal(pkd)-1);
		}
	/*
	for(i = 0; i < nSplit; ++i)
	    mdlassert(pkd->mdl,TYPEQueryACTIVE(&(pkd->pStore[i])));
	for(i = pkdActive(pkd); i < nSplitInactive; ++i)
	    mdlassert(pkd->mdl,!TYPEQueryACTIVE(&(pkd->pStore[i])));
	*/

	nSplitInactive -= pkdActive(pkd);
	/*
	 ** Now do some fancy rearrangement.
	 */
	i = nSplit;
	j = nSplit+nSplitInactive;
	while (j < pkdActive(pkd) + nSplitInactive) {
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
		++i; ++j;
		}
	pkd->nRejects = pkdLocal(pkd) - nSplit - nSplitInactive;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nActive = nSplit;
	pkd->nLocal = nSplit + nSplitInactive;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int pkdColRejects(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide)
{
	int nSplit,iRejects,i;

	mdlassert(pkd->mdl,pkd->nRejects == 0);
	if (!iSplitSide) {
		nSplit = pkdLowerPartWrap(pkd,d,fSplitInactive,fSplit,0,pkdLocal(pkd)-1);
		}
	else {
		nSplit = pkdUpperPartWrap(pkd,d,fSplitInactive,fSplit,0,pkdLocal(pkd)-1);
		}

	pkd->nRejects = pkdLocal(pkd) - nSplit;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int pkdSwapRejects(PKD pkd,int idSwap)
{
	size_t nBuf;
	size_t nOutBytes,nSndBytes,nRcvBytes;

	if (idSwap != -1) {
		nBuf = (pkdSwapSpace(pkd))*sizeof(PARTICLE);
		nOutBytes = pkd->nRejects*sizeof(PARTICLE);
		mdlassert(pkd->mdl,pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
		mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkdLocal(pkd)],
				nOutBytes,&nSndBytes,&nRcvBytes);
		pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
		pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
		}
	return(pkd->nRejects);
	}

void pkdSwapAll(PKD pkd, int idSwap)
{
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;
    int bGood;
    
    /*
     ** Move particles to High memory.
     */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    bGood = mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
			&nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,bGood);
    mdlassert(pkd->mdl,nSndBytes/sizeof(PARTICLE) == (size_t) pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/sizeof(PARTICLE);
    }

int pkdSwapSpace(PKD pkd)
{
	return(pkdFreeStore(pkd) - pkdLocal(pkd));
	}


int pkdFreeStore(PKD pkd)
{
	return(pkd->nStore);
	}

int pkdActive(PKD pkd)
{
	return(pkd->nActive);
	}

int pkdTreeActive(PKD pkd)
{
	return(pkd->nTreeActive);
	}

int pkdSmoothActive(PKD pkd)
{
	return(pkd->nSmoothActive);
	}

int pkdInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nActive);
	}

int pkdTreeInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nTreeActive);
	}

#ifndef pkdLocal
int pkdLocal(PKD pkd)
{
	return(pkd->nLocal);
	}
#endif

void pkdTotals(PKD pkd, int *nDark, int *nGas, int *nStar)
{
    int i,nd,ng,ns;
    nd=ng=ns=0;
    for (i=0;i<pkd->nLocal;++i) {
        if (pkdIsDark(pkd,&pkd->pStore[i])) nd++;
        if (pkdIsGas(pkd,&pkd->pStore[i])) ng++;
        if (pkdIsStar(pkd,&pkd->pStore[i])) ns++;
        }
    *nDark=nd;
    *nGas=ng;
    *nStar=ns;
    }

int pkdNodes(PKD pkd)
{
	return(pkd->nNodes);
	}


void pkdDomainColor(PKD pkd)
{
#ifdef COLORCODE
	int i;
	
	for (i=0;i<pkd->nLocal;++i) {
		pkd->pStore[i].fColor = (FLOAT)pkd->idSelf;
		}
#endif
	}


int pkdColOrdRejects(PKD pkd,int nOrdSplit,int iSplitSide)
{
	int nSplit,iRejects,i;

	if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	pkd->nRejects = pkdLocal(pkd) - nSplit;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}

int cmpParticles(const void *pva,const void *pvb)
{
	const PARTICLE *pa = (const PARTICLE *)pva;
	const PARTICLE *pb = (const PARTICLE *)pvb;

	return(pa->iOrder - pb->iOrder);
	}

void pkdLocalOrder(PKD pkd)
{
	qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
	}

#ifdef FAST_AGGS
int cmpAggs(const void *pva,const void *pvb)	//JVD
{
	const PARTICLE *pa = (const PARTICLE *)pva;
	const PARTICLE *pb = (const PARTICLE *)pvb;

	return(pb->iOrgIdx - pa->iOrgIdx);	/* Reversed for decreasing order */
	}

void pkdLocalAggOrder(PKD pkd)	//JVD
{
	qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpAggs);

	}
#endif /* FAST_AGGS */

void pkdWriteTipsy(PKD pkd,char *pszFileName,int nStart,
				   int bStandard,double dvFac,double duTFac,int iGasModel)
{
	PARTICLE *p;
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	struct gas_particle gp;
	struct star_particle sp;
	int nout;
	float fTmp;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	mdlassert(pkd->mdl,fp != NULL);
	pkdSeek(pkd,fp,nStart,bStandard);
	if (bStandard) {
		FLOAT vTemp;
		XDR xdrs;
		/* 
		 ** Write Stuff!
		 */
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (pkdIsDark(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
#ifdef CHANGESOFT
				fTmp = p->fSoft0;
#else
				fTmp = p->fSoft;
#endif
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else if (pkdIsGas(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
				fTmp = p->fDensity;
				xdr_float(&xdrs,&fTmp);
#ifdef GASOLINE
				/*
				 ** Convert thermal energy to tempertature.
				 */
				switch (iGasModel) {
				case GASMODEL_COOLING:
#ifndef NOCOOLING
					vTemp = CoolCodeEnergyToTemperature( pkd->Cool, &p->CoolParticle, p->u );
#else
					mdlassert(pkd->mdl,0);
#endif
					break;
				default:
					vTemp = duTFac*p->u;
					}
				fTmp = vTemp;
				xdr_float(&xdrs,&fTmp);
				/* fTmp = sqrt(0.25*p->fBall2);  Write softening in tipsy outputs */
#ifdef CHANGESOFT
				fTmp = p->fSoft0;
#else
				fTmp = p->fSoft;
#endif
				xdr_float(&xdrs,&fTmp);
#ifdef DEBUG
				/* Store divv in metals for now */
				fTmp = p->divv;
#else
				fTmp = p->fMetals;
#endif
				xdr_float(&xdrs,&fTmp);
#else /* not gasoline */
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
#ifdef CHANGESOFT
				fTmp = p->fSoft0;
#else
				fTmp = p->fSoft;
#endif
				xdr_float(&xdrs,&fTmp);
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
#endif
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else if (pkdIsStar(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
#ifdef GASOLINE
				fTmp = p->fMetals;
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fTimeForm;
				xdr_float(&xdrs,&fTmp);
#else
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);			
#endif
#ifdef CHANGESOFT
				fTmp = p->fSoft0;
#else
				fTmp = p->fSoft;
#endif
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else mdlassert(pkd->mdl,0);
			}
		xdr_destroy(&xdrs);
		}
	else {
		/* 
		 ** Write Stuff!
		 */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (pkdIsDark(pkd,p)) {
				for (j=0;j<3;++j) {
					dp.pos[j] = p->r[j];
					dp.vel[j] = dvFac*p->v[j];
					}
				dp.mass = p->fMass;
#ifdef CHANGESOFT
 				dp.eps = p->fSoft0;
#else
 				dp.eps = p->fSoft;
#endif
				dp.phi = p->fPot;
				nout = fwrite(&dp,sizeof(struct dark_particle),1,fp);
				mdlassert(pkd->mdl,nout == 1);
				}
			else if (pkdIsGas(pkd,p)) {
				for (j=0;j<3;++j) {
					gp.pos[j] = p->r[j];
					gp.vel[j] = dvFac*p->v[j];
					}
				gp.mass = p->fMass;
#ifdef CHANGESOFT
				gp.hsmooth = p->fSoft0;
#else
				gp.hsmooth = p->fSoft;
#endif
				gp.phi = p->fPot;
				gp.rho = p->fDensity;
#ifdef GASOLINE
				gp.temp = duTFac*p->u;
				gp.metals = p->fMetals;
#else
				gp.temp = 0.0;
				gp.metals = 0.0;
#endif
				nout = fwrite(&gp,sizeof(struct gas_particle),1,fp);
				mdlassert(pkd->mdl,nout == 1);
				}
			else if (pkdIsStar(pkd,p)) {
				for (j=0;j<3;++j) {
					sp.pos[j] = p->r[j];
					sp.vel[j] = dvFac*p->v[j];
					}
				sp.mass = p->fMass;
#ifdef CHANGESOFT
				sp.eps = p->fSoft0;
#else
				sp.eps = p->fSoft;
#endif
				sp.phi = p->fPot;
#ifdef GASOLINE
				sp.metals = p->fMetals;
				sp.tform = p->fTimeForm;
#else
				sp.metals = 0.0;
				sp.tform = 0.0;
#endif
				nout = fwrite(&sp,sizeof(struct star_particle),1,fp);
				mdlassert(pkd->mdl,nout == 1);
				}
			else mdlassert(pkd->mdl,0);
			}
		}
	nout = fclose(fp);
	mdlassert(pkd->mdl,nout == 0);
	}


void pkdSetSoft(PKD pkd,double dSoft)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef CHANGESOFT
		p[i].fSoft0 = dSoft;
#endif
		p[i].fSoft = dSoft;
		}
	}

#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	assert(dFac > 0.0);
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) { 
		  /*		  assert(p[i].fSoft0 > 0.);
				  printf("Compton:  %f\n",p[i].fSoft0);*/
                  
		        assert(p[i].fSoft0 > 0.0);
	                p[i].fSoft = p[i].fSoft0*dFac;
			assert(p[i].fSoft > 0.0);
		        }
	        }
	else {
 	        assert(dSoftMax > 0.0);
	        for (i=0;i<n;++i) {
	                p[i].fSoft = p[i].fSoft0*dFac;
	                if (p[i].fSoft > dSoftMax) p[i].fSoft = dSoftMax;
		        }
	        }
	}

void pkdPreVariableSoft(PKD pkd,int iVariableSoftType)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);

#ifndef DENSSOFT
	for (i=0;i<n;++i) {
	        if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
			  p[i].fSoft = p[i].fBall2;
			  }
	        }
#endif
	}

void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul,int iVariableSoftType)
{
	PARTICLE *p;
	int i,n;
	double dTmp;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	
	if (bSoftMaxMul) {
	        for (i=0;i<n;++i) {
	                if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
#ifdef DENSSOFT
	                          p[i].fSoft = pow((p[i].fMass*1.90986/p[i].fDensity),.3333333333);
#else
		                  dTmp = sqrt(p[i].fBall2*.25);
	                          p[i].fBall2 = p[i].fSoft;
	                          p[i].fSoft = (dTmp <= p[i].fSoft0*dSoftMax ? dTmp : p[i].fSoft0*dSoftMax);
#endif
                                  }
		        }
	        }
	else {
	        for (i=0;i<n;++i) {
#ifdef CHECKSOFT			  
			  if (p[i].iOrder == CHECKSOFT) fprintf(stderr,"Particle %i: %g %g %g %i %g\n",p[i].iOrder,p[i].fDensity,sqrt(p[i].fBall2*.25),p[i].fSoft,p[i].iActive,dSoftMax);
#endif
			  if (TYPETest(&(p[i]),iVariableSoftType) && TYPEQueryACTIVE(&(p[i]))) {
#ifdef DENSSOFT
				  p[i].fSoft = pow((p[i].fMass*1.90986/p[i].fDensity),.3333333333);
#else
				  dTmp = sqrt(p[i].fBall2*.25);
				  p[i].fBall2 = p[i].fSoft;
				  p[i].fSoft = (dTmp <= dSoftMax ? dTmp : dSoftMax);
#endif
			          }
#ifdef CHECKSOFT			  
			  if (p[i].iOrder == CHECKSOFT) fprintf(stderr,"Particle %iA: %g %g %g %i %g\n",p[i].iOrder,p[i].fDensity,sqrt(p[i].fBall2*.25),p[i].fSoft,p[i].iActive,dSoftMax);
#endif
		        }
	        }
	}
#endif

void pkdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];

		if (p2->bndBall.fMin[j] < p1->bndBall.fMin[j])
			pOut->bndBall.fMin[j] = p2->bndBall.fMin[j];
		else
			pOut->bndBall.fMin[j] = p1->bndBall.fMin[j];
		if (p2->bndBall.fMax[j] > p1->bndBall.fMax[j])
			pOut->bndBall.fMax[j] = p2->bndBall.fMax[j];
		else
			pOut->bndBall.fMax[j] = p1->bndBall.fMax[j];
		}
	/*
	 ** Find the center of mass and mass weighted softening.
	 */
        pOut->fMass = p1->fMass + p2->fMass;
	pOut->fSoft = p1->fMass*p1->fSoft + p2->fMass*p2->fSoft;
	for (j=0;j<3;++j) {
		pOut->r[j] = p1->fMass*p1->r[j] + p2->fMass*p2->r[j];
		}
	if (pOut->fMass > 0) {
		pOut->fSoft /= pOut->fMass;
		for (j=0;j<3;++j) {
			pOut->r[j] /= pOut->fMass;
			}
		}
	}


void pkdCalcCell(PKD pkd,KDN *pkdn,FLOAT *rcm,int iOrder,
				 struct pkdCalcCellStruct *pcc)
{
	int pj;
	double m,dx,dy,dz,d2,d1;
	struct pkdCalcCellStruct cc;

	/*
	 ** Initialize moments.
	 ** Initialize various B numbers.
	 */
	switch (iOrder) {	
	case 4:
		cc.Hxxxx = 0.0;
		cc.Hxyyy = 0.0;
		cc.Hxxxy = 0.0;
		cc.Hyyyy = 0.0;
		cc.Hxxxz = 0.0;
		cc.Hyyyz = 0.0;
		cc.Hxxyy = 0.0;
		cc.Hxxyz = 0.0;
		cc.Hxyyz = 0.0;
		cc.Hxxzz = 0.0;
		cc.Hxyzz = 0.0;
		cc.Hxzzz = 0.0;
		cc.Hyyzz = 0.0;
		cc.Hyzzz = 0.0;
		cc.Hzzzz = 0.0;
		cc.B6 = 0.0;
		/*FALLTHRU*/
	case 3:
		cc.Oxxx = 0.0;
		cc.Oxyy = 0.0;
		cc.Oxxy = 0.0;
		cc.Oyyy = 0.0;
		cc.Oxxz = 0.0;
		cc.Oyyz = 0.0;
		cc.Oxyz = 0.0;
		cc.Oxzz = 0.0;
		cc.Oyzz = 0.0;
		cc.Ozzz = 0.0;
		cc.B5 = 0.0;
		/*FALLTHRU*/
	default:
		cc.Qxx = 0.0;
		cc.Qyy = 0.0;
		cc.Qzz = 0.0;
		cc.Qxy = 0.0;
		cc.Qxz = 0.0;
		cc.Qyz = 0.0;
		cc.B2 = 0.0;
		cc.B3 = 0.0;
		cc.B4 = 0.0;
		}
	cc.Bmax = 0.0;
	/*
	 ** Calculate moments and B numbers about center-of-mass.
	 */
	if (pkdn == NULL) pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - rcm[0];
		dy = pkd->pStore[pj].r[1] - rcm[1];
		dz = pkd->pStore[pj].r[2] - rcm[2];
		d2 = dx*dx + dy*dy + dz*dz;
		d1 = sqrt(d2);
		if (d1 > cc.Bmax) cc.Bmax = d1;
		cc.B2 += m*d2;
		cc.B3 += m*d2*d1;
		cc.B4 += m*d2*d2;
		cc.B5 += m*d2*d2*d1;
		cc.B6 += m*d2*d2*d2;
#ifdef COMPLETE_LOCAL
		d2 = 0.0;
#endif
		switch (iOrder) {
		case 4:
			/*
			 ** Calculate reduced hexadecapole moment...
			 */
			cc.Hxxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
			cc.Hxyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
			cc.Hxxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
			cc.Hyyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
			cc.Hxxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
			cc.Hyyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
			cc.Hxxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
			cc.Hxxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
			cc.Hxyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
			cc.Hxxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
			cc.Hxyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
			cc.Hxzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
			cc.Hyyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
			cc.Hyzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
			cc.Hzzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
			/*FALLTHRU*/
		case 3:
			/*
			 ** Calculate reduced octopole moment...
			 */
			cc.Oxxx += m*(dx*dx*dx - 0.6*d2*dx);
			cc.Oxyy += m*(dx*dy*dy - 0.2*d2*dx);
			cc.Oxxy += m*(dx*dx*dy - 0.2*d2*dy);
			cc.Oyyy += m*(dy*dy*dy - 0.6*d2*dy);
			cc.Oxxz += m*(dx*dx*dz - 0.2*d2*dz);
			cc.Oyyz += m*(dy*dy*dz - 0.2*d2*dz);
			cc.Oxyz += m*dx*dy*dz;
			cc.Oxzz += m*(dx*dz*dz - 0.2*d2*dx);
			cc.Oyzz += m*(dy*dz*dz - 0.2*d2*dy);
			cc.Ozzz += m*(dz*dz*dz - 0.6*d2*dz);;
			/*FALLTHRU*/
		default:
			/*
			 ** Calculate quadrupole moment...
			 */
			cc.Qxx += m*dx*dx;
			cc.Qyy += m*dy*dy;
			cc.Qzz += m*dz*dz;
			cc.Qxy += m*dx*dy;
			cc.Qxz += m*dx*dz;
			cc.Qyz += m*dy*dz;
			}
		}
	*pcc = cc;
	}


double fcnAbsMono(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (3.0*pkdn->mom.B2 - 2.0*pkdn->mom.B3/r)/t;
	return t;
	}


double fcnAbsQuad(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (4.0*pkdn->mom.B3 - 3.0*pkdn->mom.B4/r)/t;
	return t;
	}


double fcnAbsOct(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (5.0*pkdn->mom.B4 - 4.0*pkdn->mom.B5/r)/t;
	return t;
	}


double fcnAbsHex(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (6.0*pkdn->mom.B5 - 5.0*pkdn->mom.B6/r)/t;
	return t;
	}


double dRootBracket(KDN *pkdn,double dErrBnd,double (*fcn)(KDN *,double))
{
	double dLower,dUpper,dMid;
	double dLowerErr,dErrDif,dErrCrit;
	int iter = 0;

	dErrCrit = 1e-6*dErrBnd;
	dLower = (1.0 + 1e-6)*pkdn->mom.Bmax;
	dLowerErr = (*fcn)(pkdn,dLower);
	if (dLowerErr < dErrBnd) {
/*
		printf("iter:%d %g\n",iter,dLower);
*/
		return dLower;
		}
	/*
	 ** Hunt phase for upper bound.
	 */
	dUpper = 2*dLower;
	while (1) {
		dErrDif = dErrBnd - (*fcn)(pkdn,dUpper);
		if (dErrDif > dErrCrit) break;
		dUpper = 2*dUpper;
		}
	/*
	 ** Bracket root.
	 */
	while (1) {
		dMid = 0.5*(dLower + dUpper);
		dErrDif = dErrBnd - (*fcn)(pkdn,dMid);
		if (dErrDif < 0) {
			dLower = dMid;
			}
		else {
			dUpper = dMid;
			if (dErrDif < dErrCrit) break;
			}
		if (++iter > 32) break;
		}
/*
	printf("iter:%d %g\n",iter,dMid);
*/
	return dMid;
	}


double pkdCalcOpen(KDN *pkdn,int iOpenType,double dCrit,int iOrder)
{
	double dOpen=0;

	if (iOpenType == OPEN_ABSPAR) {
		switch (iOrder) {
		case 1:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsMono);
			break;
		case 2:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsQuad);
			break;
		case 3:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsOct);
			break;
		case 4:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsHex);
			break;
		default:
			assert(0);
			}
		}
	else if (iOpenType == OPEN_JOSH) {
		/*
		 ** Set opening criterion to an approximation of Josh's theta.
		 */
		dOpen = 2/sqrt(3.0)*pkdn->mom.Bmax/dCrit;
		if (dOpen < pkdn->mom.Bmax) dOpen = pkdn->mom.Bmax;
		}
	else {
		/*
		 ** Set opening criterion to the minimal, i.e., Bmax.
		 */
		dOpen = pkdn->mom.Bmax;
		}
	return(dOpen);
	}


void pkdUpPass(PKD pkd,int iCell,int iOpenType,double dCrit,
			   int iOrder, int bGravity)
{
	KDN *c;
	PARTICLE *p;
	int l,u,pj,j;
	double dOpen;

	c = pkd->kdNodes;
	p = pkd->pStore;
	l = c[iCell].pLower;
	u = c[iCell].pUpper;
	if (c[iCell].iDim >= 0) {
		pkdUpPass(pkd,LOWER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdUpPass(pkd,UPPER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdCombine(&c[LOWER(iCell)],&c[UPPER(iCell)],&c[iCell]);
		}
	else {
		c[iCell].fMass = 0.0;
		c[iCell].fSoft = 0.0;
		for (j=0;j<3;++j) {
			/*
			 ** We initialize the bounds to these extreme crazy values so
			 ** that any particle will set the bounds correctly, if there 
			 ** are no particles in the loop then these remain the bounds.
			 */
			c[iCell].bnd.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bnd.fMax[j] = -FLOAT_MAXVAL;
			c[iCell].bndBall.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bndBall.fMax[j] = -FLOAT_MAXVAL;
			c[iCell].r[j] = 0.0;
			}
		for (pj=l;pj<=u;++pj) {
			for (j=0;j<3;++j) {
				if (p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = p[pj].r[j];
				if (p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = p[pj].r[j];

				if (p[pj].r[j]-p[pj].fBallMax < c[iCell].bndBall.fMin[j])
					c[iCell].bndBall.fMin[j] = p[pj].r[j]-p[pj].fBallMax;
				if (p[pj].r[j]+p[pj].fBallMax > c[iCell].bndBall.fMax[j])
					c[iCell].bndBall.fMax[j] = p[pj].r[j]+p[pj].fBallMax;
				}
			/*
			 ** Find center of mass and total mass and mass weighted softening.
			 */
			c[iCell].fMass += p[pj].fMass;
			c[iCell].fSoft += p[pj].fMass*p[pj].fSoft;
			for (j=0;j<3;++j) {
				c[iCell].r[j] += p[pj].fMass*p[pj].r[j];
				}
			}
		if (c[iCell].fMass > 0) {
			for (j=0;j<3;++j) {
				c[iCell].r[j] /= c[iCell].fMass;
				}
			c[iCell].fSoft /= c[iCell].fMass;
			}
		}
	/*
	 ** Calculate multipole moments.
	 */
	if(bGravity) {
	    pkdCalcCell(pkd,&c[iCell],c[iCell].r,iOrder,&c[iCell].mom);
	    dOpen = pkdCalcOpen(&c[iCell],iOpenType,dCrit,iOrder);
	    c[iCell].fOpen2 = dOpen*dOpen;
	    }
	}


/*
 ** JST's Select
 */
void pkdSelect(PKD pkd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	FLOAT v;
	int i,j;

	p = pkd->pStore;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


int NumBinaryNodes(PKD pkd,int nBucket,int pLower,int pUpper)
{
	BND bnd;
	int i,j,m,d,l,u;
	FLOAT fSplit;
	int bGoodBounds;	/* is the cell a finite size? */
	PARTICLE *pStore = pkd->pStore;

	if (pLower > pUpper) return(0);
	else {
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pStore[pLower].r[j];
			bnd.fMax[j] = pStore[pLower].r[j];
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pStore[i].r[j];
				else if (pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pStore[i].r[j];
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (bnd.fMax[j] > bnd.fMin[j])
			    bGoodBounds = 1;
			}
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (bnd.fMax[j]-bnd.fMin[j] > bnd.fMax[d]-bnd.fMin[d]) {
					d = j;
					}
				}
			/*
			 ** Now we do the split.
			 */
			fSplit = 0.5*(bnd.fMin[d]+bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,fSplit,pLower,pUpper);
			/*
			 ** Recursive node count.
			 */
			l = NumBinaryNodes(pkd,nBucket,pLower,m-1);
			u = NumBinaryNodes(pkd,nBucket,m,pUpper);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			mdlassert(pkd->mdl,l > 0 && u > 0);
			return(l+u+1);
			}
		else {
			return(1);
			}
		}
	}


int BuildBinary(PKD pkd,int nBucket,int pLower,int pUpper,int iOpenType,
		double dCrit,int iOrder, int bGravity)
{
	KDN *pkdn;
	int i,j,m,d,c;
	FLOAT fm;
	double dOpen;
	int bGoodBounds;	/* Is the cell a finite size? */
	BND bnd, bndBall;
	PARTICLE *pStore;

	if (pLower > pUpper) return(-1);
	else {
		/*
		 ** Grab a cell from the cell storage!
		 */
		mdlassert(pkd->mdl,pkd->iFreeCell < pkd->nNodes);
		c = pkd->iFreeCell++;
		pkdn = &pkd->kdNodes[c];
		pkdn->pLower = pLower;
		pkdn->pUpper = pUpper;
		pStore = pkd->pStore;
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pStore[pLower].r[j];
			bnd.fMax[j] = pStore[pLower].r[j];
			bndBall.fMin[j] = pStore[pLower].r[j]-pStore[pLower].fBallMax;
			bndBall.fMax[j] = pStore[pLower].r[j]+pStore[pLower].fBallMax;
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pStore[i].r[j];
				else if (pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pStore[i].r[j];
				if (pStore[i].r[j]-pStore[i].fBallMax < bndBall.fMin[j]) 
					bndBall.fMin[j] = pStore[i].r[j]-pStore[i].fBallMax;
				if (pStore[i].r[j]+pStore[i].fBallMax > bndBall.fMax[j])
					bndBall.fMax[j] = pStore[i].r[j]+pStore[i].fBallMax;
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (bnd.fMax[j] > bnd.fMin[j])
			    bGoodBounds = 1;
			}

		pkdn->bnd = bnd;
		pkdn->bndBall = bndBall;
		
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (bnd.fMax[j]-bnd.fMin[j] > 
					bnd.fMax[d]-bnd.fMin[d]) {
					d = j;
					}
				}
			pkdn->iDim = d;
			/*
			 ** Now we do the split.
			 */
			pkdn->fSplit = 0.5*(bnd.fMin[d]+bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,pkdn->fSplit,pLower,pUpper);
			pkdn->iLower = BuildBinary(pkd,nBucket,pLower,m-1,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			pkdn->iUpper = BuildBinary(pkd,nBucket,m,pUpper,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			mdlassert(pkd->mdl,pkdn->iLower != -1 && pkdn->iUpper != -1);
			/*
			 ** Now calculate the mass, center of mass and mass weighted
			 ** softening radius.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			if (pkdn->iLower != -1) {
				fm = pkd->kdNodes[pkdn->iLower].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iLower].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iLower].r[j];
					}
				}
			if (pkdn->iUpper != -1) {
				fm = pkd->kdNodes[pkdn->iUpper].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iUpper].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iUpper].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		else {
			pkdn->iDim = -1; /* it is a bucket! */
			pkdn->fSplit = 0.0;
			pkdn->iLower = -1;
			pkdn->iUpper = -1;
			/*
			 ** Calculate the bucket quantities.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
				fm = pStore[i].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pStore[i].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pStore[i].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		/*
		 ** Calculate multipole moments.
		 */
		if(bGravity) {
		    pkdCalcCell(pkd,pkdn,pkdn->r,iOrder,&pkdn->mom);
		    dOpen = pkdCalcOpen(pkdn,iOpenType,dCrit,iOrder);
		    pkdn->fOpen2 = dOpen*dOpen;
		    }
		return(c);
		}
	}


void pkdThreadTree(PKD pkd,int iCell,int iNext)
{
	int l,u;

	if (iCell == -1) return;
	else if (pkd->kdNodes[iCell].iDim != -1) {
		l = pkd->kdNodes[iCell].iLower;
		u = pkd->kdNodes[iCell].iUpper;
		if (u == -1) {
			mdlassert(pkd->mdl,l != -1);
			pkdThreadTree(pkd,l,iNext);
			}
		else if (l == -1) {
			mdlassert(pkd->mdl,u != -1);
			pkdThreadTree(pkd,u,iNext);
			/* 
			 ** It is convenient to change the "down" pointer in this case.
			 */
			pkd->kdNodes[iCell].iLower = u;
			}
		else {
			pkdThreadTree(pkd,l,u);
			pkdThreadTree(pkd,u,iNext);
			}
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	else {
		pkd->kdNodes[iCell].iLower = -1;	/* Just make sure! */
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	}

/*
 * Builds a binary tree by splitting the largest spatial axis in half
 * The bounds for each subcell are recalculated (squeezed) to be exact
 *
 */
void pkdBuildBinary(PKD pkd,int nBucket,int iOpenType,double dCrit,
		    int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	int bEmpty = 0;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_TREEACTIVE );
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	/*
	 ** First problem is to figure out how many cells we need to
	 ** allocate. We need to do this in advance because we need to
	 ** synchronize to allocate the memory with mdlMalloc().
	 */
	if (bTreeActiveOnly) {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nTreeActive-1);
		}
	else {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nLocal-1);
		}
	/*
	 ** We need at least one particle per processor.
	mdlassert(pkd->mdl,pkd->nNodes > 0);
	 ** Let's relax this --trq
	 */
	if(pkd->nNodes == 0) {
	    pkd->nNodes = 1;
	    bEmpty = 1;
	    printf("id:%d has an empty local tree\n",pkd->idSelf);
	}

	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdNodes != NULL);
	/*
	 ** Now we really build the tree.
	 */
	pkd->iFreeCell = 0;
	if (bTreeActiveOnly) {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nTreeActive-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	else {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nLocal-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	mdlassert(pkd->mdl,bEmpty || pkd->iFreeCell == pkd->nNodes);
	if(bEmpty) {
	    KDN *pkdn;
	    int j;
	    
	    /*
	     * Set up an empty bucket.
	     */
	    pkd->iRoot = 0;
	    pkdn = &pkd->kdNodes[0];
	    pkdn->pLower = 0;
	    pkdn->pUpper = -1;
	    for (j=0;j<3;++j) {
	        pkdn->bnd.fMin[j] = FLOAT_MAXVAL;
		pkdn->bnd.fMax[j] = -FLOAT_MAXVAL;
		pkdn->bndBall.fMin[j] = FLOAT_MAXVAL;
		pkdn->bndBall.fMax[j] = -FLOAT_MAXVAL;
	        }
	    pkdn->iDim = -1; /* it is a bucket! */
	    pkdn->fSplit = 0.0;
	    pkdn->iLower = -1;
	    pkdn->iUpper = -1;
	    pkdn->fMass = 0.0;
	    pkdn->fSoft = 0.0;
	    for (j=0;j<3;++j) {
	        pkdn->r[j] = 0.0;
		}
	    if(bGravity) {
	        pkdCalcCell(pkd,pkdn,pkdn->r,iOrder,&pkdn->mom);
		pkdn->fOpen2 = 0.0;
	      }
        /*Determine whether the bounds are larger than half the periodic boundary --RB debug*/
	    }
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = pkd->kdNodes[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}


/*
 * Builds a binary tree by splitting on the largest spatial axis
 * so as to equal divide the particle by number
 * After the tree is built the bounds for each subcell are later recalculated (squeezed) to be exact
 *
 */
void pkdBuildLocal(PKD pkd,size_t nBucket,int iOpenType,double dCrit,
		   int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	char ach[256];
	BND bndDum;
	int bGoodBounds;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_TREEACTIVE );
	pkd->nBucket = nBucket;
	if (bTreeActiveOnly) n = pkd->nTreeActive;
	else n = pkd->nLocal;
	if (n == 0) {
		printf("id:%d has an empty local tree\n",pkd->idSelf);
		}
	pkd->nLevels = 1;
	l = 1;
	while ((size_t) n > nBucket) {
		n = n>>1;
		l = l<<1;
		++pkd->nLevels;
		}
	pkd->nSplit = l;
	pkd->nNodes = l<<1;
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdNodes != NULL);
	pkd->iFreeCell = pkd->nNodes;
	sprintf(ach,"nNodes:%d nSplit:%d nLevels:%d nBucket:%lu\n",
			pkd->nNodes,pkd->nSplit,pkd->nLevels,nBucket);
	mdlDiag(pkd->mdl,ach);
	/*
	 ** Set up ROOT node
	 */
	c = pkd->kdNodes;
	pkd->iRoot = ROOT;
	c[pkd->iRoot].pLower = 0;
	if (bTreeActiveOnly) {
		c[pkd->iRoot].pUpper = pkd->nTreeActive-1;
		}
	else {
		c[pkd->iRoot].pUpper = pkd->nLocal-1;
		}
	/*
	 ** determine the local bound of the particles.
	 */
	if (bTreeActiveOnly) {
		pkdCalcBound(pkd,&bndDum,&bndDum,&c[pkd->iRoot].bnd,&bndDum);
		}
	else {
		pkdCalcBound(pkd,&c[pkd->iRoot].bnd,&bndDum,&bndDum,&bndDum);
		}
	i = pkd->iRoot;
	while (1) {
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (c[i].bnd.fMax[j] > c[i].bnd.fMin[j])
				bGoodBounds = 1;
			}
		if (i < pkd->nSplit && (c[i].pUpper - c[i].pLower) > 0 && bGoodBounds) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;
			m = (c[i].pLower + c[i].pUpper)/2;
			pkdSelect(pkd,d,m,c[i].pLower,c[i].pUpper);
			c[i].fSplit = pkd->pStore[m].r[d];
			c[i].iLower = LOWER(i);
			c[i].iUpper = UPPER(i);
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			/* Unknown here; fix in pkdUpPass
			c[LOWER(i)].bndBall = c[i].bndBall; */
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			/* Unknown here; fix in pkdUpPass
			c[UPPER(i)].bndBall = c[i].bndBall; */
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			mdlassert(pkd->mdl,diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;		/* to indicate a bucket! */
			c[i].iLower = -1;
			c[i].iUpper = -1;
			SETNEXT(i);
			if (i == pkd->iRoot) break;
			}
		}
	pkdUpPass(pkd,pkd->iRoot,iOpenType,dCrit,iOrder, bGravity); 
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = c[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}


void pkdBucketWeight(PKD pkd,int iBucket,FLOAT fWeight)
{
	KDN *pbuc;
	int pj;
	
	pbuc = &pkd->kdNodes[iBucket];
	for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
		if (TYPEQueryACTIVE(&(pkd->pStore[pj])))
			pkd->pStore[pj].fWeight = fWeight;
		}
	}


void pkdColorCell(PKD pkd,int iCell,FLOAT fColor)
{
#ifdef COLORCODE
	KDN *pkdn;
	int pj;
	
	pkdn = &pkd->kdNodes[iCell];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		pkd->pStore[pj].fColor = fColor;
		}
#endif	
	}

#ifdef GR_DRAG
#define C2 (10065.3*10065.3)

void pkdGRDragGetSunAccel(PKD pkd,double dSunMass,double aSun[3])
{
	/*
	** This function resembles the GR_DRAG portion of
	** pkdSunIndirect(), but in this case computes the acceleration,
	** including GR terms, of the "Sun" itself due to all the other
	** particles in the system.
	*/

	PARTICLE *p;
	double r, r2,ir,r_hat[3];
	int n,i,k;

	for (k=0;k<3;k++)
		aSun[k] = 0.0;

	n = pkdLocal(pkd);

	for (i=0;i<n;i++) {
	  p = &pkd->pStore[i];
	  r2 = 0.0;
	  for (k=0;k<3;k++)
	    r2 += p->r[k]*p->r[k];
	  assert(r2 > 0.0);
	  r = sqrt(r2);
	  ir = 1.0 / (r - 3.0 * p->fMass / C2); /* We use this modified potential to replicate pericenter precession */
	  for (k=0;k<3;k++) {
	    r_hat[k] = p->r[k] / r;
	    aSun[k] -= p->fMass*ir*ir*r_hat[k]; /* This *should* be a + sign but aSun is added backwards in normal vanilla pkdgrav */
	  }
	}
}


#define GR_RADIUS_SQUARED (1.0*1.0) /* in AU (make this a parameter?) */

void pkdGRIntegrateCloseParticles(PKD pkd,double dSunMass,double dDelta,double dTime,int *nMerged, double *dMergerMassLost, double aSunCorrection[3])
{
	/*
	** Checks if a particle is too close to the "Sun" (central object).
	** If so, sets it inactive and calls external routine to compute
	** change in orbit.
	*/

	/* NOTE: ignores multistepping (active rungs---same likely true
	   for other GR functions...) */

	void grdrag(double px0, double py0, double pz0, double vx0, double
				vy0, double vz0, double m, double M, double interval, double *pxfin,
		    double *pyfin, double *pzfin, double *vxfin, double *vyfin, double *vzfin, int *merge, 
                    double *apocenter, double *pericenter, double subcode_r2, double *t_final, int returnr_p);

	PARTICLE *p;
	double r2,r2new,v2,x,y,z,vx,vy,vz,etot,lx,ly,lz,ltot2,semi,ecc, r_a,fVnew[3], r;
	double dFinalSubcodeTime,dDeltaV[3], dOldr[3], ir, r_hat[3];
	int i,n,bMerged,bWasOutsideSubcode, j;
        FILE *fp;

	*nMerged = 0;  
	*dMergerMassLost = 0;
	for (j=0;j<3;j++)
	  aSunCorrection[j] = 0.0; 

	n = pkdLocal(pkd);
	/* loop over all particles on this processor */
	for (i=0;i<n;i++) { 
	  p = &pkd->pStore[i];
	  mdlassert(pkd->mdl,TYPEQueryACTIVE(p));
	  r2 = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2]; /* This remains unchanged throughout this whole loop */
	  
	  if (r2 < C2 * dDelta * dDelta) {  /* We only need to worry about all of this if the particle is within dDelta*c of 0,0,0, which many aren't */
	   
	    /* We calculate where it will be to see if that will be inside our GR_RADIUS_SQUARED */
	    r2new = 0;
	    for (j=0;j<3;j++) {
	      fVnew[j] = p->v[j] + 0.5 * dDelta * p->a[j]; /* Velocity halfway through the step */
	      r2new += (p->r[j] + fVnew[j] * dDelta)*(p->r[j] + fVnew[j]*dDelta);
	    }
	    
	    bWasOutsideSubcode = p->bNoKickNoDrift ? 0 : 1; /*Becomes 1 if this the particle was integrated the last time using pkdgrav */
	    p->bNoKickNoDrift = 0;

	    /* Now we test to see if the particle will be inside GR_RADIUS; if it was before, we test to see if  it's still inside its particular entry radius */
	    
	    if (((bWasOutsideSubcode) && (r2new < GR_RADIUS_SQUARED)) ||
		((!bWasOutsideSubcode) && (r2 < p->dEntryR2)))
	      {
		/* This is the first entry so we save our radius */
		if (bWasOutsideSubcode) 
		    p->dEntryR2 = r2;
		
		p->bNoKickNoDrift = 1;
		grdrag(p->r[0],p->r[1],p->r[2],p->v[0],p->v[1],p->v[2],
		       p->fMass,dSunMass,dDelta,&x,&y,&z,&vx,&vy,&vz,&bMerged,&r_a,&p->dLastPericenter,p->dEntryR2,&dFinalSubcodeTime,bWasOutsideSubcode);
				
		if (bMerged) /* If there was a plunge, or if the particle didn't come out of the subcode, we merge.  */
		  {
		    ++(*nMerged);
		    fp = fopen("mergers.dat", "a");
		    if (r_a)
		      fprintf(fp, "Particle %d stayed in subcode_r at time %f with ap %f and peri %f\n", p->iOrder, dTime, r_a, p->dLastPericenter);
		    else
		      fprintf(fp, "Particle %d plunged into the SMBH at time %f\n", p->iOrder, dTime);
		    fclose(fp);
		    *dMergerMassLost += p->fMass;
		    pkdDeleteParticle(pkd, p);
		  }
		else 
		  {
		    /* Before we copy over the coordinates, we compute the average acceleration it experienced during the step.  We will
                       use this in the calculation of aSunCorrection[] */
		    
		    dDeltaV[0] = (vx - p->v[0]);
		    dDeltaV[1] = (vy - p->v[1]);
		    dDeltaV[2] = (vz - p->v[2]);
		    
		    if (!dFinalSubcodeTime) /* If we did a full subcode step */
		      {
			/* We first add the negative of what was added erroneously in GRDragGetSunAccel */
			r = sqrt(r2);
			ir = 1.0 / (r - 3.0 * p->fMass / C2);
			for (j=0;j<3;j++) 
			  r_hat[j] = p->r[j] / r;
			for (j=0;j<3;j++)
			  aSunCorrection[j] += p->fMass*ir*ir*r_hat[j]; /* We use a + sign here where a - sign was used in pkdGRDragGetSunAccel */
			/* Now we add the new contribution from our particle */
			for (j=0;j<3;j++)
			  aSunCorrection[j] += dDeltaV[j]*p->fMass/(dSunMass*dDelta);
		      }
		    
		    /* Before we copy, we may need to save the original coordinates if we have a time fragment */
		    if (dFinalSubcodeTime) {
		      for (j=0;j<3;j++)
			dOldr[j] = p->r[j];
		    }

		    /* We copy over the new coordinates */
		    p->r[0] = x;
		    p->r[1] = y;
		    p->r[2] = z;
		    p->v[0] = vx;
		    p->v[1] = vy;
		    p->v[2] = vz;

		    if (bWasOutsideSubcode)
		      {
			/* We print when we're just entering, and use the subcode to get the *actual* pericenter */
			r2new = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2]; /* We're reusing this variable so that r2 remains unchanged */
			v2 = p->v[0]*p->v[0] + p->v[1]*p->v[1] + p->v[2]*p->v[2];
			etot = 0.5*p->fMass*v2 - dSunMass*p->fMass/(sqrt(r2new)-3.0*dSunMass/C2);
			semi = -1.0*dSunMass*p->fMass/(2.0*etot);
			fp = fopen("module.dat", "a");
			lx = p->r[1]*p->v[2] - p->r[2]*p->v[1];
			ly = p->r[2]*p->v[0] - p->r[0]*p->v[2];
			lz = p->r[0]*p->v[1] - p->r[1]*p->v[0];
			ltot2 = lx*lx + ly*ly + lz*lz;
			fprintf(fp, "%8.5f %4d %5.2f %8.6f %8.5f ", dTime, p->iOrder, p->fMass, sqrt(r2new), sqrt(v2));
			fprintf(fp, "%10.3f %8.6f %6.5f %8.6f %8.5f\n", etot, semi, ecc = sqrt(1.0 - ltot2/((dSunMass+p->fMass)*semi)), sqrt(4*M_PI*M_PI*semi*semi*semi/(dSunMass+p->fMass)), p->dLastPericenter);
			fclose(fp);		  
			
		      }

		    if (dFinalSubcodeTime) /* If we're exiting the subcode, there will be a time fragment...*/
		      {
			p->dEntryR2 = 0.0;
			/* Now we decide which point we're at, based on whether or not we got more than 1/2 way through the step integrating */
			if (dFinalSubcodeTime < 0.5 * dDelta)
			  {
			    p->bNoKickNoDrift = 0;  /* We let PKDGRAV do its thing */
			    /* ...After we reset the particle's acceleration at its new location */
			    /* First we undo the old acceleration */
			    r = sqrt(r2); /* This is still the old r2 we calculated at the very beginning! */
			    ir = 1.0/(r - 3.0*dSunMass/C2);
			    for (j=0;j<3;j++) {
			      r_hat[j] = dOldr[j] / r;
			      p->a[j] += r_hat[j]*dSunMass*ir*ir; /* This subtracts out the old acceleration due to the SMBH */
			    }

			    /* this adds in the acceleration due to the "Sun" at the new position */
			    r = sqrt(p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2]);
			    ir = 1.0/(r-3.0*dSunMass/C2);
			    for (j=0;j<3;j++){
			      r_hat[j] = p->r[j]/r;
			      p->a[j] -= r_hat[j]*dSunMass*ir*ir; 
			    }
			  }
			else
			  {
			    /* We are pretending that we've done a full step and we need to add to aSunCorrection*/
			    
			    r = sqrt(r2); /* Still the old r2 from the beginning of the step*/
			    ir = 1.0 / (r - 3.0 * p->fMass / C2);
			    for (j=0;j<3;j++) 
			      r_hat[j] = p->r[j] / r;
			    for (j=0;j<3;j++)
			      aSunCorrection[j] += p->fMass*ir*ir*r_hat[j]; /* We use a + sign here where a - sign was used in pkdGRDragGetSunAccel */
			    /* Now we add the new contribution from our particle */
			    for (j=0;j<3;j++)
			       aSunCorrection[j] += dDeltaV[j]*p->fMass/(dSunMass*dFinalSubcodeTime);
			  }
		      } /* End if exiting subcode */
		  } /* Endelse didn't merge during subcode run */
	      } /* Endif uses subcode */
	  } /* Endif close enough to calculate stuff*/
	} /* Loop over particles */
}


#undef GR_RADIUS_SQUARED

void pkdGRCorrectaSun(PKD pkd, double aSunCorrection[3])
{
  /* When we are running the GR_DRAG module, the particles inside still give their contributions to aSun, which can be 
   * ridiculously large.  We correct for this in pkdGRIntegrateCloseParticles, but we need to pass it all back up to
   * the master level, then add it to each particle's acceration, as they have all already been erroneously calculated with 
   * this gigantic acceleration spike.  That is done in this simple function.
   */
  int n, i, j;
  PARTICLE *p;

  n = pkdLocal(pkd);
  for (i=0;i<n;i++) {
    p = &pkd->pStore[i];
    mdlassert(pkd->mdl,TYPEQueryACTIVE(p));
    for (j = 0 ; j < 3 ; ++j)
      p->a[j] -= aSunCorrection[j];
  }
}

#undef C2
#endif /*GR_DRAG*/

void pkdGravAll(PKD pkd,int nReps,int bPeriodic,int iOrder,int bEwald,int iEwOrder,
		   double fEwCut,double fEwhCut,int bDoSun,double dSunSoft,double *aSun,int *nActive,
		   double *pdPartSum,double *pdCellSum,double *pdSoftSum,CASTAT *pcs,
		   double *pdFlop)
{
	KDN *c = pkd->kdNodes;
	int iCell,n;
	FLOAT fWeight;
	double dFlopI, dFlopE;
	int i,j;
	BND bndActive;
	BND bndTmp;
	char achDiag[256];
	const int bDetailTimer = 0; /* Change to 1 to get detailed timing on
				     * interacts, walks, and Ewald
				     */

	*pdFlop = 0.0;
	pkdClearTimer(pkd,1);
	pkdClearTimer(pkd,2);
	pkdClearTimer(pkd,3);
	/*
	 ** Set up Ewald tables and stuff.
	 */
	if (bPeriodic && bEwald) {
	    pkdEwaldInit(pkd,fEwhCut,iEwOrder);	/* ignored in Flop count! */
		}
	/*
	 ** Start particle caching space (cell cache already active).
	 */
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
			   pkdLocal(pkd));
	/*
	 ** Walk over the local buckets!
	 */
	*nActive = 0;
	*pdPartSum = 0.0;
	*pdCellSum = 0.0;
	*pdSoftSum = 0.0;
	iCell = pkd->iRoot;
	if(!bDetailTimer) pkdStartTimer(pkd,2);

	while (iCell != -1) {
		if (c[iCell].iLower != -1) {
			iCell = c[iCell].iLower;
			continue;
			}
		n = 0;
		for(j = 0; j < 3; j++) {
		    bndActive.fMin[j] = FLOAT_MAXVAL;
		    bndActive.fMax[j] = -FLOAT_MAXVAL;
		    }
		for (i=c[iCell].pLower;i<=c[iCell].pUpper;++i) {
			if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			    ++n;
			    for (j=0;j<3;++j) {
					if (pkd->pStore[i].r[j] < bndActive.fMin[j]) 
						bndActive.fMin[j] = pkd->pStore[i].r[j];
					if (pkd->pStore[i].r[j] > bndActive.fMax[j])
						bndActive.fMax[j] = pkd->pStore[i].r[j];
					}
			    }
			}
		if (n > 0) {
		
			/* 
			 * set bounds to bounds of active particles.
			 */
			bndTmp = c[iCell].bnd;
			c[iCell].bnd = bndActive;
			/*
			 ** Calculate gravity on this bucket.
			 */
			if(bDetailTimer) pkdStartTimer(pkd,1);
			pkdBucketWalk(pkd,iCell,nReps,iOrder); /* ignored in Flop count! */
			if(bDetailTimer) pkdStopTimer(pkd,1);
			c[iCell].bnd = bndTmp;
			*nActive += n;
			*pdPartSum += n*pkd->nPart + 
				n*(2*(c[iCell].pUpper-c[iCell].pLower) - n + 1)/2;
			*pdCellSum += n*pkd->nCellNewt;
			*pdSoftSum += n*pkd->nCellSoft;
			if(bDetailTimer) pkdStartTimer(pkd,2);
			dFlopI = pkdBucketInteract(pkd,iCell,iOrder);

			*pdFlop += dFlopI;
			if(bDetailTimer) pkdStopTimer(pkd,2);
			/*
			 * Now do Ewald part.
			 */

			if (bPeriodic && bEwald) {
				if(bDetailTimer) pkdStartTimer(pkd,3);
				dFlopE = pkdBucketEwald(pkd,iCell,nReps,fEwCut,iEwOrder);
				*pdFlop += dFlopE;
				if(bDetailTimer) pkdStopTimer(pkd,3);
				}
			else {
			    dFlopE = 0.0;
			    }
			fWeight = dFlopI + dFlopE;
			/*
			 ** pkdBucketWeight, only updates the weights of the active
			 ** particles. Although this isn't really a requirement it
			 ** might be a good idea, if weights correspond to different 
			 ** tasks at different times.
			 */
			pkdBucketWeight(pkd,iCell,fWeight);
			}
		iCell = c[iCell].iUpper;
		}
	if(!bDetailTimer) pkdStopTimer(pkd,2);
	if (bDoSun) {
		const double dTinyBox = 1e-14;
		int iDummy = pkd->nStore;
		/*
		 ** Calculate the indirect interaction for solar system problems.
		 ** we need a "dummy" particle that is pointed to by the cell,
		 ** pkd->iFreeCell which also needs its bounds set correctly for 
		 ** this to work.
		 ** Don't allow periodic BCs at the moment.
		 */
		mdlassert(pkd->mdl,nReps == 0);
		mdlassert(pkd->mdl,bPeriodic == 0);
		for (j=0;j<3;++j) {
			pkd->pStore[iDummy].r[j] = 0.0;
			pkd->pStore[iDummy].a[j] = 0.0;
			c[pkd->iFreeCell].bnd.fMin[j] = -dTinyBox;
			c[pkd->iFreeCell].bnd.fMax[j] = dTinyBox;
			}
		TYPESet(&(pkd->pStore[iDummy]),TYPE_ACTIVE);
		pkd->pStore[iDummy].fPot = 0;
		pkd->pStore[iDummy].fSoft = dSunSoft;
		c[pkd->iFreeCell].fSoft = dSunSoft;
		c[pkd->iFreeCell].pLower = iDummy;
		c[pkd->iFreeCell].pUpper = iDummy;
		c[pkd->iFreeCell].iDim = -1;	/* It is really a bucket! */

		pkdBucketWalk(pkd,pkd->iFreeCell,nReps,iOrder);
		pkdBucketInteract(pkd,pkd->iFreeCell,iOrder);

		/*
		 ** Now we should have the indirect acceleration contained in the 
		 ** particle pointed to by iDummy. Now we have to bring this acceleration
		 ** back up to the master level. We can then take care of it with an
		 ** external potential call.
		 */
		for (j=0;j<3;++j) {
			aSun[j] = pkd->pStore[iDummy].a[j];
			}
		}
	/*
	 ** Get caching statistics.
	 */
	pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
	pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
	pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
	pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
	pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
	pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
	/*
	 ** Stop particle caching space.
	 */
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	sprintf(achDiag, "nMaxPart: %d, nMaxSoftCell: %d, nMaxNewtCell: %d\n",
		pkd->nMaxPart, pkd->nMaxCellSoft, pkd->nMaxCellNewt);
	mdlDiag(pkd->mdl, achDiag);

	}


void pkdSunIndirect(PKD pkd,double *aSun,int bDoSun,double dSunMass,double dSunSoft)
{
	PARTICLE *p;
	double t,idt2,a,b;
	int i,j,n;

#ifdef GR_DRAG
	const double c2 = 1.00653e4*1.00653e4; /* c = 1.00653040e4 in pkdgrav units */
	double r_hat[3];
#endif

#ifdef COLLISIONS
	pkd->bRepel = 0; /* initialization needed because of call to SPLINE() below */
#endif
#ifdef GR_DRAG
	assert(bDoSun);
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			t = 0;
			for (j=0;j<3;++j) t += p[i].r[j]*p[i].r[j];
			if(t == 0)
				a = b = 0;
			else
				SPLINE(t, (p[i].fSoft + dSunSoft), a, b);
			p[i].fPot -= dSunMass*a;
			idt2 = (p[i].fMass + dSunMass)*b;
			if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;

			/*
			 ** The way that aSun[j] gets added in here is confusing, but this
			 ** is the correct way! (aSun[] is the acceleration on the Sun).
			 */
			if (bDoSun) {
#ifdef GR_DRAG
				for (j=0;j<3;j++) 
					r_hat[j] = p[i].r[j]*a;
      				for (j=0;j<3;++j)
					p[i].a[j] -= aSun[j] + dSunMass*(1.0 / (1.0/a - 3.0 * dSunMass / c2))*(1.0 / (1.0/a - 3.0 * dSunMass / c2)) * r_hat[j];
	
#else /* GR_DRAG */
				b *= dSunMass;
				for (j=0;j<3;++j) {					
					p[i].a[j] -= (aSun[j] + p[i].r[j]*b);
					}				
#endif /* !GR_DRAG */
				}
			else {
				b *= p[i].fMass;
				for (j=0;j<3;++j) {
					p[i].a[j] -= (aSun[j] - p[i].r[j]*b);
					}
				}
			}
		}
	}


void pkdLogHalo(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	/* v in (kpc)/(4.691822922e16s) 128 km/s x 1.5233 */
	const double v = 194.9848486;
	const double d = 12.0;	/* in kpc */
	const double fl = 1.0;	/* flattening */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the logarithmic halo potential.
			 */
			double r = sqrt(x*x + y*y + z*z/(fl*fl));
			double C = 1.0/(r*r + d*d);	
			p[i].a[0] -= 2*v*v*x*C;
			p[i].a[1] -= 2*v*v*y*C;
			p[i].a[2] -= 2*v*v*z*C/(fl*fl);
			p[i].fPot += v*v*log(r*r + d*d);
			}
		}
	}


void pkdHernquistSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_s = 3.4e5;	/* in 10^5 M_sun */
	const double c = 0.7;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double A = 1.0/(r + c);	
			p[i].a[0] -= M_s*A*A*x/r;
			p[i].a[1] -= M_s*A*A*y/r;
			p[i].a[2] -= M_s*A*A*z/r;
			p[i].fPot -= M_s/(r + c);
			}
		}
	}


void pkdNFWSpheroid(PKD pkd, double M_200, double r_200, double c, double dSoft)
{
	PARTICLE *p;
	int i,n;
        /* M_200, r_200, c, and dSoft are now set in param file as
         * dNFWm200, dNFWr200, dNFWconc, and dNFWsoft */
	const double G = 1;
	/* Assuming G = 1 (this sets a timescale) */
	/* TimeUnit = sqrt(kpc^3/(G*1e12 Msun)) = 1.1285945e+09 yr */
        /* Vunit = 2073.8081 km/s */

	const double eps = c*dSoft/r_200; 

	/* r=r_200 (x=1, cx=c), M=M_200 */
	const double M_const = M_200/
	  ( (1./3.)*eps*eps/(1+eps)/(1+eps) 
	    - 1/(1+eps) + log((1+c)/(1+eps)) + 1/(1+c) );
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double A, fPot;
			double cx = r*(c/r_200);
			
			if (cx < eps) {
			  fPot = G*M_const*c/r_200* ((1./6.)*( cx*cx - eps*eps )/(eps*(1+eps)*(1+eps)) 
			    - (1./3.)*eps/(1+eps)/(1+eps) - 1/(1+eps));

			  A = G*M_const* 		
			    (1./3.)/(eps*(1+eps)*(1+eps))
			    *(c*c*c)/(r_200*r_200*r_200);
			}			  
			else {
			  fPot = G*M_const*c/r_200 * 
			    (( -(1./3.) *eps*eps/(1+eps)/(1+eps) - eps/(1+eps) 
				   - log((1+cx)/(1+eps)) ) /cx);

			  A = G*M_const* 
			    ( (1./3.)*eps*eps/(1+eps)/(1+eps) 
			       - 1/(1+eps) + log((1+cx)/(1+eps)) + 1/(1+cx) )
			    /(r*r*r);
			}

			/*fprintf(stderr,"%i: %f %f %f  %f %f %f %f\n",p[i].iOrder,x,y,z,r,fPot,-A*r,A); */
			p[i].a[0] -= A*x;
			p[i].a[1] -= A*y;
			p[i].a[2] -= A*z;
			p[i].fPot += fPot;
			}
		}
	}

/** JH Feb 4 2004 
elliptical potential using Hernquist model for luminous matter
and Hernquist or NFW model for dark matter */ 

void pkdElliptical(PKD pkd, int bEllipticalDarkNFW)
{
	PARTICLE *p;
	int i,n;
/* if bEllipticalDarkNFW=0 in param file uses hernquist dark matter otherwise uses NFW dark matter */
	const double M_l =1.45e11/6.45e9 ;	/* mass of luminous matter in 6.45e9 M_sun */
	const double a = 2.856/3.5;     /* hernquist scale length for light matter in 6 kpc */
	const double pi=  3.1415; 
	const double M_d= 7.2e11/6.45e9;     /* mass of dark matter in 6.45e9 M_sun */
        const double d= 31.7/3.5;     /* hernquist scale length for dark matter in 6 kpc */
        const double r_s=14.4/3.5 ;     /* NFW scale radius in 6 kpc  */
        const double rho_s= (0.010e9/6.45e9)*3.5*3.5*3.5 ;   /* NFW density at scale radius M_sun/kpc^3 converted to dMsolUnits*/
        const double eps1= 1e-5 ;   /* small softening for hernquist model */
        const double dsoft=1e-5 ;
	p = pkd->pStore;
	n = pkdLocal(pkd);
	/*	printf("bEllipticalDarkNFW= %d", bEllipticalDarkNFW); */
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
	                 */
			double r = sqrt(x*x + y*y + z*z);
			double A = 1.0/(r + a)  ;	/*light matter component */
			double Ah=  1.0/(r+ d) ;        /*hernquist dark matter component */
			double Anfw= 4*pi*rho_s*r_s*r_s*r_s ;     /* nfw dark matter component */
			if (!bEllipticalDarkNFW) {
			  p[i].a[0] -= M_l*A*A*x/(r + eps1) + M_d*Ah*Ah*x/(r+ eps1);
			  p[i].a[1] -= M_l*A*A*y/(r +eps1)+ M_d*Ah*Ah*y/(r + eps1);
			  p[i].a[2] -= M_l*A*A*z/(r + eps1)+ M_d*Ah*Ah*z/(r +eps1);
			  p[i].fPot -= M_l*A + M_d*Ah ;
			} 
			else {
			      if (r>=dsoft) {
                              p[i].a[0]-= M_l*A*A*x/(r+ eps1) +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*x/r ; 
			    p[i].a[1]-= M_l*A*A*y/(r+eps1)  +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*y/r ;   
			    p[i].a[2]-= M_l*A*A*z/(r+ eps1) +Anfw*(r_s/(r+r_s)/(r*r) +log(r+r_s)/(r*r) -1/(r*r) -log(r_s)/(r*r))*z/r ;   
			    p[i].fPot-= M_l*A -Anfw*(-log(r+r_s) + log(r_s))/r;
			    }
			      else {
			    p[i].a[0]-= M_l*A*A*x/(r+ eps1) +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*x/dsoft ; 
			    p[i].a[1]-= M_l*A*A*y/(r+eps1)  +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*y/dsoft ;   
			    p[i].a[2]-= M_l*A*A*z/(r+ eps1) +Anfw*(r_s/(dsoft+r_s)/(dsoft*dsoft) +log(dsoft+r_s)/(dsoft*dsoft) -1/(dsoft*dsoft) -log(r_s)/(dsoft*dsoft))*z/dsoft ;   
			    p[i].fPot-= M_l*A -Anfw*(-log(dsoft+r_s) + log(r_s))/dsoft; 
			    }
			  }
		}
	}
}


void pkdHomogSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_s = 5;	
	const double r_s = 10;
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double Mr,A;
			Mr = (r < r_s ? pow(r/r_s,3.)*M_s : M_s);
			A = Mr/(r*r*r);
			p[i].a[0] -= A*x;
			p[i].a[1] -= A*y;
			p[i].a[2] -= A*z;
			p[i].fPot += 2*( r < r_s ? 0.5*(Mr/r-3*M_s/r_s) : -M_s/r );
			}
		}
	}

void pkdBodyForce(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			if (p[i].r[2]>0) {
				p[i].a[2] -= 1;
				p[i].fPot += p[i].r[2];
				}
			else {
				p[i].a[2] += 1;
				p[i].fPot -= p[i].r[2];
				}
			}
		}
	}

void pkdMiyamotoDisk(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_d = 1.0e6;	/* in 10^5 M_sun */
	const double aa = 6.5;		/* in kpc */
	const double b = 0.26;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the Miyamoto Disk potential
			 */
			double U = sqrt(z*z + b*b);
			double W = aa + U;
			double A = sqrt(x*x + y*y + W*W);
			p[i].a[0] -= M_d*x/(A*A*A);
			p[i].a[1] -= M_d*y/(A*A*A);
			p[i].a[2] -= M_d*z*W/(A*A*A*U);
			p[i].fPot -= M_d/A;
			}
		}
	}

/*
 ** This is a copy of the Miyamoto function above, but change it as you
 ** like.
 */
void pkdTimeVarying(PKD pkd,double dTime)
{
	const double M_d0 = 1.0e6;	/* in 10^5 M_sun */
	const double aa = 6.5;		/* in kpc */
	const double b = 0.26;		/* in kpc */
	PARTICLE *p;
	int i,n;
	double M_d;
	
	p = pkd->pStore;
	n = pkdLocal(pkd);
	M_d = M_d0 * dTime;
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the Miyamoto Disk potential
			 */
			double U = sqrt(z*z + b*b);
			double W = aa + U;
			double A = sqrt(x*x + y*y + W*W);
			p[i].a[0] -= M_d*x/(A*A*A);
			p[i].a[1] -= M_d*y/(A*A*A);
			p[i].a[2] -= M_d*z*W/(A*A*A*U);
			p[i].fPot -= M_d/A;
			}
		}
	}

#ifdef ROT_FRAME

/*#define RF_SPINUP_ONLY*/

void
pkdRotFrame(PKD pkd,double dOmega,double dOmegaDot)
{
	/* WARNING: p[i].fPot not updated */

#ifdef RF_SPINUP_ONLY
	PARTICLE *p = pkd->pStore;
	int i;
	for (i=0;i<pkd->nLocal;i++) {
		p[i].a[0] -= dOmegaDot*p[i].r[1];
		p[i].a[1] += dOmegaDot*p[i].r[0];
		}
#else
	/* note Omega & dOmega/dt are assumed to be in the +z direction */

	PARTICLE *p;
	double w2,w2r[2],wxv[2],dwxr[2];
	int i,k;

	p = pkd->pStore;
	w2 = dOmega*dOmega;
	for (i=0;i<pkd->nLocal;i++) {
		if (!TYPEQueryACTIVE(&p[i])) continue;
		w2r[0] = -w2*p[i].r[0];
		w2r[1] = -w2*p[i].r[1];
		wxv[0] = -dOmega*p[i].vPred[1];
		wxv[1] =  dOmega*p[i].vPred[0];
		dwxr[0] = -dOmegaDot*p[i].r[1];
		dwxr[1] =  dOmegaDot*p[i].r[0];
		for (k=0;k<2;k++)
			p[i].a[k] -= (w2r[k] + 2*wxv[k] + dwxr[k]);
		}
#endif
	}

#endif /* ROT_FRAME */

void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double L[])
{
	/* L is calculated with respect to the origin (0,0,0) */

	PARTICLE *p;
	FLOAT rx,ry,rz,vx,vy,vz;
	int i,n;

#ifdef COLLISIONS
	FLOAT wx,wy,wz,moi;
#endif

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*T = 0.0;
	*U = 0.0;
	*Eth = 0.0;
	L[0] = L[1] = L[2] = 0;
	for (i=0;i<n;++i) {
		rx = p[i].r[0]; ry = p[i].r[1]; rz = p[i].r[2];
		vx = p[i].v[0]; vy = p[i].v[1]; vz = p[i].v[2];
		*T += 0.5*p[i].fMass*(vx*vx + vy*vy + vz*vz);
		*U += 0.5*p[i].fMass*p[i].fPot;
#ifdef GASOLINE
		if (pkdIsGas(pkd,&p[i]))
			*Eth += p[i].fMass*p[i].u;
#endif
		L[0] += p[i].fMass*(ry*vz - rz*vy);
		L[1] += p[i].fMass*(rz*vx - rx*vz);
		L[2] += p[i].fMass*(rx*vy - ry*vx);
#ifdef COLLISIONS
		wx = p[i].w[0]; wy = p[i].w[1]; wz = p[i].w[2];
		moi = 0.4*p[i].fMass*RADIUS(&p[i])*RADIUS(&p[i]); /* 2/5 MR^2: uniform sphere */
		*T += 0.5*moi*(wx*wx + wy*wy + wz*wz);
		L[0] += moi*wx;
		L[1] += moi*wy;
		L[2] += moi*wz;
#endif
		}
	}


void pkdCalcEandLExt(PKD pkd,double *dMass,double dSumMR[],double dSumMV[],
					 double *dPot)
{
	PARTICLE *p;
	FLOAT m,r2;
	int i,k,n;

	/* Currently this is for the heliocentric reference frame only */

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*dMass = *dPot = 0;
	for (k=0;k<3;k++) dSumMR[k] = dSumMV[k] = 0;
	for (i=0;i<n;i++) {
		*dMass += (m = p[i].fMass);
		r2 = 0;
		for (k=0;k<3;k++) {
			dSumMR[k] += m*p[i].r[k];
			dSumMV[k] += m*p[i].v[k];
			r2 += p[i].r[k]*p[i].r[k];
			}
		if (r2 > 0) *dPot += m/sqrt(r2);
		}
	}

void pkdMassInR(PKD pkd, double R, double *pdMass, FLOAT *com)
{
    PARTICLE *p;
    int i;
    int nLocal;
    double dMass = 0.0;
    FLOAT fCom[3] = {0.0,0.0,0.0};

    p = pkd->pStore;
    nLocal = pkdLocal(pkd);
    for(i = 0; i < nLocal; i++) {
	double r2 = 0.0;
	int k;
	
	for(k = 0; k < 3; k++) {
	    r2 += p[i].r[k]*p[i].r[k];
	    }
	
	if(r2 <= R*R) {
	    dMass += p[i].fMass;
	    for(k = 0; k < 3; k++)
		fCom[k] += p[i].fMass*p[i].r[k];
	    }
	}
    *pdMass = dMass;
    for(i = 0; i < 3; i++)
	com[i] = fCom[i];
    }

/* 
 * Use the f and g functions to advance an unperturbed orbit.
 */
void fg(MDL mdl,double mu,FLOAT *x,FLOAT *v,double dt) {
	double f,g,fd,gd;			/* Gauss's f, g, fdot and gdot */
	double r,vsq;
	double u;					/* r v cos(phi) */
	double a;					/* semi-major axis */
	double e;					/* eccentricity */
	double ec,es;				/* e cos(E), e sin(E) */
	double en;					/* mean motion */
	double nf;					/* orbital frequency */
	double dec;					/* delta E */
	double dm;					/* delta mean anomoly */
	double lo = -2*M_PI;
	double up = 2*M_PI;
	double w;					/* function to zero */
	double wp;					/* first derivative */
	double wpp;					/* second derivative */
	double wppp;				/* third derivative */
	double dx,s,c,sh;
	double next;
	double converge;			/* converge criterion */
	int iter,i,j;
	const double DOUBLE_EPS = 1.2e-16;
	const int MAX_ITER = 256;

	/* 
	 * Evaluate some orbital quantites.
	 */
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	u = x[0]*v[0] + x[1]*v[1] + x[2]*v[2];
	a = 1/(2/r-vsq/mu);
	en = sqrt(mu/(a*a*a));
	ec = 1-r/a;
	es = u/(en*a*a);
	e = sqrt(ec*ec + es*es);
#if (0)
	/*
	 ** Only for GR!
	 */
	dt *= 1 - 1.5*mu/(csq*a);
#endif
	nf = en/(2*M_PI);
	j = (int)(nf*dt);
	dt -= j/nf;					/* reduce to single orbit */
	dm = en*dt - es;
	if ((es*cos(dm)+ec*sin(dm)) > 0) dec = dm + 0.85*e;
	else dec = dm - 0.85*e;
	dm = en*dt;					/* reset mean anomoly */
	converge = fabs(dm*DOUBLE_EPS);
	/*
	 * First double precision iterations for dec solution.
	 * This solves Kepler's equation in difference form:
	 * dM = dE - e cos(E_0) sin(dE) + e sin(E_0)(1 - cos(dE))
	 * Note that (1 - cos(dE)) == 2*sin(dE/2)*sin(dE/2).  The latter is
	 * probably more accurate for small dE.
	 */
	for(iter=1;iter<=MAX_ITER;iter++) {
		s = sin(dec);
		c = cos(dec);
		sh = sin(0.5*dec);
		w = dec - ec*s + es*2*sh*sh - dm;
		if (w > 0) up = dec;
		else lo = dec;
		wp = 1 - ec*c + es*s;
		wpp = ec*s + es*c;
		wppp = ec*c - es*s;
		dx = -w/wp;
		dx = -w/(wp + 0.5*dx*wpp);
		dx = -w/(wp + 0.5*dx*wpp + dx*dx*wppp/6);
		next = dec + dx;
		if (fabs(dx) <= converge) break;
		if (next > lo && next < up) dec = next;
		else dec = 0.5*(lo + up);
		if (dec==lo || dec==up) break;
		}
	if (iter>MAX_ITER) {
		char ach[256];

		sprintf(ach,"dec soln failed, dconverge: %g dx: %g dec: %g\n",
				converge, dx, dec);
		mdlDiag(mdl,ach);
		mdlassert(mdl, 0);
		}
	/*
	 * Update the orbit.
	 * See Danby, eq. 6.8.11 and 6.8.12, for expressions for f, g, fdot,
	 * and gdot.
	 * JS: changed f and gd expressions to take advantage of the half angle 
	 * formula given above.
	 */
	f = 1 - (a/r)*2*sh*sh;
	g = dt + (s-dec)/en;
	fd = -(a/(r*wp))*en*s;
	gd = 1 - 2*sh*sh/wp;
	for (i=0;i<3;i++) {
		s = f*x[i]+g*v[i];
		v[i] = fd*x[i]+gd*v[i];
		x[i] = s;
		}
	}

void pkdDrift(PKD pkd,double dDelta,FLOAT fCenter[3],int bPeriodic,int bFandG,FLOAT fCentMass)
{
	PARTICLE *p;
	int i,j,n;
	int bInBox = 1;
#ifdef SLIDING_PATCH
	FLOAT fShear;
#endif

	mdlDiag(pkd->mdl, "Into pkddrift\n");
	n = pkdLocal(pkd);
	if (bFandG)
		for (i=0;i<n;++i) {
			p = &pkd->pStore[i];
			fg(pkd->mdl,fCentMass + p->fMass,p->r,p->v,dDelta);
			}
	else {
		p = pkd->pStore;

		for (i=0;i<n;++i,++p) {
#ifdef AGGS
			if (IS_AGG(p))
				continue; /* skip aggregate particles; handled in msrAggsAdvance() */
#endif
#ifdef WALLS
			if (PARTICLE_STUCK(p))
				continue;
#endif
#ifdef GR_DRAG
			if (p->bNoKickNoDrift)
				continue;
#endif
#ifdef SLIDING_PATCH
			fShear = 0.0;
			p->bAzWrap = 0; /* reset azimuthal wrap flag */
#endif
#ifdef COLLMOD_ZML
			if (p->ctimer > 0) {
			    --p->ctimer;
			    if (p->ctimer == 0)
			    	p->iColor = 3;
			}
#endif /* COLLMOD_ZML */
			for (j=0;j<3;++j) {
				p->r[j] += dDelta*p->v[j];
				if (bPeriodic) {
					while (p->r[j] >= fCenter[j] + 0.5*pkd->fPeriod[j]) {/*DEBUG! was "if"!*/
						p->r[j] -= pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (pkd->PP->bPatch) {
							if (j == 0) { /* radial wrap */
								fShear = 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth; /* (dWidth is same as fPeriod[0]) */
								p->r[1] += SHEAR(-1,pkd->dTime + dDelta,pkd->PP);
								}
							/* apply randomization only if there was no radial wrap (fShear=0) */
							if (pkd->PP->bRandAzWrap == 1 && j == 1 && fShear == 0.0) /* azimuthal wrap */
								p->bAzWrap = 1;
							}
#endif /* SLIDING_PATCH */
						}
					while (p->r[j] < fCenter[j] - 0.5*pkd->fPeriod[j]) { /*DEBUG! was "if"!*/
						p->r[j] += pkd->fPeriod[j];
#ifdef SLIDING_PATCH
						if (pkd->PP->bPatch) {
							if (j == 0) {
								fShear = - 1.5*pkd->PP->dOrbFreq*pkd->PP->dWidth;
								p->r[1] += SHEAR(1,pkd->dTime + dDelta,pkd->PP);
								}
							/* apply randomization only if there was no radial wrap (fShear=0) */
							if (pkd->PP->bRandAzWrap == 1 && j == 1 && fShear == 0.0)
								p->bAzWrap = 1;
							}
#endif /* SLIDING_PATCH */
						}
					bInBox = bInBox && (p->r[j] >= fCenter[j]-0.5*pkd->fPeriod[j]);
					bInBox = bInBox && (p->r[j] <  fCenter[j]+0.5*pkd->fPeriod[j]);
					}
				}
#ifdef DEM_TRACK_ORIENT
			demUpdateOrient(p,dDelta);
#endif /* DEM_TRACK_ORIENT */
#ifdef SLIDING_PATCH
			p->v[1] += fShear;
#ifdef NEED_VPRED
			p->vPred[1] += fShear;
#endif /* NEED_VPRED */
			p->dPy -= fShear/3.0; /* angular momentum is also changed */
#endif /* SLIDING_PATCH */
			mdlassert(pkd->mdl, bInBox);
			}
		}
	mdlDiag(pkd->mdl, "Out of pkddrift\n");
	}

void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
			 double dvPredFacTwo, double duDelta, double duPredDelta,
			 int iGasModel, double z, double duDotLimit)
{
	/*
	** For non-cosmology cannonical KDK scheme,
	** dvFacOne: 1 on opening kick, 1 on closing kick
	** dvFacTwo: 0.5*dDelta on opening kick, 0.5*dDelta on closing kick
	** dvPredFacOne: 1 on opening kick, 1 on closing kick
	** dvPredFacTwo: 0 on opening kick, 0.5*dDelta on closing kick
	*/
	PARTICLE *p;
	int i,j,n;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef AGGS
		if (IS_AGG(p)) /* skip aggregate particles; handled in msrAggsKick() */
			continue;
#endif
#ifdef WALLS
		if (PARTICLE_STUCK(p))
			continue;
#endif
#ifdef GR_DRAG
		if (p->bNoKickNoDrift)
			continue;
#endif
		if (TYPEQueryACTIVE(p)) {
#ifdef NEED_VPRED
#ifdef GASOLINE
			if (pkdIsGas(pkd, p)) {
				for (j=0;j<3;++j) {
					p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
				}
				if (iGasModel != GASMODEL_ISOTHERMAL) {
#ifndef NOCOOLING				
				  p->uPred = p->u + p->uDot*duPredDelta;
				  p->u = p->u + p->uDot*duDelta;
#else
				  p->uPred = p->u + p->PdV*duPredDelta;
				  p->u = p->u + p->PdV*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) 
				  if (p->uPred < 0) p->uPred = 0;
				  if (p->u < 0) p->u = 0;
#endif
				  }
				}
#else /* GASOLINE */
			for (j=0;j<3;++j) {
				p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
#ifdef DEM
				p->wPred[j] = p->w[j]*dvPredFacOne + p->wDot[j]*dvPredFacTwo;
#endif /* DEM */
				}
#endif /* !GASOLINE */
#endif /* NEED_VPRED */
			for (j=0;j<3;++j) {
				p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
#ifdef DEM
				p->w[j] = p->w[j]*dvFacOne + p->wDot[j]*dvFacTwo;
#endif /* DEM */
				}
			}
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done pkdkick\n");
	}

#ifdef SLIDING_PATCH
void pkdKickPatch(PKD pkd, double dvFacOne, double dvFacTwo,
#ifdef NEED_VPRED
				  double dvPredFacOne, double dvPredFacTwo,
#endif
				  double dOrbFreq, /* Orbital Frequency of Patch */
				  int bOpen)	/* Is this an opening Kick? */

{
	PARTICLE *p;
	int i,j,n;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef AGGS
		if (IS_AGG(p)) { /* skip aggregate particles; COMs are handled in msrAggsKick()... */
			/* 
			   10/17/09: an aggregated particle's Py is ignored (and incorrect) until 
			   liberation from the agg, at which time the free particle recalculates Py 
			   (see pkdAggsDoCollision(), pkdAggsDelete(), pkdAggsCheckStress()...). 
			*/
			continue;
			}		
#endif /* AGGS */
		if (TYPEQueryACTIVE(p)) {
		    if(!bOpen) {
				/* perform Cross Hamiltonian */
				p->v[0] += 2.0*dvFacTwo*dOrbFreq*p->dPy;
				p->v[1] = p->dPy - 2*dOrbFreq*p->r[0];
				}
#ifdef NEED_VPRED
			for (j=0;j<3;++j) {
				p->vPred[j] = p->v[j]*dvPredFacOne + p->a[j]*dvPredFacTwo;
#ifdef DEM
				p->wPred[j] = p->w[j]*dvPredFacOne + p->wDot[j]*dvPredFacTwo;
#endif
				}
#endif /* NEED_VPRED */
			for (j=0;j<3;++j) {
				p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
#ifdef DEM
				p->w[j] = p->w[j]*dvFacOne + p->wDot[j]*dvFacTwo;
#endif
				}
		    if(bOpen) {
				p->dPy = p->v[1] + 2.0*dOrbFreq*p->r[0];
				/* perform Cross Hamiltonian */
				p->v[0] += 2.0*dvFacTwo*dOrbFreq*p->dPy;
				p->v[1] = p->dPy - dOrbFreq*p->r[0]
				  - dOrbFreq*(p->r[0] + 2.0*dvFacTwo*p->v[0]);
				}
			}
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done pkdkickpatch\n");
	}
#endif /* SLIDING_PATCH */

void pkdReadCheck(PKD pkd,char *pszFileName,int iVersion,int iOffset,
				  int nStart,int nLocal)
{
	FILE *fp;
	CHKPART cp;
	long lStart;
	int i,j,iRtnVal;

	pkd->nLocal = nLocal;
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	mdlassert(pkd->mdl,fp != NULL);
	lStart = iOffset+nStart*sizeof(CHKPART);
	iRtnVal = fseek(fp,lStart,SEEK_SET);
	mdlassert(pkd->mdl,iRtnVal == 0);
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		PARTICLE *p = &pkd->pStore[i];
	    
		iRtnVal = fread(&cp,sizeof(CHKPART),1,fp);
		mdlassert(pkd->mdl,iRtnVal == 1);
		p->iOrder = cp.iOrder;
		mdlassert(pkd->mdl,p->iOrder >= 0); /* sanity check */
		p->fMass = cp.fMass;
		mdlassert(pkd->mdl,p->fMass > 0.0);
#ifdef CHANGESOFT
		p->fSoft0 = cp.fSoft;
#endif
		p->fSoft = cp.fSoft;
		for (j=0;j<3;++j) {
			p->r[j] = cp.r[j];
			p->v[j] = cp.v[j];
#ifdef NEED_VPRED
			/*DEBUG! If accelerations depend on vPred, restarts will
			  change the integration in the KDK scheme unless vPred is
			  stored explicitly in the checkpoints...*/
			p->vPred[j] = cp.v[j];
#endif
			}
#ifdef GASOLINE
		p->u = cp.u;
		p->uPred = cp.u;
		assert(p->u >= 0);
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
#endif
#ifdef STARFORM 
		p->fESNrate = 0.0;
		p->fTimeForm = cp.fTimeForm;
		p->fMassForm = cp.fMassForm;
		for (j=0;j<3;++j) {
			p->rForm[j] = cp.rForm[j];
			p->vForm[j] = cp.vForm[j];
			}
		p->fDensity = cp.fDenForm;
		p->iGasOrder = cp.iGasOrder;
                p->fTimeCoolIsOffUntil = cp.fTimeCoolIsOffUntil;
                p->fNSN = 0.0;
                p->fNSNtot = 0.0;
                p->fMFracOxygen = cp.fMFracOxygen;
                p->fMFracIron = cp.fMFracIron;
#endif
#ifdef SIMPLESF
		p->fMassStar = cp.fMassStar;
		p->fTimeForm = cp.fTimeForm;
		for (j=0;j<3;++j) {
			p->rForm[j] = cp.rForm[j];
			p->vForm[j] = cp.vForm[j];
			}
		p->fDensity = cp.fDenForm;
		p->iGasOrder = cp.iGasOrder;
#endif
		p->fMetals = cp.fMetals;
#ifndef NOCOOLING		
		p->CoolParticle = cp.CoolParticle;
#endif
#endif
#ifdef COLLISIONS
		p->iOrgIdx = cp.iOrgIdx;
		for (j=0;j<3;++j)
			p->w[j] = cp.w[j];
		p->iColor = cp.iColor;
#endif /* COLLISIONS */
#ifdef ORIGIN_HISTOGRAM
		for(j = 0; j < NUM_ORIGIN_BINS; ++j)
			p->origin_bins[j] = cp.origin_bins[j];
#endif /* ORIGIN_HISTOGRAM */ 
		TYPEClear(p);
		TYPESet(p,TYPE_ACTIVE);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
#ifdef GR_DRAG
		p->bNoKickNoDrift = cp.bNoKickNoDrift;
		p->dLastPericenter = cp.dLastPericenter;
		p->dEntryR2 = cp.dEntryR2;
#endif
#ifdef JOHNNY
		p->bInKepler = cp.bInKepler;
#endif /* JOHNNY */
	}
	fclose(fp);
}


void pkdWriteCheck(PKD pkd,char *pszFileName,int iOffset,int nStart)
{
	FILE *fp;
	CHKPART cp;
	long lStart;
	int i,j,nLocal;
	int nout;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	mdlassert(pkd->mdl,fp != NULL);
	lStart = iOffset+nStart*sizeof(CHKPART);
	fseek(fp,lStart,0);
	/* 
	 ** Write Stuff!
	 */
	nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;++i) {
		PARTICLE *p = &pkd->pStore[i];

		cp.iOrder = p->iOrder;
		cp.fMass = p->fMass;
#ifdef CHANGESOFT
		cp.fSoft = p->fSoft0;
#else
		cp.fSoft = p->fSoft;
#endif
		for (j=0;j<3;++j) {
			cp.r[j] = p->r[j];
			cp.v[j] = p->v[j];
			}
#ifdef GASOLINE
		cp.u = p->u;
		cp.fMetals = p->fMetals;
#ifndef NOCOOLING		
		cp.CoolParticle = p->CoolParticle;
#endif
#ifdef STARFORM
		cp.fTimeForm = p->fTimeForm;
		cp.fMassForm = p->fMassForm;
		for (j=0;j<3;++j) {
			cp.rForm[j] = p->rForm[j];
			cp.vForm[j] = p->vForm[j];
			}
		cp.fDenForm = p->fDensity;
		cp.iGasOrder = p->iGasOrder;
                cp.fTimeCoolIsOffUntil = p->fTimeCoolIsOffUntil;
                cp.fMFracOxygen = p->fMFracOxygen;
                cp.fMFracIron = p->fMFracIron;
#endif
#ifdef SIMPLESF
		cp.fMassStar = p->fMassStar;
		cp.fTimeForm = p->fTimeForm;
		for (j=0;j<3;++j) {
			cp.rForm[j] = p->rForm[j];
			cp.vForm[j] = p->vForm[j];
			}
		cp.fDenForm = p->fDensity;
		cp.iGasOrder = p->iGasOrder;
#endif
#endif
#ifdef COLLISIONS
		cp.iOrgIdx = p->iOrgIdx;
		for (j=0;j<3;++j)
			cp.w[j] = p->w[j];
		cp.iColor = p->iColor;
#endif /* COLLISIONS */
#ifdef ORIGIN_HISTOGRAM
		for(j = 0; j < NUM_ORIGIN_BINS; ++j)
			cp.origin_bins[j] = p->origin_bins[j];
#endif /* ORIGIN_HISTOGRAM */
#ifdef GR_DRAG
		cp.bNoKickNoDrift = p->bNoKickNoDrift;
		cp.dLastPericenter = p->dLastPericenter;
		cp.dEntryR2 = p->dEntryR2;
#endif
#ifdef JOHNNY
		cp.bInKepler = p->bInKepler;
#endif /* JOHNNY */
		nout = fwrite(&cp,sizeof(CHKPART),1,fp);
		mdlassert(pkd->mdl,nout == 1);
		}
	nout = fclose(fp);
	mdlassert(pkd->mdl,nout == 0);
}


void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn)
{
	int i;

	if (pkd->kdTop != NULL) free(pkd->kdTop);
	if (pkd->piLeaf != NULL) free(pkd->piLeaf);
	pkd->kdTop = malloc(nCell*sizeof(KDN));
	mdlassert(pkd->mdl,pkd->kdTop != NULL);
	pkd->piLeaf = malloc(pkd->nThreads*sizeof(int));
	mdlassert(pkd->mdl,pkd->piLeaf != NULL);
	for (i=1;i<nCell;++i) {
		if (pkdn[i].pUpper) {
			pkd->kdTop[i] = pkdn[i];
			if (pkdn[i].pLower >= 0) pkd->piLeaf[pkdn[i].pLower] = i;
			}
		}
	}


void pkdCalcRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	int pj;
	double m,dx,dy,dz;
	double d2;

	/*
	 ** Initialize moments.
	 */
	pcc->xxxx = 0.0;
	pcc->xyyy = 0.0;
	pcc->xxxy = 0.0;
	pcc->yyyy = 0.0;
	pcc->xxxz = 0.0;
	pcc->yyyz = 0.0;
	pcc->xxyy = 0.0;
	pcc->xxyz = 0.0;
	pcc->xyyz = 0.0;
	pcc->xxzz = 0.0;
	pcc->xyzz = 0.0;
	pcc->xzzz = 0.0;
	pcc->yyzz = 0.0;
	pcc->yzzz = 0.0;
	pcc->zzzz = 0.0;
	pcc->xxx = 0.0;
	pcc->xyy = 0.0;
	pcc->xxy = 0.0;
	pcc->yyy = 0.0;
	pcc->xxz = 0.0;
	pcc->yyz = 0.0;
	pcc->xyz = 0.0;
	pcc->xzz = 0.0;
	pcc->yzz = 0.0;
	pcc->zzz = 0.0;
	/*
	 ** Calculate moments.
	 */
	pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - pkdn->r[0];
		dy = pkd->pStore[pj].r[1] - pkdn->r[1];
		dz = pkd->pStore[pj].r[2] - pkdn->r[2];

#ifdef REDUCED_EWALD
		d2 = dx*dx + dy*dy + dz*dz;
#else
		d2 = 0.0;
#endif
		pcc->xxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
		pcc->xyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
		pcc->xxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
		pcc->yyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
		pcc->xxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
		pcc->xxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
		pcc->xxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
		pcc->xyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
		pcc->xxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
		pcc->xyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
		pcc->xzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
		pcc->yzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
		pcc->zzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
		/*
		 ** Calculate octopole moment...
		 */
		pcc->xxx += m*(dx*dx*dx - 0.6*d2*dx);
		pcc->xyy += m*(dx*dy*dy - 0.2*d2*dx);
		pcc->xxy += m*(dx*dx*dy - 0.2*d2*dy);
		pcc->yyy += m*(dy*dy*dy - 0.6*d2*dy);
		pcc->xxz += m*(dx*dx*dz - 0.2*d2*dz);
		pcc->yyz += m*(dy*dy*dz - 0.2*d2*dz);
		pcc->xyz += m*dx*dy*dz;
		pcc->xzz += m*(dx*dz*dz - 0.2*d2*dx);
		pcc->yzz += m*(dy*dz*dz - 0.2*d2*dy);
		pcc->zzz += m*(dz*dz*dz - 0.6*d2*dz);
		}
	}


void pkdDistribRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	double tr;

	pkd->ilcnRoot = *pcc;
	/*
	 ** Must set the quadrupole, mass and cm.
	 */
	pkdn = &pkd->kdTop[ROOT];
#ifdef REDUCED_EWALD
	tr = pkdn->mom.Qxx + pkdn->mom.Qyy + pkdn->mom.Qzz;
#else
	tr = 0.0;
#endif
	pkd->ilcnRoot.m = pkdn->fMass;
	pkd->ilcnRoot.x = pkdn->r[0];
	pkd->ilcnRoot.y = pkdn->r[1];
	pkd->ilcnRoot.z = pkdn->r[2];
	pkd->ilcnRoot.xx = pkdn->mom.Qxx - tr/3.0;
	pkd->ilcnRoot.xy = pkdn->mom.Qxy;
	pkd->ilcnRoot.xz = pkdn->mom.Qxz;
	pkd->ilcnRoot.yy = pkdn->mom.Qyy - tr/3.0;
	pkd->ilcnRoot.yz = pkdn->mom.Qyz;
	pkd->ilcnRoot.zz = pkdn->mom.Qzz - tr/3.0;
	}


double pkdMassCheck(PKD pkd) 
{
	double dMass=0.0;
	int i,iRej;
	int nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;++i) {
		dMass += pkd->pStore[i].fMass;
		}
	iRej = pkdFreeStore(pkd) - pkd->nRejects;
	for (i=0;i<pkd->nRejects;++i) {
		dMass += pkd->pStore[iRej+i].fMass;
		}
	return(dMass);
	}

double pkdGetMassChange(PKD pkd) 
{
	/*
	** Returns total mass of particles added or deleted this step.
	** Must be called BEFORE msrAddDelParticles().
	*/
	double dMassChange = 0.0;
	int nLocal = pkdLocal(pkd);
	int i;

	for (i=0;i<nLocal;i++) {
		if (pkd->pStore[i].iOrder == -1)
			dMassChange += pkd->pStore[i].fMass; /* mass gain */
		else if (pkd->pStore[i].iOrder < -1)
			dMassChange -= pkd->pStore[i].fMass; /* mass loss */
		}
	return dMassChange;
	}

void pkdMassMetalsEnergyCheck(PKD pkd, double *dTotMass, double *dTotMetals, 
                    double *dTotOx, double *dTotFe, double *dTotEnergy) 
{
	int i;
	*dTotMass=0.0;
	*dTotMetals=0.0;
	*dTotOx=0.0;
	*dTotFe=0.0;
	*dTotEnergy=0.0;

	for (i=0;i<pkdLocal(pkd);++i) {
		*dTotMass += pkd->pStore[i].fMass;
#ifdef GASOLINE 
                *dTotMetals += pkd->pStore[i].fMass*pkd->pStore[i].fMetals;
#ifdef STARFORM
                *dTotOx += pkd->pStore[i].fMass*pkd->pStore[i].fMFracOxygen;
                *dTotFe += pkd->pStore[i].fMass*pkd->pStore[i].fMFracIron;
                if ( TYPETest(&pkd->pStore[i], TYPE_GAS) ){
                    *dTotEnergy += pkd->pStore[i].fMass*pkd->pStore[i].fESNrate;
                    }
#endif
#endif
		}
	}


void
pkdSetRung(PKD pkd, int iRung)
{
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		pkd->pStore[i].iRung = iRung;
		}
    }

void
pkdBallMax(PKD pkd, int iRung, int bGreater, double dhFac)
{
    int i;
    
    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
			(pkd->pStore[i].iRung == iRung ||
			 (bGreater && pkd->pStore[i].iRung > iRung))) {
			pkd->pStore[i].fBallMax = sqrt(pkd->pStore[i].fBall2)*dhFac;
			}
		}
    return;
    }

int
pkdActiveRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int nActive;
    char out[128];
    int nLocal = pkdLocal(pkd);
    
    nActive = 0;
    for (i=0;i<nLocal;++i) {
		if(pkd->pStore[i].iRung == iRung ||
		   (bGreater && pkd->pStore[i].iRung > iRung)) {
			TYPESet(&(pkd->pStore[i]),TYPE_ACTIVE);
			++nActive;
			}
		else
			TYPEReset(&(pkd->pStore[i]),TYPE_ACTIVE);
		}
    sprintf(out, "nActive: %d\n", nActive);
    mdlDiag(pkd->mdl, out);
    pkd->nActive = nActive;
    return nActive;
    }

int
pkdCurrRung(PKD pkd, int iRung)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
		if(pkd->pStore[i].iRung == iRung) {
			iCurrent = 1;
			break;
			}
		}
    return iCurrent;
    }

void
pkdGravStep(PKD pkd,double dEta)
{
	double dt;
    int i;

    for (i=0;i<pkdLocal(pkd);i++) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			mdlassert(pkd->mdl,pkd->pStore[i].dtGrav > 0.0);
			dt = dEta/sqrt(pkd->pStore[i].dtGrav);
			if (dt < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dt;
			}
		}
    }

void
pkdAccelStep(PKD pkd,double dEta,double dVelFac,double dAccFac,int bDoGravity,
             int bEpsAcc,int bSqrtPhi,double dhMinOverSoft)
{
    int i;
    double vel;
    double acc;
    int j;
    double dT;

    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			vel = 0;
			acc = 0;
			for (j=0;j<3;j++) {
				vel += pkd->pStore[i].v[j]*pkd->pStore[i].v[j];
				acc += pkd->pStore[i].a[j]*pkd->pStore[i].a[j];
				}
			mdlassert(pkd->mdl,vel >= 0.0);
			vel = sqrt(vel)*dVelFac;
			acc = sqrt(acc)*dAccFac;
			dT = FLOAT_MAXVAL;
			if (bEpsAcc && acc>0) {
#ifdef GASOLINE
			     if (pkdIsGas(pkd, &(pkd->pStore[i])) && dhMinOverSoft < 1 && pkd->pStore[i].fBall2<4.0*pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) {
			        if (pkd->pStore[i].fBall2 > 4.0*dhMinOverSoft
				    *pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) 
				   dT = dEta*sqrt(sqrt(0.25*pkd->pStore[i].fBall2)/acc);
			        else 
				   dT = dEta*sqrt((dhMinOverSoft*pkd->pStore[i].fSoft)/acc);
			        }
		             else 
#endif
			        dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
			        }
			if (bSqrtPhi && acc>0) {
				/*
				 ** NOTE: The factor of 3.5 keeps this criterion in sync
				 ** with DensityStep. The nominal value of dEta for both
				 ** cases is then 0.02-0.03.
				 */
				double dtemp =
					dEta*3.5*sqrt(dAccFac*fabs(pkd->pStore[i].fPot))/acc;
				if (dtemp < dT)
					dT = dtemp;
				}
			if (dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }

void
pkdDensityStep(PKD pkd,double dEta,double dRhoFac)
{
    int i;
    double dT;
    
    for (i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			dT = dEta/sqrt(pkd->pStore[i].fDensity*dRhoFac);
			if (dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }

int
pkdDtToRung(PKD pkd,int iRung,double dDelta,int iMaxRung,
	    int bAll, /* 0 => symplectic case */
	    int *pnMaxRung,	/* number of particles on MaxRung */
	    int *piMaxRungIdeal) /* preferred max rung */
{
    int i;
    int iMaxRungOut;
    int iTempRung;
    int iSteps;
    int nMaxRung;
    int iMaxRungIdeal;
    
    iMaxRungOut = 0;
    iMaxRungIdeal = 0;
    nMaxRung = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
		if(pkd->pStore[i].iRung >= iRung) {
			mdlassert(pkd->mdl,TYPEQueryACTIVE(&(pkd->pStore[i])));
			if(bAll) {          /* Assign all rungs at iRung and above */
			        assert(pkd->pStore[i].fSoft > 0.0);
			        assert(pkd->pStore[i].dt > 0.0);
				iSteps = floor(dDelta/pkd->pStore[i].dt);
				/* insure that integer boundary goes
				   to the lower rung. */
				if(fmod(dDelta,pkd->pStore[i].dt) == 0.0)
				    iSteps--;
				iTempRung = iRung;
				if(iSteps < 0)
				    iSteps = 0;
				while(iSteps) {
					++iTempRung;
					iSteps >>= 1;
					}
				if(iTempRung >= iMaxRungIdeal)
					iMaxRungIdeal = iTempRung+1;
				if(iTempRung >= iMaxRung)
					iTempRung = iMaxRung-1;
				pkd->pStore[i].iRung = iTempRung;
				}
			else {
				if(dDelta <= pkd->pStore[i].dt) {
					pkd->pStore[i].iRung = iRung;
					}
				else {
					pkd->pStore[i].iRung = iRung+1;
					}
                }
			}
		if(pkd->pStore[i].iRung > iMaxRungOut) {
			iMaxRungOut = pkd->pStore[i].iRung;
		        nMaxRung = 1;
                        }
		else if (pkd->pStore[i].iRung == iMaxRungOut) 
		        nMaxRung ++;
		        
		}

    *pnMaxRung = nMaxRung;
    *piMaxRungIdeal = iMaxRungIdeal;
    return iMaxRungOut;
    }

void
pkdInitDt(PKD pkd,double dDelta)
{
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		if(TYPEQueryACTIVE(&(pkd->pStore[i])))
			pkd->pStore[i].dt = dDelta;
		}
    }
    
int pkdRungParticles(PKD pkd,int iRung)
{
	int i,n;

	n = 0;
	for (i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iRung == iRung) ++n;
		}
	return n;
	}

void 
pkdMoveParticle(PKD pkd, double *xcenter,double *xoffset,int iOrder)
{
    int i,j;
    
    for (i=0;i<pkdLocal(pkd);i++) {
	if (pkd->pStore[i].iOrder == iOrder) {
	    
#ifdef SLIDING_PATCH
	    pkd->pStore[i].v[1]+=1.5*pkd->PP->dOrbFreq*pkd->pStore[i].r[0];
#endif	
	    for (j=0;j<3;j++)
		pkd->pStore[i].r[j]=xcenter[j]+xoffset[j];
#ifdef SLIDING_PATCH
	    pkd->pStore[i].v[1]-=1.5*pkd->PP->dOrbFreq*pkd->pStore[i].r[0];
#endif
	    }
	}
    }

#ifdef STARFORM
void
pkdUnDeleteParticle(PKD pkd, PARTICLE *p)
{
    assert(TYPETest(p, TYPE_DELETED)); 

    TYPEReset(p, TYPE_DELETED); 
    p->iOrder = -2 - p->iOrder; 
    }
#endif /* STARFORM */

void
pkdDeleteParticle(PKD pkd, PARTICLE *p)
{
    p->iOrder = -2 - p->iOrder;
    TYPEClearACTIVE(p); 
    TYPESet(p, TYPE_DELETED); 
    }

void
pkdNewParticle(PKD pkd, PARTICLE p)
{
    mdlassert(pkd->mdl,pkd->nLocal < pkd->nStore);
    pkd->pStore[pkd->nLocal] = p;
    pkd->pStore[pkd->nLocal].iOrder = -1;
    pkd->nLocal++;
    }

void
pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
	     int *nDeltaStar)
{
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;
    
    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for(pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if(pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if(p->iOrder == -1) {
	    ++pj;
	    ++nNew;
#ifdef GASOLINE
	    ++ndStar;
#else
	    ++ndDark;
#endif
	    if(TYPEQueryACTIVE(p))
		++pkd->nActive;
	    continue;
	    }
	else if(p->iOrder < -1){
	    --newnLocal;
	    p->iOrder = -2 - p->iOrder;
	    if(pkdIsGas(pkd, p))
		--ndGas;
	    else if(pkdIsDark(pkd, p))
		--ndDark;
	    else if(pkdIsStar(pkd, p))
		--ndStar;
	    else
		mdlassert(pkd->mdl,0);
	    if(TYPEQueryACTIVE(p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->nLocal = newnLocal;
    }

void
pkdNewOrder(PKD pkd,int nStart)
{
    int pi;
    
    for(pi=0;pi<pkdLocal(pkd);pi++) {
		if(pkd->pStore[pi].iOrder == -1) {
			pkd->pStore[pi].iOrder = nStart++;
			}
		}
    }

void
pkdGetNParts(PKD pkd, struct outGetNParts *out )
{
    int pi;
    int n;
    int nGas;
    int nDark;
    int nStar;
	int iMaxOrderGas;
	int iMaxOrderDark;
	int iMaxOrderStar;
    PARTICLE *p;
    
    n = 0;
    nGas = 0;
    nDark = 0;
    nStar = 0;
	iMaxOrderGas = -1;
	iMaxOrderDark = -1;
	iMaxOrderStar = -1;
    for(pi = 0; pi < pkdLocal(pkd); pi++) {
		p = &pkd->pStore[pi];
		n++;
	    if(pkdIsGas(pkd, p)) {
			++nGas;
			if (p->iOrder > iMaxOrderGas) iMaxOrderGas = p->iOrder;
			}
	    else if(pkdIsDark(pkd, p)) {
			++nDark;
			if (p->iOrder > iMaxOrderDark) iMaxOrderDark = p->iOrder;
			}
	    else if(pkdIsStar(pkd, p)) {
			++nStar;
			if (p->iOrder > iMaxOrderStar) iMaxOrderStar = p->iOrder;
			}
		}

	out->n  = n;
    out->nGas = nGas;
    out->nDark = nDark;
    out->nStar = nStar;
	out->iMaxOrderGas = iMaxOrderGas;
	out->iMaxOrderDark = iMaxOrderDark;
	out->iMaxOrderStar = iMaxOrderStar;
    }

void
pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar,int nMaxOrder, int nMaxOrderGas,
	     int nMaxOrderDark)
{
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
/*  pkd->nMaxOrder = nMaxOrder; */
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;
}

void pkdCoolVelocity(PKD pkd,int nSuperCool,double dCoolFac,
					 double dCoolDens,double dCoolMaxDens)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (p[i].iOrder < nSuperCool && p[i].fDensity >= dCoolDens &&
			p[i].fDensity < dCoolMaxDens) {
			for (j=0;j<3;++j) {
#ifdef SUPERCOOL
				p[i].v[j] = p[i].vMean[j] + (p[i].v[j]-p[i].vMean[j])*dCoolFac;
#else
				p[i].v[j] *= dCoolFac;
#endif
				}
			}
		}
	}

int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) { 
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPEFilter(p,iFilterMask,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		else {
			TYPEReset(p,iSetMask);
			}
		}
    if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
    if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
#ifdef COOLDEBUG
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
		assert(p->u >= 0.0);
		assert(p->uPred >= 0.0);
#endif
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPETest(p,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		}
    /*
	   Need to fix this up:
	   if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
	   if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
	   if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
	   */
    return nActive;
    }

int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if (TYPETest(p,iTestMask)) {
			TYPEReset(p,iSetMask);
			nActive++;
			}
		}
    /*
	   if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
	   if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
	   if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
	   */
    return nActive;
    }

int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		if (TYPEFilter(p,iFilterMask,iTestMask)) {
			nActive++;
			}
		}
    {
	char debug[100];
	sprintf(debug, "Filter %d:%d, Counted: %d\n",iFilterMask,iTestMask,nActive);
	mdlDiag(pkd->mdl,debug);
	}
    return nActive;
    }

int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		/* DEBUG: Paranoia check */
		/* mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL)); */
		if (TYPETest(p,iTestMask)) {
			TYPESet(p,iSetMask);
			nActive++;
			}
		else {
			TYPEReset(p,iSetMask);
			}
		}
    if (iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveMaskRung(PKD pkd, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    
    nActive = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
		if(p->iRung == iRung || (bGreater && p->iRung > iRung)) {
			TYPESet(p,iSetMask);
			++nActive;
			}
		else
			TYPEReset( p, iSetMask );
		}
    sprintf(out,"nActive: %d\n",nActive);
    mdlDiag(pkd->mdl,out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveTypeRung(PKD pkd, unsigned iTestMask, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    int nLocal = pkdLocal(pkd);
    
    nActive = 0;
    for(i=0;i<nLocal;++i) {
        p = &pkd->pStore[i];
        /* DEBUG: Paranoia check */
        mdlassert(pkd->mdl,TYPETest(p,TYPE_ALL));
		if(TYPETest(p,iTestMask) && 
           (p->iRung == iRung || (bGreater && p->iRung > iRung))) {
			TYPESet(p,iSetMask);
			++nActive;
			}
		else
			TYPEReset( p, iSetMask );
		}
    sprintf(out,"nActive: %d\n",nActive);
    mdlDiag(pkd->mdl,out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

struct SortStruct {
  int iOrder;
  int iStore;
};

int CompSortStruct(const void * a, const void * b) {
  return ( ( ((const struct SortStruct *) a)->iOrder < ((const struct SortStruct *) b)->iOrder ? -1 : 1 ) );
}

void pkdSetTypeFromFileSweep(PKD pkd, int iSetMask, char *file, 
	   struct SortStruct *ss, int nss, int *pniOrder, int *pnSet) {
  int niOrder = 0, nSet = 0;
  int iOrder, iOrderOld, nRet;
  FILE *fp;
  int iss;

  fp = fopen( file, "r" );
  assert( fp != NULL );

  iss = 0;
  iOrderOld = -1;
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	while (ss[iss].iOrder < iOrder) {
	  iss++;
	  if (iss >= nss) goto DoneSS;
	}
	if (iOrder == ss[iss].iOrder) {
	  TYPESet(&(pkd->pStore[ss[iss].iStore]),iSetMask);
	  nSet++;
	}
  }
 
DoneSS:
  /* Finish reading file to verify consistency across processors */
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	}
  fclose(fp);

  *pniOrder += niOrder;
  *pnSet += nSet;

  return;
}


int pkdSetTypeFromFile(PKD pkd, int iSetMask, int biGasOrder, char *file, int *pniOrder, int *pnSet, int *pnSetiGasOrder)
{
  struct SortStruct *ss;
  int i,nss;

  *pniOrder = 0;
  *pnSet = 0;
  *pnSetiGasOrder = 0;

  nss = pkdLocal(pkd);
  ss = malloc(sizeof(*ss)*nss);
  assert( ss != NULL );

  for(i=0;i<pkdLocal(pkd);++i) {
	ss[i].iOrder = 	pkd->pStore[i].iOrder;
	ss[i].iStore = i;
  }

  qsort( ss, nss, sizeof(*ss), CompSortStruct );

  pkdSetTypeFromFileSweep(pkd, iSetMask, file, ss, nss, pniOrder, pnSet);

#if defined(STARFORM) || defined(SIMPLESF)
  if (biGasOrder) {
	for(i=0;i<pkdLocal(pkd);++i) {
	  ss[i].iOrder =   ( TYPETest(&pkd->pStore[i], TYPE_STAR) ? 
						 pkd->pStore[i].iGasOrder : -1 );
	  ss[i].iStore = i;
	}
	
	qsort( ss, nss, sizeof(*ss), CompSortStruct );
	
	pkdSetTypeFromFileSweep(pkd, iSetMask, file, ss, nss, pniOrder, pnSetiGasOrder);
  }
#endif

  free( ss );

  return *pnSet+*pnSetiGasOrder;
}


void pkdSetParticleTypes(PKD pkd, int nSuperCool)
{
    PARTICLE *p;
    int i, iSetMask;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		iSetMask = 0;
		if (pkdIsGasByOrder(pkd,p)) iSetMask |= TYPE_GAS;
		if (pkdIsDarkByOrder(pkd,p)) iSetMask |= TYPE_DARK;
		if (pkdIsStarByOrder(pkd,p)) iSetMask |= TYPE_STAR;
		if (p->iOrder < nSuperCool) iSetMask |= TYPE_SUPERCOOL;

		TYPESet(p,iSetMask);
		}
    }

int
pkdSoughtParticleList(PKD pkd, int iTypeSought, int nMax, int *n, struct SoughtParticle *sp)
{
    int i, nFound;
	PARTICLE *p = pkd->pStore;

	nFound = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
 	  if (TYPETest((p+i),iTypeSought)) {
		/* Don't halt on overflow -- just stop copying */
		if (nFound < nMax) {
		  sp[nFound].iOrder = p[i].iOrder;
		  sp[nFound].iActive = p[i].iActive;
		  sp[nFound].x = p[i].r[0];
		  sp[nFound].y = p[i].r[1];
		  sp[nFound].z = p[i].r[2];
		}
		nFound ++;
	  }
	}
	*n = nFound;
	for (i=0;i<nFound;i++) {
	  printf("pkd sub star %d: %g %g %g\n",i,sp[i].x,sp[i].y,sp[i].z);
	}

	return (nFound);
}

void
pkdCoolUsingParticleList(PKD pkd, int nList, struct SoughtParticle *l)
{
#ifndef NOCOOLING
    int i,j;
	double r2,r2min,dx;
	PARTICLE *p = pkd->pStore;

	assert(nList > 0);
    for(i=0;i<pkdLocal(pkd);++i) {
 	  if (TYPEQueryACTIVE(&p[i]) && TYPETest(&p[i],TYPE_GAS)) {
		dx = p[i].r[0]-l[0].x;
		r2 = dx*dx;
		dx = p[i].r[1]-l[0].y;
		r2 += dx*dx;
		dx = p[i].r[2]-l[0].z;
		r2 += dx*dx;
		r2min = r2;
		for (j=1;j<nList;j++) {
		  dx = p[i].r[0]-l[j].x;
		  r2 = dx*dx;
		  dx = p[i].r[1]-l[j].y;
		  r2 += dx*dx;
		  dx = p[i].r[2]-l[j].z;
		  r2 += dx*dx;
		  if (r2 < r2min) r2min = r2;
		}
#ifdef COOLING_DISK
		p[i].CoolParticle.r = sqrt(r2min);
#endif
	  }
	}
#endif
}

void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM)
{
    int i;

    for(i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iOrder < nGrowMass) {
			pkd->pStore[i].fMass += dDeltaM;
			}
		}
    }

void pkdInitAccel(PKD pkd)
{
    int i,j;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
  		        pkd->pStore[i].fPot = 0;
			for (j=0;j<3;++j) {
				pkd->pStore[i].a[j] = 0;
#if defined(GASOLINE) && defined(SHOCKTRACK)
				pkd->pStore[i].aPres[j] = 0;
#endif
#ifdef DEM
				pkd->pStore[i].wDot[j] = 0;
#endif /* DEM */
				pkd->pStore[i].dtGrav = 0; /* why within this for-loop? */
				}
			}
		}
    }

int pkdIsGasByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder <= pkd->nMaxOrderGas) return 1;
	else return 0;
	}

int pkdIsDarkByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderGas && p->iOrder <= pkd->nMaxOrderDark)
	    return 1;
	else return 0;
	}

int pkdIsStarByOrder(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderDark) return 1;
	else return 0;
	}

#ifdef GASOLINE

void pkdUpdateuDot(PKD pkd, double duDelta, double dTime, double z, int iGasModel, int bUpdateY )
{
#ifndef NOCOOLING	
	PARTICLE *p;
	int i,n;
	int bCool = 0;
	COOL *cl = NULL;
	COOLPARTICLE cp;
	double E,dt = 0;
#endif

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

#ifndef NOCOOLING
	switch (iGasModel) {
	case GASMODEL_COOLING:
		bCool = 1;
		cl = pkd->Cool;
		CoolSetTime( cl, dTime, z  );
		dt = CoolCodeTimeToSeconds( cl, duDelta );
		break;
		}

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) {
			if ( bCool 
#ifdef STARFORM
                            &&  ((dTime >= p->fTimeCoolIsOffUntil) 
                               || ((p->fESNrate + p->PdV)*duDelta + p->u < 0 ))  
#endif                                
			  ) {
				cp = p->CoolParticle;
				E = p->u;
#ifdef STARFORM
				cl->p = p;
                                CoolIntegrateEnergyEPDRCode(cl, &cp, &E, p->fESNrate + p->PdV, p->fDensity, p->r, dt);
#else
				CoolIntegrateEnergyEPDRCode(cl, &cp, &E, p->PdV, p->fDensity, p->r, dt);
#endif
				mdlassert(pkd->mdl,E > 0);

				p->uDot = (E - p->u)/duDelta;
#ifdef STARFORM
                                if ((p->uDot < (p->fESNrate + p->PdV)) 
                                    && (dTime < p->fTimeCoolIsOffUntil)) {
                                        p->uDot = p->fESNrate + p->PdV;
                                        }
#endif
#ifdef SIMPLESF 
				if (dTime < p->fTimeForm) {
					if (p->uDot<p->PdV) p->uDot = p->PdV;
					}
#endif

				if (bUpdateY) p->CoolParticle = cp;
				}
			else { 
#ifdef STARFORM
				p->uDot = p->PdV + p->fESNrate;
#else
				p->uDot = p->PdV;
#endif
				}
			}
		}
#endif
	pkdStopTimer(pkd,1);
	}

void pkdUpdateShockTracker(PKD pkd, double dDelta, double dShockTrackerA, double dShockTrackerB )
{
#ifdef SHOCKTRACK
	PARTICLE *p;
	int i,n;
	double conh,factor;
	double dv2,a2,ap2,adotap=1;

        if (dShockTrackerB == 0) printf("Doing cheap shock tracking\n");
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) {
		        p->ShockTracker = 0;

			a2 = ((p->a[0]*p->a[0])+(p->a[1]*p->a[1])+(p->a[2]*p->a[2]));
			ap2 = ((p->aPres[0]*p->aPres[0])+(p->aPres[1]*p->aPres[1])+(p->aPres[2]*p->aPres[2]));
			adotap = ((p->a[0]*p->aPres[0])+(p->a[1]*p->aPres[1])+(p->a[2]*p->aPres[2]));

#if 0
			if (!(i%2000)) {
			  /*
printf("aP %d: %g %g %g %g %g %g\n",i,p->a[0],p->a[1],p->a[2],p->aPres[0],p->aPres[1],p->aPres[2]);
			  */
printf("r %g PdV %g a %g %g %g SW %g %g %g\n",sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2]),
       p->PdV,((p->a[0]*p->aPres[0])+(p->a[1]*p->aPres[1])+(p->a[2]*p->aPres[2])),sqrt(a2),sqrt(ap2),
       sqrt(a2*0.25*p->fBall2)/(p->c*p->c),
       (p->a[0]*p->gradrho[0]+p->a[1]*p->gradrho[1]+p->a[2]*p->gradrho[2])/
       (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
	0.4/p->fBall2)/(p->c*p->c),p->BalsaraSwitch );
			}
#endif

			/* Rarefaction or Gravity dominated compression */
			if ( p->PdV < 0 || adotap < 0.5*a2) p->ShockTracker = 0;
			else {
			 if (dShockTrackerB == 0) {
			  if (a2 < ap2) 
			    dv2 = sqrt(a2*0.25*p->fBall2);
			  else
			    dv2 = sqrt(ap2*0.25*p->fBall2);
			 }
			 else {
			  if (a2 < ap2) 
			    dv2 = -p->fDensity*(p->a[0]*p->gradrho[0]+p->a[1]*p->gradrho[1]+p->a[2]*p->gradrho[2])/
			      (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
			       p->fDensity*p->fDensity*dShockTrackerB/p->fBall2);
			  else
			    dv2 = -p->fDensity*(p->aPres[0]*p->gradrho[0]+p->aPres[1]*p->gradrho[1]+p->aPres[2]*p->gradrho[2])/
			      (p->gradrho[0]*p->gradrho[0]+p->gradrho[1]*p->gradrho[1]+p->gradrho[2]*p->gradrho[2]+
			       p->fDensity*p->fDensity*dShockTrackerB/p->fBall2);

			  }
			 /*
			  if (dv2 < dShockTrackerA*p->c*p->c) p->ShockTracker = 0;
			 */
			 p->ShockTracker = dv2/(dShockTrackerA*p->c*p->c);
			 if (p->ShockTracker > 1) p->ShockTracker = 1;
			 }
			}
		}

#endif
        }

/* Note: Uses uPred */
void pkdAdiabaticGasPressure(PKD pkd, double gammam1, double gamma)
{
    PARTICLE *p;
    double PoverRho;
    int i;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
			PoverRho = gammam1*p->uPred;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(gamma*PoverRho);
			}
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) {
			printf("Pressure %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,
			       p->PoverRho2*p->fDensity*p->fDensity,p->c);
			}
#endif            
		}
    }

#ifndef NOCOOLING
/* Note: Uses uPred */
void pkdCoolingGasPressure(PKD pkd, double gammam1, double gamma)
{
    PARTICLE *p;
	COOL *cl = pkd->Cool;
    double PoverRho;
    int i;

	/* Coolng Setup function required? 
	   Abundances may be lagging by half a step ... */

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
		    CoolCodePressureOnDensitySoundSpeed( cl, &p->CoolParticle, p->uPred, p->fDensity, gamma, gammam1, &PoverRho, &(p->c) );
			p->PoverRho2 = PoverRho/p->fDensity;
			}
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) {
			printf("Pressure %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,
			       p->PoverRho2*p->fDensity*p->fDensity,p->c);
			}
#endif            
		}
    }
#endif

void pkdLowerSoundSpeed(PKD pkd, double dhMinOverSoft)
{
    PARTICLE *p;
    double dfBall2MinOverSoft,dfBall2Min,ratio;
    int i;

    dfBall2MinOverSoft = 4*dhMinOverSoft*dhMinOverSoft;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (pkdIsGas(pkd,p)) {
		        dfBall2Min = dfBall2MinOverSoft*p->fSoft*p->fSoft;
			if (p->fBall2 < dfBall2Min) {
			  ratio = p->fBall2/dfBall2Min;
			  p->PoverRho2 *= ratio;
			  p->c *= sqrt(ratio);
#ifdef DEBUG
		if ((p->iOrder % 100)==0) {
			printf("Pressure %i: %g %g, %g %g, %g %g\n",
			       p->iOrder,
			       sqrt(p->fBall2),sqrt(dfBall2Min),
			       p->PoverRho2,p->PoverRho2/ratio,
			       p->c,p->c/sqrt(ratio));
			}
#endif            
			  }
			}
		}
    }

void pkdInitEnergy(PKD pkd, double dTuFac, double z, double dTime )
{
    PARTICLE *p;
    int i;
#ifndef NOCOOLING
    COOL *cl;
	double T,E;

    cl = pkd->Cool;
	CoolSetTime( cl, dTime, z  );
#endif

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
		if (TYPEQueryTREEACTIVE(p) && pkdIsGas(pkd,p)) {
#ifndef NOCOOLING
			T = p->u / dTuFac;
			CoolInitEnergyAndParticleData( cl, &p->CoolParticle, &E, p->fDensity, T );
			p->u = E;
#endif
			p->uPred = p->u;
#ifdef DEBUG
			if ((p->iOrder % 1000)==0) {
				printf("InitEnergy %i: %f %g   %f %f %f %g\n",
					   p->iOrder,T,p->u * cl->dErgPerGmUnit,
					   p->CoolParticle.HI,p->CoolParticle.HeI,p->CoolParticle.HeII,p->fDensity*cl->dComovingGmPerCcUnit);
				} 
#endif            
		        }
                }
    }

#ifdef GLASS
/* Currently wired to have no more than two regions with different
   Pressures (densities) split by x=0 with a linear connection */
void pkdGlassGasPressure(PKD pkd, void *vin)
{
    PARTICLE *p;
    double PoverRho,xx,nsp=2.5;
    int i;
    struct inGetGasPressure *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                if (TYPEQueryTREEACTIVE(p)) {
    		        if (p->r[0] < -nsp*in->dGlassxL) {
			  if (p->r[0] > in->dxBoundL + nsp*in->dGlassxL)
			       PoverRho=in->dGlassPoverRhoL;  
			  else {
			       xx = ( p->r[0] - in->dxBoundL + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
              		       }
			  }
                        else if (p->r[0] > nsp*in->dGlassxR) {
			  if (p->r[0] < in->dxBoundR - nsp*in->dGlassxR)
			       PoverRho=in->dGlassPoverRhoR;  
			  else {
			       xx = ( p->r[0] - in->dxBoundR + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
			       }
			  }
			else {
			       xx = ( p->r[0] + nsp*in->dGlassxL )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoL + xx
				 *(in->dGlassPoverRhoR - in->dGlassPoverRhoL);
			       }

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->gamma*PoverRho);
		        }
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
		        printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
                }
    }

#endif

void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta)
{
	PARTICLE *p;
	int i,n;

#ifdef DEBUG
	printf("pkdKickRhopred: %g %g\n",dHubbFac,dDelta);
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef DEBUG
		if (pkdIsGas(pkd,(p+i)) && ((p+i)->iOrder % 3000)==0) {
			printf("Rhopreding %i: %i %i %f %f %f %f   %f %f %f %f\n",
			       (p+i)->iOrder,TYPEQueryACTIVE(p+i),
				   TYPEQueryTREEACTIVE(p+i),
			       sqrt(0.25*(p+i)->fBall2),(p+i)->fDensity,(p+i)->u,(p+i)->uPred,
			       (p+i)->a[0],(p+i)->a[1],(p+i)->a[2],(p+i)->PdV);
			}
#endif
		if(TYPEFilter( &p[i], TYPE_GAS|TYPE_ACTIVE, TYPE_GAS )) {
			p[i].fDensity = p[i].fDensity*(1 + dDelta*(dHubbFac - p[i].divv));
			}
		}
	}

int pkdSphCurrRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(pkdIsGas(pkd, &pkd->pStore[i]) &&
		   (pkd->pStore[i].iRung == iRung ||
			(pkd->pStore[i].iRung > iRung && bGreater))) {
			iCurrent = 1;
			break;
			}
		}
    return iCurrent;
    }

void
pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, int bViscosityLimitdt)
{
    int i;
    PARTICLE *p;    
    double dT,dTu;

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if(pkdIsGas(pkd, p) && TYPEQueryACTIVE(p)) {
			/*
			 * Courant condition goes here.
			 */
	       if (p->mumax>0.0) {
				if (bViscosityLimitdt) 
				  dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(p->c + 0.6*(p->c + 2*p->BalsaraSwitch*p->mumax)));
				else
				  dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(p->c + 0.6*(p->c + 2*p->mumax)));
	                   }
	       else
#if defined(PRES_HK) || defined(PRES_MONAGHAN) || defined(SIMPLESF)
			   dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(1.6*p->c+(10./FLT_MAX)));
#else
			   dT = dEtaCourant*dCosmoFac*(sqrt(0.25*p->fBall2)/(1.6*p->c));
#endif

	       if (dEtauDot > 0.0 && p->PdV < 0.0) { /* Prevent rapid adiabatic cooling */
#ifdef SSFDEBUG
			    if (p->u<=0) printf("p->iOrder: %i\n",p->iOrder);
#endif
#ifdef SSFDEBUG
			    if (p->iOrder==5514) printf("dT1 %i: %f %f %f %f\n",p->iOrder,p->u,p->uPred,p->uDot,dT);
#endif
			    assert(p->u > 0.0);
				dTu = dEtauDot*p->u/fabs(p->PdV);
				if (dTu < dT) 
					dT = dTu;
				}
#ifdef SSFDEBUG
		   if (p->iOrder==5514) printf("dT2 %i: %f %f %f %f %f %f\n",p->iOrder,p->u,p->uPred,p->uDot,dT,dTu,p->PdV);
#endif

		if(dT < p->dt)
				p->dt = dT;
			}
		}
    }

void 
pkdSphViscosityLimiter(PKD pkd, int bOn, int bShockTracker)
{
    int i;
    PARTICLE *p;    

    if (bOn) {
        for(i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];

			if(pkdIsGas(pkd, p)) {
				if (p->divv!=0.0) {         	 
					p->BalsaraSwitch = fabs(p->divv)/
						(fabs(p->divv)+sqrt(p->curlv[0]*p->curlv[0]+
											p->curlv[1]*p->curlv[1]+
											p->curlv[2]*p->curlv[2]));
					}
				else { 
					p->BalsaraSwitch = 0;
					}
				}
			}
        }
    else {
        for(i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if(pkdIsGas(pkd, p)) {
			        p->BalsaraSwitch = 1;
				}
			}
        }
    }

void pkdDensCheck(PKD pkd, int iRung, int bGreater, int iMeasure, void *data) {
#if (0)
    int i;
    struct {     
		double dMaxDensError;
		double dAvgDensError;
		int nError;
		int nTotal;
    } *tmp=data;
    double error;

    char ach[256];

    tmp->dMaxDensError=0;
    tmp->dAvgDensError=0;
    tmp->nError=0;
    tmp->nTotal=0;

    if (!iMeasure) {
        for(i = 0; i < pkdLocal(pkd); ++i) {
			if(TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
			   (pkd->pStore[i].iRung == iRung ||
				(bGreater && pkd->pStore[i].iRung > iRung))) {
				if (pkd->pStore[i].fDensity == 0) {
					sprintf(ach, "dens zero i: %d dens %g iAct %d\n",
							pkd->pStore[i].iOrder,pkd->pStore[i].fDensity,
							pkd->pStore[i].iActive);
					mdlDiag(pkd->mdl, ach);
					}
				pkd->pStore[i].fDensSave = pkd->pStore[i].fDensity;
				}
			}
		return;
        }

    for(i=0;i<pkdLocal(pkd);++i) {
		if(TYPETest(&(pkd->pStore[i]),TYPE_GAS) &&
		   (pkd->pStore[i].iRung == iRung ||
			(bGreater && pkd->pStore[i].iRung > iRung))) {
			error = abs((pkd->pStore[i].fDensSave - pkd->pStore[i].fDensity)/pkd->pStore[i].fDensity);
			tmp->dAvgDensError += error;
			tmp->nTotal++;
			if (error>tmp->dMaxDensError) tmp->dMaxDensError=error;
			if (error>1e-5) {
				tmp->nError++;
				sprintf(ach, "dens error i: %d save %g dens %g  iAct %d\n",
						pkd->pStore[i].iOrder,pkd->pStore[i].fDensSave,
						pkd->pStore[i].fDensity,pkd->pStore[i].iActive);
				mdlDiag(pkd->mdl, ach);
				}
			}
		}

    tmp->dAvgDensError/=tmp->nTotal; 
    return;
#endif
	}

#endif /* GASOLINE */

#ifdef GLASS
/* Currently wired to have no more than two regions with different
   Pressures (densities) split by x=0 with a linear connection */
void pkdGlassGasPressure(PKD pkd, void *vin)
{
    PARTICLE *p;
    double PoverRho,xx,nsp=2.5;
    int i;
    struct inGetGasPressure *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (TYPEQueryTREEACTIVE(p)) {
			if (p->r[0] < -nsp*in->dGlassxL) {
				if (p->r[0] > in->dxBoundL + nsp*in->dGlassxL)
					PoverRho=in->dGlassPoverRhoL;  
				else {
					xx =	 (p->r[0] - in->dxBoundL + nsp*in->dGlassxR)/
						(nsp*in->dGlassxL+nsp*in->dGlassxR);
			       xx = xx*xx*(-2*xx + 3);
					PoverRho = in->dGlassPoverRhoR +
						xx*(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
					}
				}
			else if (p->r[0] > nsp*in->dGlassxR) {
				if (p->r[0] < in->dxBoundR - nsp*in->dGlassxR)
					PoverRho=in->dGlassPoverRhoR;  
				else {
					xx = (p->r[0] - in->dxBoundR + nsp*in->dGlassxR)/
						(nsp*in->dGlassxL+nsp*in->dGlassxR);
					xx = xx*xx*(-2*xx + 3);
					PoverRho = in->dGlassPoverRhoR +
						xx*(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
					}
				}
			else {
				xx = (p->r[0] + nsp*in->dGlassxL)/
					(nsp*in->dGlassxL+nsp*in->dGlassxR);
				xx = xx*xx* ( -2*xx + 3 );
				PoverRho = in->dGlassPoverRhoL +
					xx*(in->dGlassPoverRhoR - in->dGlassPoverRhoL);
				}

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->gamma*PoverRho);
			}
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
			printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
		}
    }

void
pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR)
{
    int i,j;
    PARTICLE *p;  
    double v;

    for(i=0;i<pkdLocal(pkd);++i) {
        p = &pkd->pStore[i];
        if (p->r[0]<0.0) {
			for (j=0;j<3;j++) {
				p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVL); /* should reconcile with randUniform() */
				p->vPred[j]+=v;
				}
			}
        else {
			for (j=0;j<3;j++) {
				p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVR); /* should reconcile with randUniform() */
				p->vPred[j]+=v;
				}
			}
		}
    }

#endif /* GLASS */

#ifdef COLLISIONS

int
pkdNumRejects(PKD pkd)
{
	int i,nRej = 0;

	for (i=0;i<pkd->nLocal;i++)
		if (REJECT(&pkd->pStore[i])) ++nRej;

	return nRej;
	}

void
pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	SSIO ssio;
	SSDATA data;
	PARTICLE *p;
	int i,j, iSetMask;

#ifdef SSIO_USE_MPI
    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        fprintf(stdout, "Reading %d particles\n", nLocal), fflush(stdout);
    }
#endif

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization (modeled after pkdReadTipsy()).
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fBallMax = 0.0;
		}
#ifdef SSIO_USE_MPI
    /*
     *  XXX ssio group change
     *  Load it as a group file here - don't need to seek ahead as the
     *  reading is all done on a shared file pointer so the same bit
     *  shouldn't be read twice
     */
    ssioSetGroupFile();
	if (ssioOpen(pszFileName,&ssio,SSIO_READ))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
#else
	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_READ))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
#endif
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		p->iOrder = nStart + i;
		if (!pkdIsDarkByOrder(pkd,p)) /* determined by p->iOrder */
			mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
		iSetMask = TYPE_DARK;
		if (ssioData(&ssio,&data))
			mdlassert(pkd->mdl,0); /* error during read in ss file */
		p->iOrgIdx = data.org_idx;
		p->fMass = data.mass;
		p->fSoft = SOFT_FROM_SSDATA(&data);
#ifdef CHANGESOFT 
 		p->fSoft0 = p->fSoft;
#endif
		for (j=0;j<3;++j) p->r[j] = data.pos[j];
		for (j=0;j<3;++j) p->v[j] = data.vel[j];
		for (j=0;j<3;++j) p->w[j] = data.spin[j];
		p->iColor = data.color;
#ifdef NEED_VPRED
		/*
		** Note: if accelerations depend on vPred, restarting
		** from an ss file will change the integration in the
		** KDK scheme unless vPred is stored in a separate
		** file for reading.
		*/
		for (j=0;j<3;++j) p->vPred[j] = p->v[j];
#ifdef DEM
		for (j=0;j<3;++j) p->wPred[j] = p->w[j];
#endif
#endif
		TYPESet(p,iSetMask);

#ifdef GR_DRAG
		p->bNoKickNoDrift = 0;
		p->dEntryR2 = 0.0;
#endif 
#ifdef JOHNNY
		p->bInKepler = 0;
		p->bRepeat = 0;  /*Flag for repeat entries into restricted ares */
		p->sOrbElem.ecc = -1.0; /* for consistency check */
		p->dMotionPerStep = 0.0; /* don't really need to initialize this */
#endif
	}
	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */

#ifdef SSIO_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank)
    {
        fprintf(stdout, "Done\n"), fflush(stdout);
    }
    // set back to single
    ssioSetSingleFile();
#endif

#ifndef AGGS
    for (i=0;i<nLocal;++i) {
        p = &pkd->pStore[i];
	if (p->iOrgIdx < 0) {
	    fprintf(stderr, "WARNING: particle %i has negative original index (%i)\n", p->iOrder, p->iOrgIdx);
	    fprintf(stderr, "   (skipping check of other particles on this processor).\n");
	    break;
	}
    }
#endif /* !AGGS */

#ifdef SPRINGS
	/* read springs data from a separate file, if desired */
	if (pkd->bReadSpringsData) {
		char achSpringsFilename[MAXPATHLEN];
		XDR xdrs;
		FILE *fp;
		snprintf(achSpringsFilename,MAXPATHLEN,"%s.spr",pszFileName);
		if (pkd->idSelf == 0)
			printf("Reading springs data from \"%s\"...\n",achSpringsFilename);
		fp = fopen(achSpringsFilename,"r");
		mdlassert(pkd->mdl,fp != NULL);
		fseek(fp,SPRHEAD_SIZE + nStart*SPRDATA_SIZE,SEEK_SET);
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->vPred[j]);
			for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
				xdr_int(&xdrs,&p->springs[j].iOrder);
				xdr_float(&xdrs,&p->springs[j].fZeroStrainLength);
				xdr_float(&xdrs,&p->springs[j].fYoungsModulus);
				xdr_float(&xdrs,&p->springs[j].fStressLimit);
			}
		}
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* SPRINGS */

#ifdef DEM
	/* read DEM data from a separate file, if desired */
	if (pkd->bReadDEMData) {
		char achDEMFilename[MAXPATHLEN];
		XDR xdrs;
		FILE *fp;
		int k;
		snprintf(achDEMFilename,MAXPATHLEN,"%s.dem",pszFileName);
		if (pkd->idSelf == 0)
			printf("Reading DEM data from \"%s\"...\n",achDEMFilename);
		fp = fopen(achDEMFilename,"r");
		mdlassert(pkd->mdl,fp != NULL);
		fseek(fp,DEMHEAD_SIZE + nStart*DEMDATA_SIZE,SEEK_SET);
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->vPred[j]);
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->wPred[j]);
			for (j=0;j<MAX_NUM_OVERLAPS_PER_PARTICLE;j++) {
				xdr_int(&xdrs,&p->overlaps[j].iOrder);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vShear[k]);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vnOld[k]);
#ifdef DEM_ROTATION_DASHPOT
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vRoll[k]);
				xdr_double(&xdrs,&p->overlaps[j].dTwist);
#endif /* DEM_ROTATION_DASHPOT */
				xdr_long(&xdrs,&p->overlaps[j].liOverlapCounter);
				}
#ifdef WALLS
			for (j=0;j<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;j++) {
				xdr_int(&xdrs,&p->walloverlaps[j].iOrder);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vShear[k]);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vnOld[k]);
#ifdef DEM_ROTATION_DASHPOT
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vRoll[k]);
				xdr_double(&xdrs,&p->walloverlaps[j].dTwist);
#endif /* DEM_ROTATION_DASHPOT */
				xdr_long(&xdrs,&p->walloverlaps[j].liOverlapCounter);
				}
#endif /* WALLS */
			}
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#ifdef DEM_TRACK_ORIENT
	/* read orientation data from a separate file, if desired */
	if (pkd->bReadOrient) {
		ORIDATA ori;
		Vector vTmp;
		Matrix mTmp;
		char achOriFilename[MAXPATHLEN];
		snprintf(achOriFilename,MAXPATHLEN,"%s.ori",pszFileName);
		if (ssioOpen(achOriFilename,&ssio,SSIO_READ))
			mdlassert(pkd->mdl,0); /* unable to open .ori file */
		if (ssioSetPos(&ssio,ORIHEAD_SIZE + nStart*ORIDATA_SIZE))
			mdlassert(pkd->mdl,0); /* unable to seek in .ori file */
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			if (ssioOriData(&ssio,&ori))
				mdlassert(pkd->mdl,0); /* error during read in .ori file */
			vectorCross(ori.p1,ori.p2,vTmp); /* p3 = p1 x p2 */
			vectorCopy(ori.p1,mTmp[0]); /* load temporary matrix by rows */
			vectorCopy(ori.p2,mTmp[1]);
			vectorCopy(vTmp,mTmp[2]);
			matrixTranspose(mTmp,p->orient); /* principal axes in columns */
			}
		}
	/* otherwise initialize orientations to identity matrix */
	else {
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			matrixIdentity(p->orient);
			}
		}
#endif /* DEM_TRACK_ORIENT */
#endif /* DEM */

#ifdef CHARGE
	pkdChargeInit(pkd); /* assign charges now!  (this is a bit of a hack...) */
#endif
	}

void
pkdWriteSS(PKD pkd,char *pszFileName,int nStart,int bReduced)
{
	SSIO ssio;
	PARTICLE *p;
	int i,j;

#ifdef SSIO_USE_MPI
	assert(bReduced == 0); /* reduced files parallel I/O not supported yet */
    ssioSetGroupFile();
	if (ssioOpen(pszFileName,&ssio,SSIO_WRITE))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
#else
	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_UPDATE))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
#endif /* !SSIO_USE_MPI */

	if (bReduced) {
		SSRDATA data;
#ifndef SSIO_USE_MPI
		/* seek */
		if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSRDATA_SIZE))
			mdlassert(pkd->mdl,0); /* unable to seek in reduced ss file */
#endif
		/* write */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (!pkdIsDark(pkd,p))
				mdlassert(pkd->mdl,0); /* only dark particles allowed in reduced ss file */
			data.fMass = p->fMass;
			data.fRadius = RADIUS(p);
			for (j=0;j<3;j++)
				data.vPos[j] = p->r[j];
			data.iColor = p->iColor;
			data.iOrgIdx = p->iOrgIdx;
			if (ssioDataReduced(&ssio,&data))
				mdlassert(pkd->mdl,0); /* unable to write in reduced ss file */
			}
		}
	else {
		SSDATA data;
#ifndef SSIO_USE_MPI
		/* seek */
		if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
			mdlassert(pkd->mdl,0); /* unable to seek in ss file */
#endif
		/* write */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (!pkdIsDark(pkd,p))
				mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
			data.org_idx = p->iOrgIdx;
			data.mass = p->fMass;
			data.radius = RADIUS(p);
			for (j=0;j<3;++j) data.pos[j]  = p->r[j];
			for (j=0;j<3;++j) data.vel[j]  = p->v[j];
			for (j=0;j<3;++j) data.spin[j] = p->w[j];
			data.color = p->iColor;
			if (ssioData(&ssio,&data))
				mdlassert(pkd->mdl,0); /* unable to write in ss file */
			}
		}

	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */

#ifdef SSIO_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    // set back to single
    ssioSetSingleFile();
#endif

#ifdef SPRINGS
	/* write springs data to a separate file */
	if (!bReduced && !pkd->bNoSpringsOutput) {
		char achSpringsFilename[256];
		XDR xdrs;
		FILE *fp;
		snprintf(achSpringsFilename,256,"%s.spr",pszFileName);
		fp = fopen(achSpringsFilename,"r+");
		mdlassert(pkd->mdl,fp != NULL);
		fseek(fp,SPRHEAD_SIZE + nStart*SPRDATA_SIZE,SEEK_SET);
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->vPred[j]);
			for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
				xdr_int(&xdrs,&p->springs[j].iOrder);
				xdr_float(&xdrs,&p->springs[j].fZeroStrainLength);
				xdr_float(&xdrs,&p->springs[j].fYoungsModulus);
				xdr_float(&xdrs,&p->springs[j].fStressLimit);
				}
			}
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* SPRINGS */

#ifdef DEM
	/* write DEM data to a separate file  */
	if (!bReduced && pkd->nDEMOutputs != 0) {
		char achDEMFilename[256];
		int k;
		XDR xdrs;
		FILE *fp;
		snprintf(achDEMFilename,256,"%s.dem",pszFileName);
		fp = fopen(achDEMFilename,"r+");
		mdlassert(pkd->mdl,fp != NULL);
		fseek(fp,DEMHEAD_SIZE + nStart*DEMDATA_SIZE,SEEK_SET);
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		for (i=0;i<pkdLocal(pkd);i++) {
			p = &pkd->pStore[i];
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->vPred[j]);
			for (j=0;j<3;++j) xdr_double(&xdrs,&p->wPred[j]);
			for (j=0;j<MAX_NUM_OVERLAPS_PER_PARTICLE;j++) {
				xdr_int(&xdrs,&p->overlaps[j].iOrder);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vShear[k]);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vnOld[k]);
#ifdef DEM_ROTATION_DASHPOT
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->overlaps[j].vRoll[k]);
                                xdr_double(&xdrs,&p->overlaps[j].dTwist);
#endif /* DEM_ROTATION_DASHPOT */
				xdr_long(&xdrs,&p->overlaps[j].liOverlapCounter);
				}
#ifdef WALLS
			for (j=0;j<MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;j++) {
				xdr_int(&xdrs,&p->walloverlaps[j].iOrder);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vShear[k]);
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vnOld[k]);
#ifdef DEM_ROTATION_DASHPOT
				for (k=0;k<3;k++) xdr_double(&xdrs,&p->walloverlaps[j].vRoll[k]);
                                xdr_double(&xdrs,&p->walloverlaps[j].dTwist);
#endif /* DEM_ROTATION_DASHPOT */
				xdr_long(&xdrs,&p->walloverlaps[j].liOverlapCounter);
				}
#endif /* WALLS */
			}
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* DEM */

#ifdef DEM_TRACK_ORIENT
	/* write orientation data to a separate file */
	{
		union {
		  ORIDATA ori;
		  ORIRDATA orir;
		  } oridata;
		char achOriFilename[512/*PST_FILENAME_SIZE*/];
		int k;
		snprintf(achOriFilename,512,"%s.ori",pszFileName);
		if (ssioOpen(achOriFilename,&ssio,SSIO_UPDATE))
			mdlassert(pkd->mdl,0); /* unable to open .ori file */
		if (bReduced) {
			if (ssioSetPos(&ssio,ORIHEAD_SIZE + nStart*ORIRDATA_SIZE))
				mdlassert(pkd->mdl,0); /* unable to seek in reduced .ori file */
			}
		else {
			if (ssioSetPos(&ssio,ORIHEAD_SIZE + nStart*ORIDATA_SIZE))
				mdlassert(pkd->mdl,0); /* unable to seek in .ori file */
			}
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (bReduced) {
				for (k=0;k<3;k++) {
					oridata.orir.p1[k] = p->orient[k][0]; /* axes in columns */
					oridata.orir.p2[k] = p->orient[k][1];
					}
				if (ssioOriDataReduced(&ssio,&oridata.orir))
					mdlassert(pkd->mdl,0); /* unable to write in reduced .ori file */
				}
			else {
				for (k=0;k<3;k++) {
					oridata.ori.p1[k] = p->orient[k][0]; /* axes in columns */
					oridata.ori.p2[k] = p->orient[k][1];
					}
				if (ssioOriData(&ssio,&oridata.ori))
					mdlassert(pkd->mdl,0); /* unable to write in .ori file */
				}
			}
		if (ssioClose(&ssio))
			mdlassert(pkd->mdl,0); /* unable to close .ori file */
		}
#endif /* DEM_TRACK_ORIENT */
	}

void pkdAddUnifGrav(PKD pkd,double dgx,double dgy,double dgz)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
		if (TYPEQueryACTIVE(&p[i])) {
			p[i].a[0] += dgx;
			p[i].a[1] += dgy;
			p[i].a[2] += dgz;
			}
		}
	}

#ifdef JOHNNY /*DEBUG will want to make this standard at some point*/

void _pkdReplaceEscaper(PKD pkd,double dCentMass,double dStarMass,PARTICLE *p)
{
	/*
	** Replaces escaped star with star from a set distribution.
	** We generate orbital elements for the star around an SMBH, 
	** drawing the semimajor axis from a power law distribution with
	** index -alpha, the eccentricity from a thermal distribution P(e)=2e, 
	** the true anomaly using a proper time weighting, and the other 
	** angles from an isotropic distribution.  These are then changed to 
	** Cartesian coordinates.
	*/

	double Minternal,avar,prob,probmax,xrand;
	double a,e,inc,Om,om,nu;
	double E,h,R,param;
	double amin = 100;   /* Minimum allowed semimajor axis in AU */
	double amax;   /* Maximum allowed semimajor axis in AU */
	double numStars, rinfl;

	rinfl = 2.063e5*sqrt(dCentMass/1.0e6);
	numStars = dStarMass;
	/* amax = rinfl*pow(numStars/dCentMass,1.0/(3.0-alpha)); */
	amax = rinfl*(numStars/dCentMass);/* stream-lined version (alpha = 2) */

	avar=randUniform();
	do {
		/*a=pow((pow(amin,3.0-alpha)+avar*(pow(amax,3.0-alpha)-pow(amin,3.0-alpha))),1.0/(3.0-alpha));*/ /* general version */
		a=amin+avar*(amax-amin); /* stream-lined version (alpha = 2) */
		Minternal=dCentMass+dStarMass*avar;
		e=sqrt(randUniform());
		} while (a*(1.0-e)<10.0); /* pericenter must be > 10 AU */
	probmax=(1.0+e)*(1.0+e);
	do {
		nu=2.0*M_PI*randUniform();
		prob=(1.0-e*e)/(1.0+e*cos(nu));
		prob=prob*prob;
		xrand=randUniform();
		} while (xrand>prob/probmax);
	Om=2.0*M_PI*randUniform();
	om=2.0*M_PI*randUniform();
	inc=acos(1.0-2.0*randUniform());
	if (randUniform()<0.5) inc*=-1.0;

	E=acos((e+cos(nu))/(1.0+e*cos(nu)));
	if (nu>M_PI) E=2.0*M_PI-E;

	R=a*(1.0-e*cos(E));
	h=sqrt(Minternal*a*(1.0-e*e)); /* G = 1 */
	param=(1.0+e*cos(nu))*a*(1.0-e*cos(E));
	p->r[0]=R*(cos(Om)*cos(om+nu)-sin(Om)*sin(om+nu)*cos(inc));
	p->r[1]=R*(sin(Om)*cos(om+nu)+cos(Om)*sin(om+nu)*cos(inc));
	p->r[2]=R*(sin(inc)*sin(om+nu));

	p->v[0]= p->r[0]*h*e*sin(nu)/(R*param);
	p->v[0]-=(h/R)*(cos(Om)*sin(om+nu)+sin(Om)*cos(om+nu)*cos(inc));
	p->v[1]= p->r[1]*h*e*sin(nu)/(R*param);
	p->v[1]-=(h/R)*(sin(Om)*sin(om+nu)-cos(Om)*cos(om+nu)*cos(inc));
	p->v[2]= p->r[2]*h*e*sin(nu)/(R*param);
	p->v[2]+=(h/R)*sin(inc)*cos(om+nu);

	}

int pkdFindEscapers(PKD pkd,double dCentMass,double dStarMass)
{
	/*
	** Marks for deletion any particle with non-negative energy relative
	** to the central mass (intended for heliocentric coordinate system).
	** Here the central mass is the mass of the SMBH; the "star mass" is
	** the mass of everything else in the system, needed for the
	** distribution function (cf. _pkdReplaceEscaper()).
	*/

	PARTICLE *p;
	double r2, rinfl, numStars, amax;
	int i,n,nEsc = 0;

	rinfl = 2.063e5*sqrt(dCentMass/1.0e6);
	numStars = dStarMass;
	/* amax = rinfl*pow(numStars/dCentMass,1.0/(3.0-alpha)); */
	amax = rinfl*(numStars/dCentMass);/* stream-lined version (alpha = 2) */

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
		r2 = p[i].r[0]*p[i].r[0] + p[i].r[1]*p[i].r[1] + p[i].r[2]*p[i].r[2];
		mdlassert(pkd->mdl,r2 > 0.0);
		if (sqrt(r2) >= 6*amax) {
			printf("Particle %i (processor %i, index %i) escaped!\n",
				   p[i].iOrder,pkd->idSelf,i);
			_pkdReplaceEscaper(pkd,dCentMass,dStarMass,&p[i]);
			++nEsc;
			}
		}

	return nEsc;
	}

int pkdCheckForKepler(PKD pkd,double dCentMass,double dDelta)
{
	/*
	** Checks to see whether any particles satisfy the criterion for
	** switching to (or from) a Kepler orbit relative to the central
	** mass.  Here we are assuming the central mass is a supermassive
	** black hole located at the origin (0,0,0), i.e., we are using
	** "heliocentric" coordinates.  Particles switching to a Kepler
	** orbit have their orbital parameters computed and a flag set to
	** disable drifts and kicks.  NOTE: this function should be called
	** at the end of the closing kick.  Particles in a Kepler orbit
	** are moved forward by the interval dDelta and have their
	** heliocentric coordinates computed and stored as r & v in the
	** particle structure.  Those leaving a Kepler orbit have their
	** flag reset so drifts and kicks resume normally.  Returns the
	** number of particles entering the Kepler zone on this processor.
	*/

	const double r2max = 2500.0; /* 50 AU squared: hardcoded limit */

	struct helio sCoords;
	PARTICLE *p;
	double r2;
	int i,n,nKep = 0;

	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
		p = &pkd->pStore[i];
		if (!p->bInKepler) {
		        r2 = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2];
			if (r2 < r2max) { /* particle was outside critical radius on previous step */
				/* fill struct and get orbital elements */
				sCoords.x = p->r[0];
				sCoords.y = p->r[1];
				sCoords.z = p->r[2];
				sCoords.vx = p->v[0];
				sCoords.vy = p->v[1];
				sCoords.vz = p->v[2];
				heltodel(dCentMass + p->fMass,&sCoords,&p->sOrbElem);
				assert(p->sOrbElem.sma > 0.0); /* only closed orbits allowed! */
				p->dMotionPerStep = sqrt((dCentMass + p->fMass)/(p->sOrbElem.sma*p->sOrbElem.sma*p->sOrbElem.sma))*dDelta;
				p->bInKepler = 1; /* set flag */
				if (!p->bRepeat) { /* First time particle has entered critical radius */
				  ++nKep;
				}
				p->bRepeat += 1; /* count number of entries into critical radius*/
				/* printf("Particle %i has entered critical radius %i time(s).\n",
				   p->iOrder,p->bRepeat); */
			}
			else { p->bInKepler = 0; }
		}
		else {
				assert(p->sOrbElem.ecc >= 0.0); /* sanity check */
				/* advance mean anomaly */
				p->sOrbElem.mea += p->dMotionPerStep;
				p->sOrbElem.mea = fmod(p->sOrbElem.mea,2*M_PI); /* Make sure between 0 and 2pi */
				/* get new heliocentric coords and store in particle struct */
				deltohel(dCentMass + p->fMass,&p->sOrbElem,&sCoords);
				p->r[0] = sCoords.x;
				p->r[1] = sCoords.y;
				p->r[2] = sCoords.z;
				p->v[0] = sCoords.vx;
				p->v[1] = sCoords.vy;
				p->v[2] = sCoords.vz;
		                r2 = p->r[0]*p->r[0] + p->r[1]*p->r[1] + p->r[2]*p->r[2];
				 if (r2 > r2max) {
				  p->bInKepler = 0;
				  }
		}
	}
	return nKep;
	}

#endif /* JOHNNY */

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

void pkdPatch(PKD pkd)
{
	PARTICLE *p;
	PATCH_PARAMS *PP;
	int i;

	PP = pkd->PP;

	assert(PP->bPatch);

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		/*
		** Apply the velocity-independent part of Hill's Equations
		** (Kick2 in the symplectic integration paper)
		*/
		p->a[0] -= PP->dOrbFreq*PP->dOrbFreq*p->r[0];
		p->a[2] -= PP->dOrbFreqZ2*p->r[2];
		}

	if (PP->bExtPert) {

		/* following only for external perturbers in the patch model */

		double rPert[3];      /* position of perturber relative to saturn */
		double rpp[3];        /* position of patch relative to perturber */
		double R[3];          /* position of a patch particle relative to perturber (goal) */
		double Rmag,Rmag3inv; /* magnitude of R, and inverse cube */

		/*
		** NOTE: positive x: away from saturn, y: direction of orbital velocity, z: RH rule
		** NOTE: All vectors are calculated in the *patch* frame (rotating wrt saturn)
		*/

		/* calculate position of perturber relative to saturn, in the patch frame,
		** by using a rotation matrix:
		**       cos(omega_patch * t)        sin(omega_patch * t)
		**      -sin(omega_patch * t)        cos(omega_patch * t)
		*/

		rPert[0] = (PP->dPertOrbDist*cos((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*cos(PP->dOrbFreq*pkd->dTime))
		  + (PP->dPertOrbDist*sin((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*sin(PP->dOrbFreq*pkd->dTime));

		rPert[1] = (PP->dPertOrbDist*cos((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*(-sin(PP->dOrbFreq*pkd->dTime)))
		  + (PP->dPertOrbDist*sin((PP->dPertOrbFreq*pkd->dTime) + PP->dPertPhase)*cos(PP->dOrbFreq*pkd->dTime));

		rPert[2] = PP->dPertMaxZ*cos((PP->dPertOrbFreqZ*pkd->dTime) + PP->dPertPhaseZ);

		/* calculate vector pointing from perturber (center) to patch center */
		rpp[0] = PP->dOrbDist - rPert[0];
		rpp[1] = -rPert[1];
		rpp[2] = -rPert[2];

		for (i=0;i<pkd->nLocal;i++) {
			p = &pkd->pStore[i];
			if (!TYPEQueryACTIVE(p)) continue;

			/* calculate force/mass of an external perturber on sliding-patch particles */

			/* calculate position vector */
			R[0] = rpp[0] + p->r[0];
			R[1] = rpp[1] + p->r[1];
			R[2] = rpp[2] + p->r[2];

			Rmag = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
			Rmag3inv = 1.0/(Rmag*Rmag*Rmag);

			p->a[0] -= PP->dPertMass*R[0]*Rmag3inv;
			p->a[1] -= PP->dPertMass*R[1]*Rmag3inv;
			p->a[2] -= PP->dPertMass*R[2]*Rmag3inv;
			}
		}
	}

int pkdRandAzWrap(PKD pkd)
{
	/*
	** This routine probably only valid for low surface densities
	** that do not lead to strong wake/aggregate formation.
	*/

	/* NOTE: Currently ignoring spin... */

	const double Sqrt2OverPi = M_2_SQRTPI/M_SQRT2; /* sqrt(2/PI) */

	PATCH_PARAMS *PP;
	PARTICLE *p;
	double amp,phase;
	int i,iSide,nRand = 0;

	assert(pkd->nLocal > 0);

	PP = pkd->PP;

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];
		TYPESet(p,TYPE_TREEACTIVE); /* all particles go into density tree */
		TYPEReset(p,TYPE_SMOOTHACTIVE); /* not all are centers for overlap searches */
		if (p->bAzWrap) {
			/* assign uniform random radial position depending on side */
			iSide = (p->r[0] < 0.0 ? -1 : p->r[0] > 0.0 ? 1 : 0); /* -1 ==> left side of patch; +1 ==> right side */
			if (iSide == 0)
				iSide = (randUniform() < 0.5 ? -1 : 1); /* handle rare case */
			/* flip side if equilibrium stream on opposite side */
			if ((iSide == -1 && PP->iStripOption == STRIP_RIGHT_ONLY) ||
				(iSide ==  1 && PP->iStripOption == STRIP_LEFT_ONLY))
				iSide *= -1;
			/*
			** sqrt() below is to account for the fact that particles
			** at the outer edge stream by faster and should therefore
			** wrap more often. If we don't account for this, particle
			** density at the inner edge grows.  Uses probability
			** function from origin to patch edge, but constructs the
			** transformation on the strip only
			*/
#ifndef SQ
#define SQ(x) ((x)*(x))
#endif
			if (!PP->bNoRandomX) {
				/* Randomize radial location */
				p->r[0] = iSide*sqrt(SQ(PP->dStripInner) +
									 randUniform()*(SQ(PP->dStripOuter) - SQ(PP->dStripInner)));
				}
			else {
				/* Assign new radial location based on particle's guiding center */
				p->r[0] = 4.0*p->r[0] + ( 2.0*(p->v[1])/(PP->dOrbFreq));
				// OLD VERSION:   p->r[0] = p->r[0] + ( (p->v[1] /*+ 1.5*PP->dOrbFreq*p->r[0]*/)/(2.0*PP->dOrbFreq));
				/* check if particle now outside patch; if so, fix it */
				if (p->r[0] < -0.5*PP->dWidth) {
					p->r[0] += PP->dWidth;
					printf("DEBUG-rp: had to wrap radially during guiding center calculation (+)\n");
					}
				else if (p->r[0] > 0.5*PP->dWidth) {
					p->r[0] -= PP->dWidth;
					printf("DEBUG-rp: had to wrap radially during guiding center calculation (-)\n");
					}
				assert(p->r[0] >= -0.5*PP->dWidth && p->r[0] <= 0.5*PP->dWidth);
				} /* radial randomization option */
#undef SQ
			/* old version...
			   p->r[0] = iSide*(PP->dStripInner + (PP->dStripOuter - PP->dStripInner)*randUniform());
			*/
			/* assign random velocity relative to shear */
			p->v[0] = randGaussian()*PP->dVelDispX*sqrt(PP->dAvgMass/p->fMass);
			p->v[1] = randGaussian()*PP->dVelDispY*sqrt(PP->dAvgMass/p->fMass);
			p->v[1] -= 1.5*PP->dOrbFreq*p->r[0]; /* add shear */
			/* deal with the third dimension */
			amp = PP->dAvgVertAmp*Sqrt2OverPi*sqrt(-2.0*log(randUniform()))*sqrt(PP->dAvgMass/p->fMass); /* Rayleigh deviate */
			phase = 2.0*M_PI*randUniform();
			p->r[2] = amp*cos(phase);
			p->v[2] = - amp*PP->dOrbFreq*sin(phase);
			/* set cannonical momentum for new symplectic patch integrator */
			p->dPy = p->v[1] + 2.0*PP->dOrbFreq*p->r[0];
			/* set flag to do an overlap search around this particle */
			TYPESet(p,TYPE_SMOOTHACTIVE);
			++nRand;
			if ((iSide == -1 && p->r[1] > -0.4*PP->dLength) ||
				(iSide == 1 && p->r[1] < 0.4*PP->dLength))
				(void) fprintf(stderr,"DEBUG iOrd=%i iSide=%i y=%g (Ly=%g)\n",p->iOrder,iSide,p->r[1],PP->dLength);
			}
		}

	return nRand;
	}

#endif /* SLIDING_PATCH */

#ifdef SIMPLE_GAS_DRAG
void
_pkdSimpleGasVel(double r[],int iFlowOpt,double dTime,double u[])
{
	/* From Paolo Tanga's e-mail Feb 19, 2001: M, f_1 and f_2 are
	   parameters (double precision) that can be considered internal
	   to _pkdGasVel, since they define the streamfunction. It define
	   the attached velocity field. f_1 and f_2 are scaling parameters
	   that allow to control the size of the vortex and its
	   width/length ratio.  Note that I have centered a vortex in 0,0,
	   which should be good in the local coordinates of a single
	   patch. Anyway, the given function is not suitable for the
	   simulation of an entire disk...! */

	switch (iFlowOpt) {
	case 1: {
		const double M	= 1.0e-7;
		const double f1	= 1.0e7;
		const double f2	= 1.0e7;
#ifndef SLIDING_PATCH
		assert(0); /* this formula only valid in patch frame */
#endif

		u[0] = M*cos(f2*r[1] - M_PI_2)*cos(f1*r[0] - M_PI);
		u[1] = M*sin(f2*r[1] - M_PI_2)*sin(f1*r[0] - M_PI);
		u[2] = 0.0;
		break;
		}
	default:
		assert(0); /* no other cases defined */
		}
}

void pkdSimpleGasDrag(PKD pkd,int iFlowOpt,int bEpstein,double dGamma,
					  double dTime)
{
	PARTICLE *p;
	double u[3],dCoef;
	int i;

	for (i=0;i<pkd->nLocal;i++) {
		p = &pkd->pStore[i];
		if (!TYPEQueryACTIVE(p)) continue;
		_pkdSimpleGasVel(p->r,iFlowOpt,dTime,u); /* get gas velocity */
		if (bEpstein)
			dCoef = -dGamma/RADIUS(p); /* i.e. -gamma/R */
		else
			dCoef = -dGamma/(RADIUS(p)*RADIUS(p)); /* i.e. -gamma/R^2 */
		/* note: following assumes dust density << gas density */
		p->a[0] += dCoef*(p->vPred[0] - u[0]);
		p->a[1] += dCoef*(p->vPred[1] - u[1]);
		p->a[2] += dCoef*(p->vPred[2] - u[2]);
		}
	}
#endif /* SIMPLE_GAS_DRAG */

#ifdef NEED_VPRED
#ifdef GASOLINE
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo,double duDelta,
			 int iGasModel,double z,double duDotLimit)
{
	PARTICLE *p;
	int i,j,n;

	mdlDiag(pkd->mdl, "Into Vpred\n");

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
		if (pkdIsGas(pkd,p)) {
			for (j=0;j<3;++j) {
				p->vPred[j] = p->vPred[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
			if (iGasModel != GASMODEL_ISOTHERMAL) {
#ifndef NOCOOLING
#ifdef COOLDEBUG
				if (p->uPred+p->uDot*duDelta < 0) 
					fprintf(stderr,"upred error %i: %g %g %g -> %g %i\n",p->iOrder,p->uPred,p->uDot,duDelta,p->uPred + p->uDot*duDelta,p->iRung);
#endif
#ifdef SSFDEBUG
			    if (p->iOrder==5514) printf("%i: %f %f %f %f\n",p->iOrder,p->u,p->uPred,p->uDot,duDelta);
#endif
			  p->uPred = p->uPred + p->uDot*duDelta;
#else
			  p->uPred = p->uPred + p->PdV*duDelta;
#endif
#if defined(PRES_HK) || defined(PRES_MONAGHAN) || defined(SIMPLESF)
			  if (p->uPred < 0) p->uPred = 0;
#endif
			  mdlassert(pkd->mdl,p->uPred >= 0.0);
			  }
			}
		}

	pkdStopTimer(pkd,1);
	mdlDiag(pkd->mdl, "Done Vpred\n");
	}
#else
void
pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;i++) {
#ifdef WALLS
		if (PARTICLE_STUCK(&p[i]))
			continue;
#endif
		for (j=0;j<3;j++) {
			p[i].vPred[j] = p[i].vPred[j]*dvFacOne + p[i].a[j]*dvFacTwo;
#ifdef DEM
			p[i].wPred[j] = p[i].wPred[j]*dvFacOne + p[i].wDot[j]*dvFacTwo;
#endif /* DEM */
			}
		}
	}
#endif
#endif /* NEED_VPRED */

void
pkdCOM(PKD pkd, double *com)
{
    int i;
    int nLocal = pkdLocal(pkd);
	double m;
    
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = 0;
	com[4] = 0;
	com[5] = 0;
	com[6] = 0;
	com[7] = 0;
	com[8] = 0;
	com[9] = 0;
	com[10] = 0;
	com[11] = 0;

    for (i=0;i<nLocal;++i) {
	  m = pkd->pStore[i].fMass;
	  if ( TYPETest(&pkd->pStore[i], TYPE_GAS) ) {
		com[0] += m*pkd->pStore[i].r[0];
		com[1] += m*pkd->pStore[i].r[1];
		com[2] += m*pkd->pStore[i].r[2];
		com[3] += m;
		}
	  else if ( TYPETest(&pkd->pStore[i], TYPE_DARK) ) {
		com[4] += m*pkd->pStore[i].r[0];
		com[5] += m*pkd->pStore[i].r[1];
		com[6] += m*pkd->pStore[i].r[2];
		com[7] += m;
		}
	  else if ( TYPETest(&pkd->pStore[i], TYPE_STAR) ) {
		com[8] += m*pkd->pStore[i].r[0];
		com[9] += m*pkd->pStore[i].r[1];
		com[10] += m*pkd->pStore[i].r[2];
		com[11] += m;
		}
  	  }
    }

void
pkdCOMByType(PKD pkd, int type, double *com)
{
    int i;
    int nLocal = pkdLocal(pkd);
	double m;
    
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = 0;

    for (i=0;i<nLocal;++i) {
	  if ( TYPETest(&pkd->pStore[i], type) ) {
		m = pkd->pStore[i].fMass;
		com[0] += m*pkd->pStore[i].r[0];
		com[1] += m*pkd->pStore[i].r[1];
		com[2] += m*pkd->pStore[i].r[2];
		com[3] += m;
		}
  	  }
    }

void
pkdOldestStar(PKD pkd, double *com)
{
#if defined(STARFORM) || defined (SIMPLESF)
    int nLocal = pkdLocal(pkd);
#endif

	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	com[3] = FLT_MAX;

#if defined(STARFORM) || defined (SIMPLESF)
    for (i=0;i<nLocal;++i) {
	  if ( TYPETest(&pkd->pStore[i], TYPE_STAR) && pkd->pStore[i].fTimeForm < com[3]) {
		com[0] = pkd->pStore[i].r[0];
		com[1] = pkd->pStore[i].r[1];
		com[2] = pkd->pStore[i].r[2];
		com[3] = pkd->pStore[i].fTimeForm;
		}
  	  }
#endif
   }

int pkdSetSink(PKD pkd, double dSinkMassMin)
{
#ifdef GASOLINE
    PARTICLE *p;
    int i,nSink = 0;
    int nLocal = pkdLocal(pkd);

    for(i=0;i<nLocal;++i) { 
		p = &pkd->pStore[i];
		if ((TYPETest(p,TYPE_STAR) && p->fTimeForm < 0)) {
			TYPESet(p,TYPE_SINK);
			nSink++;
			}
		}

    return nSink;
#else
    assert(0);	/* GASOLINE needed for sink particles to work */
#endif
	return -1; /* to keep compiler happy */
    }

#ifdef CHARGE

void pkdChargeInit(PKD pkd)
{
	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		/* for now, assign charge between 0 & 1 billion e- based on color */
		p->dCharge = ((double) (p->iColor - 16)/(255 - 16))*1.0e9;
		}
	}

void pkdChargeZGetMoments(PKD pkd,CHARGE_PARAMS *CP)
{
	PARTICLE *p;
	int i,nLocal = pkdLocal(pkd);

	CP->dQ = CP->dZ = CP->dSigZ = 0.0;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		CP->dQ += p->dCharge;
		CP->dZ += p->dCharge*p->r[2];
		CP->dSigZ += p->dCharge*p->r[2]*p->r[2];
		}
	}

void pkdChargeZApplyMoments(PKD pkd,const CHARGE_PARAMS *CP)
{
	/* NOTE: units of permittivity = C^2 s^2/(m^3 kg) */

	const double eps_conv = 9.107e76; /* permittivity in pkdgrav units! */

	PARTICLE *p;
	double dPreFac;
	int i,nLocal = pkdLocal(pkd);

	double dLim = 1.5*CP->dSigZ;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		dPreFac = CP->dQ*p->dCharge/(eps_conv*CP->dArea*p->fMass);
		if (p->r[2] > CP->dZ + dLim)
			p->a[2] += dPreFac;
		else if (p->r[2] < CP->dZ - dLim)
			p->a[2] -= dPreFac;
		else
			/* may be unstable if dLim small, but this can be fixed if necessary */
			p->a[2] += dPreFac*(p->r[2] - CP->dZ)/dLim;
		}
	}

#endif /* CHARGE */

#ifdef DEM_TIDAL_SPACE

void pkdDEMTidalFindPlanet(PKD pkd,int iColor,int *bFound,double dPos[],double dVel[])
{
	PARTICLE *p;
	int i,k,nLocal = pkdLocal(pkd);

	*bFound = 0;
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->iColor == iColor) {
			for (k=0;k<3;k++) {
				dPos[k] = p->r[k];
				dVel[k] = p->v[k];
				}
			*bFound = 1;
			return; /* assuming only 1 planet particle! */
			}
		}
	}

void pkdDEMTidalFindMarker(PKD pkd,int iOrder,int *bFound,double dPos[],double dVel[])
{
	PARTICLE *p;
	int i,k,nLocal = pkdLocal(pkd);

	*bFound = 0;
	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		if (p->iOrder == iOrder) {
			for (k=0;k<3;k++) {
				dPos[k] = p->r[k];
				dVel[k] = p->v[k];
				}
			*bFound = 1;
			return; /* assuming iOrder is unique! */
			}
		}
	}

#endif /* DEM_TIDAL_SPACE */

#ifdef DEM_TIDAL_LOCAL

/* frames: SPT - space translating frame, where all Agg and Planet vectors are represented.
           MCT - mass center translating frame, parallel with SPT, origin at the mass center of rubble pile.
           BDY - body fixed frame, origin at mass center, initially coincided with MCT.
           LOC - local frame at the surface, origin at the position of local simulation. */

void getMatrixBDY2LOC(const Vector vlocalOrig0,Matrix MatrixBDY2LOC)
{
    Vector vlocalOrient;
    double phi,theta,Cth,Sth,Cph,Sph;
    
	/*generate orientation vector of the local subflat here due to the topography data, here a simple strategy is used below, working as a test*/
    vectorCopy(vlocalOrig0,vlocalOrient);
    vectorNorm(vlocalOrient);
	/*calculate 1- 2- rotation from the body fixed frame to local frame*/
    phi=asin(vlocalOrient[0]);//phi is ranged in [-pi/2,pi/2]
    if (vlocalOrient[1]==0.0&&vlocalOrient[2]==0.0) theta=0.0;
    else theta=atan2(vlocalOrient[2],-vlocalOrient[1]);//theta is ranged in [-pi,pi]
    /*generate transforming matrix from body frame to local frame*/
    Cth=cos(theta);
    Sth=sin(theta);
    Cph=cos(phi);
    Sph=sin(phi);
    MatrixBDY2LOC[0][0]=Cph;
    MatrixBDY2LOC[0][1]=Sph*Sth;
    MatrixBDY2LOC[0][2]=-Sph*Cth;
    MatrixBDY2LOC[1][0]=0.0;
    MatrixBDY2LOC[1][1]=Cth;
    MatrixBDY2LOC[1][2]=Sth;
    MatrixBDY2LOC[2][0]=Sph;
    MatrixBDY2LOC[2][1]=-Cph*Sth;
    MatrixBDY2LOC[2][2]=Cph*Sth;
}

void getMatrixOrient(const DEM_TIDAL *d,Matrix m)
{
    Vector r1,r2,r3;
    int i;

    vectorSub(d->vMarker1Pos,d->vAggPos,r1);
	vectorSub(d->vMarker2Pos,d->vAggPos,r2);
	vectorSub(d->vMarker3Pos,d->vAggPos,r3);
    
    for (i=0;i<3;i++)
    {
        m[i][0]=r1[i];
        m[i][1]=r2[i];
        m[i][2]=r3[i];
    }
}

void getPlntGrav(const Vector r,Vector plntGrav)
{
    double rNorm;

    rNorm = vectorMag(r);
    vectorScale(r,PLANET_MASS/(rNorm*rNorm*rNorm),plntGrav);
}

void pkdDEMTidal(PKD pkd,const DEM_TIDAL *d0,const DEM_TIDAL *d)
{
	PARTICLE *p;
	Vector vLocalOrig,vLocalGrav,vLocalSurf,vPlanetPos,plntGrav,tmpVector,r,v,a,vAggSpin,vAggSpinDot,vAggAcc;
	Matrix MatrixBDY2LOC,MatrixMCT2BDY,MatrixOrient0,MatrixOrient,tmpMatrix,MatrixMCT2LOC;
	int i,nLocal = pkdLocal(pkd);

    /* set the value of initial local origin and local gravity, can be out of the loop*/
    vectorSet(vLocalOrig,LOCAL_ORIG_X,LOCAL_ORIG_Y,LOCAL_ORIG_Z);
	vectorSet(vLocalGrav,LOCAL_GRAV_X,LOCAL_GRAV_Y,LOCAL_GRAV_Z);
    vectorSet(vLocalSurf,LOCAL_SURF_X,LOCAL_SURF_Y,LOCAL_SURF_Z);
    vectorSub(vLocalOrig,d0->vAggPos,vLocalOrig);
    /* get Matrix from body frame to local frame, can be out of the loop*/
    getMatrixBDY2LOC(vLocalSurf,MatrixBDY2LOC);
    /* get initial orientation of the rubble pile, can be out of the loop*/
    getMatrixOrient(d0,MatrixOrient0);
    /* get current orientation of the rubble pile*/
    getMatrixOrient(d,MatrixOrient);
	/* calculate the current transforming matrix from mass center translating frame to body frame*/
    matrixCopy(MatrixOrient,tmpMatrix);
    matrixInvert(tmpMatrix);
    matrixMultiply(MatrixOrient0,tmpMatrix,MatrixMCT2BDY);
	/* construct MCT -> LOC transform */
	matrixMultiply(MatrixBDY2LOC,MatrixMCT2BDY,MatrixMCT2LOC);
    /* get attraction from the planet, and tranform it to local frame*/
	vectorSub(d->vPlanetPos,d->vAggPos,vPlanetPos);
    matrixTransform(MatrixMCT2BDY,vPlanetPos,tmpVector);
	vectorSub(tmpVector,vLocalOrig,tmpVector);
    getPlntGrav(tmpVector,plntGrav);
	vectorCopy(plntGrav,tmpVector);
    matrixTransform(MatrixBDY2LOC,tmpVector,plntGrav);
    /* transform local position to local frame, can be out of loop*/
	vectorCopy(vLocalOrig,tmpVector);
    matrixTransform(MatrixBDY2LOC,tmpVector,vLocalOrig);
    /* transform local gravity to local frame, can be out of loop*/
	vectorCopy(vLocalGrav,tmpVector);
    matrixTransform(MatrixBDY2LOC,tmpVector,vLocalGrav);
    /* transform vAggSpin to local frame*/
    matrixTransform(MatrixMCT2LOC,d->vAggSpin,vAggSpin);
    /* transform dAggSpinDot to local frame*/
    matrixTransform(MatrixMCT2LOC,d->vAggSpinDot,vAggSpinDot);
    /* transform dAggAcc to local frame*/
    matrixTransform(MatrixMCT2LOC,d->vAggAcc,vAggAcc);

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		vectorAdd(vLocalOrig,p->r,r);
		vectorCross(vAggSpin,r,v);
		vectorCross(vAggSpin,v,a); /* centrifugal force (per unit mass) */
		vectorSub(p->a,a,p->a);
		vectorCross(vAggSpinDot,r,a); /* torque (per unit mass) */
		vectorSub(p->a,a,p->a);
		vectorCross(vAggSpin,p->v,a);
		vectorScale(a,2.0,a); /* coriolis force (per unit mass) */
		vectorSub(p->a,a,p->a);
		vectorAdd(p->a,vLocalGrav,p->a);
		vectorAdd(p->a,plntGrav,p->a);
		vectorSub(p->a,vAggAcc,p->a);
		}
	}
    
#endif /* DEM_TIDAL_LOCAL */

#ifdef SPINUP

void pkdSpinupGetAgg(PKD pkd, double dRInner, double dROuter, int *bSuspend)
{
	PARTICLE *p;
	double r;
	int i, nLocal = pkdLocal(pkd);

	*bSuspend = 0;
	for (i = 0; i < nLocal; i++) {
		p = &pkd->pStore[i];
#ifdef SPINUP_WITH_AGGS
        p->iActualAgg = p->iOrgIdx;
#endif
		r = vectorMag(p->r); /* simple radial distance from origin */
		if (r < dRInner)
#ifdef SPINUP_WITH_AGGS
			p->iOrgIdx = -1 - INT_MAX; /* tag this particle as belonging to agg INT_MAX */
#else
			p->iOrgIdx = 0; /* tag this particle as belonging to agg -1 */
#endif /* SPINUP_WITH_AGGS */
		else if (r < dROuter)
			*bSuspend = 1;
		else
			p->iOrgIdx = 1 + p->iOrder; /* do not spinup this particle */
		}
	}

void pkdSpinupAlign(PKD pkd, Vector r_com, Vector vAxis, double dAngle)
{
	PARTICLE *p;
	Vector r;
	double sa, ca;
	int i, nLocal = pkdLocal(pkd);

	sa = sin(dAngle);
	ca = cos(dAngle);
	for (i = 0; i < nLocal; i++) {
		p = &pkd->pStore[i];
		if (p->iOrgIdx == 0) {
			vectorSub(p->r, r_com, r);
			vectorRotate2(r, vAxis, sa, ca);
			vectorAdd(r, r_com, p->r);
			}
		}
	}

void pkdSpinup(PKD pkd, Vector r_com, Vector spin)
{
	PARTICLE *p;
	Vector r, wxr;
	int i, nLocal = pkdLocal(pkd);
	
	for (i = 0; i < nLocal; i++) {
		p = &pkd->pStore[i];
		if (p->iOrgIdx == 0) {
			vectorSub(p->r, r_com, r);
			vectorCross(spin, r, wxr);
			vectorAdd(p->v, wxr, p->v);
			vectorCopy(spin, p->w);
			}
		}
	}


#ifdef SPINUP_WITH_AGGS
void pkdSpinupAccel(PKD pkd, Vector r_com, Vector omegadot, int bSpinupAlignAxes)
{
	PARTICLE *p;
	Vector r, axr;
	int i, nLocal = pkdLocal(pkd);

	for (i = 0; i < nLocal; i++) {
		p = &pkd->pStore[i];
		if (p->iOrgIdx == -1 - INT_MAX) {
			vectorSub(p->r, r_com, r);
            vectorCross(omegadot, r, axr);
            vectorAdd(p->a, axr, p->a);
            if (bSpinupAlignAxes != 0)
                vectorSub(p->r, r_com, p->r);
			}
		}
	}

void pkdSpinupResetAggs(PKD pkd)
{
    PARTICLE *p;
    int i, nLocal = pkdLocal(pkd);

    for (i = 0; i < nLocal; i++) {
        p = &pkd->pStore[i];
        p->iOrgIdx = p->iActualAgg;
        }
    }

#endif /* SPINUP_WITH_AGGS */

#endif /* SPINUP */
