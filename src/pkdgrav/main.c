#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "mdl.h"
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"

#ifdef COLLISIONS
#include <rpc/xdr.h> /* needed for time stamping terse collision log on restart */
#endif

/* 
** FPE trapping for Linux/glibc 2.2 and later.  Some compilers (e.g.,
** gcc, icc) recognize feenableexcept() without fully conforming to
** the C99 standard.  Since FPE trapping is expensive, we keep it off
** by default; set TRAP_FPE to 1 during compilation to turn it on.
** WARNING! The Intel compiler (icc) will still ignore some FPEs
** unless the "-mp" compile flag is also used!
*/
#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
/* feenableexcept() is part of the C99 standard */
# if !__APPLE__ /* FPE handing not supported under macOS */
#  define _GNU_SOURCE
#  include <fenv.h>
int feenableexcept(int);
# endif /* !__APPLE__ */
#endif

#ifdef DEM_TIDAL_SPACE
#include <unistd.h> /* for unlink() */
#include "dem.h"
void DEMTidalWriteData();
#endif /* DEM_TIDAL_SPACE */

void main_ch(MDL mdl)
{
	PST pst;
	LCL lcl;

	lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
	lcl.pkd = NULL;
	pstInitialize(&pst,mdl,&lcl);

	pstAddServices(pst,mdl);

	mdlHandler(mdl);

	pstFinish(pst);
	}

#ifdef AMPI
#define printf CmiPrintf
/* Charm MPI requires this name as "main" */
int AMPI_Main(int argc,char **argv)
#else
int main(int argc,char **argv)
#endif
{
	MDL mdl;
	MSR msr;
	FILE *fpLog = NULL;
	FILE *fpLogTiming = NULL;
	char achFile[MAXPATHLEN];
	double dTime;
	double E=0,T=0,U=0,Eth=0,L[3]={0,0,0};
	double dWMax=0,dIMax=0,dEMax=0,dMass=0,dMultiEff=0;
	long lSec=0,lStart;
	int i,iStep,iSec=0,iInterrupt,iStop=0,nActive,iNumOutputs, OutputList[NUMOUTPUTS];
	char achBaseMask[MAXPATHLEN];

#ifdef COLLISIONS
#ifdef RORY_EXTENSION
	double sec,dsec;
#endif
	int bReordered = 0;
#endif

#ifdef DEM_TIDAL_SPACE
	FILE *fpDEMTidal = NULL;
#endif

#ifdef TINY_PTHREAD_STACK
	static int first = 1;
	static char **save_argv;

	/*
	 * Hackery to get around SGI's tiny pthread stack.
	 * Main will be called twice.  The second time, argc and argv
	 * will be garbage, so we have to save them from the first.
	 * Another way to do this would involve changing the interface
	 * to mdlInitialize(), so that this hackery could be hidden
	 * down there.
	 */
	if(first) {
	    save_argv = malloc((argc+1)*sizeof(*save_argv));
	    for(i = 0; i < argc; i++)
	        save_argv[i] = strdup(argv[i]);
	    save_argv[argc] = NULL;
	    }
	else {
	    argv = save_argv;
	    }
	first = 0;
#endif /* TINY_PTHREAD_STACK */

#ifndef CCC
	/* no stdout buffering */
	setbuf(stdout,(char *) NULL);
#endif

#ifdef __FAST_MATH__
    // catch at compile time
    #error Fast math cannot be enabled
#endif

#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
	/* enable explicit FPE trapping -- see #include above */
# if __APPLE__
	#error FPE handling not supported under macOS
	/* it is possible to do, e.g., http://philbull.wordpress.com/2012/12/09/update-floating-point-exception-handling-on-mac-os-x/ -- also see README_FPE in the docs directory of this package */
# else
	feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
# endif
#endif

	lStart=time(0);
	mdlInitialize(&mdl,argv,main_ch);
	for(argc = 0; argv[argc]; argc++); /* some MDLs can trash argv */
	msrInitialize(&msr,mdl,argc,argv);

#if !(TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/)
	(void) fprintf(stderr,"WARNING: No floating point exception trapping enabled\n");
#endif

#ifdef SPRINGS
	printf("NOTE: MAX_NUM_SPRINGS_PER_PARTICLE = %i\n",MAX_NUM_SPRINGS_PER_PARTICLE);
	printf("==> .spr head size = %lu, data size = %lu per particle\n",SPRHEAD_SIZE,SPRDATA_SIZE);
#endif /* SPRINGS */

#ifdef DEM
	
	printf("NOTE: MAX_NUM_OVERLAPS_PER_PARTICLE = %i\n",MAX_NUM_OVERLAPS_PER_PARTICLE);
#ifdef WALLS
	printf("NOTE: MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS = %i\n",MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS);
#endif /* WALLS */
	printf("==> .dem head size = %lu, data size = %lu per particle\n",DEMHEAD_SIZE,DEMDATA_SIZE);

#ifdef DEM_FIXED_BALL
	printf("DEM_FIXED_BALL: Particles should be in ");
	switch (msr->param.FO.DP.iFixedBallOption) {
	case 0: /* Maximum Radius */
		printf("INCREASING");
		break;
	case 1: /* Mean */
		printf("DECREASING");
		break;
	case 2: /* Next Largest */
		printf("DECREASING");
		break;
  default:
    assert(0); /* shouldn't get here */
		}
	printf(" radius order.\n");
#else
	printf("NOTE: Particles should be in INCREASING radius order.\n");
#endif /* DEM_FIXED_BALL */
	
#endif /* DEM */

	(void) strncpy(achBaseMask,msr->param.achDigitMask,MAXPATHLEN);

	/*
	 Look for checkpoint files.  If not found, we start as normal.
	 If found, msrFindCheck() will move most recent to .chk, and 
	 we restart. bOverwrite means start from beginning, even if 
	 checkpoints exist.
	 */
	if(!msr->param.bOverwrite && msrFindCheck(msr)) {
#if defined(SPRINGS) || defined(DEM)
		msrRestartsNotSupportedMessage(msr);
#endif
		msr->param.bRestart = 1;
		dTime = msrReadCheck(msr,&iStep);
#ifdef COLLISIONS
		if (msr->param.nSmooth > msr->N) {
			msr->param.nSmooth = msr->N;
			if (msr->param.bVWarnings)
				printf("WARNING: main() restart: nSmooth reduced to %i\n",msr->N);
			}
#endif /* COLLISIONS */

#ifdef AGGS
		/*
		 ** Aggregate info not currently stored in checkpoints, so
		 ** reconstruct now.
		 */
#ifdef SLIDING_PATCH
		msr->dTime = dTime;
#endif
		msrAggsFind(msr);
#endif /* AGGS */

#ifdef WALLS
		/*
		** Positions and velocities of moving walls are not currently
		** stored in checkpoints, but positions and velocities of any
		** stuck particles do account for wall motion.  The
		** msrWallsInitPosAndVel() function needs to distinguish
		** between a restart from a checkpoint and one from an output
		** (with iStartStep > 0).  Since here we have a restart, reset
		** iStartStep just in case.
		**DEBUG! fails if restarting a run for which iStartStep > 0...
		*/
		msr->param.iStartStep = 0; /* just in case */
#endif /* WALLS */

#ifdef GASOLINE
#ifndef NOCOOLING
		if (msr->param.iGasModel == GASMODEL_COOLING
			|| msr->param.bStarForm) 
		    msrInitCooling(msr);
#endif
#endif
		if(msr->param.bRotatingBar) {
		    msrInitRotatingBar(msr, dTime);
		    }
		msrInitStep(msr);
		dMass = msrMassCheck(msr,-2.0,"Initial (Restart)");
		if (msr->param.bVStart) printf("Restart Step:%d\n",iStep);
		if (msrLogInterval(msr)) {
			snprintf(achFile,MAXPATHLEN,"%s.log",msrOutName(msr));
			fpLog = fopen(achFile,"a");
			assert(fpLog != NULL);
			setbuf(fpLog,(char *) NULL); /* no buffering */
			fprintf(fpLog,"# RESTART (dTime = %g)\n# ",dTime);
			for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
			fprintf(fpLog,"\n");
			msrLogParams(msr,fpLog);
			/* Timing data, if requested */
			if ((fpLogTiming = LogTimingInit( msr, "a" ))) {
			    fprintf(fpLogTiming,"# RESTART (dTime = %g)\n# ",dTime);
			    }
		        }
#ifdef COLLISIONS
		if (msr->param.iCollLogOption != COLL_LOG_NONE) {
			FILE *fp = fopen(msr->param.achCollLog,"r");
			if (fp) { /* add RESTART tag only if file already exists */
				fclose(fp);
				fp = fopen(msr->param.achCollLog,"a");
				assert(fp != NULL);
				switch (msr->param.iCollLogOption) {
				case COLL_LOG_VERBOSE:
					fprintf(fp,"RESTART:T=%e\n",dTime);
					break;
				case COLL_LOG_TERSE:
					{
					XDR xdrs;
					int dum = -1;
					xdrstdio_create(&xdrs,fp,XDR_ENCODE);
					(void) xdr_double(&xdrs,&dTime);
					(void) xdr_int(&xdrs,&dum);
					(void) xdr_int(&xdrs,&dum);
					(void) xdr_int(&xdrs,&dum);
					xdr_destroy(&xdrs);
					break;
					}
				default:
					assert(0); /* should never happen */
					}
				fclose(fp);
				}
			}
#endif
/* problem for restarts here, need to find second largest before following calls --RB*/
		if(msrKDK(msr) || msr->param.bGravStep || msr->param.bAccelStep) {
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
			msrInitAccel(msr);

			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
			msrUpdateSoft(msr,dTime);
			msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
			msrBuildTree(msr,0,dMass,0);
			msrMassCheck(msr,dMass,"After msrBuildTree");
			if (msrDoGravity(msr)) {
				msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
				}
			}
#ifdef GASOLINE
		msrInitSph(msr,dTime);
#endif
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr,0.0);
		/* 
		 ** Dump Frame Initialization
		 */
		/* Bring frame count up to correct place for restart. */
		if( msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 1 )
                    && msr->df->dDumpFrameStep > 0) {
			while(msr->df->dStep + msr->df->dDumpFrameStep < iStep) {
				msr->df->dStep += msr->df->dDumpFrameStep;
				msr->df->nFrame++;
			}
		}

		if (msrSteps(msr) == 0) goto CheckForDiagnosticOutput;
		goto Restart;
		}
	if(msr->param.bRestart) {
	    printf("Error: restart requested and no checkpoint file found\n");
	    msrFinish(msr);
	    mdlFinish(mdl);
	    return 1;
	    }
	    
	/*
	 ** Read in the binary file, this may set the number of timesteps or
	 ** the size of the timestep when the zto parameter is used.
	 */

#ifndef COLLISIONS
	
	dTime = msrReadTipsy(msr);

#else

	dTime = msrReadSS(msr); /* must use "Solar System" (SS) I/O format... */
	assert(dTime >= 0.0);
	if (dTime == 0.0 && msr->param.iStartStep > 0) {
		fprintf(stderr,"WARNING: time stamp in initial conditions is zero\n==> assuming dTime = iStartStep*dDelta = %g instead.\n",dTime);
		dTime = msr->param.iStartStep*msr->param.dDelta;
		}
	else if (dTime > 0.0 && msr->param.iStartStep == 0) {
		fprintf(stderr,"WARNING: time stamp in initial conditions is non-zero\n==> assuming dTime = 0.0 instead.\n");
		dTime = 0.0;
		}
	else if (dTime != msr->param.iStartStep*msr->param.dDelta) {
		/* due to possible round-off error, we can only warn here, not abort */
		fprintf(stderr,"WARNING: time stamp in initial conditions and start step may be inconsistent\nfor requested timestep (%.16e != %.16e).\n",dTime,msr->param.iStartStep*msr->param.dDelta);
		dTime = msr->param.iStartStep*msr->param.dDelta; /* force the issue! */
		fprintf(stderr,"==> time reset to %.16e.\n",dTime);
		/*assert(0);*/ /* inconsistent clocks */
		}
#ifdef DEM
	msrDEMCheckFile(msr);
#endif /* DEM */
#ifndef DEM_FIXED_BALL
	if (msr->param.nSmooth > msr->N) {
		msr->param.nSmooth = msr->N;
		if (msr->param.bVWarnings)
			printf("WARNING: main(): nSmooth reduced to %i\n",msr->N);
		}
#endif /* !DEM_FIXED_BALL */
	if (msr->param.iCollLogOption != COLL_LOG_NONE) {
		FILE *fp;
		if (msr->param.iStartStep > 0) { /* append if non-zero start step */
			fp = fopen(msr->param.achCollLog,"a");
			assert(fp != NULL);
			fprintf(fp,"START:T=%e\n",dTime);
			}
		else { /* otherwise erase any old log */
			fp = fopen(msr->param.achCollLog,"w");
			assert(fp != NULL);
			}
		fclose(fp);
		}
#ifdef RORY_EXTENSION
#ifdef SLIDING_PATCH
	if (msr->param.iRandStep) {
	    FILE *rfp = fopen("random.log","w");
	    assert(rfp);
	    fclose(rfp);
	    msr->param.iNextRandomization = msrGetNextRandomTime(msr->param.iRandStep,msr->param.iStartStep + 1);
	    }	
#endif /* SLIDING_PATCH */
#endif /* RORY_EXTENSION */

#endif /* COLLISIONS */

#ifdef GASOLINE
#ifndef NOCOOLING
	if (msr->param.iGasModel == GASMODEL_COOLING ||
		msr->param.bStarForm)
		msrInitCooling(msr);
#endif /* NOCOOLING */
#endif /* GASOLINE */

	if(msr->param.bRotatingBar)
	    msrInitRotatingBar(msr, dTime);

	msrInitStep(msr);

#ifdef GLASS
	msrInitGlass(msr);
#endif

	dMass = msrMassCheck(msr,-2.0,"Initial");

	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	msrMassCheck(msr,dMass,"After msrSetSoft");

	msrSetSink(msr);

#ifdef COLLISIONS
	if (msr->param.bFindRejects)
		msrFindRejects(msr);
#endif

#ifdef WALLS
	msrWallsInitPosAndVel(msr,dTime);
#endif

#ifdef SPRINGS
	if (msr->param.FO.iForceOverrideOption == FO_STRENGTH && !msr->param.FO.SP.bReadSpringsData)
		msrAssignSprings(msr);
#endif

#ifdef DEM
#ifdef DEM_FIXED_BALL /*This has to happen before msrAssignDEM*/
  if (msr->param.FO.DP.iFixedBallOption == 2) {
    msrDEMFixedBallNextLargest(msr);
  }
  /*Look through processors and get size of next largest particles (last particle uses itself)
   populate p->dNextLargestRadius */
#endif
	if (msr->param.FO.iForceOverrideOption == FO_STRENGTH && !msr->param.FO.DP.bReadDEMData)
		msrAssignDEM(msr);
	if (msr->param.FO.DP.iDEMStatsInterval > 0)
		msrDEMStats(msr,dTime,msr->param.iStartStep);
#endif

#ifdef AGGS
	/* find and initialize any aggregates */
#ifdef SLIDING_PATCH
	msr->dTime = dTime;
#endif
	msrAggsFind(msr);
	msrMassCheck(msr,dMass,"After msrAggsFind");
#endif /* AGGS */

#ifdef ORIGIN_HISTOGRAM
	msrInitializeOriginBins(msr);
#endif /* ORIGIN_HISTOGRAM */
	/*
	 ** If the simulation is periodic make sure to wrap all particles into
	 ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
	 */
	msrDrift(msr,dTime,0.0); /* also finds initial overlaps if applicable */
	msrMassCheck(msr,dMass,"After initial msrDrift");

 CheckForDiagnosticOutput:
	if (msrSteps(msr) > 0) {
		if (msrComove(msr)) {
			msrSwitchTheta(msr,dTime);
			}
		/*
		 ** Now we have all the parameters for the simulation we can make a 
		 ** log file entry.
		 */
		if (msrLogInterval(msr)) {
			snprintf(achFile,MAXPATHLEN,"%s.log",msrOutName(msr));
			fpLog = fopen(achFile,"w");
			assert(fpLog != NULL);
			setbuf(fpLog,(char *) NULL); /* no buffering */
			/*
			 ** Include a comment at the start of the log file showing the
			 ** command line options.
			 */
			fprintf(fpLog,"# ");
			for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
			fprintf(fpLog,"\n");
			msrLogParams(msr,fpLog);
			fprintf(fpLog,"#\n# Time Redshift TotalEnergy Kinetic Potential Thermal Lx Ly Lz GravSec WalkSecMax InteractSecMax EwaldSecMax MultistepEff\n");
			/* Timing data, if requested */
			fpLogTiming = LogTimingInit( msr, "w" );
			}

		/*
		 ** Build tree, activating all particles first (just in case).
		 */
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
		msrDomainDecomp(msr,0,1);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
		msrInitAccel(msr);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrUpdateSoft(msr,dTime);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrBuildTree(msr,0,dMass,0);
		msrMassCheck(msr,dMass,"After msrBuildTree");
		if (msrDoGravity(msr)) {
			msrGravity(msr,msr->param.iStartStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			msrMassCheck(msr,dMass,"After msrGravity");
			msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
			msrMassCheck(msr,dMass,"After msrCalcEandL");
			dMultiEff = 1.0;
			if (msrLogInterval(msr)) {
				(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
							   "%i %e %e %e %e\n",dTime,
							   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
							   E,T,U,Eth,L[0],L[1],L[2],iSec,dWMax,dIMax,dEMax,
							   dMultiEff);
				}
			/* LogTimingOutput( msr, fpLogTiming, dTime, 0 ); */
			}
#ifdef DEM_TIDAL_SPACE
		/* prepare file for recording space-frame accelerations */
		unlink(DEM_TIDAL_FILENAME); /* remove any existing file first */
		fpDEMTidal = fopen(DEM_TIDAL_FILENAME,"a");
		assert(fpDEMTidal != NULL);
		setbuf(fpDEMTidal,(char *) NULL); /* in case we forget to close it! */
		DEMTidalWriteData(msr,0.0,fpDEMTidal); /* write out data for time 0 */
#endif /* DEM_TIDAL_SPACE */

#ifdef GASOLINE
		msrInitSph(msr,dTime);
#endif
		if (msr->param.bDoSinksAtStart) msrDoSinks(msr,0.0);
		/* 
		 ** Dump Frame Initialization
		 */
		msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 0);

		LogTimingZeroCounters( msr );

#ifdef COLLISIONS
		if (msr->param.iRedOutInterval > 0) {
			/* write out a "reduced" version of the initial conditions */
			/* (includes any boundary wrapping that may have occurred) */
			int iRetVal;
			msrReorder(msr);
			iRetVal = snprintf(achFile,MAXPATHLEN,msr->param.achDigitMask,msrOutName(msr),msr->param.iStartStep);
			assert(iRetVal > 0 && iRetVal < MAXPATHLEN);
			strncat(achFile,".r",MAXPATHLEN - 1);
			msrWriteSS(msr,achFile,dTime,1/*reduced*/);
#ifdef CHARGE
			/* also output the particle charge array */
			iRetVal = snprintf(achFile,MAXPATHLEN,msr->param.achDigitMask,msrOutName(msr),msr->param.iStartStep);
			assert(iRetVal > 0 && iRetVal < MAXPATHLEN - 5);
			strncat(achFile,".charge",MAXPATHLEN - 6);
			msrOutArray(msr,achFile,OUT_CHARGE_ARRAY);
#endif /* CHARGE */
			}
#endif /* COLLISIONS */
		
		for (iStep=msr->param.iStartStep+1;iStep<=msrSteps(msr);++iStep) {
#ifdef TIMESTAMP /*DEBUG rationalize this*/
        struct tm* tm_info;
        time_t timer;
        static char buffer[100];
        time(&timer);
        tm_info = localtime(&timer);
        strftime(buffer, 25, "%Y/%m/%d %H:%M:%S", tm_info);
        fprintf(stdout, "@@ STEP %d at %s @@\n", iStep, buffer);
#endif
			if (msrComove(msr)) {
				msrSwitchTheta(msr,dTime);
				}
			if (msrKDK(msr)) {
				dMultiEff = 0.0;
				lSec = time(0);

#ifdef RORY_EXTENSION
				if (msr->param.iMinBinaryRung > 0 && 
					msr->iCurrMaxRung >= msr->param.iMinBinaryRung) {
					if (msr->param.bVDetails) {
					  sec = msrTime(msr);
					  printf("\nSearching for binaries...\n");
					  msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					  msrDomainDecomp(msr,0,1);
					  msrBuildTree(msr,0,dMass,1);
					  msrCheckForBinary(msr,dTime);
					  dsec=msrTime(msr) - sec;
					  printf("Binary search complete, Wallclock: %f sec\n\n",dsec);
					} else {
					  msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					  msrDomainDecomp(msr,0,1);
					  msrBuildTree(msr,0,dMass,1);

					  msrCheckForBinary(msr,dTime);
					}
				}
#ifdef SLIDING_PATCH
				if (msr->param.iRandStep) {
					if (iStep >= msr->param.iNextRandomization) {
					    msrRandomizeLargeMasses(msr,iStep,dTime);
					}
				}
#endif /* SLIDING_PATCH */
#endif /* RORY_EXTENSION */

#ifdef RUBBLE_ZML
				{
				  int j;
#ifdef ORIGIN_HISTOGRAM
				  int k;
#endif /* ORIGIN_HISTOGRAM */
				  msrMassCheck(msr,dMass,"Before msrRubCleanup");

				  /*
				  ** Are there any dust particles that need to be added 
				  ** to the dust bins?
				  ** Skip step 0 so initial conditions are preserved.
				  */

				  if (iStep > 0)
					  msrRubCleanup(msr,dTime);
				  msrMassCheck(msr,dMass,"After msrRubCleanup");

				  if (msr->param.CP.DB.nDustBins > 0) {
					  /*
					  ** Is it time to add dust to planetesimals?
					  */
					  if (iStep%msr->param.CP.DB.iDustBinsApplyInt == 0)
						  msrDustBinsApply(msr);
					  /* time for a dust output? */
					  if (iStep%msr->param.iOutInterval == 0) {
						  FILE *fpDust;
						  snprintf(achFile,MAXPATHLEN - 5 /* for ".dust" extension */,msr->param.achDigitMask,msrOutName(msr),iStep);
						  strcat(achFile,".dust");
						  fpDust = fopen(achFile,"w");
						  assert(fpDust != NULL);
						  for (j=0;j<msr->param.CP.DB.nDustBins;j++) {
							  fprintf(fpDust,"%i %e",j,msr->aDustBins[j].dMass);
#ifdef ORIGIN_HISTOGRAM
							  for (k=0;k<NUM_ORIGIN_BINS;k++)
								  fprintf(fpDust," %e",msr->aDustBins[j].origin_bins[k]);
#endif /* ORIGIN_HISTOGRAM */
							  fprintf(fpDust,"\n");
							  }
						  fprintf(fpDust,"Trash %e",msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
						  for (j=0;j<NUM_ORIGIN_BINS;j++)
							  fprintf(fpDust," %e",msr->DustBinsTrash.origin_bins[j]);
#endif /* ORIGIN_HISTOGRAM */
						  fprintf(fpDust,"\n");
						  (void) fclose(fpDust);
						  }
					  }
				  msrMassCheck(msr,dMass,"After dust applied to planetesimals");
				  /*
				  ** The rubble routines need to know if two
				  ** planetesimals will collide during the drift
				  ** interval so that they can be forced to the
				  ** smallest rung. But this may actually result in
				  ** the two planetesimals *not* colliding (since
				  ** their orbits will be better integrated), so
				  ** it's necessary before each top step to reset
				  ** the flags warning of imminent collision.
				  */
				  msrRubbleResetColFlag(msr);
				  msrMassCheck(msr,dMass,"Before msrTopStepKDK");
				}
#endif /* RUBBLE_ZML */

#ifdef COLLMOD_ZML /* I don't think that we need this for CollMod */
				{
				  int j;
#ifdef ORIGIN_HISTOGRAM
				  int k;
#endif /* ORIGIN_HISTOGRAM */
				  if (msr->param.CP.DB.nDustBins > 0) {
					  /*
					  ** Is it time to add dust to planetesimals?
					  */
					  if (iStep%msr->param.CP.DB.iDustBinsApplyInt == 0)
						  msrDustBinsApply(msr);
					  /* time for a dust output? */
					  if (iStep%msr->param.iOutInterval == 0) {
						  FILE *fpDust;
						  snprintf(achFile,MAXPATHLEN - 5 /* for ".dust" extension */,msr->param.achDigitMask,msrOutName(msr),iStep);
						  strcat(achFile,".dust");
						  fpDust = fopen(achFile,"w");
						  assert(fpDust != NULL);
						  for (j=0;j<msr->param.CP.DB.nDustBins;j++) {
							  fprintf(fpDust,"%i %e",j,msr->aDustBins[j].dMass);
#ifdef ORIGIN_HISTOGRAM
							  for (k=0;k<NUM_ORIGIN_BINS;k++)
								  fprintf(fpDust," %e",msr->aDustBins[j].origin_bins[k]);
#endif /* ORIGIN_HISTOGRAM */
							  fprintf(fpDust,"\n");
							  }
						  fprintf(fpDust,"Trash %e",msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
						  for (j=0;j<NUM_ORIGIN_BINS;j++)
							  fprintf(fpDust," %e",msr->DustBinsTrash.origin_bins[j]);
#endif /* ORIGIN_HISTOGRAM */
						  fprintf(fpDust,"\n");
						  (void) fclose(fpDust);
						  }
					  }
				  msrMassCheck(msr,dMass,"After dust applied to planetesimals"); /*DEBUG*/
				  /*
				  ** The collision model routines need to know if two
				  ** planetesimals will collide during the drift
				  ** interval so that they can be forced to the
				  ** smallest rung. But this may actually result in
				  ** the two planetesimals *not* colliding (since
				  ** their orbits will be better integrated), so
				  ** it's necessary before each top step to reset
				  ** the flags warning of imminent collision.
				  */
				  msrCollModResetColFlag(msr); //this can break PLANETESIMAL-PLANETESIMAL collisions
				  msrMassCheck(msr,dMass,"Before msrTopStepKDK"); /*DEBUG*/
				}
#endif /* COLLMOD_ZML */

				msrTopStepKDK(msr,iStep-1,dTime,msrDelta(msr),0,0,1,
							  &dMultiEff,&dWMax,&dIMax,&dEMax,&iSec);
				
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in KDK");
				dTime += msrDelta(msr);
				/*
				** Output a log file line if requested.
				** Note: no extra gravity calculation required.
				*/
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
				  msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
				  msrMassCheck(msr,dMass,"After msrCalcEandL in KDK");
				  lSec = time(0) - lSec;
				  (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
								 "%.16e %li %e %e %e %e\n",dTime,
								 1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
								 E,T,U,Eth,L[0],L[1],L[2],lSec,dWMax,dIMax,
								 dEMax,dMultiEff);
				}
				LogTimingOutput( msr, fpLogTiming, dTime, 0 );
			}
			else {
				lSec = time(0);
				msr->bDoneDomainDecomp = 0;
				msrTopStepDKD(msr,iStep-1,dTime,msrDelta(msr),&dMultiEff);
				msrRungStats(msr);
				msrCoolVelocity(msr,dTime,dMass); /* Supercooling if specified */
				msrMassCheck(msr,dMass,"After CoolVelocity in DKD");
				msrGrowMass(msr,dTime,msrDelta(msr)); /* Grow Masses if specified */
				dTime += msrDelta(msr);
				if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
					/*
					 ** Output a log file line.
					 ** Reactivate all particles.
					 */
					if (msrDoGravity(msr)) {
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrDomainDecomp(msr,0,1);
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrUpdateSoft(msr,dTime);
						msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
						msrBuildTree(msr,0,dMass,0);
						msrMassCheck(msr,dMass,"After msrBuildTree in DKD-log");
						msrInitAccel(msr);
						msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
						msrMassCheck(msr,dMass,"After msrGravity in DKD-log");
						}
					msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
					msrMassCheck(msr,dMass,"After msrCalcEandL in DKD-log");
					(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
								   "%.16e %li %e %e %e %e\n",dTime,
								   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
								   E,T,U,Eth,L[0],L[1],L[2],time(0)-lSec,dWMax,
								   dIMax,dEMax,dMultiEff);

					}
				LogTimingOutput( msr, fpLogTiming, dTime, 0 );
				lSec = time(0) - lSec;
				}
#ifdef DEM_TIDAL_SPACE
			if ((msr->param.iRedOutInterval > 0 && iStep%msr->param.iRedOutInterval == 0))
				DEMTidalWriteData(msr,dTime,fpDEMTidal);
#endif /* DEM_TIDAL_SPACE */
			/*
			 ** Check for user interrupt.
			 */
			iInterrupt = msrCheckForInterrupt(msr);
			switch (iInterrupt) {
			case INTR_NONE:
				break;
			case INTR_STOP:
				fprintf(stderr,"FORCED STOP: iStep=%i\n",iStep);
				iStop = 1;
				break;
			case INTR_STATUS:
				fprintf(stderr,"STATUS: iStep=%i dTime=%g N=%i\n",iStep,dTime,msr->N);
				break;
			case INTR_CHECK:
				fprintf(stderr,"FORCED CHECKPOINT: iStep=%i\n",iStep);
				break;
			case INTR_OUTPUT:
				fprintf(stderr,"FORCED OUTPUT: iStep=%i\n",iStep);
				break;
			default:
				assert(0); /* shouldn't get here */
				}
			/*
			** Output if 1) we've hit an output time
			**           2) we are stopping
			**           3) we're at an output interval
			**           4) user interrupt request
			*/
#ifndef BENCHMARK
			if (msrOutTime(msr,dTime) || iStep == msrSteps(msr) || iStop || iInterrupt == INTR_OUTPUT ||
				(msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0)) {
				if (msr->nGas && !msr->param.bKDK) {
					msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrBuildTree(msr,1,-1.0,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
					}
				msrReorder(msr);
#ifdef COLLISIONS
				bReordered = 1;
#endif
				msrMassCheck(msr,dMass,"After msrReorder in OutTime");
				snprintf(achFile,MAXPATHLEN,msr->param.achDigitMask,msrOutName(msr),iStep);
				msrCreateGasOutputList(msr, &iNumOutputs, OutputList);
				msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);

#ifdef DEM
				if (msr->param.FO.DP.nDEMOutputs > 0 || msr->param.FO.DP.nDEMOutputs < -1)
					msrDEMRemoveOldOutput(msr,iStep);
#endif /* DEM */

				if (msrDoDensity(msr) || msr->param.bDohOutput) {
					msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
					msrDomainDecomp(msr,0,1);
					msrBuildTree(msr,0,dMass,1);
					msrSmooth(msr,dTime,SMX_DENSITY,1);
				    msrReorder(msr);
					}
                msrCreateAllOutputList(msr, &iNumOutputs, OutputList);
                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);

				/*
				 ** Don't allow duplicate outputs.
				 */
				while (msrOutTime(msr,dTime));
				}
#ifdef COLLISIONS
			if ((msr->param.iRedOutInterval > 0 && iStep%msr->param.iRedOutInterval == 0)) {
				/* do a "reduced" ss output (good for movie making) */
				int iRetVal;
				msrReorder(msr);
				bReordered = 1;
				iRetVal = snprintf(achFile,MAXPATHLEN,msr->param.achDigitMask,msrOutName(msr),iStep);
				assert(iRetVal > 0 && iRetVal < MAXPATHLEN);
				strncat(achFile,".r",MAXPATHLEN - 1);
				msrWriteSS(msr,achFile,dTime,1/*reduced*/);
#ifdef CHARGE
				/* also output the particle charge array */
				iRetVal = snprintf(achFile,MAXPATHLEN,msr->param.achDigitMask,msrOutName(msr),iStep);
				assert(iRetVal > 0 && iRetVal < MAXPATHLEN - 5);
				strncat(achFile,".charge",MAXPATHLEN - 6);
				msrOutArray(msr,achFile,OUT_CHARGE_ARRAY);
#endif /* CHARGE */
				}
			if (bReordered) {
				 /*
				 ** Need to do domain decomposition now because the
				 ** smooth done as part of the collision search during
				 ** the drift assumes domain order!  The gravity step
				 ** does a decomposition, but only AFTER the drift.
				 ** Since outputs are probably infrequent (fewer than
				 ** 1 per step), it doesn't hurt much to force the
				 ** decomposition now, after every output.  This will
				 ** be inefficient if outputting every step!
				 */
				msrDomainDecomp(msr,0,1);
				bReordered = 0;
				}
#endif /* COLLISIONS */
#ifdef DEM
			if (msr->param.FO.DP.iDEMStatsInterval > 0 && iStep%msr->param.FO.DP.iDEMStatsInterval == 0)
			  msrDEMStats(msr,dTime,iStep);
#endif /* DEM */
#endif /* !BENCHMARK */
			if (!iStop && msr->param.iWallRunTime > 0) {
			    if (msr->param.iWallRunTime*60 - (time(0)-lStart) < ((int) (lSec*1.5)) ) {
					printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
					printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
						   msr->param.iWallRunTime*60,time(0)-lStart,lSec);
					iStop = 1;
					}
				}
			if (iStop || iInterrupt == INTR_CHECK || iStep == msrSteps(msr) ||
				(msrCheckInterval(msr) && iStep%msrCheckInterval(msr) == 0)) {
				/*
				 ** Write a checkpoint.
				 */
#ifndef BENCHMARK
#if !defined(SPRINGS) && !defined(DEM)
				msrWriteCheck(msr,dTime,iStep);
				msrMassCheck(msr,dMass,"After msrWriteCheck");
#endif
#endif /* !BENCHMARK */
			Restart:
				;
				}
			if (iStop) break;
			}
		if (msrLogInterval(msr)) {
		    (void) fclose(fpLog);
		    LogTimingFinish( msr, fpLogTiming, dTime );
		    }
#ifdef DEM_TIDAL_SPACE
		fclose(fpDEMTidal);
#endif /* DEM_TIDAL_SPACE */
		if (msr->param.bVStart) printf("Integration complete\n");
		}
	else {
            /* Do DiagnosticOutput */
            struct inInitDt in;
            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);

            in.dDelta = 1e37; /* large number */
            pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
            msrInitAccel(msr);
            puts("Initialized Accel and dt\n");
            snprintf(achFile,MAXPATHLEN,"%s",msrOutName(msr));

            if (msrRestart(msr)) {
                msrReorder(msr);
                snprintf(achFile,MAXPATHLEN,"%s",msrOutName(msr));
#ifndef COLLISIONS
                msrWriteTipsy(msr,achFile,dTime);
#else
                msrWriteSS(msr,achFile,dTime,0/*standard*/);
#endif
                }

            if (msrDoGravity(msr)) {
                msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
                msrDomainDecomp(msr,0,1);
                msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
                msrUpdateSoft(msr,dTime);
                msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
                msrBuildTree(msr,0,dMass,0);
                msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Gravity");
                msrGravity(msr,msr->param.iStartStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
                msrMassCheck(msr,dMass,"After msrGravity in OutSingle Gravity");
                msrReorder(msr);
                msrMassCheck(msr,dMass,"After msrReorder in OutSingle Gravity");
                iNumOutputs = 0;
                OutputList[iNumOutputs++]=OUT_ACCELG_VECTOR;
                OutputList[iNumOutputs++]=OUT_POT_ARRAY;
                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
                msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Gravity");
                }
#ifdef GASOLINE
            if (msr->nGas > 0) {
                msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                msrDomainDecomp(msr,0,1);
                msrBuildTree(msr,1,-1.0,1);
                msrSmooth(msr,dTime,SMX_DENSITY,1);
                if (msr->param.bBulkViscosity) {
                    msrReSmooth(msr,dTime,SMX_DIVVORT,1);
                    msrGetGasPressure(msr);
                    msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
                    } 
                else {
                    if (msr->param.bViscosityLimiter || msr->param.bShockTracker)
                    msrReSmooth(msr,dTime,SMX_DIVVORT,1);

                    msrSphViscosityLimiter(msr, dTime);
                    msrGetGasPressure(msr);
                    /*
                    msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
                    */
                    msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
                    msrUpdateShockTracker(msr, 0.0);
                    msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
                    /*
                    if (msr->param.bShockTracker)
                    msrReSmooth(msr,dTime,SMX_SHOCKTRACK,1);
                    */
                    /*
                    msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
                    */
                    }
                msrReorder(msr);
                if (msr->param.bSphStep) {
                    fprintf(stderr,"Adding SphStep dt\n");
                    msrSphStep(msr,dTime);
                    }
                msrCreateGasStepZeroOutputList(msr, &iNumOutputs,OutputList);
                msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
                }
#endif
            /*
             ** Build tree, activating all particles first (just in case).
             */
            if (msrDoDensity(msr) || msr->param.bDensityStep) {
                    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    msrDomainDecomp(msr,0,1);
                    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    msrBuildTree(msr,0,-1.0,1);
                    msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Density");
                    msrSmooth(msr,dTime,SMX_DENSITY,1);
                    msrMassCheck(msr,dMass,"After msrSmooth in OutSingle Density");
                    } 
            if (msrDoDensity(msr)) {
                    msrReorder(msr);
                    msrMassCheck(msr,dMass,"After msrReorder in OutSingle Density");
                    iNumOutputs = 0;
                    OutputList[iNumOutputs++]=OUT_DENSITY_ARRAY;
                    msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
            /*	snprintf(achFile,MAXPATHLEN,"%s.den",msrOutName(msr));
                    msrReorder(msr);
                    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);*/
                    msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Density");
                    }

            if (msrDoGravity(msr)) {
                    if (msr->param.bGravStep) {
                            fprintf(stderr,"Adding GravStep dt\n");
                            msrGravStep(msr,dTime);
                            }
                    if (msr->param.bAccelStep) {
                            fprintf(stderr,"Adding AccelStep dt\n");
                            msrAccelStep(msr,dTime);
                            }
                    }

            if (msr->param.bDensityStep) {
                fprintf(stderr,"Adding DensStep dt\n");
                msrDensityStep(msr,dTime);
                    }

            if (msr->param.bDeltaAccelStep) {
                fprintf(stderr,"Adding DeltaAccelStep dt\n");
                    if (!msr->param.bDeltaAccelStepGasTree) {
                        msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
                        msrBuildTree(msr,0,-1.0,1);
                        }
                    else {
                        msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
                        msrBuildTree(msr,0,-1.0,1);
                    }
                    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    msrBuildTree(msr,0,-1.0,1);
                    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
                    /* This smooth sets dt directly -- hardwired coefficient */
                    msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
                }
            /*
            msrDtToRung(msr,0,msrDelta(msr),1);
            msrRungStats(msr);
            */
            msrReorder(msr);
            iNumOutputs = 0;
            OutputList[iNumOutputs++]=OUT_DT_ARRAY;
            msrWriteOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
            /*snprintf(achFile,MAXPATHLEN,"%s.dt",msrOutName(msr));
            msrOutArray(msr,achFile,OUT_DT_ARRAY);*/
            if(msr->param.iMaxRung > 1 && (msr->param.bDensityStep || msrDoGravity(msr))) {
                msrDtToRung(msr,0,msrDelta(msr),1);
                msrRungStats(msr);
                }
		}
	
	dfFinalize( msr->df );
	msrFinish(msr);
	mdlFinish(mdl);
	return 0;
	}

#ifdef DEM_TIDAL_SPACE

void DEMTidalWriteData(MSR msr,double dTime,FILE *fp)
{
	/* output acceleration data at reduced output interval */

	DEM_TIDAL d;

	d.dTime = dTime;
	msrDEMTidalGetData(msr,&d);
	fprintf(fp,"%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
			d.dTime,
			d.vPlanetPos[0],d.vPlanetPos[1],d.vPlanetPos[2],
			d.vPlanetVel[0],d.vPlanetVel[1],d.vPlanetVel[2],
			d.vAggPos[0],d.vAggPos[1],d.vAggPos[2],
			d.vAggVel[0],d.vAggVel[1],d.vAggVel[2],
			d.vAggSpin[0],d.vAggSpin[1],d.vAggSpin[2],
			d.vAggAcc[0],d.vAggAcc[1],d.vAggAcc[2],
			d.vAggSpinDot[0],d.vAggSpinDot[1],d.vAggSpinDot[2],
			d.vMarker1Pos[0],d.vMarker1Pos[1],d.vMarker1Pos[2],
			d.vMarker1Vel[0],d.vMarker1Vel[1],d.vMarker1Vel[2],
			d.vMarker2Pos[0],d.vMarker2Pos[1],d.vMarker2Pos[2],
			d.vMarker2Vel[0],d.vMarker2Vel[1],d.vMarker2Vel[2],
			d.vMarker3Pos[0],d.vMarker3Pos[1],d.vMarker3Pos[2],
			d.vMarker3Vel[0],d.vMarker3Vel[1],d.vMarker3Vel[2]);
	}

#endif /* DEM_TIDAL_SPACE */
