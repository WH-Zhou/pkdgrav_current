#include <stdio.h>
#include <stdlib.h> /* includes malloc() macros */
#include <unistd.h> /* for unlink() */
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <sys/stat.h>

#define max(A,B) ((A) > (B) ? (A) : (B))

#include <sys/param.h> /* for MAXPATHLEN, if available */
#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

#ifdef CRAY_XT3
#include "../xdr/types.h"
#include "../xdr/xdr.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif

#include "master.h"
#include "tipsydefs.h"
#include "opentype.h"
#include "fdl.h"
#include "outtype.h"
#include "smoothfcn.h"

#ifdef AMPI
#include "charm.h"
#define printf CmiPrintf
#endif

#ifdef COLLISIONS
#include "ssdefs.h" /* in turn includes ssio.h */
#include "collision.h"
#include "linalg.h"
#endif

#ifdef AGGS
#include "matrix.h"
#endif

#ifdef RUBBLE_ZML
#include "rubble.h"
#elif defined(COLLMOD_ZML)
#include "collmod.h"
#endif

#ifdef DEM_TIDAL_SPACE
#include "dem.h"
#endif

#ifdef SLIDING_PATCH
#include "random.h"
#endif

#ifdef DEM_WALLS_REACT
#include "diffeq.h"
#endif

#define DEN_CGS_SYS 1.6831e6 /* multiply density in cgs by this to get
								density in system units (AU, M_Sun) */

void _msrLeader(void)
{
#ifdef GASOLINE
    puts("USAGE: gasoline [SETTINGS | FLAGS] [SIM_FILE]");
#else
    puts("USAGE: pkdgrav [SETTINGS | FLAGS] [SIM_FILE]");
#endif
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
	}

void _msrTrailer(void)
{
	puts("(see man page for more information)");
	}


void _msrExit(MSR msr,int status)
{
	MDL mdl=msr->mdl;
	msrFinish(msr);
	mdlFinish(mdl);
	exit(status);
	}

void _msrMakePath(const char achDir[],const char achBase[],char achPath[])
{
	/*
	** Prepends "achDir" to "achBase" and returns the result in
	** "achPath".  It is the caller's responsibility to ensure enough
	** memory has been allocated for "achPath".
	*/

	assert(achPath != NULL);
	if (achBase != NULL && achBase[0] == '/') { /* ignore achDir if already in achBase */
		strcpy(achPath,achBase);
		return;
	}
	achPath[0] = '\0';
	if (achDir != NULL) {
		strcat(achPath,achDir);
		strcat(achPath,"/");
		}
	if (achBase == NULL)
		return;
	strcat(achPath,achBase);
	}

void _msrStripQuotes(const char achIn[],char achOut[])
{
	/* Removes leading and trailing double quotes from string */

	assert(achIn != NULL && achOut != NULL);
	if (strlen(achIn) > 0 && achIn[0] == '"')
		(void) strcpy(achOut,achIn + 1);
	else
		(void) strcpy(achOut,achIn);
	if (strlen(achOut) > 0 && achOut[strlen(achOut) - 1] == '"')
		achOut[strlen(achOut) - 1] = '\0';
	}

#if defined(SPRINGS) || defined(DEM)
void msrRestartsNotSupportedMessage(MSR msr)
{
#ifdef SPRINGS
	fprintf(stderr,"ERROR: Restarts from checkpoints not supported in SPRINGS.\n");
#endif /* SPRINGS */
#ifdef DEM
	fprintf(stderr,"ERROR: Restarts from checkpoints not supported in DEM.\n");
#endif /* DEM */
	fprintf(stderr,"Instead, do the following:\n");
	fprintf(stderr,"1. Set achInFile to the latest full output.\n");
	fprintf(stderr,"2. Change iStartStep to the timestep of the output.\n");
	fprintf(stderr,"3. Make sure iCheckInterval is 0.\n");
	fprintf(stderr,"4. Remove any checkpoint files.\n");
	fprintf(stderr,"5. Rerun with the +overwrite option enabled.\n");
#ifdef SPRINGS
	fprintf(stderr,"NOTE: a corresponding .spr file is needed for spring data.\n");
	fprintf(stderr,"In that case, make sure bReadSpringsData is 1.\n");
#endif /* SPRINGS */
#ifdef DEM
	fprintf(stderr,"NOTE: a corresponding .dem file is needed to preserve contacts.\n");
	fprintf(stderr,"In that case, make sure bReadDEMData is 1.\n");
	fprintf(stderr,"IMPORTANT: DEM file format depends on certain compile options!\n");
#endif /* DEM */
#ifdef DEM_TRACK_ORIENT
	fprintf(stderr,"To read particle orientations, set bReadOrient to 1.\n");
#endif /* DEM_TRACK_ORIENT */
	_msrExit(msr,1);
	}
#endif /* SPRINGS || DEM */

void msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv)
{
	MSR msr;
	int j,ret;
	int id,nDigits;
	struct inSetAdd inAdd;
	struct inLevelize inLvl;
	struct inGetMap inGM;

#ifdef COLLISIONS
	char achUnifGravFile[MAXPATHLEN];
#endif
#ifdef WALLS
	char achWallsFile[MAXPATHLEN];
#endif
#ifdef SPECIAL_PARTICLES
	char achSpecialFile[MAXPATHLEN];
#endif
#ifdef DEM_COHESION
    char achCohesionFile[MAXPATHLEN];
#endif
    
	msr = (MSR)malloc(sizeof(struct msrContext));
	assert(msr != NULL);

	msr->bDumpFrame = 0;
	msr->df = NULL;

	msr->mdl = mdl;
	msr->pst = NULL;
	msr->lcl.pkd = NULL;
	msr->pMap = NULL;
	*pmsr = msr;
	/*
	 ** default dTimeOld value
	 */
	msr->dTimeOld = 1e-20;
	csmInitialize(&msr->param.csm);
	/*
	 ** Now setup for the input parameters.
	 **
	 ** NOTE: nThreads & bDiag are parsed here, but the actual values are
	 ** read from the command line via mdlInitialize(). This means the
	 ** values of nThreads & bDiag read by prmAddParam() are ignored!
	 */
	prmInitialize(&msr->prm,_msrLeader,_msrTrailer);
	msr->param.nThreads = 1;
	prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,sizeof(int),"sz",
				"<nThreads>");
	msr->param.bDiag = 0;
	prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,sizeof(int),"d",
				"enable/disable per thread diagnostic output");
	msr->param.bOverwrite = 0;
	prmAddParam(msr->prm,"bOverwrite",0,&msr->param.bOverwrite,sizeof(int),
				"overwrite","enable/disable checkpoint overwrite = +overwrite");
	msr->param.bVWarnings = 1;
	prmAddParam(msr->prm,"bVWarnings",0,&msr->param.bVWarnings,sizeof(int),
				"vwarnings","enable/disable warnings = +vwarnings");
	msr->param.bVStart = 1;
	prmAddParam(msr->prm,"bVStart",0,&msr->param.bVStart,sizeof(int),
				"vstart","enable/disable verbose start = +vstart");
	msr->param.bVStep = 1;
	prmAddParam(msr->prm,"bVStep",0,&msr->param.bVStep,sizeof(int),
				"vstep","enable/disable verbose step = +vstep");
	msr->param.bVRungStat = 1;
	prmAddParam(msr->prm,"bVRungStat",0,&msr->param.bVRungStat,sizeof(int),
				"vrungstat","enable/disable rung statistics = +vrungstat");
	msr->param.bVDetails = 0;
	prmAddParam(msr->prm,"bVDetails",0,&msr->param.bVDetails,sizeof(int),
				"vdetails","enable/disable verbose details = -vdetails");
	msr->param.bLogTiming = 0;
	prmAddParam(msr->prm,"bLogTiming",0,&msr->param.bLogTiming,sizeof(int),
				"logtiming","enable/disable log of timing data = -timing");
	msr->param.bLogTimingStep = 0;
	prmAddParam(msr->prm,"bLogTimingStep",0,&msr->param.bLogTimingStep,sizeof(int),
				"logtimingstep","log of timing data per step = +timings");
	msr->param.bLogTimingSubStep = 0;
	prmAddParam(msr->prm,"bLogTimingSubStep",0,&msr->param.bLogTimingSubStep,sizeof(int),
				"logtimingstep","log of timing data per sub step = +timingss");
	msr->param.bLogTimingStepTot = 0;
	prmAddParam(msr->prm,"bLogTimingStepTot",0,&msr->param.bLogTimingStepTot,sizeof(int),
				"logtimingstep","log of total timing data per step = +timingst");
	msr->param.bLogTimingSubStepTot = 0;
	prmAddParam(msr->prm,"bLogTimingSubStepTot",0,&msr->param.bLogTimingSubStepTot,sizeof(int),
				"logtimingstep","log of total timing data per substep = +timingsst");
	nDigits = 5;
	prmAddParam(msr->prm,"nDigits",1,&nDigits,sizeof(int),"nd",
				"<number of digits to use in output filenames> = 5");
	msr->param.bPeriodic = 0;
	prmAddParam(msr->prm,"bPeriodic",0,&msr->param.bPeriodic,sizeof(int),"p",
				"periodic/non-periodic = -p");
	msr->param.bRestart = 0;
	prmAddParam(msr->prm,"bRestart",0,&msr->param.bRestart,sizeof(int),"restart",
				"restart from checkpoint");
	msr->param.bParaRead = 1;
	prmAddParam(msr->prm,"bParaRead",0,&msr->param.bParaRead,sizeof(int),"par",
				"enable/disable parallel reading of files = +par");
	msr->param.bParaWrite = 1;
	prmAddParam(msr->prm,"bParaWrite",0,&msr->param.bParaWrite,sizeof(int),"paw",
				"enable/disable parallel writing of files = +paw");
	msr->param.bUseWallClock = 1;
	prmAddParam(msr->prm,"bUseWallClock",0,&msr->param.bUseWallClock,sizeof(int),"wallclock",
				"enable/disable use of wall clock for timing = +wallclock");
	msr->param.bCannonical = 1;
	prmAddParam(msr->prm,"bCannonical",0,&msr->param.bCannonical,sizeof(int),"can",
				"enable/disable use of cannonical momentum = +can");
	msr->param.bKDK = 1;
	prmAddParam(msr->prm,"bKDK",0,&msr->param.bKDK,sizeof(int),"kdk",
				"enable/disable use of kick-drift-kick integration = +kdk");
	msr->param.bBinary = 1;
	prmAddParam(msr->prm,"bBinary",0,&msr->param.bBinary,sizeof(int),"bb",
				"spatial/density binary trees = +bb");
	msr->param.iBinaryOutput = 0;
	prmAddParam(msr->prm,"iBinaryOutput",1,&msr->param.iBinaryOutput,sizeof(int),
				"binout","<array outputs 0 ascii, 1 float, 2 double, 3 FLOAT(internal)> = 0");
	msr->param.bPackedVector = 0;
	prmAddParam(msr->prm,"bPackedVector",0,&msr->param.bPackedVector,sizeof(int),
				"pvec","enable/disable packed vector outputs = +pvec");
	msr->param.bDoDensity = 1;
	prmAddParam(msr->prm,"bDoDensity",0,&msr->param.bDoDensity,sizeof(int),
				"den","enable/disable density outputs = +den");
	msr->param.iReadIOrder = 0;
	prmAddParam(msr->prm,"iReadIOrder",1,&msr->param.iReadIOrder,sizeof(int),
				"iordin","<array outputs 0 NO, 1 int, 2 long, 3 int (internal)> = 0");
	msr->param.bDoIOrderOutput = 0;
	prmAddParam(msr->prm,"bDoIOrderOutput",0,&msr->param.bDoIOrderOutput,
		    sizeof(int), "iordout","enable/disable iOrder outputs = -iordout");
	msr->param.bDohOutput = 0;
	prmAddParam(msr->prm,"bDohOutput",0,&msr->param.bDohOutput,sizeof(int),
				"hout","enable/disable h outputs = -hout");
	msr->param.bDoSphhOutput = 0;
	prmAddParam(msr->prm,"bDoSphhOutput",0,&msr->param.bDoSphhOutput,sizeof(int),
				"sphhout","enable/disable Sph h outputs = -sphhout");
	msr->param.bDodtOutput = 0;
	prmAddParam(msr->prm,"bDodtOutput",0,&msr->param.bDodtOutput,sizeof(int),
				"dtout","enable/disable dt outputs = -dtout");
	msr->param.nBucket = 8;
	prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,sizeof(int),"b",
				"<max number of particles in a bucket> = 8");
	msr->param.iStartStep = 0;
	prmAddParam(msr->prm,"iStartStep",1,&msr->param.iStartStep,
				sizeof(int),"nstart","<initial step numbering> = 0");
	msr->param.nSteps = 0;
	prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,sizeof(int),"n",
				"<number of timesteps> = 0");
	msr->param.iOutInterval = 0;
	prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,sizeof(int),
				"oi","<number of timesteps between snapshots> = 0");
#ifdef COLLISIONS /* "reduced" ss output (good for movie making) */
	msr->param.iRedOutInterval = 0;
	prmAddParam(msr->prm,"iRedOutInterval",1,&msr->param.iRedOutInterval,sizeof(int),
				"oir","<number of timesteps between reduced snapshots> = 0");
#endif
	msr->param.dDumpFrameStep = 0;
	prmAddParam(msr->prm,"dDumpFrameStep",2,&msr->param.dDumpFrameStep,sizeof(double),
				"dfi","<number of steps between dumped frames> = 0");
	msr->param.dDumpFrameTime = 0;
	prmAddParam(msr->prm,"dDumpFrameTime",2,&msr->param.dDumpFrameTime,sizeof(double),
				"dft","<number of timesteps between dumped frames> = 0");
	msr->param.iLogInterval = 10;
	prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,sizeof(int),
				"ol","<number of timesteps between logfile outputs> = 10");
	msr->param.iCheckInterval = 10;
	prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,sizeof(int),
				"oc","<number of timesteps between checkpoints> = 10");
	msr->param.iOrder = 4;
	prmAddParam(msr->prm,"iOrder",1,&msr->param.iOrder,sizeof(int),"or",
				"<multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.bEwald = 1;
	prmAddParam(msr->prm,"bEwald",0,&msr->param.bEwald,sizeof(int),"ewald",
				"enable/disable Ewald correction = +ewald");
	msr->param.iEwOrder = 4;
	prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,sizeof(int),"ewo",
				"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
	msr->param.nReplicas = 0;
	prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,sizeof(int),"nrep",
				"<nReplicas> = 0 for -p, or 1 for +p");
	msr->param.dSoft = 0.0;
	prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,sizeof(double),"e",
				"<gravitational softening length> = 0.0");
	msr->param.dSoftMax = 0.0;
	prmAddParam(msr->prm,"dSoftMax",2,&msr->param.dSoftMax,sizeof(double),"eMax",
				"<maximum comoving gravitational softening length (abs or multiplier)> = 0.0");
	msr->param.bPhysicalSoft = 0;
	prmAddParam(msr->prm,"bPhysicalSoft",0,&msr->param.bPhysicalSoft,sizeof(int),"PhysSoft",
				"<Physical gravitational softening length> -PhysSoft");
	msr->param.bSoftMaxMul = 1;
	prmAddParam(msr->prm,"bSoftMaxMul",0,&msr->param.bSoftMaxMul,sizeof(int),"SMM",
				"<Use maximum comoving gravitational softening length as a multiplier> +SMM");
	msr->param.bVariableSoft = 0;
	prmAddParam(msr->prm,"bVariableSoft",0,&msr->param.bVariableSoft,sizeof(int),"VarSoft",
				"<Variable gravitational softening length> -VarSoft");
	msr->param.bVariableSoftStar = 1;
	prmAddParam(msr->prm,"bVariableSoftStar",0,&msr->param.bVariableSoftStar,sizeof(int),"VarSoftStar",
				"<Variable gravitational softening length stars> -VarSoftStar");
	msr->param.bVariableSoftGas = 1;
	prmAddParam(msr->prm,"bVariableSoftGas",0,&msr->param.bVariableSoftGas,sizeof(int),"VarSoftGas",
				"<Variable gravitational softening length gas> -VarSoftGas");
	msr->param.bVariableSoftDark = 1;
	prmAddParam(msr->prm,"bVariableSoftDark",0,&msr->param.bVariableSoftDark,sizeof(int),"VarSoft",
				"<Variable gravitational softening length dark> -VarSoftDark");
	msr->param.nSoftNbr = 32;
	prmAddParam(msr->prm,"nSoftNbr",1,&msr->param.nSoftNbr,sizeof(int),"VarSoft",
				"<Neighbours for Variable gravitational softening length> 32");
	msr->param.bSoftByType = 1;
	prmAddParam(msr->prm,"bSoftByType",0,&msr->param.bSoftByType,sizeof(int),"SBT",
				"<Variable gravitational softening length by Type> +SBT");
	msr->param.bDoSoftOutput = 0;
	prmAddParam(msr->prm,"bDoSoftOutput",0,&msr->param.bDoSoftOutput,sizeof(int),
				"softout","enable/disable soft outputs = -softout");

	msr->param.dDelta = 0.0;
	prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,sizeof(double),"dt",
				"<time step>");
	msr->param.dEta = 0.1;
	prmAddParam(msr->prm,"dEta",2,&msr->param.dEta,sizeof(double),"eta",
				"<time step criterion> = 0.1");
	msr->param.dEtaDeltaAccel = 0.2;
	prmAddParam(msr->prm,"dEtaDeltaAccel",2,&msr->param.dEtaDeltaAccel,sizeof(double),"etadrda",
				"<drda time step criterion> = 0.2");
	msr->param.dEtaCourant = 0.4;
	prmAddParam(msr->prm,"dEtaCourant",2,&msr->param.dEtaCourant,sizeof(double),"etaC",
				"<Courant criterion> = 0.4");
	msr->param.dEtauDot = 0.25;
	prmAddParam(msr->prm,"dEtauDot",2,&msr->param.dEtauDot,sizeof(double),"etau",
				"<uDot criterion> = 0.25");
	msr->param.duDotLimit = -0.2;
	prmAddParam(msr->prm,"duDotLimit",2,&msr->param.duDotLimit,sizeof(double),"uDL",
				"<uDotLimit:  Treat udot/u < duDotLimit specially> = -0.2 < 0");
	msr->param.bGravStep = 0;
	prmAddParam(msr->prm,"bGravStep",0,&msr->param.bGravStep,sizeof(int),
				"gs","<Symmetric gravity timestepping (sqrt(r^3/(mi + mj)))>");
	msr->param.bEpsAccStep = 1;
	prmAddParam(msr->prm,"bEpsAccStep",0,&msr->param.bEpsAccStep,sizeof(int),
				"ea", "<Sqrt(Epsilon on a) timestepping>");
	msr->param.bSqrtPhiStep = 0;
	prmAddParam(msr->prm,"bSqrtPhiStep",0,&msr->param.bSqrtPhiStep,sizeof(int),
				"sphi", "<Sqrt(Phi) on a timestepping>");
	msr->param.bDensityStep = 0;
	prmAddParam(msr->prm,"bDensityStep",0,&msr->param.bDensityStep,sizeof(int),
				"isrho", "<Sqrt(1/Rho) timestepping>");
	msr->param.bDeltaAccelStep = 0;
	prmAddParam(msr->prm,"bDeltaAccelStep",0,&msr->param.bDeltaAccelStep,sizeof(int),
				"isdrda", "<Sqrt(dr/da) timestepping>");
	msr->param.bDeltaAccelStepGasTree = 0;
	prmAddParam(msr->prm,"bDeltaAccelStepGasTree",0,&msr->param.bDeltaAccelStepGasTree,sizeof(int),
				"isdrdagt", "<Sqrt(dr/da) timestepping via gas tree>");
	msr->param.nTruncateRung = 0;
	prmAddParam(msr->prm,"nTruncateRung",1,&msr->param.nTruncateRung,sizeof(int),"nTR",
				"<number of MaxRung particles to delete MaxRung> = 0");
	msr->param.bNonSymp = 1;
	prmAddParam(msr->prm,"bNonSymp",0,&msr->param.bNonSymp,sizeof(int),
				"ns", "<Non-symplectic density stepping>");
	msr->param.iMaxRung = 1;
	prmAddParam(msr->prm,"iMaxRung",1,&msr->param.iMaxRung,sizeof(int),
				"mrung", "<maximum timestep rung>");
#ifdef COLLISIONS
	/* Set to < 0 to disable */
	msr->param.iMinBinaryRung=-1;
	prmAddParam(msr->prm,"iMinBinaryRung",1,&msr->param.iMinBinaryRung,
				sizeof(int),"MinBinaryRung",
				"<minimum rung to search for binaries>");
	msr->param.dBallVelFact=0;
	prmAddParam(msr->prm,"dBallVelFact",2,&msr->param.dBallVelFact,
		    sizeof(double),"BallVelFact",
		    "<radius of collision search: r=dBVF*dDelta*velocity>");
	msr->param.dMaxBinaryEcc=1;
	prmAddParam(msr->prm,"dMaxBinaryEcc",2,&msr->param.dMaxBinaryEcc,
		    sizeof(double),"MaxBinaryEcc",
		    "<maximum binary eccentricity>");
	assert(msr->param.dMaxBinaryEcc <= 1 && msr->param.dMaxBinaryEcc >= 0);
	msr->param.iMinCollRung=0;
	prmAddParam(msr->prm,"iMinCollRung",1,&msr->param.iMinCollRung,
		    sizeof(int),"MinCollisionRung",
		    "<minimum rung to search for collisions");
	if (msr->param.iMinCollRung)
	  assert(msr->param.iMaxRung > 1);
#ifdef SLIDING_PATCH
	msr->param.iRandStep=0;
	prmAddParam(msr->prm,"iRandStep",1,&msr->param.iRandStep,
		    sizeof(int),"StepsForRandomizingLargeMasses",
		    "<whole steps between randomizing big masses>");
	msr->param.dLargeMass=0;
	prmAddParam(msr->prm,"dLargeMass",2,&msr->param.dLargeMass,
		    sizeof(double),"MinMassForRandomization",
		    "<minimum mass to randomly move>");
	if (msr->param.iRandStep)
	    assert(msr->param.dLargeMass > 0);
	msr->param.dRandBall=0;
	prmAddParam(msr->prm,"dRandBall",2,&msr->param.dRandBall,
		    sizeof(double),"SizeOfBallDuringRandomMove",
		    "<size of sphere to rondomly move>");
	if (msr->param.dRandBall > 0) {
	    assert(msr->param.dLargeMass != 0);
	    assert(msr->param.iRandStep);
	    }
#endif /* SLIDING_PATCH */
#endif /* COLLISIONS */
	msr->param.dEwCut = 2.6;
	prmAddParam(msr->prm,"dEwCut",2,&msr->param.dEwCut,sizeof(double),"ew",
				"<dEwCut> = 2.6");
	msr->param.dEwhCut = 2.8;
	prmAddParam(msr->prm,"dEwhCut",2,&msr->param.dEwhCut,sizeof(double),"ewh",
				"<dEwhCut> = 2.8");
	msr->param.dTheta = 0.8;
	msr->param.dTheta2 = msr->param.dTheta;
	prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,sizeof(double),"theta",
				"<Barnes opening criterion> = 0.8");
	prmAddParam(msr->prm,"dTheta2",2,&msr->param.dTheta2,sizeof(double),
				"theta2","<Barnes opening criterion after a < daSwitchTheta> = 0.8");
	msr->param.daSwitchTheta = 1./3.;
	prmAddParam(msr->prm,"daSwitchTheta",2,&msr->param.daSwitchTheta,sizeof(double),"aSwitchTheta",
				"<a to switch theta at> = 1./3.");
	msr->param.dAbsPartial = 0.0;
	prmAddParam(msr->prm,"dAbsPartial",2,&msr->param.dAbsPartial,sizeof(double),"ap",
				"<absolute partial error opening criterion>");
	msr->param.dRelPartial = 0.0;
	prmAddParam(msr->prm,"dRelPartial",2,&msr->param.dRelPartial,sizeof(double),"rp",
				"<relative partial error opening criterion>");
	msr->param.dAbsTotal = 0.0;
	prmAddParam(msr->prm,"dAbsTotal",2,&msr->param.dAbsTotal,sizeof(double),"at",
				"<absolute total error opening criterion>");
	msr->param.dRelTotal = 0.0;
	prmAddParam(msr->prm,"dRelTotal",2,&msr->param.dRelTotal,sizeof(double),"rt",
				"<relative total error opening criterion>");

	msr->param.bDoSinks = 0;
	prmAddParam(msr->prm,"bDoSinks",0,&msr->param.bDoSinks,sizeof(int),
				"sinks","enable/disable sinks = -sinks");
	msr->param.bBHSink = 0;
	prmAddParam(msr->prm,"bBHSink",0,&msr->param.bBHSink,sizeof(int),
				"bhsink","Bondi-Hoyle type sink = -bhsink");
	msr->param.dBHSinkEddEff = 0.1;
	prmAddParam(msr->prm,"dBHSinkEddEff",2,&msr->param.dBHSinkEddEff,sizeof(double),"bhsinkeddeff",
				"<BHSink Eddington Efficiency>");
	msr->param.dBHSinkFeedbackEff = 0.05;
	prmAddParam(msr->prm,"dBHSinkFeedbackEff",2,&msr->param.dBHSinkFeedbackEff,sizeof(double),"bhsinkfbeff",
				"<BHSink Feedback Efficiency>");
	msr->param.dBHSinkAlpha = 1;
	prmAddParam(msr->prm,"dBHSinkAlpha",2,&msr->param.dBHSinkAlpha,sizeof(double),"bhsinkalpha",
				"<BHSink Alpha>");
	msr->param.bDoSinksAtStart = 0;
	prmAddParam(msr->prm,"bDoSinksAtStart",0,&msr->param.bDoSinksAtStart,sizeof(int),
				"sinksas","enable/disable sinks at start = -sinksas");
	msr->param.bSinkThermal = 0;
	prmAddParam(msr->prm,"bSinkThermal",0,&msr->param.bSinkThermal,sizeof(int),
				"tsinks","enable/disable thermal energy in sink calcs = -tsinks");
	msr->param.dSinkRadius = 0.0;
	prmAddParam(msr->prm,"dSinkRadius",2,&msr->param.dSinkRadius,sizeof(double),"sinkr",
				"<Sink Radius>");
	msr->param.dSinkBoundOrbitRadius = 0.0;
	prmAddParam(msr->prm,"dSinkBoundOrbitRadius",2,&msr->param.dSinkBoundOrbitRadius,sizeof(double),"sinkbor",
				"<Sink Bound Orbit Radius>");
	msr->param.dDeltaSink = msr->param.dDelta;
	prmAddParam(msr->prm,"dDeltaSink", 2, &msr->param.dDeltaSink,
		    sizeof(double), "dDeltaSink",
		    "<Minimum sink timestep in years> = dDelta");
	msr->param.dSinkMassMin = 0;  /* Default reset to FLT_MAX for BH sink */
	prmAddParam(msr->prm,"dSinkMassMin", 2, &msr->param.dSinkMassMin,
		    sizeof(double), "dSinkMassMin", "<Minimum Mass to act as a sink> = 0" );
	msr->param.iSinkRung = 0; 
	prmAddParam(msr->prm,"iSinkRung", 2, &msr->param.iSinkRung,
		    sizeof(double), "iSinkRung",
		    "<Sink Rung> = 0");
	msr->param.dPeriod = 1.0;
	prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,sizeof(double),"L",
				"<periodic box length> = 1.0");
	msr->param.dxPeriod = 1.0;
	prmAddParam(msr->prm,"dxPeriod",2,&msr->param.dxPeriod,sizeof(double),"Lx",
				"<periodic box length in x-dimension> = 1.0");
	msr->param.dyPeriod = 1.0;
	prmAddParam(msr->prm,"dyPeriod",2,&msr->param.dyPeriod,sizeof(double),"Ly",
				"<periodic box length in y-dimension> = 1.0");
	msr->param.dzPeriod = 1.0;
	prmAddParam(msr->prm,"dzPeriod",2,&msr->param.dzPeriod,sizeof(double),"Lz",
				"<periodic box length in z-dimension> = 1.0");
	msr->param.achInFile[0] = '\0';
	prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,256,"I",
				"<input file name> (file in TIPSY binary format)");
#ifdef GASOLINE
	strcpy(msr->param.achOutName,"gasoline");
	prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"gasoline\"");
#else
	strcpy(msr->param.achOutName,"pkdgrav");
	prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
				"<output name for snapshots and logfile> = \"pkdgrav\"");
#endif
	msr->param.csm->bComove = 0;
	prmAddParam(msr->prm,"bComove",0,&msr->param.csm->bComove,sizeof(int),
				"cm", "enable/disable comoving coordinates = -cm");
	msr->param.csm->dHubble0 = 0.0;
	prmAddParam(msr->prm,"dHubble0",2,&msr->param.csm->dHubble0, 
				sizeof(double),"Hub", "<dHubble0> = 0.0");
	msr->param.csm->dOmega0 = 1.0;
	prmAddParam(msr->prm,"dOmega0",2,&msr->param.csm->dOmega0,
				sizeof(double),"Om", "<dOmega0> = 1.0");
	msr->param.csm->dLambda = 0.0;
	prmAddParam(msr->prm,"dLambda",2,&msr->param.csm->dLambda,
				sizeof(double),"Lambda", "<dLambda> = 0.0");
	msr->param.csm->dOmegaRad = 0.0;
	prmAddParam(msr->prm,"dOmegaRad",2,&msr->param.csm->dOmegaRad,
				sizeof(double),"Omrad", "<dOmegaRad> = 0.0");
	msr->param.csm->dOmegab = 0.0;
	prmAddParam(msr->prm,"dOmegab",2,&msr->param.csm->dOmegab,
				sizeof(double),"Omb", "<dOmegab> = 0.0");
	msr->param.csm->dQuintess = 0.0;
	prmAddParam(msr->prm,"dQuintess",2,&msr->param.csm->dQuintess,
				sizeof(double),"Quint",
		    "<dQuintessence (constant w = -1/2) > = 0.0");
	strcpy(msr->param.achDataSubPath,".");
	prmAddParam(msr->prm,"achDataSubPath",3,msr->param.achDataSubPath,256,
				NULL,NULL);
	msr->param.dExtraStore = 0.1;
	prmAddParam(msr->prm,"dExtraStore",2,&msr->param.dExtraStore,
				sizeof(double),NULL,NULL);
	msr->param.nSmooth = 64;
	prmAddParam(msr->prm,"nSmooth",1,&msr->param.nSmooth,sizeof(int),"s",
				"<number of particles to smooth over> = 64");
	msr->param.bStandard = 0;
	prmAddParam(msr->prm,"bStandard",0,&msr->param.bStandard,sizeof(int),"std",
				"output in standard TIPSY binary format = -std");
	msr->param.dRedTo = 0.0;	
	prmAddParam(msr->prm,"dRedTo",2,&msr->param.dRedTo,sizeof(double),"zto",
				"specifies final redshift for the simulation");
	msr->param.nSuperCool = 0;
	prmAddParam(msr->prm,"nSuperCool",1,&msr->param.nSuperCool,sizeof(int),
				"scn","<number of supercooled particles> = 0");
	msr->param.dCoolFac = 0.95;
	prmAddParam(msr->prm,"dCoolFac",2,&msr->param.dCoolFac,sizeof(double),
				"scf","<Velocity Cooling factor> = 0.95 (no cool = 1.0)");
	msr->param.dCoolDens = 50.0;
	prmAddParam(msr->prm,"dCoolDens",2,&msr->param.dCoolDens,sizeof(double),
				"scd","<Velocity Cooling Critical Density> = 50");
	msr->param.dCoolMaxDens = 1e8;
	prmAddParam(msr->prm,"dCoolMaxDens",2,&msr->param.dCoolMaxDens,
				sizeof(double),"scmd",
				"<Velocity Cooling Maximum Density> = 1e8");
	msr->param.bSymCool = 0;
	prmAddParam(msr->prm,"bSymCool",0,&msr->param.bSymCool,sizeof(int),NULL,
				NULL);
	msr->param.nGrowMass = 0;
	prmAddParam(msr->prm,"nGrowMass",1,&msr->param.nGrowMass,sizeof(int),
				"gmn","<number of particles to increase mass> = 0");
	msr->param.dGrowDeltaM = 0.0;
	prmAddParam(msr->prm,"dGrowDeltaM",2,&msr->param.dGrowDeltaM,
				sizeof(double),"gmdm","<Total growth in mass/particle> = 0.0");
	msr->param.dGrowStartT = 0.0;
	prmAddParam(msr->prm,"dGrowStartT",2,&msr->param.dGrowStartT,
				sizeof(double),"gmst","<Start time for growing mass> = 0.0");
	msr->param.dGrowEndT = 1.0;
	prmAddParam(msr->prm,"dGrowEndT",2,&msr->param.dGrowEndT,
				sizeof(double),"gmet","<End time for growing mass> = 1.0");
	msr->param.dFracNoDomainDecomp = 0.002;
	prmAddParam(msr->prm,"dFracNoDomainDecomp",2,&msr->param.dFracNoDomainDecomp,
				sizeof(double),"fndd",
				"<Fraction of Active Particles for no new DD> = 0.002");
	msr->param.dFracNoDomainDimChoice = 0.1;
	prmAddParam(msr->prm,"dFracNoDomainDimChoice",2,&msr->param.dFracNoDomainDimChoice,
				sizeof(double),"fnddc",
				"<Fraction of Active Particles for no new DD dimension choice> = 0.1");
	msr->param.dFracFastGas = 0.2;
	prmAddParam(msr->prm,"dFracFastGas",2,&msr->param.dFracFastGas,
				sizeof(double),"fndd",
				"<Fraction of Active Particles for Fast Gas> = 0.01");
	msr->param.dhMinOverSoft = 0.0;
	prmAddParam(msr->prm,"dhMinOverSoft",2,&msr->param.dhMinOverSoft,
				sizeof(double),"hmin",
				"<Minimum h as a fraction of Softening> = 0.0");
	msr->param.bDoGravity = 1;
	prmAddParam(msr->prm,"bDoGravity",0,&msr->param.bDoGravity,sizeof(int),"g",
				"enable/disable gravity (interparticle and external potentials) = +g");
	msr->param.bDoSelfGravity = 1;
	prmAddParam(msr->prm,"bDoSelfGravity",0,&msr->param.bDoSelfGravity,sizeof(int),"sg",
				"enable/disable interparticle self gravity = +sg");
	msr->param.bRungDD = 0;
	prmAddParam(msr->prm,"bRungDomainDecomp",0,&msr->param.bRungDD,sizeof(int),
				"RungDD","<Rung Domain Decomp> = 0");
	msr->param.dRungDDWeight = 1.0;
	prmAddParam(msr->prm,"dRungDDWeight",2,&msr->param.dRungDDWeight,sizeof(int),
				"RungDDWeight","<Rung Domain Decomp Weight> = 1.0");
	msr->param.bFandG = 0;
	prmAddParam(msr->prm,"bFandG",0,&msr->param.bFandG,sizeof(int),"fg",
				"use/don't use Kepler orbit drifts = -fg");
	msr->param.bHeliocentric = 0;
	prmAddParam(msr->prm,"bHeliocentric",0,&msr->param.bHeliocentric,
				sizeof(int),"hc","use/don't use Heliocentric coordinates = -hc");
	msr->param.dCentMass = 0.0;
	prmAddParam(msr->prm,"dCentMass",2,&msr->param.dCentMass,sizeof(double),
				"fgm","specifies the central mass for Keplerian orbits");
	msr->param.bLogHalo = 0;
	prmAddParam(msr->prm,"bLogHalo",0,&msr->param.bLogHalo,
				sizeof(int),"halo","use/don't use galaxy halo = -halo");
	msr->param.bHernquistSpheroid = 0;
	prmAddParam(msr->prm,"bHernquistSpheroid",0,&msr->param.bHernquistSpheroid,
				sizeof(int),"hspher","use/don't use galaxy Hernquist Spheroid = -hspher");
	msr->param.bNFWSpheroid = 0;
	prmAddParam(msr->prm,"bNFWSpheroid",0,&msr->param.bNFWSpheroid,
				sizeof(int),"NFWspher","use/don't use galaxy NFW Spheroid = -NFWspher");
        prmAddParam(msr->prm,"dNFWm200",2,&msr->param.dNFWm200,sizeof(double),
                    "dNFWm200","Mass inside rho/rho_c = 200");
        prmAddParam(msr->prm,"dNFWr200",2,&msr->param.dNFWr200,sizeof(double),
                    "dNFWr200","Radius of rho/rho_c = 200");
        prmAddParam(msr->prm,"dNFWconc",2,&msr->param.dNFWconc,sizeof(double),
                    "dNFWconc","NFW concentration");
        prmAddParam(msr->prm,"dNFWsoft",2,&msr->param.dNFWsoft,sizeof(double),
                    "dNFWsoft","Fixed potential softening length");
	msr->param.bElliptical=0;
	prmAddParam(msr->prm,"bElliptical",0,&msr->param.bElliptical,
				sizeof(int),"elliptical","use/don't");
	msr->param.bEllipticalDarkNFW=0;
	prmAddParam(msr->prm,"bEllipticalDarkNFW",0,&msr->param.bEllipticalDarkNFW,
		    sizeof(int),"ellipticaldarknfw","use/dont");
	msr->param.bHomogSpheroid = 0;
	prmAddParam(msr->prm,"bHomogSpheroid",0,&msr->param.bHomogSpheroid,
				sizeof(int),"hspher","use/don't use galaxy Homog Spheroid = -homogspher");
	msr->param.bBodyForce = 0;
	prmAddParam(msr->prm,"bBodyForce",0,&msr->param.bBodyForce,
				sizeof(int),"bodyforce","use/don't use body force = -bf");
	msr->param.bMiyamotoDisk = 0;
	prmAddParam(msr->prm,"bMiyamotoDisk",0,&msr->param.bMiyamotoDisk,
				sizeof(int),"mdisk","use/don't use galaxy Miyamoto Disk = -mdisk");
	msr->param.bTimeVarying = 0;
	prmAddParam(msr->prm,"bTimeVarying",0,&msr->param.bTimeVarying,
				sizeof(int),"tvar","use/don't use the time varying externalpotential = -tvar");
	rotbarInitialize(&msr->param.rotbar);
	msr->param.bRotatingBar = 0;
	prmAddParam(msr->prm,"bRotatingBar",0,&msr->param.bRotatingBar,
		    sizeof(int),"rotbar",
		    "use/don't use rotating bar = -rotbar");
	rotbarAddParams(msr->param.rotbar, msr->prm);
#ifdef ROT_FRAME
	msr->param.bRotFrame = 0;
	prmAddParam(msr->prm,"bRotFrame",0,&msr->param.bRotFrame,
				sizeof(int),"rframe","use/don't use rotating frame = -rframe");
	msr->param.dOmega = 0;
	prmAddParam(msr->prm,"dOmega",2,&msr->param.dOmega,sizeof(double),
				"omega","rotating frame <dOmega> = 0");
	msr->param.dOmegaDot = 0;
	prmAddParam(msr->prm,"dOmegaDot",2,&msr->param.dOmegaDot,sizeof(double),
				"omegadot","<dOmegaDot> = 0");
#endif
	msr->param.iWallRunTime = 0;
	prmAddParam(msr->prm,"iWallRunTime",1,&msr->param.iWallRunTime,
				sizeof(int),"wall",
				"<Maximum Wallclock time (in minutes) to run> = 0 = infinite");
	msr->param.dSunSoft = 0.0;
	prmAddParam(msr->prm,"dSunSoft", 2, &msr->param.dSunSoft,
		    sizeof(double), "sunSoft",
		    "<Softening length of the Sun in heliocentric coordinates> = 0.0");
#ifdef GASOLINE
	msr->param.bSphStep = 1;
	prmAddParam(msr->prm,"bSphStep",0,&msr->param.bSphStep,sizeof(int),
				"ss","<SPH timestepping>");
	msr->param.bDoGas = 1;
	prmAddParam(msr->prm,"bDoGas",0,&msr->param.bDoGas,sizeof(int),"gas",
				"calculate gas/don't calculate gas = +gas");
#ifndef GASOLINE
	msr->param.bDoGas = 0;
#endif
	msr->param.bGeometric = 0;
	prmAddParam(msr->prm,"bGeometric",0,&msr->param.bGeometric,sizeof(int),
				"geo","geometric/arithmetic mean to calc Grad(P/rho) = +geo");
	msr->param.bGasAdiabatic = 0;
	prmAddParam(msr->prm,"bGasAdiabatic",0,&msr->param.bGasAdiabatic,
				sizeof(int),"GasAdiabatic",
				"<Gas is Adiabatic> = +GasAdiabatic");
	msr->param.bGasIsothermal = 0;
	prmAddParam(msr->prm,"bGasIsothermal",0,&msr->param.bGasIsothermal,
				sizeof(int),"GasIsothermal",
				"<Gas is Isothermal> = +GasIsothermal");
	msr->param.bGasCooling = 0;
	prmAddParam(msr->prm,"bGasCooling",0,&msr->param.bGasCooling,
				sizeof(int),"GasCooling",
				"<Gas is Cooling> = +GasCooling");
	msr->param.iGasModel = GASMODEL_UNSET; /* Deprecated in for backwards compatibility */
	prmAddParam(msr->prm,"iGasModel",0,&msr->param.iGasModel,
				sizeof(int),"GasModel",
				"<Gas model employed> = 0 (Adiabatic)");
#ifndef NOCOOLING
	CoolAddParams( &msr->param.CoolParam, msr->prm );
#endif
	msr->param.dShockTrackerA = 0.16; 
	prmAddParam(msr->prm,"dShockTrackerA",2,&msr->param.dShockTrackerA,
				sizeof(double),"STA",
				"<Shock Tracker A constant> = 0.16");
	msr->param.dShockTrackerB = 0.4; 
	prmAddParam(msr->prm,"dShockTrackerB",2,&msr->param.dShockTrackerB,
				sizeof(double),"STB",
				"<Shock Tracker B constant> = 0.4");
	msr->param.dConstAlpha = 1.0; 	/* Default changed to 0.5 later if bBulkViscosity */
	prmAddParam(msr->prm,"dConstAlpha",2,&msr->param.dConstAlpha,
				sizeof(double),"alpha",
				"<Alpha constant in viscosity> = 1.0 or 0.5 (bBulkViscosity)");
	msr->param.dConstBeta = 2.0; 	/* Default changed to 0.5 later if bBulkViscosity */
	prmAddParam(msr->prm,"dConstBeta",2,&msr->param.dConstBeta,
				sizeof(double),"beta",
				"<Beta constant in viscosity> = 2.0 or 0.5 (bBulkViscosity)");
	msr->param.dConstGamma = 5.0/3.0;
	prmAddParam(msr->prm,"dConstGamma",2,&msr->param.dConstGamma,
				sizeof(double),"gamma",
				"<Ratio of specific heats> = 5/3");
	msr->param.dMeanMolWeight = 1.0;
	prmAddParam(msr->prm,"dMeanMolWeight",2,&msr->param.dMeanMolWeight,
				sizeof(double),"mmw",
				"<Mean molecular weight in amu> = 1.0");
	msr->param.dGasConst = 1.0;
	prmAddParam(msr->prm,"dGasConst",2,&msr->param.dGasConst,
				sizeof(double),"gcnst",
				"<Gas Constant>");
	msr->param.dKBoltzUnit = 1.0;
	prmAddParam(msr->prm,"dKBoltzUnit",2,&msr->param.dKBoltzUnit,
				sizeof(double),"gcnst",
				"<Boltzmann Constant in System Units>");
	msr->param.dMsolUnit = 1.0;
	prmAddParam(msr->prm,"dMsolUnit",2,&msr->param.dMsolUnit,
				sizeof(double),"msu",
				"<Solar mass/system mass unit>");
	msr->param.dKpcUnit = 1000.0;
	prmAddParam(msr->prm,"dKpcUnit",2,&msr->param.dKpcUnit,
				sizeof(double),"kpcu",
				"<Kiloparsec/system length unit>");
	msr->param.ddHonHLimit = 0.1;
	prmAddParam(msr->prm,"ddHonHLimit",2,&msr->param.ddHonHLimit,
				sizeof(double),"dhonh",
				"<|dH|/H Limiter> = 0.1");
	msr->param.bViscosityLimiter = 0;
	prmAddParam(msr->prm,"bViscosityLimiter",0,&msr->param.bViscosityLimiter,sizeof(int),
				"vlim","<Balsara Viscosity Limiter> = 0");
	msr->param.bViscosityLimitdt = 0;
	prmAddParam(msr->prm,"bViscosityLimitdt",0,&msr->param.bViscosityLimitdt,sizeof(int),
				"vlim","<Balsara Viscosity Limit dt> = 0");
	msr->param.bShockTracker = 0;
	prmAddParam(msr->prm,"bShockTracker",0,&msr->param.bShockTracker,sizeof(int),
				"st","<Shock Tracker> = 0");
	msr->param.bBulkViscosity = 0;
	prmAddParam(msr->prm,"bBulkViscosity",0,&msr->param.bBulkViscosity,sizeof(int),
				"bulk","<Bulk Viscosity> = 0");
	msr->param.bGasDomainDecomp = 0;
	prmAddParam(msr->prm,"bGasDomainDecomp",0,&msr->param.bGasDomainDecomp,sizeof(int),
				"gasDD","<Gas Domain Decomp> = 0");
	msr->param.bLowerSoundSpeed = 0;
	prmAddParam(msr->prm,"bLowerSoundSpeed",0,&msr->param.bLowerSoundSpeed,sizeof(int),
				"bLowc","<Lower Sound Speed> = 0");
	msr->param.bFastGas = 1;
	prmAddParam(msr->prm,"bFastGas",0,&msr->param.bFastGas,sizeof(int),
				"Fgas","<Fast Gas Method> = 1");
	msr->param.bStarForm = 0;
	prmAddParam(msr->prm,"bStarForm",0,&msr->param.bStarForm,sizeof(int),
				"stfm","<Star Forming> = 0");
	msr->param.bFeedBack = 0;
	prmAddParam(msr->prm,"bFeedBack",0,&msr->param.bFeedBack,sizeof(int),
				"fdbk","<Stars provide feedback> = 0");
	msr->param.bFormOutputs = 1;
	prmAddParam(msr->prm,"bFormOutputs",0,&msr->param.bFormOutputs,sizeof(int),
				"fdbk","<Write *form files?> = 0");
#ifdef SIMPLESF
	msr->param.SSF_dComovingDenMin = 200.0;
	prmAddParam(msr->prm,"SSF_dComovingDenMin", 2, &msr->param.SSF_dComovingDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
	msr->param.SSF_dPhysDenMin =  7e-26;
	prmAddParam(msr->prm,"SSF_dPhysDenMin", 2, &msr->param.SSF_dPhysDenMin,
		    sizeof(double), "stPDmin",
		    "<Minimum physical density for forming stars (gm/cc)> =  7e-26");
    msr->param.SSF_dInitStarMass = 0;
    prmAddParam(msr->prm,"SSF_dInitStarMass", 2, &msr->param.SSF_dInitStarMass,
				sizeof(double), "stm0",
				"<Initial star mass> = 0");
	msr->param.SSF_dESNPerStarMass = 1.25e16;
	prmAddParam(msr->prm,"SSF_dESNPerStarMass", 2, &msr->param.SSF_dESNPerStarMass,
		    sizeof(double), "ESNPerStarMass",
		    "<ESN per star mass, erg per g of stars> = 1.25e16");
	msr->param.SSF_dTMax = 3e4;
	prmAddParam(msr->prm,"SSF_dTMax", 2, &msr->param.SSF_dTMax,
		    sizeof(double), "SSF_dTMax",
		    "<Maximum temperature for forming stars, K> = 3e4");
	msr->param.SSF_dEfficiency = 0.1;
	prmAddParam(msr->prm,"SSF_dEfficiency", 2, &msr->param.SSF_dEfficiency,
		    sizeof(double), "SSF_dEfficiency",
		    "<SF Efficiency> = 0.1");
	msr->param.SSF_dtCoolingShutoff = 30e6;
	prmAddParam(msr->prm,"SSF_dtCoolingShutoff", 2, &msr->param.SSF_dtCoolingShutoff,
		    sizeof(double), "SSF_dtCoolingShutoff",
		    "<SF Cooling Shutoff duration> = 30e6");
	msr->param.SSF_bdivv = 1;
	prmAddParam(msr->prm,"SSF_bdivv", 0, &msr->param.SSF_bdivv,
		    sizeof(int), "SSF_bdivv",
		    "<SF Use div v for star formation> = 1");
#endif /* SIMPLESF */

#ifdef STARFORM
	stfmInitialize(&msr->param.stfm);
	msr->param.stfm->dOverDenMin = 2.0;
	prmAddParam(msr->prm,"dOverDenMin", 2, &msr->param.stfm->dOverDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
	msr->param.stfm->dPhysDenMin = 0.1;
	prmAddParam(msr->prm,"dPhysDenMin", 2, &msr->param.stfm->dPhysDenMin,
		    sizeof(double), "stPDmin",
		    "<Minimum physical density for forming stars (atoms/cc)> = .1");
	msr->param.stfm->dStarEff = .3333;
	prmAddParam(msr->prm,"dStarEff", 2, &msr->param.stfm->dStarEff,
		    sizeof(double), "stEff",
		    "<Fraction of gas converted into stars per timestep> = .3333");
    msr->param.stfm->dInitStarMass = 0;
    prmAddParam(msr->prm,"dInitStarMass", 2, &msr->param.stfm->dInitStarMass,
				sizeof(double), "stm0",
				"<Initial star mass> = 0");
	msr->param.stfm->dMinMassFrac = .1;
	prmAddParam(msr->prm,"dMinMassFrac", 2, &msr->param.stfm->dMinMassFrac,
		    sizeof(double), "stMinFrac",
		    "<Minimum fraction of average mass of neighbour particles required for gas particles to avoid deletion> = .1");
	msr->param.stfm->dMinGasMass = 0.0;
	prmAddParam(msr->prm,"dMinGasMass", 2, &msr->param.stfm->dMinGasMass,
		    sizeof(double), "stMinGas",
		    "<Minimum mass of a gas particle> = 0.0");
	msr->param.stfm->dMaxStarMass = 0.0;
	prmAddParam(msr->prm,"dMaxStarMass", 2, &msr->param.stfm->dMaxStarMass,
		    sizeof(double), "stMaxStarMass",
		    "<Maximum amount of star mass a hybrid particle can contain = 0.0");
	msr->param.stfm->dCStar = 0.05;
	prmAddParam(msr->prm,"dCStar", 2, &msr->param.stfm->dCStar,
		    sizeof(double), "stCStar",
		    "<Star formation coefficient> = 0.1");
	msr->param.stfm->dTempMax = 1.5e4;
	prmAddParam(msr->prm,"dTempMax", 2, &msr->param.stfm->dTempMax,
		    sizeof(double), "stTempMax",
		    "<Maximum temperature at which star formation occurs> = 0.0");
	msr->param.stfm->dSoftMin = 1.0;
	prmAddParam(msr->prm,"dSoftMin", 2, &msr->param.stfm->dSoftMin,
		    sizeof(double), "stSoftMin",
		    "<Minimum softening for star formation> = 0.0");
	msr->param.dDeltaStarForm = 1e6;
	prmAddParam(msr->prm,"dDeltaStarForm", 2, &msr->param.dDeltaStarForm,
		    sizeof(double), "dDeltaStarForm",
		    "<Minimum SF timestep in years> = 1e6");
	msr->param.bShortCoolShutoff = 0;
	prmAddParam(msr->prm,"bShortCoolShutoff", 0, &msr->param.bShortCoolShutoff,
		    sizeof(int), "bShortCoolShutoff",
		    "<Which cooling shutoff time to use> = long one");
	msr->param.iStarFormRung = 0;
	prmAddParam(msr->prm,"iStarFormRung", 2, &msr->param.iStarFormRung,
		    sizeof(double), "iStarFormRung",
		    "<Star Formation Rung> = 0");

/* supernova constants */
	fbInitialize(&msr->param.fb);
        snInitialize(&msr->param.sn);
        snInitConstants(msr->param.sn);
	msr->param.sn->dESN = 0.1e51;
	prmAddParam(msr->prm,"dESN", 2, &msr->param.sn->dESN,
		    sizeof(double), "snESN",
		    "<Energy of supernova in ergs> = 0.1e51");
	if (msr->param.sn->dESN > 0.0) msr->param.bSmallSNSmooth = 1;
        else msr->param.bSmallSNSmooth = 0;
	prmAddParam(msr->prm,"bSmallSNSmooth", 0, &msr->param.bSmallSNSmooth,
		    sizeof(int), "bSmallSNSmooth",
		    "<smooth SN ejecta over blast or smoothing radius> = blast radius");
        msr->param.bSNTurnOffCooling = 1;
	prmAddParam(msr->prm,"bSNTurnOffCooling", 0, &msr->param.bSNTurnOffCooling,
		    sizeof(int), "bSNTurnOffCooling",
		    "<Do SN turn off cooling> = 1");

#endif /* STARFORM */
#endif /* GASOLINE */
#ifdef GLASS
	msr->param.dGlassDamper = 0.0;
	prmAddParam(msr->prm,"dGlassDamper",2,&msr->param.dGlassDamper,
		    sizeof(double),"dGlassDamper",
		    "Lose 0.0 dt velocity per step (no damping)");
	msr->param.dGlassPoverRhoL = 1.0;
	prmAddParam(msr->prm,"dGlassPoverRhoL",2,&msr->param.dGlassPoverRhoL,
				sizeof(double),"dGlassPoverRhoL","Left P on Rho = 1.0");
	msr->param.dGlassPoverRhoL = 1.0;
	prmAddParam(msr->prm,"dGlassPoverRhoR",2,&msr->param.dGlassPoverRhoR,
				sizeof(double),"dGlassPoverRhoR","Right P on Rho = 1.0");
	msr->param.dGlassxL = 0.0;
	prmAddParam(msr->prm,"dGlassxL",2,&msr->param.dGlassxL,sizeof(double),
				"dGlassxL","Left smoothing range = 0.0");
	msr->param.dGlassxR = 0.0;
	prmAddParam(msr->prm,"dGlassxR",2,&msr->param.dGlassxR,sizeof(double),
				"dGlassxR","Right smoothing range = 0.0");
	msr->param.dGlassVL = 0.0;
	prmAddParam(msr->prm,"dGlassVL",2,&msr->param.dGlassVL,sizeof(double),
				"dGlassVL","Left Max Random Velocity = 0.0");
	msr->param.dGlassVR = 0.0;
	prmAddParam(msr->prm,"dGlassVR",2,&msr->param.dGlassVR,sizeof(double),
				"dGlassVR","Right Max Random Velocity = 0.0");
#endif /* GLASS */

#ifdef SLIDING_PATCH
	msr->param.PP.bPatch = 0;
	prmAddParam(msr->prm,"bPatch",0,&msr->param.PP.bPatch,sizeof(int),
				"patch","enable/disable patch reference frame = -patch");
	msr->param.PP.dOrbDist = 0.0;
	prmAddParam(msr->prm,"dOrbDist",2,&msr->param.PP.dOrbDist,
				sizeof(double),"orbdist","<Patch orbital distance>");
	msr->param.PP.bExtPert = 0;
	prmAddParam(msr->prm,"bExtPert",0,&msr->param.PP.bExtPert,
				sizeof(int),"extpert","enable/disable patch perturber = -extpert");
	msr->param.PP.dPertOrbDist = 0.0;
	prmAddParam(msr->prm,"dPertOrbDist",2,&msr->param.PP.dPertOrbDist,
				sizeof(double),"pertorbdist","<Perturber orbital distance>");
	msr->param.PP.dPertMass = 0.0;
	prmAddParam(msr->prm,"dPertMass",2,&msr->param.PP.dPertMass,
				sizeof(double),"pertmass","<Perturber mass>");
	msr->param.PP.dPertMaxZ = 0.0;
	prmAddParam(msr->prm,"dPertMaxZ",2,&msr->param.PP.dPertMaxZ,
				sizeof(double),"pertmaxz","<Perturber maximum z height>");
	msr->param.PP.dPertOrbFreqZ = 0.0;
	prmAddParam(msr->prm,"dPertOrbFreqZ",2,&msr->param.PP.dPertOrbFreqZ,
				sizeof(double),"pertorbfreqz","<Perturber vertical orbital frequency>");
	msr->param.PP.dPertPhase = 0.0;
	prmAddParam(msr->prm,"dPertPhase",2,&msr->param.PP.dPertPhase,
				sizeof(double),"pertphase","<Perturber orbital phase>");
	msr->param.PP.dPertPhaseZ = 0.0;
	prmAddParam(msr->prm,"dPertPhaseZ",2,&msr->param.PP.dPertPhaseZ,
				sizeof(double),"pertphasez","<Perturber vertical orbital phase>");
	msr->param.PP.bRandAzWrap = 0;
	prmAddParam(msr->prm,"bRandAzWrap",0,&msr->param.PP.bRandAzWrap,
				sizeof(int),"randazwrap","enable/disable random azimuthal wrap = -randazwrap");
	msr->param.PP.bNoRandomX = 0;
	prmAddParam(msr->prm,"bNoRandomX",0,&msr->param.PP.bNoRandomX,
				sizeof(int),"norandx","enable/disable no randomization of x component = -norandx");
	msr->param.PP.nWrapAttempts = 1;
	prmAddParam(msr->prm,"nWrapAttempts",1,&msr->param.PP.nWrapAttempts,
				sizeof(int),"nwraptries","<Number of wrap attempts (without overlap) before giving up> = 1");
	msr->param.PP.iStripOption = 2;
	prmAddParam(msr->prm,"iStripOption",1,&msr->param.PP.iStripOption,
				sizeof(int),"stripoption","<Strip option> = 2");
	msr->param.PP.dStripInner = 0.0;
	prmAddParam(msr->prm,"dStripInner",2,&msr->param.PP.dStripInner ,
				sizeof(double),"stripinner","<Inner edge of strip in patch width units for random azimuthal wrap> = 0.0");
	msr->param.PP.dStripOuter = 0.5;
	prmAddParam(msr->prm,"dStripOuter",2,&msr->param.PP.dStripOuter,
				sizeof(double),"stripouter","<Outer edge of strip in patch width units for random azimuthal wrap> = 0.5");
	msr->param.PP.dVelDispX = 0.0;
	prmAddParam(msr->prm,"dVelDispX",2,&msr->param.PP.dVelDispX,
				sizeof(double),"veldispx","<Radial velocity dispersion for random azimuthal wrap>");
	msr->param.PP.dVelDispY = 0.0;
	prmAddParam(msr->prm,"dVelDispY",2,&msr->param.PP.dVelDispY,
				sizeof(double),"veldispy","<Azimuthal velocity dispersion for random azimuthal wrap>");
	msr->param.PP.dAvgVertAmp = 0.0;
	prmAddParam(msr->prm,"dAvgVertAmp",2,&msr->param.PP.dAvgVertAmp,
				sizeof(double),"avgvertamp","<Mean vertical oscillation amplitude for random azimuthal wrap>");	
	msr->param.PP.dAvgMass = 0.0;
	prmAddParam(msr->prm,"dAvgMass",2,&msr->param.PP.dAvgMass,
				sizeof(double),"avgmass","<Average particle mass for random azimuthal wrap>");
#endif /* SLIDING_PATCH */

#ifdef SIMPLE_GAS_DRAG
	msr->param.bSimpleGasDrag = 0;
	prmAddParam(msr->prm,"bSimpleGasDrag",0,&msr->param.bSimpleGasDrag,
				sizeof(int),"sgd","enable/disable simple gas drag = -sgd");
	msr->param.bEpstein = 1;
	prmAddParam(msr->prm,"bEpstein",0,&msr->param.bEpstein,sizeof(int),
				"epstein","enable/disable Epstein regime for gas drag");
	msr->param.dGamma = 1.0e-11;
	prmAddParam(msr->prm,"dGamma",2,&msr->param.dGamma,sizeof(double),
				"gamma","coefficient for inverse stopping time = 1.0e-11");
#endif /* SIMPLE_GAS_DRAG */

#ifdef SPINUP
	msr->param.bUseAngMomSpinup = 0;
	prmAddParam(msr->prm,"bUseAngMomSpinup",0,&msr->param.bUseAngMomSpinup,
				sizeof(int),"spinupm","use constant angular momentum increment = 0");
	msr->param.bSpinupAlignAxes = 1;
	prmAddParam(msr->prm,"bSpinupAlignAxes",0,&msr->param.bSpinupAlignAxes,
				sizeof(int),"spinupalign","align spinup axes = True");
	msr->param.dSpinupPeriodStart = 6.0;
	prmAddParam(msr->prm,"dSpinupPeriodStart",2,&msr->param.dSpinupPeriodStart,
				sizeof(double),"spinupps","start spin period in hours = 6.0");
	msr->param.dSpinupPeriodInter = 3.0;
	prmAddParam(msr->prm,"dSpinupPeriodInter",2,&msr->param.dSpinupPeriodInter,
				sizeof(double),"spinuppi","intermediate spin period in hours = 3.0");
	msr->param.dSpinupPeriodFinal = 2.26;
	prmAddParam(msr->prm,"dSpinupPeriodFinal",2,&msr->param.dSpinupPeriodFinal,
				sizeof(double),"spinuppf","final spin period in hours = 2.26");
	msr->param.nSpinupStepsStart = 0;
	prmAddParam(msr->prm,"nSpinupStepsStart",1,&msr->param.nSpinupStepsStart,
				sizeof(int),"spinupss","number of steps before starting spinup = 0");
	msr->param.nSpinupStepsInter = 0;
	prmAddParam(msr->prm,"nSpinupStepsInter",1,&msr->param.nSpinupStepsInter,
				sizeof(int),"spinupsi","number of steps from start to intermediate spinup = 0");
	msr->param.nSpinupStepsFinal = 0;
	prmAddParam(msr->prm,"nSpinupStepsFinal",1,&msr->param.nSpinupStepsFinal,
				sizeof(int),"spinupsf","number of steps from intermediate to final spinup = n");
	msr->param.nSpinupStepsHold = 0;
	prmAddParam(msr->prm,"nSpinupStepsHold",1,&msr->param.nSpinupStepsHold,
				sizeof(int),"spinupsh","number of steps to hold at final spin = 0");
	msr->param.nSpinupStepsRestNum = 0;
	prmAddParam(msr->prm,"nSpinupStepsRestNum",1,&msr->param.nSpinupStepsRestNum,
				sizeof(int),"spinupsrn","number of rest to settle down");
	msr->param.nSpinupStepsRestTime = 0;
	prmAddParam(msr->prm,"nSpinupStepsRestTime",1,&msr->param.nSpinupStepsRestTime,
				sizeof(int),"spinupsrt","number of steps to settle down for each rest");
	msr->param.dSpinupRadiusInner = 0.4;
	prmAddParam(msr->prm,"dSpinupRadiusInner",2,&msr->param.dSpinupRadiusInner,
				sizeof(double),"spinupri","inner radius for spinup in km = 0.4");
	msr->param.dSpinupRadiusOuter = 0.8;
	prmAddParam(msr->prm,"dSpinupRadiusOuter",2,&msr->param.dSpinupRadiusOuter,
				sizeof(double),"spinupro","outer radius for spinup in km = 0.8");
        msr->param.dSpinupCoheCoeff = 0;
	prmAddParam(msr->prm,"dSpinupCoheCoeff",2,&msr->param.dSpinupCoheCoeff,
				sizeof(double),"spinupcc","increment of cohesion coefficient");
	msr->param.nSpinupCoheTime = 0;
	prmAddParam(msr->prm,"nSpinupCoheTime",1,&msr->param.nSpinupCoheTime,
				sizeof(int),"spinupct","number of steps to hold a certain cohesion coefficient");
#endif /* SPINUP */
	
#ifdef COLLISIONS
	msr->param.bAllowSimulColl = 1;
	prmAddParam(msr->prm,"bAllowSimulColl",0,&msr->param.bAllowSimulColl,
				sizeof(int),"simulcoll","allow/disallow simultaneous collisions = -simulcoll");
	msr->param.bFindRejects = 0;
	prmAddParam(msr->prm,"bFindRejects",0,&msr->param.bFindRejects,
				sizeof(int),"rejects","enable/disable check for rejected ICs = -rejects");
	msr->param.iCollLogOption = 0;
	prmAddParam(msr->prm,"iCollLogOption",1,&msr->param.iCollLogOption,
				sizeof(int),"clog","<Collision log option> = 0");
	msr->param.dSmallStep = 0.0;
	prmAddParam(msr->prm,"dSmallStep",2,&msr->param.dSmallStep,
				sizeof(double),"sstep","<Rectilinear time-step>");
	msr->param.dxUnifGrav = msr->param.dyUnifGrav = msr->param.dzUnifGrav = 0;
	prmAddParam(msr->prm,"dxUnifGrav",2,&msr->param.dxUnifGrav,sizeof(double),
				"gx","<x component of uniform gravity field> = 0");
	prmAddParam(msr->prm,"dyUnifGrav",2,&msr->param.dyUnifGrav,sizeof(double),
				"gy","<y component of uniform gravity field> = 0");
	prmAddParam(msr->prm,"dzUnifGrav",2,&msr->param.dzUnifGrav,sizeof(double),
				"gz","<z component of uniform gravity field> = 0");
	msr->param.CP.iOutcomes = BOUNCE;
	prmAddParam(msr->prm,"iOutcomes",1,&msr->param.CP.iOutcomes,
				sizeof(int),"outcomes","<Allowed collision outcomes> = 0");
	msr->param.CP.dMergeLimit = 1.0;
	prmAddParam(msr->prm,"dMergeLimit",2,&msr->param.CP.dMergeLimit,
				sizeof(double),"mlim","<Merge limit> = 1.0");
	msr->param.CP.dDensity = 0.0;
	prmAddParam(msr->prm,"dDensity",2,&msr->param.CP.dDensity,
				sizeof(double),"density","<Merged particle density> = 0");
	msr->param.CP.iDensityAltCol = 0;
	prmAddParam(msr->prm,"iDensityAltCol",1,&msr->param.CP.iDensityAltCol,
				sizeof(int),"denaltcol","<Color of alternate density particle> = 0");
	msr->param.CP.dDensityAltVal = 0.0;
	prmAddParam(msr->prm,"dDensityAltVal",2,&msr->param.CP.dDensityAltVal,
				sizeof(double),"denaltval","<Alternate merged particle density> = 0");
	msr->param.CP.iEpsNOption = ConstEps;
	prmAddParam(msr->prm,"iEpsNOption",1,&msr->param.CP.iEpsNOption,
				sizeof(int),"epsnopt","<EpsN option> = 0");
	msr->param.CP.dEpsN = 1.0;
	prmAddParam(msr->prm,"dEpsN",2,&msr->param.CP.dEpsN,
				sizeof(double),"epsn","<Coefficient of restitution> = 1");
	msr->param.CP.dEpsT = 1.0;
	prmAddParam(msr->prm,"dEpsT",2,&msr->param.CP.dEpsT,
				sizeof(double),"epst","<Coefficient of surface friction> = 1");
	msr->param.CP.dEpsNCoef = 0.52;
	prmAddParam(msr->prm,"dEpsNCoef",2,&msr->param.CP.dEpsNCoef,
				sizeof(double),"epsncoef","<Coefficient of EpsN power law> = 0.52");
	msr->param.CP.dEpsNExp = -0.14;
	prmAddParam(msr->prm,"dEpsNExp",2,&msr->param.CP.dEpsNExp,
				sizeof(double),"epsnexp","<Exponent of EpsN power law> = -0.14");
	msr->param.CP.dEpsNVStar = 0.01;
	prmAddParam(msr->prm,"dEpsNVStar",2,&msr->param.CP.dEpsNVStar,
				sizeof(double),"epsnvstar","<VStar in Borderies EpsN law> = 0.01");
	msr->param.CP.dEpsNMin = 0.01;
	prmAddParam(msr->prm,"dEpsNMin",2,&msr->param.CP.dEpsNMin,
				sizeof(double),"epsnmind","<Minimum EpsN for variable EpsN laws> = 0.01");
	msr->param.CP.iSlideOption = EscVel;
	prmAddParam(msr->prm,"iSlideOption",1,&msr->param.CP.iSlideOption,
				sizeof(int),"sopt","<Slide option> = 0");
	msr->param.CP.dSlideLimit = 0.0;
	prmAddParam(msr->prm,"dSlideLimit",2,&msr->param.CP.dSlideLimit,
				sizeof(double),"slide","<Sliding motion parameter> = 0");
	msr->param.CP.dSlideEpsN = 1.0;
	prmAddParam(msr->prm,"dSlideEpsN",2,&msr->param.CP.dSlideEpsN,
				sizeof(double),"sepsn","<epsn if speed less than minimum> = 1");
	msr->param.CP.dSlideEpsT = 1.0;
	prmAddParam(msr->prm,"dSlideEpsT",2,&msr->param.CP.dSlideEpsT,
				sizeof(double),"sepst","<epst if speed less than minimum> = 1");
	msr->param.CP.dCollapseLimit = 0.0;
	prmAddParam(msr->prm,"dCollapseLimit",2,&msr->param.CP.dCollapseLimit,
				sizeof(double),"collapse","<Inelastic collapse parameter> = 0");
	msr->param.CP.dCollapseEpsN = 1.0;
	prmAddParam(msr->prm,"dCollapseEpsN",2,&msr->param.CP.dCollapseEpsN,
				sizeof(double),"cepsn","<epsn for inelastic collapse> = 1");
	msr->param.CP.dCollapseEpsT = 1.0;
	prmAddParam(msr->prm,"dCollapseEpsT",2,&msr->param.CP.dCollapseEpsT,
				sizeof(double),"cepst","<epst for inelastic collapse> = 1");
	msr->param.CP.dCrushLimit = 0.0; /* i.e. no limit, set to DBL_MAX below */
	prmAddParam(msr->prm,"dCrushLimit",2,&msr->param.CP.dCrushLimit,
				sizeof(double),"crush","<Maximum impact speed squared> = 0");
	msr->param.CP.dCrushEpsN = 1.0;
	prmAddParam(msr->prm,"dCrushEpsN",2,&msr->param.CP.dCrushEpsN,
				sizeof(double),"cepsn","<epsn if speed greater than maximum> = 1");
	msr->param.CP.dCrushEpsT = 1.0;
	prmAddParam(msr->prm,"dCrushEpsT",2,&msr->param.CP.dCrushEpsT,
				sizeof(double),"cepst","<epst if speed greater than maximum> = 1");
	msr->param.CP.iOverlapOption = OverlapIsError;
	prmAddParam(msr->prm,"iOverlapOption",1,&msr->param.CP.iOverlapOption,
				sizeof(int),"overlap","<Overlap option> = 0");
	msr->param.CP.bStrictOverlap = FALSE;
	prmAddParam(msr->prm,"bStrictOverlap",0,&msr->param.CP.bStrictOverlap,
				sizeof(int),"strictoverlap","enable/disable strict interpretation of overlaps");
	msr->param.CP.dBackstepLimit = 0.0;
	prmAddParam(msr->prm,"dBackstepLimit",2,&msr->param.CP.dBackstepLimit,
				sizeof(double),"bsteplim","<Backstep limit> = 0.0");
	msr->param.CP.dAdjPosLimit = 0.0;
	prmAddParam(msr->prm,"dAdjPosLimit",2,&msr->param.CP.dAdjPosLimit,
				sizeof(double),"adjposlim","<Adjust position limit> = 0.0");
	msr->param.CP.dRepelFac = 0.0;
	prmAddParam(msr->prm,"dRepelFac",2,&msr->param.CP.dRepelFac,
				sizeof(double),"repelfac","<Repel factor> = 0.0");
	msr->param.CP.dFragLimit = 0.0;
	prmAddParam(msr->prm,"dFragLimit",2,&msr->param.CP.dFragLimit,
				sizeof(double),"fraglim","<Frag limit> = 0.0");

	/* optional time-variable parameters */

	strcpy(achUnifGravFile,""); /* same as *achUnifGravFile = '\0'; */
	prmAddParam(msr->prm,"achUnifGravFile",3,achUnifGravFile,MAXPATHLEN,
				"unifgravfile","<name of uniform gravity data file> = \"\"");

	/*
	 ** Added to ss.par by ZML 08.20.03 - for resolved rubble collisions
	 */

#ifdef RUBBLE_ZML
	msr->param.CP.dRubbleMinFracMass = 0.2;
	prmAddParam(msr->prm,"dRubbleMinFracMass",2,&msr->param.CP.dRubbleMinFracMass,
				sizeof(double),"rubminfrac","<rubminfrac collision outcome criterion> = 0.2");
	msr->param.CP.DB.dRubMinMass = 1.5e-10;
	prmAddParam(msr->prm,"dRubMinMass",2,&msr->param.CP.DB.dRubMinMass,sizeof(double),
				"rubmin","<rubmin max size put into dust> = 1.5e-10");
	msr->param.CP.DB.iRubNumDynToBounce = 10;
	prmAddParam(msr->prm,"iRubNumDynToBounce",1,&msr->param.CP.DB.iRubNumDynToBounce,
				sizeof(int),"rndtb","<rndtb number of dyn times bouncing> = 10");
	msr->param.CP.DB.iRubNumDynToMerge = 10;
	prmAddParam(msr->prm,"iRubNumDynToMerge",1,&msr->param.CP.DB.iRubNumDynToMerge,
				sizeof(int),"rndtm","<rndtm number of dyn times merging> = 10");
	msr->param.CP.DB.nDustBins = 10;
	prmAddParam(msr->prm,"nDustBins",1,&msr->param.CP.DB.nDustBins,sizeof(int),
				"ndb","<Number of dust bins (0 disables)> = 0");
	msr->param.CP.DB.iDustBinsApplyInt = 10; /* should be < # steps/orbits */
	prmAddParam(msr->prm,"iDustBinsApplyInt",1,&msr->param.CP.DB.iDustBinsApplyInt,
				sizeof(int),"dba","<Dust bins application interval (> 0)> = 10");
	msr->param.CP.DB.iDustBinsVelDispOpt = 0;
	prmAddParam(msr->prm,"iDustBinsVelDispOpt",1,&msr->param.CP.DB.iDustBinsVelDispOpt,
				sizeof(int),"dbv","<Dust bins velocity dispersion (0,mean,max)> = 0");
	msr->param.CP.DB.dDustBinsInner = 0.0; /* in AU */
	prmAddParam(msr->prm,"dDustBinsInner",2,&msr->param.CP.DB.dDustBinsInner,
				sizeof(double),"dbi","<Inner dust bin radius (AU)> = 0.5");
	msr->param.CP.DB.dDustBinsOuter = 2.0; /* in AU */
	prmAddParam(msr->prm,"dDustBinsOuter",2,&msr->param.CP.DB.dDustBinsOuter,
				sizeof(double),"dbo","<Outer dust bin radius (AU)> = 2.0");
	msr->param.CP.DB.dDustBinsScaleHeight = 1.0e-4; /* in AU */
	prmAddParam(msr->prm,"dDustBinsScaleHeight",2,&msr->param.CP.DB.dDustBinsScaleHeight,
				sizeof(double),"dbsh","<Dust bins scale height (AU)> = 1.0e-4");
	msr->param.CP.DB.dDustBinsInitSigma = 0.1; /* in g/cm^2 */
	prmAddParam(msr->prm,"dDustBinsInitSigma",2,&msr->param.CP.DB.dDustBinsInitSigma,
				sizeof(double),"dbis","<Init. dust surf. mass den. at 1 AU (g/cm^2)> = 0.1");
	msr->param.CP.DB.dDustBinsInitAlpha = -1.5;
	prmAddParam(msr->prm,"dDustBinsInitAlpha",2,&msr->param.CP.DB.dDustBinsInitAlpha,
				sizeof(double),"dbia","<Init. dust surf. mass den. power law exponent> = -1.5");
#elif defined(COLLMOD_ZML)
	msr->param.bCollDelay = 0;
	prmAddParam(msr->prm,"bCollDelay",0,&msr->param.bCollDelay,
				sizeof(int),"colldelay","enable/disable collision delay = -delay");
	msr->param.dCollTimer = 6.2;
	prmAddParam(msr->prm,"dCollTimer",2,&msr->param.dCollTimer,
				sizeof(double),"colltimer","<Collision delay timer> = 6.2");
	msr->param.CP.DB.dCollMinMass = 8.4e-12;
	prmAddParam(msr->prm,"dCollMinMass",2,&msr->param.CP.DB.dCollMinMass,
				sizeof(double),"collmodmin","<Collmodmin max size put into dust> = 1.5e-10");
	msr->param.CP.DB.nDustBins = 10;
	prmAddParam(msr->prm,"nDustBins",1,&msr->param.CP.DB.nDustBins,
				sizeof(int),"ndb","<Number of dust bins (0 disables)> = 10");
	msr->param.CP.DB.iDustBinsApplyInt = 10; /* should be < # steps/orbits */
	prmAddParam(msr->prm,"iDustBinsApplyInt",1,&msr->param.CP.DB.iDustBinsApplyInt,
				sizeof(int),"dba","<Dust bins application interval (> 0)> = 10");
	msr->param.CP.DB.iDustBinsVelDispOpt = 0; 
	prmAddParam(msr->prm,"iDustBinsVelDispOpt",1,&msr->param.CP.DB.iDustBinsVelDispOpt,
				sizeof(int),"dbv","<Dust bins velocity dispersion (0,mean,max)> = 0");
	msr->param.CP.DB.dDustBinsInner = 0.0; /* in AU */
	prmAddParam(msr->prm,"dDustBinsInner",2,&msr->param.CP.DB.dDustBinsInner,
				sizeof(double),"dbi","<Inner dust bin radius (AU)> = 0.5");
	msr->param.CP.DB.dDustBinsOuter = 2.0; /* in AU */
	prmAddParam(msr->prm,"dDustBinsOuter",2,&msr->param.CP.DB.dDustBinsOuter,
				sizeof(double),"dbo","<Outer dust bin radius (AU)> = 2.0");
	msr->param.CP.DB.dDustBinsScaleHeight = 1.0e-4; /* in AU */
	prmAddParam(msr->prm,"dDustBinsScaleHeight",2,&msr->param.CP.DB.dDustBinsScaleHeight,
				sizeof(double),"dbsh","<Dust bins scale height (AU)> = 1.0e-4");
	msr->param.CP.DB.dDustBinsInitSigma = 0.1; /* in g/cm^2 */
	prmAddParam(msr->prm,"dDustBinsInitSigma",2,&msr->param.CP.DB.dDustBinsInitSigma,
				sizeof(double),"dbis","<Init. dust surf. mass den. at 1 AU (g/cm^2)> = 0.1");
	msr->param.CP.DB.dDustBinsInitAlpha = -1.5;
	prmAddParam(msr->prm,"dDustBinsInitAlpha",2,&msr->param.CP.DB.dDustBinsInitAlpha,
				sizeof(double),"dbia","<Init. dust surf. mass den. power law exponent> = -1.5");
#endif /* RUBBLE_ZML, COLLMOD_ZML */

#ifdef AGGS
	msr->param.CP.bAggsSolveQuartic = 0;
	prmAddParam(msr->prm,"bAggsSolveQuartic",0,&msr->param.CP.bAggsSolveQuartic,sizeof(int),
				"aggssolvequartic","enable/disable quartic solver for aggs = -aggssolvequartic");	
    msr->param.CP.SP.dTensileCoef = -1.0; /* i.e. infinitely rigid */
    prmAddParam(msr->prm,"dTensileCoef",2,&msr->param.CP.SP.dTensileCoef,
                sizeof(double),"tensilecoef","<tensile strength coefficient> = -1");
    msr->param.CP.SP.dTensileExp = 0.0;
    prmAddParam(msr->prm,"dTensileExp",2,&msr->param.CP.SP.dTensileExp,
                sizeof(double),"tensileexp","<tensile strength exponent> = 0");
    msr->param.CP.SP.dShearCoef = -1.0;
    prmAddParam(msr->prm,"dShearCoef",2,&msr->param.CP.SP.dShearCoef,
                sizeof(double),"shearcoef","<shear strength coefficient> = -1");
    msr->param.CP.SP.dShearExp = 0.0;
    prmAddParam(msr->prm,"dShearExp",2,&msr->param.CP.SP.dShearExp,
                sizeof(double),"shearexp","<shear strength exponent> = 0");
#endif /* AGGS */

#ifdef WALLS
	strcpy(achWallsFile,"");
	prmAddParam(msr->prm,"achWallsFile",3,achWallsFile,MAXPATHLEN,
				"wallsfile","<name of wall data file> = \"\"");
	msr->param.CP.WP.bWallsEdgeDetect = 1;
	prmAddParam(msr->prm,"bWallsEdgeDetect",0,&msr->param.CP.WP.bWallsEdgeDetect,sizeof(int),
				"wallsedgedetect","enable/disable finite wall edge detection = +wallsedgedetect");		
	msr->param.CP.WP.bWallsSolveQuartic = 0;
	prmAddParam(msr->prm,"bWallsSolveQuartic",0,&msr->param.CP.WP.bWallsSolveQuartic,sizeof(int),
				"wallssolvequartic","enable/disable quartic solver for particles stuck on rotating walls = -wallssolvequartic");	
#endif /* WALLS */

#ifdef SPECIAL_PARTICLES
	strcpy(achSpecialFile,"");
	prmAddParam(msr->prm,"achSpecialFile",3,achSpecialFile,MAXPATHLEN,
				"specialfile","<name of special particle file> = \"\"");
#endif

#endif /* COLLISIONS */

	/* force override options */

	msr->param.FO.iForceOverrideOption = 0;
	prmAddParam(msr->prm,"iForceOverrideOption",1,&msr->param.FO.iForceOverrideOption,
				sizeof(int),"overrideopt","<Force override option> = 0");

#ifdef SPRINGS
	{
		SPRING_PARAMS *SP = &msr->param.FO.SP; /* shorthand */
		SP->dMeanYoungsModulus = 0.0;
		prmAddParam(msr->prm,"dMeanYoungsModulus",2,&SP->dMeanYoungsModulus,
					sizeof(double),"youngsmod","<Young's modulus>");
		SP->dMeanStrainLimit = 0.0;
		prmAddParam(msr->prm,"dMeanStrainLimit",2,&SP->dMeanStrainLimit,
					sizeof(double),"strainlim","<Strain limit>");
		SP->dYoungsStdDev = 0.0;
		prmAddParam(msr->prm,"dYoungsStdDev",2,&SP->dYoungsStdDev,
					sizeof(double),"youngsstddev","<Standard deviation of Young's modulus>");
		SP->dStrainStdDev = 0.0;
		prmAddParam(msr->prm,"dStrainStdDev",2,&SP->dStrainStdDev,
					sizeof(double),"strainstddev","<Standard deviation of Strain limit>");
		SP->dMaxStrainLimit = 1.0;
		prmAddParam(msr->prm,"dMaxStrainLimit",2,&SP->dMaxStrainLimit,
					sizeof(double),"maxstrainlim","<Maximum allowed Strain limit>");
		SP->dLinkageLength = 2.5;
		prmAddParam(msr->prm,"dLinkageLength",2,&SP->dLinkageLength,
					sizeof(double),"linkagelength","<Linkage length> = 2.5 effective radii");
		SP->dZeroStrainLength = 0.0;
		prmAddParam(msr->prm,"dZeroStrainLength",2,&SP->dZeroStrainLength,
					sizeof(double),"zerostrainlength","<Zero strain length>");
		SP->dDamp = 0.0;
		prmAddParam(msr->prm,"dDamp",2,&SP->dDamp,
					sizeof(double),"dDamp","<Spring Dampening Factor>");
		SP->dPorosityFactor = 1.0;
		prmAddParam(msr->prm,"dPorosityFactor",2,&SP->dPorosityFactor,
					sizeof(double),"porfac","<Porosity factor> = 1");
		SP->bReadSpringsData = 0;
		prmAddParam(msr->prm,"bReadSpringsData",0,&SP->bReadSpringsData,
					sizeof(int),"readsprings","enable/disable read springs data = -readsprings");
		SP->bNoSpringsOutput = 0;
		prmAddParam(msr->prm,"bNoSpringsOutput",0,&SP->bNoSpringsOutput,
					sizeof(int),"nospringsoutput","enable/disable write springs data = -nospringsoutput");
		}
#endif /* SPRINGS */

#ifdef DEM
	{
		DEM_PARAMS *DP = &msr->param.FO.DP; /* shorthand */
		DP->dKn = 0.0; /* (m)(v_max^2)/(x_max^2) */
		prmAddParam(msr->prm,"dKn",2,&DP->dKn,
					sizeof(double),"kn","<kn> = 0.0 (no force)");
		DP->dKt = -1; /* Kt default is (2/7) kn */
		prmAddParam(msr->prm,"dKt",2,&DP->dKt,
					sizeof(double),"kt","<kt> = -1");
		DP->dMuS = 0.0; /* coefficient of static friction */
		prmAddParam(msr->prm,"dMuS",2,&DP->dMuS,
					sizeof(double),"mus","<mus> = 0.0");
		DP->dMuR = 0.0; /* coefficient of rolling friction */
		prmAddParam(msr->prm,"dMuR",2,&DP->dMuR,
					sizeof(double),"mur","<mur> = 0.0");
		DP->dMuT = 0.0; /* coefficient of twisting friction */
		prmAddParam(msr->prm,"dMuT",2,&DP->dMuT,
					sizeof(double),"mut","<mut> = 0.0");
		DP->dAccCrit = 0.0; /* minimum net acceleration */
		prmAddParam(msr->prm,"dAccCrit",2,&DP->dAccCrit,
					sizeof(double),"acccrit","<acccrit> = 0.0");		
		DP->dTangentialSpringDrag = 0.0; /* fraction of spring remaining after exceeding static friction */
		/*DEBUG: may want to make this rate of "spring loss", e.g., how long before spring is 1/2 or 1/e */
		prmAddParam(msr->prm,"dTangentialSpringDrag",2,&DP->dTangentialSpringDrag,
					sizeof(double),"dtangspringdrag","<Tangential spring drag> = 0.0");
		DP->dMinorFrac = 0.03;
		prmAddParam(msr->prm,"dMinorFrac",2,&DP->dMinorFrac,
                    sizeof(double),"demminor","<demminor> = 0.03");
		DP->dMajorFrac = 0.15;
		prmAddParam(msr->prm,"dMajorFrac",2,&DP->dMajorFrac,
                    sizeof(double),"demmajor","<demmajor> = 0.15");
		DP->dErrorFrac = 0.50;
		prmAddParam(msr->prm,"dErrorFrac",2,&DP->dErrorFrac,
                    sizeof(double),"demerror","<demerror> = 0.50");
		DP->bReadDEMData = 0;
		prmAddParam(msr->prm,"bReadDEMData",0,&DP->bReadDEMData,
					sizeof(int),"readdem","enable/disable read DEM data = -readdem");
#ifdef DEM_TRACK_ORIENT
		DP->bReadOrient = 0;
		prmAddParam(msr->prm,"bReadOrient",0,&DP->bReadOrient,
					sizeof(int),"readori","enable/disable read .ori data = -readori");
#endif /* DEM_TRACK_ORIENT */
		DP->nDEMOutputs = 1;
		prmAddParam(msr->prm,"nDEMOutputs",1,&DP->nDEMOutputs,
					sizeof(int),"ndemoutputs","<number of DEM outputs to keep> = 1");
		DP->iDEMStatsInterval = 0;
		prmAddParam(msr->prm,"iDEMStatsInterval",1,&DP->iDEMStatsInterval,sizeof(int),
					"si","<number of timesteps between DEM stats> = 0");
#ifdef DEM_FIXED_BALL
		DP->iFixedBallOption = 0;
		prmAddParam(msr->prm,"iFixedBallOption",1,&DP->iFixedBallOption,
					sizeof(int),"fboption","<DEM Fixed Ball Option> = 0");
#endif /* DEM_FIXED_BALL */
#ifdef DEM_ROTATION_DASHPOT
		DP->dBetaS = 0.5; /* shape parameter used in spring-dashpot-slider rotation model */
		prmAddParam(msr->prm,"dBetaS",2,&DP->dBetaS,
					sizeof(double),"betas","<betas> = 0.5");
#endif /* DEM_ROTATION_DASHPOT */
		DP->bUseLegacyRoll = 0;
		prmAddParam(msr->prm,"bUseLegacyRoll",0,&DP->bUseLegacyRoll,
					sizeof(int),"legacyroll","<use legacy rolling friction model (nonequal torques)> = 0");
		DP->bUseContactRadius = 0;
		prmAddParam(msr->prm,"bUseContactRadius",0,&DP->bUseContactRadius,
					sizeof(int),"conrad","<use contact radius instead of reduced radius for rolling friction> = 0");
#ifdef DEM_COHESION
		DP->iCohesionModel = 0; /* 0=microscopic contact, 1=granular bridge */
		prmAddParam(msr->prm,"iCohesionModel",1,&DP->iCohesionModel,
					sizeof(int),"CohesionModel","<CohesionModel> = 0");
		DP->dCohesiveCoeff = 0.0;  /* the cohesive coefficient */
		prmAddParam(msr->prm,"dCohesiveCoeff",2,&DP->dCohesiveCoeff,
					sizeof(double),"CohesiveCoeff","<CohesiveCoeff> = 0.0");
        /* optional particle-variable cohesive coefficient */
        strcpy(achCohesionFile,"");
        prmAddParam(msr->prm,"achCohesionFile",3,achCohesionFile,MAXPATHLEN,
                    "cohesionfile","<name of cohesion data file> = \"\"");
#endif /* DEM_COHESION */
#ifdef DEM_TWOLAYERS
		{
			DP->dInnerOverlapBndry = 0.0; /* InnerOverlapBndry is at particles surfaces by default */
			prmAddParam(msr->prm,"dInnerOverlapBndry",2,&DP->dInnerOverlapBndry,
						sizeof(double),"InnerOverlapBndry","<InnerOverlapBndry> = 0");
			DP->dKnOuter = -1; /* (m)(v_max^2)/(x_max^2) */
			prmAddParam(msr->prm,"dKnOuter",2,&DP->dKnOuter,
						sizeof(double),"knOuter","<knOuter> = -1");
			DP->dKtOuter = -1; /* KtOuter default is (2/7) knOuter */
			prmAddParam(msr->prm,"dKtOuter",2,&DP->dKtOuter,
						sizeof(double),"ktOuter","<ktOuter> = -1");
			DP->dKnInner = -1; /* (m)(v_max^2)/(x_max^2) */
			prmAddParam(msr->prm,"dKnInner",2,&DP->dKnInner,
						sizeof(double),"knInner","<knInner> = -1");
			DP->dKtInner = -1; /* KtInner default is (2/7) knInner */
			prmAddParam(msr->prm,"dKtInner",2,&DP->dKtInner,
						sizeof(double),"ktInner","<ktInner> = -1");
			}
#endif /* DEM_TWOLAYERS */
		}
#endif /* DEM */

#if (INTERNAL_WARNINGS == 0)
	fprintf(stderr,"WARNING: INTERNAL_WARNINGS disabled\n");
#endif

	/*
	 ** Set the box center to (0,0,0) for now!
	 */
	for (j=0;j<3;++j) msr->fCenter[j] = 0.0;
	/*
	 ** Define any "LOCAL" parameters (LCL)
	 */
	msr->lcl.pszDataPath = getenv("PTOOLS_DATA_PATH");
	/*
	 ** Process command line arguments.
	 */
	ret = prmArgProc(msr->prm,argc,argv);
	if (!ret) {
		_msrExit(msr,1);
		}

#if defined(GASOLINE)
	assert(msr->param.bKDK); /* DKD broken at the moment... */
#elif defined(COLLISIONS)
	assert(msr->param.bKDK); /* ditto */
#endif

	if (nDigits < 1 || nDigits > MAXPATHLEN/2) {
		fprintf(stderr,"Invalid number of filename digits.\n");
		_msrExit(msr,1);
		}

	snprintf(msr->param.achDigitMask,MAXPATHLEN,"%%s.%%0%ii",nDigits);

#ifndef BENCHMARK
#if !defined(SPRINGS) && !defined(DEM)
	if ( getenv("PKDGRAV_CHECKPOINT_FDL") == NULL ) {
          fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set\n");
	}
#endif /* !SPRINGS && !DEM */
#endif /* !BENCHMARK */
	/*
	 ** Don't allow periodic BC's for Kepler orbital problems.
	 ** It just doesn't make sense, does it?
	 */
	if (msr->param.bFandG) msr->param.bPeriodic = 0;
	/*
	 ** Make sure that we have some setting for nReplicas if bPeriodic is set.
	 */
	if (msr->param.bPeriodic && !prmSpecified(msr->prm,"nReplicas")) {
		msr->param.nReplicas = 1;
		}
	/*
	 ** Warn that we have a setting for nReplicas if bPeriodic NOT set.
	 */
	if (!msr->param.bPeriodic && msr->param.nReplicas != 0)
		fprintf(stderr,"WARNING: nReplicas set to non-zero value for non-periodic!\n");

	if (!msr->param.achInFile[0] && !msr->param.bRestart) {
		puts("ERROR: no input file specified");
		_msrExit(msr,1);
		}

	if (msr->param.dTheta <= 0) {
		if (msr->param.dTheta == 0 && msr->param.bVWarnings)
			fprintf(stderr,"WARNING: Zero opening angle may cause numerical problems\n");
		else if (msr->param.dTheta < 0) {
			puts("ERROR: Opening angle must be non-negative");
			_msrExit(msr,1);
			}
		}

	msr->nThreads = mdlThreads(mdl);

	if (msr->param.bDoSelfGravity && !msr->param.bDoGravity)
		fprintf(stderr,"WARNING: May need bDoGravity on for bDoSelfGravity to work\n");

	/*
	 ** Always set bCannonical = 1 if bComove == 0
	 */
	if (!msrComove(msr)) {
		if (!msr->param.bCannonical)
			fprintf(stderr,"WARNING: bCannonical reset to 1 for non-comoving (bComove == 0)\n");
		msr->param.bCannonical = 1;
		} 
	/* 
	 * Softening 
	 */
	
	if (msr->param.bPhysicalSoft || msr->param.bVariableSoft) {
	  if (msr->param.bPhysicalSoft && !msrComove(msr)) {
	    fprintf(stderr,"WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
	    msr->param.bPhysicalSoft = 0;
	  }
#ifndef CHANGESOFT
	  puts("ERROR: You must compile with -DCHANGESOFT to use changing softening options");
	  _msrExit(msr,1);
#endif
	  if (msr->param.bVariableSoft && !prmSpecified(msr->prm,"bDoSoftOutput")) msr->param.bDoSoftOutput=1;
  
	  if (msr->param.bPhysicalSoft && msr->param.bVariableSoft) {
	    puts("ERROR: You may only choose one of Physical or Variable softening");
	    _msrExit(msr,1);
	  }
	}

	/*
	 ** Determine the period of the box that we are using.
	 ** Set the new d[xyz]Period parameters which are now used instead
	 ** of a single dPeriod, but we still want to have compatibility
	 ** with the old method of setting dPeriod.
	 */
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dxPeriod")) {
		msr->param.dxPeriod = msr->param.dPeriod;
		}
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dyPeriod")) {
		msr->param.dyPeriod = msr->param.dPeriod;
		}
	if (prmSpecified(msr->prm,"dPeriod") && 
		!prmSpecified(msr->prm,"dzPeriod")) {
		msr->param.dzPeriod = msr->param.dPeriod;
		}
	/*
	 ** Periodic boundary conditions can be disabled along any of the
	 ** x,y,z axes by specifying a period of zero for the given axis.
	 ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
	 ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
	 */
	if (msr->param.dPeriod  == 0) msr->param.dPeriod  = FLOAT_MAXVAL;
	if (msr->param.dxPeriod == 0) msr->param.dxPeriod = FLOAT_MAXVAL;
	if (msr->param.dyPeriod == 0) msr->param.dyPeriod = FLOAT_MAXVAL;
	if (msr->param.dzPeriod == 0) msr->param.dzPeriod = FLOAT_MAXVAL;
#ifdef GASOLINE
	assert(msr->param.duDotLimit <= 0.0);
	if (msr->param.bBulkViscosity) {
		if (!prmSpecified(msr->prm,"dConstAlpha"))
			msr->param.dConstAlpha=0.5;
		if (!prmSpecified(msr->prm,"dConstBeta"))
			msr->param.dConstBeta=0.5;
		}
#ifndef SHOCKTRACK
	if (msr->param.bShockTracker != 0) {
	        fprintf(stderr,"Compile with -DSHOCKTRACK for Shock Tracking.\n");
			assert(0);
	        }
#endif
#endif
	/*
	 ** Determine opening type.
	 */
	msr->iOpenType = 0;
	msr->bOpenSpec = 1;
	if (prmFileSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmFileSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmFileSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmFileSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
	if (prmArgSpecified(msr->prm,"dTheta")) msr->iOpenType = OPEN_JOSH;
	if (prmArgSpecified(msr->prm,"dAbsPartial")) msr->iOpenType = OPEN_ABSPAR;
	if (prmArgSpecified(msr->prm,"dRelPartial")) msr->iOpenType = OPEN_RELPAR;
	if (prmArgSpecified(msr->prm,"dAbsTotal")) msr->iOpenType = OPEN_ABSTOT;
	if (prmArgSpecified(msr->prm,"dRelTotal")) msr->iOpenType = OPEN_RELTOT;
	if (!msr->iOpenType) {
		msr->iOpenType = OPEN_JOSH;
		msr->bOpenSpec = 0;
		}
	switch (msr->iOpenType) {
	case OPEN_JOSH:
		msr->dCrit = msr->param.dTheta;
		if (!prmSpecified(msr->prm,"dTheta2")) 
			msr->param.dTheta2 = msr->param.dTheta;
		break;
	case OPEN_ABSPAR:
		msr->dCrit = msr->param.dAbsPartial;
		break;
	case OPEN_RELPAR:
		msr->dCrit = msr->param.dRelPartial;
		break;
	case OPEN_ABSTOT:
		msr->dCrit = msr->param.dAbsTotal;
		break;
	case OPEN_RELTOT:
		msr->dCrit = msr->param.dRelTotal;
		break;
	default:
		msr->dCrit = msr->param.dTheta;
		if (!prmSpecified(msr->prm,"dTheta2")) 
			msr->param.dTheta2 = msr->param.dTheta;
		}
	/*
	 ** Initialize comove variables.
	 */
	msr->nMaxOuts = 100;
	msr->pdOutTime = malloc(msr->nMaxOuts*sizeof(double));
	assert(msr->pdOutTime != NULL);
	msr->nOuts = 0;

	/*
	 ** Check timestepping.
	 */

	if (msr->param.iMaxRung < 1) {
		msr->param.iMaxRung = 1;
		if (msr->param.bVWarnings)
			fprintf(stderr,"WARNING: iMaxRung set to 1\n");
		}
	if(msr->param.iMaxRung > 29) {
		fprintf(stderr, "WARNING: Cannot have %d rungs, reducing to 29 (the maximum)\n", msr->param.iMaxRung);
		msr->param.iMaxRung = 29;
		}

	if (msr->param.bGravStep && !msr->param.bDoGravity) {
		puts("ERROR: need gravity to use gravity stepping...");
		_msrExit(msr,1);
		}
	if (msr->param.bEpsAccStep || msr->param.bSqrtPhiStep) {
		msr->param.bAccelStep = 1;
		}
#ifdef GASOLINE
	assert(msr->param.bSphStep);
#endif

#ifdef DEM
	msr->param.FO.DP.dDelta = msrDelta(msr); /* adjustable timesteps?? -srs */
	msr->param.FO.DP.dEpsN = msr->param.CP.dEpsN;
	msr->param.FO.DP.dEpsT = msr->param.CP.dEpsT;
	if (msr->param.FO.DP.dKt == -1)
		msr->param.FO.DP.dKt = 0.2857142857142857143*msr->param.FO.DP.dKn; /* (2/7) * kn */
#ifdef DEM_TWOLAYERS /* could add that if kn_outer = kn_inner & kt_outer = kt_inner to just set and use kn, but leaving for now for debugging and baselining purposes */
	if (msr->param.FO.DP.dKnOuter == -1 && msr->param.FO.DP.dKnInner == -1) msr->param.FO.DP.dKnOuter = msr->param.FO.DP.dKnInner = msr->param.FO.DP.dKn; /* copy kn */
	if (msr->param.FO.DP.dKtOuter == -1) msr->param.FO.DP.dKtOuter = 0.2857142857142857143*msr->param.FO.DP.dKnOuter; /* (2/7) * kn_outer */
	if (msr->param.FO.DP.dKtInner == -1) msr->param.FO.DP.dKtInner = 0.2857142857142857143*msr->param.FO.DP.dKnInner; /* (2/7) * kn_inner */
#endif /* DEM_TWOLAYERS */

#endif /* DEM */

#define SECONDSPERYEAR   31557600.
#ifdef GASOLINE
#define msrSetGasModel( iModel ) { \
  if (msr->param.iGasModel == GASMODEL_UNSET) msr->param.iGasModel = iModel; \
  else { \
    fprintf( stderr, "More than one Gas Model specified [%i and %i]", msr->param.iGasModel, iModel ); \
    assert(0); \
  } \
}
	if (msr->param.bGasAdiabatic    ) msrSetGasModel( GASMODEL_ADIABATIC );
	if (msr->param.bGasIsothermal   ) msrSetGasModel( GASMODEL_ISOTHERMAL );
	if (msr->param.bGasCooling      ) msrSetGasModel( GASMODEL_COOLING );

	if (msr->param.iGasModel == GASMODEL_UNSET) {
		msr->param.iGasModel = GASMODEL_ADIABATIC;
		fprintf( stderr, "Defaulting to Adiabatic Gas Model."); 
		}	

	switch(msr->param.iGasModel) {
	case GASMODEL_ADIABATIC:
		msr->param.bGasAdiabatic = 1;	
		break;
	case GASMODEL_ISOTHERMAL:
		msr->param.bGasIsothermal = 1;
		break;
	case GASMODEL_COOLING:
		msr->param.bGasCooling = 1;
		break;
		}

        if( msr->param.bNFWSpheroid ){
            assert (prmSpecified(msr->prm, "dNFWm200") &&
                         prmSpecified(msr->prm, "dNFWr200")&&
                         prmSpecified(msr->prm, "dNFWconc")&&
                         prmSpecified(msr->prm, "dNFWsoft"));
            }
	if(msr->param.iGasModel == GASMODEL_COOLING) {
		assert (prmSpecified(msr->prm, "dMsolUnit") &&
				prmSpecified(msr->prm, "dKpcUnit"));
		}
	/* bolzman constant in cgs */
#define KBOLTZ	1.38e-16
	/* mass of hydrogen atom in grams */
#define MHYDR 1.67e-24
	/* solar mass in grams */
#define MSOLG 1.99e33
	/* G in cgs */
#define GCGS 6.67e-8
	/* kiloparsec in centimeters */
#define KPCCM 3.085678e21
	/* Thompson cross-section (cm^2) */
#define SIGMAT 6.6524e-25
	/* Speed of Light cm/s */
#define LIGHTSPEED 2.9979e10
	/*
	 ** Convert kboltz/mhydrogen to system units, assuming that
	 ** G == 1.
	 */
	if(prmSpecified(msr->prm, "dMsolUnit") &&
	   prmSpecified(msr->prm, "dKpcUnit")) {
		msr->param.dGasConst = msr->param.dKpcUnit*KPCCM*KBOLTZ
			/MHYDR/GCGS/msr->param.dMsolUnit/MSOLG;
		/* code energy per unit mass --> erg per g */
		msr->param.dErgPerGmUnit = GCGS*msr->param.dMsolUnit*MSOLG/(msr->param.dKpcUnit*KPCCM);
		/* code density --> g per cc */
		msr->param.dGmPerCcUnit = (msr->param.dMsolUnit*MSOLG)/pow(msr->param.dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		msr->param.dSecUnit = sqrt(1/(msr->param.dGmPerCcUnit*GCGS));
		/* code comoving density --> g per cc = msr->param.dGmPerCcUnit (1+z)^3 */
		msr->param.dComovingGmPerCcUnit = msr->param.dGmPerCcUnit;
		}

	if(msr->param.bStarForm || msr->param.bFeedBack) {
	    assert (prmSpecified(msr->prm, "dMsolUnit") &&
		    prmSpecified(msr->prm, "dKpcUnit"));

		if (!(msr->param.iGasModel == GASMODEL_COOLING)) {
			fprintf(stderr,"Warning: You are not running a cooling"
					"EOS with starformation\n");
			}

	if (msr->param.bBHSink) {
	    /* For BH sinks -- default to metallicity as sink indicator */
	    if(!prmSpecified(msr->prm, "dSinkMassMin")) msr->param.dSinkMassMin = FLT_MAX;
#ifndef GASOLINE
	    fprintf(stderr, "Gas required for BH Sinks\n");
	    assert(0);
#endif
	    msr->param.bDoSinks = 1;
            /* Units of inverse time -- code units */
	    msr->param.dBHSinkEddFactor = GCGS*4*M_PI*MHYDR/
		(SIGMAT*LIGHTSPEED*msr->param.dBHSinkEddEff)/msr->param.dSecUnit;
	    /* c^2 times efficiency factor (ergs per g) -- code units */
	    msr->param.dBHSinkFeedbackFactor = msr->param.dBHSinkFeedbackEff*msr->param.dBHSinkEddEff*(LIGHTSPEED*LIGHTSPEED)/msr->param.dErgPerGmUnit;
	    }

#ifdef STARFORM
	    assert((msr->param.stfm->dStarEff > 0.0 && 
                msr->param.stfm->dStarEff < 1.0) ||
			   msr->param.stfm->dInitStarMass > 0.0);
	    assert(msr->param.stfm->dMinMassFrac > 0.0 && 
                msr->param.stfm->dMinMassFrac < 1.0);
		if (msr->param.stfm->dInitStarMass > 0) {
/*			if (msr->param.stfm->dMinGasMass <= 0) */
 			  /* Only allow 10% underweight star particles */
				msr->param.stfm->dMinGasMass = 0.9*msr->param.stfm->dInitStarMass;
			}
		else
			assert(msr->param.stfm->dMinGasMass > 0.0);

	    msr->param.stfm->dSecUnit = msr->param.dSecUnit;
	    msr->param.stfm->dGmPerCcUnit = msr->param.dGmPerCcUnit;
	    msr->param.stfm->dGmUnit = msr->param.dMsolUnit*MSOLG;
	    msr->param.stfm->dErgUnit =
		GCGS*pow(msr->param.dMsolUnit*MSOLG, 2.0)
		/(msr->param.dKpcUnit*KPCCM);
            msr->param.dKBoltzUnit = msr->param.dKpcUnit*KPCCM*KBOLTZ
			/GCGS/msr->param.dMsolUnit/MSOLG/msr->param.dMsolUnit/MSOLG;
	    /* convert to system units */
	    msr->param.stfm->dPhysDenMin *= MHYDR/msr->param.stfm->dGmPerCcUnit;
            msr->param.dDeltaStarForm *= SECONDSPERYEAR/msr->param.dSecUnit;
            msr->param.stfm->dDeltaT = msr->param.dDeltaStarForm;

	    msr->param.fb->dSecUnit = msr->param.dSecUnit;
	    msr->param.fb->dGmUnit = msr->param.dMsolUnit*MSOLG;
	    msr->param.fb->dErgPerGmUnit = msr->param.dErgPerGmUnit;
	    msr->param.fb->dInitStarMass = msr->param.stfm->dInitStarMass;
#endif /* STARFORM */
#ifdef SIMPLESF		
		assert(msr->param.SSF_dInitStarMass > 0.0);
#endif
	    }
	    
#endif /* GASOLINE */

#ifdef ROT_FRAME
	puts("ERROR: rotating frame no longer supported.");
	_msrExit(msr,1);
	if (msr->param.bVWarnings && !msr->param.bRotFrame) {
		fprintf(stderr,"WARNING: ROT_FRAME set without bRotFrame\n");
		}
#endif /* ROT_FRAME */

	/*
	 ** Parameter checks and initialization.
	 */

	if (msr->param.bHeliocentric) {
		if (!msr->param.bDoGravity) {
			puts("ERROR: must enable gravity (bDoGravity) to use heliocentric frame.\n");
			_msrExit(msr,1);
			}
#ifdef SIMPLE_GAS_DRAG
		if (msr->param.bSimpleGasDrag) {
			puts("ERROR: cannot combine heliocentric frame with simple gas drag (yet!).\n");
			_msrExit(msr,1);
			}
#endif
		}

#ifdef SPINUP
	/* units converted in msrSpinup() */
	if (msr->param.dSpinupPeriodStart <= 0.0 ||
		msr->param.dSpinupPeriodInter <= 0.0 ||
		msr->param.dSpinupPeriodFinal <= 0.0 ||
		msr->param.nSpinupStepsStart < 0 ||
		msr->param.nSpinupStepsInter < 0 ||
		msr->param.nSpinupStepsFinal < 0 ||
		msr->param.nSpinupStepsHold < 0 ||
		msr->param.nSpinupStepsRestNum < 0.0 ||
		msr->param.nSpinupStepsRestTime < 0.0 ||
		msr->param.dSpinupRadiusInner < 0.0 ||
		msr->param.dSpinupRadiusOuter < 0.0 ||
		msr->param.dSpinupRadiusInner > msr->param.dSpinupRadiusOuter ||
		msr->param.nSpinupCoheTime < 0) {
		puts("ERROR: invalid SPINUP parameter(s).\n");
		_msrExit(msr,1);
		}
#endif /* SPINUP */

#ifdef SPINUP_WITH_AGGS
    puts("WARNING: SPINUP_WITH_AGGS set -- aggregates leaving spinup region may cause errors.\n");
#endif /* SPINUP_WITH_AGGS */
	
#ifdef SLIDING_PATCH
	assert(msr->param.bDoGravity);
	assert(!msr->param.bHeliocentric);
	assert(!msr->param.bFandG);
	if (msr->param.PP.bPatch) {
		PATCH_PARAMS *PP = &msr->param.PP;
		if (!msr->param.bPeriodic) {
			puts("ERROR: must use periodic BCs for patch model");
			_msrExit(msr,1);
			}
		assert(msr->param.dxPeriod > 0.0);
		assert(msr->param.dyPeriod > 0.0);
		assert(msr->param.dzPeriod > 0.0);
		if (msr->param.dyPeriod == FLOAT_MAXVAL) {
			puts("ERROR: must specify positive y period for patch");
			_msrExit(msr,1);
			}
		assert(msr->param.nReplicas > 0);
		if (msr->param.bEwald) {
			puts("ERROR: cannot use Ewald correction in patch model");
			_msrExit(msr,1);
			}
		if (msr->param.dCentMass <= 0.0) {
			puts("ERROR: must specify positive central mass for patch model");
			_msrExit(msr,1);
			}
		if (PP->dOrbDist <= 0.0) {
			puts("ERROR: must specify positive orbital distance for patch model");
			_msrExit(msr,1);
			}
		if (PP->bExtPert) {
			if (PP->dPertOrbDist <= 0.0) {
				puts("ERROR: must specify positive perturber orbital distance");
				_msrExit(msr,1);
				}
			if (PP->dPertMass <= 0.0) {
				puts("ERROR: must specify positive perturber mass");
				_msrExit(msr,1);
				}
			if (PP->dPertMaxZ < 0.0) {
				puts("ERROR: perturber maximum vertical displacement must be non-negative");
				_msrExit(msr,1);
				}
			if (PP->dPertOrbFreqZ < 0.0) {
				puts("ERROR: perturber vertical orbital frequency must be non-negative");
				_msrExit(msr,1);
				}
			}
		if (PP->bRandAzWrap) {
			if (!PP->bNoRandomX) {
				if (PP->nWrapAttempts <= 0) {
					puts("ERROR: number of wrap attempts must be positive");
					_msrExit(msr,1);
					}
				if (PP->iStripOption != STRIP_LEFT_ONLY && PP->iStripOption !=
					STRIP_RIGHT_ONLY && PP->iStripOption != STRIP_BOTH) {
					puts("ERROR: invalid strip option");
					_msrExit(msr,1);
					}
				if (PP->dStripInner < 0.0 || PP->dStripInner >= 0.5) {
					puts("ERROR: strip inner edge must be >= 0.0 and < 0.5");
					_msrExit(msr,1);
					}
				if (PP->dStripOuter <= 0.0 || PP->dStripOuter > 0.5) {
					puts("ERROR: strip outer edge must be > 0.0 and <= 0.5");
					_msrExit(msr,1);
					}
				if (PP->dStripInner >= PP->dStripOuter) {
					puts("ERROR: strip inner edge must be < strip outer edge");
					_msrExit(msr,1);
					}
				}
			if (PP->dVelDispX < 0.0) {
				puts("ERROR: radial velocity dispersion must be non-negative");
				_msrExit(msr,1);
				}
			if (PP->dVelDispY < 0.0) {
				puts("ERROR: azimuthal velocity dispersion must be non-negative");
				_msrExit(msr,1);
				}
			if (PP->dAvgVertAmp < 0.0) {
				puts("ERROR: vertical amplitude must be non-negative");
				_msrExit(msr,1);
				}
			if (PP->dAvgMass <= 0.0) {
				puts("ERROR: average particle mass must be positive");
				_msrExit(msr,1);
				}
			}
		/* load up PATCH_PARAMS struct */
		PP->dCentMass = msr->param.dCentMass;
		PP->dWidth = msr->param.dxPeriod;
		PP->dLength = msr->param.dyPeriod;
		PP->dOrbFreq = sqrt(PP->dCentMass/(PP->dOrbDist*PP->dOrbDist*PP->dOrbDist));
		PP->dDelta = msr->param.dDelta;
		if (msr->param.bVStart)
			(void) printf("Patch: orb freq = %g = %g rad/s <==> P = %g h = %g d\n",
						  PP->dOrbFreq,2*M_PI*PP->dOrbFreq/(365.25*24*3600),
						  365.25*24/PP->dOrbFreq,365.25/PP->dOrbFreq);
		if (PP->bExtPert) {
			PP->dPertOrbFreq = sqrt(PP->dCentMass/(PP->dPertOrbDist*PP->dPertOrbDist*PP->dPertOrbDist));
			if (msr->param.bVStart)
				(void) printf("Patch: perturber orb freq = %g = %g rad/s <==> P = %g h = %g d\n",
							  PP->dPertOrbFreq,2*M_PI*PP->dPertOrbFreq/(365.25*24*3600),
							  365.25*24/PP->dPertOrbFreq,365.25/PP->dPertOrbFreq);
			}
		if (PP->bRandAzWrap) {
			PP->dStripInner *= PP->dWidth;
			PP->dStripOuter *= PP->dWidth;
			if (msr->param.bVStart)
				(void) printf("Patch: strip limits = %g to %g AU\n",PP->dStripInner,PP->dStripOuter);
			}
		}
	else
		fprintf(stderr,"WARNING: SLIDING_PATCH set without bPatch\n");
#endif /* SLIDING_PATCH */

#ifdef COLLISIONS
#ifdef GASOLINE
	puts("ERROR: can't mix COLLISIONS and GASOLINE!");
	_msrExit(msr,1);
#endif
	assert(MAX_NUM_FRAG > 1); /* must allow at least 2 particles */
	if (!msr->param.bKDK) {
		puts("ERROR: must use KDK scheme for collisions");
		_msrExit(msr,1);
		}
#ifdef DEM
#ifndef DEM_FIXED_BALL
	if (msr->param.bVWarnings && msr->param.nSmooth < 2)
		fprintf(stderr,"WARNING: neighbor detection disabled (nSmooth < 2)\n");
#endif /* !DEM_FIXED_BALL */
#else
#ifdef WALLS
	if (msr->param.bVWarnings && msr->param.nSmooth < 1)
		fprintf(stderr,"WARNING: collision detection disabled (nSmooth < 1)\n");
#else
	if (msr->param.nSmooth < 1) {
		puts("ERROR: nSmooth must be positive");
		_msrExit(msr,1);
		}
	if (msr->param.bVWarnings && msr->param.nSmooth < 2)
		fprintf(stderr,"WARNING: collision detection disabled (nSmooth < 2)\n");
#endif /* !WALLS */
#endif /* !DEM */
	assert(msr->param.dCentMass >= 0.0);
	if (msr->param.bFandG) {
#ifdef WALLS
		assert(0);
#endif
		assert(msr->param.bHeliocentric);
		if (!msr->param.bCannonical) {
			puts("ERROR: must use cannonical momentum in FandG collision model");
			_msrExit(msr,1);
			}
		if (msr->param.dSmallStep < 0) {
			puts("ERROR: dSmallStep cannot be negative");
			_msrExit(msr,1);
			}
		if (msr->param.bVWarnings && msr->param.dSmallStep == 0)
			fprintf(stderr,"WARNING: encounter detection disabled (dSmallStep = 0)\n");
		if (msr->param.dSmallStep > msr->param.dDelta) {
			puts("ERROR: inner step must be less than dDelta");
			_msrExit(msr,1);
			}
		}
	if (msr->param.bFindRejects && msr->param.nSmooth < 2) {
		puts("ERROR: nSmooth must be > 1 to find rejects");
		_msrExit(msr,1);
		}
	switch (msr->param.iCollLogOption) {
	case COLL_LOG_NONE:
		break;
	case COLL_LOG_VERBOSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_TXT);
		break;
	case COLL_LOG_TERSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_BIN);
		break;
	default:
		puts("ERROR: Invalid collision log option");
		_msrExit(msr,1);
		}
	if (msr->param.iCollLogOption && msr->param.bVStart)
		printf("Collision log: \"%s\"\n",msr->param.achCollLog);
#ifndef DEM
#ifdef WALLS
	if (msr->param.nSmooth >= 1)
#else
	if (msr->param.nSmooth >= 2)
#endif
	{
		COLLISION_PARAMS *CP = &msr->param.CP;
		if (!(CP->iOutcomes & (MERGE | BOUNCE | FRAG))) {
			puts("ERROR: must specify one of MERGE/BOUNCE/FRAG");
			_msrExit(msr,1);
			}
		if ((CP->iOutcomes & MERGE)) {
			if (CP->iOutcomes != MERGE && !msr->param.bDoSelfGravity) {
				puts("ERROR: need interparticle gravity for conditional merging");
				_msrExit(msr,1);
				}
			if (CP->dDensity < 0.0) {
				puts("ERROR: invalid merge density");
				_msrExit(msr,1);
				}
			if (CP->iDensityAltCol < 0 || CP->iDensityAltCol == PLANETESIMAL) {
				puts("ERROR: invalid alternate density color");
				_msrExit(msr,1);
				}
			if (CP->iDensityAltCol > 0 && CP->dDensityAltVal <= 0.0) {
				puts("ERROR: invalid alternate density");
				_msrExit(msr,1);
				}
			}
		else if (CP->dMergeLimit != 0.0 || CP->dDensity != 0.0 ||
				 CP->iDensityAltCol != 0 || CP->dDensityAltVal != 0.0) /* check against defaults */
				fprintf(stderr,"WARNING: merge parameters ignored (no merging)\n");
		if (CP->iOutcomes & BOUNCE) {
			switch (CP->iEpsNOption) {
			case ConstEps:
				if (CP->dEpsN <= 0.0 || CP->dEpsN > 1.0) {
					puts("ERROR: coef. of restitution must be > 0 and <= 1");
					_msrExit(msr,1);
					}
				break;
			case PowerLaw:
				if (CP->dEpsNCoef <= 0.0) {
					puts("ERROR: invalid EpsN power law coefficient");
					_msrExit(msr,1);
					}
				break;
			case Compacted:
				break;
			case Borderies:
				if (CP->dEpsNVStar < 0.0) {
					puts("ERROR: invalid VStar for Borderies law");
					_msrExit(msr,1);
					}
				break;
			default:
				puts("ERROR: invalid EpsN option");
				_msrExit(msr,1);
				}
			if (CP->iEpsNOption != ConstEps && (CP->dEpsNMin <= 0.0 || CP->dEpsNMin > 1.0)) {
				puts("ERROR: minimum EpsN must be positive and less than 1");
				_msrExit(msr,1);
				}
			if (CP->dSlideLimit < 0.0 || CP->dCollapseLimit < 0.0 || CP->dCrushLimit < 0.0) {
				puts("ERROR: invalid slide, collapse, and/or crush limit");
				_msrExit(msr,1);
				}
			if (CP->dCrushLimit == 0.0)
				CP->dCrushLimit = DBL_MAX;
			if (CP->dSlideLimit > 0.0) {
				if (CP->dSlideEpsN <= 0.0 || CP->dSlideEpsN > 1.0) {
					puts("ERROR: slide limit coef of restitution must be > 0 and <= 1");
					_msrExit(msr,1);
					}
				if (CP->dSlideEpsT < -1.0 || CP->dSlideEpsT > 1.0) {
					puts("ERROR: slide limit coef of surface friction must be >= -1 and <= 1");
					_msrExit(msr,1);
					}
				}
			if (CP->dCrushLimit < DBL_MAX) {
				if (CP->dCrushEpsN <= 0.0 || CP->dCrushEpsN > 1.0) {
					puts("ERROR: crush limit coef of restitution must be > 0 and <= 1");
					_msrExit(msr,1);
					}
				if (CP->dCrushEpsT < -1.0 || CP->dCrushEpsT > 1.0) {
					puts("ERROR: crush limit coef of surface friction must be >= -1 and <= 1");
					_msrExit(msr,1);
					}
				}
			}
		else if (CP->iEpsNOption != ConstEps || CP->dEpsN != 0.0 || CP->dEpsT != 0.0 ||
				 CP->iSlideOption != EscVel || CP->dSlideLimit != 0.0 || CP->dSlideEpsN != 1.0 || CP->dSlideEpsT != 1.0 ||
				 CP->dCollapseLimit != 0.0 || CP->dCollapseEpsN != 1.0 || CP->dCollapseEpsT != 1.0 ||
				 CP->dCrushLimit != DBL_MAX || CP->dCrushEpsN != 1.0 || CP->dCrushEpsT != 1.0) {
			fprintf(stderr,"WARNING: bounce parameters ignored (no bouncing)\n");
			CP->dCollapseLimit = 0.0; /* just in case */
			}
		if (CP->iOutcomes & FRAG) {
#ifndef AGGS
#ifndef SPRINGS
			puts("ERROR: FRAG not available with this compile option");
			_msrExit(msr,1);
#endif
#endif
			}
		else if (CP->dFragLimit != 0.0)
			fprintf(stderr,"WARNING: frag limit ignored\n");
		if (CP->iOutcomes & (MERGE | BOUNCE | FRAG)) {
			switch (CP->iOverlapOption) {
			case OverlapIgnore:
				puts("WARNING: overlap ignore option selected");
				break;
			case OverlapIsError:
			case OverlapBackstep:
			case OverlapAdjPos:
			case OverlapRepel:
			case OverlapMerge:
				break;
			default:
				puts("ERROR: invalid overlap option");
				_msrExit(msr,1);
				}
			if (CP->iOverlapOption == OverlapRepel && !msr->param.bDoSelfGravity) {
				puts("ERROR: repulsion force overlap option requires interparticle forces");
				_msrExit(msr,1);
				}
			if (CP->iOverlapOption == OverlapMerge && !(CP->iOutcomes & MERGE)) {
				puts("ERROR: merge overlap option requires merge outcome enabled");
				_msrExit(msr,1);
				}
			if (CP->iOverlapOption != OverlapIgnore && CP->iOverlapOption != OverlapIsError) {
				(void) printf("WARNING: ");
				switch (CP->iOverlapOption) {
				case OverlapBackstep:
					(void) printf("backstep");
					break;
				case OverlapAdjPos:
					(void) printf("adjust position");
					break;
				case OverlapRepel:
					(void) printf("repulsion force");
					break;
				case OverlapMerge:
					(void) printf("merge");
					break;
				default:
					assert(0); /* this shouldn't happen */
					}
				(void) printf(" overlap option selected.\n");
				puts("***");
				puts("Often an overlap is an indication of a problem, but");
				puts("with this overlap option turned on you will only");
				puts("receive a warning (depending on the values of the");
				puts("warning macros found in pkd.h).  Sometimes, usually");
				puts("when there are many bouncing particles confined to a");
				puts("small space, a \"real\" collapse will occur that can");
				puts("only be avoided be tweaking dDelta, nSmooth, dEpsN,");
				puts("dSlideLimit, or dCollapseLimit, or by selecting an");
				puts("overlap option.  However, in extreme circumstances,");
				puts("it may still not be possible to prevent large-scale");
				puts("unphysical anomalies.  Also, when particle mergers");
				puts("happen close to other particles, an unavoidable");
				puts("overlap may occur that can only be fixed by using an");
				puts("option such as \"merge on overlap\".  Finally, note");
				puts("near misses will also be tolerated.  Simultaneous");
				puts("collisions, however, are still not allowed, because");
				puts("of the risk of entering an infinite loop.");
				if (CP->iOverlapOption != OverlapRepel) {
					puts("NOTE: overlapping particles that are not approaching");
					puts("one another are ignored with this option.");
					}
				puts("***");
				}
			if (CP->iOverlapOption == OverlapBackstep) {
				if (CP->dBackstepLimit < 0.0)
					CP->dBackstepLimit = -CP->dBackstepLimit*msr->param.dDelta;
				else if (CP->dBackstepLimit == 0.0)
					CP->dBackstepLimit = DBL_MAX;
				if (CP->dBackstepLimit > msr->param.dDelta)
					fprintf(stderr,"WARNING: backstep limit exceeds timestep\n");
#ifdef WALLS
				if (!CP->bStrictOverlap) {
					puts("ERROR: recommend turning on strict overlap for WALLS with backstep"); /*DEBUG! why? maybe make this a warning instead?*/
					_msrExit(msr,1);
					}
#endif
				}
			if (CP->iOverlapOption == OverlapAdjPos) {
				if (CP->dAdjPosLimit <= 0.0 || CP->dAdjPosLimit >= 1.0) {
					puts("ERROR: adjust position limit must be > 0 and < 1");
					_msrExit(msr,1);
					}
				}
			if (CP->iOverlapOption == OverlapRepel) {
				/*if (CP->dRepelFac >= 0.0) {
				  puts("ERROR: repulsion factor must be negative");
				  _msrExit(msr,1);
				  }*/
				  /* RP-REWRITE 6/19/09: Changing 'repel' to a linear force */
				if (CP->dRepelFac <= 0.0) {
					puts("ERROR: as of 6/19/09, repulsion factor must be *positive*");
					_msrExit(msr,1);
					}
				else 
					/* k = user_constant * dt^2 */
					CP->dRepelFac /= (msr->param.dDelta * msr->param.dDelta);
				}
			if (CP->dBackstepLimit != 0.0 && CP->iOverlapOption != OverlapBackstep)
				fprintf(stderr,"WARNING: backstep limit ignored\n");
			if (CP->dAdjPosLimit != 0.0 && CP->iOverlapOption != OverlapAdjPos)
				fprintf(stderr,"WARNING: adjust position limit ignored\n");
			if (CP->dRepelFac != 0.0 && CP->iOverlapOption != OverlapRepel)
				fprintf(stderr,"WARNING: repulsion factor ignored\n");
			CP->dDensity *= DEN_CGS_SYS; /* convert: cgs to system units */
			CP->dDensityAltVal *= DEN_CGS_SYS;
			if (CP->iSlideOption == MaxTrv) {
				/*DEBUG note: variable gravity vector ignored... */
				/* (best to use representative value here) */
				double g = sqrt(msr->param.dxUnifGrav*msr->param.dxUnifGrav +
								msr->param.dyUnifGrav*msr->param.dyUnifGrav +
								msr->param.dzUnifGrav*msr->param.dzUnifGrav);
				CP->dSlideLimit = sqrt(2.0*g*CP->dSlideLimit);
				}
			CP->dSlideLimit2 = CP->dSlideLimit*CP->dSlideLimit;
			}
		} /* collisions enabled */
	else if (msr->param.CP.iOverlapOption == OverlapRepel) {
		puts("[RP] ERROR: repel overlap option invoked with collisions disabled!");
		_msrExit(msr,1);
		/* k = user_constant * dt^2 */
		//CP->dRepelFac /= (msr->param.dDelta * msr->param.dDelta); NOTE: to use this here, you must move CP's declaration to this scope!
		//fprintf(stderr,"WARNING: repel overlap option invoked with collisions disabled\n");
		}
#endif /* !DEM */
#endif /* COLLISIONS */

#ifdef AGGS
#ifndef COLLISIONS
	puts("ERROR: For aggregates, collisions must be defined.");
	_msrExit(msr,1);
#endif
	if (!msrKDK(msr)) {
		puts("ERROR: For aggregates, only KDK scheme may be used.");
		_msrExit(msr,1);
		}
	if (!msrDoGravity(msr))
	    puts("WARNING: AGGS accel & torque calcs expect bDoGravity = 1.");
	if (msr->param.CP.iOutcomes & FRAG) {
		if (!(msr->param.CP.iOutcomes & BOUNCE)) {
			puts("ERROR: FRAG requires BOUNCE with AGGS");
			_msrExit(msr,1);
			}
		if (msr->param.CP.iOutcomes & MERGE) {
			/* RP-DEBUG: error if these have the same sign (ie both relative or both absolute) and frag is lower than merge */
			if((msr->param.CP.dFragLimit * msr->param.CP.dMergeLimit) < 0.0)
		        puts("[RP] WARNING: mixed units for mergeLimit and fragLimit!");
			else if(fabs(msr->param.CP.dFragLimit) < fabs(msr->param.CP.dMergeLimit)) {
				puts("ERROR: frag limit must be at least as large in magnitude as merge limit");
				_msrExit(msr,1);
				}
			}
		}
	if (msr->param.CP.iEpsNOption != ConstEps)
		puts("[RP] WARNING: variable dEpsN with AGGS is hacked!");
	if (msr->param.CP.dEpsT != 1.0)
		fprintf(stderr,"WARNING: tangential restitution coefficient ignored for AGGS\n");
	if (msr->param.CP.SP.dTensileCoef < 0.0)
		msr->param.CP.SP.dTensileCoef = DBL_MAX;
	if (msr->param.CP.SP.dTensileExp > 0.0)
		fprintf(stderr,"WARNING: positive tensile strength exponent\n");
	if (msr->param.CP.SP.dTensileCoef == DBL_MAX && msr->param.CP.SP.dTensileExp != 0.0)
		fprintf(stderr,"WARNING: tensile strength exponent ignored\n");
	if (msr->param.CP.SP.dShearCoef < 0.0)
		msr->param.CP.SP.dShearCoef = DBL_MAX;
	if (msr->param.CP.SP.dShearExp > 0.0)
		fprintf(stderr,"WARNING: positive shear strength exponent\n");
	if (msr->param.CP.SP.dShearCoef == DBL_MAX && msr->param.CP.SP.dShearExp != 0.0)
		fprintf(stderr,"WARNING: shear strength exponent ignored\n");
#endif /* AGGS */

#ifdef AGGS_IN_PATCH
	if (!msr->param.bPeriodic) {
		puts("ERROR: must use periodicBCs for AGGS_IN_PATCH");
		_msrExit(msr,1);
		}
	if (!msr->param.PP.bPatch) {
		puts("ERROR: must set bPatch for AGGS_IN_PATCH");
		_msrExit(msr,1);
		}
#endif /* AGGS_IN_PATCH */

#if defined(DEM) && defined(AGGS)
	if (msr->param.FO.DP.dMuR > 0.0 || msr->param.FO.DP.dMuT > 0.0) {
	    puts("Warning: When using aggregates, dMuR and dMuT are forced to zero to prevent unphysical behavior.");
	    }
	puts("WARNING: DEM & AGGS both set -- not fully implemented.");
	puts("         Currently merging and collisionally induced fragmentation");
	puts("         are ignored (because these are only implemented in the");
	puts("         hard-sphere code).");
#endif

#ifdef COLLISIONS
	msrUnifGravGetData(msr,achUnifGravFile);
#endif
#ifdef WALLS
	msrWallsGetData(msr,achWallsFile);
	if (msr->param.bVWarnings && msr->param.CP.WP.nWalls == 0)
		fprintf(stderr,"WARNING: no walls data specified.\n");
#ifdef DEM_WALLS_REACT
	if (!msrDoGravity(msr))
	    puts("WARNING: DEM_WALLS_REACT accel & torque calcs expect bDoGravity = 1.");
#endif /* DEM_WALLS_REACT */
#endif /* WALLS */
#ifdef SPECIAL_PARTICLES
	msrSpecialGetData(msr,achSpecialFile);
#endif
#ifdef DEM_COHESION
    msrCohesionGetData(msr,achCohesionFile);
#endif

#ifdef RUBBLE_ZML

	assert(msr->param.bHeliocentric);
	assert(msr->param.iMaxRung > 1);
	assert(msr->param.CP.dDensity > 0);
	if (msr->param.CP.DB.nDustBins > 0) {
		DUST_BINS_PARAMS *p = &msr->param.CP.DB;
		double ri,ro,a;
		int i;
#ifdef ORIGIN_HISTOGRAM
		int jj;
#endif /* ORIGIN_HISTOGRAM */
		/* checks */
		assert(p->nDustBins <= DUST_BINS_MAX_NUM);
#ifdef ORIGIN_HISTOGRAM
		assert(p->nDustBins == NUM_ORIGIN_BINS);
#endif /* ORIGIN_HISTOGRAM */
		assert(p->iDustBinsApplyInt > 0);
		assert(p->dDustBinsInner < p->dDustBinsOuter);
		assert(p->dDustBinsScaleHeight > 0.0);
		assert(p->dDustBinsInitSigma >= 0.0);
		if (p->dDustBinsInitSigma > 0.0) {
			p->dDustBinsInitSigma *= 10*SQ(AU)/M_SUN; /* from g/cm^2 to system units */
			assert(p->dDustBinsInitAlpha != -2);
			}
		/* initialization */
		p->dDustBinsWidth = (p->dDustBinsOuter - p->dDustBinsInner)/p->nDustBins;
		a = p->dDustBinsInitAlpha + 2; /* from integration of surface mass density */
		for (i=0;i<p->nDustBins;i++) {
			ri = p->dDustBinsInner + i*p->dDustBinsWidth;
			ro = ri + p->dDustBinsWidth;
			msr->aDustBins[i].dMass = TWO_PI*p->dDustBinsInitSigma*(pow(ro,a) - pow(ri,a))/a;
			printf("Mass in DustBin = %e\n", msr->aDustBins[i].dMass);
			msr->aDustBins[i].dVolume = M_PI*(SQ(ro) - SQ(ri))*p->dDustBinsScaleHeight;
#ifdef ORIGIN_HISTOGRAM
			for (jj=0;jj<NUM_ORIGIN_BINS;jj++)
				msr->aDustBins[i].origin_bins[jj] = 0.0;
			msr->aDustBins[i].origin_bins[i] = 1.0; /* dust originates in the bin it's in -- THIS REQUIRES nDustBins == NUM_ORIGINS_BINS! */
#endif /* ORIGIN_HISTOGRAM */
			}
		}
	msr->DustBinsTrash.dMass = 0.0; /* no dust in the overflow variable */
#ifdef ORIGIN_HISTOGRAM
	{
		int i;
		for (i=0;i<NUM_ORIGIN_BINS;i++)
			msr->DustBinsTrash.origin_bins[i] = 0.0;
		}
#endif /* ORIGIN_HISTOGRAM */
	msr->re.nEvents = 0; /* no events for rubble clocks yet */	

#elif defined(COLLMOD_ZML) /* TIDY -- This section of code is almost identical to that above may want to tidy this up */

	//assert(msr->param.bHeliocentric);
	assert(msr->param.iMaxRung > 1);
	if (msr->param.CP.DB.nDustBins > 0) {
		DUST_BINS_PARAMS *p = &msr->param.CP.DB;
		double ri,ro,a;
		int i;
#ifdef ORIGIN_HISTOGRAM
		int jj;
#endif /* ORIGIN_HISTOGRAM */
		/* checks */
		assert(p->nDustBins <= DUST_BINS_MAX_NUM);
#ifdef ORIGIN_HISTOGRAM
		assert(p->nDustBins == NUM_ORIGIN_BINS);
#endif /* ORIGIN_HISTOGRAM */
		assert(p->iDustBinsApplyInt > 0);
		assert(p->dDustBinsInner < p->dDustBinsOuter);
		assert(p->dDustBinsScaleHeight > 0.0);
		assert(p->dDustBinsInitSigma >= 0.0);
		if (p->dDustBinsInitSigma > 0.0) {
			p->dDustBinsInitSigma *= 10*SQ(AU)/M_SUN; /* from g/cm^2 to system units */
			assert(p->dDustBinsInitAlpha != -2);
			}
		/* initialization */
		p->dDustBinsWidth = (p->dDustBinsOuter - p->dDustBinsInner)/p->nDustBins;
		a = p->dDustBinsInitAlpha + 2; /* from integration of surface mass density */
		for (i=0;i<p->nDustBins;i++) {
			ri = p->dDustBinsInner + i*p->dDustBinsWidth;
			ro = ri + p->dDustBinsWidth;
			msr->aDustBins[i].dMass = TWO_PI*p->dDustBinsInitSigma*(pow(ro,a) - pow(ri,a))/a;
			printf("Mass in DustBin = %e\n", msr->aDustBins[i].dMass);
			msr->aDustBins[i].dVolume = M_PI*(SQ(ro) - SQ(ri))*p->dDustBinsScaleHeight;
#ifdef ORIGIN_HISTOGRAM
			for (jj=0;jj<NUM_ORIGIN_BINS;jj++)
				msr->aDustBins[i].origin_bins[jj] = 0.0;
			msr->aDustBins[i].origin_bins[i] = 1.0; /* dust originates in the bin it's in -- THIS REQUIRES nDustBins == NUM_ORIGINS_BINS! */
#endif /* ORIGIN_HISTOGRAM */
			}
		}
	msr->DustBinsTrash.dMass = 0.0; /* no dust in the overflow variable */
#ifdef ORIGIN_HISTOGRAM
	{
		int i;
		for (i=0;i<NUM_ORIGIN_BINS;i++)
			msr->DustBinsTrash.origin_bins[i] = 0.0;
		}
#endif /* ORIGIN_HISTOGRAM */

#endif /* RUBBLE_ZML, COLLMOD_ZML */

#ifdef CHARGE
#if !defined(DEM) || !defined(WALLS)
	puts("ERROR: CHARGE must be compiled with DEM and WALLS currently.");
	_msrExit(msr,1);
#endif /* !DEM || !WALLS */
#endif /* CHARGE */

#ifdef DEM_TIDAL_LOCAL
	msrDEMTidalReadData(msr);
#endif /* DEM_TIDAL_LOCAL */

	/* check force override options */

	switch (msr->param.FO.iForceOverrideOption) {
	case FO_NONE:
#ifdef SPRINGS
		puts("ERROR: code compiled with SPRINGS but force override option not set.");
		_msrExit(msr,1);
#endif
#ifdef DEM
		puts("ERROR: code compiled with DEM but force override option not set.");
		_msrExit(msr,1);
#endif
		break;
	case FO_VANDERWAALS:
		puts("ERROR: van der Walls override option not implemented.");
		_msrExit(msr,1); /*FALLTHRU*/
	case FO_STRENGTH:
#ifdef SPRINGS
		if (msr->param.iCheckInterval > 0)
			msrRestartsNotSupportedMessage(msr);			
		if (msr->param.CP.iOutcomes & MERGE) {
			puts("ERROR: merging not supported with SPRINGS");
			_msrExit(msr,1);
			}
		if (!(msr->param.CP.iOutcomes & BOUNCE)) {
			puts("ERROR: bouncing must be used with SPRINGS");
			_msrExit(msr,1);
			}
#ifdef AGGS
		puts("WARNING: SPRINGS: material strength force override option invoked with AGGS.");
#endif
		{
			SPRING_PARAMS *SP = &msr->param.FO.SP;
			if (SP->dMeanYoungsModulus < 0.0) {
				puts("ERROR: mean Young's modulus must be positive.");
				_msrExit(msr,1);
				}
			if (SP->dMeanStrainLimit < 0.0) {
				puts("ERROR: mean strain limit must be positive.");
				_msrExit(msr,1);
				}
			if (SP->dMeanStrainLimit > SP->dMaxStrainLimit) {
				puts("ERROR: mean strain limit must not exceed the maximum allowed strain limit.");
				_msrExit(msr,1);
			}
			if (SP->dYoungsStdDev < 0.0) {
		        puts("WARNING: Young's modulus standard deviation negative... continuing...");
				SP->dYoungsStdDev *= -1.;
		        }
			if (SP->dStrainStdDev < 0.0) {
		        puts("WARNING: strain limit standard deviation negative... continuing...");
				SP->dStrainStdDev *= -1.;
		        }
			if ((SP->dYoungsStdDev) && !(SP->dMeanYoungsModulus))
		        puts("WARNING: dMeanYoungsModulus set to zero, but dYoungsStdDev is nonzero.");
			if ((SP->dStrainStdDev) && !(SP->dMeanStrainLimit))
		        puts("WARNING: dMeanStrainLimit set to zero, but dStrainStdDev is nonzero.");
			if (SP->dMaxStrainLimit > 1.)
		        puts("WARNING: strain limits could be large ( >1.0 ).");
			if (SP->dLinkageLength <= 0.0) {
				puts("ERROR: linkage length must be positive.");
				_msrExit(msr,1);
				}
			if (SP->dLinkageLength <= 2.0)
				puts("WARNING: linkage length <= 2.0 only particles in contact or overlapping will form links.");
			if (SP->dLinkageLength > 2.5)
				puts("WARNING: linkage length > 2.5 may be unphysical for HCP rubble pile");
			if (SP->dZeroStrainLength < 0.0) {
				puts("ERROR: zero strength length must be non-negative.");
				_msrExit(msr,1);
				}
			if (SP->dPorosityFactor <= 0.0 || SP->dPorosityFactor > 1.0) {
				puts("ERROR: porosity factor must be > 0 and <= 1.");
				_msrExit(msr,1);
				}
			if (SP->dDamp < 0.0)
				puts("WARNING: UNSTABLE! spring damping negative... oscillations will grow.");
			if (SP->dDamp >= 3.0)
		        puts("WARNING: spring very overdamped... may want to ensure damping acceleration is reasonable for timestep.");
			if (msr->param.nSmooth > MAX_NUM_SPRINGS_PER_PARTICLE) {
				puts("ERROR: nSmooth cannot exceed MAX_NUM_SPRINGS_PER_PARTICLE.");
				_msrExit(msr,1);
				}
			{
				/* the following few estimates approximate zero strain length = 2 R_eff, dYoungsStdDev = dStrainStdDev = 0 */

		        if (SP->dMeanYoungsModulus && SP->dMeanStrainLimit) {

			        const double dConvert = (1.49597870e11/1.9891e30*(365.25*24*3600)*(365.25*24*3600)/(2.0*M_PI)/(2.0*M_PI)); /* Pascals --> pkdgrav units (multiply) */
					const double dSQRT_kg_m = 3646245620.668213553; /* sqrt[kg/m] --> sqrt[M_sun/AU], mks --> pkdgrav units (multiply) */

					double YM = SP->dMeanYoungsModulus;
					double MaxStrain = SP->dMeanStrainLimit;
					double MaxMax = SP->dMaxStrainLimit;
					double L = YM*MaxStrain;
					double x = SP->dPorosityFactor;
					double dDamp = SP->dDamp;

					double k = 0.5*((YM*L)/(YM + L))*dConvert*M_PI*x;

					(void) printf("SPRINGS: recommended maximum dDelta = %g sqrt([m/kg]/[R_eff/m])\n",0.02*M_PI/sqrt(k)*dSQRT_kg_m); /* may want to check distribution of springs for minimum */

					if (dDamp) {
				        double k2 = sqrt(2.*YM*M_PI*x*dConvert)*dDamp;
						(void) printf("SPRINGS: damping unstable for dDelta > %g sqrt([m/kg]/[R_eff/m])  [for strictly equal mass particles, safe to add factor of sqrt(2)]\n",dSQRT_kg_m/k2); /* may want to check distribution of springs for minimum */
						}

					(void) printf("SPRINGS: cutoff distance ~ %g [R_eff/m] \n",2.0*(1 + L/YM)*149597892000.);

					if (SP->dYoungsStdDev)
				        (void) printf("SPRINGS: Range in Young's modulus [Pa] --> (0 < Y < %e)   with a mean value of %e and a standard deviation of %e.\n",YM+YM,YM,SP->dYoungsStdDev);
					if (SP->dStrainStdDev && SP->dStrainStdDev <= 0.5*MaxMax)
				        (void) printf("SPRINGS: Range in strain limit --> (0 < max_strain < %g)   with a mean value of %g and a standard deviation of %g.\n",MaxStrain+MaxStrain,MaxStrain,SP->dStrainStdDev); /* no StdDev for REFORM_SPRINGS */
					if (SP->dStrainStdDev > 0.5*MaxMax)
				        (void) printf("SPRINGS: Range in strain limit --> (%g <= max_strain <= %g)   with a mean value of %g and a standard deviation of %g.\n",MaxStrain+MaxStrain-MaxMax,MaxMax,MaxStrain,SP->dStrainStdDev); /* no StdDev for REFORM_SPRINGS */
					}
				}
			}
#endif /* SPRINGS */
#ifdef DEM
		if (msr->param.iCheckInterval > 0)
			msrRestartsNotSupportedMessage(msr);			
		if (msr->param.bVWarnings && (msr->param.CP.iOutcomes & MERGE || msr->param.CP.iOutcomes & BOUNCE || msr->param.CP.iOutcomes & FRAG))
			fprintf(stderr,"WARNING: iOutcomes ignored with soft sphere DEM.\n");
		{
			DEM_PARAMS *DP = &msr->param.FO.DP;
			if (DP->dKn < 0.0) {
				puts("ERROR: DEM normal spring constant must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dKn == 0.0)
				puts("WARNING: DEM normal spring constant set to zero.");
			if (DP->dKt < 0.0 && DP->dKt != -1.0 /* special case */) {
				puts("ERROR: DEM tangential spring constant must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dKt == 0.0)
				puts("WARNING: DEM tangential spring constant set to zero.");
			if (DP->dMuS < 0.0) {
				puts("ERROR: DEM static friction coefficient must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dMuR < 0.0) {
				puts("ERROR: DEM rolling friction coefficient must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dMuT < 0.0) {
				puts("ERROR: DEM twisting friction coefficient must not be negative.");
				_msrExit(msr,1);
				}
			if (msr->param.CP.dEpsT < 1.0 && DP->dMuS == 0.0)
				puts("WARNING: currently dEpsT has no effect if dMuS is zero.");
			if (DP->dAccCrit < 0.0) {
				puts("ERROR: DEM critical acceleration must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dTangentialSpringDrag < 0.0 || DP->dTangentialSpringDrag > 1.) {
				puts("ERROR: DEM tangential spring drag coefficient must be >= 0 and <= 1.");
				_msrExit(msr,1);
				}
			if (DP->dTangentialSpringDrag != 0.0 && DP->dTangentialSpringDrag != 1.)
				puts("WARNING: DEM tangential spring drag coefficient not 0 or 1...effect may be timestep dependent!");
			if (DP->dMinorFrac < 0.0) {
				puts("ERROR: Overlap minor warning fraction must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dMajorFrac < DP->dMinorFrac) {
				puts("ERROR: Overlap major warning fraction must not be less than minor warning fraction.");
				_msrExit(msr,1);
				}
			if (DP->dErrorFrac < DP->dMajorFrac) {
				puts("ERROR: Overlap error fraction must not be less than major warning fraction.");
				_msrExit(msr,1);
				}
			if (DP->iDEMStatsInterval < 0) {
				puts("ERROR: DEM stats output interval must not be negative.");
				_msrExit(msr,1);
				}
#ifdef DEM_FIXED_BALL
			if (DP->iFixedBallOption != 0 && DP->iFixedBallOption != 1 && DP->iFixedBallOption != 2) {
				puts("ERROR: DEM Fixed Ball Option needs to be set to 0, 1, or 2");
				_msrExit(msr,1);
				}
#endif /* DEM_FIXED_BALL */
#ifdef DEM_ROTATION_DASHPOT
			if (DP->dBetaS < 0.0) {
				puts("ERROR: DEM shape parameter must not be negative.");
				_msrExit(msr,1);
				}
#endif /* DEM_ROTATION_DASHPOT */
#ifdef DEM_COHESION
			if (DP->iCohesionModel != 0 && DP->iCohesionModel != 1) {
				puts("ERROR: iCohesionModel must take a value of 0 or 1.");
				_msrExit(msr,1);
				}
			if (DP->dCohesiveCoeff < 0.0) {
				puts("ERROR: DEM cohesive coefficient must not be negative.");
				_msrExit(msr,1);
				}
#endif /* DEM_COHESION */
#ifdef DEM_TWOLAYERS
			if (DP->dInnerOverlapBndry < 0.0) {
				puts("ERROR: DEM inner overlap boundary must not be negative.");
				_msrExit(msr,1);
				}
			if (DP->dKnOuter <= 0.0 || DP->dKnInner <= 0.0 || DP->dKtOuter <= 0.0 || DP->dKtInner <= 0.0) {
				puts("ERROR: At least one two-layer spring component not positive.");
				_msrExit(msr,1);
				}
			if (DP->dKnOuter >= DP->dKnInner)
				puts("WARNING: DEM normal outer spring constant equal to or stiffer than inner.");
#endif /* DEM_TWOLAYERS */
			}
#endif /* DEM */
#if !defined(SPRINGS) && !defined(DEM)
		puts("ERROR: Material strength force override not supported with this compile option.");
		_msrExit(msr,1);
#endif /* !SPRINGS && !DEM */
		break;
	default:
		puts("ERROR: Unrecognized force override option.");
		_msrExit(msr,1);
		}

#ifdef GR_DRAG
	msr->dMergerMassLost = 0.0; /* initialization */
#endif

#ifdef WALLS
	msr->dDeathWallsMassLost = 0.0;
#endif

	/* end of options checking */

	pstInitialize(&msr->pst,msr->mdl,&msr->lcl);

	pstAddServices(msr->pst,msr->mdl);
	/*
	 ** Create the processor subset tree.
	 */
	for (id=1;id<msr->nThreads;++id) {
		if (msr->param.bVDetails) printf("Adding %d to the pst\n",id);
		inAdd.id = id;
		pstSetAdd(msr->pst,&inAdd,sizeof(inAdd),NULL,NULL);
		}
	if (msr->param.bVDetails) printf("\n");

	/*
	** Make an easy-to-grep-for strings for the number of
	** processors we are running on.
	*/
#ifdef GASOLINE

#ifdef BENCHMARK
	printf("GASOLINE (BENCHMARK) running on %d processor",msr->nThreads);
#else
	printf("GASOLINE running on %d processor",msr->nThreads);
#endif

#else

#ifdef BENCHMARK
	printf("PKDGRAV (BENCHMARK) running on %d processors",msr->nThreads);
#else
	printf("PKDGRAV running on %d processor",msr->nThreads);
#endif

#endif /* GASOLINE */
	if (msr->nThreads>1) printf("s");
	printf(".\n");

	/*
	 ** Levelize the PST.
	 */
	inLvl.iLvl = 0;
	pstLevelize(msr->pst,&inLvl,sizeof(inLvl),NULL,NULL);
	/*
	 ** Create the processor mapping array for the one-node output
	 ** routines.
	 */
	msr->pMap = malloc(msr->nThreads*sizeof(int));
	assert(msr->pMap != NULL);
	inGM.nStart = 0;
	pstGetMap(msr->pst,&inGM,sizeof(inGM),msr->pMap,NULL);
	/*
	 ** Initialize tree type to none.
	 */
	msr->iTreeType = MSR_TREE_NONE;
	msr->iCurrMaxRung = 0;
	/*
	 ** Mark the Domain Decompositon as not done
	 */
	msr->bDoneDomainDecomp = 0;
	msr->iLastRungDomainDecomp = 0;
	msr->nRung = (int *) malloc( (msr->param.iMaxRung+1)*sizeof(int) );
	/*
	** Seed the random number generator on each processor.
	**
	** NOTE: although the seed chosen for each processor is displayed
	** on stdout, there is currently no way to run the code with a
	** pre-determined seed.  In an emergency, pstRandSeedGenerator()
	** in pst.c can be hacked.
	*/
	pstRandSeedGenerator(msr->pst,NULL,0,NULL,NULL);
	/* handy reminder */
	if (msr->param.bVWarnings && !msr->param.bUseWallClock)
		puts("NOTE: Times reported as \"Wallclock\" are in fact CPU times.");
	}

double msrTime(MSR msr) {
	if (msr->param.bUseWallClock) {
		struct timeval tv;
		gettimeofday(&tv,NULL);
		return tv.tv_sec + tv.tv_usec*1.0e-6;
		}
	else {
		struct rusage ru;
		getrusage(RUSAGE_SELF,&ru);
		return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec*1.0e-6;
		}
	}

void msrLogParams(MSR msr,FILE *fp)
{
	double z, testDelta;
	int i;

#ifdef __DATE__
#ifdef __TIME__
	fprintf(fp,"# Code compiled: %s %s\n",__DATE__,__TIME__);
#endif
#endif
	fprintf(fp,"# Preprocessor macros:");
#ifdef GASOLINE
	fprintf(fp," GASOLINE");
#endif
#ifdef STARFORM
	fprintf(fp," STARFORM");
#endif
#ifdef KROUPA
	fprintf(fp," KROUPA");
#endif
#ifdef SIMPLESF
	fprintf(fp," SIMPLESF");
#endif
#ifdef LARGEFBALL
	fprintf(fp," LARGEFBALL");
#endif
#ifdef SHOCKTRACK
	fprintf(fp," SHOCKTRACK");
#endif
#ifdef PEAKEDKERNEL
	fprintf(fp," PEAKEDKERNEL");
#endif
#ifdef CHANGESOFT
 	fprintf(fp," CHANGESOFT");
#endif
#ifdef NOCOOLING
 	fprintf(fp," NOCOOLING");
#endif
#ifdef COOLING_COSMO
 	fprintf(fp," COOLING_COSMO");
#endif
#ifdef COOLING_PLANET
 	fprintf(fp," COOLING_PLANET");
#endif
#ifdef COOLING_BATE
 	fprintf(fp," COOLING_BATE");
#endif
#ifdef COOLING_DISK
 	fprintf(fp," COOLING_DISK");
#endif
#ifdef GLASS
	fprintf(fp," GLASS");
#endif
#ifdef HSHRINK
	fprintf(fp," HSHRINK");
#endif
#ifdef DEBUG
	fprintf(fp," DEBUG");
#endif
#ifdef ALTSPH
	fprintf(fp," ALTSPH");
#endif
#ifdef SUPERCOOL
	fprintf(fp," SUPERCOOL");
#endif
#ifdef ROT_FRAME
	fprintf(fp," ROT_FRAME");
#endif
/* following entries should match Makefile.in order */
#ifdef COLLISIONS
	fprintf(fp," COLLISIONS");
#endif
#ifdef WALLS
	fprintf(fp," WALLS");
#endif
#ifdef AGGS
	fprintf(fp," AGGS");
#endif
#ifdef FAST_AGGS
	fprintf(fp," FAST_AGGS");
#endif
#ifdef SPRINGS
	fprintf(fp," SPRINGS");
#endif
#ifdef SLIDING_PATCH
	fprintf(fp," SLIDING_PATCH");
#endif
#ifdef DEM
	fprintf(fp," DEM");
#endif
#ifdef DEM_FIXED_BALL
	fprintf(fp," DEM_FIXED_BALL");
#endif
#ifdef DEM_ROTATION_DASHPOT
	fprintf(fp," DEM_ROTATION_DASHPOT");
#endif
#ifdef DEM_TRACK_ORIENT
	fprintf(fp," DEM_TRACK_ORIENT");
#endif
#ifdef DEM_COHESION
	fprintf(fp," DEM_COHESION");
#endif
#ifdef DEM_TWOLAYERS
	fprintf(fp," DEM_TWOLAYERS");
#endif
#ifdef DEM_WALLS_REACT
	fprintf(fp," DEM_WALLS_REACT");
#endif
#ifdef DEM_OVERLAP_OUTPUT
	fprintf(fp," DEM_OVERLAP_OUTPUT");
#endif
#ifdef DEM_DIAG
	fprintf(fp," DEM_DIAG");
#endif
#ifdef SPINUP
	fprintf(fp," SPINUP");
#endif
#ifdef SEISMIC_WAVE
	fprintf(fp," SEISMIC_WAVE");
#endif
#ifdef RUBBLE_ZML
	fprintf(fp," RUBBLE_ZML");
#elif defined(COLLMOD_ZML)
	fprintf(fp," COLLMOD_ZML");
#endif
#ifdef ORIGIN_HISTOGRAM
	fprintf(fp," ORIGIN_HISTOGRAM");
#endif
#ifdef NUMREC
	fprintf(fp," NUMREC");
#endif
#ifdef GSL
	fprintf(fp," GSL");
#endif
#ifdef SSIO_USE_MPI
	fprintf(fp," SSIO_USE_MPI");
#endif
/* end of Makefile.in order */
/* deprecated? */
#ifdef CHARGE
	fprintf(fp," CHARGE");
#endif
#ifdef SPECIAL_PARTICLES
	fprintf(fp," SPECIAL_PARTICLES");
#endif
#ifdef SIMPLE_GAS_DRAG
	fprintf(fp," SIMPLE_GAS_DRAG");
#endif
#ifdef GR_DRAG
	fprintf(fp," GR_DRAG");
#endif
/* end of deprecated? */
#ifdef _REENTRANT
	fprintf(fp," _REENTRANT");
#endif
#ifdef CRAY_T3D
	fprintf(fp," CRAY_T3D");
#endif
#ifdef PRES_MONAGHAN
	fprintf(fp," PRES_MONAGHAN");
#endif 
#ifdef PRES_HK
	fprintf(fp," PRES_HK");
#endif 
	{
	time_t timep;

	(void) time(&timep);
	fprintf(fp,"\n# Run started: %s",ctime(&timep));
	}
	{
	char ach[MAXPATHLEN];
	fprintf(fp,"# Master host: ");
	if (gethostname(ach,MAXPATHLEN))
		fprintf(fp,"unknown");
	else
		fprintf(fp,"%s",ach);
	fprintf(fp,"\n# Current working directory: ");
	if (getcwd(ach,MAXPATHLEN) == NULL)
		fprintf(fp,"unknown");
	else
		fprintf(fp,"%s",ach);
	}
	fprintf(fp,"\n# N: %d",msr->N);
	fprintf(fp," nThreads: %d",msr->param.nThreads);
	fprintf(fp," bDiag: %d",msr->param.bDiag);
	fprintf(fp," Verbosity flags: (%d,%d,%d,%d,%d)",msr->param.bVWarnings,
			msr->param.bVStart,msr->param.bVStep,msr->param.bVRungStat,
			msr->param.bVDetails);
	fprintf(fp,"\n# bPeriodic: %d",msr->param.bPeriodic);
	fprintf(fp," bRestart: %d",msr->param.bRestart);
	fprintf(fp," bComove: %d",msr->param.csm->bComove);
	fprintf(fp,"\n# bParaRead: %d",msr->param.bParaRead);
	fprintf(fp," bParaWrite: %d",msr->param.bParaWrite);
	fprintf(fp," bUseWallClock: %d",msr->param.bUseWallClock);
	fprintf(fp," bCannonical: %d",msr->param.bCannonical);
	fprintf(fp," bStandard: %d",msr->param.bStandard);
	fprintf(fp,"\n# bKDK: %d",msr->param.bKDK);
	fprintf(fp," nBucket: %d",msr->param.nBucket);
	fprintf(fp," iOutInterval(%d,%d): %d",msr->param.iBinaryOutput,msr->param.bPackedVector,msr->param.iOutInterval);
	fprintf(fp," dDumpFrameStep: %g",msr->param.dDumpFrameStep);
	fprintf(fp," dDumpFrameTime: %g",msr->param.dDumpFrameTime);
	fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
	fprintf(fp," bLogTiming: %d (%d%d%d%d)",msr->param.bLogTiming,msr->param.bLogTimingSubStep,msr->param.bLogTimingStep,msr->param.bLogTimingSubStepTot,msr->param.bLogTimingStepTot);
	fprintf(fp,"\n# iCheckInterval: %d",msr->param.iCheckInterval);
	fprintf(fp," iOrder: %d",msr->param.iOrder);
	fprintf(fp," iEwOrder: %d",msr->param.iEwOrder);
	fprintf(fp," nReplicas: %d",msr->param.nReplicas);
	fprintf(fp,"\n# dEwCut: %f",msr->param.dEwCut);
	fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
	fprintf(fp,"\n# iStartStep: %d",msr->param.iStartStep);
	fprintf(fp," nSteps: %d",msr->param.nSteps);
	fprintf(fp," nSmooth: %d",msr->param.nSmooth);
	fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
	if (prmSpecified(msr->prm,"dSoft"))
		fprintf(fp," dSoft: %g",msr->param.dSoft);
	else
		fprintf(fp," dSoft: input");
	fprintf(fp,"\n# bPhysicalSoft: %d",msr->param.bPhysicalSoft);
	fprintf(fp," bVariableSoft: %d (%d %d %d)",msr->param.bVariableSoft,msr->param.bVariableSoftStar,msr->param.bVariableSoftGas,msr->param.bVariableSoftDark);
	fprintf(fp," nSoftNbr: %d",msr->param.nSoftNbr);
	fprintf(fp," bSoftByType: %d",msr->param.bSoftByType);
	fprintf(fp," bSoftMaxMul: %d",msr->param.bSoftMaxMul);
	fprintf(fp," dSoftMax: %g",msr->param.dSoftMax);
	fprintf(fp," bDoSoftOutput: %d",msr->param.bDoSoftOutput);
	fprintf(fp,"\n# dDelta: %g",msr->param.dDelta);
	fprintf(fp," dEta: %g",msr->param.dEta);
	fprintf(fp," dEtaDeltaAccel: %g",msr->param.dEtaDeltaAccel);
	fprintf(fp," dEtaCourant: %g",msr->param.dEtaCourant);
	fprintf(fp," iMaxRung: %d",msr->param.iMaxRung);
	fprintf(fp,"\n# bGravStep: %d",msr->param.bGravStep);
	fprintf(fp," bEpsAccStep: %d",msr->param.bEpsAccStep);
	fprintf(fp," bSqrtPhiStep: %d",msr->param.bSqrtPhiStep);
	fprintf(fp," bDensityStep: %d",msr->param.bDensityStep);
	fprintf(fp," bDeltaAccelStep: %d",msr->param.bDeltaAccelStep);
	fprintf(fp," (gt): %d",msr->param.bDeltaAccelStepGasTree);
	fprintf(fp," nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp," bNonSymp: %d",msr->param.bNonSymp);
	fprintf(fp,"\n# bDoGravity: %d",msr->param.bDoGravity);
	fprintf(fp," bDoSelfGravity: %d",msr->param.bDoSelfGravity);
	fprintf(fp," bFandG: %d",msr->param.bFandG);
	fprintf(fp," bHeliocentric: %d",msr->param.bHeliocentric);
	fprintf(fp," dCentMass: %g",msr->param.dCentMass);
	fprintf(fp," dSunSoft: %g",msr->param.dSunSoft);
	fprintf(fp,"\n# bLogHalo: %d",msr->param.bLogHalo );
	fprintf(fp," bHernquistSpheroid: %d",msr->param.bHernquistSpheroid );
	fprintf(fp," bNFWSpheroid: %d",msr->param.bNFWSpheroid );
        if( msr->param.bNFWSpheroid ){
            fprintf(fp," dNFWm200: %g",msr->param.dNFWm200 );
            fprintf(fp," dNFWr200: %g",msr->param.dNFWr200 );
            fprintf(fp," dNFWsoft: %g",msr->param.dNFWsoft );
            fprintf(fp," dNFWconc: %g",msr->param.dNFWconc );
            }
	fprintf(fp," bHomogSpheroid: %d",msr->param.bHomogSpheroid );
	fprintf(fp," bBodyForce: %d",msr->param.bBodyForce );
	fprintf(fp," bMiyamotoDisk: %d",msr->param.bMiyamotoDisk );
	fprintf(fp," bTimeVarying: %d",msr->param.bTimeVarying );
	fprintf(fp,"\n# bRotatingBar: %d",msr->param.bRotatingBar);
	if (msr->param.bRotatingBar != 0) {
	    rotbarLogParams( msr->param.rotbar, fp );
	    }
#ifdef ROT_FRAME
	fprintf(fp,"\n# bRotFrame: %d",msr->param.bRotFrame);
	fprintf(fp," dOmega: %g",msr->param.dOmega);
	fprintf(fp," dOmegaDot: %g",msr->param.dOmegaDot);
#endif /* ROT_FRAME */
	fprintf(fp,"\n# bDoSinks: %d",msr->param.bDoSinks );
	fprintf(fp," bBHSink: %d",msr->param.bBHSink );
	fprintf(fp," dBHSinkEddEff: %g",msr->param.dBHSinkEddEff);
	fprintf(fp," dBHSinkFeedbackEff: %g",msr->param.dBHSinkFeedbackEff);
	fprintf(fp," dBHSinkAlpha: %g",msr->param.dBHSinkAlpha);
	fprintf(fp," bDoSinksAtStart: %d",msr->param.bDoSinksAtStart );
	fprintf(fp," bSinksThermal: %d",msr->param.bSinkThermal );
	fprintf(fp," dSinkRadius: %g",msr->param.dSinkRadius);
	fprintf(fp," dSinkBoundOrbitRadius: %g",msr->param.dSinkBoundOrbitRadius);
	fprintf(fp," dSinkMassMin: %g",msr->param.dSinkMassMin);
	fprintf(fp,"\n# dFracNoDomainDecomp: %g",msr->param.dFracNoDomainDecomp);
	fprintf(fp," dFracNoDomainDimChoice: %g",msr->param.dFracNoDomainDimChoice);
	fprintf(fp," bFastGas: %d",msr->param.bFastGas);
	fprintf(fp," dFracFastGas: %g",msr->param.dFracFastGas);
	fprintf(fp," dhMinOverSoft: %g",msr->param.dhMinOverSoft);
	fprintf(fp," bRungDD: %d",msr->param.bRungDD);
	fprintf(fp," dRungDDWeight: %g ",msr->param.dRungDDWeight);
	fprintf(fp,"\n# nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp," bLowerSoundSpeed: %d",msr->param.bLowerSoundSpeed);
	fprintf(fp," bShockTracker: %d",msr->param.bShockTracker);
	fprintf(fp," dShockTrackerA: %f",msr->param.dShockTrackerA);
	fprintf(fp," dShockTrackerB: %f",msr->param.dShockTrackerB);
	fprintf(fp,"\n# GROWMASS: nGrowMass: %d",msr->param.nGrowMass);
	fprintf(fp," dGrowDeltaM: %g",msr->param.dGrowDeltaM);
	fprintf(fp," dGrowStartT: %g",msr->param.dGrowStartT);
	fprintf(fp," dGrowEndT: %g",msr->param.dGrowEndT);
#ifdef GASOLINE
	fprintf(fp,"\n# SPH: bDoGas: %d",msr->param.bDoGas);	
	fprintf(fp," bGeometric: %d",msr->param.bGeometric);
	/* fprintf(fp," iGasModel: %d",msr->param.iGasModel); // Deprecated usage */
	fprintf(fp," bGasAdiabatic: %d",msr->param.bGasAdiabatic);	
	fprintf(fp," bGasIsothermal: %d",msr->param.bGasIsothermal);	
	fprintf(fp," bGasCooling: %d",msr->param.bGasCooling);	
	fprintf(fp," dConstAlpha: %g",msr->param.dConstAlpha);
	fprintf(fp," dConstBeta: %g",msr->param.dConstBeta);
	fprintf(fp,"\n# dConstGamma: %g",msr->param.dConstGamma);
	fprintf(fp," dMeanMolWeight: %g",msr->param.dMeanMolWeight);
	fprintf(fp," dGasConst: %g",msr->param.dGasConst);
	fprintf(fp," dKBoltzUnit: %g",msr->param.dKBoltzUnit);
	fprintf(fp," dMsolUnit: %g",msr->param.dMsolUnit);
	fprintf(fp," dKpcUnit: %g",msr->param.dKpcUnit);
	fprintf(fp," ddHonHLimit: %g",msr->param.ddHonHLimit);
	fprintf(fp,"\n# bViscosityLimiter: %d",msr->param.bViscosityLimiter);
	fprintf(fp," bBulkViscosity: %d",msr->param.bBulkViscosity);
	fprintf(fp," bGasDomainDecomp: %d",msr->param.bGasDomainDecomp);
	fprintf(fp," bSphStep: %d",msr->param.bSphStep);
	fprintf(fp,"\n#bSN: %d",msr->param.bSN);
	fprintf(fp," dSNRhoCut: %g",msr->param.dSNRhoCut);
 	fprintf(fp," dSNTMin: %g",msr->param.dSNTMin);
        fprintf(fp," dSNTMax: %g",msr->param.dSNTMax);
	fprintf(fp," dSNMetalCut: %g",msr->param.dSNMetalCut);
	fprintf(fp," dSNHeatFraction: %g",msr->param.dSNHeatFraction);
#ifndef NOCOOLING
	CoolLogParams( &msr->param.CoolParam, fp );
#endif
#endif
#ifdef STARFORM
	fprintf(fp,"\n# Star Formation: bStarForm: %d",msr->param.bStarForm);
	fprintf(fp," bFormOutputs: %d",msr->param.bFormOutputs);
	fprintf(fp," bFeedBack: %d",msr->param.bFeedBack);	
	fprintf(fp," dOverDenMin: %g",msr->param.stfm->dOverDenMin);
	fprintf(fp," dPhysDenMin: %g",msr->param.stfm->dPhysDenMin);
	fprintf(fp," dStarEff: %g",msr->param.stfm->dStarEff);
	fprintf(fp," dCStar: %g",msr->param.stfm->dCStar);
	fprintf(fp," dTempMax: %g",msr->param.stfm->dTempMax);
	fprintf(fp," dSoftMin: %g",msr->param.stfm->dSoftMin);
	fprintf(fp," dMinMassFrac: %g",msr->param.stfm->dMinMassFrac);
	fprintf(fp," dMinGasMass: %g",msr->param.stfm->dMinGasMass);
	fprintf(fp," dMaxStarMass: %g",msr->param.stfm->dMaxStarMass);
	fprintf(fp," dESN: %g",msr->param.sn->dESN);
	fprintf(fp," bSNTurnOffCooling: %i",msr->param.bSNTurnOffCooling);
	fprintf(fp," bShortCoolShutoff: %i",msr->param.bShortCoolShutoff);
	fprintf(fp," bSmallSNSmooth: %i",msr->param.bSmallSNSmooth);

        for ( testDelta = msr->param.dDelta; 
            testDelta >= msr->param.dDeltaStarForm && 
            msr->param.dDeltaStarForm > 0.0; testDelta *= 0.5 ){
                    if ( !(prmSpecified(msr->prm,"iStarFormRung")) )
                        msr->param.iStarFormRung++;
                    }
        if ( testDelta <= msr->param.dDelta ){ 
            fprintf(fp," dDeltaStarForm (set): %g, effectively: %g = %g yrs, iStarFormRung: %i",
                    msr->param.dDeltaStarForm, testDelta,
                    testDelta*msr->param.dSecUnit/SECONDSPERYEAR,
                    msr->param.iStarFormRung );
            msr->param.stfm->dDeltaT = msr->param.dDeltaStarForm = testDelta;
            }
        else if ( msr->param.dDeltaStarForm == 0.0 ) {
            fprintf(fp," dDeltaStarForm (set): %g, effectively: 0.0 = 0.0 yrs, iStarFormRung: maxRung",
                    msr->param.dDeltaStarForm );
            msr->param.iStarFormRung = msr->param.iMaxRung;
            }
        else {
            fprintf(fp," dDeltaStarForm (set): %g, effectively:  NO STARS WILL FORM", msr->param.dDeltaStarForm);
            }

#endif
  	  if (msr->param.bDoSinks) {
	      if (prmSpecified(msr->prm,"iSinkRung")) {
		  /* Find associated timestep for iSinkRung */
		  int iRung;
		  if (msr->param.iSinkRung > msr->param.iMaxRung) 
		      msr->param.iSinkRung = msr->param.iMaxRung;

		  testDelta = msr->param.dDelta;
		  for ( iRung = 0; iRung < msr->param.iSinkRung ; iRung++ ) 
		      testDelta *= 0.5;
		  }
	      else {
		  /* Find associate Rung for dDeltaSink */
		  /* NB: dDeltaSink is always in code units */
		  msr->param.iSinkRung = 0;
		  for ( testDelta = msr->param.dDelta; 
			testDelta > msr->param.dDeltaSink ; testDelta *= 0.5 ) {
                        msr->param.iSinkRung++;
			if (msr->param.iSinkRung >= msr->param.iMaxRung)
				exit(-1);
		      }
		  }		  
	      fprintf(fp," dDeltaSink (set): %g, effectively: %g = %g yrs, iSinkRung: %i",
		      msr->param.dDeltaSink, testDelta,
		      testDelta*msr->param.dSecUnit/SECONDSPERYEAR,
		      msr->param.iSinkRung );
	      msr->param.dDeltaSink = testDelta;
	      }


#ifdef SIMPLESF
	fprintf(fp,"\n# SSF: bStarForm: %d",msr->param.bStarForm);
	fprintf(fp," bFeedBack: %d",msr->param.bFeedBack);	
	fprintf(fp," SSF_dEfficiency: %g",msr->param.SSF_dEfficiency);
	fprintf(fp," SSF_dTMax: %g",msr->param.SSF_dTMax);
	fprintf(fp," SSF_dPhysDenMin: %g",msr->param.SSF_dPhysDenMin);
	fprintf(fp," SSF_dComovingDenMin: %g",msr->param.SSF_dComovingDenMin);
	fprintf(fp," SSF_dESNPerStarMass: %g",msr->param.SSF_dESNPerStarMass);
	fprintf(fp," SSF_dInitStarMass: %g",msr->param.SSF_dInitStarMass);
	fprintf(fp," SSF_dtCoolingShutoff: %g",msr->param.SSF_dtCoolingShutoff);
	fprintf(fp," SSF_bdivv: %d",msr->param.SSF_bdivv);
#endif
#ifdef SLIDING_PATCH
	fprintf(fp,"\n# bPatch: %d",msr->param.PP.bPatch);
	if (msr->param.PP.bPatch) {
		PATCH_PARAMS *PP = &msr->param.PP;
		fprintf(fp," dOrbDist: %g",PP->dOrbDist);
		fprintf(fp," (dOrbFreq: %g)",PP->dOrbFreq);
		fprintf(fp," bExtPert: %d",PP->bExtPert);
		if (PP->bExtPert) {
			fprintf(fp,"\n# dPertOrbDist: %g",PP->dPertOrbDist);
			fprintf(fp," dPertMass: %g",PP->dPertMass);
			fprintf(fp," dPertMaxZ: %g",PP->dPertMaxZ);
			fprintf(fp," dPertOrbFreqZ: %g",PP->dPertOrbFreqZ);
			fprintf(fp,"\n# dPertPhase: %g",PP->dPertPhase);
			fprintf(fp," dPertPhaseZ: %g",PP->dPertPhaseZ);
			fprintf(fp," (dPertOrbFreq: %g)",PP->dPertOrbFreq);
			}
		fprintf(fp,"\n# bRandAzWrap: %d",PP->bRandAzWrap);
		if (PP->bRandAzWrap) {
			fprintf(fp," bNoRandomX: %d",PP->bNoRandomX);
			fprintf(fp," nWrapAttempts: %d",PP->nWrapAttempts);
			fprintf(fp," dStripInner: %g",PP->dStripInner/PP->dWidth);
			fprintf(fp," dStripOuter: %g",PP->dStripOuter/PP->dWidth);
			fprintf(fp," dVelDispX: %g",PP->dVelDispX);
			fprintf(fp," dVelDispY: %g",PP->dVelDispY);
			fprintf(fp," dAvgVertAmp: %g",PP->dAvgVertAmp);
			fprintf(fp," dAvgMass: %g",PP->dAvgMass);
			}
		}
#endif /* SLIDING_PATCH */
#ifdef SIMPLE_GAS_DRAG
	fprintf(fp,"\n# bSimpleGasDrag: %d",msr->param.bSimpleGasDrag);
	fprintf(fp," bEpstein: %d",msr->param.bEpstein);
	fprintf(fp," dGamma: %g",msr->param.dGamma);
#endif /* SIMPLE_GAS_DRAG */
#ifdef SPINUP
	fprintf(fp,"\n# bUseAngMomSpinup: %d",msr->param.bUseAngMomSpinup);
	fprintf(fp," bSpinupAlignAxes: %d",msr->param.bSpinupAlignAxes);
	fprintf(fp," dSpinupPeriodStart: %g",msr->param.dSpinupPeriodStart);
	fprintf(fp," dSpinupPeriodInter: %g",msr->param.dSpinupPeriodInter);
	fprintf(fp," dSpinupPeriodFinal: %g",msr->param.dSpinupPeriodFinal);
	fprintf(fp," nSpinupStepsStart: %i",msr->param.nSpinupStepsStart);
	fprintf(fp," nSpinupStepsInter: %i",msr->param.nSpinupStepsInter);
	fprintf(fp," nSpinupStepsFinal: %i",msr->param.nSpinupStepsFinal);
	fprintf(fp," nSpinupStepsHold: %i",msr->param.nSpinupStepsHold);
	fprintf(fp," nSpinupStepsRestNum: %i",msr->param.nSpinupStepsRestNum);
	fprintf(fp," nSpinupStepsRestTime: %i",msr->param.nSpinupStepsRestTime);
	fprintf(fp," dSpinupRadiusInner: %g",msr->param.dSpinupRadiusInner);
	fprintf(fp," dSpinupRadiusOuter: %g",msr->param.dSpinupRadiusOuter);	
	fprintf(fp," dSpinupCoheCoeff: %g",msr->param.dSpinupCoheCoeff);
	fprintf(fp," nSpinupCoheTime: %i",msr->param.nSpinupCoheTime);
#endif /* SPINUP */
#ifdef COLLISIONS
	fprintf(fp,"\n# Collisions...");
	fprintf(fp," bAllowSimulColl: %i",msr->param.bAllowSimulColl);
	fprintf(fp," bFindRejects: %i",msr->param.bFindRejects);
	fprintf(fp," iCollLogOption: %i",msr->param.iCollLogOption);
	fprintf(fp," dSmallStep: %g",msr->param.dSmallStep);
	fprintf(fp," iMinCollRung: %i",msr->param.iMinCollRung);
	fprintf(fp,"\n# dxUnifGrav: %g",msr->param.dxUnifGrav);
	fprintf(fp," dyUnifGrav: %g",msr->param.dyUnifGrav);
	fprintf(fp," dzUnifGrav: %g",msr->param.dzUnifGrav);
    fprintf(fp,"\n# iOutcomes: %i",msr->param.CP.iOutcomes);
	fprintf(fp," dMergeLimit: %g",msr->param.CP.dMergeLimit);
	fprintf(fp," dDensity: %g",msr->param.CP.dDensity/DEN_CGS_SYS);
	fprintf(fp," iDensityAltCol: %i",msr->param.CP.iDensityAltCol);
	fprintf(fp," dDensityAltVal: %g",msr->param.CP.dDensityAltVal/DEN_CGS_SYS);
	fprintf(fp,"\n# iEpsNOption: %i",msr->param.CP.iEpsNOption);
    fprintf(fp," dEpsN: %g",msr->param.CP.dEpsN);
    fprintf(fp," dEpsT: %g",msr->param.CP.dEpsT);
	fprintf(fp," dEpsNCoef: %g",msr->param.CP.dEpsNCoef);
	fprintf(fp," dEpsNExp: %g",msr->param.CP.dEpsNExp);
	fprintf(fp," dEpsNVStar: %g",msr->param.CP.dEpsNVStar);
	fprintf(fp," dEpsNMin: %g",msr->param.CP.dEpsNMin);
    fprintf(fp,"\n# iMinBinaryRung: %i",msr->param.iMinBinaryRung);
    fprintf(fp," dBallVelFact: %g",msr->param.dBallVelFact);
    fprintf(fp," dMaxBinaryEcc: %g",msr->param.dMaxBinaryEcc);
    fprintf(fp,"\n# iSlideOption: %i",msr->param.CP.iSlideOption);
	fprintf(fp," dSlideLimit: %g",msr->param.CP.dSlideLimit);
    fprintf(fp," dSlideEpsN: %g",msr->param.CP.dSlideEpsN);
    fprintf(fp," dSlideEpsT: %g",msr->param.CP.dSlideEpsT);
    fprintf(fp,"\n# dCollapseLimit: %g",msr->param.CP.dCollapseLimit);
    fprintf(fp," dCollapseEpsN: %g",msr->param.CP.dCollapseEpsN);
    fprintf(fp," dCollapseEpsT: %g",msr->param.CP.dCollapseEpsT);
    fprintf(fp,"\n# dCrushLimit: %g",msr->param.CP.dCrushLimit);
    fprintf(fp," dCrushEpsN: %g",msr->param.CP.dCrushEpsN);
    fprintf(fp," dCrushEpsT: %g",msr->param.CP.dCrushEpsT);
	fprintf(fp,"\n# iOverlapOption: %i",msr->param.CP.iOverlapOption);
	fprintf(fp," bStrictOverlap: %i",msr->param.CP.bStrictOverlap);
	fprintf(fp," dBackstepLimit: %g",msr->param.CP.dBackstepLimit);
	fprintf(fp," dAdjPosLimit: %g",msr->param.CP.dAdjPosLimit);
	fprintf(fp," dRepelFac: %g",msr->param.CP.dRepelFac);
	fprintf(fp,"\n# dFragLimit: %g",msr->param.CP.dFragLimit);
	fprintf(fp,"\n# nUnifGrav: %i",msr->param.sTimeVarParams.nUnifGrav);
#endif
#ifdef RUBBLE_ZML
	fprintf(fp,"\n# Rubble...");
	fprintf(fp," dRubbleMinFracMass: %g",msr->param.CP.dRubbleMinFracMass);
	fprintf(fp," dRubMinMass: %g",msr->param.CP.DB.dRubMinMass);
	fprintf(fp," iRubNumDynToBounce: %i",msr->param.CP.DB.iRubNumDynToBounce);
	fprintf(fp," iRubNumDynToMerge: %i",msr->param.CP.DB.iRubNumDynToMerge);
	fprintf(fp,"\n# nDustBinsVelDispOpt: %i",msr->param.CP.DB.iDustBinsVelDispOpt);
	fprintf(fp," nDustBins: %i",msr->param.CP.DB.nDustBins);
	fprintf(fp," iDustBinsApplyInt: %i",msr->param.CP.DB.iDustBinsApplyInt);
	fprintf(fp," dDustBinsInner: %g",msr->param.CP.DB.dDustBinsInner);
	fprintf(fp," dDustBinsOuter: %g",msr->param.CP.DB.dDustBinsOuter);
	fprintf(fp,"\n# dDustBinsScaleHeight: %g",msr->param.CP.DB.dDustBinsScaleHeight);
	fprintf(fp," dDustBinsInitSigma: %g",msr->param.CP.DB.dDustBinsInitSigma);
	fprintf(fp," dDustBinsInitAlpha: %g",msr->param.CP.DB.dDustBinsInitAlpha);
	fprintf(fp," (dDustBinsWidth: %g)",msr->param.CP.DB.dDustBinsWidth);
#elif defined(COLLMOD_ZML)
	fprintf(fp,"\n# Collision model...");
	fprintf(fp," dColModMinMass: %g",msr->param.CP.DB.dCollMinMass);
	fprintf(fp,"\n# nDustBinsVelDispOpt: %i",msr->param.CP.DB.iDustBinsVelDispOpt);
	fprintf(fp," nDustBins: %i",msr->param.CP.DB.nDustBins);
	fprintf(fp," iDustBinsApplyInt: %i",msr->param.CP.DB.iDustBinsApplyInt);
	fprintf(fp," dDustBinsInner: %g",msr->param.CP.DB.dDustBinsInner);
	fprintf(fp," dDustBinsOuter: %g",msr->param.CP.DB.dDustBinsOuter);
	fprintf(fp,"\n# dDustBinsScaleHeight: %g",msr->param.CP.DB.dDustBinsScaleHeight);
	fprintf(fp," dDustBinsInitSigma: %g",msr->param.CP.DB.dDustBinsInitSigma);
	fprintf(fp," dDustBinsInitAlpha: %g",msr->param.CP.DB.dDustBinsInitAlpha);
	fprintf(fp," (dDustBinsWidth: %g)",msr->param.CP.DB.dDustBinsWidth);
#endif /* RUBBLE_ZML, COLLMOD_ZML */
#ifdef ORIGIN_HISTOGRAM
	fprintf(fp," NUM_ORIGIN_BINS: %i", NUM_ORIGIN_BINS);
#endif /* ORIGIN_HISTOGRAM */ 
#ifdef AGGS
	fprintf(fp,"\n# Aggs...");
	fprintf(fp," bAggsSolveQuartic: %i",msr->param.CP.bAggsSolveQuartic);
	fprintf(fp," dTensileCoef: %g",msr->param.CP.SP.dTensileCoef);
	fprintf(fp," dTensileExp: %g",msr->param.CP.SP.dTensileExp);
	fprintf(fp," dShearCoef: %g",msr->param.CP.SP.dShearCoef);
	fprintf(fp," dShearExp: %g",msr->param.CP.SP.dShearExp);
#endif
#ifdef WALLS
	{
		const WALL_PARAMS *WP = &msr->param.CP.WP;
		const WALL_DATA *w;
		int j;
		fprintf(fp,"\n# nWalls: %i",WP->nWalls);
		for (j=0;j<WP->nWalls;j++) {
			assert(j == WP->pWalls[j].iWallID);
			w = &WP->pWalls[j].wd;
			fprintf(fp,"\n# Wall %i -- type: %i origin: %g,%g,%g orient: %g,%g,%g vertex1: %g,%g,%g vertex2: %g,%g,%g velocity: %g,%g,%g step: %g osc-ampl: %g osc-freq: %g osc-phase: %gosc-vec: %g %g %g radius: %g hole-radius: %g length: %g taper: %g open-angle: %g ang-speed: %g epsn: %g epst: %g color: %i transparency: %g",j,w->iType,w->vOrigin[0],w->vOrigin[1],w->vOrigin[2],w->vOrient[0],w->vOrient[1],w->vOrient[2],w->vVertex1[0],w->vVertex1[1],w->vVertex1[2],w->vVertex2[0],w->vVertex2[1],w->vVertex2[2],w->vVel[0],w->vVel[1],w->vVel[2],w->dStep,w->dOscAmp,w->dOscFreq,w->dOscPhase,w->vOscVec[0],w->vOscVec[1],w->vOscVec[2],w->dRadius,w->dHoleRadius,w->dLength,w->dTaper,w->dOpenAngle,w->dAngSpeed,w->dEpsN,w->dEpsT,w->iColor,w->dTrans);
#ifdef DEM
			fprintf(fp," k_n: %g k_t: %g mu_s: %g mu_r: %g mu_t: %g",w->dKn,w->dKt,w->dMuS,w->dMuR,w->dMuT);
#ifdef DEM_ROTATION_DASHPOT
			fprintf(fp," beta_s: %g",w->dBetaS);
#endif /* DEM_ROTATION_DASHPOT */
#ifdef DEM_COHESION
			fprintf(fp," cohe_model: %i cohe_coeff: %g",w->iCohesionModel,w->dCohesiveCoeff);
#endif /* DEM_COHESION */
#ifndef DEM_TWOLAYERS
			fprintf(fp," inner-overlap-boundary: %g k_n_outer: %g k_t_outer: %g k_n_inner: %g k_t_inner: %g",
					w->dInnerOverlapBndry,w->dKnOuter,w->dKtOuter,w->dKnInner,w->dKtInner);
#endif /* DEM_TWOLAYERS */
#ifdef DEM_WALLS_REACT
		fprintf(fp,"mass: %g cog: %g,%g,%g rot-vec: %g,%g,%g Ixj: %g,%g,%g Iyj: %g,%g,%g Izj: %g,%g,%g Pxj: %g,%g,%g Pyj: %g,%g,%g Pzj: %g,%g,%g assembly: %i",w->dMass,w->vCoG[0],w->vCoG[1],w->vCoG[2],w->vRot[0],w->vRot[1],w->vRot[2],w->vIxj[0],w->vIxj[1],w->vIxj[2],w->vIyj[0],w->vIyj[1],w->vIyj[2],w->vIzj[0],w->vIzj[1],w->vIzj[2],w->vPxj[0],w->vPxj[1],w->vPxj[2],w->vPyj[0],w->vPyj[1],w->vPyj[2],w->vPzj[0],w->vPzj[1],w->vPzj[2],w->iAssembly);
		if (w->dMass > 0.0 && w->dAngSpeed != 0.0) {
		    puts("ERROR: angular speed of reactive walls must be zero.");
		    _msrExit(msr,1);
		    }
#endif /* DEM_WALLS_REACT */
#endif /* DEM */
			}
		fprintf(fp,"\n# bWallsEdgeDetect: %i bWallsSolveQuartic: %i",WP->bWallsEdgeDetect,WP->bWallsSolveQuartic);
		}
#endif /* WALLS */
#ifdef SPECIAL_PARTICLES
	{
	SPECIAL_PARTICLE_DATA *s;
	int j;
	fprintf(fp,"\n# nSpecial: %i",msr->param.nSpecial);
	for (j=0;j<msr->param.nSpecial;j++) {
		s = &msr->param.sSpecialData[j];
		fprintf(fp,"\n# Special %i: org_idx: %i type: %i ",j,
				msr->param.iSpecialId[j],s->iType);
		if (s->iType & SPECIAL_OBLATE)
			fprintf(fp,"dRadEq: %g J2: %g J4: %g p[0]: %g p[1]: %g p[2]: %g",
					s->oblate.dRadEq,s->oblate.J2,s->oblate.J4,
					s->oblate.p[0],s->oblate.p[1],s->oblate.p[2]);
		if (s->iType & SPECIAL_GR)
			fprintf(fp,"UNSUPPORTED");
		if (s->iType & SPECIAL_FORCE)
			fprintf(fp,"dMag: %g",s->force.dMag);
	        if (s->iType & SPECIAL_NOGHOSTPERT) 
		        fprintf(fp,"Excluded from ghost cell gravity calculations");
		}
	}
#endif
	fprintf(fp,"\n# iForceOverrideOption: %i",msr->param.FO.iForceOverrideOption);
	switch (msr->param.FO.iForceOverrideOption) {
	case FO_NONE:
		break;
	case FO_VANDERWAALS:
		assert(0);
	case FO_STRENGTH:
#ifdef SPRINGS
		fprintf(fp,"\n# dMeanYoungsModulus: %g",msr->param.FO.SP.dMeanYoungsModulus);
		fprintf(fp," dMeanStrainLimit: %g",msr->param.FO.SP.dMeanStrainLimit);
		fprintf(fp," dYoungsStdDev: %g",msr->param.FO.SP.dYoungsStdDev);
		fprintf(fp," dStrainStdDev: %g",msr->param.FO.SP.dStrainStdDev);
		fprintf(fp," dMaxStrainLimit: %g",msr->param.FO.SP.dMaxStrainLimit);
		fprintf(fp," dLinkageLength: %g",msr->param.FO.SP.dLinkageLength);
		fprintf(fp," dZeroStrainLength: %g",msr->param.FO.SP.dZeroStrainLength);
		fprintf(fp," dDamp: %g",msr->param.FO.SP.dDamp);
		fprintf(fp," dPorosityFactor: %g",msr->param.FO.SP.dPorosityFactor);
		fprintf(fp," bReadSpringsData: %i",msr->param.FO.SP.bReadSpringsData);
		fprintf(fp," bNoSpringsOutput: %i",msr->param.FO.SP.bNoSpringsOutput);
#endif /* SPRINGS */
#ifdef DEM
		fprintf(fp,"\n# dKn: %g",msr->param.FO.DP.dKn);
		fprintf(fp," dKt: %g",msr->param.FO.DP.dKt);
		fprintf(fp," dMuS: %g",msr->param.FO.DP.dMuS);
		fprintf(fp," dMuR: %g",msr->param.FO.DP.dMuR);
		fprintf(fp," dMuT: %g",msr->param.FO.DP.dMuT);
		fprintf(fp," dAccCrit: %g",msr->param.FO.DP.dAccCrit);
		fprintf(fp,"\n# dMinorFrac: %g",msr->param.FO.DP.dMinorFrac);
		fprintf(fp," dMajorFrac: %g",msr->param.FO.DP.dMajorFrac);
		fprintf(fp," dErrorFrac: %g",msr->param.FO.DP.dErrorFrac);
		fprintf(fp," bReadDEMData: %i",msr->param.FO.DP.bReadDEMData);
#ifdef DEM_TRACK_ORIENT
		fprintf(fp," bReadOrient: %i",msr->param.FO.DP.bReadOrient);
#endif /* DEM_TRACK_ORIENT */
		fprintf(fp," nDEMOutputs: %i",msr->param.FO.DP.nDEMOutputs);
		fprintf(fp," iDEMStatsInterval: %i",msr->param.FO.DP.iDEMStatsInterval);
#ifdef DEM_FIXED_BALL
		fprintf(fp," iFixedBallOption: %i",msr->param.FO.DP.iFixedBallOption);
#endif /* DEM_FIXED_BALL */
#ifdef DEM_ROTATION_DASHPOT
		fprintf(fp," dBetaS: %g",msr->param.FO.DP.dBetaS);
		fprintf(fp," bUseLegacyRoll: %i",msr->param.FO.DP.bUseLegacyRoll);
		fprintf(fp," bUseContactRadius: %i",msr->param.FO.DP.bUseContactRadius);
#endif /* DEM_ROTATION_DASHPOT */
#ifdef DEM_COHESION
		fprintf(fp," iCohesionModel: %i",msr->param.FO.DP.iCohesionModel);
		fprintf(fp," dCohesiveCoeff: %g",msr->param.FO.DP.dCohesiveCoeff);
#endif /* DEM_COHESION */
#ifdef DEM_TWOLAYERS
		fprintf(fp," dInnerOverlapBndry: %g",msr->param.FO.DP.dInnerOverlapBndry);
		fprintf(fp," dKnOuter: %g",msr->param.FO.DP.dKnOuter);
		fprintf(fp," dKtOuter: %g",msr->param.FO.DP.dKtOuter);
		fprintf(fp," dKnInner: %g",msr->param.FO.DP.dKnInner);
		fprintf(fp," dKtInner: %g",msr->param.FO.DP.dKtInner);
#endif /* DEM_TWOLAYERS */
#endif /* DEM */
#if !defined(SPRINGS) && !defined(DEM)
		assert(0);
#endif /* !SPRINGS && !DEM */
		break;
	default:
		assert(0);
		}
	switch (msr->iOpenType) {
	case OPEN_JOSH:
		fprintf(fp,"\n# iOpenType: JOSH");
		break;
	case OPEN_ABSPAR:
		fprintf(fp,"\n# iOpenType: ABSPAR");
		break;
	case OPEN_RELPAR:
		fprintf(fp,"\n# iOpenType: RELPAR");
		break;
	case OPEN_ABSTOT:
		fprintf(fp,"\n# iOpenType: ABSTOT");
		break;
	case OPEN_RELTOT:
		fprintf(fp,"\n# iOpenType: RELTOT");
		break;
	default:
		fprintf(fp,"\n# iOpenType: NONE?");
		}
	fprintf(fp," dTheta: %f",msr->param.dTheta);
	fprintf(fp,"\n# dAbsPartial: %g",msr->param.dAbsPartial);
	fprintf(fp," dRealPartial: %g",msr->param.dRelPartial);
	fprintf(fp," dAbsTotal: %g",msr->param.dAbsTotal);
	fprintf(fp," dRelTotal: %g",msr->param.dRelTotal);
	fprintf(fp,"\n# dPeriod: %g",msr->param.dPeriod);
	fprintf(fp," dxPeriod: %g",
			msr->param.dxPeriod == FLOAT_MAXVAL ? 0.0 : msr->param.dxPeriod);
	fprintf(fp," dyPeriod: %g",
			msr->param.dyPeriod == FLOAT_MAXVAL ? 0.0 : msr->param.dyPeriod);
	fprintf(fp," dzPeriod: %g",
			msr->param.dzPeriod == FLOAT_MAXVAL ? 0.0 : msr->param.dzPeriod);
	fprintf(fp,"\n# dHubble0: %g",msr->param.csm->dHubble0);
	fprintf(fp," dOmega0: %g",msr->param.csm->dOmega0);
	fprintf(fp," dLambda: %g",msr->param.csm->dLambda);
	fprintf(fp," dOmegaRad: %g",msr->param.csm->dOmegaRad);
	fprintf(fp," dOmegab: %g",msr->param.csm->dOmegab);
	fprintf(fp," dQuintess: %g",msr->param.csm->dQuintess);
	fprintf(fp,"\n# achInFile: %s",msr->param.achInFile);
	fprintf(fp,"\n# achOutName: %s",msr->param.achOutName); 
	fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
	if (msrComove(msr)) {
		fprintf(fp,"\n# RedOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			z = 1.0/csmTime2Exp(msr->param.csm, msr->pdOutTime[i]) - 1.0;
			fprintf(fp," %f",z);
			}
		fprintf(fp,"\n");
		}
	else {
		fprintf(fp,"\n# TimeOut:");
		if (msr->nOuts == 0) fprintf(fp," none");
		for (i=0;i<msr->nOuts;i++) {
			if (i%5 == 0) fprintf(fp,"\n#   ");
			fprintf(fp," %f",msr->pdOutTime[i]);
			}
		fprintf(fp,"\n");
	    }
    }

int msrCheckForInterrupt(MSR msr)
{
	/*
	 ** Checks for existence of command in INTR_FILE in run directory.
	 ** If found, the file is removed and the return status is set to
	 ** the interrupt token, otherwise INTR_NONE.  For backwards
	 ** compatability, if the file is empty, INTR_STOP is assumed.
	 ** The tokens are defined in master.h (search for "interrupt").
	 */
	
	static char INTERRUPT[INTERRUPTS][256] = {
	  /* these must match tokens in master.h */
	  "NONE",
	  "STOP",
	  "STATUS",
	  "CHECK",
	  "OUTPUT"
	};

	FILE *fp = NULL;
	char achFile[PST_FILENAME_SIZE],achTmp[PST_FILENAME_SIZE],*rv;
	int i,iReturn=INTR_NONE;

	_msrMakePath(msr->param.achDataSubPath,INTR_FILE,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFile);
	
	if ((fp = fopen(achFile,"r")) != NULL) {
		rv = fgets(achTmp,256,fp);
		if (rv == NULL) { /* empty (or unreadable) file */
			printf("User stop request detected.\n");
			iReturn = INTR_STOP;
			}
		else {
			if (strlen(achTmp) > 0)
				achTmp[strlen(achTmp) - 1] = '\0'; /* chomp '\n' */
			for (i=1;i<INTERRUPTS;i++) /* 0 is INT_NONE */	
				if (strcasecmp(achTmp,INTERRUPT[i]) == 0) {
					printf("User interrupt detected: %s\n",INTERRUPT[i]);
					iReturn = i;
					break;
					}
			if (iReturn == INTR_NONE)
				fprintf(stderr,"User interrupt ERROR: unrecognized command \"%s\".\n",achTmp);
			}
		fclose(fp);
		unlink(achFile);
		}
	
	return iReturn;
	}

void msrFinish(MSR msr)
{
	int id;

	for (id=1;id<msr->mdl->nThreads;++id) {
		if (msr->param.bVDetails) printf("Stopping thread %d\n",id);		
		mdlReqService(msr->mdl,id,SRV_STOP,NULL,0);
		mdlGetReply(msr->mdl,id,NULL,NULL);
		}
	pstFinish(msr->pst);
#ifdef AGGS
	if (msr->pAggs != NULL)
		free((void *) msr->pAggs);
#endif
	/*
	 ** finish with parameter stuff, deallocate and exit.
	 */
	prmFinish(msr->prm);
	if (msr->pMap != NULL)
		free(msr->pMap);
	free(msr);
	}

void msrOneNodeReadTipsy(MSR msr, struct inReadTipsy *in)
{
    int i,id;
    int *nParts;				/* number of particles for each processor */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    /*struct inSetParticleTypes intype; -- not used: see JW's comment below -- DCR 12/19/02*/

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert(nid == (int) (msr->nThreads*sizeof(*nParts)));
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadTipsy(plcl->pkd,achInFile,nStart,nParts[id],
					 in->bStandard,in->iReadIOrder,in->dvFac,in->dTuFac);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadTipsy(plcl->pkd,achInFile,0,nParts[0],in->bStandard,in->iReadIOrder,in->dvFac,
				 in->dTuFac);

/* I think this code can be removed -- this is done later  JW Sept 2002 */
/*
	if (msr->param.iReadIOrder) {
		struct outGetNParts outget;
		struct inSetNParts inset;
		
		pstGetNParts(msr->pst,NULL,0,&outget,NULL);
		assert(outget.nGas == msr->nGas);
		assert(outget.nDark == msr->nDark);
		assert(outget.nStar == msr->nStar);
		inset.nGas = outget.nGas;
		inset.nDark = outget.nDark;
		inset.nStar = outget.nStar;
		msr->nMaxOrderGas = inset.nMaxOrderGas = outget.iMaxOrderGas;
		msr->nMaxOrderDark = inset.nMaxOrderDark = outget.iMaxOrderDark;
        msr->nMaxOrder = inset.nMaxOrder     = outget.iMaxOrderStar;
		pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
		}

    intype.nSuperCool = msr->param.nSuperCool;
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
*/
    }

int xdrHeader(XDR *pxdrs,struct dump *ph)
{
	int pad = 0;
	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
	return 1;
	}


double msrReadTipsy(MSR msr)
{
	FILE *fp;
	struct dump h;
	struct inReadTipsy in;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double dTime,aTo,tTo,z;
	struct inSetParticleTypes intype;
	double sec=0.0,dsec=0.0;
	
	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,in.achInFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achInFile,achInFile);

		fp = fopen(achInFile,"r");
		if (!fp) {
			printf("Could not open InFile:%s\n",achInFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No input file specified\n");
		_msrExit(msr,1);
		return -1.0;
		}
	/*
	 ** Assume tipsy format for now, and dark matter only.
	 */
	if (msr->param.bStandard) {
		XDR xdrs;

		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		xdrHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
		}
	else {
		fread(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);

	msr->N = h.nbodies;
	msr->nDark = h.ndark;
	msr->nGas = h.nsph;
	msr->nStar = h.nstar;
	msr->nMaxOrder = msr->N - 1;
	msr->nMaxOrderGas = msr->nGas - 1;
	msr->nMaxOrderDark = msr->nGas + msr->nDark - 1;

	assert(msr->N == msr->nDark+msr->nGas+msr->nStar);
#ifndef GASOLINE
	if (msr->nGas != 0) fprintf(stderr,"GASOLINE compile flag not set:  Treating %d Gas particles as Dark\n",msr->nGas);
#endif
	if (msrComove(msr)) {
		if(msr->param.csm->dHubble0 == 0.0) {
			printf("No hubble constant specified\n");
			_msrExit(msr,1);
			}
		dTime = csmExp2Time(msr->param.csm,h.time);
		z = 1.0/h.time - 1.0;
		if (msr->param.bVStart)
			printf("Input file, Time:%g Redshift:%g Expansion factor:%g\n",
				   dTime,z,h.time);
		if (prmSpecified(msr->prm,"dRedTo")) {
			if (!prmArgSpecified(msr->prm,"nSteps") &&
				prmArgSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr,1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmArgSpecified(msr->prm,"dDelta") &&
					 prmArgSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr,1);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta =
					(tTo-dTime)/(msr->param.nSteps -
						     msr->param.iStartStep);
				
				else
				    msr->param.dDelta = 0.0;
				}
			else if (!prmSpecified(msr->prm,"nSteps") &&
					 prmFileSpecified(msr->prm,"dDelta")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");
					_msrExit(msr,1);
					}
				msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
				}
			else if (!prmSpecified(msr->prm,"dDelta") &&
					 prmFileSpecified(msr->prm,"nSteps")) {
				aTo = 1.0/(msr->param.dRedTo + 1.0);
				tTo = csmExp2Time(msr->param.csm,aTo);
				if (msr->param.bVStart)
					printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
						   tTo,1.0/aTo-1.0,aTo);
				if (tTo < dTime) {
					printf("Badly specified final redshift, check -zto parameter.\n");	
					_msrExit(msr,1);
					}
				if(msr->param.nSteps != 0)
				    msr->param.dDelta =	(tTo-dTime)/(msr->param.nSteps
													 - msr->param.iStartStep);
				else
				    msr->param.dDelta = 0.0;
				}
			}
		else {
			tTo = dTime + msr->param.nSteps*msr->param.dDelta;
			aTo = csmTime2Exp(msr->param.csm,tTo);
			if (msr->param.bVStart)
				printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
					   tTo,1.0/aTo-1.0,aTo);
			}
		if (msr->param.bVStart)
			printf("Reading file...\nN:%d nDark:%d nGas:%d nStar:%d\n",msr->N,
				   msr->nDark,msr->nGas,msr->nStar);
		if (msr->param.bCannonical) {
			in.dvFac = h.time*h.time;
			}
		else {
			in.dvFac = 1.0;
			}
		}
	else {
		dTime = h.time;
		if (msr->param.bVStart) printf("Input file, Time:%g\n",dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		if (msr->param.bVStart) {
			printf("Simulation to Time:%g\n",tTo);
			printf("Reading file...\nN:%d nDark:%d nGas:%d nStar:%d Time:%g\n",
				   msr->N,msr->nDark,msr->nGas,msr->nStar,dTime);
			}
		in.dvFac = 1.0;
		}
	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;
	in.nStar = msr->nStar;
	in.iOrder = msr->param.iOrder;
	in.bStandard = msr->param.bStandard;
	in.iReadIOrder = msr->param.iReadIOrder;
#ifdef GASOLINE
	in.dTuFac = msr->param.dGasConst/(msr->param.dConstGamma - 1)/
		msr->param.dMeanMolWeight;
#else
	in.dTuFac = 1.0;
#endif
	/*
	 ** Since pstReadTipsy causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;
	/*
	 ** Provide the period.
	 */
	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;

	/* Read Timings --JPG */
	if (msr->param.bVDetails) {
	     printf("Reading input file data...\n");
	     sec = msrTime(msr);
	}

	if(msr->param.bParaRead)
	    pstReadTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadTipsy(msr, &in);

	if (msr->param.bVDetails) {
	     dsec = msrTime(msr) - sec;
		 printf("Data read complete, Wallclock: %f secs\n",dsec);
	}

	if (msr->param.iReadIOrder) {
		struct outGetNParts outget;
		struct inSetNParts inset;
		
		pstGetNParts(msr->pst,NULL,0,&outget,NULL);
		assert(outget.nGas == msr->nGas);
		assert(outget.nDark == msr->nDark);
		assert(outget.nStar == msr->nStar);
		inset.nGas = outget.nGas;
		inset.nDark = outget.nDark;
		inset.nStar = outget.nStar;
		msr->nMaxOrderGas = inset.nMaxOrderGas = outget.iMaxOrderGas;
		msr->nMaxOrderDark = inset.nMaxOrderDark = outget.iMaxOrderDark;
        msr->nMaxOrder = inset.nMaxOrder     = outget.iMaxOrderStar;
		pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
		if (msr->param.bVDetails) puts("IOrder file has been successfully read.");
		}

	intype.nSuperCool = msr->param.nSuperCool;
	pstSetParticleTypes(msr->pst, &intype, sizeof(intype), NULL, NULL);
	if (msr->param.bVDetails) puts("Input file has been successfully read.");
	/*
	 ** Now read in the output points, passing the initial time.
	 ** We do this only if nSteps is not equal to zero.
	 */
	if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}


/*
 ** This function makes some DANGEROUS assumptions!!!
 ** Main problem is that it calls pkd level routines, bypassing the
 ** pst level. It uses plcl pointer which is not desirable.
 */
void msrOneNodeWriteTipsy(MSR msr, struct inWriteTipsy *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,in->bStandard,
				  in->dvFac,in->duTFac,in->iGasModel); 
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteTipsy(plcl->pkd,achOutFile,nStart,
					  in->bStandard, in->dvFac, in->duTFac,in->iGasModel); 
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }


void msrCalcWriteStart(MSR msr) 
{
	struct outSetTotal out;
	struct inSetWriteStart in;

	pstSetTotal(msr->pst,NULL,0,&out,NULL);
	assert(out.nTotal == msr->N);
	in.nWriteStart = 0;
	pstSetWriteStart(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrCalcNCWriteStart(MSR msr) 
{
	struct outSetTotals out;
	struct inSetNCWriteStart in;

	pstSetTotals(msr->pst,NULL,0,&out,NULL);
	assert(out.nGas == msr->nGas);
	assert(out.nDark == msr->nDark);
	assert(out.nStar == msr->nStar);
	in.nGasWriteStart = 0;
	in.nDarkWriteStart = 0;
	in.nStarWriteStart = 0;
	pstSetNCWriteStart(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrWriteTipsy(MSR msr,char *pszFileName,double dTime)
{
	FILE *fp;
	struct dump h;
	struct inWriteTipsy in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double sec=0.0,dsec=0.0;

	/*
	 ** Calculate where to start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	_msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	
	fp = fopen(achOutFile,"w");
	if (!fp) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}
	in.bStandard = msr->param.bStandard;
#ifdef GASOLINE
	in.duTFac = (msr->param.dConstGamma - 1)*msr->param.dMeanMolWeight/
		msr->param.dGasConst;
#else
	in.duTFac = 1.0;
#endif
	in.iGasModel = msr->param.iGasModel;
	/*
	 ** Assume tipsy format for now.
	 */
	h.nbodies = msr->N;
	h.ndark = msr->nDark;
	h.nsph = msr->nGas;
	h.nstar = msr->nStar;
	if (msrComove(msr)) {
		h.time = csmTime2Exp(msr->param.csm,dTime);
		if (msr->param.bCannonical) {
			in.dvFac = 1.0/(h.time*h.time);
			}
		else {
			in.dvFac = 1.0;
			}
		}
	else {
		h.time = dTime;
		in.dvFac = 1.0;
		}
	h.ndim = 3;
	if (msr->param.bVDetails) {
		if (msrComove(msr)) {
			printf("Writing file...\nTime:%g Redshift:%g\n",
				   dTime,(1.0/h.time - 1.0));
			}
		else {
			printf("Writing file...\nTime:%g\n",dTime);
			}
		}
	if (in.bStandard) {
		XDR xdrs;

		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdrHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
		}
	else {
		fwrite(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);

	/* Write Timings --JPG */
	if (msr->param.bVDetails) {
	     puts("Writing output file data...");
	     sec = msrTime(msr);
	     }
	if(msr->param.bParaWrite)
	    pstWriteTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteTipsy(msr, &in);
	if (msr->param.bVDetails) {
	     dsec = msrTime(msr) - sec;
		 printf("Data write complete, Wallclock: %f secs\n",dsec);
		 puts("Output file has been successfully written, Wallclock");
	     }
	}

void msrWriteTipsyHead(MSR msr,char *achOutFile,double dTime, struct inWriteTipsy *in)
{
	FILE *fp;
	struct dump h;
	
	fp = fopen(achOutFile,"w");
        assert(fp);
/*	if (!fp) {
            printf("Could not open OutFile:%s\n",achOutFile);
            _msrExit(msr,1);
            }*/
	in->bStandard = msr->param.bStandard;
#ifdef GASOLINE
	in->duTFac = (msr->param.dConstGamma - 1)*msr->param.dMeanMolWeight/
		msr->param.dGasConst;
#else
	in->duTFac = 1.0;
#endif
	in->iGasModel = msr->param.iGasModel;
	/*
	 ** Assume tipsy format for now.
	 */
	h.nbodies = msr->N;
	h.ndark = msr->nDark;
	h.nsph = msr->nGas;
	h.nstar = msr->nStar;
	if (msrComove(msr)) {
		h.time = csmTime2Exp(msr->param.csm,dTime);
		if (msr->param.bCannonical) {
			in->dvFac = 1.0/(h.time*h.time);
			}
		else {
			in->dvFac = 1.0;
			}
		}
	else {
		h.time = dTime;
		in->dvFac = 1.0;
		}
	h.ndim = 3;
	if (msr->param.bVDetails) {
		if (msrComove(msr)) {
			printf("Writing file...\nTime:%g Redshift:%g\n",
				   dTime,(1.0/h.time - 1.0));
			}
		else {
			printf("Writing file...\nTime:%g\n",dTime);
			}
		}
	if (in->bStandard) {
		XDR xdrs;

		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdrHeader(&xdrs,&h);
		xdr_destroy(&xdrs);
		}
	else {
		fwrite(&h,sizeof(struct dump),1,fp);
		}
	fclose(fp);
    }
    
void msrWriteTipsyBody(MSR msr,char *pszFileName,double dTime, struct inWriteTipsy *in)
{
	FILE *fp;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double sec=0.0,dsec=0.0;

	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);
	
	fp = fopen(achOutFile,"w");
        assert(fp);
/*	if (!fp) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}*/
	/* Write Timings --JPG */
	if (msr->param.bVDetails) {
	     puts("Writing output file data...");
	     sec = msrTime(msr);
	     }
	if(msr->param.bParaWrite)
	    pstWriteTipsy(msr->pst,in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteTipsy(msr, in);
	if (msr->param.bVDetails) {
	     dsec = msrTime(msr) - sec;
	     printf("Data write complete, Wallclock: %f secs\n",dsec);
		 puts("Output file has been successfully written, Wallclock");
	     }
	}

void msrSetSoft(MSR msr,double dSoft)
{
	struct inSetSoft in;
  
	if (msr->param.bVDetails) printf("Set Softening...\n");
	in.dSoft = dSoft;
	pstSetSoft(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrSetSink(MSR msr) 
{
    struct inSetSink in;
    struct outSetSink out;
	
	if (msr->param.bDoSinks) {
	  in.dSinkMassMin = msr->param.dSinkMassMin;
	  pstSetSink(msr->pst,&in,sizeof(in),&out,NULL);
	  if (msr->param.bVDetails) printf("Identified %d sink particles\n",out.nSink);
	  msr->nSink = out.nSink;
	  }
    }


void msrDomainDecomp(MSR msr, int iRung, int bGreater)
{
	struct inDomainDecomp in;
	int iRungDD,iRungSD,nActive;

#ifdef GASOLINE
	/* Sanity check on Gas particles being present */
	if (msr->nGas==0 && (msr->param.bDoGas==1 || msr->param.bGasDomainDecomp)) {
		if (msr->param.bGasDomainDecomp) {
			printf("MDD: switching bGasDomainDecomp off.\n");
			msr->param.bGasDomainDecomp=0;
			msr->bDoneDomainDecomp=0;
			}
		if (msr->param.bDoGas==1) {
			printf("MDD: switching bDoGas off.\n");
			msr->param.bDoGas=0;
			}
		}
#endif

	in.bDoRootFind = 1;
	in.bDoSplitDimFind = 1;
	
	nActive=0;
	if (bGreater) {
		iRungDD=msr->iCurrMaxRung+1; 
		while (iRungDD > iRung) {
			iRungDD--;
			nActive+=msr->nRung[iRungDD];
			}
		while(iRungDD > 0 && nActive < msr->N*msr->param.dFracNoDomainDecomp) {
			iRungDD--;
			nActive+=msr->nRung[iRungDD];
			}
		iRungSD = iRungDD;
		while(iRungSD > 0 && nActive < msr->N*msr->param.dFracNoDomainDimChoice) {
			iRungSD--;
			nActive+=msr->nRung[iRungSD];
			}
		}
	else {
		iRungDD = iRung;
		while(iRungDD > 0 && msr->nRung[iRungDD] < msr->N*msr->param.dFracNoDomainDecomp) {
			iRungDD--;
			}
		iRungSD = iRungDD;
		while(iRungSD > 0 && msr->nRung[iRungSD] < msr->N*msr->param.dFracNoDomainDimChoice) {
			iRungSD--;
			}
		}

	if (msr->nActive < msr->N*msr->param.dFracNoDomainDecomp) {
		if (msr->bDoneDomainDecomp && msr->iLastRungDomainDecomp >= iRungDD) {
			if (msr->param.bVRungStat) printf("Skipping Root Finder (nActive = %d/%d, iRung %d/%d/%d)\n",msr->nActive,msr->N,iRung,iRungDD,msr->iLastRungDomainDecomp);
			in.bDoRootFind = 0;
			in.bDoSplitDimFind = 0;
			}
		else if (iRungDD < iRung) {
			/* Set up the DD for the highest rung that still gets one */
			msrActiveRung(msr,iRungDD,bGreater);
			}
		}
	else iRungDD = iRung;

	if (in.bDoRootFind && msr->bDoneDomainDecomp && iRungDD > iRungSD && msr->iLastRungDomainDecomp >= iRungSD) {
		if (msr->param.bVRungStat) printf("Skipping Split Dim Finding (nDDActive = %d/%d, iRung %d/%d/%d/%d)\n",msr->nActive,msr->N,iRung,iRungDD,iRungSD,msr->iLastRungDomainDecomp);
		in.bDoSplitDimFind = 0;
		}

	if (msr->param.bGasDomainDecomp) {
		pstGasWeight(msr->pst,NULL,0,NULL,NULL);
		}

	if (msr->param.bRungDD) {
	        struct inRungDDWeight inRDD;
		inRDD.iMaxRung = msr->iCurrMaxRung;
                inRDD.dWeight = msr->param.dRungDDWeight;
		pstRungDDWeight(msr->pst,&inRDD,sizeof(struct inRungDDWeight),NULL,NULL);
		}

	if (msr->param.bVDetails) printf("DD: nActive %d nTreeActive %d nSmoothActive %d\n",msr->nActive,msr->nTreeActive,msr->nSmoothActive);
	LOGTIME( pstDomainDecomp(msr->pst,&in,sizeof(in),NULL,NULL), "Domain Decomposition", TIMING_DD );
	msr->bDoneDomainDecomp = 1; 

	msr->iLastRungDomainDecomp = iRungDD;
	if (iRungDD < iRung) {
	        /* Restore Active data */
		msrActiveRung(msr,iRung,bGreater);
		}
	}

void msrBuildTree(MSR msr,int bTreeActiveOnly, double dMass,int bSmooth)
{
	struct inBuildTree in;
	struct outBuildTree out;
	struct inColCells inc;
	struct ioCalcRoot root;
	KDN *pkdn;
	int iDum,nCell;

	if (msr->param.bVDetails) printf("Building local trees...\n");

	/*
	 ** First make sure the particles are in (Tree) Active/Inactive order.
	 */
	msrActiveTypeOrder(msr, TYPE_ACTIVE|TYPE_TREEACTIVE );
	in.nBucket = msr->param.nBucket;
	in.iOpenType = msr->iOpenType;
	in.iOrder = (msr->param.iOrder >= msr->param.iEwOrder)?
		msr->param.iOrder:msr->param.iEwOrder;
	in.dCrit = msr->dCrit;

	in.bActiveOnly = bTreeActiveOnly;
	in.bTreeActiveOnly = bTreeActiveOnly;
	if (bSmooth) {
		in.bBinary = 0;
		in.bGravity = 0;
		msr->iTreeType = MSR_TREE_DENSITY;
		msr->bGravityTree = 0;
		LOGTIME( pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum), "Tree built", TIMING_SPHTree );
		}
	else {
		in.bBinary = msr->param.bBinary;
		in.bGravity = 1;
		msr->bGravityTree = 1;
		if (msr->param.bBinary) {
			msr->iTreeType = MSR_TREE_SPATIAL;
			}
		else {
			msr->iTreeType = MSR_TREE_DENSITY;
			}
		LOGTIME( pstBuildTree(msr->pst,&in,sizeof(in),&out,&iDum), "Tree built", TIMING_GravTree );
		}
	msrMassCheck(msr,dMass,"After pstBuildTree in msrBuildTree");

	nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
	pkdn = malloc(nCell*sizeof(KDN));
	assert(pkdn != NULL);
	inc.iCell = ROOT;
	inc.nCell = nCell;
	pstColCells(msr->pst,&inc,sizeof(inc),pkdn,NULL);
	msrMassCheck(msr,dMass,"After pstColCells in msrBuildTree");

	pstDistribCells(msr->pst,pkdn,nCell*sizeof(KDN),NULL,NULL);
	msrMassCheck(msr,dMass,"After pstDistribCells in msrBuildTree");
	free(pkdn);
	if (!bSmooth) {
		pstCalcRoot(msr->pst,NULL,0,&root,&iDum);
		msrMassCheck(msr,dMass,"After pstCalcRoot in msrBuildTree");
		pstDistribRoot(msr->pst,&root,sizeof(struct ioCalcRoot),NULL,NULL);
		msrMassCheck(msr,dMass,"After pstDistribRoot in msrBuildTree");
	    }
    }


void msrDomainColor(MSR msr)
{
	pstDomainColor(msr->pst,NULL,0,NULL,NULL);
	}

#ifdef FAST_AGGS
void msrAggsCacheInit(MSR msr)	//JVD Local reorder rewrite
{
	if (msr->param.bVDetails) {
		printf("Initializing FAST_AGGS cache line ints...\n");
		pstAggsCacheInit(msr->pst,NULL,0,NULL,NULL);
		printf("FAST_AGGS cache line initialized\n");
		}
	else {
		pstAggsCacheInit(msr->pst,NULL,0,NULL,NULL);
		}
	}

void msrAggsReorder(MSR msr)
{
	/*
	** Reorder by decreasing original index. Needed for FAST_AGGS
	** compile option to allow aggs of varying sphere size.
	** Patterned after msrReorder().
	*/
	if (msr->param.bVDetails) {
		double sec,dsec;
		printf("Ordering aggs...\n");
		sec = msrTime(msr);
		pstLocalAggOrder(msr->pst,NULL,0,NULL,NULL);
		dsec = msrTime(msr) - sec;
		printf("Local agg orders established, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstLocalAggOrder(msr->pst,NULL,0,NULL,NULL);
		}

	/* Tree is void & domain decomp needs to be redone */
	msr->iTreeType = MSR_TREE_NONE;
	msr->bDoneDomainDecomp = 0;
	}
#endif /* FAST_AGGS */


void msrReorder(MSR msr)
{
	struct inDomainOrder in;

	in.iMaxOrder = msrMaxOrder(msr);
	if (msr->param.bVDetails) {
		double sec,dsec;
		printf("Ordering...\n");
		sec = msrTime(msr);
		pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
		pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
		dsec = msrTime(msr) - sec;
		printf("Order established, Wallclock: %f secs\n\n",dsec);
		}
	else {
		pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
		pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
		}
 	/*
	 ** Mark tree as void.
	 */
	msr->iTreeType = MSR_TREE_NONE;
 	/*
	 ** Mark domain decomp as not done.
	 */
	msr->bDoneDomainDecomp = 0;
	}


void msrCreateAllStepZeroOutputList(MSR msr, int *iNumOutputs, int OutputList[])
{
    /* Do all the stuff smoothed over all particles. */
    *iNumOutputs = 0;
    OutputList[(*iNumOutputs)++]=OUT_ACCELG_VECTOR;
    OutputList[(*iNumOutputs)++]=OUT_POT_ARRAY;
    OutputList[(*iNumOutputs)++]=OUT_DT_ARRAY;
    if (msrDoDensity(msr)) OutputList[(*iNumOutputs)++]=OUT_DENSITY_ARRAY;
}

void msrCreateGasStepZeroOutputList(MSR msr, int *iNumOutputs, int OutputList[])
{
    /* Do all the stuff smoothed over all particles. */
    *iNumOutputs = 0;
#ifdef GASOLINE				
    if (msr->param.bDoSphhOutput) OutputList[(*iNumOutputs)++]=OUT_SPHH_ARRAY;
    if (msr->param.bSphStep) OutputList[(*iNumOutputs)++]=OUT_SPHDT_ARRAY;
    if (!msr->param.bBulkViscosity){
        OutputList[(*iNumOutputs)++]=OUT_BALSARASWITCH_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_DIVV_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_MUMAX_ARRAY;
        if (msr->param.bShockTracker) {
            OutputList[(*iNumOutputs)++]=OUT_SHOCKTRACKER_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_DIVONCONH_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_DIVONCONX_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_DIVRHOV_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_GRADRHO_VECTOR;
            OutputList[(*iNumOutputs)++]=OUT_ACCELPRES_VECTOR;
        }

        OutputList[(*iNumOutputs)++]=OUT_ACCEL_VECTOR;
        OutputList[(*iNumOutputs)++]=OUT_PDV_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_PDVPRES_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_PDVVISC_ARRAY;
        }
#ifndef NOCOOLING				
    {
    int ArrayCnt = 0;
    char OutSuffix[20];
    int OutType;
    
    for (;;) {	
        CoolOutputArray( &msr->param.CoolParam, ArrayCnt, &OutType, OutSuffix );
        if (OutType == OUT_NULL) break;
        OutputList[(*iNumOutputs)++]=OutType;
        ArrayCnt++;
        }
    }
#endif
#endif
    
}

void msrCreateAllOutputList(MSR msr, int (*iNumOutputs), int OutputList[])
{
    /* Do all the stuff smoothed over all particles. */
    (*iNumOutputs) = 0;
    if (msrDoDensity(msr))  OutputList[(*iNumOutputs)++]=OUT_DENSITY_ARRAY;
    if (msr->param.bDoSoftOutput) OutputList[(*iNumOutputs)++]=OUT_SOFT_ARRAY;
    if (msr->param.bDohOutput) OutputList[(*iNumOutputs)++]=OUT_H_ARRAY;
}

void msrCreateGasOutputList(MSR msr, int (*iNumOutputs), int OutputList[])
{
    /* Add your new output file to the list after you've added
     * your item to the enumerated list in outtype.h, what the
     * value is in outtype.c ArrType or VecType and what the 
     * postfix is in outtype.c ArrFilename or VecFilename
     */
    (*iNumOutputs) = 0;

    if(msr->param.iBinaryOutput ==6)  {
        OutputList[(*iNumOutputs)++]=OUT_POS_VECTOR;
        OutputList[(*iNumOutputs)++]=OUT_VEL_VECTOR;
        OutputList[(*iNumOutputs)++]=OUT_MASS_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_POT_ARRAY;
#ifdef GASOLINE				
        OutputList[(*iNumOutputs)++]=OUT_GASDENSITY_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_TEMP_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_H_ARRAY;
#endif
        }
    else OutputList[(*iNumOutputs)++]=BIG_FILE; /*Tipsy, SS or whatever*/
    if(msr->param.bDoIOrderOutput) OutputList[(*iNumOutputs)++]=OUT_IORDER_ARRAY;
    if (msr->param.bDodtOutput) OutputList[(*iNumOutputs)++]=OUT_DT_ARRAY;
#ifdef GASOLINE				
    if (msr->param.bDoSphhOutput) OutputList[(*iNumOutputs)++]=OUT_SPHH_ARRAY;
#ifdef PDVDEBUG
    OutputList[(*iNumOutputs)++]=OUT_PDVPRES_ARRAY;
    OutputList[(*iNumOutputs)++]=OUT_PDVVISC_ARRAY;
#endif
    if (msr->param.bShockTracker) {
        OutputList[(*iNumOutputs)++]=OUT_SHOCKTRACKER_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_BALSARASWITCH_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_SPHH_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_DIVV_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_DIVRHOV_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_GRADRHO_VECTOR;
    }
#ifndef NOCOOLING				
    {
    int ArrayCnt = 0;
    char OutSuffix[20];
    int OutType;
    
    for (;;) {	
        CoolOutputArray( &msr->param.CoolParam, ArrayCnt, &OutType, OutSuffix );
        if (OutType == OUT_NULL) break;
        OutputList[(*iNumOutputs)++]=OutType;
        ArrayCnt++;
        }
    }
#endif

#ifdef STARFORM
    if(msr->param.bStarForm || msr->param.bFeedBack) {
        OutputList[(*iNumOutputs)++]=OUT_IGASORDER_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_COOLTURNONTIME_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_OXYGENMASSFRAC_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_IRONMASSFRAC_ARRAY;
        if(msr->param.bFormOutputs){
            OutputList[(*iNumOutputs)++]=OUT_TIMEFORM_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_MASSFORM_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_DENSITYFORM_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_TEMPFORM_ARRAY;
            OutputList[(*iNumOutputs)++]=OUT_RFORM_VECTOR;
            OutputList[(*iNumOutputs)++]=OUT_VFORM_VECTOR;
            }
#ifdef SIMPLESF
        OutputList[(*iNumOutputs)++]=OUT_DIVV_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_TCOOLAGAIN_ARRAY;
        OutputList[(*iNumOutputs)++]=OUT_MSTAR_ARRAY;
#endif
        }
#endif
#endif
#ifdef DEM_DIAG
	OutputList[(*iNumOutputs)++] = OUT_DEM_DIAG;
#endif /* DEM_DIAG */
#ifdef ORIGIN_HISTOGRAM
	OutputList[(*iNumOutputs)++] = OUT_ORIGIN_BINS;
#endif /* ORIGIN_HISTOGRAM */ 
}

void msrWriteNCOutputs(MSR msr, char *achFile, int OutputList[], int iNumOutputs, double dTime)
{
    FILE *fp;
    char dirname[256];
    int i, k, iDim, nDim, code, magic=1062053;
    LCL *plcl = msr->pst->plcl;
    char achOutFile[PST_FILENAME_SIZE];
    const char *typenames[3];
    int nTypes[3];
    struct inOutput inOut;
    struct outNC out;
    XDR xdrs;

#ifdef GASOLINE
    inOut.duTFac = (msr->param.dConstGamma - 1)*msr->param.dMeanMolWeight/
		msr->param.dGasConst;
#else
    inOut.duTFac = 1.0;
#endif
    /*DEBUG not used preminmax = 4*sizeof(int)+sizeof(double);*/
    typenames[0]="gas";
    typenames[1]="dark";
    typenames[2]="star";
    /*
     ** Calculate where to start writing.
     ** This sets plcl->n(Gas|Dark|Start)WriteStart.
     */
    msrCalcNCWriteStart(msr);
    /*
     ** Add Data Subpath for local and non-local names.
     */
    _msrMakePath(msr->param.achDataSubPath,achFile,inOut.achOutFile);
    /*
     ** Add local Data Path only for writing from this function.
     ** Other nodes may have different directory structures.
     */
    _msrMakePath(plcl->pszDataPath,inOut.achOutFile,achOutFile);
    
    assert(mkdir(achOutFile, 0775)<1);
/*    sprintf(xmlFile,"%s/description.xml",achOutFile);
    xmlfp = fopen(xmlFile,"w");
    fprintf(xmlfp,"<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n<simulation>\n");*/
    sprintf(dirname,"%s/gas",achOutFile);
    if (msr->nGas) assert(mkdir(dirname,0775)<1);
    sprintf(dirname,"%s/dark",achOutFile);
    if (msr->nDark) assert(mkdir(dirname,0775)<1);
    sprintf(dirname,"%s/star",achOutFile);
    if (msr->nStar) assert(mkdir(dirname,0775)<1);

    for (i=0; i<iNumOutputs;i++){
        code = FLOAT32;
        nTypes[0] = msr->nGas;nTypes[1] = msr->nDark;nTypes[2] = msr->nStar;
        switch (OutputList[i]){
            case OUT_TIMEFORM_ARRAY:
            case OUT_MASSFORM_ARRAY:
            case OUT_RFORM_VECTOR:
            case OUT_VFORM_VECTOR:
            case OUT_DENSITYFORM_ARRAY:
            case OUT_TEMPFORM_ARRAY:
                nTypes[0]=nTypes[1]=0;
                break;
            case OUT_IGASORDER_ARRAY:
                nTypes[0]=nTypes[1]=0; /*FALLTHRU*/
            case OUT_IORDER_ARRAY:
                code=INT32;
                break;
            /* Gas only floats*/
            case OUT_COOLTURNONTIME_ARRAY:
            case OUT_COOL_ARRAY0:
            case OUT_COOL_ARRAY1:
            case OUT_COOL_ARRAY2:
            case OUT_SPHH_ARRAY:
            case OUT_TEMP_ARRAY:
            case OUT_GASDENSITY_ARRAY:
            case OUT_PDVPRES_ARRAY:
            case OUT_PDVVISC_ARRAY:
                nTypes[1]=nTypes[2]=0;
                break;
            case OUT_OXYGENMASSFRAC_ARRAY:
            case OUT_IRONMASSFRAC_ARRAY:
            case OUT_METALS_ARRAY:
                nTypes[1]=0;
                break;
            }
            
        nDim = (OutputList[i] > OUT_1D3DSPLIT) ? 3 : 1;
        inOut.iBinaryOutput = msr->param.iBinaryOutput;
        inOut.N = msr->N;
        inOut.iType=OutputList[i];
        pstOutNCVector(msr->pst,&inOut,sizeof(inOut),&out,NULL);
        for (k=0;k<3;k++){
            _msrMakePath(plcl->pszDataPath,inOut.achOutFile,achOutFile);
            if (nTypes[k]) {
                sprintf(achOutFile,"%s/%s/",achOutFile,typenames[k]);
                VecFilename(achOutFile,OutputList[i]);
                fp = fopen(achOutFile,"r+");
                assert(fp != NULL);
                xdrstdio_create(&xdrs,fp,XDR_ENCODE);
                xdr_int(&xdrs,&magic);
                xdr_double(&xdrs,&dTime);
                xdr_int(&xdrs,&nTypes[k]);
                xdr_int(&xdrs,&nDim);
                xdr_int(&xdrs,&code);
                for (iDim=0; iDim<nDim; iDim++) 
                    xdr_float(&xdrs,&out.min[k][iDim]);
                for (iDim=0; iDim<nDim; iDim++) 
                    xdr_float(&xdrs,&out.max[k][iDim]);
                xdr_destroy(&xdrs);
                fclose(fp);
                }
            }
        }

    }

void msrWriteOutputs(MSR msr, char *achFile, int OutputList[], int iNumOutputs, double dTime)
{
    FILE *fp;
    int i, iDim, nDim;
    LCL *plcl = msr->pst->plcl;
    char achOutFile[PST_FILENAME_SIZE];
    struct inOutput inOut;
#ifdef COLLISIONS
    struct inWriteSS in;
#else
    struct inWriteTipsy in;
#endif

    if (msr->param.iBinaryOutput == 6) {
        msrWriteNCOutputs(msr, achFile, OutputList, iNumOutputs, dTime);
        return;
        }

    /*
     ** Calculate where to start writing.
     ** This sets plcl->nWriteStart.
     */
    msrCalcWriteStart(msr);
    /*
     ** Add Data Subpath for local and non-local names.
     */
    _msrMakePath(msr->param.achDataSubPath,achFile,in.achOutFile);
    /*
     ** Add local Data Path only for writing from this function.
     ** Other nodes may have different directory structures.
     */
    _msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

#ifdef COLLISIONS
	in.bReduced = 0; /* standard */
#endif

#ifdef SPRINGS
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.bNoSpringsOutput = msr->param.FO.SP.bNoSpringsOutput;
#endif /* SPRINGS */

#ifdef DEM
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.nDEMOutputs = msr->param.FO.DP.nDEMOutputs;
#endif /* DEM */

    /* Write Headers */

    snprintf(inOut.achOutFile,PST_FILENAME_SIZE,"%s.",in.achOutFile);
    for (i=0;i<iNumOutputs;i++){
        if ( OutputList[i] == BIG_FILE ){
#ifdef COLLISIONS
            msrWriteSSHead(msr,achOutFile,dTime,in.bReduced);
#else
            msrWriteTipsyHead(msr,achOutFile,dTime,&in);
#endif
            } 
        else {
            _msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
			strcat(achOutFile,".");
            VecFilename(achOutFile,OutputList[i]);
            fp = fopen(achOutFile,"w");
            assert(fp != NULL);
            if (msr->param.iBinaryOutput && msr->param.bStandard) {
                XDR xdrs;
                xdrstdio_create(&xdrs,fp,XDR_ENCODE);
                xdr_int(&xdrs,&msr->N);
                }
            else if (msr->param.iBinaryOutput) {
				fwrite(&msr->N,sizeof(int),1,fp);
				}
            else {
				fprintf(fp,"%d\n",msr->N);
#ifdef DEM_DIAG /*DEBUG a bit of a hack unfortunately*/
				fprintf(fp, "%i", MAX_NUM_OVERLAPS_PER_PARTICLE);
#ifdef WALLS
				fprintf(fp, " %i", MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS);
#endif /* WALLS */
				fprintf(fp, "\n");
#endif /* DEM_DIAG */
				}
            fclose(fp);
            }
        }

    /* Write Data */
    inOut.iBinaryOutput = msr->param.iBinaryOutput;
    inOut.N = msr->N;
#ifndef SSIO_USE_MPI
    if (msr->param.iBinaryOutput) {
#endif
        if(msr->param.bParaWrite) {
            for (i=0; i<iNumOutputs;i++){
                if ( OutputList[i] == BIG_FILE ){
#ifdef COLLISIONS
                    pstWriteSS(msr->pst,&in,sizeof(in),NULL,NULL);
#else
                    pstWriteTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
#endif
                    } 
                else {
                    inOut.iType=OutputList[i];
                    if ((OutputList[i] > OUT_1D3DSPLIT)&& msr->param.bPackedVector) {
                        inOut.iDim = -3;
                        pstOutVector(msr->pst,&inOut,sizeof(inOut),NULL,NULL);
                    }
                    else
                    {
                        nDim=(OutputList[i] > OUT_1D3DSPLIT) ? 3 : 1;
                        for (iDim=0; iDim<nDim; iDim++) {
                            inOut.iDim = iDim;
                            pstOutVector(msr->pst,&inOut,sizeof(inOut),NULL,NULL);
                        }
                    }
                }
            }
        }
        else /* Serial Binary */
        {
            msrOneNodeWriteOutputs(msr, OutputList, iNumOutputs, &in);
        }
#ifndef SSIO_USE_MPI
    }
    else  /* ASCII:  NO PARALLEL OPTION! Only packed vectors supported. */
    {
        msrOneNodeWriteOutputs(msr, OutputList, iNumOutputs, &in);
    }
#endif
}
    
void msrOneNodeWriteOutputs(MSR msr, int OutputList[], int iNumOutputs,
#ifdef COLLISIONS
							struct inWriteSS *in
#else
							struct inWriteTipsy *in
#endif
							)
{
    int i,id,iDim,nDim;
    int iOut;  /* iterator through OutputList */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    char achOutFileVec[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    _msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);
	snprintf(achOutFileVec, PST_FILENAME_SIZE, "%s.", achOutFile);

#ifdef SPRINGS
	plcl->pkd->bNoSpringsOutput = in->bNoSpringsOutput;
#endif
#ifdef DEM
	plcl->pkd->nDEMOutputs = in->nDEMOutputs;
#endif

    /* 
     * First write our own particles.
     */
    assert(msr->pMap[0] == 0);
    for (iOut=0; iOut<iNumOutputs;iOut++){
        if( OutputList[iOut]== BIG_FILE){

#ifdef COLLISIONS
            pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart,in->bReduced);
#else
            pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,in->bStandard,
                                      in->dvFac,in->duTFac,in->iGasModel); 
#endif
            } else {
            /* Only packed ASCII format supported!
             * Use readpackedvector in Tipsy */
			  if ((OutputList[iOut] > OUT_1D3DSPLIT) && (!msr->param.iBinaryOutput || msr->param.bPackedVector)) {
                pkdOutVector(plcl->pkd,achOutFileVec,plcl->nWriteStart, -3,
							 OutputList[iOut], msr->param.iBinaryOutput,msr->N,
							 msr->param.bStandard);
                } else {
                nDim=(OutputList[iOut] > OUT_1D3DSPLIT) ? 3 : 1;
                for (iDim=0; iDim<nDim; iDim++) 
				  pkdOutVector(plcl->pkd,achOutFileVec,plcl->nWriteStart,
							   iDim, OutputList[iOut],msr->param.iBinaryOutput,
							   msr->N,msr->param.bStandard);
                }
            }
        }

	nStart = plcl->pkd->nLocal;
    /* Write out the particles on all the other nodes */
    for (i=1;i<msr->nThreads;++i) {
        id = msr->pMap[i];
        /* 
         * Swap particles with the remote processor.
         */
        inswap = 0;
        mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
        pkdSwapAll(plcl->pkd, id);
        mdlGetReply(pst0->mdl,id,NULL,NULL);
        /* 
         * Write the swapped particles.
         */
        for (iOut=0; iOut<iNumOutputs;iOut++){
            if( OutputList[iOut]== BIG_FILE){
#ifdef COLLISIONS
                pkdWriteSS(plcl->pkd,achOutFile,nStart,in->bReduced);
#else
                pkdWriteTipsy(plcl->pkd,achOutFile,nStart,in->bStandard,
							  in->dvFac,in->duTFac,in->iGasModel); 
#endif
                } else {
                /* Only packed ASCII format supported!
                 * Use readpackedvector in Tipsy */
                if ((OutputList[iOut] > OUT_1D3DSPLIT) && (!msr->param.iBinaryOutput || msr->param.bPackedVector)) {
				  pkdOutVector(plcl->pkd,achOutFileVec, nStart,
							   -3, OutputList[iOut], msr->param.iBinaryOutput,
							   msr->N,msr->param.bStandard);
				  } else {
                  nDim=(OutputList[iOut] > OUT_1D3DSPLIT) ? 3 : 1;
                  for (iDim=0; iDim<nDim; iDim++) {
#ifdef SIMPLESF
					  /* The SF variables are written in binary no
						 matter what */
					  pkdOutVector(plcl->pkd,achOutFileVec,
								   nStart, iDim, OutputList[iOut],
								   ((OutputList[iOut]==OUT_TCOOLAGAIN_ARRAY || OutputList[iOut]==OUT_MSTAR_ARRAY) ? 1 : msr->param.iBinaryOutput),
								   msr->N,msr->param.bStandard);
#else
					  pkdOutVector(plcl->pkd,achOutFileVec, nStart,
								   iDim, OutputList[iOut],
								   msr->param.iBinaryOutput, msr->N,
								   msr->param.bStandard); 
#endif
                      }
                   }
                }
            }
        nStart += plcl->pkd->nLocal;
        /* 
         * Swap them back again.
         */
        inswap = 0;
        mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
        pkdSwapAll(plcl->pkd, id);
        mdlGetReply(pst0->mdl,id,NULL,NULL);
        }
    assert(nStart == msr->N);

    }
    
void msrOneNodeOutArray(MSR msr, struct inOutput *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    _msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    assert(msr->pMap[0] == 0);
    nStart = plcl->pkd->nLocal;
    pkdOutVector(plcl->pkd,in->achOutFile,nStart, 0, in->iType,in->iBinaryOutput, msr->N,in->bStandard); 
    for (i=1;i<msr->nThreads;++i) {
            id = msr->pMap[i];
        /* 
         * Swap particles with the remote processor.
         */
        inswap = 0;
        mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
        pkdSwapAll(plcl->pkd, id);
        mdlGetReply(pst0->mdl,id,NULL,NULL);
        /* 
         * Write the swapped particles.
         */
#ifdef SIMPLESF
        pkdOutVector(plcl->pkd,in->achOutFile, nStart, 0, in->iType, ((in->iType==OUT_TCOOLAGAIN_ARRAY || in->iType==OUT_MSTAR_ARRAY) ? 1 : in->iBinaryOutput), msr->N,in->bStandard);
#else
        pkdOutVector(plcl->pkd,in->achOutFile,nStart, 0, in->iType,in->iBinaryOutput, msr->N,in->bStandard); 
#endif
        nStart += plcl->pkd->nLocal;
        /* 
         * Swap them back again.
         */
        inswap = 0;
        mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
        pkdSwapAll(plcl->pkd, id);
        mdlGetReply(pst0->mdl,id,NULL,NULL);
        }
    assert(nStart == msr->N);
    }

void msrOutArray(MSR msr,const char *pszFile,int iType)
{
	FILE *fp;
	struct inOutput in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;

	/*
	 ** Calculate where to start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);

	if (pszFile) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Array Output File:%s\n",achOutFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No Array Output File specified\n");
		_msrExit(msr,1);
		return;
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	in.iType = iType;
	in.iBinaryOutput = msr->param.iBinaryOutput;
        in.iDim=1;
        in.N = msr->N;
	if (msr->param.iBinaryOutput) {
		fwrite(&msr->N,sizeof(int),1,fp);
                fclose(fp);
                if(msr->param.bParaWrite)
                    pstOutArray(msr->pst,&in,sizeof(in),NULL,NULL);
                else
                    msrOneNodeOutArray(msr, &in);
		}
	else {
		fprintf(fp,"%d\n",msr->N);
                fclose(fp);
                msrOneNodeOutArray(msr, &in);
		}
	}


void msrOneNodeOutVector(MSR msr, struct inOutput *in)
{
    int i,id, iDim;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
    _msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    assert(msr->pMap[0] == 0);
    if (in->iBinaryOutput || msr->param.bPackedVector) {
        nStart = plcl->pkd->nLocal;
        /* 
         * First write our own particles.
         */
        if (msr->param.bPackedVector) {
            pkdOutVector(plcl->pkd,in->achOutFile,nStart, -3, in->iType,in->iBinaryOutput, msr->N,in->bStandard); 
            } else {
            for (iDim=0;iDim<3;++iDim) {
                pkdOutVector(plcl->pkd,achOutFile, nStart, iDim, in->iType,msr->param.iBinaryOutput, msr->N,in->bStandard);
                }
            }
            
        for (i=1;i<msr->nThreads;++i) {
            id = msr->pMap[i];
            /* 
             * Swap particles with the remote processor.
             */
            inswap = 0;
            mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
            pkdSwapAll(plcl->pkd, id);
            mdlGetReply(pst0->mdl,id,NULL,NULL);
            /* 
             * Write the swapped particles.
             */
           if (msr->param.bPackedVector) {
                pkdOutVector(plcl->pkd,in->achOutFile,nStart,-3,in->iType,in->iBinaryOutput, msr->N,in->bStandard); 
                } else {
                for (iDim=0;iDim<3;++iDim) {
                    pkdOutVector(plcl->pkd,achOutFile,nStart, iDim,in->iType,msr->param.iBinaryOutput, msr->N,in->bStandard);
                    }
                }
            nStart += plcl->pkd->nLocal;
            /* 
             * Swap them back again.
             */
            inswap = 0;
            mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
            pkdSwapAll(plcl->pkd, id);
            mdlGetReply(pst0->mdl,id,NULL,NULL);
            }
        assert(nStart == msr->N);
        } else { /* ASCII, non packed vectors are a pain! 
                  * Note all the swap all's.  You should
                  * definitely use packed vectors or binary
                  * format for large simulations!
                  */
            nStart = 0;
            for (iDim=0;iDim<3;++iDim) {
                pkdOutVector(plcl->pkd,achOutFile,nStart, iDim,in->iType,msr->param.iBinaryOutput, msr->N,in->bStandard);
                nStart += plcl->pkd->nLocal;
                for (i=1;i<msr->nThreads;++i) {
                    id = msr->pMap[i];
                    inswap = 0;
                    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
                    pkdSwapAll(plcl->pkd, id);
                    mdlGetReply(pst0->mdl,id,NULL,NULL);

                    pkdOutVector(plcl->pkd,achOutFile,nStart, iDim,in->iType,msr->param.iBinaryOutput, msr->N,in->bStandard);
                    
                    nStart += plcl->pkd->nLocal;
                    /* 
                     * Swap them back again.
                     */
                    inswap = 0;
                    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
                    pkdSwapAll(plcl->pkd, id);
                    mdlGetReply(pst0->mdl,id,NULL,NULL);
                    }
                }
            assert(nStart == 3*msr->N);
            }
    }

void msrOutVector(MSR msr,char *pszFile,int iType)
{
	struct inOutput in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	FILE *fp;
	int iDim;

        msrCalcWriteStart(msr);

	if (pszFile) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

		fp = fopen(achOutFile,"w");
		if (!fp) {
			printf("Could not open Vector Output File:%s\n",achOutFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No Vector Output File specified\n");
		_msrExit(msr,1);
		return;
		}
	/*
	 ** Write the Header information and close the file again.
	 */
	in.iType = iType;
	in.iBinaryOutput = msr->param.iBinaryOutput;
	if (msr->param.iBinaryOutput) {
            fwrite(&msr->N,sizeof(int),1,fp);
            fclose(fp);
            in.N=msr->N;
            if(msr->param.bParaWrite){
                if (msr->param.bPackedVector) {
                    in.iDim = -3;
                    pstOutVector(msr->pst,&in,sizeof(in),NULL,NULL);
                    } else {
                    for (iDim=0;iDim<3;++iDim) {
                        in.iDim = iDim;
                        pstOutVector(msr->pst,&in,sizeof(in),NULL,NULL);
                        }
                    }
                } else msrOneNodeOutVector(msr, &in);
            } else {
		fprintf(fp,"%d\n",msr->N);
                fclose(fp);
                msrOneNodeOutVector(msr, &in);
		}
	}


void msrSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric)
{
	struct inSmooth in;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
						   4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
	{
	    double dAccFac = 1.0/(in.smf.a*in.smf.a*in.smf.a);
	    in.smf.dDeltaAccelFac = msr->param.dEtaDeltaAccel/sqrt(dAccFac);
	    }
	in.smf.dBHSinkAlphaFactor = msr->param.dBHSinkAlpha*4*M_PI;
	in.smf.dBHSinkEddFactor = msr->param.dBHSinkEddEff;
	in.smf.dBHSinkFeedbackFactor = msr->param.dBHSinkFeedbackFactor;
	in.smf.dSinkCurrentDelta = msr->param.dSinkCurrentDelta;
	in.smf.bSinkThermal = msr->param.bSinkThermal;
	in.smf.dSinkRadius = msr->param.dSinkRadius;
	in.smf.dSinkBoundOrbitRadius = msr->param.dSinkBoundOrbitRadius;
	in.smf.iSmoothFlags = 0; /* Initial value, return value in outSmooth */
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
#if defined(STARFORM) || defined(CHECKSOFT)
	in.smf.dTime = dTime;
#endif
#ifdef STARFORM
        in.smf.dSecUnit = msr->param.dSecUnit;  /*if you want to output feedback shutoff time in years*/
        in.smf.dGmUnit = msr->param.dMsolUnit*MSOLG;  /*if you want to use snCalcSNIIFeedback to calculate feedback*/
        in.smf.sn = *msr->param.sn;
        in.smf.dMinMassFrac = msr->param.stfm->dMinMassFrac;
	in.smf.dTime = dTime;
	in.smf.bSNTurnOffCooling = msr->param.bSNTurnOffCooling;
	in.smf.bSmallSNSmooth = msr->param.bSmallSNSmooth;
	in.smf.bShortCoolShutoff = msr->param.bShortCoolShutoff;
        /* from McKee and Ostriker (1977) ApJ 218 148 */
        in.smf.dRadPreFactor = pow(10,1.74)/(msr->param.dKpcUnit*1000.0)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),-0.16)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.2);
	if (msr->param.bShortCoolShutoff){        /* end of snowplow */
        	in.smf.dTimePreFactor = SECONDSPERYEAR*pow(10,5.92)/(msr->param.dSecUnit)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),0.27)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.64);
		} else {       /* t_{max}*/
        	in.smf.dTimePreFactor = SECONDSPERYEAR*pow(10,6.85)/(msr->param.dSecUnit)*
    pow(MSOLG*msr->param.dMsolUnit/(msr->param.dMeanMolWeight*MHYDR*pow(KPCCM*msr->param.dKpcUnit,3)),0.32)*
    pow(0.0001*GCGS*pow(MSOLG*msr->param.dMsolUnit,2)/(pow(KPCCM*msr->param.dKpcUnit,4)*KBOLTZ),-0.70);
		}
#endif /*STARFORM*/
#ifdef COLLISIONS
	/* for Hill sphere checks in FindRejects() */
	in.smf.dCentMass = msr->param.dCentMass;
#endif
	if (msr->param.bVStep) {
	    struct outSmooth out;
	    LOGTIME( pstSmooth(msr->pst,&in,sizeof(in),&out,NULL), "Smooth Calculated", TIMING_Smooth );
	    if (msr->nThreads > 1) {
		double iP = 1.0/msr->nThreads;
		printf("Particle Cache Statistics (average per processor):\n");
		printf("    Accesses:    %10g\n",out.dpASum*iP);
		printf("    Miss Ratio:  %10g\n",out.dpMSum*iP);
		printf("    Min Ratio:   %10g\n",out.dpTSum*iP);
		printf("    Coll Ratio:  %10g\n",out.dpCSum*iP);
		printf("Cell Cache Statistics (average per processor):\n");
		printf("    Accesses:    %10g\n",out.dcASum*iP);
		printf("    Miss Ratio:  %10g\n",out.dcMSum*iP);
		printf("    Min Ratio:   %10g\n",out.dcTSum*iP);
		printf("    Coll Ratio:  %10g\n",out.dcCSum*iP);
		printf("\n");
		}
	    }
	else {
	    LOGTIME( pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL), "Smooth Calculated", TIMING_Smooth );
	    }
	}


void msrReSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric)
{
	struct inReSmooth in;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = iSmoothType;
	in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
			       4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
	{
	    double dAccFac = 1.0/(in.smf.a*in.smf.a*in.smf.a);
	    in.smf.dDeltaAccelFac = msr->param.dEtaDeltaAccel/sqrt(dAccFac);
	    }
	in.smf.dBHSinkAlphaFactor = msr->param.dBHSinkAlpha*4*M_PI;
	in.smf.dBHSinkEddFactor = msr->param.dBHSinkEddEff;
	in.smf.dBHSinkFeedbackFactor = msr->param.dBHSinkFeedbackFactor;
	in.smf.dSinkCurrentDelta = msr->param.dSinkCurrentDelta;
	in.smf.bSinkThermal = msr->param.bSinkThermal;
	in.smf.dSinkRadius = msr->param.dSinkRadius;
	in.smf.dSinkBoundOrbitRadius = msr->param.dSinkBoundOrbitRadius;
	in.smf.iSmoothFlags = 0; /* Initial value, return value in outSmooth */
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
	if (msr->param.bVStep) {
		struct outSmooth out;

		LOGTIME( pstReSmooth(msr->pst,&in,sizeof(in),&out,NULL), "ReSmooth Calculated", TIMING_ReSmooth );
		if (msr->nThreads > 1) {
		    double iP = 1.0/msr->nThreads;
		    printf("Particle Cache Statistics (average per processor):\n");
		    printf("    Accesses:    %10g\n",out.dpASum*iP);
		    printf("    Miss Ratio:  %10g\n",out.dpMSum*iP);
		    printf("    Min Ratio:   %10g\n",out.dpTSum*iP);
		    printf("    Coll Ratio:  %10g\n",out.dpCSum*iP);
		    printf("Cell Cache Statistics (average per processor):\n");
		    printf("    Accesses:    %10g\n",out.dcASum*iP);
		    printf("    Miss Ratio:  %10g\n",out.dcMSum*iP);
		    printf("    Min Ratio:   %10g\n",out.dcTSum*iP);
		    printf("    Coll Ratio:  %10g\n",out.dcCSum*iP);
		    printf("\n");
		    }
		}
	else {
		LOGTIME( pstReSmooth(msr->pst,&in,sizeof(in),NULL,NULL), "ReSmooth Calculated", TIMING_ReSmooth );
		}
	}

void msrMarkSmooth(MSR msr,double dTime,int bSymmetric,int iMarkType)
{
	struct inMarkSmooth in;

	/*
	 ** Make sure that the type of tree is a density binary tree!
	 */
	assert(msr->iTreeType == MSR_TREE_DENSITY);
	in.nSmooth = msr->param.nSmooth;
	in.bPeriodic = msr->param.bPeriodic;
	in.bSymmetric = bSymmetric;
	in.iSmoothType = SMX_MARK; 
	in.iMarkType = iMarkType;
	if (msrComove(msr)) {
		in.smf.H = csmTime2Hub(msr->param.csm,dTime);
		in.smf.a = csmTime2Exp(msr->param.csm,dTime);
		}
	else {
		in.smf.H = 0.0;
		in.smf.a = 1.0;
		}
#ifdef GASOLINE
	in.smf.alpha = msr->param.dConstAlpha;
	in.smf.beta = msr->param.dConstBeta;
	in.smf.gamma = msr->param.dConstGamma;
	in.smf.algam = in.smf.alpha*sqrt(in.smf.gamma*(in.smf.gamma - 1));
	in.smf.bGeometric = msr->param.bGeometric;
	in.smf.bCannonical = msr->param.bCannonical;
	in.smf.bGrowSmoothList = 0;
#endif
	LOGTIME( pstMarkSmooth(msr->pst,&in,sizeof(in),NULL,NULL), "MarkSmooth Calculated", TIMING_MarkSmooth );
	}

void msrUpdateSoft(MSR msr,double dTime) {
#ifdef CHANGESOFT
       if (!(msr->param.bPhysicalSoft || msr->param.bVariableSoft)) return;
       if (msr->param.bPhysicalSoft) {
	 struct inPhysicalSoft in;

	 in.dFac = 1./csmTime2Exp(msr->param.csm,dTime);
	 in.bSoftMaxMul = msr->param.bSoftMaxMul;
	 in.dSoftMax = msr->param.dSoftMax;

	 if (msr->param.bSoftMaxMul && in.dFac > in.dSoftMax) in.dFac = in.dSoftMax;

	 pstPhysicalSoft(msr->pst,&in,sizeof(in),NULL,NULL);
       }
       else {
		 int type;
		 struct inPreVariableSoft inPre;
		 struct inPostVariableSoft inPost;
		 type =
		   (msr->param.bVariableSoftDark ? TYPE_DARK : 0) |
		   (msr->param.bVariableSoftStar ? TYPE_STAR : 0) |
		   (msr->param.bVariableSoftGas ? TYPE_GAS : 0);
		 inPre.iVariableSoftType = type;
		 pstPreVariableSoft(msr->pst,&inPre,sizeof(inPre),NULL,NULL);
		 
		 if (msr->param.bSoftByType) {
		   if (msr->nDark && msr->param.bVariableSoftDark) {
			 msrActiveType(msr,TYPE_DARK,TYPE_TREEACTIVE);
			 msrBuildTree(msr,1,-1.0,1);
			 msrActiveExactType(msr,TYPE_ACTIVE|TYPE_DARK,TYPE_ACTIVE|TYPE_DARK,TYPE_SMOOTHACTIVE);
			 msrSmooth(msr,dTime,SMX_NULL,0);

		   }
		   if (msr->nGas && msr->param.bVariableSoftGas) {
#ifdef DENSSOFT
			 msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE|TYPE_DensZeroed);
			 msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			 msrBuildTree(msr,1,-1.0,1);
			 msrActiveType(msr,TYPE_ACTIVE,TYPE_DensACTIVE );
			 msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1);
#else
			 msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
			 msrBuildTree(msr,1,-1.0,1);
			 msrActiveExactType(msr,TYPE_ACTIVE|TYPE_GAS,TYPE_ACTIVE|TYPE_GAS,TYPE_SMOOTHACTIVE);
			 msrSmooth(msr,dTime,SMX_NULL,0);
#endif
		   }
		   if (msr->nStar && msr->param.bVariableSoftStar) {
			 msrActiveType(msr,TYPE_STAR,TYPE_TREEACTIVE);
			 msrBuildTree(msr,1,-1.0,1);
			 msrActiveExactType(msr,TYPE_ACTIVE|TYPE_STAR,TYPE_ACTIVE|TYPE_STAR,TYPE_SMOOTHACTIVE);
			 msrSmooth(msr,dTime,SMX_NULL,0);
		   }
		 }
		 else {
		   msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
		   msrBuildTree(msr,1,-1.0,1);
		   msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE);
		   msrSmooth(msr,dTime,SMX_NULL,0);
		 }

		 inPost.dSoftMax = msr->param.dSoftMax;
		 inPost.bSoftMaxMul = msr->param.bSoftMaxMul;
		 inPost.iVariableSoftType = type;
		 pstPostVariableSoft(msr->pst,&inPost,sizeof(inPost),NULL,NULL);
       }
#endif
}

#ifdef GR_DRAG

void msrGRDragGetSunAccel(MSR msr,double aSun[3])
{
	/*
	** Brute force calculation of the acceleration on the "Sun"
	** (central massive object) with GR terms included.  We do it this
	** way (overwriting aSun computed in msrGravity()), because the
	** tree does not know about GR.  This is only an order(N)
	** operation, so it should not slow things down too much.
	*/

	struct inGRDragGetSunAccel in;
	struct outGRDragGetSunAccel out;
	int k;

	in.dSunMass = msr->param.dCentMass;
	pstGRDragGetSunAccel(msr->pst,&in,sizeof(in),&out,NULL);
	for (k=0;k<3;k++)
		aSun[k] = out.aSun[k];
	}

void msrGRIntegrateCloseParticles(MSR msr,double dDelta,double dTime)
{
	/*
	** Identifies particles within msr->param.dGRRadius of the central
	** body and integrates their orbits forward an interval dDelta
	** using explicit treatment of GR.  The particles are set to
	** inactive so they are neither kicked nor drifted.
	*/

	/* Also corrects the aSun[] vector so that it is not affected by
	** a particle that's at a very low pericenter */

	struct inGRIntegrateCloseParticles in;
	struct outGRIntegrateCloseParticles out;
	struct inGRCorrectaSun in2;
	int j;

	in.dSunMass = msr->param.dCentMass;
	in.dDelta = dDelta;
	in.dTime = dTime;
	pstGRIntegrateCloseParticles(msr->pst,&in,sizeof(in),&out,NULL);
	if (out.nMerged > 0)
	  msrAddDelParticles(msr);
	msr->dMergerMassLost += out.dMergerMassLost;
	for (j = 0 ; j < 3 ; ++j)
	  in2.aSunCorrection[j] = out.aSunCorrection[j];
	pstGRCorrectaSun(msr->pst,&in2,sizeof(in2),NULL,NULL); 
}

#endif /*GR_DRAG*/

void msrGravity(MSR msr,double dStep,int bDoSun,
				int *piSec,double *pdWMax,double *pdIMax,double *pdEMax,
				int *pnActive)
{
	struct inGravity in;
	struct outGravity out;
	struct inGravExternal inExt;
	int iDum,j;
	double sec,dsec;

#ifdef AGGS_IN_PATCH
	msr->dTime = dStep*msr->param.dDelta;
#endif
	if (msr->param.bDoSelfGravity) {
		assert(msr->bGravityTree == 1);
		assert(msr->iTreeType == MSR_TREE_SPATIAL || 
			   msr->iTreeType == MSR_TREE_DENSITY);
		if (msr->param.bVStep) printf("Calculating Gravity, Step:%f\n",dStep);
		in.nReps = msr->param.nReplicas;
		in.bPeriodic = msr->param.bPeriodic;
		in.iOrder = msr->param.iOrder;
		in.bEwald = msr->param.bEwald;
		in.iEwOrder = msr->param.iEwOrder;
#ifdef COLLISIONS
		if (msr->param.CP.iOverlapOption == OverlapRepel) {
			in.bRepel = 1;
			in.dRepelFac = msr->param.CP.dRepelFac;
			}
		else {
			in.bRepel = 0;
			in.dRepelFac = 0.0;
			}
#endif
#ifdef SLIDING_PATCH
		in.dTime = dStep*msr->param.dDelta;
		in.PP = msr->param.PP; /* struct copy */
#endif
		/*
		 ** The meaning of 'bDoSun' here is that we want the accel on (0,0,0)
		 ** to use in creating the indirect acceleration on each particle. This
		 ** is why it looks inconsistent with call to pstGravExternal() below.
		 */
		in.bDoSun = msr->param.bHeliocentric;
		in.dEwCut = msr->param.dEwCut;
		in.dEwhCut = msr->param.dEwhCut;
		in.dSunSoft = msr->param.dSunSoft;
		sec = msrTime(msr);
		pstGravity(msr->pst,&in,sizeof(in),&out,&iDum);
		dsec = msrTime(msr) - sec;
		LOGTIMINGUPDATE( dsec, TIMING_Gravity );
#ifdef SPECIAL_PARTICLES
		if (msr->param.nSpecial) {
			/*
			 ** Handle contributions from "special" particles now,
			 ** e.g. oblateness effects, GR, etc.
			 */
			struct inGetSpecial inGetSpec;
			struct outGetSpecial outGetSpec;
			struct inDoSpecial inDoSpec;
			struct outDoSpecial outDoSpec;
			inGetSpec.nSpecial = msr->param.nSpecial;
			for (j=0;j<msr->param.nSpecial;j++)
				inGetSpec.iId[j] = msr->param.iSpecialId[j]; /* orig indices */
			inGetSpec.mInfo.dCentMass = msr->param.dCentMass; /* if special frame */
			pstGetSpecialParticles(msr->pst,&inGetSpec,sizeof(inGetSpec),
								   &outGetSpec,NULL);
			inDoSpec.nSpecial = msr->param.nSpecial;
			for (j=0;j<msr->param.nSpecial;j++) {
				inDoSpec.sData[j] = msr->param.sSpecialData[j]; /* struct cp */
				inDoSpec.sInfo[j] = outGetSpec.sInfo[j]; /* ditto */
				}
			inDoSpec.mInfo.bNonInertial = msr->param.bHeliocentric;
			inDoSpec.mInfo.dxPeriod = msr->param.dxPeriod;
			inDoSpec.mInfo.dyPeriod = msr->param.dyPeriod;
			inDoSpec.mInfo.dzPeriod = msr->param.dzPeriod;
			inDoSpec.mInfo.dOmega = msr->param.dOmega;
			inDoSpec.mInfo.dTime =  dStep*msr->param.dDelta;
			inDoSpec.mInfo.nReplicas = msr->param.nReplicas;
			if (inDoSpec.mInfo.bNonInertial)
				for (j=0;j<3;j++)
					outDoSpec.aFrame[j] = 0.0; /* initialize */
			pstDoSpecialParticles(msr->pst,&inDoSpec,sizeof(inDoSpec),
								  &outDoSpec,NULL);
			if (inDoSpec.mInfo.bNonInertial)
				for (j=0;j<3;j++)
					out.aSun[j] += outDoSpec.aFrame[j]; /* add non-inertial term */
			}
#endif
#ifdef GR_DRAG
		msrGRDragGetSunAccel(msr,out.aSun);
#endif /* GR_DRAG */
	    } /* if including particle self-gravity */
	else { /* pstGravity() not called, so zero all output parameters */
		dsec = out.nActive = out.nTreeActive = out.aSun[0] =
		out.aSun[1] = out.aSun[2] = out.dPartSum = out.dCellSum =
		out.dSoftSum = out.dFlop = out.dWSum = out.dWMax = out.dWMin =
		out.dISum = out.dIMax = out.dIMin = out.dESum = out.dEMax =
		out.dEMin = out.dpASum = out.dpMSum = out.dpCSum = out.dpTSum
		= out.dcASum = out.dcMSum = out.dcCSum = out.dcTSum = 0;
		}
	/* enforced initialization */
	inExt.bIndirect = 0;
	inExt.bDoSun = 0;
	inExt.bLogHalo = 0;
	inExt.bHernquistSpheroid = 0;
	inExt.bNFWSpheroid = 0;
	inExt.bElliptical= 0;
	inExt.bEllipticalDarkNFW = 0;
	inExt.bHomogSpheroid = 0;
	inExt.bBodyForce = 0;
	inExt.bMiyamotoDisk = 0;
	inExt.bTimeVarying = 0;
	inExt.bRotatingBar = 0;
#ifdef ROT_FRAME
	inExt.bRotFrame = 0;
#endif
#ifdef SLIDING_PATCH
/*	inExt.PP = NULL; --RB Debugging*/
#endif
#ifdef SIMPLE_GAS_DRAG
	inExt.bSimpleGasDrag = 0;
#endif
	/*
	 ** Calculate any external potential stuff.
	 ** This may contain a huge list of flags in the future, so we may want
	 ** to replace this test with something like bAnyExternal.
	 */
	if (msr->param.bHeliocentric || msr->param.bLogHalo ||
		msr->param.bHernquistSpheroid || msr->param.bNFWSpheroid ||
		msr->param.bElliptical ||
		msr->param.bHomogSpheroid || msr->param.bBodyForce ||
        	msr->param.bMiyamotoDisk || msr->param.bTimeVarying ||
	    	msr->param.bRotatingBar) {
	        struct outGravExternal outExt;
		/*
		 ** Provide the time.
		 */
		inExt.dTime = dStep*msr->param.dDelta;
		/*
		 ** Only allow inclusion of solar terms if we are in Heliocentric 
		 ** coordinates.
		 */
		inExt.bIndirect = msr->param.bHeliocentric;
		inExt.bDoSun = bDoSun;  /* Treat the Sun explicitly. */
		inExt.dSunMass = msr->param.dCentMass;
		inExt.dSunSoft = msr->param.dSunSoft;
		if (inExt.bIndirect)
		    for (j=0;j<3;++j) inExt.aSun[j] = out.aSun[j];
		inExt.bLogHalo = msr->param.bLogHalo;
		inExt.bHernquistSpheroid = msr->param.bHernquistSpheroid;
		if (inExt.bNFWSpheroid == msr->param.bNFWSpheroid) {
		  inExt.dNFWm200= msr->param.dNFWm200;
		  inExt.dNFWr200= msr->param.dNFWr200;
		  inExt.dNFWconc= msr->param.dNFWconc;
		  inExt.dNFWsoft= msr->param.dNFWsoft;
		}
		inExt.bElliptical = msr->param.bElliptical;
		inExt.bEllipticalDarkNFW = msr->param.bEllipticalDarkNFW;
		inExt.bHomogSpheroid = msr->param.bHomogSpheroid;
		inExt.bBodyForce = msr->param.bBodyForce;
		inExt.bMiyamotoDisk = msr->param.bMiyamotoDisk;
		inExt.bTimeVarying = msr->param.bTimeVarying;
		inExt.bRotatingBar = msr->param.bRotatingBar;
		if(msr->param.bRotatingBar) {
		    for (j=0;j<3;++j)
			inExt.aCom[j] = msr->param.rotbar->dPos[j];
		    
		    inExt.dRotBarAmp = msr->param.rotbar->amplitude;
		    inExt.dRotBarPosAng = msr->param.rotbar->dPosAng;
		    inExt.dRotBarB5 = msr->param.rotbar->dB5;
		    }
		pstGravExternal(msr->pst,&inExt,sizeof(inExt),&outExt, NULL);
		if(msr->param.bRotatingBar) {
		    struct outCalcEandL outL;
		    int iIgnore;
		    
		    pstCalcEandL(msr->pst, NULL, 0, &outL, &iIgnore);
		    
		    msr->param.rotbar->dLzPart = outL.L[2];
		    for (j=0;j<3;++j) {
			msr->param.rotbar->dAcc[j] = outExt.dAcc[j];
			msr->param.rotbar->dTorque[j] = outExt.dTorque[j];
			}
		    }
	    }

	/* if patch, do aggs next because the agg calculations are in a special frame */
#ifdef AGGS_IN_PATCH /*DEBUG but how do torques from Hill's eqns apply?*/
	msrAggsGravity(msr);
#endif

	/*DEBUG! extra calls to pstGravExternal() follow depending on
	  compile options, but there should only be *one* call per
	  step -- currently there is no check to prevent this!*/
	
#ifdef ROT_FRAME
	if (msr->param.bRotFrame) { /* general rotating frame */
	  inExt.bRotFrame = msr->param.bRotFrame;
	  inExt.dOmega = msr->param.dOmega +
	    dStep*msr->param.dDelta*msr->param.dOmegaDot;
	  inExt.dOmegaDot = msr->param.dOmegaDot;
	  pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
	}
#endif /* ROT_FRAME */

#ifdef SLIDING_PATCH
	if (msr->param.PP.bPatch) {
	  static int bFirstCall = 1;
	  if (bFirstCall) {
	    PATCH_PARAMS *PP = &msr->param.PP;
	    /* need to compute vertical frequency enhancement */
	    PP->dOrbFreqZ2 = PP->dOrbFreq*PP->dOrbFreq;
	    if (msr->param.bDoSelfGravity) {
	      if (PP->dWidth == FLOAT_MAXVAL || PP->dLength == FLOAT_MAXVAL)
		(void) printf("WARNING: Vert. freq. enhancement disabled\n"
			      "(no boundary condition in either x or y, or both!)\n");
	      else {
		double m_total = msrMassCheck(msr,-2.0,"");
		PP->dOrbFreqZ2 += 2*M_PI*m_total/msr->param.nReplicas/
		  pow(PP->dWidth*PP->dLength,1.5);
	      }
	    }
	    if (msr->param.bVStart)
	      (void) printf("Patch: vert. freq. enhancement factor = %g\n",
			    sqrt(PP->dOrbFreqZ2)/PP->dOrbFreq);
	    bFirstCall = 0;
	  }
	  inExt.PP = msr->param.PP; /* struct copy */
	  pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
	}
#endif /* SLIDING_PATCH */

#ifdef SIMPLE_GAS_DRAG
	if (msr->param.bSimpleGasDrag) {
	  inExt.bSimpleGasDrag = msr->param.bSimpleGasDrag;
	  inExt.iFlowOpt	= 1; /* temporary */
	  inExt.bEpstein	= msr->param.bEpstein;
	  inExt.dGamma	= msr->param.dGamma;
	  inExt.dTime		= dStep*msr->param.dDelta;
	  pstGravExternal(msr->pst,&inExt,sizeof(inExt),NULL,NULL);
	}
#endif /* SIMPLE_GAS_DRAG */

#ifdef COLLISIONS
	{
		struct inAddUnifGrav inaug;
		msrSetUnifGrav(msr,dStep*msr->param.dDelta); /* load values */
		inaug.dgx = msr->param.dxUnifGrav;
		inaug.dgy = msr->param.dyUnifGrav;
		inaug.dgz = msr->param.dzUnifGrav;
		pstAddUnifGrav(msr->pst,&inaug,sizeof(inaug),NULL,NULL);
	}
#endif

#ifdef CHARGE
	msrChargeZ(msr);
#endif

#ifdef DEM_TIDAL_LOCAL
	msrDEMTidal(msr,dStep*msr->param.dDelta);
#endif

	/* do springs/soft-sphere after all external potentials resolved */
	if (msr->param.FO.iForceOverrideOption == FO_STRENGTH) {
#ifdef COLLISIONS
#ifdef SPRINGS
		msrDoSprings(msr,dStep*msr->param.dDelta);
#endif
#ifdef DEM
		msrDoDEM(msr,dStep*msr->param.dDelta);
#ifdef DEM_WALLS_REACT
		msrDEMWallsReact(msr,msrDelta(msr));
#endif
		msrDEMZeroSmallMotions(msr);
#endif
#if !defined(SPRINGS) && !defined(DEM)
		assert(0);
#endif /* !SPRINGS && !DEM */
#endif /* COLLISIONS */	
		}

	/* finally, now that all particle accelerations are known, compute any aggregate center-of-mass accelerations and torques */

#ifdef SPINUP_WITH_AGGS
    msrSpinupAccel(msr, dStep);
#endif

#ifndef AGGS_IN_PATCH /* (msrAggsGravity(msr) called earlier in this case) */
#ifdef AGGS
	msrAggsGravity(msr);
#endif /* AGGS */
#endif /* !AGGS_IN_PATCH */

	/*
	 ** Output.
	 */
	*piSec = dsec;
	*pnActive = out.nActive;
	*pdWMax = out.dWMax;
	*pdIMax = out.dIMax;
	*pdEMax = out.dEMax;
	if (msr->param.bVStep) {
		double dPartAvg,dCellAvg,dSoftAvg;
		double dWAvg,dWMax,dWMin;
		double dIAvg,dIMax,dIMin;
		double dEAvg,dEMax,dEMin;
		double iP;

		/*
		 ** Output some info...
		 */
		if (dsec > 0.0) {
			double dMFlops = out.dFlop/dsec*1e-6;
			printf("Gravity Calculated, Wallclock: %f secs, MFlops:%.1f, Flop:%.3g\n",
				   dsec,dMFlops,out.dFlop);
			}
		else {
			printf("Gravity Calculated, Wallclock: %f secs, MFlops:unknown, Flop:%.3g\n",
				   dsec,out.dFlop);
			}
		if (out.nActive > 0) {
			dPartAvg = out.dPartSum/out.nActive;
			dCellAvg = out.dCellSum/out.nActive;
			dSoftAvg = out.dSoftSum/out.nActive;
			}
		else {
			dPartAvg = dCellAvg = dSoftAvg = 0;
#ifdef COLLISIONS
			/*
			 ** This is allowed to happen in the collision model because a
			 ** time-step rung may be left vacant following a merger event.
			 */
#else
			if (msr->param.bVWarnings)
				printf("WARNING: no particles found!\n");
#endif
			}
		iP = 1.0/msr->nThreads;
		dWAvg = out.dWSum*iP;
		dIAvg = out.dISum*iP;
		dEAvg = out.dESum*iP;
		dWMax = out.dWMax;
		dIMax = out.dIMax;
		dEMax = out.dEMax;
		dWMin = out.dWMin;
		dIMin = out.dIMin;
		dEMin = out.dEMin;
		printf("dPartAvg:%f dCellAvg:%f dSoftAvg:%f\n",
			   dPartAvg,dCellAvg,dSoftAvg);
		printf("Walk CPU     Avg:%10f Max:%10f Min:%10f\n",dWAvg,dWMax,dWMin);
		printf("Interact CPU Avg:%10f Max:%10f Min:%10f\n",dIAvg,dIMax,dIMin);
		if (msr->param.bEwald) printf("Ewald CPU    Avg:%10f Max:%10f Min:%10f\n",dEAvg,dEMax,dEMin);
		if (msr->nThreads > 1) {
			printf("Particle Cache Statistics (average per processor):\n");
			printf("    Accesses:    %10g\n",out.dpASum*iP);
			printf("    Miss Ratio:  %10g\n",out.dpMSum*iP);
			printf("    Min Ratio:   %10g\n",out.dpTSum*iP);
			printf("    Coll Ratio:  %10g\n",out.dpCSum*iP);
			printf("Cell Cache Statistics (average per processor):\n");
			printf("    Accesses:    %10g\n",out.dcASum*iP);
			printf("    Miss Ratio:  %10g\n",out.dcMSum*iP);
			printf("    Min Ratio:   %10g\n",out.dcTSum*iP);
			printf("    Coll Ratio:  %10g\n",out.dcCSum*iP);
			printf("\n");
			}
		}
	}


void msrCalcEandL(MSR msr,int bFirst,double dTime,double *E,double *T,
				  double *U,double *Eth,double L[])
{
	struct outCalcEandL out;
	struct inCalcEandLExt inExt;
	struct outCalcEandLExt outExt;
	double a;
	int k;

	pstCalcEandL(msr->pst,NULL,0,&out,NULL);
	*T = out.T;
	*U = out.U;
	*Eth = out.Eth;
	for (k=0;k<3;k++) L[k] = out.L[k];
	/*
	 ** Calculate E & L contribution from external potential and/or
	 ** reference frame. Currently only heliocentric frame implemented.
	 ** NOTE: In some cases the correction terms may be comparable to
	 ** the original values, resulting in precision errors.
	 */
	if (msr->param.bHeliocentric) {
		double dSunPos[3],dSunVel[3],dRxV[3],dTotMass,dSunVel2;
		inExt.bHeliocentric = msr->param.bHeliocentric;
		pstCalcEandLExt(msr->pst,&inExt,sizeof(inExt),&outExt,NULL);
		dTotMass = outExt.dMass + msr->param.dCentMass;
		assert(dTotMass > 0.0);
		dSunVel2 = 0;
		for (k=0;k<3;k++) {
			dSunPos[k] = - outExt.dSumMR[k]/dTotMass;
			dSunVel[k] = - outExt.dSumMV[k]/dTotMass;
			dSunVel2 += dSunVel[k]*dSunVel[k];
			}
		dRxV[0] = dSunPos[1]*dSunVel[2] - dSunPos[2]*dSunVel[1];
		dRxV[1] = dSunPos[2]*dSunVel[0] - dSunPos[0]*dSunVel[2];
		dRxV[2] = dSunPos[0]*dSunVel[1] - dSunPos[1]*dSunVel[0];
		*T -= 0.5*dTotMass*dSunVel2;
		*U -= 0.5*msr->param.dCentMass*outExt.dPot;
		for (k=0;k<3;k++) L[k] -= dTotMass*dRxV[k];
		}
	/*
	 ** Do the comoving coordinates stuff.
	 ** Currently L is not adjusted for this. Should it be?
	 */
	a = csmTime2Exp(msr->param.csm,dTime);
	if (!msr->param.bCannonical) *T *= pow(a,4.0);
	/*
	 * Estimate integral (\dot a*U*dt) over the interval.
	 * Note that this is equal to integral (W*da) and the latter
	 * is more accurate when a is changing rapidly.
	 */
	if (msrComove(msr) && !bFirst) {
		msr->dEcosmo += 0.5*(a - csmTime2Exp(msr->param.csm, msr->dTimeOld))
			*((*U) + msr->dUOld);
		}
	else {
		msr->dEcosmo = 0.0;
		}
	msr->dTimeOld = dTime;
	msr->dUOld = *U;
	*U *= a;
#ifdef COLLISIONS
	if (bFirst)
		msr->dTcoll = 0;
	else
		*T -= msr->dTcoll;
#endif
	*E = (*T) + (*U) - msr->dEcosmo + a*a*(*Eth);
	}


void msrDrift(MSR msr,double dTime,double dDelta)
{

	struct inDrift in;
	int j;

#ifdef COLLMOD_ZML
	if (msr->param.bCollDelay == 1) {
		if (dTime > msr->param.dCollTimer)
			msr->param.bCollDelay = 0; /*DEBUG not pretty to change input parameter*/
	}
#endif

#ifdef NEED_VPRED
	struct inKickVpred invpr;
	double a;
#endif

	/*
	 ** This only gets done if growmass parameters have been set!
	 */
	msrGrowMass(msr,dTime,dDelta);

#ifdef AGGS
	msrAggsAdvanceOpen(msr);
#endif

#ifdef COLLISIONS
#ifdef COLLMOD_ZML
	if (msr->param.bCollDelay == 0)
#endif
	{
		if (msr->param.iMinCollRung) {
			if (msr->iCurrMaxRung >= msr->param.iMinCollRung)
				msrDoCollisions(msr,dTime,dDelta);
			} else {
			msrDoCollisions(msr,dTime,dDelta);
			}
#ifdef RUBBLE_ZML
		if (msr->param.CP.bDoRubbleKDKRestart)
			return;
#elif defined(COLLMOD_ZML)
		if (msr->param.CP.bDoCollModKDKRestart)
			return;
#endif
		}
#endif

	if (msr->param.bCannonical) {
		in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
		}
	else {
		in.dDelta = dDelta;
		}
	for (j=0;j<3;++j) {
		in.fCenter[j] = msr->fCenter[j];
		}
	in.bPeriodic = msr->param.bPeriodic;
	in.bFandG = msr->param.bFandG;
	in.fCentMass = msr->param.dCentMass;
#ifdef SLIDING_PATCH
	in.dTime = dTime;
	in.PP = msr->param.PP; /* struct copy */
#endif
	pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);

	/*
	 ** Once we move the particles the tree should no longer be considered 
	 ** valid.
	 */
	msr->iTreeType = MSR_TREE_NONE;

	if(msr->param.bRotatingBar) {
	    rotbarDrift(msr->param.rotbar, dTime, dDelta);
	    }
#ifdef NEED_VPRED

#ifdef PREDRHO
	if (msr->param.bPredRho == 2) {
		struct inKickRhopred inRhop;
		if (msrComove(msr)) 
			inRhop.dHubbFac = 3*csmTime2Hub(msr->param.csm, dTime + dDelta/2.0);
		else
			inRhop.dHubbFac = 0.0;
		inRhop.dDelta = dDelta;
		/* Non Active Gas particles need to have updated densities */
		pstKickRhopred(msr->pst,&inRhop,sizeof(inRhop),NULL,NULL);
		}
#endif

	if (msr->param.bCannonical) {
#ifdef GLASS	  
		invpr.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		invpr.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
		invpr.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
		invpr.duDelta  = dDelta;
		invpr.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		invpr.z = 1/a - 1;
		invpr.duDotLimit = msr->param.duDotLimit;
		}
	else {
		double H;
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		invpr.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		invpr.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
		invpr.duDelta  = dDelta;
		invpr.iGasModel = msr->param.iGasModel;
		invpr.z = 1/a - 1;
		invpr.duDotLimit = msr->param.duDotLimit;
		}
	if (dDelta != 0.0) {
		struct outKick out;
		int nout;
		pstKickVpred(msr->pst,&invpr,sizeof(invpr),&out,&nout);
		if (nout) printf("Drift (Vpred): Avg Wallclock %f, Max Wallclock %f\n",
						 out.SumTime/out.nSum,out.MaxTime);
		LOGTIMINGUPDATE( out.MaxTime, TIMING_Drift );
	    }
#endif /* NEED_VPRED */

#ifdef SLIDING_PATCH
	if (msr->param.PP.bRandAzWrap == 1) {
		struct inRandAzWrap inrand;
		struct outRandAzWrap outrand;
		int n = 0;

		inrand.PP = msr->param.PP; /* struct copy */
		do {
			pstRandAzWrap(msr->pst,&inrand,sizeof(inrand),&outrand,NULL);
			/* (no point checking for wrap overlaps if collisions disabled) */
			if (inrand.PP.bNoRandomX || outrand.nRandomized == 0 || msr->param.nSmooth == 1)
				break;
			msrBuildTree(msr,0,-1.0,1); /* force rebuild of density tree since particles have moved */
			msrSmooth(msr,0.0,SMX_FIND_OVERLAPS,1); /* check for any overlaps, keep looping until none */
			} while (++n < inrand.PP.nWrapAttempts);
		if (n == inrand.PP.nWrapAttempts) {
			/*
			** We limit the number of attempts to reposition a
			** particle after azimuthal wrap in order to avoid a
			** possible infinite loop.  If all attempts fail, the
			** problem will be handled like a normal overlap during
			** the collision search (cf. CheckForCollision() in
			** smoothfcn.c).
			*/
			(void) fprintf(stderr,"WARNING: Unable to reposition particle after azimuthal boundary wrap.\n");
			}
		}
#endif

#ifdef AGGS
	msrAggsAdvanceClose(msr,dDelta);
#endif

#ifdef WALLS
	msrWallsMove(msr,dTime,dDelta);
#endif
	}

/* Untested */
void msrCoolVelocity(MSR msr,double dTime,double dMass)
{
#ifdef SUPERCOOL
	struct inCoolVelocity in;
	
	if (msr->param.nSuperCool > 0) {
		/*
		 ** Activate all for densities if bSymCool == 0
		 */
		if (msr->param.bSymCool) {
			msrActiveType(msr,TYPE_SUPERCOOL,TYPE_ACTIVE);
			msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveType(msr,TYPE_SUPERCOOL,TYPE_ACTIVE);
			/* Unsure what is desired here -- assuming all particles are in tree 
			   as per above setting of TREEACTIVE-- JW */
			msrBuildTree(msr,0,dMass,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrReSmooth(msr,dTime,SMX_MEANVEL,1);
			}
		else {
			/*
			 ** Note, here we need to calculate the densities of all
			 ** the particles so that the mean velocity can be 
			 ** calculated.
			 */
			/* activate all */
			msrActiveTypeRung(msr,TYPE_SUPERCOOL,TYPE_ACTIVE,0,1); 
			msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrDomainDecomp(msr,0,1);
			msrActiveTypeRung(msr,TYPE_SUPERCOOL,TYPE_ACTIVE,0,1); 
			msrBuildTree(msr,0,SMX_DENSITY,1);
			msrSmooth(msr,dTime,SMX_DENSITY,0);
			msrReSmooth(msr,dTime,SMX_MEANVEL,0);
			}
		/*
		 ** Now cool them.
		 */
		in.nSuperCool = msr->param.nSuperCool;
		in.dCoolFac = msr->param.dCoolFac;
		in.dCoolDens = msr->param.dCoolDens;
		in.dCoolMaxDens = msr->param.dCoolMaxDens;
		pstCoolVelocity(msr->pst,&in,sizeof(in),NULL,NULL);
		}
#endif
	}

void msrGrowMass(MSR msr, double dTime, double dDelta)
{
    struct inGrowMass in;
    
    if (msr->param.nGrowMass > 0 && dTime > msr->param.dGrowStartT &&
		dTime <= msr->param.dGrowEndT) {
		in.nGrowMass = msr->param.nGrowMass;
		in.dDeltaM = msr->param.dGrowDeltaM*dDelta/
			(msr->param.dGrowEndT - msr->param.dGrowStartT);
		pstGrowMass(msr->pst, &in, sizeof(in), NULL, NULL);
		}
    }

/*
 * For gasoline, updates predicted velocities to middle of timestep.
 */
void msrKickDKD(MSR msr,double dTime,double dDelta)
{
	double H,a;
	struct inKick in;
	struct outKick out;
	int j;

#ifndef NEED_VPRED
	in.dvPredFacOne = in.dvPredFacTwo = in.duDelta = in.duPredDelta = in.duDotLimit =
		in.iGasModel = in.z = 0; /* to prevent FPE from uninitialized values */
#endif

	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = csmComoveKickFac(msr->param.csm,dTime,0.5*dDelta);
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.5*dDelta;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		dTime -= dDelta/4.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvPredFacOne = (1.0 - 0.5*H*dDelta)/(1.0 + 0.5*H*dDelta);
		in.dvPredFacTwo = 0.5*dDelta/pow(a,3.0)/(1.0 + 0.5*H*dDelta);
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.5*dDelta;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	printf("Kick: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	if(msr->param.bRotatingBar && msr->param.dDelta == dDelta) {
	    for (j=0;j<3;++j) {
		msr->param.rotbar->dVel[j]
		    = msr->param.rotbar->dVel[j]*in.dvFacOne
		    + msr->param.rotbar->dAcc[j]*in.dvFacTwo;
		}
	    }
	LOGTIMINGUPDATE( out.MaxTime, TIMING_Kick );
	}

/*
 * For gasoline, updates predicted velocities to beginning of timestep.
 */
void msrKickKDKOpen(MSR msr,double dTime,double dDelta)
{
	/* NOTE: dDelta should be 1/2 the drift step here... */

	double H,a;
	struct inKick in;
	struct outKick out;

#ifndef NEED_VPRED
	in.dvPredFacOne = in.dvPredFacTwo = in.duDelta = in.duPredDelta = in.duDotLimit =
		in.iGasModel = in.z = 0; /* to prevent FPE from uninitialized values */
#endif

	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper;  /* Damp velocities */
#else
		in.dvFacOne = 1.0;		/* no hubble drag, man! */
#endif
		in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper;  /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = 0.0;
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.0;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		in.dvPredFacOne = 1.0;
		in.dvPredFacTwo = 0.0;
		in.duDelta      = dDelta;
		in.duPredDelta  = 0.0;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
#ifndef SLIDING_PATCH
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
#else
	if(!msr->param.PP.bPatch) {
	    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	    }
	else {
	    struct inKickPatch inp;

	    inp.bOpen = 1;
	    inp.dOrbFreq = msr->param.PP.dOrbFreq;
	    inp.dvFacOne = 1.0;
	    inp.dvFacTwo = in.dvFacTwo;
#ifdef NEED_VPRED
		inp.dvPredFacOne = 1.0;
		inp.dvPredFacTwo = 0.0;
#endif /* NEED_VPRED */
	    pstKickPatch(msr->pst,&inp,sizeof(inp),&out,NULL);
	    }
#endif /* SLIDING_PATCH */
	    
	if (msr->param.bVDetails) 
		printf("KickOpen: Avg Wallclock %f, Max Wallclock %f\n",
			   out.SumTime/out.nSum,out.MaxTime);
	if(msr->param.bRotatingBar && msr->param.dDelta*.5 == dDelta) {
	    rotbarKick(msr->param.rotbar, in.dvFacOne, in.dvFacTwo);
	    }
	
	LOGTIMINGUPDATE( out.MaxTime, TIMING_Kick );

#ifdef DEM_WALLS_REACT
	msrWallsKick(msr, dDelta);
#endif /* DEM_WALLS_REACT */

#ifdef AGGS
#ifdef AGGS_IN_PATCH
	msrAggsInPatchKick(msr,dDelta,1/*open*/);
#else
	msrAggsKick(msr,dDelta);
#endif
#endif /* AGGS */
	}

/*
 * For gasoline, updates predicted velocities to end of timestep.
 */
void msrKickKDKClose(MSR msr,double dTime,double dDelta)
{
	/* NOTE: dDelta should be 1/2 the drift step here... */

	double H,a;
	struct inKick in;
	struct outKick out;

#ifndef NEED_VPRED
	in.dvPredFacOne = in.dvPredFacTwo = in.duDelta = in.duPredDelta = in.duDotLimit =
		in.iGasModel = in.z = 0; /* to prevent FPE from uninitialized values */
#endif
	
	if (msr->param.bCannonical) {
#ifdef GLASS	  
		in.dvFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvFacOne = 1.0; /* no hubble drag, man! */
#endif
		in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
#ifdef NEED_VPRED
#ifdef GLASS	  
		in.dvPredFacOne = 1.0 - fabs(dDelta)*msr->param.dGlassDamper; /* Damp velocities */
#else
		in.dvPredFacOne = 1.0;
#endif
		in.dvPredFacTwo = in.dvFacTwo;
		in.duDelta      = dDelta;
		in.duPredDelta  = dDelta;
		in.iGasModel = msr->param.iGasModel;
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
	else {
		/*
		 ** Careful! For non-cannonical we want H and a at the 
		 ** HALF-STEP! This is a bit messy but has to be special
		 ** cased in some way.
		 */
		dTime += dDelta/2.0;
		a = csmTime2Exp(msr->param.csm,dTime);
		H = csmTime2Hub(msr->param.csm,dTime);
		in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
		in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
#ifdef NEED_VPRED
		in.dvPredFacOne = in.dvFacOne;
		in.dvPredFacTwo = in.dvFacTwo;
		in.duDelta      = dDelta;
		in.duPredDelta  = dDelta;
		in.iGasModel = msr->param.iGasModel;
		in.z = 1/a - 1;
		in.duDotLimit = msr->param.duDotLimit;
#endif /* NEED_VPRED */
		}
#ifndef SLIDING_PATCH
	pstKick(msr->pst,&in,sizeof(in),&out,NULL);
#else
	if(!msr->param.PP.bPatch) {
	    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
	    }
	else {
	    struct inKickPatch inp;

	    inp.bOpen = 0;
	    inp.dOrbFreq = msr->param.PP.dOrbFreq;
	    inp.dvFacOne = 1.0;
	    inp.dvFacTwo = in.dvFacTwo;
#ifdef NEED_VPRED
		inp.dvPredFacOne = 1.0;
		inp.dvPredFacTwo = in.dvFacTwo;
#endif /* NEED_VPRED */
	    pstKickPatch(msr->pst,&inp,sizeof(inp),&out,NULL);
	    }
#endif /* SLIDING_PATCH*/
	if (msr->param.bVDetails)
		printf("KickClose: Avg Wallclock %f, Max Wallclock %f\n",
			   out.SumTime/out.nSum,out.MaxTime);
	if(msr->param.bRotatingBar && msr->param.dDelta*.5 == dDelta) {
	    rotbarKick(msr->param.rotbar, in.dvFacOne, in.dvFacTwo);
	    }
	LOGTIMINGUPDATE( out.MaxTime, TIMING_Kick );

#ifdef DEM_WALLS_REACT
	msrWallsKick(msr, dDelta);
#endif /* DEM_WALLS_REACT */

#ifdef AGGS
#ifdef AGGS_IN_PATCH
	msrAggsInPatchKick(msr,dDelta,0/*close*/);
#else
	msrAggsKick(msr,dDelta);
#endif
#endif /* AGGS */
	}

void msrOneNodeReadCheck(MSR msr, struct inReadCheck *in)
{
    int i,id;
    int *nParts;				/* number of particles for each processor */
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
#ifdef COLLISIONS
	struct inReadSS initin;
#else
    struct inReadTipsy initin;
#endif
    int j;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    initin.nFileStart = in->nFileStart;
    initin.nFileEnd = in->nFileEnd;
    initin.iOrder = in->iOrder;
    initin.fExtraStore = in->fExtraStore;
    initin.nDark = in->nDark;
    initin.nGas = in->nGas;
    initin.nStar = in->nStar;
#ifndef COLLISIONS /* extra stuff for tipsy format */
    initin.bStandard = 0;
	initin.iReadIOrder = 0;
    initin.dvFac = 1;
    initin.dTuFac = 1;
#endif

    for(j = 0; j < 3; j++)
		initin.fPeriod[j] = in->fPeriod[j];

    pstOneNodeReadInit(msr->pst, &initin, sizeof(initin), nParts, &nid);
    assert(nid == (int) (msr->nThreads*sizeof(*nParts)));
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
    assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadCheck(plcl->pkd,achInFile,in->iVersion,
					 in->iOffset, nStart, nParts[id]);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadCheck(plcl->pkd,achInFile,in->iVersion,in->iOffset,0,nParts[0]);
    }

int msrFindCheck(MSR msr) {
	char achCheckFile[PST_FILENAME_SIZE];
	char achLatestCheckFile[PST_FILENAME_SIZE];
	time_t latestTime;
	char extraString[PST_FILENAME_SIZE];
	int iNotCorrupt;
	FDL_CTX* fdl;
	struct stat statbuf;
	const char* suffixes[3] = {"chk0", "chk1", "chk"};
	int i;
	
	latestTime = 0;
	
	for(i = 0; i < 3; ++i) {
		/* Get the filename of the checkpoint file formatted properly. */
		snprintf(achCheckFile, PST_FILENAME_SIZE, "%s/%s.%s", msr->param.achDataSubPath, msr->param.achOutName, suffixes[i]);
		if(msr->pst->plcl->pszDataPath) {
			strcpy(extraString, achCheckFile);
			snprintf(achCheckFile, PST_FILENAME_SIZE, "%s/%s", msr->pst->plcl->pszDataPath, extraString);
		}

		/* Check if the checkpoint file exists. */
		if(!stat(achCheckFile, &statbuf)) {
			/* Check if the checkpoint file is corrupt. */
			fdl = FDL_open(achCheckFile);
			/* Read and ignore the version number. */
			FDL_read(fdl, "version" ,&iNotCorrupt);
			FDL_read(fdl, "not_corrupt_flag", &iNotCorrupt);
			if(iNotCorrupt == 1) {
				/* Get the file creation time. */
				if(statbuf.st_mtime > latestTime) {
					/* The checkpoint file exists, isn't corrupt, and is more recent, so mark it. */
					latestTime = statbuf.st_mtime;
					strcpy(achLatestCheckFile, achCheckFile);
				}
			}
			FDL_finish(fdl);
		}
	}
	
	/* Did we find any checkpoint files? */
	if(latestTime == 0)
		return 0;
	
	/* On the Cray XT3, renaming a file to the same name will not return success,
	 * so check to make sure LatestCheckFile and CheckFile are not the same name. */
	if (strcmp(achLatestCheckFile, achCheckFile) == 0)
	     return 1;	     /* They are the same name, so no need to rename */

	/* Rename latest valid checkpoint file to .chk, return success. */
	if(!rename(achLatestCheckFile, achCheckFile)) /* The last filename checked is the one we want. */
		return 1;
	else
		return 0;
}

double msrReadCheck(MSR msr,int *piStep)
{
	struct inReadCheck in;
	struct inSetNParts inset;
	struct inSetParticleTypes intype;
	char achInFile[PST_FILENAME_SIZE];
	int i;
	LCL *plcl = msr->pst->plcl;
	double dTime;
	int iVersion,iNotCorrupt;
	FDL_CTX *fdl;
	int nMaxOutsTmp,nOutsTmp;
	double *pdOutTimeTmp;
	double sec=0.0,dsec=0.0;
#ifdef ORIGIN_HISTOGRAM
	int j;
#endif /* ORIGIN_HISTOGRAM */

	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	snprintf(achInFile,PST_FILENAME_SIZE,"%s/%s.chk",
			 msr->param.achDataSubPath,msr->param.achOutName);
	strcpy(in.achInFile,achInFile);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		snprintf(achInFile,PST_FILENAME_SIZE,"%s/%s",plcl->pszDataPath,in.achInFile);
		}
	fdl = FDL_open(achInFile);
	FDL_read(fdl,"version",&iVersion);
	if (msr->param.bVDetails)
		printf("Reading Version-%d Checkpoint file.\n",iVersion);
	FDL_read(fdl,"not_corrupt_flag",&iNotCorrupt);
	if (iNotCorrupt != 1) {
		printf("Sorry the checkpoint file is corrupted.\n");
		_msrExit(msr,1);
		}
	FDL_read(fdl,"number_of_particles",&msr->N);
	/*
	 ** As of checkpoint version 3 we include numbers of dark gas and star 
	 ** particles to support GASOLINE.
	 */
	if (iVersion > 2) {
		FDL_read(fdl,"number_of_gas_particles",&msr->nGas);
		FDL_read(fdl,"number_of_dark_particles",&msr->nDark);
		FDL_read(fdl,"number_of_star_particles",&msr->nStar);
		assert(msr->N == msr->nDark+msr->nGas+msr->nStar);
		}
	else {
		msr->nDark = msr->N;
		msr->nGas = 0;
		msr->nStar = 0;
		}
	FDL_read(fdl,"current_timestep",piStep);
	FDL_read(fdl,"current_time",&dTime);
	FDL_read(fdl,"current_ecosmo",&msr->dEcosmo);
	FDL_read(fdl,"old_time",&msr->dTimeOld);
	FDL_read(fdl,"old_potentiale",&msr->dUOld);
#ifdef COLLISIONS
	msr->dTcoll = 0; /* collisional energy loss not currently stored in checkpoint file... */
#endif
	if (!msr->bOpenSpec) {
		FDL_read(fdl,"opening_type",&msr->iOpenType);
		FDL_read(fdl,"opening_criterion",&msr->dCrit);
		}
	FDL_read(fdl,"number_of_outs",&msr->nOuts);
	if (msr->nOuts > msr->nMaxOuts) {
		msr->nMaxOuts = msr->nOuts;
		msr->pdOutTime = realloc(msr->pdOutTime,msr->nMaxOuts*sizeof(double));
		assert(msr->pdOutTime != NULL);
		}
	for (i=0;i<msr->nOuts;++i) {
		FDL_index(fdl,"out_time_index",i);
		FDL_read(fdl,"out_time",&msr->pdOutTime[i]);
		}
#ifdef RUBBLE_ZML
	if (iVersion > 80) {
		printf("checkpoint version = %d\n", iVersion);
		FDL_read(fdl,"num_rub_events", &msr->re.nEvents);
		for (i=0;i<msr->re.nEvents;i++) {
			FDL_index(fdl,"rub_event_index",i);
			FDL_read(fdl,"rub_iColor",&msr->re.rc[i].iColor);
			FDL_read(fdl,"rub_dTStartMergePhase",
					 &msr->re.rc[i].dTStartMergePhase);
			FDL_read(fdl,"rub_dTEndRubblePhase",
					 &msr->re.rc[i].dTEndRubblePhase);
			}
		FDL_read(fdl,"num_dust_bins", &msr->param.CP.DB.nDustBins);
		/*Do I need to add iDustBinsVelDispOpt 05.24.07*/
#ifdef ORIGIN_HISTOGRAM
		assert(msr->param.CP.DB.nDustBins == NUM_ORIGIN_BINS);
#endif /* ORIGIN_HISTOGRAM */
		for (i=0;i<msr->param.CP.DB.nDustBins;i++) {
			FDL_index(fdl,"dust_bins_index",i);
			FDL_read(fdl,"dust_bins_mass",&msr->aDustBins[i].dMass);
			FDL_read(fdl,"dust_bins_volume",&msr->aDustBins[i].dVolume);
#ifdef ORIGIN_HISTOGRAM
			if (iVersion > 81) {
				for (j=0;j<NUM_ORIGIN_BINS;j++) {
					FDL_index(fdl,"dust_bins_origin_bins_index",j);
					FDL_read(fdl,"dust_bins_origin_bins_value",&msr->aDustBins[i].origin_bins[j]);
					}
				}
#endif /* ORIGIN_HISTOGRAM */
			}
		FDL_read(fdl,"dust_bins_trash",&msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
		if (iVersion > 81) {
			for (i=0;i<NUM_ORIGIN_BINS;i++) {
				FDL_index(fdl,"trash_bin_origin_bins_index",i);
				FDL_read(fdl,"trash_bin_origin_bins_value",&msr->DustBinsTrash.origin_bins[i]);
				}
			}
#endif /* ORIGIN_HISTOGRAM */
	}
#elif defined(COLLMOD_ZML) /* TIDY */
	if (iVersion > 90) {
		printf("checkpoint version = %d\n", iVersion);
		FDL_read(fdl,"num_dust_bins", &msr->param.CP.DB.nDustBins);
		printf("num dust bins %d\n",msr->param.CP.DB.nDustBins);
		/*Do I need to add iDustBinsVelDispOpt 05.24.07*/
#ifdef ORIGIN_HISTOGRAM
		assert(msr->param.CP.DB.nDustBins == NUM_ORIGIN_BINS);
#endif /* ORIGIN_HISTOGRAM */
		for (i=0;i<msr->param.CP.DB.nDustBins;i++) {
			FDL_index(fdl,"dust_bins_index",i);
			FDL_read(fdl,"dust_bins_mass",&msr->aDustBins[i].dMass);
			FDL_read(fdl,"dust_bins_volume",&msr->aDustBins[i].dVolume);
#ifdef ORIGIN_HISTOGRAM
			if (iVersion > 91) {
				for (j=0;j<NUM_ORIGIN_BINS;j++) {
					FDL_index(fdl,"dust_bins_origin_bins_index",j);
					FDL_read(fdl,"dust_bins_origin_bins_value",&msr->aDustBins[i].origin_bins[j]);
					}
				}
#endif /* ORIGIN_HISTOGRAM */
			}
		FDL_read(fdl,"dust_bins_trash",&msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
		if (iVersion > 91) {
			for (i=0;i<NUM_ORIGIN_BINS;i++) {
				FDL_index(fdl,"trash_bin_origin_bins_index",i);
				FDL_read(fdl,"trash_bin_origin_bins_value",&msr->DustBinsTrash.origin_bins[i]);
				}
			}
#endif /* ORIGIN_HISTOGRAM */
	}
#endif /* RUBBLE_ZML, COLLMOD_ZML */

	/*
	 ** Read the old parameters.
	 */
	FDL_read(fdl,"bPeriodic",&msr->param.bPeriodic);
	FDL_read(fdl,"bComove",&msr->param.csm->bComove);
	if (!prmSpecified(msr->prm,"bParaRead"))
		FDL_read(fdl,"bParaRead",&msr->param.bParaRead);
	if (!prmSpecified(msr->prm,"bParaWrite"))
		FDL_read(fdl,"bParaWrite",&msr->param.bParaWrite);
	/*
	 ** Checkpoints can NOT be switched to a different coordinate system!
	 */
	FDL_read(fdl,"bCannonical",&msr->param.bCannonical);
	if (!prmSpecified(msr->prm,"bStandard"))
		FDL_read(fdl,"bStandard",&msr->param.bStandard);
	FDL_read(fdl,"bKDK",&msr->param.bKDK);
	/*
	 ** bBinary somehow got left out of the previous versions of the
	 ** checkpoint file. We fix this as of checkpoint version 5.
	 */
	if (iVersion > 4) {
		if (!prmSpecified(msr->prm,"bBinary")) 
			FDL_read(fdl,"bBinary",&msr->param.bBinary);
		}
	if (!prmSpecified(msr->prm,"nBucket"))
		FDL_read(fdl,"nBucket",&msr->param.nBucket);
	if (!prmSpecified(msr->prm,"iOutInterval"))
		FDL_read(fdl,"iOutInterval",&msr->param.iOutInterval);
	if (!prmSpecified(msr->prm,"iLogInterval"))
		FDL_read(fdl,"iLogInterval",&msr->param.iLogInterval);
	if (!prmSpecified(msr->prm,"iCheckInterval"))
		FDL_read(fdl,"iCheckInterval",&msr->param.iCheckInterval);
	if (!prmSpecified(msr->prm,"iOrder"))
		FDL_read(fdl,"iExpOrder",&msr->param.iOrder);
	if (!prmSpecified(msr->prm,"iEwOrder"))
		FDL_read(fdl,"iEwOrder",&msr->param.iEwOrder);
	if (!prmSpecified(msr->prm,"nReplicas"))
		FDL_read(fdl,"nReplicas",&msr->param.nReplicas);
	if (!prmSpecified(msr->prm,"nSteps"))
		FDL_read(fdl,"nSteps",&msr->param.nSteps);
	if (!prmSpecified(msr->prm,"dExtraStore"))
		FDL_read(fdl,"dExtraStore",&msr->param.dExtraStore);
	if (!prmSpecified(msr->prm,"dDelta"))
		FDL_read(fdl,"dDelta",&msr->param.dDelta);
	if (!prmSpecified(msr->prm,"dEta"))
		FDL_read(fdl,"dEta",&msr->param.dEta);
	if (iVersion > 4) {
	    if (!prmSpecified(msr->prm,"dEtaCourant"))
		FDL_read(fdl,"dEtaCourant",&msr->param.dEtaCourant);
	    }
	if (!prmSpecified(msr->prm,"bEpsAccStep"))
		FDL_read(fdl,"bEpsAccStep",&msr->param.bEpsAccStep);
	if (!prmSpecified(msr->prm,"bNonSymp"))
		FDL_read(fdl,"bNonSymp",&msr->param.bNonSymp);
	if (!prmSpecified(msr->prm,"iMaxRung"))
		FDL_read(fdl,"iMaxRung",&msr->param.iMaxRung);
	if (!prmSpecified(msr->prm,"dEwCut"))
		FDL_read(fdl,"dEwCut",&msr->param.dEwCut);
	if (!prmSpecified(msr->prm,"dEwhCut"))
		FDL_read(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_read(fdl,"dPeriod",&msr->param.dPeriod);
	if (iVersion > 3) {
		FDL_read(fdl,"dxPeriod",&msr->param.dxPeriod);
		FDL_read(fdl,"dyPeriod",&msr->param.dyPeriod);
		FDL_read(fdl,"dzPeriod",&msr->param.dzPeriod);
		}
	else {
		msr->param.dxPeriod = msr->param.dPeriod;
		msr->param.dyPeriod = msr->param.dPeriod;
		msr->param.dzPeriod = msr->param.dPeriod;
		}
	FDL_read(fdl,"dHubble0",&msr->param.csm->dHubble0);
	FDL_read(fdl,"dOmega0",&msr->param.csm->dOmega0);
	if (iVersion > 4) {
	    FDL_read(fdl,"dLambda",&msr->param.csm->dLambda);
	    FDL_read(fdl,"dOmegaRad",&msr->param.csm->dOmegaRad);
	    FDL_read(fdl,"dQuintess",&msr->param.csm->dQuintess);
	    }
	if (iVersion > 3) {
		if (!prmSpecified(msr->prm,"dTheta"))
			FDL_read(fdl,"dTheta",&msr->param.dTheta);
		if (!prmSpecified(msr->prm,"dTheta2"))
			FDL_read(fdl,"dTheta2",&msr->param.dTheta2);
		}
	else {
		if (!prmSpecified(msr->prm,"dTheta"))
			msr->param.dTheta = msr->dCrit;
		if (!prmSpecified(msr->prm,"dTheta2"))
			msr->param.dTheta2 = msr->dCrit;
		}
	if (iVersion > 5) {
		if (!prmSpecified(msr->prm,"bDoGravity"))
			FDL_read(fdl,"bDoGravity",&msr->param.bDoGravity);
		if (!prmSpecified(msr->prm,"bDoGas"))
			FDL_read(fdl,"bDoGas",&msr->param.bDoGas);
		if (!prmSpecified(msr->prm,"dEtaCourant"))
			FDL_read(fdl,"dEtaCourant",&msr->param.dEtaCourant);
		if (!prmSpecified(msr->prm,"iStartStep"))
			FDL_read(fdl,"iStartStep",&msr->param.iStartStep);
		if (!prmSpecified(msr->prm,"dFracNoDomainDecomp"))
			FDL_read(fdl,"dFracNoDomainDecomp",&msr->param.dFracNoDomainDecomp);
#ifndef NOCOOLING
#if defined(COOLING_COSMO) && defined(GASOLINE)
		if (!prmSpecified(msr->prm,"dMassFracHelium"))
			FDL_read(fdl,"dMassFracHelium",&msr->param.CoolParam.dMassFracHelium);
#endif
#endif
		}
	if (iVersion > 3) {
	    FDL_read(fdl, "max_order", &msr->nMaxOrder);
	    FDL_read(fdl, "max_order_gas", &msr->nMaxOrderGas);
	    FDL_read(fdl, "max_order_dark", &msr->nMaxOrderDark);
	    }
	else {
	    msr->nMaxOrder = msr->N - 1;
	    msr->nMaxOrderGas = -1;
	    msr->nMaxOrderDark = msr->nMaxOrder;
	    }
	if (iVersion > 4) {
		FDL_read(fdl,"bFandG",&msr->param.bFandG);
		FDL_read(fdl,"bHeliocentric",&msr->param.bHeliocentric);
		FDL_read(fdl,"dCentMass",&msr->param.dCentMass);
		FDL_read(fdl,"bLogHalo",&msr->param.bLogHalo);
		FDL_read(fdl,"bHernquistSpheroid",&msr->param.bHernquistSpheroid);
		FDL_read(fdl,"bMiyamotoDisk",&msr->param.bMiyamotoDisk);
		}
	else {
		msr->param.bFandG = 0;
		msr->param.bHeliocentric = 0;
		msr->param.dCentMass = 0.0;
		msr->param.bLogHalo = 0;
		msr->param.bHernquistSpheroid = 0;
		msr->param.bMiyamotoDisk = 0;
		}

	if (iVersion > 6) {
	    FDL_read(fdl,"bRotatingBar", &msr->param.bRotatingBar);
	    rotbarCheckRead(msr->param.rotbar, fdl);
	    }
	else {
	    msr->param.bRotatingBar = 0;
	    }
	
	/*
	 * Check if redshift file is present, and if so reread it --JPG
	 */
	/* Store old values in temporary space */
	pdOutTimeTmp = malloc(msr->nMaxOuts*sizeof(double));
	nMaxOutsTmp = msr->nMaxOuts;
	nOutsTmp = msr->nOuts;
	for (i=0;i<msr->nOuts;++i) {
	    pdOutTimeTmp[i] = msr->pdOutTime[i];
	}
	/* Calculate initial time to give to msrReadOuts*/
	msrReadOuts(msr,dTime - (msr->param.nSteps*msr->param.dDelta));
	if (msr->nOuts == 0) { /* No redshift file...use old values */
	    free(msr->pdOutTime);
	    msr->pdOutTime = pdOutTimeTmp;
	    msr->nMaxOuts = nMaxOutsTmp;
	    msr->nOuts = nOutsTmp;
	}

	if (msr->param.bVDetails) {
		printf("Reading checkpoint file...\nN:%d Time:%g\n",msr->N,dTime);
		sec = msrTime(msr);
	        }
	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;
	in.nStar = msr->nStar;
	in.iOrder = msr->param.iOrder;
	/*
	 ** Since pstReadCheck causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;
	/*
	 ** Provide the period.
	 */
	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;
	in.iVersion = iVersion;
	in.iOffset = FDL_offset(fdl,"particle_array");
	FDL_finish(fdl);
	if(msr->param.bParaRead)
	    pstReadCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadCheck(msr,&in);
	if (msr->param.bVDetails) {
		puts("Checkpoint file has been successfully read.");
		dsec = msrTime(msr) - sec;
		printf("Data read complete, Wallclock: %f secs\n",dsec);
	        }
	inset.nGas = msr->nGas;
	inset.nDark = msr->nDark;
	inset.nStar = msr->nStar;
	inset.nMaxOrderGas = msr->nMaxOrderGas;
	inset.nMaxOrderDark = msr->nMaxOrderDark;
	inset.nMaxOrder = msr->nMaxOrder;
	pstSetNParts(msr->pst,&inset,sizeof(inset),NULL,NULL);
	intype.nSuperCool = msr->param.nSuperCool;
	pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
	
	i = msrCountType(msr, TYPE_GAS, TYPE_GAS);
	assert(i == msr->nGas);
	i = msrCountType(msr, TYPE_DARK, TYPE_DARK);
	assert(i == msr->nDark);
	i = msrCountType(msr, TYPE_STAR, TYPE_STAR);
	assert(i == msr->nStar);
	
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}

void msrOneNodeWriteCheck(MSR msr, struct inWriteCheck *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset, 0); 
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset, nStart); 
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }

void msrWriteCheck(MSR msr,double dTime,int iStep)
{
	struct inWriteCheck in;
	char achOutFile[PST_FILENAME_SIZE],achTmp[PST_FILENAME_SIZE];
	int i;
	LCL *plcl = msr->pst->plcl;
	FDL_CTX *fdl;
	char *pszFdl;
	int iVersion,iNotCorrupt;
	static int first = 1;
	double sec=0.0,dsec=0.0;
#ifdef ORIGIN_HISTOGRAM
	int j;
#endif /* ORIGIN_HISTOGRAM */
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	if(first) {
	    snprintf(achTmp,PST_FILENAME_SIZE,"%s.chk0",msr->param.achOutName);
	    first = 0;
	} else {
		snprintf(achTmp,PST_FILENAME_SIZE,"%s.chk1",msr->param.achOutName);
		first = 1;
		}
	_msrMakePath(msr->param.achDataSubPath,achTmp,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	pszFdl = getenv("PKDGRAV_CHECKPOINT_FDL");
	if (pszFdl == NULL) {
          fprintf(stderr,"PKDGRAV_CHECKPOINT_FDL environment variable not set: no Checkpoint written\n");
          return;
        }	  
	fdl = FDL_create(achOutFile,pszFdl);
	iVersion = CHECKPOINT_VERSION;
	FDL_write(fdl,"version",&iVersion);
	iNotCorrupt = 0;
	FDL_write(fdl,"not_corrupt_flag",&iNotCorrupt);
	/*
	 ** Checkpoint header.
	 */
	FDL_write(fdl,"number_of_particles",&msr->N);
	FDL_write(fdl,"number_of_gas_particles",&msr->nGas);
	FDL_write(fdl,"number_of_dark_particles",&msr->nDark);
	FDL_write(fdl,"number_of_star_particles",&msr->nStar);
	FDL_write(fdl, "max_order", &msr->nMaxOrder);
	FDL_write(fdl, "max_order_gas", &msr->nMaxOrderGas);
	FDL_write(fdl, "max_order_dark", &msr->nMaxOrderDark);
	FDL_write(fdl,"current_timestep",&iStep);
	FDL_write(fdl,"current_time",&dTime);
	FDL_write(fdl,"current_ecosmo",&msr->dEcosmo);
	FDL_write(fdl,"old_time",&msr->dTimeOld);
	FDL_write(fdl,"old_potentiale",&msr->dUOld);
	FDL_write(fdl,"opening_type",&msr->iOpenType);
	FDL_write(fdl,"opening_criterion",&msr->dCrit);
	FDL_write(fdl,"number_of_outs",&msr->nOuts);
	for (i=0;i<msr->nOuts;++i) {
		FDL_index(fdl,"out_time_index",i);
		FDL_write(fdl,"out_time",&msr->pdOutTime[i]);
		}
#ifdef RUBBLE_ZML
	FDL_write(fdl,"num_rub_events", &msr->re.nEvents);
	for (i=0;i<msr->re.nEvents;i++) {
		FDL_index(fdl,"rub_event_index",i);
		FDL_write(fdl,"rub_iColor",&msr->re.rc[i].iColor);
		FDL_write(fdl,"rub_dTStartMergePhase",
				  &msr->re.rc[i].dTStartMergePhase);
		FDL_write(fdl,"rub_dTEndRubblePhase",
				  &msr->re.rc[i].dTEndRubblePhase);
		}
	FDL_write(fdl,"num_dust_bins",&msr->param.CP.DB.nDustBins);
	for (i=0;i<msr->param.CP.DB.nDustBins;i++) {
		FDL_index(fdl,"dust_bins_index",i);
		FDL_write(fdl,"dust_bins_mass",&msr->aDustBins[i].dMass);
		FDL_write(fdl,"dust_bins_volume",&msr->aDustBins[i].dVolume);
#ifdef ORIGIN_HISTOGRAM
		/* for some reason it's not possible to use num_dust_bins
		   in the FDL file for SIZE of dust_bins_origin_bins_array --
		   it must be hardcoded! */
		for (j=0;j<NUM_ORIGIN_BINS;j++) {
			FDL_index(fdl,"dust_bins_origin_bins_index",j);
			FDL_write(fdl,"dust_bins_origin_bins_value",&msr->aDustBins[i].origin_bins[j]);
			}
#endif /* ORIGIN_HISTOGRAM */
		}	
	FDL_write(fdl,"dust_bins_trash",&msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
	for (i=0;i<NUM_ORIGIN_BINS;i++) {
		FDL_index(fdl,"trash_bin_origin_bins_index",i);
		FDL_write(fdl,"trash_bin_origin_bins_value",&msr->DustBinsTrash.origin_bins[i]);
		}
#endif /* ORIGIN_HISTOGRAM */
#elif defined(COLLMOD_ZML)
	FDL_write(fdl,"num_dust_bins",&msr->param.CP.DB.nDustBins);
	for (i=0;i<msr->param.CP.DB.nDustBins;i++) {
		FDL_index(fdl,"dust_bins_index",i);
		FDL_write(fdl,"dust_bins_mass",&msr->aDustBins[i].dMass);
		FDL_write(fdl,"dust_bins_volume",&msr->aDustBins[i].dVolume);
#ifdef ORIGIN_HISTOGRAM
		/*DEBUG for some reason it's not possible to use num_dust_bins
		  in the FDL file for SIZE of dust_bins_origin_bins_array --
		  it must be hardcoded!*/
		for (j=0;j<NUM_ORIGIN_BINS;j++) {
			FDL_index(fdl,"dust_bins_origin_bins_index",j);
			FDL_write(fdl,"dust_bins_origin_bins_value",&msr->aDustBins[i].origin_bins[j]);
			}
#endif /* ORIGIN_HISTOGRAM */
		}
		
	FDL_write(fdl,"dust_bins_trash",&msr->DustBinsTrash.dMass);
#ifdef ORIGIN_HISTOGRAM
	for (i=0;i<NUM_ORIGIN_BINS;i++) {
		FDL_index(fdl,"trash_bin_origin_bins_index",i);
		FDL_write(fdl,"trash_bin_origin_bins_value",&msr->DustBinsTrash.origin_bins[i]);
		}
#endif /* ORIGIN_HISTOGRAM */
#endif /* RUBBLE_ZML, COLLMOD_ZML */
	/*
	 ** Write the old parameters.
	 */
	FDL_write(fdl,"bPeriodic",&msr->param.bPeriodic);
	FDL_write(fdl,"bComove",&msr->param.csm->bComove);
	FDL_write(fdl,"bParaRead",&msr->param.bParaRead);
	FDL_write(fdl,"bParaWrite",&msr->param.bParaWrite);
	FDL_write(fdl,"bCannonical",&msr->param.bCannonical);
	FDL_write(fdl,"bStandard",&msr->param.bStandard);
	FDL_write(fdl,"bKDK",&msr->param.bKDK);
	FDL_write(fdl,"bBinary",&msr->param.bBinary);
	FDL_write(fdl,"bDoGravity",&msr->param.bDoGravity);
	FDL_write(fdl,"bDoGas",&msr->param.bDoGas);
	FDL_write(fdl,"bFandG",&msr->param.bFandG);
	FDL_write(fdl,"bHeliocentric",&msr->param.bHeliocentric);
	FDL_write(fdl,"bLogHalo",&msr->param.bLogHalo);
	FDL_write(fdl,"bHernquistSpheroid",&msr->param.bHernquistSpheroid);
	FDL_write(fdl,"bMiyamotoDisk",&msr->param.bMiyamotoDisk);
	FDL_write(fdl,"bRotatingBar", &msr->param.bRotatingBar);
	rotbarCheckWrite(msr->param.rotbar, fdl);
	
	FDL_write(fdl,"nBucket",&msr->param.nBucket);
	FDL_write(fdl,"iOutInterval",&msr->param.iOutInterval);
	FDL_write(fdl,"iLogInterval",&msr->param.iLogInterval);
	FDL_write(fdl,"iCheckInterval",&msr->param.iCheckInterval);
	FDL_write(fdl,"iExpOrder",&msr->param.iOrder);
	FDL_write(fdl,"iEwOrder",&msr->param.iEwOrder);
	FDL_write(fdl,"nReplicas",&msr->param.nReplicas);
	FDL_write(fdl,"iStartStep",&msr->param.iStartStep);
	FDL_write(fdl,"nSteps",&msr->param.nSteps);
	FDL_write(fdl,"dExtraStore",&msr->param.dExtraStore);
	FDL_write(fdl,"dDelta",&msr->param.dDelta);
	FDL_write(fdl,"dEta",&msr->param.dEta);
	FDL_write(fdl,"dEtaCourant",&msr->param.dEtaCourant);
	FDL_write(fdl,"bEpsAccStep",&msr->param.bEpsAccStep);
	FDL_write(fdl,"bNonSymp",&msr->param.bNonSymp);
	FDL_write(fdl,"iMaxRung",&msr->param.iMaxRung);
	FDL_write(fdl,"dEwCut",&msr->param.dEwCut);
	FDL_write(fdl,"dEwhCut",&msr->param.dEwhCut);
	FDL_write(fdl,"dPeriod",&msr->param.dPeriod);
	FDL_write(fdl,"dxPeriod",&msr->param.dxPeriod);
	FDL_write(fdl,"dyPeriod",&msr->param.dyPeriod);
	FDL_write(fdl,"dzPeriod",&msr->param.dzPeriod);
	FDL_write(fdl,"dHubble0",&msr->param.csm->dHubble0);
	FDL_write(fdl,"dOmega0",&msr->param.csm->dOmega0);
	FDL_write(fdl,"dLambda",&msr->param.csm->dLambda);
	FDL_write(fdl,"dOmegaRad",&msr->param.csm->dOmegaRad);
	FDL_write(fdl,"dQuintess",&msr->param.csm->dQuintess);
	FDL_write(fdl,"dTheta",&msr->param.dTheta);
	FDL_write(fdl,"dTheta2",&msr->param.dTheta2);
	FDL_write(fdl,"dCentMass",&msr->param.dCentMass);
#if defined(GASOLINE) && !defined(NOCOOLING) && defined(COOLING_COSMO)
	FDL_write(fdl,"dMassFracHelium",&msr->param.CoolParam.dMassFracHelium);
#else
	{
	FLOAT dummy = 0.75; /* Nasty! */
	FDL_write(fdl,"dMassFracHelium",&dummy);
	}
#endif
	FDL_write(fdl,"dFracNoDomainDecomp",&msr->param.dFracNoDomainDecomp);
	if (msr->param.bVDetails) {
		printf("Writing checkpoint file...\nTime:%g\n",dTime);
		sec = msrTime(msr);
	        }
	/*
	 ** Do a parallel or serial write to the output file.
	 */
	msrCalcWriteStart(msr);
	in.iOffset = FDL_offset(fdl,"particle_array");
	if(msr->param.bParaWrite)
	    pstWriteCheck(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteCheck(msr, &in);
	if (msr->param.bVDetails) {
		puts("Checkpoint file has been successfully written.");
		dsec = msrTime(msr) - sec;
		printf("Data write complete, Wallclock: %f secs\n",dsec);
	        }
	iNotCorrupt = 1;
	FDL_write(fdl,"not_corrupt_flag",&iNotCorrupt);
	FDL_finish(fdl);
	}


int msrOutTime(MSR msr,double dTime)
{	
	if (msr->iOut < msr->nOuts) {
		if (dTime >= msr->pdOutTime[msr->iOut]) {
			++msr->iOut;
			return(1);
			}
		else return(0);
		}
	else return(0);
	}


int cmpTime(const void *v1,const void *v2) 
{
	const double *d1 = (const double *)v1;
	const double *d2 = (const double *)v2;

	if (*d1 < *d2) return(-1);
	else if (*d1 == *d2) return(0);
	else return(1);
	}

void msrReadOuts(MSR msr,double dTime)
{
	char achFile[PST_FILENAME_SIZE];
	char ach[PST_FILENAME_SIZE];
	LCL *plcl = &msr->lcl;
	FILE *fp;
	int i,ret;
	double z,a,n;
	char achIn[80];
	
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achFile[0] = '\0';
	snprintf(achFile,PST_FILENAME_SIZE,"%s/%s.red",
			 msr->param.achDataSubPath,msr->param.achOutName);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		strcpy(ach,achFile);
		snprintf(achFile,PST_FILENAME_SIZE,"%s/%s",plcl->pszDataPath,ach);
		}
	fp = fopen(achFile,"r");
	if (!fp) {
#ifndef COLLISIONS
		if (msr->param.bVWarnings)
			printf("WARNING: Could not open redshift input file:%s\n",achFile);
#endif
		msr->nOuts = 0;
		return;
		}
	i = 0;
	while (1) {
		if (!fgets(achIn,80,fp)) goto NoMoreOuts;
		switch (achIn[0]) {
		case 'z':
			ret = sscanf(&achIn[1],"%lf",&z);
			if (ret != 1) goto NoMoreOuts;
			a = 1.0/(z+1.0);
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
			break;
		case 'a':
			ret = sscanf(&achIn[1],"%lf",&a);
			if (ret != 1) goto NoMoreOuts;
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
			break;
		case 't':
			ret = sscanf(&achIn[1],"%lf",&msr->pdOutTime[i]);
			if (ret != 1) goto NoMoreOuts;
			break;
		case 'n':
			ret = sscanf(&achIn[1],"%lf",&n);
			if (ret != 1) goto NoMoreOuts;
			msr->pdOutTime[i] = dTime + (n-0.5)*msrDelta(msr);
			break;
		default:
			ret = sscanf(achIn,"%lf",&z);
			if (ret != 1) goto NoMoreOuts;
			a = 1.0/(z+1.0);
			msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
			}
		++i;
		if(i > msr->nMaxOuts) {
			msr->nMaxOuts *= 2;
			msr->pdOutTime = realloc(msr->pdOutTime,
									 msr->nMaxOuts*sizeof(double));
			assert(msr->pdOutTime != NULL);
		    }
		}
 NoMoreOuts:
	msr->nOuts = i;
	/*
	 ** Now sort the array of output times into ascending order.
	 */
	qsort(msr->pdOutTime,msr->nOuts,sizeof(double),cmpTime);
	fclose(fp);
	}


int msrSteps(MSR msr)
{
	return(msr->param.nSteps);
	}


char *msrOutName(MSR msr)
{
	return(msr->param.achOutName);
	}


double msrDelta(MSR msr)
{
	return(msr->param.dDelta);
	}


int msrCheckInterval(MSR msr)
{
	return(msr->param.iCheckInterval);
	}


int msrLogInterval(MSR msr)
{
	return(msr->param.iLogInterval);
	}


int msrOutInterval(MSR msr)
{
	return(msr->param.iOutInterval);
	}


int msrRestart(MSR msr)
{
	return(msr->param.bRestart);
	}


int msrComove(MSR msr)
{
	return(msr->param.csm->bComove);
	}


int msrKDK(MSR msr)
{
	return(msr->param.bCannonical && msr->param.bKDK);
	}


int msrDoSun(MSR msr)
{
	if (msr->param.bFandG) return(0);
	else return(1);
	}


double msrSoft(MSR msr)
{
	return(msr->param.dSoft);
	}


void msrSwitchTheta(MSR msr,double dTime)
{
	double a;

	a = csmTime2Exp(msr->param.csm,dTime);
	if (a >= msr->param.daSwitchTheta) msr->dCrit = msr->param.dTheta2; 
	}


double msrMassCheck(MSR msr,double dMass,const char *pszWhere)
{
	struct outMassCheck out;

#ifdef GROWMASS
	out.dMass = 0.0;
#else
	if (dMass == -1.0)
		return 0.0; /* nothing to compare against! (use -2.0 for init'n) */
	if (msr->param.bVDetails)
		printf("Mass check (%s)...",pszWhere);
	pstMassCheck(msr->pst,NULL,0,&out,NULL);
#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
	{
	int i;

	for (i=0;i<msr->param.CP.DB.nDustBins;i++) 
		out.dMass += msr->aDustBins[i].dMass;
	}
	out.dMass += msr->DustBinsTrash.dMass;
#endif /* RUBBLE_ZML || COLLMOD_ZML */

#ifdef GR_DRAG
	out.dMass += msr->dMergerMassLost;
#endif

#ifdef WALLS
	out.dMass += msr->dDeathWallsMassLost;
#endif

	if (dMass < 0.0) {
		if (msr->param.bVDetails)
			printf("mass = %.15e\n",out.dMass);
		return(out.dMass);
		}
	else if (fabs(dMass - out.dMass) > CONSERVE_FRAC*dMass) {
		printf("ERROR: Mass not conserved: %.15e != %.15e!",
			   dMass,out.dMass);
		if (msr->param.bVDetails)
			printf("\n");
		else
			printf(" (%s)\n",pszWhere);
		}
	else if (msr->param.bVDetails)
		printf("ok\n");
#endif /* if !(GROWMASS) */
	return(out.dMass);
	}

void msrMassMetalsEnergyCheck(MSR msr,double *dTotMass, double *dTotMetals, 
        double *dTotFe, double *dTotOx, double *dTotEnergy,char *pszWhere)
{
	struct outMassMetalsEnergyCheck out;
	
#ifdef GROWMASS
	out.dTotMass = 0.0;
	out.dTotMetals = 0.0;
	out.dTotFe = 0.0;
	out.dTotOx = 0.0;
	out.dTotEnergy = 0.0;
#else
	if (msr->param.bVDetails)
		puts("doing mass check...");
	pstMassMetalsEnergyCheck(msr->pst,NULL,0,&out,NULL);
	if (*dTotMass < 0.0) {}
	else {
		if ( fabs(*dTotMass - out.dTotMass) > CONSERVE_FRAC*(*dTotMass) ) {
			printf("ERROR: Mass not conserved (%s): %.15e != %.15e!\n",
				   pszWhere,*dTotMass,out.dTotMass);
			}
#ifdef STARFORM
            /* Commented out because metals are almost conserved.  Error
             * comes because some of the gas mass gets converted into stars
             * so conservation error is reported as a net loss in metals
             * Remaining metals are in stars
			 if ( fabs(*dTotMetals - out.dTotMetals) > CONSERVE_FRAC*(*dTotMetals) ) {
			 printf("ERROR: Metal mass not conserved (%s): %.15e != %.15e!\n",
			 pszWhere,*dTotMetals,out.dTotMetals);
			 }*/
                if ( fabs(*dTotMetals - out.dTotMetals) > CONSERVE_FRAC*(*dTotMetals) ) {
			 printf("ERROR: Metal mass not conserved (%s): %.15e != %.15e!\n",
			 pszWhere,*dTotMetals,out.dTotMetals);
			 }
                if ( fabs(*dTotFe - out.dTotFe) > CONSERVE_FRAC*(*dTotFe) ) {
			 printf("ERROR: Iron mass not conserved (%s): %.15e != %.15e!\n",
			 pszWhere,*dTotFe,out.dTotFe);
			 }
                if ( fabs(*dTotOx - out.dTotOx) > CONSERVE_FRAC*(*dTotOx) ) {
			 printf("ERROR: Oxygen mass not conserved (%s): %.15e != %.15e!\n",
			 pszWhere,*dTotOx,out.dTotOx);
			 }
                if ( fabs(out.dTotMetals - out.dTotOx - out.dTotFe) > CONSERVE_FRAC*(out.dTotMetals) ) {
			 printf("ERROR: Oxygen and Iron do not add up to total metals (%s): %.15e != %.15e!\n",
			 pszWhere,(out.dTotOx + out.dTotFe),out.dTotMetals);
			 }
		if ( fabs(*dTotEnergy - out.dTotEnergy*msr->param.dDeltaStarForm) > CONSERVE_FRAC*(*dTotEnergy) ) {
			printf("ERROR: SN Energy not conserved (%s): %.15e != %.15e!\n",
				   pszWhere,*dTotEnergy,out.dTotEnergy*msr->param.dDeltaStarForm);
			}
#endif
		}
#endif
	*dTotMass = out.dTotMass;
	*dTotMetals = out.dTotMetals;
	*dTotFe = out.dTotFe;
	*dTotOx = out.dTotOx;
	*dTotEnergy = out.dTotEnergy;
	return;
	}

void
msrInitStep(MSR msr)
{
    struct inSetRung in;

    in.iRung = msr->param.iMaxRung - 1;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.iRung;
    }


void
msrSetRung(MSR msr, int iRung)
{
    struct inSetRung in;

    in.iRung = iRung;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.iRung;
    }


int msrMaxRung(MSR msr)
{
    return msr->param.iMaxRung;
    }


int msrCurrMaxRung(MSR msr)
{
    return msr->iCurrMaxRung;
    }


int msrCurrMaxRungInclDF(MSR msr)
{
	int iRung = msr->iCurrMaxRung;
	if (msr->bDumpFrame && iRung < msr->df->iMaxRung) iRung = msr->df->iMaxRung;
    return iRung;
    }


double msrEta(MSR msr)
{
    return msr->param.dEta;
    }

double msrEtaCourant(MSR msr)
{
    return msr->param.dEtaCourant;
    }


void msrBallMax(MSR msr, int iRung, int bGreater)
{
    struct inBallMax in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    in.dhFac = 1.0+msr->param.ddHonHLimit;
    pstBallMax(msr->pst, &in, sizeof(in), NULL, NULL);
    }

/*
 * bGreater = 1 => activate all particles at this rung and greater.
 */
void msrActiveRung(MSR msr, int iRung, int bGreater)
{
    struct inActiveRung in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstActiveRung(msr->pst, &in, sizeof(in), &(msr->nActive), NULL);
    }

void msrActiveTypeOrder(MSR msr, unsigned int iTestMask )
{
    struct inActiveTypeOrder in;
    int nActive;

    in.iTestMask = iTestMask;
    pstActiveTypeOrder(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iTestMask & TYPE_ACTIVE)       msr->nActive      = nActive;
    if (iTestMask & TYPE_TREEACTIVE)   msr->nTreeActive   = nActive;
    if (iTestMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveOrder(MSR msr)
{
    pstActiveOrder(msr->pst,NULL,0,&(msr->nActive),NULL);
    }

void msrActiveExactType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iFilterMask = iFilterMask;
    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstActiveExactType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstActiveType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrSetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstSetType(msr->pst,&in,sizeof(in),&nActive,NULL);
    }

void msrResetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;

    pstResetType(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (msr->param.bVDetails) printf("nResetType: %d\n",nActive);
    }

int msrCountType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask) 
{
    struct inActiveType in;
    int nActive;

    in.iFilterMask = iFilterMask;
    in.iTestMask = iTestMask;

    pstCountType(msr->pst,&in,sizeof(in),&nActive,NULL);

    return nActive;
    }

void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = (~0);
    in.iSetMask = iSetMask;
    in.iRung = iRung;
    in.bGreater = bGreater;

    pstActiveMaskRung(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

void msrActiveTypeRung(MSR msr, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater) 
{
    struct inActiveType in;
    int nActive;

    in.iTestMask = iTestMask;
    in.iSetMask = iSetMask;
    in.iRung = iRung;
    in.bGreater = bGreater;

    pstActiveTypeRung(msr->pst,&in,sizeof(in),&nActive,NULL);

    if (iSetMask & TYPE_ACTIVE      ) msr->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) msr->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) msr->nSmoothActive = nActive;
    }

int msrCurrRung(MSR msr, int iRung)
{
    struct inCurrRung in;
    struct outCurrRung out;

    in.iRung = iRung;
    pstCurrRung(msr->pst, &in, sizeof(in), &out, NULL);
    return out.iCurrent;
    }

void
msrGravStep(MSR msr, double dTime)
{
    struct inGravStep in;
    double a;

    if (msrComove(msr)) {
        a = csmTime2Exp(msr->param.csm,dTime);
        in.dEta = msrEta(msr)*pow(a,1.5);
        }
    else {
        in.dEta = msrEta(msr);
        }
    pstGravStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrAccelStep(MSR msr,double dTime)
{
    struct inAccelStep in;
    double a;

    in.dEta = msrEta(msr);
    a = csmTime2Exp(msr->param.csm,dTime);
    if (msr->param.bCannonical) {
		in.dVelFac = 1.0/(a*a);
	}
    else {
		in.dVelFac = 1.0;
	}
    in.dAccFac = 1.0/(a*a*a);
    in.bDoGravity = msrDoGravity(msr);
    in.bEpsAcc = msr->param.bEpsAccStep;
    in.bSqrtPhi = msr->param.bSqrtPhiStep;
    in.dhMinOverSoft = (msr->param.bLowerSoundSpeed ? msr->param.dhMinOverSoft : 0);
    pstAccelStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrDensityStep(MSR msr,double dTime)
{
    struct inDensityStep in;
    double expand;

    if (msr->param.bVDetails) printf("Calculating Rung Densities...\n");
    msrSmooth(msr,dTime,SMX_DENSITY,0);
    in.dEta = msrEta(msr);
    expand = csmTime2Exp(msr->param.csm,dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    pstDensityStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrInitDt(MSR msr)
{
    struct inInitDt in;
    
    in.dDelta = msrDelta(msr);
    pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrDtToRung(MSR msr, int iRung, double dDelta, int bAll)
{
    struct inDtToRung in;
    struct outDtToRung out;

    in.iRung = iRung;
    in.dDelta = dDelta;
    in.iMaxRung = msrMaxRung(msr);
    in.bAll = bAll;

    pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);

    if(out.iMaxRungIdeal > msrMaxRung(msr)) {
	fprintf(stderr, "WARNING, TIMESTEPS TOO LARGE: nMaxRung (%d) is greater than ideal rung (%d)\n", 
		msrMaxRung(msr), out.iMaxRungIdeal);
	}
    if (out.nMaxRung <= msr->param.nTruncateRung && out.iMaxRung > iRung) {
	if (msr->param.bVDetails)
	    printf("n_CurrMaxRung = %d  (iCurrMaxRung = %d):  Promoting particles to iCurrMaxrung = %d\n",
		   out.nMaxRung,out.iMaxRung,out.iMaxRung-1);

	in.iMaxRung = out.iMaxRung; /* Note this is the forbidden rung so no -1 here */
	pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);
	}

    msr->iCurrMaxRung = out.iMaxRung;
    }

/* Not fully tested: */
void msrTopStepSym(MSR msr, double dStep, double dTime, double dDelta, 
				   int iRung, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax, dIMax, dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
	    if(iRung < msrMaxRung(msr)-1) {
			if (msr->param.bVDetails) printf("Adjust, iRung: %d\n", iRung);
			msrActiveRung(msr, iRung, 1);
			msrDrift(msr,dTime,0.5*dDelta);
			dTime += 0.5*dDelta;
			msrInitDt(msr);
			if (msr->param.bGravStep || msr->param.bAccelStep) {
			    msrInitAccel(msr);
			    msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrUpdateSoft(msr,dTime);
			    msrBuildTree(msr,0,dMass,0);
			    msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
				if (msr->param.bGravStep) {
					msrGravStep(msr,dTime);
					}
				if (msr->param.bAccelStep) {
					msrAccelStep(msr,dTime);
					}
			    }
			if (msr->param.bDensityStep) {
			    msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrBuildTree(msr,0,dMass,1);
			    msrDensityStep(msr,dTime);
			    }
			msrDtToRung(msr,iRung,dDelta,0);
			msrDrift(msr,dTime,-0.5*dDelta);
			dTime -= 0.5*dDelta;
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepSym(msr, dStep, dTime, 0.5*dDelta,iRung+1,pdActiveSum);
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr, iRung, 0);
		msrInitAccel(msr);
#ifdef GASOLINE
		if(msrSphCurrRung(msr, iRung, 0)) {
			if (msr->param.bVDetails) printf("SPH, iRung: %d\n", iRung);
           	        msrActiveTypeRung(msr, TYPE_GAS, TYPE_ACTIVE, iRung, 0 );
                        msrDomainDecomp(msr,iRung,0);
           	        msrActiveTypeRung(msr, TYPE_GAS, TYPE_ACTIVE, iRung, 0 );
           	        msrActiveType(msr, TYPE_GAS, TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
			msrBuildTree(msr,1,-1.0,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrGetGasPressure(msr);
			if (msrDoGas(msr)) {
			  if (msr->param.bBulkViscosity) {
  			    msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			    msrActiveRung(msr, iRung, 0);
			    msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
			    }
			  else {
                            if (msr->param.bViscosityLimiter
				|| msr->param.bShockTracker
				|| msr->param.bStarForm) 
			      msrReSmooth(msr, dTime, SMX_DIVVORT, 1);
			    msrUpdateShockTracker(msr, dDelta);
  			    msrSphViscosityLimiter(msr, dTime);
			    msrActiveRung(msr, iRung, 0);
			    msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
                            }
			  }
			}
#endif
		if(msrCurrRung(msr, iRung)) {
		    if(msrDoGravity(msr)) {
   		        if (msr->param.bVDetails) printf("Gravity, iRung: %d\n", iRung);
                        msrDomainDecomp(msr,iRung,0);
			msrActiveRung(msr, iRung, 0);
			msrUpdateSoft(msr,dTime);
			msrBuildTree(msr,0,dMass,0);
			msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
			*pdActiveSum += (double)nActive/msr->N;
			}
		    }
		if (msr->param.bVDetails) printf("Kick, iRung: %d\n", iRung);
		msrKickDKD(msr, dTime, dDelta);
		msrRungStats(msr);
		msrTopStepSym(msr,dStep,dTime+0.5*dDelta,0.5*dDelta,iRung+1,
					  pdActiveSum);
		}
	else {    
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n",iRung-1);
		msrDrift(msr,dTime,dDelta);
		}
	}


void msrRungStats(MSR msr)
{
	if (msr->param.bVRungStat) {
		struct inRungStats in;
		struct outRungStats out;
		int i;

		printf("Rung distribution: (");
		for (i=0;i<msr->param.iMaxRung;++i) {
			if (i != 0) printf(",");
			in.iRung = i;
			pstRungStats(msr->pst,&in,sizeof(in),&out,NULL);
			msr->nRung[i] = out.nParticles;
			printf("%d",out.nParticles);
			}
		printf(")\n");
		}
	}

void msrTopStepNS(MSR msr, double dStep, double dTime, double dDelta, int
				  iRung, int iAdjust, double *pdActiveSum)
{
    double dMass = -1.0;
    int iSec;
    double dWMax,dIMax,dEMax;
	int nActive;

	if(msrCurrMaxRung(msr) >= iRung) { /* particles to be kicked? */
		if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
			if (msr->param.bVDetails) printf("Adjust, iRung: %d\n", iRung);
			msrActiveRung(msr,iRung,1);
 			msrActiveType(msr,TYPE_ALL,TYPE_SMOOTHACTIVE|TYPE_TREEACTIVE);
			msrInitDt(msr);
			if (msr->param.bGravStep) {
				msrGravStep(msr,dTime);
				}
			if (msr->param.bAccelStep) {
			    msrAccelStep(msr,dTime);
				}
			if (msr->param.bDensityStep) {
				msrDomainDecomp(msr,iRung,1);
			    msrActiveRung(msr,iRung,1);
			    msrBuildTree(msr,0,dMass,1);
			    msrDensityStep(msr,dTime);
			    }
#ifdef GASOLINE
			if (msr->param.bSphStep) {
				msrSphStep(msr,dTime);
				}
#endif
			msrDtToRung(msr,iRung,dDelta,1);
			}
		/*
		 ** Actual Stepping.
		 */
		msrTopStepNS(msr,dStep,dTime,0.5*dDelta,iRung+1,0,pdActiveSum);
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr,iRung,0);
		msrDomainDecomp(msr,iRung,0);
		msrActiveRung(msr,iRung,0);
		msrInitAccel(msr);
#ifdef GASOLINE
		if(msrSphCurrRung(msr, iRung, 0)) {
			if (msr->param.bVDetails) printf("SPH, iRung: %d\n", iRung);
			msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iRung,0);
			msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			msrBuildTree(msr,1,-1.0,1);
			msrSmooth(msr,dTime,SMX_DENSITY,1);
			msrGetGasPressure(msr);
			if (msrDoGas(msr)) {
 			   if (msr->param.bBulkViscosity) {
  			       msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			       msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);
			       } 
			   else {
			       if (msr->param.bViscosityLimiter
				   || msr->param.bShockTracker
				   || msr->param.bStarForm)
				   msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			       msrUpdateShockTracker(msr, dDelta);
			       msrSphViscosityLimiter(msr, dTime);
			       msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1);
			       }
			   }
		}
#endif
		if(msrCurrRung(msr, iRung)) {
		    if(msrDoGravity(msr)) {
				if (msr->param.bVDetails) printf("Gravity, iRung: %d\n", iRung);
				msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
				msrUpdateSoft(msr,dTime);
				msrBuildTree(msr,0,dMass,0);
				msrGravity(msr,dStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,
						   &nActive);
				*pdActiveSum += (double)nActive/msr->N;
				}
		    }
		if (msr->param.bVDetails) printf("Kick, iRung: %d\n", iRung);
		msrKickDKD(msr,dTime,dDelta);
		msrTopStepNS(msr,dStep,dTime+0.5*dDelta,0.5*dDelta,iRung+1,1,
					 pdActiveSum);
		}
	else {    
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n", iRung-1);
		msrDrift(msr,dTime,dDelta);
		}
	}

void msrTopStepDKD(MSR msr, double dStep, double dTime, double dDelta, 
				double *pdMultiEff)
{
	int iRung = 0;

#ifdef COLLISIONS
	assert(0); /* DKD multi-stepping unsupported for COLLISIONS */
#endif

	*pdMultiEff = 0.0;
	if(msr->param.bNonSymp)
		msrTopStepNS(msr,dStep,dTime,dDelta,iRung,1,pdMultiEff);
	else
		msrTopStepSym(msr,dStep,dTime,dDelta,iRung,pdMultiEff);

	if (msr->param.bVStep)
		printf("Multistep Efficiency (average number of microsteps per step):%f\n",
			   *pdMultiEff);
	}

void msrTopStepKDK(MSR msr,
				   double dStep,	/* Current step */
				   double dTime,	/* Current time */
				   double dDelta,	/* Time step */
				   int iRung,		/* Rung level */
				   int iKickRung,	/* Gravity on all rungs from iRung
									   to iKickRung */
				   int iAdjust,		/* Do an adjust? */
				   double *pdActiveSum,
				   double *pdWMax,
				   double *pdIMax,
				   double *pdEMax,
				   int *piSec)
{
    double dMass = -1.0;
    int nActive;

    LogTimingSetRung( msr, iKickRung );
    if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
		if (msr->param.bVDetails) printf("Adjust, iRung: %d\n",iRung);
		msrActiveRung(msr, iRung, 1);
		msrInitDt(msr);
		if (msr->param.bGravStep) {
			msrGravStep(msr,dTime);
			}
		if (msr->param.bAccelStep) {
		    msrAccelStep(msr,dTime);
			}
		if (msr->param.bDensityStep) {
		    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		    msrDomainDecomp(msr,iRung,1);
		    msrActiveRung(msr,iRung,1);
		    msrBuildTree(msr,0,dMass,1);
		    msrDensityStep(msr,dTime);
		    }
		if (msr->param.bDeltaAccelStep) {
		  /* Ensure we have a Density tree ready to use */
			if (!msr->param.bDeltaAccelStepGasTree || 
				msr->iTreeType != MSR_TREE_DENSITY) {
#ifdef DELTAACCELACTIVE
			    msrActiveRung(msr,iRung,1);
			    msrActiveType(msr,TYPE_ACTIVE,TYPE_TREEACTIVE);
		        msrBuildTree(msr,1,-1.0,1);   /* bTreeActive */
#else
			    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
			    msrBuildTree(msr,1,-1.0,1);   /* bTreeActive */
#endif
			    }
		    msrActiveRung(msr,iRung,1);
			msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE);
			msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE);
			/* This smooth sets dt directly -- hardwired coefficient */
			msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
		    }

#ifdef GASOLINE
		if (msr->param.bSphStep) {
			msrSphStep(msr,dTime);
			}
#endif

#ifdef RUBBLE_ZML
		msrRubbleStep(msr);
#elif defined(COLLMOD_ZML)
		msrCollModStep(msr);
#endif

		msrDtToRung(msr,iRung,dDelta,1);
		if (iRung == 0) {
		  /*
		  msrReorder(msr);
		  msrOutArray(msr,"test.dt",OUT_DT_ARRAY);
		  msrActiveOrder(msr);
		  */
		  msrRungStats(msr);
		  }
		}
    if (msr->param.bVDetails) printf("Kick, iRung: %d\n",iRung);
    msrActiveRung(msr,iRung,0);
#ifdef GASOLINE
    msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
    msrUpdateuDot(msr,dTime,0.5*dDelta,1);
#endif
#ifdef GR_DRAG
	msrGRIntegrateCloseParticles(msr,dDelta,dTime);
#endif

#ifdef SPINUP
#ifndef SPINUP_WITH_AGGS
	msrSpinup(msr, dStep); //If SPINUP_WITH_AGGS is on, call msrSpinupAccel() instead.
#endif
#endif

    msrKickKDKOpen(msr,dTime,0.5*dDelta);
    if (msrCurrMaxRungInclDF(msr) > iRung) {
		/*
		 ** Recurse.
		 */
		msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,0,
					  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
                /* Call to TopStep can change the rung setting so redo Set Rung */
		LogTimingSetRung( msr, iKickRung );
		dStep += 1.0/(2 << iRung);
		dTime += 0.5*dDelta;
		msrActiveRung(msr,iRung,0);
#ifdef GASOLINE
		msrUpdateuDot(msr,dTime,0.5*dDelta,0); /* Need forward uDot for Upreds */
#endif
		msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,1,
					  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
		LogTimingSetRung( msr, iKickRung );
		}
    else {
		/* This Drifts everybody */
		if (msr->param.bVDetails) printf("Drift, iRung: %d\n", iRung);
#ifdef GASOLINE
		msrDrift(msr,dTime,0.5*dDelta);
		dTime += 0.5*dDelta;
		dStep += 1.0/(2 << iRung);
		msrActiveRung(msr,iRung,0);
		msrUpdateuDot(msr,dTime,0.5*dDelta,0); /* Need forward uDot for uPred */
		msrDrift(msr,dTime,0.5*dDelta);
		dTime += 0.5*dDelta;
		dStep += 1.0/(2 << iRung);
#else
		msrDrift(msr,dTime,dDelta);
#ifdef RUBBLE_ZML
		/*
		 ** Since we may need to turn planetesimals to rubble in the
		 ** middle of the step, and rubble pieces are forced to be on
		 ** the lowest rung, we may need to recompute the kick
		 ** circumstances now and proceed with any colliding
		 ** planetesimals forced to the lowest rung. This should ONLY
		 ** happen if dDelta is equal to msr->param.dDelta, so there
		 ** shouldn't be any recursion problems. However, if, in the
		 ** future, we decide to have a more continuous distribution of
		 ** rungs for the rubble problem, this quick fix will need to
		 ** be carefully reconsidered to see if there are any
		 ** unfortunate repercussions.
		 */

		if (msr->param.CP.bDoRubbleKDKRestart) {
			assert(dDelta == msr->param.dDelta);
			assert(iRung == 0);
			msr->param.CP.bDoRubbleKDKRestart = 0;
			msrKickKDKOpen(msr,dTime,-0.5*dDelta); /* go to start of kick */
			msrTopStepKDK(msr,dStep,dTime,dDelta,iRung,iKickRung,1,
						  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
			return;
			}
#elif defined(COLLMOD_ZML)
		/*
		 ** Although there are no resolved rubble piles the time step 
		 ** will still need to be reduced to the bottom rung in order
		 ** to resolve the collisions.
		 */

		if (msr->param.CP.bDoCollModKDKRestart) {
			assert(dDelta == msr->param.dDelta);
			assert(iRung == 0);
			msr->param.CP.bDoCollModKDKRestart = 0;
			msrKickKDKOpen(msr,dTime,-0.5*dDelta); /* go to start of kick */
			msrTopStepKDK(msr,dStep,dTime,dDelta,iRung,iKickRung,1,
						  pdActiveSum,pdWMax,pdIMax,pdEMax,piSec);
			return;
			}
#endif /* RUBBLE_ZML, COLLMOD_ZML */

		dTime += dDelta;
		dStep += 1.0/(1 << iRung);
#endif
#ifdef SIMPLESF
		msrSimpleStarForm(msr, dTime, dDelta);
#endif
#ifdef STARFORM
                /* only form stars at user defined intervals */
                /* JW: Is this dDelta choice correct? */
                if ( iKickRung <= msr->param.iStarFormRung )
                    msrFormStars(msr, dTime, max(dDelta,msr->param.dDeltaStarForm));
#endif
		/* 
		 ** Dump Frame
		 */
		if (msr->param.dDumpFrameTime > 0 && dTime >= msr->df->dTime)
			msrDumpFrame( msr, dTime, dStep );
		else if (msr->param.dDumpFrameStep > 0 && dStep >= msr->df->dStep) 
			msrDumpFrame( msr, dTime, dStep );

		/* 
		 ** Calculate Forces (if required)
		 */
#if 0
		msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
#endif
		msrActiveMaskRung(msr,TYPE_ACTIVE,iKickRung,1);
		LogTimingSetN( msr, msr->nActive );

		if (msr->nActive) {
			msrDomainDecomp(msr,iKickRung,1);
			msrInitAccel(msr);

			if (msr->param.bVStep) printf("Forces, Step:%f nActive %i\n",dStep,msr->nActive);
			if(msrDoGravity(msr)) {
				if (msr->param.bDoSelfGravity) {
					msrActiveRung(msr,iKickRung,1);
					msrUpdateSoft(msr,dTime);
					msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
					if (msr->param.bVDetails)
						printf("Gravity, iRung: %d to %d\n", iRung, iKickRung);
					msrBuildTree(msr,0,dMass,0);
					}
				msrGravity(msr,dStep,msrDoSun(msr),piSec,pdWMax,pdIMax,pdEMax,&nActive);
#ifdef CHECKSOFT			  
	   {
	   char achFile[256]; 

  	   fprintf(stderr,"Outputing .soft .dt .den tipsy\n");
	   msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE|TYPE_DensZeroed);
	   msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iKickRung,1);
	   msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	   msrBuildTree(msr,1,-1.0,1);
	   msrActiveType(msr,TYPE_ACTIVE,TYPE_DensACTIVE );
	   msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1);
	   /*	   msrSmooth(msr,dTime,SMX_DENSITY,1);*/
	   msrReorder(msr);
	   sprintf(achFile,"step%015.10f.soft",dTime);
	   msrOutArray(msr,achFile,OUT_SOFT_ARRAY);
	   sprintf(achFile,"step%015.10f.dt",dTime);
	   msrOutArray(msr,achFile,OUT_DT_ARRAY);
	   sprintf(achFile,"step%015.10f.pot",dTime);
	   msrOutArray(msr,achFile,OUT_POT_ARRAY);
	   sprintf(achFile,"step%015.10f.den",dTime);
	   msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
	   sprintf(achFile,"step%015.10f",dTime);
	   msrWriteTipsy(msr,achFile,dTime);
	   }
#endif
				*pdActiveSum += (double)nActive/msr->N;
				}
			
#ifdef GASOLINE
			if (msr->param.bVDetails)
				printf("SPH, iRung: %d to %d\n",iRung,iKickRung);
			msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iKickRung,1);
			msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
			printf("nActive %d nTreeActive %d nSmoothActive %d\n",msr->nActive,
				   msr->nTreeActive,msr->nSmoothActive);
			if(msrDoGas(msr) && msrSphCurrRung(msr,iKickRung,1)) {
				msrBuildTree(msr,1,-1.0,1);
 				msrActiveTypeRung(msr,TYPE_GAS,TYPE_ACTIVE,iKickRung,1);
				if (msr->param.bFastGas && msr->nActive < msr->nGas*msr->param.dFracFastGas) {
			        msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE|TYPE_Scatter|TYPE_DensZeroed );
					msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE|TYPE_DensACTIVE );
					if (msr->param.bVDetails)
						printf("Dens Active Particles: %d\n",msr->nSmoothActive );
					/* Density for Actives and mark Gather neighbours */
					msrSmooth(msr,dTime,SMX_MARKDENSITY,1); 
					/* mark Scatter Neighbours */
					msrMarkSmooth(msr,dTime,1,TYPE_Scatter); 
					/* They need density too... */
					msrActiveType(msr,TYPE_ACTIVE|TYPE_NbrOfACTIVE|TYPE_Scatter, TYPE_DensACTIVE );
					/* ...but don't redo Actives in smooth*/
					msrActiveExactType(msr,TYPE_DensACTIVE|TYPE_ACTIVE, 
									   TYPE_DensACTIVE,TYPE_SMOOTHACTIVE);
					/* Density for Neighbours and mark scatter neighbours of actives */
					msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1);
					/* mark Scatter Neighbours of Neighbours */
					msrMarkSmooth(msr,dTime,1,TYPE_Scatter);
					/* Scatter Density contribution from those particles... */
					msrActiveExactType(msr,TYPE_Scatter|TYPE_DensACTIVE, 
									   TYPE_Scatter,TYPE_SMOOTHACTIVE);
					msrSmooth(msr,dTime,SMX_MARKIIDENSITY,1); 
					
					/* We want direct neighbours of Actives only */
					msrActiveType(msr,TYPE_NbrOfACTIVE,TYPE_SMOOTHACTIVE);
					if (msr->param.bVDetails)
						printf("Density Zeroed: %d ",
							   msrCountType(msr,TYPE_DensZeroed,TYPE_DensZeroed));
					if (msr->param.bVDetails)
						printf("Neighbours: %d ",
							   msrCountType(msr,TYPE_ACTIVE|TYPE_NbrOfACTIVE,TYPE_NbrOfACTIVE));
			        }
				else {
					msrResetType(msr,TYPE_GAS,TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE );
					msrActiveType(msr,TYPE_ACTIVE,TYPE_SMOOTHACTIVE|TYPE_DensACTIVE );
					if (msr->param.bVDetails)
						printf("Dens Active Particles: %d\n",msr->nSmoothActive );
					msrActiveType(msr,TYPE_GAS,TYPE_SMOOTHACTIVE );
					
					msrSmooth(msr,dTime,SMX_MARKDENSITY,1);
					
					if (msr->param.bVDetails)
						printf("Neighbours: %d ",
							   msrCountType(msr,TYPE_NbrOfACTIVE,TYPE_NbrOfACTIVE ) );
					msrActiveType(msr,TYPE_NbrOfACTIVE,TYPE_SMOOTHACTIVE);
			        }

				if (msr->param.bVDetails)
					printf("Smooth Active Particles: %d\n",msr->nSmoothActive);

				if (msr->param.bViscosityLimiter
					|| msr->param.bBulkViscosity
					|| msr->param.bShockTracker
					|| msr->param.bStarForm) {
					msrReSmooth(msr,dTime,SMX_DIVVORT,1);
					}
				msrSphViscosityLimiter(msr, dTime);
				
				msrGetGasPressure(msr);
				
				if (msr->param.bShockTracker) { 
			        msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
					msrUpdateShockTracker(msr, dDelta);
			        if (msr->param.bBulkViscosity) 
						msrReSmooth(msr,dTime,SMX_HKVISCOSITY,1);
					else
						msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
			        }
				else {
			        if (msr->param.bBulkViscosity) 
						msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);     
					else
						msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1); 
			        }

				msrBallMax(msr,iKickRung,1);
				}
#endif /* GASOLINE */
			}
                /* only accrete onto sinks at user defined intervals 
		   Can I do this after the gas kick and save a treebuild? */
                if ( iKickRung <= msr->param.iSinkRung )
                    msrDoSinks(msr, max(dDelta,msr->param.dDeltaSink) );

		}
    if (msr->param.bVDetails) printf("Kick, iRung: %d\n",iRung);
    msrActiveRung(msr,iRung,0);
#ifdef GASOLINE
    msrUpdateuDot(msr,dTime,0.5*dDelta,1);
#endif
    msrKickKDKClose(msr,dTime,0.5*dDelta);
#ifdef JOHNNY
	msrFindEscapers(msr);
	msrCheckForKepler(msr);
#endif
	}

int
msrMaxOrder(MSR msr)
{
    return msr->nMaxOrder;
    }

void
msrAddDelParticles(MSR msr)
{
    struct outColNParts *pColNParts;
    int *pNewOrder;
    struct inSetNParts in;
    struct inSetParticleTypes intype;
    int iOut;
    int i;
    
    if (msr->param.bVDetails) printf("Changing Particle number\n");
    pColNParts = malloc(msr->nThreads*sizeof(*pColNParts));
    pstColNParts(msr->pst, NULL, 0, pColNParts, &iOut);
    /*
     * Assign starting numbers for new particles in each processor.
     */
    pNewOrder = malloc(msr->nThreads*sizeof(*pNewOrder));
    for(i=0;i<msr->nThreads;i++) {
		/*
		 * Detect any changes in particle number, and force a tree
		 * build.
		 */
		if (pColNParts[i].nNew != 0 || pColNParts[i].nDeltaGas != 0 ||
			pColNParts[i].nDeltaDark != 0 || pColNParts[i].nDeltaStar != 0)
			msr->iTreeType = MSR_TREE_NONE;
		pNewOrder[i] = msr->nMaxOrder + 1;
		msr->nMaxOrder += pColNParts[i].nNew;
		msr->nGas += pColNParts[i].nDeltaGas;
		msr->nDark += pColNParts[i].nDeltaDark;
		msr->nStar += pColNParts[i].nDeltaStar;
		}
    msr->N = msr->nGas + msr->nDark + msr->nStar;
#ifndef GASOLINE
    msr->nMaxOrderDark = msr->nMaxOrder;
#endif

    pstNewOrder(msr->pst,pNewOrder,sizeof(*pNewOrder)*msr->nThreads,NULL,NULL);

    if (msr->param.bVDetails)
	printf("New numbers of particles: %d gas %d dark %d star\n",
	       msr->nGas, msr->nDark, msr->nStar);
    
    in.nGas = msr->nGas;
    in.nDark = msr->nDark;
    in.nStar = msr->nStar;
    in.nMaxOrderGas = msr->nMaxOrderGas;
    in.nMaxOrderDark = msr->nMaxOrderDark;
    in.nMaxOrder = msr->nMaxOrder;
    pstSetNParts(msr->pst,&in,sizeof(in),NULL,NULL);
    intype.nSuperCool = msr->param.nSuperCool;
	/* This shouldn't really be necessary -- it is undesirable to do a fix-up like this */
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL); 

    i = msrCountType(msr, TYPE_GAS, TYPE_GAS);
    assert(i == msr->nGas);
    i = msrCountType(msr, TYPE_DARK, TYPE_DARK);
    assert(i == msr->nDark);
    i = msrCountType(msr, TYPE_STAR, TYPE_STAR);
    assert(i == msr->nStar);

    free(pNewOrder);
    free(pColNParts);
    }

void
msrDoSinks(MSR msr, double dDelta)
{
	double sec,sec1,dsec,dMass;
	int nAccreted;

    if(msr->param.bDoSinks == 0 || msr->nSink == 0) return;
    if (msr->param.bBHSink && dDelta <= 0.0) return;
	    
    sec = msrTime(msr);

    dMass = msrMassCheck(msr, -2.0, "Accrete onto Sinks: Initial Value");

    /* Note: Only gas particles are accreted by sinks */
    if (msr->iTreeType != MSR_TREE_DENSITY) {
	    msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
		msrBuildTree(msr,1,-1.0,1);  /* bTreeActive */
	    }

    msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
    msrActiveType(msr,TYPE_SINK,TYPE_ACTIVE|TYPE_SMOOTHACTIVE);

    msr->param.dSinkCurrentDelta = dDelta;
    if (msr->param.bBHSink) {
        /* Smooth Bondi-Hoyle Accretion: radius set by nSmooth */
	msrSmooth(msr,0.0,SMX_BHSINKACCRETE,1);
	}
    else {
	/* Fixed Radius Accretion: particle by particle (cf. Bate) */
	msrSmooth(msr,0.0,SMX_SINKACCRETE,1);
	}
    
    nAccreted = msr->nGas;

    msrMassCheck(msr, dMass, "Accrete onto Sinks: before particle adjustment");

    msrAddDelParticles(msr);
    msrMassCheck(msr, dMass, "Accrete onto Sinks: after particle adjustment");

	nAccreted -= msr->nGas;

	sec1 = msrTime(msr);
	dsec = sec1 - sec;
	printf("Sinks Done (%d accreted) Calculated, Wallclock: %f secs\n\n",nAccreted,dsec);
	LOGTIMINGUPDATE( dsec, TIMING_Sink );
    }

/* In principle this code is general for any search but for now
   it will be restricted to looking for a nearby star particle */
void
msrCoolUsingParticleList(MSR msr )
{
#ifndef NOCOOLING
    double sec,sec1,dsec;	
	struct inSoughtParticleList in;
	struct inoutParticleList list;
	int i;
  
	in.iTypeSought = TYPE_STAR;
	in.nMax = MAXSOUGHTPARTICLELIST;

	sec = msrTime(msr);
	
	/* O(N_seek N_sought) */
	pstSoughtParticleList(msr->pst,&in,sizeof(in),&list,NULL);
	if (list.n > in.nMax) {
	  fprintf(stderr," Sought Particles returned more than Max (%d > %d)\n",list.n,in.nMax);
	  assert(list.n <= in.nMax);
	}
	for (i=0;i<list.n;i++) {
	  printf("star %d: %g %g %g\n",i,list.p[i].x,list.p[i].y,list.p[i].z);
	}

	pstCoolUsingParticleList(msr->pst,&list,sizeof(list),NULL,NULL);
	

	sec1 = msrTime(msr);
	dsec = sec1 - sec;
	printf("Cooling Radii Calculated (%d Stars), Wallclock: %f secs\n\n",list.n,dsec);
	LOGTIMINGUPDATE( dsec, TIMING_Cool );
	
#endif
    }

int msrDoDensity(MSR msr)
{
	return(msr->param.bDoDensity);
	}

int msrDoGravity(MSR msr)
{
	return(msr->param.bDoGravity);
	}

int msrDoGas(MSR msr)
{
	return(msr->param.bDoGas);
	}

void msrInitAccel(MSR msr)
{
	pstInitAccel(msr->pst,NULL,0,NULL,NULL);
	}

void msrInitTimeSteps(MSR msr,double dTime,double dDelta) 
{
	double dMass = -1.0;

	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	msrInitDt(msr);
	if (msr->param.bGravStep) {
		msrGravStep(msr,dTime);
		}
	if (msr->param.bAccelStep) {
		msrAccelStep(msr,dTime);
		}
	if (msr->param.bDensityStep) {
		msrDomainDecomp(msr,0,1);
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		msrBuildTree(msr,0,dMass,1);
		msrDensityStep(msr,dTime);
		}
#ifdef GASOLINE
	if (msr->param.bSphStep) {
		msrSphStep(msr,dTime);
		}
#endif

#ifdef RUBBLE_ZML
	msrRubbleStep(msr);
#elif defined(COLLMOD_ZML)
	msrCollModStep(msr);
#endif

	msrDtToRung(msr,0,dDelta,1);
	msrRungStats(msr);
	}

#ifdef GASOLINE

void msrGetGasPressure(MSR msr)
{
	struct inGetGasPressure in;
  
	in.iGasModel = (enum GasModel) msr->param.iGasModel;

	switch (in.iGasModel) {

	case GASMODEL_ADIABATIC:
	case GASMODEL_ISOTHERMAL:
	case GASMODEL_COOLING:
		in.gamma = msr->param.dConstGamma;
		in.gammam1 = in.gamma-1;
		break;
	case GASMODEL_GLASS:
#ifdef GLASS
		in.dGlassPoverRhoL = msr->param.dGlassPoverRhoL;
		in.dGlassPoverRhoR = msr->param.dGlassPoverRhoR;
		in.dGlassxL = msr->param.dGlassxL;
		in.dGlassxR = msr->param.dGlassxR;
		in.dxBoundL = -0.5*msr->param.dxPeriod;
		in.dxBoundR = +0.5*msr->param.dxPeriod;
#else
		assert(0);
#endif
		break;
		}

	pstGetGasPressure(msr->pst,&in,sizeof(in),NULL,NULL);

	if (msr->param.bLowerSoundSpeed) msrLowerSoundSpeed(msr);
	}

void msrLowerSoundSpeed(MSR msr)
{
	struct inLowerSoundSpeed in;
  
	in.dhMinOverSoft = msr->param.dhMinOverSoft;

	pstLowerSoundSpeed(msr->pst,&in,sizeof(in),NULL,NULL);
	}


void msrUpdateuDot(MSR msr,double dTime,double dDelta,int bUpdateY)
{
	struct inUpdateuDot in;
	struct outUpdateuDot out;
	double a;
	
#if defined(COOLING_DISK) 
	switch (msr->param.iGasModel) {
	case GASMODEL_COOLING:
	  msrCoolUsingParticleList( msr );
	  break;
	}
#endif

	in.duDelta = dDelta;
	dTime += dDelta/2.0;
	a = csmTime2Exp(msr->param.csm,dTime);
	in.z = 1/a - 1;
	in.dTime = dTime;
	in.iGasModel = msr->param.iGasModel;
	in.bUpdateY = bUpdateY;

	pstUpdateuDot(msr->pst,&in,sizeof(in),&out,NULL);

	printf("UpdateUdot: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
	LOGTIMINGUPDATE( out.MaxTime, TIMING_Cool );
	}

void msrUpdateShockTracker(MSR msr,double dDelta)
{
        struct inUpdateShockTracker in;

        if (!msr->param.bShockTracker) return;

	in.dDelta = dDelta;
	in.dShockTrackerA = msr->param.dShockTrackerA;
	in.dShockTrackerB = msr->param.dShockTrackerB;

	pstUpdateShockTracker(msr->pst,&in,sizeof(struct inUpdateShockTracker),NULL,NULL);

	printf("UpdateShockTracker Done.\n");
	}

void msrInitSph(MSR msr,double dTime)
{
#ifndef NOCOOLING
	struct inInitEnergy in;
	double a;
#endif

	msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
	msrBuildTree(msr,1,-1.0,1);
	msrSmooth(msr,dTime,SMX_DENSITY,1);

#ifndef NOCOOLING
	switch (msr->param.iGasModel) {
	case GASMODEL_COOLING:
	    if(msr->param.bRestart) break;  /* Already OK from checkpoint */
            /*
            * Get a consistent initial state where energy is consistent with 
            * the initial density and input temperature and the ionization
            * fraction is the consistent equilibrium state.
            **/
            in.dTuFac = msr->param.dGasConst/(msr->param.dConstGamma - 1)/
                    msr->param.dMeanMolWeight;
            a = csmTime2Exp(msr->param.csm,dTime);
            in.z = 1/a - 1;
            in.dTime = dTime;
            pstInitEnergy(msr->pst, &in, sizeof(in), NULL, NULL);
            break;
            }
#endif

	if (msrDoGas(msr)) {
	    msrBallMax(msr, 0, 1);
	    if (msr->param.bViscosityLimiter || msr->param.bBulkViscosity
		    || msr->param.bStarForm) {
		        msrReSmooth(msr,dTime,SMX_DIVVORT,1);
			}
		msrSphViscosityLimiter(msr, dTime);

		msrGetGasPressure(msr);
			
		if (msr->param.bShockTracker) { 
			msrReSmooth(msr,dTime,SMX_SPHPRESSURE,1);
			msrUpdateShockTracker(msr, 0.0);
			if (msr->param.bBulkViscosity) 
			        msrReSmooth(msr,dTime,SMX_HKVISCOSITY,1);
			else
			        msrReSmooth(msr,dTime,SMX_SPHVISCOSITY,1);
		        }
		else {
			if (msr->param.bBulkViscosity) 
				msrReSmooth(msr,dTime,SMX_HKPRESSURETERMS,1);     
			else
				msrReSmooth(msr,dTime,SMX_SPHPRESSURETERMS,1); 
		        }

   	        msrUpdateuDot(msr,dTime,0.5*msr->param.dDelta,0);
		}

	}

#ifndef NOCOOLING
void msrInitCooling(MSR msr)
{
	int cntTable;
	int nTableRows, nTableColumns;
	void *dTableData = NULL;
	char TableFileSuffix[20];
	struct inInitCooling in;

	in.dGmPerCcUnit = msr->param.dGmPerCcUnit;
	in.dComovingGmPerCcUnit = msr->param.dComovingGmPerCcUnit;
	in.dErgPerGmUnit = msr->param.dErgPerGmUnit;
	in.dSecUnit = msr->param.dSecUnit;
	in.dKpcUnit = msr->param.dKpcUnit;
	in.z = 60.0; /*dummy value*/
	in.dTime = 0.0; /* dummy value */
	in.CoolParam = msr->param.CoolParam;

	pstInitCooling(msr->pst,&in,sizeof(struct inInitCooling),NULL,NULL);

	/* Read in tables from files as necessary */
	cntTable = 0;
	for (;;) {
		CoolTableReadInfo( &msr->param.CoolParam, cntTable, &nTableColumns, TableFileSuffix );
		if (!nTableColumns) break;

		cntTable++;
		nTableRows = msrReadASCII(msr, TableFileSuffix, nTableColumns, NULL);
		if (nTableRows) {
			assert( sizeof(double)*nTableRows*nTableColumns <= CL_NMAXBYTETABLE );
			dTableData = malloc(sizeof(double)*nTableRows*nTableColumns);
			assert( dTableData != NULL );
			nTableRows = msrReadASCII(msr,TableFileSuffix, 7, dTableData );
			
			pstCoolTableRead(msr->pst,dTableData,sizeof(double)*nTableRows*nTableColumns,NULL,NULL);
			}
		}
	}
#endif

int msrSphCurrRung(MSR msr, int iRung, int bGreater)
{
    struct inSphCurrRung in;
    struct outSphCurrRung out;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstSphCurrRung(msr->pst,&in,sizeof(in),&out,NULL);
    return out.iCurrent;
    }

void msrSphStep(MSR msr, double dTime)
{
    struct inSphStep in;
    
    if (!msrDoGas(msr)) return;

    in.dCosmoFac = csmTime2Exp(msr->param.csm,dTime);
    in.dEtaCourant = msrEtaCourant(msr);
    in.dEtauDot = msr->param.dEtauDot;
    in.bViscosityLimitdt = msr->param.bViscosityLimitdt;
    pstSphStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrSphViscosityLimiter(MSR msr, double dTime)
{
    struct inSphViscosityLimiter in;
    
    in.bOn = msr->param.bViscosityLimiter;
    in.bShockTracker = msr->param.bShockTracker;

    pstSphViscosityLimiter(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#endif /* GASOLINE */

int msrDumpFrameInit(MSR msr, double dTime, double dStep, int bRestart) {
	/*LCL *plcl = &msr->lcl; -- not used: DCR 12/19/02*/
	char achFile[160];
	
	if (msr->param.dDumpFrameStep > 0 || msr->param.dDumpFrameTime > 0) {
		msr->bDumpFrame = 1;
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		achFile[0] = '\0';
		snprintf(achFile,160,"%s/%s.director",msr->param.achDataSubPath,
				 msr->param.achOutName);
		
		dfInitialize( &msr->df, msr->param.dSecUnit/SECONDSPERYEAR, 
					 dTime, msr->param.dDumpFrameTime, dStep, 
					 msr->param.dDumpFrameStep, msr->param.dDelta, 
					 msr->param.iMaxRung, msr->param.bVDetails,
					 achFile );

		/* Read in photogenic particle list */
		if (msr->df->bGetPhotogenic) {
		  achFile[0] = 0;
		  snprintf(achFile,160,"%s/%s.photogenic",msr->param.achDataSubPath,
				   msr->param.achOutName);
		  msrSetTypeFromFile( msr, achFile, TYPE_PHOTOGENIC );
		}

		if(!bRestart)
			msrDumpFrame( msr, dTime, dStep );
                return 1;
		} else { return 0; }
	}

void msrDumpFrame(MSR msr, double dTime, double dStep)
{
	double sec,dsec1,dsec2,dExp;

	sec = msrTime(msr);

	if (msr->df->iDimension == DF_3D) {
#ifdef VOXEL
		/* 3D Voxel Projection */
		struct inDumpVoxel in;
		assert(0);

		dfSetupVoxel( msr->df, dTime, dStep, &in );

		pstDumpVoxel(msr->pst, &in, sizeof(struct inDumpVoxel), NULL, NULL );
		dsec1 = msrTime(msr) - sec;
		
		dfFinishVoxel( msr->df, dTime, dStep, &in );
		
		dsec2 = msrTime(msr) - sec;
		
		printf("DF Dumped Voxel %i at %g (Wallclock: Render %f tot %f secs).\n",
			   msr->df->nFrame-1,dTime,dsec1,dsec2);
		LOGTIMINGUPDATE( dsec2, TIMING_DumpFrame );
#endif
		}
	else {
		/* 2D Projection */
	    struct inDumpFrame in;
		void *Image; 
		int nImage;
		double com[12];

		if (msr->df->bGetCentreOfMass) {
		  pstCOM(msr->pst, NULL, 0, &com[0], NULL);
		  }

		if (msr->df->bGetPhotogenic) {
		  int type = TYPE_PHOTOGENIC;
		  pstCOMByType(msr->pst, &type, sizeof(int), &com[0], NULL);
		  }

		if (msr->df->bGetOldestStar) {
		  pstOldestStar(msr->pst, NULL, 0, &com[0], NULL);
		  }

		dExp = csmTime2Exp(msr->param.csm,dTime);
		dfSetupFrame( msr->df, dTime, dStep, dExp, &com[0], &in );

		Image = dfAllocateImage( in.nxPix, in.nyPix );
		
		pstDumpFrame(msr->pst, &in, sizeof(struct inDumpFrame), Image, &nImage );
		dsec1 = msrTime(msr) - sec;
		
		dfFinishFrame( msr->df, dTime, dStep, &in, Image );
		
		dsec2 = msrTime(msr) - sec;
		
		printf("DF Dumped Image %i at %g (Wallclock: Render %f tot %f secs).\n",
			   msr->df->nFrame-1,dTime,dsec1,dsec2);
		LOGTIMINGUPDATE( dsec2, TIMING_DumpFrame );

		dfFreeImage( Image );
		}
	}

void msrFormStars(MSR msr, double dTime, double dDelta)
{
#ifdef STARFORM
    struct inFormStars in;
    struct outFormStars outFS;
    struct inFeedback inFB;
    struct outFeedback outFB;
    double dTotMass = -1.0, dTotMetals = -1.0, dTotFe = -1.0, 
            dTotOx = -1.0, dTotEnergy = -1.0;
    double dTotSNEnergy = 0.0;
    int i;
    int iDum;

	double sec,dsec;

	sec = msrTime(msr);

        msrMassMetalsEnergyCheck(msr, &dTotMass, &dTotMetals, 
            &dTotFe, &dTotOx, &dTotEnergy, "Form Stars");
    if(msr->param.bStarForm){
/*		return;*/
    
        in.dTime = dTime;
        msr->param.stfm->dDeltaT = dDelta;
        in.stfm = *msr->param.stfm;
        
        if (msr->param.bVDetails) printf("Form Stars ... ");

        msrMassMetalsEnergyCheck(msr, &dTotMass, &dTotMetals, 
            &dTotFe, &dTotOx, &dTotEnergy, "Form Stars");
        
        msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_ACTIVE);
          /* Only worth the trouble if you're deleting the gas particles
            * as described in method 2 below.  Does not scale well otherwise.
            * msrDomainDecomp(msr, 0, 1);
            */
        msrBuildTree(msr,1,dTotMass,1);
        pstFormStars(msr->pst, &in, sizeof(in), &outFS, NULL);
        if (msr->param.bVDetails)
                    printf("%d Stars formed with mass %g, %d gas deleted\n",
                               outFS.nFormed, outFS.dMassFormed, outFS.nDeleted);
        /* there are two gas particle deletion criteria:
               
           1) in pstFormStars: gas particles with mass less than
           stfm->dMinGasMass are marked for deletion
               
           2) in DeleteGas (see smoothfcn.c): gas particles with 
           mass less than dMinMassFrac of the average mass of neighbouring
           gas particles are also marked for deletion 
               
           - eh, Feb 7/01*/

        /*
         * Find low mass gas particles and mark them for deletion.
             * For better numerical treatment
        if (msr->param.bVDetails) printf("Delete Gas ...\n");
        msrSmooth(msr, dTime, SMX_DELETE_GAS, 0);
         */
        
        /*
         * Record star formation events XXX - not done.
         * NB.  At the moment each star is a star formation event.
         */
          
        /*
         * Distribute mass, and metals of deleted gas particles.
         */
        if (msr->param.bVDetails) printf("Distribute Deleted Gas ...\n");
        msrActiveType(msr, TYPE_DELETED, TYPE_SMOOTHACTIVE);
        msrSmooth(msr, dTime, SMX_DIST_DELETED_GAS, 1);
        /*
         * adjust particle numbers
         */
        msrAddDelParticles(msr);
        msrMassCheck(msr, dTotMass, "Form stars: after particle adjustment");

            dsec = msrTime(msr) - sec;
            printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);
	    LOGTIMINGUPDATE( dsec, TIMING_StarForm );
        }
    /*
     * Calculate energy of SN for any stars for the next timestep.  This
     * requires looking at past star forming events.  Also calculate
     * mass loss.
     */
     
    if(msr->param.bFeedBack) {
        inFB.dTime = dTime;
        inFB.dDelta = dDelta;
        inFB.fb  = *msr->param.fb;
        inFB.sn  = *msr->param.sn;
		if (msr->param.bVDetails) printf("Calculate Feedback ...\n");
		sec = msrTime(msr);
		pstFeedback(msr->pst, &inFB, sizeof(inFB),
					&outFB, &iDum);
		if(msr->param.bVDetails) {
			printf("Feedback totals: mass, energy, metalicity\n");
			for(i = 0; i < FB_NFEEDBACKS; i++){
				printf("feedback %d: %g %g %g\n", i,
					   outFB.fbTotals[i].dMassLoss,
					   outFB.fbTotals[i].dEnergy,
					   outFB.fbTotals[i].dMassLoss != 0.0 ?
					   outFB.fbTotals[i].dMetals
					   /outFB.fbTotals[i].dMassLoss : 0.0);
                                dTotMetals += outFB.fbTotals[i].dMetals;
                                dTotFe += outFB.fbTotals[i].dMIron;
                                dTotOx += outFB.fbTotals[i].dMOxygen;
                                dTotSNEnergy += outFB.fbTotals[i].dEnergy;
                                }
			}


		/*
		 * spread mass lost from SN, (along with energy and metals)
		 * to neighboring gas particles.
		 */
		if (msr->param.bVDetails) printf("Distribute SN Energy ...\n");
		msrActiveType(msr, TYPE_GAS, TYPE_ACTIVE|TYPE_TREEACTIVE);
		msrBuildTree(msr,1,-1.0,1);

		msrResetType(msr, TYPE_STAR, TYPE_SMOOTHDONE);
		msrActiveType(msr, TYPE_STAR, TYPE_SMOOTHACTIVE);
		assert(msr->nSmoothActive == msr->nStar);
		msrSmooth(msr, dTime, SMX_DIST_SN_ENERGY, 1);
		msrMassMetalsEnergyCheck(msr, &dTotMass, &dTotMetals, &dTotFe, 
                    &dTotOx, &dTotSNEnergy, "Form stars: after feedback");

		dsec = msrTime(msr) - sec;
		printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
		LOGTIMINGUPDATE( dsec, TIMING_Feedback );
		}

#endif
    }

void msrSimpleStarForm(MSR msr, double dTime, double dDelta)
{
/* Note: Must be called with an SPH tree built and available */
#ifdef SIMPLESF
    struct inSimpleStarForm in;
    struct outSimpleStarForm out;
	double a,d1,d2;

    double dMass = -1.0;
	double sec,sec1,dsec;

    if(msr->param.bStarForm == 0) return;
    
	sec = msrTime(msr);

    a = csmTime2Exp(msr->param.csm,dTime);

	/* Convert input parameters to code units */
    in.dRateCoeff = msr->param.SSF_dEfficiency*sqrt(4.*M_PI/pow(a,3)); /* G=1 */
	in.dTMax = msr->param.SSF_dTMax;
	d1 = msr->param.SSF_dComovingDenMin;
	d2 = msr->param.SSF_dPhysDenMin/msr->param.dGmPerCcUnit*pow(a,3);
	in.dDenMin = (d1>d2 ? d1 : d2);
	in.dDelta = dDelta;

	in.dTime = dTime;
	in.dInitStarMass = msr->param.SSF_dInitStarMass;
	in.dESNPerStarMass = msr->param.SSF_dESNPerStarMass/msr->param.dErgPerGmUnit;
#define SECONDSPERYEAR   31557600.
	in.dtCoolingShutoff = msr->param.SSF_dtCoolingShutoff*SECONDSPERYEAR/msr->param.dSecUnit;
	in.bdivv = msr->param.SSF_bdivv;
    
    if (msr->param.bVDetails) printf("Simple Star Form ... ");

    dMass = msrMassCheck(msr, -2.0, "Form Stars");
	msrActiveType(msr,0,TYPE_SMOOTHACTIVE|TYPE_SMOOTHDONE|TYPE_NbrOfACTIVE);
	/* New stars will be set to TYPE_SMOOTHACTIVE when created */
	
    pstSimpleStarForm(msr->pst, &in, sizeof(in), &out, NULL);
    if (msr->param.bVDetails)
		printf("%d Stars formed with mass %g, %d gas deleted\n",
			   out.nFormed, out.dMassFormed, out.nDeleted);

    /*
     * adjust particle numbers
     */
    msrAddDelParticles(msr);
    msrMassCheck(msr, dMass, "Form stars: after particle adjustment");

	sec1 = msrTime(msr);
	dsec = sec1 - sec;
	printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec)
	LOGTIMINGUPDATE( dsec, TIMING_StarForm );

	if (msr->param.bFeedBack && out.nFormed) {
		/* Build a tree to distribute energy from SN, if nFormed > 0  */

		/* Any new stars have been set to SMOOTHACTIVE */
		msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE|TYPE_ACTIVE);
		msrBuildTree(msr,1,dMass,1);
		msrSmooth(msr, dTime, SMX_SIMPLESF_FEEDBACK, 1);
		dsec = msrTime(msr) - sec1;
		printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
		LOGTIMINGUPDATE( dsec, TIMING_Feedback );
		}
#endif
    }

#ifdef GLASS

void msrInitGlass(MSR msr)
{
    struct inRandomVelocities in;
    
    in.dMaxVelocityL = msr->param.dGlassVL; 
    in.dMaxVelocityR = msr->param.dGlassVR; 
    pstRandomVelocities(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#endif

#ifdef COLLISIONS

/*DEBUG IMPORTANT: many of the following routines use smooth and
  resmooth but bypass msrSmooth() and msrReSmooth() -- this means that
  changes to the input parameters for smooth/resmooth MUST be
  reflected here (cf. inSmooth in pst.h, and SMF in smoothfcn.h),
  otherwise there is the risk of passing uninitialized values with
  unpredictable consequences that are hard to debug.  Ideally we
  should use msrSmooth() and msrReSmooth(), so that everything is kept
  in the same place.  Or at least an initialization function, e.g.,
  msrInitinSmooth(), should be used to make sure everything is
  initialized to values that will trip an assert in smooth().
  (Similar arguments can be made for how external potentials are
  handled in msrGravity().)*/

int
msrNumRejects(MSR msr)
{
	struct outNumRejects out;

	pstNumRejects(msr->pst,NULL,0,&out,NULL);
	return out.nRej;
	}

void msrFindRejects(MSR msr)
{
	/*
	** Checks initial conditions for particles with overlapping
	** physical or Hill spheres.  The latter case only makes sense for
	** particles orbiting a massive central body, like the Sun, and is
	** controlled by the value of msr->param.dCentMass.  Rejects are
	** written to REJECTS_FILE (see ssdefs.h).  This procedure is
	** intended to be called iteratively from an external
	** initial-conditions program, such as ssic or patchic.  Note that
	** the time is assumed to be zero, and in the case of periodic
	** boundary conditions, the particles are assumed to be in the
	** bounding box.
	*/
	
	/*
	** IMPORTANT NOTE: the number of rejects detected by this
	** algorithm may differ when running in parallel, even when the
	** same number of processes are used (in multiple instances).
	** This is because particle (and/or Hill sphere) overlaps are only
	** considered in pairwise fashion here, but there could be chains
	** of particles that are overlapping.  In the latter case, the
	** number of particles rejected will depend on which overlaps are
	** detected first, and in parallel that order depends on the
	** number of processes and other factors.  This is NOT an error --
	** it just means that some particles will be flagged unnecessarily
	** as rejects.  This results in some inefficiency, but a better
	** approach that still relies on pairwise detection is not
	** immediately obvious.
	*/
	
	struct inSmooth smooth;
	int nRej = 0;
	
	if (msr->param.bVStart)
		puts("Checking for rejected ICs...");
	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
	msrDomainDecomp(msr,0,1);
	msrActiveType(msr,TYPE_ALL,TYPE_SMOOTHACTIVE|TYPE_TREEACTIVE);
	msrBuildTree(msr,0,-1.0,1); /* 1=density tree */
	/* can't use msrSmooth() because of special cases (compile options) */
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 1; /* use combiner cache */
	smooth.iSmoothType = SMX_REJECTS;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
	smooth.smf.dTime = 0.0; /*DEBUG assumes start time is zero...*/
#ifdef SLIDING_PATCH
	smooth.smf.PP = msr->param.PP; /* struct copy */
#endif
	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
	nRej = msrNumRejects(msr);
	if (nRej > 0) {
		printf("%i reject%s found!\n",nRej,(nRej == 1 ? "":"s"));
		msrReorder(msr);
		msrOutArray(msr,REJECTS_FILE,OUT_REJECTS_ARRAY);
		_msrExit(msr,1);
		}
	else {
		puts("No rejects found.");
		_msrExit(msr,0);
		}
	assert(0); /* unreachable statement */
	}

void msrOneNodeReadSS(MSR msr,struct inReadSS *in)
{
    int i,id;
    int *nParts;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    pstOneNodeReadInit(msr->pst,in,sizeof(*in),nParts,&nid);
    assert(nid == (int) (msr->nThreads*sizeof(*nParts)));
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
#ifdef SPRINGS
	plcl->pkd->bReadSpringsData = in->bReadSpringsData;
#endif
#ifdef DEM
	plcl->pkd->bReadDEMData = in->bReadDEMData;
#endif
#ifdef DEM_TRACK_ORIENT
	plcl->pkd->bReadOrient = in->bReadOrient;
#endif
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadSS(plcl->pkd,achInFile,nStart,nParts[id]);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadSS(plcl->pkd,achInFile,0,nParts[0]);
    }

double msrReadSS(MSR msr)
{
	SSIO ssio;
	SSHEAD head;
	struct inReadSS in;
	struct inSetParticleTypes intype;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double dTime;

	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,in.achInFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achInFile,achInFile);

		if (ssioOpen(achInFile,&ssio,SSIO_READ)) {
			fprintf(stderr,"Could not open InFile:%s\n",achInFile);
			_msrExit(msr,1);
			}
		}
	else {
		fprintf(stderr,"No input file specified\n");
		_msrExit(msr,1);
		}

	/* Read header */

	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Could not read header of InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		fprintf(stderr,"Could not close InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}

	msr->N = msr->nDark = head.n_data;
	msr->nGas = msr->nStar = 0;
	msr->nMaxOrder = msr->N - 1;
	msr->nMaxOrderGas = msr->nGas - 1; /* always -1 */
	msr->nMaxOrderDark = msr->nDark - 1;

	dTime = head.time;
	if (msr->param.bVStart) {
		double tTo;
		printf("Input file...N=%i,Time=%g\n",msr->N,dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		printf("Simulation to Time:%g\n",tTo);
		}

#ifdef SLIDING_PATCH
	if (dTime != 0.0)
		fprintf(stderr,"WARNING: SLIDING_PATCH assumes start time = 0!\n"); /* only a problem if you don't have an appropriate ss file to start from */
#endif

	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;	/* always zero */
	in.nStar = msr->nStar;	/* always zero */
	in.iOrder = msr->param.iOrder;
	/*
	 ** Since pstReadSS causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;

	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;

#ifdef SPRINGS
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.bReadSpringsData = msr->param.FO.SP.bReadSpringsData;
	if (in.bReadSpringsData) {
		/* check header for consistency */
		char achSpringsFilename[PST_FILENAME_SIZE];
		XDR xdrs;
		FILE *fp;
		double dDum;
		int iDum;
		snprintf(achSpringsFilename,PST_FILENAME_SIZE,"%s.spr",achInFile);
		fp = fopen(achSpringsFilename,"r");
		if (fp == NULL) {
			fprintf(stderr,"ERROR: could not open springs file \"%s\".\n",achSpringsFilename);
			_msrExit(msr,1);
			}
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		xdr_double(&xdrs,&dDum);
		if (dDum != head.time)
			fprintf(stderr,"WARNING: time stamps in springs file (%g) and ss file (%g) do not match exactly.\n",dDum,head.time);
		xdr_int(&xdrs,&iDum);
		if (iDum != head.n_data) {
			fprintf(stderr,"ERROR: number of particles in springs file (%i) and ss file (%i) do not match.\n",iDum,head.n_data);
			_msrExit(msr,1);
			}
		xdr_int(&xdrs,&iDum);
		if (iDum != MAX_NUM_SPRINGS_PER_PARTICLE) {
			fprintf(stderr,"ERROR: MAX_NUM_SPRINGS_PER_PARTICLE in springs file (%i) does not match compiled value (%i).\n",iDum,MAX_NUM_SPRINGS_PER_PARTICLE);
			_msrExit(msr,1);
			}
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* SPRINGS */

#ifdef DEM
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.bReadDEMData = msr->param.FO.DP.bReadDEMData;
	if (in.bReadDEMData) {
		/* check header for consistency */
		char achDEMFilename[PST_FILENAME_SIZE];
		XDR xdrs;
		FILE *fp;
		double dDum;
		int iDum;
		snprintf(achDEMFilename,PST_FILENAME_SIZE,"%s.dem",achInFile);
		fp = fopen(achDEMFilename,"r");
		if (fp == NULL) {
			fprintf(stderr,"ERROR: could not open DEM file \"%s\".\n",achDEMFilename);
			_msrExit(msr,1);
			}
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		xdr_double(&xdrs,&dDum);
		if (dDum != head.time)
			fprintf(stderr,"WARNING: time stamps in DEM file (%g) and ss file (%g) do not match exactly.\n",dDum,head.time);
		xdr_int(&xdrs,&iDum);
		if (iDum != head.n_data) {
			fprintf(stderr,"ERROR: number of particles in DEM file (%i) and ss file (%i) do not match.\n",iDum,head.n_data);
			_msrExit(msr,1);
			}
		xdr_int(&xdrs,&iDum);
		if (iDum != MAX_NUM_OVERLAPS_PER_PARTICLE) {
			fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE in DEM file (%i) does not match compiled value (%i).\n",iDum,MAX_NUM_OVERLAPS_PER_PARTICLE);
			_msrExit(msr,1);
			}
		xdr_int(&xdrs,&iDum);
#ifdef WALLS
		if (iDum != MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS) {
			fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS in DEM file (%i) does not match compiled value (%i).\n",iDum,MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS);
			_msrExit(msr,1);
			}
#else
		if (iDum != -1) {
			fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS in DEM file (%i) does not match expected value (-1).\n",iDum);
			_msrExit(msr,1);
			}
#endif /* WALLS */
		xdr_int(&xdrs,&iDum);
#ifdef DEM_ROTATION_DASHPOT
		if (iDum != 1) {
			fprintf(stderr,"ERROR: Rotation dashpot flag in DEM file not set while code compiled with DEM_ROTATION_DASHPOT.\n");
			_msrExit(msr,1);
			}
#else
		if (iDum != 0) {
			fprintf(stderr,"ERROR: Rotation dashpot flag in DEM file does not match expected value (0).\n");
			_msrExit(msr,1);
			}
#endif /* DEM_ROTATION_DASHPOT */
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* DEM */

#ifdef DEM_TRACK_ORIENT
	in.bReadOrient = msr->param.FO.DP.bReadOrient;
	if (in.bReadOrient) {
		ORIHEAD orihead;
		char achOriFilename[PST_FILENAME_SIZE];
		snprintf(achOriFilename,PST_FILENAME_SIZE,"%s.ori",achInFile);
		if (ssioOpen(achOriFilename,&ssio,SSIO_READ)) {
			fprintf(stderr,"Could not open %s\n",achOriFilename);
			_msrExit(msr,1);
			}
		if (ssioHead(&ssio,&orihead)) {
			fprintf(stderr,"Could not read header of %s\n",achOriFilename);
			_msrExit(msr,1);
			}
		if (orihead.time != head.time)
			fprintf(stderr,"WARNING: time stamps in .ori file (%g) and ss file (%g) do not match exactly.\n",orihead.time,head.time);
		if (orihead.n_data != head.n_data) {
			fprintf(stderr,"ERROR: number of particles in .ori file (%i) and ss file (%i) do not match.\n",orihead.n_data,head.n_data);
			_msrExit(msr,1);
			}
		if (orihead.iMagicNumber != SSIO_MAGIC_STANDARD) {
			fprintf(stderr,"ERROR: only standard format for orientation file supported on read.\n");
			_msrExit(msr,1);
			}
		if (ssioClose(&ssio)) {
			fprintf(stderr,"Could not close %s\n",achOriFilename);
			_msrExit(msr,1);
			}
		}
#endif /* DEM_TRACK_ORIENT */

	if (msr->param.bParaRead)
	    pstReadSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadSS(msr,&in);
	if (msr->param.bVDetails) puts("Input file successfully read.");
	/*
	 ** Set particle ACTIVE flags to correspond to appropriate type.
	 */
	intype.nSuperCool = msr->param.nSuperCool;
	assert(intype.nSuperCool == 0); /* better be zero... */
	pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
	/*
	 ** Now read in the output points, passing the initial time.
	 ** We do this only if nSteps is not equal to zero.
	 */
	if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}

void msrOneNodeWriteSS(MSR msr,struct inWriteSS *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

#ifdef SPRINGS
	plcl->pkd->bNoSpringsOutput = in->bNoSpringsOutput;
#endif
#ifdef DEM
	plcl->pkd->nDEMOutputs = in->nDEMOutputs;
#endif

    /* 
     * First write our own particles.
     */
    pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart,in->bReduced);
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd,id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteSS(plcl->pkd,achOutFile,nStart,in->bReduced);
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }

void msrWriteSSHead(MSR msr,char *achOutFile,double dTime,int bReduced)
{
	SSIO ssio;
	SSHEAD head;

	if (ssioOpen(achOutFile,&ssio,SSIO_WRITE)) {
            printf("Could not open OutFile:%s\n",achOutFile);
            _msrExit(msr,1);
            }

	/* Write header */

	head.time = dTime;
	head.n_data = msr->N;
	head.iMagicNumber = (bReduced ? SSIO_MAGIC_REDUCED : SSIO_MAGIC_STANDARD);

	if (ssioHead(&ssio,&head)) {
		printf("Could not write header of OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		printf("Could not close OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}

#ifdef SPRINGS
	/* write springs data header to a separate file */
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	if (!bReduced && !msr->param.FO.SP.bNoSpringsOutput) {
		char achSpringsFilename[PST_FILENAME_SIZE];
		XDR xdrs;
		FILE *fp;
		snprintf(achSpringsFilename,PST_FILENAME_SIZE,"%s.spr",achOutFile);
		fp = fopen(achSpringsFilename,"w");
		assert(fp != NULL);
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdr_double(&xdrs,&head.time);
		xdr_int(&xdrs,&head.n_data);
		int iDum = MAX_NUM_SPRINGS_PER_PARTICLE;
		xdr_int(&xdrs,&iDum);
		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* SPRINGS */

#ifdef DEM
	/* write DEM data header to a separate file */
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	if (!bReduced && msr->param.FO.DP.nDEMOutputs != 0) {
		char achDEMFilename[PST_FILENAME_SIZE];
		XDR xdrs;
		FILE *fp;
		snprintf(achDEMFilename,PST_FILENAME_SIZE,"%s.dem",achOutFile);
		fp = fopen(achDEMFilename,"w");
		assert(fp != NULL);
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		xdr_double(&xdrs,&head.time);
		xdr_int(&xdrs,&head.n_data);
		int iDum = MAX_NUM_OVERLAPS_PER_PARTICLE;
		xdr_int(&xdrs,&iDum);
		
#ifdef WALLS
		iDum = MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS;
#else
		iDum = -1;
#endif /* WALLS */
		xdr_int(&xdrs,&iDum);

#ifdef DEM_ROTATION_DASHPOT
		iDum = 1;
#else
		iDum = 0;
#endif /* DEM_ROTATION_DASHPOT */
		xdr_int(&xdrs,&iDum);

		xdr_destroy(&xdrs);
		fclose(fp);
		}
#endif /* DEM */

#ifdef DEM_TRACK_ORIENT
	{
		ORIHEAD ohead;
		char achOriFilename[PST_FILENAME_SIZE];
		snprintf(achOriFilename,PST_FILENAME_SIZE,"%s.ori",achOutFile);
		if (ssioOpen(achOriFilename,&ssio,SSIO_WRITE)) {
            printf("Could not open \"%s\" for writing.\n",achOriFilename);
            _msrExit(msr,1);
            }
		ohead.time = dTime;
		ohead.n_data = msr->N;
		ohead.iMagicNumber = (bReduced ? SSIO_MAGIC_REDUCED : SSIO_MAGIC_STANDARD);
		if (ssioOriHead(&ssio,&ohead)) {
			printf("Could not write header for \"%s\".\n",achOriFilename);
			_msrExit(msr,1);
			}
		if (ssioClose(&ssio)) {
			printf("Could not close \"%s\".\n",achOriFilename);
			_msrExit(msr,1);
			}
		}
#endif /* DEM_TRACK_ORIENT */
    }

void msrWriteSS(MSR msr,char *pszFileName,double dTime,int bReduced)
{
	struct inWriteSS in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;

#ifdef SEISMIC_WAVE //outputing data
	WALL_PARAMS *WP = &msr->param.CP.WP;
	WALL *w;
	WALL_DATA *wd;
	int i;
	FILE *fp;
	fp = fopen("SEISMICWAVESviz.out","a"); //output file
	assert(fp);
	 
	for (i=0;i<WP->nWalls;i++) {
		w = &WP->pWalls[i];
		wd = &w->wd;
		fprintf(fp,"%.16e %.16e %.16e\n",wd->vOrigin[0],wd->vOrigin[1],wd->vOrigin[2]);
		}
	fclose(fp);
#endif /* SEISMIC_WAVE */

#ifdef DEM_WALLS_REACT //outputing data
	WALL_PARAMS *WP = &msr->param.CP.WP;
	WALL *w;
	WALL_DATA *wd;
	int i, iAssemblySwitch[MAX_NUM_WALL_ASSEMBLIES];
	FILE *fp;
	static int iRedCounter = 0; //loopcounter.. static int, can't be re-initialized
  if (msr->param.iStartStep!=0 && iRedCounter==0){ //don't reprint end step stats on restart
    //do nothing
    }
  else {
    fp = fopen("WALLSREACTviz.out","a"); //This file will have extraneous detail for visualization (origin, vertex1, vertex 2) for each wall
    assert(fp);
    fprintf(fp,"%d",(iRedCounter*msr->param.iRedOutInterval)+msr->param.iStartStep); //loopcounter multiplied by spacing interval
    for (i=0;i<WP->nWalls;i++) {
      w = &WP->pWalls[i];
      wd = &w->wd;
      //if (wd->dMass != 0.0) {
        if (wd->iType == WallTriangle || wd->iType == WallRectangle)
          fprintf(fp," 0 %d %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e",w->iWallID,wd->vOrigin[0],wd->vOrigin[1],wd->vOrigin[2],wd->vVertex1[0],wd->vVertex1[1],wd->vVertex1[2],wd->vVertex2[0],wd->vVertex2[1],wd->vVertex2[2]);
        if ((wd->iType==WallCylinderFinite) || (wd->iType==WallDisk) || (wd->iType==WallShell))
          fprintf(fp," 1 %d %.16e %.16e %.16e %.16e %.16e %.16e",w->iWallID,wd->vOrigin[0],wd->vOrigin[1],wd->vOrigin[2],wd->vOrient[0],wd->vOrient[1],wd->vOrient[2]);
        //}
      }
    fprintf(fp,"\n"); //new line
    fclose(fp);
    
    for (i=0;i<MAX_NUM_WALL_ASSEMBLIES;i++)
      iAssemblySwitch[i]=1;
    fp = fopen("WALLSREACTdata.out","a"); // CoG, Velocities, Rotation Rates, Principal Axes
    assert(fp);
    fprintf(fp,"%d",(iRedCounter*msr->param.iRedOutInterval)+msr->param.iStartStep); //loopcounter multiplied by spacing interval
    for (i=0;i<WP->nWalls;i++) {
      w = &WP->pWalls[i];
      wd = &w->wd;
      if (wd->dMass != 0.0){
        if (iAssemblySwitch[wd->iAssembly] == 1){
          fprintf(fp," %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e",wd->vCoG[0],wd->vCoG[1],wd->vCoG[2],wd->vVel[0],wd->vVel[1],wd->vVel[2],wd->vRot[0],wd->vRot[1],wd->vRot[2],wd->vPxj[0],wd->vPxj[1],wd->vPxj[2],wd->vPyj[0],wd->vPyj[1],wd->vPyj[2],wd->vPzj[0],wd->vPzj[1],wd->vPzj[2]);
          iAssemblySwitch[wd->iAssembly]=0;
          }
        //break; /*Only need the values from the first inertial Wall, so terminate loop*/
        }
      }
    fprintf(fp,"\n"); //new line
    fclose(fp);
    }
  iRedCounter++;
#endif /* DEM_WALLS_REACT */

	/*
	 ** Calculate where each processor should start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	_msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

	msrWriteSSHead(msr,achOutFile,dTime,bReduced);

	in.bReduced = bReduced;
#ifdef SPRINGS
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.bNoSpringsOutput = msr->param.FO.SP.bNoSpringsOutput;
#endif /* SPRINGS */
#ifdef DEM
	assert(msr->param.FO.iForceOverrideOption == FO_STRENGTH);
	in.nDEMOutputs = msr->param.FO.DP.nDEMOutputs;
#endif /* DEM */
	
	if (msr->param.bParaWrite)
	        pstWriteSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
		msrOneNodeWriteSS(msr,&in);

	if (msr->param.bVDetails) puts("Output successfully written.");
	}

static const char *_msrParticleLabel(MSR msr,int iColor)
{
	/* For use with msrDoCollisions() only */

#ifdef WALLS
	if (iColor < 0) {
		static char ach[32];
		snprintf(ach,32,"WALL %i",-1 - iColor);
		return ach;
		}
#endif

#ifdef AGGS
	return "BODY"; /* this entire function needs to be reworked... */
#endif

	switch (iColor) {
	case SUN:
		return "SUN";
	case JUPITER:
		return "JUPITER";
	case SATURN:
		return "SATURN";
	case URANUS:
		return "URANUS";
	case NEPTUNE:
		return "NEPTUNE";
	case PLANETESIMAL:
		return "PLANETESIMAL";
	default:
		return "UNKNOWN";
		}
	}

void msrDoCollLog(MSR msr,COLLIDER *c1,COLLIDER *c2,struct outDoCollision *outDo,int option,double dt,double dTime) 
{
	FILE *fp;
	int i;
	XDR xdrs;
	double dDum,dtOut;
	COLLIDER *c;

	if ((dtOut = dt) < 0.0) {
		/*
		** Negative timesteps can mess up, for example, merger history
		** tracking (cf. ssg in ss_core), so we don't allow them in
		** the output.
		*/
		static int bFirstWarning = 1;
		if (bFirstWarning) {
			if (msr->param.bVWarnings)
				fprintf(stderr,"WARNING: Negative time to collision (%g) in collision log set to zero\n",dt);
			bFirstWarning = 0;
			}
		dtOut = 0.0;
		}

	if (option == COLL_LOG_VERBOSE) {
		fp = fopen(msr->param.achCollLog,"a");
		assert(fp != NULL);
#ifdef AGGS
		for (i=0;i<3;i++) {
			if (!COLLIDER_IS_AGG(c1))
				c1->r[i] += c1->v[i]*dt;
			if (!COLLIDER_IS_AGG(c2))
				c2->r[i] += c2->v[i]*dt;
			}
#else
#ifdef RORY_EXTENSION
		if (dt != -1.0) /* -1.0 => binary merge: we don't want to update pos. */
#endif
		for (i=0;i<3;i++) {
			c1->r[i] += c1->v[i]*dt;
			c2->r[i] += c2->v[i]*dt;
			}
#endif
		fprintf(fp,"%s-%s COLLISION:t=%e\n",
				_msrParticleLabel(msr,c1->iColor),
				_msrParticleLabel(msr,c2->iColor),dTime + dtOut);
#ifdef AGGS
		if (COLLIDER_IS_AGG(c1))
			fprintf(fp,"***1:p=%i,o=%i,i=%i,oi=%i,M=%e,R=*%e*,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c1->id.iPid,c1->id.iOrder,c1->id.iIndex,c1->id.iOrgIdx,
					c1->agg.mass,c1->fRadius,c1->dt,c1->iRung,
					c1->agg.r_com[0],c1->agg.r_com[1],c1->agg.r_com[2],
					c1->agg.v_com[0],c1->agg.v_com[1],c1->agg.v_com[2],
					c1->agg.omega[0],c1->agg.omega[1],c1->agg.omega[2]);
		else
#endif
			fprintf(fp,"***1:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c1->id.iPid,c1->id.iOrder,c1->id.iIndex,c1->id.iOrgIdx,
					c1->fMass,c1->fRadius,c1->dt,c1->iRung,
					c1->r[0],c1->r[1],c1->r[2],
					c1->v[0],c1->v[1],c1->v[2],
					c1->w[0],c1->w[1],c1->w[2]);
#ifdef AGGS
		if (COLLIDER_IS_AGG(c2))
			fprintf(fp,"***2:p=%i,o=%i,i=%i,oi=%i,M=%e,R=*%e*,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c2->id.iPid,c2->id.iOrder,c2->id.iIndex,c2->id.iOrgIdx,
					c2->agg.mass,c2->fRadius,c2->dt,c2->iRung,
					c2->agg.r_com[0],c2->agg.r_com[1],c2->agg.r_com[2],
					c2->agg.v_com[0],c2->agg.v_com[1],c2->agg.v_com[2],
					c2->agg.omega[0],c2->agg.omega[1],c2->agg.omega[2]);
		else
#endif
			fprintf(fp,"***2:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
					"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
					c2->id.iPid,c2->id.iOrder,c2->id.iIndex,c2->id.iOrgIdx,
					c2->fMass,c2->fRadius,c2->dt,c2->iRung,
					c2->r[0],c2->r[1],c2->r[2],
					c2->v[0],c2->v[1],c2->v[2],
					c2->w[0],c2->w[1],c2->w[2]);
		fprintf(fp,"***OUTCOME=%s dT=%e\n",
				outDo->iOutcome == MISS ? "MISS" :
				outDo->iOutcome == MERGE ? "MERGE" :
				outDo->iOutcome == BOUNCE ? "BOUNCE" :
#ifdef RORY_EXTENSION
				outDo->iOutcome == BINARY_MERGE ? "BINARY MERGE" :
#endif
				outDo->iOutcome == FRAG ? "FRAG" :
#ifdef WALLS
				outDo->iOutcome == DEATH ? "DEATH" :
#endif
				"UNKNOWN",outDo->dT);
		for (i=0;i<(outDo->nOut < MAX_NUM_FRAG ? outDo->nOut : MAX_NUM_FRAG);i++) {
			c = &outDo->Out[i];
#ifdef AGGS
			if (COLLIDER_IS_AGG(c))
				fprintf(fp,"***out%i:p=%i,o=%i,i=%i,oi=%i,M=%e,R=*%e*,rung=%i,"
						"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",i,
						c->id.iPid,c->id.iOrder,c->id.iIndex,c->id.iOrgIdx,
						c->agg.mass,c->fRadius,c->iRung,
						c->agg.r_com[0],c->agg.r_com[1],c->agg.r_com[2],
						c->agg.v_com[0],c->agg.v_com[1],c->agg.v_com[2],
						c->agg.omega[0],c->agg.omega[1],c->agg.omega[2]);
			else
#endif
				fprintf(fp,"***out%i:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,rung=%i,"
						"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",i,
						c->id.iPid,c->id.iOrder,c->id.iIndex,c->id.iOrgIdx,
						c->fMass,c->fRadius,c->iRung,
						c->r[0],c->r[1],c->r[2],
						c->v[0],c->v[1],c->v[2],
						c->w[0],c->w[1],c->w[2]);
			}
		fclose(fp);
		}
	else if (option==COLL_LOG_TERSE) {
		/*
		** FORMAT: For each event, time (double), collider 1 iOrgIdx
		** (int), collider 2 iOrgIdx (int), number of post-collision
		** particles (int), iOrgIdx for each of these (n * int).
		*/
		if (outDo->iOutcome != MERGE &&
#ifdef RORY_EXTENSION
			outDo->iOutcome != BINARY_MERGE &&
#endif
			outDo->iOutcome != FRAG)
			return; /* only care when particle indices change */
		fp = fopen(msr->param.achCollLog,"a");
		assert(fp != NULL);
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		dDum = dTime + dtOut;
		(void) xdr_double(&xdrs,&dDum);
		(void) xdr_int(&xdrs,&c1->id.iOrgIdx);
		(void) xdr_int(&xdrs,&c2->id.iOrgIdx);
		(void) xdr_int(&xdrs,&outDo->nOut);
		for (i=0;i<(outDo->nOut < MAX_NUM_FRAG ? outDo->nOut : MAX_NUM_FRAG);i++)
			(void) xdr_int(&xdrs,&outDo->Out[i].id.iOrgIdx);
		xdr_destroy(&xdrs);
		(void) fclose(fp);
		}
	else
		assert(0); /* invalid collision log option */
	}

#ifdef RORY_EXTENSION

void msrCheckForBinary(MSR msr,double dTime) 
{
  /* 
   ** Determine if any particles in planetary simulations are binaries. We
   ** only are concerned with particles on a rung higher the iMinBinaryRung.
   */

  struct inSmooth smooth;
  struct outSmooth outSm;
  struct inGetColliderInfo inGet;
  struct outGetColliderInfo outGet;
  struct inBinary inBnry;
  struct outBinary outBnry;
  struct outDoCollision outDo;
  struct inMrgBnry inMrg;
  struct outMrgBnry outMrg;
  double dt;

  msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
  /* Get all particles with iRung>iMinBinaryRung. */
  msrActiveRung(msr,msr->param.iMinBinaryRung,1);
  
  smooth.nSmooth = 2;
  smooth.bPeriodic = msr->param.bPeriodic;
  smooth.bSymmetric = 0;
  smooth.iSmoothType = SMX_FINDBINARY;
  smooth.dfBall2OverSoft2 = 0.0; /* No softening limit */
  smooth.smf.dMaxBinaryEcc = msr->param.dMaxBinaryEcc;
  smooth.smf.dTime = dTime;
#ifdef SLIDING_PATCH
  smooth.smf.PP = msr->param.PP; /* struct copy */
#endif

  /* 
   ** Determine which particle meet the criterion for merging. We will only 
   ** search for one binary per big timestep (dDelta). Although it is 
   ** possible that multiple binaries could occur in a single timestep, it is
   ** most likely less computationally expensive to only search for one binary.
   ** To maximize speed-up, we will merge the binary on the highest rung. If
   ** multiple binaries appear to be slowing down the code considerably, a
   ** pstNextCollision call should be added.
   */
  pstSmooth(msr->pst,&smooth,sizeof(smooth),&outSm,NULL);
  inBnry.n = outBnry.n = 0;
  pstFindTightestBinary(msr->pst,&inBnry,sizeof(inBnry),&outBnry,NULL);
  if (outBnry.n==2) {
    assert(outBnry.dBindEn < 0);
    inGet.iOrder=outBnry.iOrder1;
    pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
    inMrg.c1 = outGet.Collider;
    inGet.iOrder=outBnry.iOrder2;
    pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
    inMrg.c2 = outGet.Collider;
#ifdef SLIDING_PATCH
    inMrg.PP = msr->param.PP; /* struct copy */
#endif
    inMrg.dTime = dTime;
    inMrg.dDelta=msr->param.dDelta;
    inMrg.iStartStep=msr->param.iStartStep;
    inMrg.bPeriodic = msr->param.bPeriodic;
    inMrg.dDensity=msr->param.CP.dDensity;
    outMrg.n=0;
    pstMergeBinary(msr->pst,&inMrg,sizeof(inMrg),&outMrg,NULL);
    if (msr->param.iCollLogOption) {
      outDo.iOutcome=BINARY_MERGE;
      outDo.nOut=1;
      outDo.Out[0]=outMrg.cOut;
      dt=0;
	  if (msr->param.iCollLogOption != COLL_LOG_NONE)
		  msrDoCollLog(msr,&inMrg.c1,&inMrg.c2,&outDo,msr->param.iCollLogOption,dt,dTime);
    }
    if (msr->param.bVDetails) {
      printf("Binary Merge at time: %e.\n",dTime);
    }
    msrAddDelParticles(msr);
  }
}

#ifdef SLIDING_PATCH

void
msrPickNewCoordinates(PARTICLE *p,int n,double *dHill,double dxPeriod,double dyPeriod,double **new,int pick)
{
    int j,k,ok;
    double hill2,dist2;
    
    *new=malloc(3*sizeof(double));
    (*new)[0]=dxPeriod*(randUniform() - 0.5);
    (*new)[1]=dyPeriod*(randUniform() - 0.5);
    (*new)[2]=p[pick].r[2];
    for (j=0;j<n;j++) {
	if (j!=pick) {
	    hill2=(dHill[pick] + dHill[j])*(dHill[pick] + dHill[j]);
	    ok=0;
	    while (!ok) {
		dist2=0;
		for (k=0;k<3;k++)
		    dist2+=((*new)[k] - p[j].r[k])*((*new)[k] - p[j].r[k]);
		if (dist2 > hill2) ok=1;
		else {
		    (*new)[0]=dxPeriod*(randUniform() - 0.5);
		    (*new)[1]=dyPeriod*(randUniform() - 0.5);
		    }		    
		}
	    }    
	}
    }

void
msrRandomLog(int iStep,PARTICLE p,PARTICLE *pn,double *nsep2,int nn,PARTICLE *pr,double *rsep2,int nr,double *x, double dHillRadius)
{
    FILE *fp;
    int i;
    

    fp=fopen("random.log","a");
    assert(fp);
    fprintf(fp,"Radomizing Large Mass, step: %d\n",iStep);
    fprintf(fp,"Large Mass Info: iOrder=%d, Mass=%e, r=(%e,%e,%e), Hill Radius=%e, nNeighbors=%d\n",p.iOrder,p.fMass,p.r[0],p.r[1],p.r[2],dHillRadius,nn);
    if (nn == 0)
	fprintf(fp,"\tNo Neighbor Particles.\n");
    else {
	for (i=0;i<nn;i++) {
	    fprintf(fp,"\tNeighbor %d: iOrder=%d, r=(%.3e,%.3e,%.3e), v=(%.3e,%.3e,%.3e), Distance=%e\n",i,pn[i].iOrder,pn[i].r[0],pn[i].r[1],pn[i].r[2],pn[i].v[0],pn[i].v[1],pn[i].v[2],sqrt(nsep2[i]));
	    }
	
	}
    fprintf(fp,"New Position: r=(%e,%e,%e), nNeighbors=%d\n",x[0],x[1],x[2],nr);
    if (nr == 0)
	fprintf(fp,"\tNo replacement particles.\n");
    else {	
	for (i=0;i<nr;i++) {
	    fprintf(fp,"\tReplacement %d: iOrder=%d, r=(%.3e,%.3e,%.3e), v=(%.3e,%.3e,%.3e), Distance=%e\n",i,pr[i].iOrder,pr[i].r[0],pr[i].r[1],pr[i].r[2],pr[i].v[0],pr[i].v[1],pr[i].v[2],sqrt(rsep2[i]));
	    }
	}

    fclose(fp);    
    }

int 
msrCheckLargeMassOverlap(PARTICLE *p,double *hill,int n)
{
    /* If two large masses have overlapping spheres, then we will not
       perform the randomization this time. */

    int i,j,k;
    double dist2;
    
    for (i=0;i<n;i++) {
	for (j=i+1;j<n-i;j++) {
	    dist2=0;
	    for (k=0;k<3;k++) 
		dist2+=(p[i].r[k]-p[j].r[k])*(p[i].r[k]-p[j].r[k]);
	    if (dist2 < (hill[i]+hill[j])*(hill[i]+hill[j]) )
		return 0;
	}
    }
    return 1;	    	    
}

int
msrGetNextRandomTime(int iBaseTime,int iTimeNow)
{
    int dev;
    int ok=0;
    
    
    while (!ok) {
	dev=randPoisson(iBaseTime);
	if (dev > 0 && dev < 5*iBaseTime) ok=1;
	}
    
    return dev + iTimeNow;
    }
    
int
msrPickMass(PARTICLE *p,int n)
{
    int i;
    double *prob,pick,msum2=0;
    

    prob=malloc(n*sizeof(double));
    for (i=0;i<n;i++) 
	msum2+= p[i].fMass*p[i].fMass;
    prob[0]=(p[0].fMass*p[0].fMass)/msum2;    
    for (i=1;i<n;i++)
	prob[i]=prob[i-1]+(p[i].fMass*p[i].fMass)/msum2;
    assert(abs(prob[n-1] - 1) < 1e-10);
    pick = randUniform();
    for (i=0;i<n;i++) {
	if (pick < prob[i]) {
	    return i;
	    }
	}
    assert(0); /* Shouldn't get here! */
	return -1; /* to keep compiler happy */
    }

void
msrRandomizeLargeMasses(MSR msr,int iStep,double dTime)
{
    /* This picks the masses which are large enough to require
       randomization */

    int j,k,iDoRand,pick;
    double *newcoo;
    
    struct inLargeMass inLM;
    struct outLargeMass outLM;
    struct inMoveParticle inMove;
    struct inGetNeighbors inGNP;
    struct outGetNeighbors outMass,outReplace;
    
    inLM.fNumHillSphere = (FLOAT)msr->param.dRandBall;
    inLM.dMass = msr->param.dLargeMass;
    inLM.dCentMass = 1; /* for now  assume solar system */
    inLM.dOrbRad = pow((inLM.dCentMass*msr->param.PP.dOrbFreq),(-2.0/3));
    inMove.PP = msr->param.PP; /* struct copy */
    outLM.n = 0;
    inGNP.PP = msr->param.PP; /* struct copy */
    inGNP.dTime = dTime;
    
    /* How many Large Masses are eligible to be moved? Returns arrays of
       PARTICLEs and radii, and total number (scalar) */
    pstFindLargeMasses(msr->pst,&inLM,sizeof(inLM),&outLM,NULL);

    if (outLM.n) {
	/* Check for overlap between the large particles' Hill spheres */
	iDoRand=msrCheckLargeMassOverlap(outLM.p,outLM.dRadius,outLM.n);
	/* Pick a random position which will become the new location
	   of the large mass, keep z the same, ensuring there is no
	   overlap in the new coordinates.. */
	if (iDoRand) {		    
	    printf("Randomizing Large Masses...\n");
	    pick=msrPickMass(outLM.p,outLM.n);
	    
	    msrPickNewCoordinates(outLM.p,outLM.n,outLM.dRadius,msr->param.dxPeriod,msr->param.dyPeriod,&newcoo,pick);
    
	    inGNP.dDist=outLM.dRadius[pick];
	    inGNP.id=outLM.p[pick].iOrder;
	    for (j=0;j<3;j++)
		inGNP.x[j]=outLM.p[pick].r[j];
	    pstGetNeighborParticles(msr->pst,&inGNP,sizeof(inGNP),&outMass,NULL);
	    for (j=0;j<3;j++) 
		inGNP.x[j]=newcoo[j];
	
	    pstGetNeighborParticles(msr->pst,&inGNP,sizeof(inGNP),&outReplace,NULL);

	    /* Now do the switching */
	    /* Replacement Particles */
	    for (k=0;k<3;k++) 
		inMove.dOrigin[k]=outLM.p[pick].r[k];
	    for (j=0;j<outReplace.n;j++) {	    
		inMove.p=outReplace.p[j];
		for (k=0;k<3;k++) 
		    inMove.dRelx[k]=outReplace.p[j].r[k] - newcoo[k];
		pstMoveParticle(msr->pst,&inMove,sizeof(inMove),NULL,NULL);
		}
	    /* Large Mass Particle */
	    for (j=0;j<3;j++)
		inMove.dOrigin[j]=newcoo[j];
	    inMove.p = outLM.p[pick];
	    for (k=0;k<3;k++) 
		inMove.dRelx[k]=0.0;
	    pstMoveParticle(msr->pst,&inMove,sizeof(inMove),NULL,NULL);
	    /* Large Mass Neighbor Particles */
	    for (j=0;j<outMass.n;j++) {    
		inMove.p = outMass.p[j];
		for (k=0;k<3;k++)
		    inMove.dRelx[k]=outMass.p[j].r[k] - outLM.p[pick].r[k];   
		pstMoveParticle(msr->pst,&inMove,sizeof(inMove),NULL,NULL);
		}
	    msrRandomLog(iStep,outLM.p[pick],outMass.p,outMass.dSep2,outMass.n,outReplace.p,outReplace.dSep2,outReplace.n,newcoo,inGNP.dDist);	
	    }
	msr->param.iNextRandomization = msrGetNextRandomTime(msr->param.iRandStep,iStep);	
	}    
    }

#endif /* SLIDING_PATCH */

#endif /* RORY_EXTENSION */

void msrDoCollisions(MSR msr,double dTime,double dDelta)
{
	/*
	** Performs smooth operation to determine if a collision occurs
	** in the next interval. If so, the collision is processed and
	** collision flags are updated so that subsequent searches in the
	** interval can use resmooth over far fewer particles. This
	** continues until no further collisions occur in the interval.
	*/

#ifdef DEM
	return; /*DEBUG this function is ignored in DEM*/
#endif

	static double dFirstTime = DBL_MAX;

	struct inSmooth smooth;
	struct outNextCollision next;
	struct inGetColliderInfo inGet;
	struct outGetColliderInfo outGet;
	struct inDoCollision inDo;
	struct outDoCollision outDo;
	struct inResetColliders reset;
	COLLIDER *c1 = &inDo.Collider1,*c2 = &inDo.Collider2;
	double dSec;
	unsigned int nOvrlp=0,nCol=0,nMis=0,nMrg=0,nBnc=0,nFrg=0;
#ifdef COLLMOD_ZML
	unsigned int nExp=0;
#endif
	int bFirstPass;

#ifdef RORY_EXTENSION
	struct inSetBall setball; 
#endif

#ifdef AGGS_IN_PATCH
	/* RP-DEBUG: 6/5/09 hang/radial clearing errors */
	int iInterProcColl = 0; 
	int iRepeatColl = 0;
#endif

#ifdef WALLS
	unsigned int nDeath=0;
#endif

#ifdef WALLS
	if (msr->param.nSmooth < 1)
		return; /* might hit a wall */
#else
	if (msr->param.nSmooth < 2)
		return; /* won't find any colliders */
#endif
	if (dFirstTime == DBL_MAX)
		dFirstTime = dTime;
	if (msr->param.bVStep)
		printf("Start collision search (dTime=%e,dDelta=%e,run step=%.6g)...\n",dTime,dDelta,dDelta==0.0?0.0:(dTime - dFirstTime)/dDelta);
	dSec = msrTime(msr);
	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	/* following in pst.h (inSmooth) */
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0;
	smooth.iSmoothType = SMX_COLLISION;
	smooth.dfBall2OverSoft2 = 0.0; /* No softening limit */
	/* following in smoothfcn.h (SMF) */
	smooth.smf.dTime = dTime;
	smooth.smf.dStart = 0.0;
	smooth.smf.dEnd = dDelta;
	smooth.smf.bAllowSimulColl = msr->param.bAllowSimulColl;
	smooth.smf.dCollapseLimit = msr->param.CP.dCollapseLimit;
	smooth.smf.iOverlapOption = msr->param.CP.iOverlapOption;
	smooth.smf.bStrictOverlap = msr->param.CP.bStrictOverlap;
	smooth.smf.dBackstepLimit = msr->param.CP.dBackstepLimit;
	smooth.smf.dAdjPosLimit = msr->param.CP.dAdjPosLimit;
	smooth.smf.bOverlap = 0;
#ifdef AGGS
	smooth.smf.bAggsSolveQuartic = msr->param.CP.bAggsSolveQuartic;
#endif
#ifdef SLIDING_PATCH
	smooth.smf.PP = msr->param.PP; /* struct copy */
#endif
#ifdef WALLS
	smooth.smf.WP = msr->param.CP.WP; /* struct copy (needed for quartic-solver-stuck-particles) */
#endif

	inDo.bPeriodic = smooth.bPeriodic;
#ifdef SLIDING_PATCH
	inDo.dTime = dTime;
	inDo.PP = msr->param.PP;
#endif
	bFirstPass = 1;

#ifdef COLLMOD_ZML
	inDo.DB = msr->param.CP.DB; /* struct copy (global dust parameters) */	
#endif

	do {
		assert(smooth.smf.dStart == 0.0 || !smooth.smf.bOverlap); /* safety check */
		if (msr->param.dBallVelFact > 0.0) { /* default is use nSmooth */
#ifdef RORY_EXTENSION
			/* search ball determined by particle size and speed */
			if (msr->iTreeType != MSR_TREE_DENSITY) {

				/* uh oh!  Rory doesn't use gravity tree, so this is never called! */

				setball.dDelta = dDelta;
				setball.dBallVelFact = msr->param.dBallVelFact;
				/* set ball BEFORE building tree! */
				pstSetBall(msr->pst,&setball,sizeof(setball),NULL,NULL);
				msrBuildTree(msr,0,-1.0,1);
				}
			pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#else
			assert(0);
#endif /* RORY_EXTENSION */
			}
		else {
			/* search ball determined by nSmooth nearest neighbors */
			if (bFirstPass) {
				if (msr->iTreeType != MSR_TREE_DENSITY)
					msrBuildTree(msr,0,-1.0,1);
				pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
				bFirstPass = 0;
#ifdef RUBBLE_ZML
				if (dDelta == msr->param.dDelta && smooth.smf.dStart == 0.0) {
					struct outRubbleCheckForKDKRestart out;
					pstRubbleCheckForKDKRestart(msr->pst,NULL,0,&out,NULL);
					if (out.bRestart) {
						msr->param.CP.bDoRubbleKDKRestart = 1;
						/* flag reset in msrTopStepKDK() */
						return;
						}
					}
#elif defined(COLLMOD_ZML)
				if (dDelta == msr->param.dDelta && smooth.smf.dStart == 0.0) {
					struct outCollModCheckForKDKRestart out;
					pstCollModCheckForKDKRestart(msr->pst,NULL,0,&out,NULL);
					if (out.bRestart) {
						msr->param.CP.bDoCollModKDKRestart = 1;
						/* flag reset in msrTopStepKDK() */
						return;
						}
					}
#endif /* RUBBLE_ZML, COLLMOD_ZML */
				}
			else {
				assert(msr->iTreeType == MSR_TREE_DENSITY);
				/* following assumes inSmooth and inReSmooth are identical */
				/* NOTE: pstReSmooth() has caused problems in parallel
				  (MDL_CACHE_LINE): see DCR's e-mails in pkd folder
				  around end of Oct '01 & Jul '03.  Also mid Nov '09!  The assert is in mdlAquire().  */
				pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
				}
			}
		pstNextCollision(msr->pst,NULL,0,&next,NULL);
		/*
		** The following assert ensures that no two collisions occur
		** at precisely the same instant.  Physically this is possible
		** but in our case it's probably a sign of trouble, so it's
		** not allowed for now, except for certain overlap handling
		** options.  Note CheckForCollision() asserts this too, for a
		** more limited case.  CAVEAT: the tests rely on *exact*
		** double-precision equality which is almost impossible to
		** achieve -- the test should really be fabs(dt - dt0) <
		** PRECISION, but this is expensive.  Consequently, a
		** simultaneous collision may be missed, resulting in an
		** overlap during the next step!  The practical reason for not
		** allowing simultaneous collisions is that it's hard (without
		** bookkeeping) to distinguish between truly simultaneous
		** collisions involving different particles and repeated
		** "collisions" between the same particles that are just
		** touching and not moving relative to one another, resulting
		** in an infinite loop.  It's for this reason that we also
		** require msr->param.CP.dEpsN > 0.
		**
		** UPDATE: A new run-time option, bAllowSimulColl, has been
		** added to allow simultaneous collisions to occur.  But we
		** log each instance and give a warning if it's happening a
		** lot...  For walls, this can be legitimate (round-off).
		*/
		if (next.dt == smooth.smf.dStart &&
			(msr->param.CP.iOverlapOption == OverlapIsError ||
			 msr->param.CP.iOverlapOption == OverlapBackstep ||
			 msr->param.CP.iOverlapOption == OverlapRepel)) {
			if (msr->param.bAllowSimulColl || msr->param.CP.iOverlapOption == OverlapBackstep) {
#if (INTERNAL_WARNINGS != 0)
				static int nWarn = 1;
				if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
					fprintf(stderr,"SIMULTANEOUS COLLISION/ZERO STEP WARNING #%i [t=%e,run step=%.16g]: %i & %i, next.dt = %g\n",
							nWarn,dTime,dDelta==0.0?0.0:(dTime - dFirstTime)/dDelta,next.iOrder1,next.iOrder2,next.dt);
				++nWarn;
#endif /* INTERNAL_WARNINGS */
				}
			else
				assert(0); /* simultaneous collisions not allowed */
			}
		/* overlaps are indicated by negative or (possibly) zero next.dt */
		if (next.dt <= 0.0)
			switch (msr->param.CP.iOverlapOption) {
			case OverlapIgnore:
				break;
			case OverlapIsError:
				if (msr->param.bAllowSimulColl && next.dt == smooth.smf.dStart)
					break;
				assert(0); /* shouldn't be here */
			case OverlapRepel:
				assert(next.dt == 0.0); /* simultaneous collision may be allowed (see above) */
				break;
			case OverlapBackstep:
				++nOvrlp;
				break;
			case OverlapAdjPos:
			case OverlapMerge:
				assert(smooth.smf.dStart == 0.0); /* (next.dt encoded with degree of overlap) */
				smooth.smf.bOverlap = 1; /* use this flag to keep checking for overlaps */
				++nOvrlp;
				break;
			default:
				assert(0); /* unrecognized overlap option */
				}
		else
			smooth.smf.bOverlap = 0; /* (can only be set if AdjPos or Merge selected) */
		/* process the collision */
		if (COLLISION(next.dt)) {
			assert(next.iOrder1 >= 0);
#ifndef WALLS
			assert(next.iOrder2 >= 0);
#endif
			inDo.dt = next.dt;
			inGet.iOrder = next.iOrder1;
			pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
			*c1 = outGet.Collider; /* struct copy */
			assert(c1->id.iOrder == inGet.iOrder);
#ifdef WALLS
			if (next.iOrder2 < 0) { /* wall collision */
				c2->id.iOrder = next.iOrder2;
				/* the rest is mostly for logging purposes, but also sanity checks */
				c2->id.iPid = -1; /* must be < 0: cf. pstDoCollision() */
				c2->id.iIndex = -1;
				c2->fMass = DBL_MAX;
				c2->fRadius = DBL_MAX;
				c2->r[0] = c2->r[1] = c2->r[2] =
				  c2->v[0] = c2->v[1] = c2->v[2] =
				  c2->w[0] = c2->w[1] = c2->w[2] = 0.0;
				c2->dt = next.dt;
				c2->iColor = next.iOrder2;
				}
			else {
#endif
			inGet.iOrder = next.iOrder2;
			pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
			*c2 = outGet.Collider;
			assert(c2->id.iOrder == inGet.iOrder);
#ifdef WALLS
				}
#endif
			inDo.CP = msr->param.CP;
			inDo.dTime = dTime;
			inDo.dDelta = msr->param.dDelta; /* this stuff for Rory's binary merging... */
			inDo.iStartStep = msr->param.iStartStep;
#ifdef AGGS
			c1->agg.bAssigned = c2->agg.bAssigned = 0;
			/*
			** Must advance any aggregates to impact time now by
			** integrating the Euler equations.  Otherwise integrating
			** forward from the START of the step after a collision
			** would mean that the "ghost" orientation of the
			** aggregate after backdrifting would need to be computed,
			** requiring yet more expensive Runge-Kutta integrations.
			** We compromise by simply integrating forward as we go,
			** keeping track of the most recent update.  Note: if the
			** overlap flag is set, the aggregate data is simply
			** copied without forward integration.
			*/
			if (COLLIDER_IS_AGG(c1)) {
				int iAggIdx;
				Aggregate *agg;
				iAggIdx = COLLIDER_AGG_IDX(c1);
				agg = &msr->pAggs[iAggIdx];
				msrAggsAdvance(msr,iAggIdx,agg,smooth.smf.bOverlap ? 0.0 : next.dt);
				c1->agg = *agg; /* struct copy */
			}
			if (COLLIDER_IS_AGG(c2)) {
				int iAggIdx;
				Aggregate *agg;
				iAggIdx = COLLIDER_AGG_IDX(c2);
				agg = &msr->pAggs[iAggIdx];
				msrAggsAdvance(msr,iAggIdx,agg,smooth.smf.bOverlap ? 0.0 : next.dt);
				c2->agg = *agg; /* struct copy */
			}
			inDo.iAggNewIdx = msr->iAggNewIdx;
#endif /* AGGS */
#ifdef RUBBLE_ZML
			/* decide now if collision outcome is forced */
			/* ultimately want to check CP.iOutcomes to ensure rubble outcome OK */
			if (c1->iColor != PLANETESIMAL && c2->iColor != PLANETESIMAL)
				; /* don't output if rubble - rubble collision */
			else
				printf("Time = %e, dt = %e\n",dTime,next.dt);
			if (c1->iColor == PLANETESIMAL && c2->iColor == PLANETESIMAL) {
				printf("COLLISION: planetesimal - planetesimal\n");
				inDo.CP.iRubForcedOutcome = RUB_FORCED_NONE;
				}
			else if (c1->iColor != PLANETESIMAL && c2->iColor != PLANETESIMAL) {
				RUB_CLOCK *prc1=NULL,*prc2=NULL;
				double dt;
				int i;
				for (i=0;i<msr->re.nEvents;i++) {
					if (msr->re.rc[i].iColor == c1->iColor)
						prc1 = &msr->re.rc[i];
					if (msr->re.rc[i].iColor == c2->iColor)
						prc2 = &msr->re.rc[i];
					}
				assert(prc1 && prc2);
				if (prc1->dTStartMergePhase >= prc2->dTStartMergePhase)
					dt = prc1->dTStartMergePhase;
				else
					dt = prc2->dTStartMergePhase;
				if (dTime + next.dt >= dt) 
					inDo.CP.iRubForcedOutcome = RUB_FORCED_MERGE;
				else 
					inDo.CP.iRubForcedOutcome = RUB_FORCED_BOUNCE;
				}
			/* May not want a perfect merge if planetesimal hits rubble 10.27.03 */
			else { /* one collider is planetesimal */
				inDo.CP.iRubForcedOutcome = RUB_FORCED_MERGE;
				printf("COLLISION: planetesimal - rubble\n");
				}

			/* determine next available rubble color if needed */
/*			inDo.CP.iRubColor = RUB_BASE_COLOR + msr->re.nEvents;*/
			{
			int i;

			inDo.CP.iRubColor = RUB_BASE_COLOR;
			for (i=0;i<msr->re.nEvents;i++)
				if (msr->re.rc[i].iColor >= inDo.CP.iRubColor)
					inDo.CP.iRubColor = msr->re.rc[i].iColor + 1;
			}
			assert(inDo.CP.iRubColor != PLANETESIMAL); /* just in case */
#endif /* RUBBLE_ZML */

#ifdef COLLMOD_ZML
			puts(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>> COLLISION DETECTED <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
			printf("In msrDoCollision: time = %.9e, dt = %e\n",dTime,next.dt);/*DEBUG*/
#endif

#ifdef AGGS_IN_PATCH
			/* RP-DEBUG */
			/* Do some statistics on this collision */

			/* if pId for each collider is different, they
			   are on separate procs -> iInterProcColl++
			   ONCE */
			if(c1->id.iPid != c2->id.iPid) iInterProcColl++;
			/*if iPrevColl is NOT infinite -> iRepeatColl
			  (do for each collider)*/
			/* Note that this doesn't work accurately in
			   parallel!!!!! */
			if(msr->pst->plcl->pkd->pStore[c1->id.iIndex].iPrevCol < INT_MAX) 
			  iRepeatColl++;
			if(msr->pst->plcl->pkd->pStore[c2->id.iIndex].iPrevCol < INT_MAX)  
			  iRepeatColl++;
#endif

			/* do the collision! */
			pstDoCollision(msr->pst,&inDo,sizeof(inDo),&outDo,NULL);

#ifdef COLLMOD_ZML
			printf("Collision Outcome N = %i\n", outDo.nOut);
			assert(outDo.nOut >= 0);
//			ZML DEBUG 01.07.14 - looking for denisty issues
			{
			  int z;
			  double density_postcoll;
			  for (z=0;z<outDo.nOut;z++){
			    density_postcoll = (outDo.Out[z].fMass*2.0e33)/((4.0/3.0)*M_PI*(outDo.Out[z].fRadius*outDo.Out[z].fRadius*outDo.Out[z].fRadius)*(1.5e13*1.5e13*1.5e13));
			    printf("Density = %e\n", density_postcoll);
			    assert(density_postcoll > 0.009);
			  }
			}
#else
			assert(outDo.nOut > 0);
#endif

#ifdef AGGS_IN_PATCH
			/* [OBSOLETE]
			** In the particular case of 2 aggs merging, need to
			** account for possibility second collider is a ghost, in
			** which case the entire aggregate needs to be offset
			** prior to the merge.  NOTE: center of mass position and
			** velocity already adjusted in pkdAggsDoCollision().

			** RP: There's no reason to adjust the stored
			** wrapped position and velocity for a
			** particle joining an aggregate.  Wrapped
			** positions MUST remain in-patch, and
			** unwrapped positions will be taken care of
			** in msrAggsUpdate().
			*/
		/* 	if (outDo.dx != 0.0 || outDo.dy != 0.0) { */
/* 				struct inAggsInPatchOffset inOffset; */
/* 				assert(outDo.iOutcome == MERGE); */
/* 				assert(!msr->param.PP.bRandAzWrap); */
/* 				inOffset.iAggIdx = COLLIDER_AGG_IDX(c2); */
/* 				inOffset.dx = -outDo.dx; */
/* 				inOffset.dy = -outDo.dy; */
/* 				inOffset.dvy = -outDo.dvy; */
/* 				inOffset.bDoUnwrapped = 0; */
/* 				//pstAggsInPatchOffset(msr->pst,&inOffset,sizeof(inOffset),NULL,NULL); */
/* 				msrAggsInPatchGetUnwrapped(msr,COLLIDER_AGG_IDX(c2),dTime + next.dt); */
/* 				} */
#endif /* AGGS_IN_PATCH */
			msr->dTcoll += outDo.dT; /* account for kinetic energy loss */
			++nCol;
#ifdef COLLMOD_ZML
			// If nOut = 0 then we have an unresolved collision with no (planetesimal) remnants
			if (outDo.nOut == 0)
				outDo.iOutcome = EXPLODE;
#endif
			switch (outDo.iOutcome) {
			case MISS:
				++nMis;
				--nCol;
				break;
			case MERGE:
				++nMrg;
				break;
			case BOUNCE:
				++nBnc;
				break;
			case FRAG:
				++nFrg;
				break;
#ifdef COLLMOD_ZML
			case EXPLODE:
				++nExp;
				break;
#endif
#ifdef WALLS
			case DEATH:
				++nDeath;
				break;
#endif /* WALLS */
			default:
				assert(0); /* unknown outcome */
				}

			/*
			 ** Set SMOOTHACTIVE for colliders and any particles
			 ** predicted to collide with them in this interval to
			 ** force recomputation of their (and only their) collision
			 ** circumstances.  If there's a merger, the deleted particle
			 ** has ACTIVE set to zero, so we can still use the old tree
			 ** and ReSmooth().  Deleted particles are cleaned up outside
			 ** this loop.
			 */
			reset.iOrder1 = c1->id.iOrder;
			reset.iOrder2 = c2->id.iOrder;
			pstResetColliders(msr->pst,&reset,sizeof(reset),NULL,NULL);
#ifdef AGGS
			/*
			 ** For aggs, need to also handle effects on the center of
			 ** mass and/or spin of each agg.  This can change future
			 ** collision circumstances, so the following routines,
			 ** which all call msrAggsBackDrift(), have the side
			 ** effect of setting the consitutent particles in any
			 ** aggs involved in the collision to SMOOTHACTIVE,
			 ** forcing recomputation of possible future collision
			 ** circumstances during the remainder of the step.
			 */
			switch (outDo.iOutcome) {
			case MERGE:
				assert(outDo.nOut == 1);
				/* merge and store result in outDo for logging */
				/* (outDo.Out[0] is garbage at the moment) */
				/*
				** NOTE: overlaps are only detected at start of step;
				** force update times to be zero in this case, even for backstep
				** (cf. calls to msrAggsAdvance() above).
				*/

				msrAggsMerge(msr,c1,c2,smooth.smf.bOverlap ? 0.0 : next.dt,&outDo.Out[0]);
				
				break;
			case BOUNCE:
				assert(outDo.nOut == 2);
				/* copy any output agg structs to master storage and update */
				/*
				** NOTE: if an overlap occurred that resulted in
				** particle position adjustment, the agg(s) involved
				** must be updated.
				*/
				msrAggsBounce(msr,&outDo.Out[0],&outDo.Out[1],smooth.smf.bOverlap ? 0.0 : next.dt,
							  smooth.smf.bOverlap && msr->param.CP.iOverlapOption == OverlapAdjPos);
				
				break;
			case FRAG:
			    assert(outDo.nOut == 2); /* only 2 particles released at a time */
				/* need to identify and update leftover agg(s) */
				msrAggsFrag(msr,c1,c2,smooth.smf.bOverlap ? 0.0 : next.dt);
				break;
			default:
			    assert(0);
			}
#endif /* AGGS */
#ifdef SPRINGS
			if (outDo.iOutcome == FRAG)
				msrBreakSprings(msr,c1->id.iOrder,c2->id.iOrder);
#endif /* SPRINGS */
#ifndef AGGS
#ifndef SPRINGS
			if (outDo.iOutcome == FRAG) {
				/* note this sets msr->iTreeType to MSR_TREE_NONE */
				msrAddDelParticles(msr);
#ifdef RUBBLE_ZML
				printf("N = %i\n", msr->nDark); /* DEBUG */
#elif defined(COLLMOD_ZML)
				printf("System Total N = %i\n", msr->nDark); /* DEBUG */
#endif
				/* have to rebuild tree here and do new smooth... */
				bFirstPass = 1;
				/* also have to reactivate all particles */
				msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
				/* if nSmooth is too small, need to grow it -- ugly... */
				if (msr->param.nSmooth < 32) {
				  smooth.nSmooth = msr->param.nSmooth = (outDo.nOut < 32 ? outDo.nOut : 32);
				  if (msr->param.bVWarnings)
				    printf("WARNING: nSmooth set to %i\n",msr->param.nSmooth);				  
				}
#ifdef RUBBLE_ZML
				assert(inDo.CP.iRubForcedOutcome == RUB_FORCED_NONE);
				/* also need to enforce new timestep rungs for rubble */
				msrRubbleStep(msr);
				msrDtToRung(msr,0,msr->param.dDelta,1);
				/* and store event in array and set clocks */
				{
				RUB_CLOCK *rc;
				double dyn_time;

				assert(msr->re.nEvents < RUB_MAX_NUM_EVENTS);
				rc = &msr->re.rc[msr->re.nEvents++];
				rc->iColor = inDo.CP.iRubColor;
				dyn_time = 
					3.0/sqrt(c1->fMass/(4.0/3.0*M_PI*c1->fRadius*c1->fRadius*c1->fRadius));
				/* make sure dynamical time is consistent with minimum step */
				if (msr->param.dDelta/(1 << (msr->param.iMaxRung - 1)) >= 
					   0.05*dyn_time) /* DEBUG*/
					printf("min step = %e dyn_time = %e dens = %e dens = %e\n", 
						   msr->param.dDelta/(1 << (msr->param.iMaxRung - 1)), dyn_time,
						  c1->fMass/(4.0/3.0*M_PI*c1->fRadius*c1->fRadius*c1->fRadius),
						    c2->fMass/(4.0/3.0*M_PI*c2->fRadius*c2->fRadius*c2->fRadius));
				assert(msr->param.dDelta/(1 << (msr->param.iMaxRung - 1)) <= 
					   0.05*dyn_time);
				rc->dTStartMergePhase = dTime + next.dt + 
					msr->param.CP.DB.iRubNumDynToBounce*dyn_time;
				rc->dTEndRubblePhase = rc->dTStartMergePhase +
					msr->param.CP.DB.iRubNumDynToMerge*dyn_time;
				}
#endif /* RUBBLE_ZML */
				}
#endif /* !SPRINGS */
#endif /* !AGGS */
			/*
			** Collisions are expected to occur in monotonic sequence
			** after the start of the drift interval, unless we're
			** still looking for/dealing with overlap conditions (or
			** bAllowSimulColl is set -- see above).
			*/
			if (!smooth.smf.bOverlap) {
				assert(next.dt > smooth.smf.dStart || (next.dt == smooth.smf.dStart && msr->param.bAllowSimulColl) || msr->param.CP.iOverlapOption == OverlapBackstep);
				if (next.dt > 0.0)
					smooth.smf.dStart = next.dt;
				else
					smooth.smf.dStart = 0.0; /* in case of overlap (with backstep), force clock back to start of step -- keep doing this until everything's fixed */
				}
#ifdef RUBBLE_ZML
			if (outDo.iOutcome == MERGE && 
				c1->iColor == PLANETESIMAL && c2->iColor == PLANETESIMAL) {

				/*
				 ** dust from interpolated collision result is added
				 ** to DustBins here.  Note: dust is also generated in
				 ** resolved (rubble pile collisions) -- see Cleanup
				 ** routines.
				 */
				struct inRubInterpCleanup in;
				struct outRubInterpCleanup out;

				assert(outDo.nOut == 1); /* paranoid! */
				assert(inDo.CP.iRubForcedOutcome == RUB_FORCED_NONE);
				in.DB = msr->param.CP.DB; /* struct copy */
				in.iOrder = outDo.Out[0].id.iOrder;
				in.dCentMass = msr->param.dCentMass;
				/*printf("Before pstRubInterpCleanup\n");*/
				pstRubInterpCleanup(msr->pst,&in,sizeof(in),&out,NULL);
				if (out.DustBinsInterp.dMass > 0.0) {
					if (out.iBin >= 0 && out.iBin < msr->param.CP.DB.nDustBins) {
#ifdef ORIGIN_HISTOGRAM
						MergeHistograms(msr->aDustBins[out.iBin].dMass, msr->aDustBins[out.iBin].origin_bins, out.DustBinsInterp.dMass, out.DustBinsInterp.origin_bins);
#endif /* ORIGIN_HISTOGRAM */
						msr->aDustBins[out.iBin].dMass += out.DustBinsInterp.dMass;
						}
					else {
#ifdef ORIGIN_HISTOGRAM
						MergeHistograms(msr->DustBinsTrash.dMass, msr->DustBinsTrash.origin_bins, out.DustBinsInterp.dMass, out.DustBinsInterp.origin_bins);
#endif /* ORIGIN_HISTOGRAM */
						msr->DustBinsTrash.dMass += out.DustBinsInterp.dMass; /* put dust outside bin range in trash */
						}
					}
				}

#elif defined(COLLMOD_ZML)

			{
				/*
				 ** dust from collision model result is added to DustBins here 
				 */

				if (outDo.DustBin.dMass > 0.0) {
					if (outDo.iDustBin >= 0 && outDo.iDustBin < msr->param.CP.DB.nDustBins) {
#ifdef ORIGIN_HISTOGRAM
						MergeHistograms(msr->aDustBins[outDo.iDustBin].dMass, msr->aDustBins[outDo.iDustBin].origin_bins, outDo.DustBin.dMass, outDo.DustBin.origin_bins);
#endif /* ORIGIN_HISTOGRAM */
						msr->aDustBins[outDo.iDustBin].dMass += outDo.DustBin.dMass;
					}
					else {
#ifdef ORIGIN_HISTOGRAM
						MergeHistograms(msr->DustBinsTrash.dMass, msr->DustBinsTrash.origin_bins, outDo.DustBin.dMass, outDo.DustBin.origin_bins);
#endif /* ORIGIN_HISTOGRAM */
						msr->DustBinsTrash.dMass += outDo.DustBin.dMass; /* put dust outside bin range in trash */
					}
				}
			}

#endif /* RUBBLE_ZML, COLLMOD_ZML */

			if (msr->param.iCollLogOption != COLL_LOG_NONE)
				msrDoCollLog(msr,c1,c2,&outDo,msr->param.iCollLogOption,next.dt,dTime);
#ifdef WALLS
			/* account for mass lost for msrMassCheck() */
			if (outDo.iOutcome == DEATH)
				msr->dDeathWallsMassLost += c1->fMass;
#endif
			}
		} while (COLLISION(next.dt) && (smooth.smf.dStart < smooth.smf.dEnd || smooth.smf.bOverlap));
	msrAddDelParticles(msr); /* clean up any deletions -- not always needed */
	if (msr->param.nSmooth > msr->N) {
		msr->param.nSmooth = msr->N;
		if (msr->param.bVWarnings)
			(void) fprintf(stderr,"WARNING: msrDoCollisions(): nSmooth reduced to %i\n",msr->param.nSmooth);
		}

	if (msr->param.bVStep) {
		(void) printf("%i collision%s: %i overlap%s, %i miss%s, %i merger%s, %i bounce%s, %i frag%s\n",
					  nCol,nCol==1?"":"s",nOvrlp,nOvrlp==1?"":"s",nMis,nMis==1?"":"es",
					  nMrg,nMrg==1?"":"s",nBnc,nBnc==1?"":"s",nFrg,nFrg==1?"":"s");

#ifdef AGGS_IN_PATCH
		/* RP-DEBUG */
		(void) printf(" %i collisions between procs, %i repeating colliders.\n",
			      iInterProcColl/2, iRepeatColl);
#endif

#ifdef COLLMOD_ZML
		if (nExp > 0)
			printf("%i particle%s in unresolved collisions.\n",nExp,nExp==1?"":"s");
#endif

#ifdef WALLS
		if (nDeath > 0)
			(void) printf("%i particle%s removed.\n",nDeath,nDeath==1?"":"s");
#endif

		(void) printf("Collision search completed, time = %g sec\n",msrTime(msr) - dSec);
		}
	}

void msrUnifGravGetData(MSR msr,const char achFilenameWithQuotes[])
{
	FILE *fp;
	struct parameters *p;
	TIME_VAR_PARAMS *tvp;
	char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
	int i;

	assert(msr != NULL && achFilenameWithQuotes != NULL);
	p = &msr->param;
	tvp = &p->sTimeVarParams;
	_msrStripQuotes(achFilenameWithQuotes,achFilename);
	if (!strlen(achFilename)) {
		tvp->nUnifGrav = 0;
		return;
		}
	_msrMakePath(p->achDataSubPath,achFilename,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
	if (!(fp = fopen(achFilename,"r"))) {
		fprintf(stderr,"Unable to open \"%s\"\n",achFilename);
		goto abort;
		}
	if (fscanf(fp,"%i",&tvp->nUnifGrav) != 1) {
		fprintf(stderr,"Expected number of uniform gravity vectors in \"%s\"\n",achFilename);
		goto abort;
		}
	if (tvp->nUnifGrav <= 0) {
		fprintf(stderr,"Invalid number of uniform gravity vectors in \"%s\"\n",achFilename);
		goto abort;
		}
	tvp->dtUnifGrav = (double *) malloc(tvp->nUnifGrav*sizeof(double));
	assert(tvp->dtUnifGrav != NULL);
	tvp->dxUnifGrav = (double *) malloc(tvp->nUnifGrav*sizeof(double));
	assert(tvp->dxUnifGrav != NULL);
	tvp->dyUnifGrav = (double *) malloc(tvp->nUnifGrav*sizeof(double));
	assert(tvp->dyUnifGrav != NULL);
	tvp->dzUnifGrav = (double *) malloc(tvp->nUnifGrav*sizeof(double));
	assert(tvp->dzUnifGrav != NULL);
	for (i=0;i<tvp->nUnifGrav;i++) {
		if (fscanf(fp,"%lf%lf%lf%lf",&tvp->dtUnifGrav[i],&tvp->dxUnifGrav[i],
				   &tvp->dyUnifGrav[i],&tvp->dzUnifGrav[i]) != 4) {
			fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
					achFilename,i);
			fprintf(stderr,"Make sure the first line gives the correct number of vectors to use.\n");
			goto abort;
			}
		if ((i == 0 && tvp->dtUnifGrav[0] < 0.0) ||
			(i > 0 && tvp->dtUnifGrav[i] <= tvp->dtUnifGrav[i - 1])) {
			fprintf(stderr,"Invalid time in \"%s\", entry %i: %g\n"
					"(times must be non-zero and increasing)\n",
					achFilename,i,tvp->dtUnifGrav[i]);
			goto abort;
			}
		}
	fclose(fp);
	return;
 abort:
	if (fp) fclose(fp);
	_msrExit(msr,1);
	}

void msrSetUnifGrav(MSR msr,double dTime)
{
	static int iCurIdx = 0;
	struct parameters *p = &msr->param;
	const TIME_VAR_PARAMS *tvp = &p->sTimeVarParams;

	if (tvp->nUnifGrav == 0)
		return; /* not using variable uniform gravity */

	for (;iCurIdx<tvp->nUnifGrav;iCurIdx++) {
		/* find where we are in the vector list */
		if (tvp->dtUnifGrav[iCurIdx] >= dTime)
			break;
		}

	if (iCurIdx == 0) {
		if (tvp->dtUnifGrav[iCurIdx] > dTime) {
			static int bFirstWarning = 1;
			if (bFirstWarning) {
				if (p->bVWarnings)
					fprintf(stderr,"WARNING: Time %g before uniform gravity list (using first value)\n",dTime);
				bFirstWarning = 0;
				}
			}
		p->dxUnifGrav = tvp->dxUnifGrav[0]; /* time index, not spatial index! */
		p->dyUnifGrav = tvp->dyUnifGrav[0];
		p->dzUnifGrav = tvp->dzUnifGrav[0];
		}
	else if (iCurIdx == tvp->nUnifGrav) {
		if (dTime > tvp->dtUnifGrav[iCurIdx - 1]) { /* should always be true */
			static int bFirstWarning = 1;
			if (bFirstWarning) {
				if (p->bVWarnings)
					fprintf(stderr,"WARNING: Time %g beyond uniform gravity list (using last value)\n",dTime);
				bFirstWarning = 0;
				}
			}
		p->dxUnifGrav = tvp->dxUnifGrav[iCurIdx - 1];
		p->dyUnifGrav = tvp->dyUnifGrav[iCurIdx - 1];
		p->dzUnifGrav = tvp->dzUnifGrav[iCurIdx - 1];
		}
	else {
		/* interpolate between values using shortest arc */
		Vector g0,g1,vAxis;
		Vector g1sub1,g1sub2;
		double t0,t1,dt,f,dNorm,dAngle;
		double g0mag, g1mag, g0scale, dNormsub1, dNormsub2, g01;
		t0 = tvp->dtUnifGrav[iCurIdx - 1]; /* start time of interval */
		t1 = tvp->dtUnifGrav[iCurIdx]; /* end time of interval */
		dt = t1 - t0;
		assert(dt > 0.0);
		f = (dTime - t0)/dt; /* current fraction of interval */
		g0[0] = tvp->dxUnifGrav[iCurIdx - 1]; /* starting gravity */
		g0[1] = tvp->dyUnifGrav[iCurIdx - 1];
		g0[2] = tvp->dzUnifGrav[iCurIdx - 1];
		g1[0] = tvp->dxUnifGrav[iCurIdx]; /* ending gravity */
		g1[1] = tvp->dyUnifGrav[iCurIdx];
		g1[2] = tvp->dzUnifGrav[iCurIdx];
		vectorCross(g0,g1,vAxis); /* get axis of rotation */
		dNorm = vectorMag(vAxis);		
/*		assert(dNorm > 0.0); */ /* reject invalid rotation */ 
		/* Added by Soko */
		g0mag = vectorMag(g0);
		g1mag = vectorMag(g1);
		g0scale = 1.0 + (g1mag - g0mag)/g0mag*f;
		if (dNorm > 0.0) {
			dAngle = atan2(dNorm,vectorDot(g0,g1)); /* angle between endpoint vectors */
			vectorScale(vAxis,1.0/dNorm,vAxis); /* make rotation axis into unit vector */
			vectorRotate(g0,vAxis,f*dAngle); /* rotate by needed fraction of angle */
			vectorScale(g0,g0scale,g0); /* scale g0 toward g1 */ 
			/* printf("T= %g yr: g0mag is %g after scaling and g0scale is %g.\n", dTime, vectorMag(g0), g0scale); */
		}
		else if (dNorm <= 0.0) {
		        dNormsub1 = vectorDot(g0,g1);
			dNormsub2 = fabs(dNormsub1);  
			assert(dNormsub2 > 0.0);  /* reject if one of the vector is null */   
			if ( dNormsub1 > 0.0) { /* no rotation needed */
				/* printf("g0mag is %g before scaling.\n", g0mag); */
				vectorScale(g0,g0scale,g0); /* scale g0 toward g1 */
				/* printf("T= %g yr: g0mag is %g after scaling and g0scale is %g.\n", dTime, vectorMag(g0), g0scale); */
			} 
			else if ( dNormsub1 < 0.0) { /* rotation by 180 degrees. */
				g01 = g0mag + g1mag;
				vectorScale(g1,1.0/g1mag,g1sub1); /* scale g1 to unit vector */ 
				vectorScale(g1sub1,g01*f,g1sub2); /* scale the unit vector by a fraction of the total length of two vectors */
				vectorAdd(g0,g1sub2,g0); /* change g0 toward g1 */
				/* printf("T= %g yr: g0mag is %g after scaling.\n", dTime, vectorMag(g0)); */
			}  
		}
		/* END: Added by Soko */ 
		p->dxUnifGrav = g0[0];
		p->dyUnifGrav = g0[1];
		p->dzUnifGrav = g0[2];
		}
	}

#ifdef JOHNNY

void msrFindEscapers(MSR msr)
{
	static double dStarMass = 0.0;
	static int nEsc = 0;

	struct inFindEscapers in;
	struct outFindEscapers out;

	if (dStarMass == 0.0) {
		dStarMass = msrMassCheck(msr,-2.0,"In msrFindEscapers");
		assert(dStarMass > 0.0);
		}

	in.dCentMass = msr->param.dCentMass;
	in.dStarMass = dStarMass;
	pstFindEscapers(msr->pst,&in,sizeof(in),&out,NULL);
	if (out.nEsc > 0) {
		nEsc += out.nEsc;
		printf("%i particle%s replaced...running total now %i.\n",out.nEsc,out.nEsc == 1 ? "" : "s",nEsc);
		/* particles moved, so need to update processors & tree */
		msrDomainDecomp(msr,0,1); /* redo domain decomposition, all rungs */
		msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
		msrBuildTree(msr,0,-1.0,0); /* rebuild gravity tree, all particles, ignore mass check */
		}
	}

void msrCheckForKepler(MSR msr)
{
	static int nKep = 0;

	struct inCheckForKepler in;
	struct outCheckForKepler out;

	in.dCentMass = msr->param.dCentMass;
	in.dDelta = msr->param.dDelta;
	assert(msr->param.iMaxRung == 1); /* no multistepping, for now */
	pstCheckForKepler(msr->pst,&in,sizeof(in),&out,NULL);
	if (out.nKep > 0) {
		nKep += out.nKep;
		printf("%i new particle%s entered critical radius...running total now %i.\n",out.nKep,out.nKep == 1 ? "" : "s",nKep);
		}
	}	

#endif /* JOHNNY */

#ifdef SPRINGS

void msrAssignSprings(MSR msr)
{
	/* initalizes springs structs and assigns springs to sufficiently proximal particles */

	struct inSmooth smooth;
	double dSec;

	if (msr->param.nSmooth < 2)
		return; /* need at least 1 neighbour per particle */

	if (msr->param.bVStart)
		(void) printf("Assigning springs...\n");

	dSec = msrTime(msr);

	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	msrDomainDecomp(msr,0,1);

	if (msr->iTreeType != MSR_TREE_DENSITY)
		msrBuildTree(msr,0,-1.0,1);

	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0; /* use read-only (not combiner) cache */
	smooth.iSmoothType = SMX_ASSIGN_SPRINGS;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
	smooth.smf.FO = msr->param.FO; /* struct copy */

	/* can't set searchball because size changes depending on particle size */

	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
	smooth.iSmoothType = SMX_COLOR_SPRINGS;
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);

	if (msr->param.bVStart)
		(void) printf("Spring assignment complete, time = %g sec.\n",msrTime(msr) - dSec);
	}

void msrDoSprings(MSR msr,double dTime)
{
	/* apply spring forces to sufficiently proximal particles */

	struct inSmooth smooth;
	double dSec;

	if (msr->param.nSmooth < 2)
		return; /* need at least 1 neighbour per particle */

	if (msr->param.bVStep)
		(void) printf("Applying spring forces...\n");

	dSec = msrTime(msr);

	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0; /* use read-only (not combiner) cache */
	smooth.iSmoothType = SMX_MOVESPECIALS_SPRINGS;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
	smooth.smf.dTime = dTime;
	smooth.smf.FO = msr->param.FO; /* struct copy */

	/*
	** Could set search ball here and do resmooth instead of smooth,
	** except we would still have to do a full smooth at the start of
	** the collision search anyway.
	*/

	if (msr->iTreeType != MSR_TREE_DENSITY)
		msrBuildTree(msr,0,-1.0,1);

	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
	smooth.bSymmetric = 1; /* use combiner cache */
	smooth.iSmoothType = SMX_DO_SPRINGS;
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
	smooth.bSymmetric = 0; /* use read-only (not combiner) cache */
	smooth.iSmoothType = SMX_COLOR_SPRINGS;
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);

	if (msr->param.bVStep)
		(void) printf("Spring force calculation complete, time = %g sec.\n",msrTime(msr) - dSec);
	}

void msrBreakSprings(MSR msr,int iOrder1,int iOrder2)
{
	/*
	** Removes all springs from particles iOrder1 and iOrder2, and all
	** springs from all other particles that are attached to these
	** two.  Called as a result of a fragmentation event.
	*/

	struct inBreakSprings in;

	in.iOrder1 = iOrder1;
	in.iOrder2 = iOrder2;
	in.bSpringsChanged = 0;
	pstBreakSprings(msr->pst,&in,sizeof(in),NULL,NULL);

        /* Recolor springs with a ReSmooth... */

	if (in.bSpringsChanged) {
	        struct inSmooth smooth;

		msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
		smooth.nSmooth = msr->param.nSmooth;
		smooth.bPeriodic = msr->param.bPeriodic;
		smooth.bSymmetric = 0; /* use read-only (not combiner) cache */
		smooth.iSmoothType = SMX_COLOR_SPRINGS;
		smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
		smooth.smf.FO = msr->param.FO; /* struct copy */
		pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
	        }
        }

#endif /* SPRINGS */

#ifdef CHARGE

void msrChargeZ(MSR msr)
{
	/*
	** Computes total charge and z-components of 1st & 2nd charge
	** distribution moments and applies them to particle z
	** accelerations.  This is intended for an axisymmetric geometry,
	** such as particles in a cylindrical container aligned with the z
	** axis and centered at x=0,y=0.
	*/

	struct ioChargeZMoments io;
	CHARGE_PARAMS *CP = &msr->param.CP.CP;

	pstChargeZGetMoments(msr->pst,NULL,0,&io,NULL);
	CP->dQ = io.CP.dQ; /* don't do struct copy else dArea overwritten... */
	CP->dZ = io.CP.dZ;
	CP->dSigZ = io.CP.dSigZ;
	assert(CP->dQ != 0.0); /* otherwise no net charge! */
	CP->dSigZ = CP->dSigZ/CP->dQ - (CP->dZ*CP->dZ)/(CP->dQ*CP->dQ);
	assert(CP->dSigZ > 0.0);
	CP->dSigZ = sqrt(CP->dSigZ);
	CP->dZ = CP->dZ/CP->dQ;
	io.CP.dArea = CP->dArea;
	pstChargeZApplyMoments(msr->pst,&io,sizeof(io),NULL,NULL);
	}

#endif /* CHARGE */

#ifdef DEM

void msrDEMCheckFile(MSR msr) {
	FILE *fp;
	char achTmp[PST_FILENAME_SIZE],achInFile[PST_FILENAME_SIZE];
	_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,achTmp);
	_msrMakePath(msr->pst->plcl->pszDataPath,achTmp,achInFile);
	if (msr->param.FO.DP.bReadDEMData == 0) {
		/* check if .dem file corresponding to starting file exists, and warn if so */
		snprintf(achTmp,PST_FILENAME_SIZE,"%s.dem",achInFile);
		fp = fopen(achTmp,"r");
		if (fp != NULL) {
			fprintf(stderr,"%s exists but is not read (bReadDEMData = 0).\n",achTmp);
			fclose(fp);
			}
		}
#ifdef DEM_TRACK_ORIENT
	if (msr->param.FO.DP.bReadOrient == 0) {
		/* repeat check for .ori file */
		snprintf(achTmp,PST_FILENAME_SIZE,"%s.ori",achInFile);
		fp = fopen(achTmp,"r");
		if (fp != NULL) {
			fprintf(stderr,"%s exists but is not read (bReadOrient = 0).\n",achTmp);
			fclose(fp);
			}
		}
#endif /* DEM_TRACK_ORIENT */
	}

void msrDEMRemoveOldOutput(MSR msr,int iStep) {
	/* Remove unwanted DEM files */
	int iRemStep, bRemove;	/* Should maybe change 'iOldStep' to be more general*/

	bRemove = FALSE;

	if (msr->param.FO.DP.nDEMOutputs > 0) {	/* Remove oldest DEM file not being kept */
		iRemStep = iStep - msr->param.FO.DP.nDEMOutputs*msrOutInterval(msr);
		if (iRemStep >= 0) {
			bRemove = TRUE;
			}
		}

	else if (msr->param.FO.DP.nDEMOutputs < -1) { /* Keep only every Nth DEM file */
		iRemStep = iStep % (msr->param.FO.DP.nDEMOutputs*msrOutInterval(msr));
		if (iRemStep != 0) {
			iRemStep = iStep;
			bRemove = TRUE;
			}
		}

	if (bRemove == TRUE) {
		char achTmp[PST_FILENAME_SIZE],achOldFile[PST_FILENAME_SIZE];
		snprintf(achOldFile,PST_FILENAME_SIZE - 4/*(for ".dem")*/,msr->param.achDigitMask,msrOutName(msr),iRemStep);
		_msrMakePath(msr->param.achDataSubPath,achOldFile,achTmp);
		_msrMakePath(msr->pst->plcl->pszDataPath,achTmp,achOldFile);
		snprintf(achTmp,PST_FILENAME_SIZE,"%s.dem",achOldFile);
		unlink(achTmp);
		}

	}

#ifdef DEM_FIXED_BALL

void msrDEMSetBall(MSR msr)
{
	/* IMPORTANT: this routine assumes particle sizes do not change
	   during the run... */
	struct inDEMSetBall setball;
	static double dFixedRadius = 0.0;
	static int bFirstCall = 1;

	int iFixedBallOption = msr->param.FO.DP.iFixedBallOption;
	
	setball.iFixedBallOption = iFixedBallOption;
	if (bFirstCall) {
		struct inDEMGetBallRadius in;
		struct outDEMGetBallRadius out;
		in.iFixedBallOption = iFixedBallOption;
		switch (iFixedBallOption) {
		case 0:
			pstDEMGetBallRadius(msr->pst,&in,sizeof(in),&out,NULL);
			dFixedRadius = out.dFixedRadius;
			printf("DEM_FIXED_BALL: set to maximum radius: %.8g\n",dFixedRadius);
			if (msr->param.bVWarnings) printf("(use caution if any particle is a LOT bigger than the others)\n");
			break;
		case 1:
			pstDEMGetBallRadius(msr->pst,&in,sizeof(in),&out,NULL);
			assert(msr->N > 0);
			dFixedRadius = out.dFixedRadius / msr->N;
			printf("DEM_FIXED_BALL: set to particle radius plus mean radius: %.8g\n",dFixedRadius);
			if (msr->param.bVWarnings) printf("(use caution if any particle is a LOT bigger than the others)\n");
			break;
		case 2:
			dFixedRadius = 0.0;
			printf("DEM_FIXED_BALL: set to particle radius plus next largest\n");
			break;
		default:
			assert(0); /* shouldn't get here */
			}
		bFirstCall = 0;
		}
	/* if particle search balls are not reset, could put following above */
	setball.dFixedRadius = dFixedRadius;
	pstDEMSetBall(msr->pst,&setball,sizeof(setball),NULL,NULL);
	}

void msrDEMFixedBallNextLargest(MSR msr) {
    pstDEMFixedBallNextLargest(msr->pst,NULL,0,NULL,NULL);
	}

#endif /* DEM_FIXED_BALL */

void msrAssignDEM(MSR msr)
{
	/* initializes DEM structs */
	struct inSmooth smooth;
	double dSec;

	if (msr->param.bVStart)
		(void) printf("Initializing DEM parameters...\n");

	dSec = msrTime(msr);

	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	msrDomainDecomp(msr,0,1);

	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0; /* use read-only (not combiner) cache */
	smooth.iSmoothType = SMX_ASSIGN_DEM;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
	smooth.smf.FO = msr->param.FO; /* struct copy */

	/*DEBUG! need to initialize these to prevent garbage being sent to smooth
	 (which begs the question why are we using a smooth function here??)*/
	smooth.smf.dTime = 0.0;
#ifdef AGGS
	smooth.smf.bAggsSolveQuartic = msr->param.CP.bAggsSolveQuartic;
#endif
#ifdef SLIDING_PATCH
	smooth.smf.PP = msr->param.PP; /* struct copy */
#endif
#ifdef WALLS
	smooth.smf.WP = msr->param.CP.WP; /* struct copy */
#endif

	if (msr->iTreeType != MSR_TREE_DENSITY)
		msrBuildTree(msr,0,-1.0,1);

#ifdef DEM_FIXED_BALL
	msrDEMSetBall(msr);
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#else
	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#endif /* DEM_FIXED_BALL */

	if (msr->param.bVStart)
		printf("DEM initialization complete, time = %g sec.\n",msrTime(msr) - dSec);

	}

void msrDoDEM(MSR msr,double dTime)
{
	/* apply DEM forces to overlapping particles */
	struct inSmooth smooth;

#ifdef WALLS
	struct outGetMassChange out;
#endif

	double dSec;

	/*DEBUG: should revisit and ensure no extra full smooths if also doing springs*/

	if (msr->param.bVStep)
		(void) printf("Applying DEM forces...\n");

	dSec = msrTime(msr);

	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 1; /* use combiner cache */
	smooth.iSmoothType = SMX_DO_DEM;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */
	smooth.smf.dTime = dTime;
	smooth.smf.FO = msr->param.FO; /* struct copy */
#ifdef AGGS
	smooth.smf.bAggsSolveQuartic = msr->param.CP.bAggsSolveQuartic;
#endif
#ifdef SLIDING_PATCH
	smooth.smf.PP = msr->param.PP; /* struct copy */
#endif
#ifdef WALLS
	smooth.smf.WP = msr->param.CP.WP; /* struct copy */
#endif

	if (msr->iTreeType != MSR_TREE_DENSITY)
		msrBuildTree(msr,0,-1.0,1);
	
#ifdef DEM_FIXED_BALL
	msrDEMSetBall(msr);
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#else
	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#endif /* DEM_FIXED_BALL */
#ifdef WALLS
	/*
	** Currently, the only implemented mass change is particle loss
	** due to "death" walls.  We account for that here.
	*/
	pstGetMassChange(msr->pst,NULL,0,&out,NULL);
	/*assert(out.dMassChange <= 0.0);*/ /* only mass loss expected */ /*DEBUG!*/
	msr->dDeathWallsMassLost -= out.dMassChange;
#endif
	msrAddDelParticles(msr);

	if (msr->param.bVStep)
		printf("DEM force calculation complete, time = %g sec.\n",msrTime(msr) - dSec);
	}

#ifdef DEM_WALLS_REACT

void msrWallsKick(MSR msr, double dt)
{
    /*
    ** Only applies to reactive walls since reactive wall velocities
    ** are found by integrating the equations of motion.
    */

    /* dt is one half current drift interval */

    WALL_PARAMS *WP = &msr->param.CP.WP;
    WALL *w;
    WALL_DATA *wd;
    int i, k;

    for (i = 0; i < WP->nWalls; i++) {
        w = &WP->pWalls[i];
        wd = &w->wd;
        
        if (wd->dMass == 0.0)
            continue; /* not a reactive wall */

	for (k = 0; k < 3; k++)
	    wd->vVel[k] += wd->accel[k] * dt;
	}
    }

void msrDEMWallsReact(MSR msr, double dDelta)
{
    /*
    ** Called during the gravity step to get acclerations and torques
    ** on reactive walls/assemblies.  Wall orientations are updated in
    ** msrWallsMove().
    */

    /* allow non-infinite inertia for walls */
    
    struct outDEMWallsReact out;//defined in pst.h
    WALL_PARAMS *WP = &msr->param.CP.WP;
    WALL *w;
    WALL_DATA *wd;
    double dWallsReactMass[MAX_NUM_WALL_ASSEMBLIES];
    double dWallsXAccel[MAX_NUM_WALL_ASSEMBLIES];
    double dWallsYAccel[MAX_NUM_WALL_ASSEMBLIES];
    double dWallsZAccel[MAX_NUM_WALL_ASSEMBLIES];
    double dTotalWallMass = 0;
    Matrix A, A_T;
    Vector vTmp;
    int i;
    
    /*
     ** Collect particle-level net force.
     **
     ** After the following call, "out" contains a double specifying the total
     ** force that the particles are exerting on the inertial wall assemblage.
     */
    
    pstDEMWallsReact(msr->pst,NULL,0,&out,NULL);
    
    for (i=0;i<MAX_NUM_WALL_ASSEMBLIES;i++)
        dWallsReactMass[i]=0;
    
    for (i=0;i<WP->nWalls;i++) {//Getting the wall mass
        w = &WP->pWalls[i];
        wd = &w->wd;
        
        if (wd->dMass == 0.0)
            continue;
        
        dWallsReactMass[wd->iAssembly] += wd->dMass;
    }
    
    /*Check if there are any inertial walls*/
    for (i=0; i<MAX_NUM_WALL_ASSEMBLIES;i++)
        dTotalWallMass += dWallsReactMass[i];
    
    if (dTotalWallMass == 0.0)
        return;
    
    double axConst[6] = {0, 0, 0, 0, 0, 0};
    double ayConst[6] = {0, 0, 0, 0, 0, 0};
    double azConst[6] = {0, 0, 0, 0, 0, 0};
    
    /*Find accelerations and torques based on Reaction force and gravity*/
	  
    for (i=0; i<MAX_NUM_WALL_ASSEMBLIES;i++) {
        if (dWallsReactMass[i]==0)
            continue;
        dWallsXAccel[i] = out.dTotalXForceFromParticles[i]/dWallsReactMass[i] + axConst[i] + msr->param.dxUnifGrav;
        dWallsYAccel[i] = out.dTotalYForceFromParticles[i]/dWallsReactMass[i] + ayConst[i] + msr->param.dyUnifGrav;
        dWallsZAccel[i] = out.dTotalZForceFromParticles[i]/dWallsReactMass[i] + azConst[i] + msr->param.dzUnifGrav;
    }

    for (i=0;i<WP->nWalls;i++) {
        w = &WP->pWalls[i];
        wd = &w->wd;
        if (wd->dMass == 0.0)
            continue;

	wd->accel[0] = dWallsXAccel[wd->iAssembly];
	wd->accel[1] = dWallsYAccel[wd->iAssembly];
	wd->accel[2] = dWallsZAccel[wd->iAssembly];
        
        A[0][0]=wd->vPxj[0];
        A[1][0]=wd->vPxj[1];
        A[2][0]=wd->vPxj[2];
        A[0][1]=wd->vPyj[0];
        A[1][1]=wd->vPyj[1];
        A[2][1]=wd->vPyj[2];
        A[0][2]=wd->vPzj[0];
        A[1][2]=wd->vPzj[1];
        A[2][2]=wd->vPzj[2];
        
        wd->torque[0]=out.dTotalXTorqueFromParticles[wd->iAssembly];
        wd->torque[1]=out.dTotalYTorqueFromParticles[wd->iAssembly];
        wd->torque[2]=out.dTotalZTorqueFromParticles[wd->iAssembly];
        wd->moments[0]=wd->vIxj[0];
        wd->moments[1]=wd->vIyj[1];
        wd->moments[2]=wd->vIzj[2];
        //if (wd->torque[0]!=0 || wd->torque[1]!=0 || wd->torque[2]!=0)
        //printf("%.8g %.8g %.8g\n",wd->torque[0],wd->torque[1],wd->torque[2]);
        /*converted torque from space to body frame*/
        matrixTranspose(A, A_T);
        matrixTransform(A_T, wd->torque, vTmp);/* x_body=A_T*x_space */
        vectorCopy(vTmp, wd->torque);
    }
}

#endif /* DEM_WALLS_REACT */

void msrDEMZeroSmallMotions(MSR msr)
{
	static int bFirstCall = 1;
	static struct inDEMZeroSmallMotions in;

	if (bFirstCall == 1) {
		in.dAccCritSq = msr->param.FO.DP.dAccCrit*msr->param.FO.DP.dAccCrit;
		in.dDeltaSq = msr->param.dDelta*msr->param.dDelta; /* change for multistepping?... */
		bFirstCall = 0;
		}
	pstDEMZeroSmallMotions(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrDEMStats(MSR msr,double dTime,int iStep)
{
	/* collect DEM statistics! */

	static int cumul_ohist[DEM_NUM_OVERLAP_BINS];
	static int cumul_ahist[DEM_NUM_COS_A_BINS];
	static int cumul_shist[DEM_NUM_S_BINS];
	static int bFirstCall = TRUE;

	struct inSmooth smooth;
	struct outDEMStats out;
	char achFile[MAXPATHLEN],achTemp[MAXPATHLEN];
	FILE *fp;
	double dSec = 0.0;
	int i;

	if (msr->param.bVStep) {
		printf("Collecting DEM statistics...\n");
		dSec = msrTime(msr);
		}

	/* initialize accumulators if first call (data lost if restart) */

	if (bFirstCall) {
		for (i=0;i<DEM_NUM_OVERLAP_BINS;i++)
			cumul_ohist[i] = 0.0;
		for (i=0;i<DEM_NUM_COS_A_BINS;i++)
			cumul_ahist[i] = 0.0;
		for (i=0;i<DEM_NUM_S_BINS;i++)
			cumul_shist[i] = 0.0;
		bFirstCall = FALSE;
		}

	/* first collect particle-level stats using smooth */

	msrActiveType(msr,TYPE_ALL,TYPE_ALLACTIVE);
	smooth.nSmooth = msr->param.nSmooth;
	smooth.bPeriodic = msr->param.bPeriodic;
	smooth.bSymmetric = 0; /* don't need combiner cache */
	smooth.iSmoothType = SMX_DEM_STATS;
	smooth.dfBall2OverSoft2 = 0.0; /* ignore softening */

	smooth.smf.dTime = dTime;
	smooth.smf.FO = msr->param.FO; /* struct copy (not needed currently) */
#ifdef AGGS
	smooth.smf.bAggsSolveQuartic = msr->param.CP.bAggsSolveQuartic;
#endif
#ifdef SLIDING_PATCH
	smooth.smf.PP = msr->param.PP; /* struct copy */
#endif
#ifdef WALLS
	smooth.smf.WP = msr->param.CP.WP; /* struct copy */
#endif

	if (msr->iTreeType != MSR_TREE_DENSITY)
		msrBuildTree(msr,0,-1.0,1);

#ifdef DEM_FIXED_BALL
	msrDEMSetBall(msr);
	pstReSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#else
	pstSmooth(msr->pst,&smooth,sizeof(smooth),NULL,NULL);
#endif /* DEM_FIXED_BALL */

	/* now collect aggregate statistics */

	pstDEMStats(msr->pst,NULL,0,&out,NULL);

	/* accumulate */

	for (i=0;i<DEM_NUM_OVERLAP_BINS;i++)
		cumul_ohist[i] += out.DEMStats.pOverlapHist[i];
	for (i=0;i<DEM_NUM_COS_A_BINS;i++)
		cumul_ahist[i] += out.DEMStats.pCosAHist[i];
	for (i=0;i<DEM_NUM_S_BINS;i++)
		cumul_shist[i] += out.DEMStats.pSHist[i];

	/* write data to file */

	snprintf(achFile,MAXPATHLEN - 10,msr->param.achDigitMask,msrOutName(msr),iStep);
	strncat(achFile,".demstats",MAXPATHLEN - 1);
	_msrMakePath(msr->param.achDataSubPath,achFile,achTemp);
	_msrMakePath(msr->lcl.pszDataPath,achTemp,achFile);
	fp = fopen(achFile,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",achFile);
		return;
		}

	fprintf(fp,"# DEM Stats file\n");
	fprintf(fp,"# line 11 = current time (pkdgrav units)\n");
	fprintf(fp,"# line 12 = number of particles\n");
	fprintf(fp,"# line 13 = global minimum distance & maximum relative neighbor speed (pkdgrav units)\n");
	fprintf(fp,"# line 14 = number of overlap bins\n");
	fprintf(fp,"# lines 15-%i = overlap bins (see dem.c for min and max)\n",14 + DEM_NUM_OVERLAP_BINS);
	fprintf(fp,"# line %i = number of cos(alpha) bins\n",15 + DEM_NUM_OVERLAP_BINS);
	fprintf(fp,"# lines %i-%i = cos(alpha) bins, min = -1, max = +1\n",16 + DEM_NUM_OVERLAP_BINS,15 + DEM_NUM_OVERLAP_BINS + DEM_NUM_COS_A_BINS);
	fprintf(fp,"# line %i = number of S/R bins\n",16 + DEM_NUM_OVERLAP_BINS + DEM_NUM_COS_A_BINS);
	fprintf(fp,"# lines %i-%i = S/R bins (see dem.c for min and max)\n",17 + DEM_NUM_OVERLAP_BINS + DEM_NUM_COS_A_BINS,16 + DEM_NUM_OVERLAP_BINS + DEM_NUM_COS_A_BINS + DEM_NUM_S_BINS);
	fprintf(fp,"%.16e\n",dTime);
	fprintf(fp,"%i\n",msr->N);
	fprintf(fp,"%.16e %.16e\n",out.DEMStats.fDistMin,out.DEMStats.fSpeedMax);
	fprintf(fp,"%i\n",DEM_NUM_OVERLAP_BINS);
	for (i=0;i<DEM_NUM_OVERLAP_BINS;i++)
		fprintf(fp,"%i %i\n",out.DEMStats.pOverlapHist[i],cumul_ohist[i]);
	fprintf(fp,"%i\n",DEM_NUM_COS_A_BINS);
	for (i=0;i<DEM_NUM_COS_A_BINS;i++)
		fprintf(fp,"%i %i\n",out.DEMStats.pCosAHist[i],cumul_ahist[i]);
	fprintf(fp,"%i\n",DEM_NUM_S_BINS);
	for (i=0;i<DEM_NUM_S_BINS;i++)
		fprintf(fp,"%i %i\n",out.DEMStats.pSHist[i],cumul_shist[i]);

	fclose(fp);

	if (msr->param.bVStep)
		printf("DEM stats calculation complete, time = %g sec.\n",msrTime(msr) - dSec);
	}

#endif /* DEM */

#ifdef SPINUP

void msrSpinupAlignAxes(MSR msr, Vector r_com, Vector vVector)
{
	/*
	** Helper routine for msrSpinup() (see below).  Rotates the
	** aggregate so vVector (either the spin or angular momentum
	** vector) aligns with the positive z axis.  Only particle
	** positions are affected because the bulk spin should already be
	** removed: any non-zero particle velocities are assumed to
	** reflect the center-of-mass motion (which is NOT rotated), and
	** any particle spins are overwritten anyway (which is something
	** that may change in the future!).
	*/
	
	const Vector z = {0.0, 0.0, 1.0}; /* align to this axis */
	struct inSpinupAlign in;
	double v2;

	vectorCross(vVector, z, in.vAxis); /* get rotation axis */
	v2 = vectorMagSq(in.vAxis);
	if (v2 == 0.0)
		return; /* already aligned, or input vector has zero magnitude */
	vectorScale(in.vAxis, 1.0 / sqrt(v2), in.vAxis); /* unit vector */
	in.dAngle = acos(vVector[2] / vectorMag(vVector)); /* from dot product */
	if (in.dAngle < 0.0) /* shouldn't really happen -- could warn/assert */
		in.dAngle += M_PI; /* rotation angle between 0 and pi */
	vectorCopy(r_com, in.r_com);
	
	pstSpinupAlign(msr->pst, &in, sizeof(in), NULL, NULL);
	}

/*
** Routine to gradually increase the bulk rotation of all particles
** within a fixed distance of the origin at a specified rate.  We do
** this by treating the inner particles as an aggregate.  In addition,
** any particles between the inner and outer boundary are not spun up,
** and if present, prevent spinup of the inner particles.  Any
** particles beyond the outer boundary are not spun up either, but do
** not prevent spinup of the inner particles.  The idea is to
** temporarily suspend spinup when there are loose particles in the
** vicinity of the main body, but to also allow for the possibility of
** a satellite forming in orbit.  Currently we assume the main body
** will remain quite close (say within 1 particle radius) of the
** origin.  (If that's not the case, new code would be needed to
** identify the center of mass for the purpose of finding the
** particles belonging to the main body, distinct from those escaping
** and/or in orbit.)
*/

#ifndef SPINUP_WITH_AGGS

void msrSpinup(MSR msr, double dStep)
{
	static int bFirstCall = 1;
	static int bFirstInter = 1;
	static int nSuspendedSteps = 0; /*DEBUG this will not work with restarts!*/

	/* define following for convenience */
	static double dStart, dInter, dFinal, dInner, dOuter, dCoheCoe, dAngMomInc, dSpin_old;
	static int bMode, nStart, nInter, nFinal, nHold, nSpinup, nRestNum, nRestTime, nCoheTime;
	static int nRest = 1; /* record the current rest times */
	
	Aggregate agg;
	double dSpin;
	int nSteps;

	if (bFirstCall == 1) {
		/* convert parameters to pkdgrav units... */
		bMode = msr->param.bUseAngMomSpinup;
		dStart = msr->param.dSpinupPeriodStart * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		dInter = msr->param.dSpinupPeriodInter * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		dFinal = msr->param.dSpinupPeriodFinal * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		nStart = msr->param.nSpinupStepsStart;
		nInter = msr->param.nSpinupStepsInter;
		nFinal = msr->param.nSpinupStepsFinal;
		nHold = msr->param.nSpinupStepsHold;
		if (nFinal == 0)
			nFinal = msr->param.nSteps;
		nRestNum = msr->param.nSpinupStepsRestNum;
		nSpinup = nFinal/(nRestNum + 1);
		nRestTime = msr->param.nSpinupStepsRestTime;
		dInner = msr->param.dSpinupRadiusInner / 1.495978707e8; /* km -> AU */
		dOuter = msr->param.dSpinupRadiusOuter / 1.495978707e8; /* km -> AU */
#ifdef DEM_COHESION
		dCoheCoe = msr->param.dSpinupCoheCoeff;
		nCoheTime = msr->param.nSpinupCoheTime;
		if (nCoheTime == 0 || dCoheCoe == 0.0)
			nCoheTime = msr->param.nSteps;
#endif
		bFirstCall = 0;
		}

	nSteps = (int) dStep - nSuspendedSteps; /* number of steps so far */

#ifdef DEM_COHESION
	if((int) dStep > 0 && (int) dStep%nCoheTime == 0) {
		msr->param.FO.DP.dCohesiveCoeff += dCoheCoe;
		printf("%i %g\n",(int) dStep,msr->param.FO.DP.dCohesiveCoeff); /*DEBUG!!! YUN*/
		}
#endif

	/* stop spin-up during rest periods */
	if (nSteps >= nStart + nInter + nRest * nSpinup && nRest <= nRestNum) {
		nSuspendedSteps++;
		nRestTime--;
		if ( nRestTime >= 0)
			return;
		else {
			nRest++;
			nRestTime = msr->param.nSpinupStepsRestTime;
			}
		}

	/* stop spin up after rubble pile achieves final spin */
	if (nSteps >= nStart + nInter + nFinal + nHold)
		return;

	/* identify particles in the main body to make them into an aggregate */

	struct inSpinupGetAgg in_get;
	struct outSpinupGetAgg out_get;

	in_get.dRInner = dInner;
	in_get.dROuter = dOuter;
	pstSpinupGetAgg(msr->pst, &in_get, sizeof(in_get), &out_get, NULL);
	if (out_get.bSuspend == 1) {
		++nSuspendedSteps;/* at least one stray particle ==> no spinup (could add diagnostic message here) */
		return; 
		}
	
	/* get center of mass of aggregate (should be close to 0,0,0) */
	
	struct inAggsGetCOM in_com;
	struct outAggsGetCOM out_com;

	in_com.iAggIdx = -1; /* requires all particles' org idx = 0 */
	pstAggsGetCOM(msr->pst, &in_com, sizeof(in_com), &out_com, NULL);
	assert(out_com.m > 0.0);

	agg.mass = out_com.m;
	vectorScale(out_com.mr, 1.0 / out_com.m, agg.r_com);
	vectorScale(out_com.mv, 1.0 / out_com.m, agg.v_com);

	/* (could check to see how far r_com & v_com are from 0,0,0 here) */

	/* now get inertia tensor and angular momentum */
	
	struct inAggsGetAxesAndSpin in_axes;
	struct outAggsGetAxesAndSpin out_axes;

	in_axes.iAggIdx = -1;
	vectorCopy(agg.r_com, in_axes.r_com); /* could also copy from out_com */
	vectorCopy(agg.v_com, in_axes.v_com);
	pstAggsGetAxesAndSpin(msr->pst, &in_axes, sizeof(in_axes), &out_axes, NULL);

	/* fill out matrix (see msrAggsGetAxesAndSpin()) */

	out_axes.I[1][0] = out_axes.I[0][1];
	out_axes.I[2][0] = out_axes.I[0][2];
	out_axes.I[2][1] = out_axes.I[1][2];

	matrixCopy(out_axes.I, agg.lambda);
	matrixInvert(agg.lambda); /* lambda = inertia tensor inverse */
	matrixTransform(agg.lambda, out_axes.L, agg.omega); /* omega current spin */

	/* remove bulk spin component */

	struct inSpinup in_spin;

	vectorCopy(agg.r_com, in_spin.r_com);
	vectorScale(agg.omega, -1.0, in_spin.spin); /* for subtraction */
	pstSpinup(msr->pst, &in_spin, sizeof(in_spin), NULL, NULL);

	/* optionally re-align old aggregate spin/angmom vector with z-axis */

	if (msr->param.bSpinupAlignAxes != 0) {
		/*
		** NOTE: if spinup is suspended, this realignment will not
		** occur. If/when spinup resumes, there will likely be a
		** "glitch" as the body will have drifted away from
		** alignment. However, if we were to force realignment at all
		** times, (1) we would need to remove the aggregate spin and
		** put it back each time, even if no spinup is happening
		** (which is fine, but it requires altering the order of steps
		** in this routine) and, more seriously, (2) it interferes
		** with feedback from the nearby orbiting debris that may be
		** naturally causing the misalignment. Because of (2) for now
		** we suspend realignment when spinup is suspended, i.e., this
		** code segment is not reached because of the return earlier.
		*/
		if (bMode != 0) /* align with ang mom vector */
			msrSpinupAlignAxes(msr, agg.r_com, out_axes.L);
		else /* otherwise with spin vector */
			msrSpinupAlignAxes(msr, agg.r_com, agg.omega);
		}
	
	/* set new bulk spin */

	if (nSteps <= nStart)
		dSpin = dStart; /* keep forcing start spin */
	else if (nSteps <= nStart + nInter)
		dSpin = dStart + (nSteps - nStart) * (dInter - dStart) / nInter; /* interpolate between start & inter */
	else if (nSteps <= nStart + nInter + nFinal) {
		if (bMode == 0)
			dSpin = dInter + (nSteps - (nStart + nInter)) * (dFinal - dInter) / nFinal; /* interpolate between inter & final */
		else { /* using the angular momentum increment */
			if (msr->param.bSpinupAlignAxes) {
				/* unfortunately, need to recompute inertia tensor and inverse */
				pstAggsGetAxesAndSpin(msr->pst, &in_axes, sizeof(in_axes), &out_axes, NULL);
				out_axes.I[1][0] = out_axes.I[0][1];
				out_axes.I[2][0] = out_axes.I[0][2];
				out_axes.I[2][1] = out_axes.I[1][2];
				matrixCopy(out_axes.I, agg.lambda);
				matrixInvert(agg.lambda);
				}
			if (bFirstInter == 1) {
				dAngMomInc = out_axes.I[2][2] / nFinal * (1.0/dFinal - 1.0/dInter);
				bFirstInter = 0;
				dSpin_old = dInter;
				}
			//dSpin = 2.0 * M_PI / (agg.omega[2] + dAngMomInc * agg.lambda[2][2] ); /* depend on current spin state */
			dSpin = 1.0 / ( 1.0 / dSpin_old + dAngMomInc * agg.lambda[2][2] );
			dSpin_old = dSpin;
			//printf("%.10e %.10e %.10e %.10e %.10e\n",dSpin,dSpin_old,agg.lambda[2][2],out_axes.I[2][2],dAngMomInc); /*DEBUG!!! YUN*/
			}
		}
	else
		dSpin = dFinal; /* keep forcing final spin */
	dSpin = 2.0 * M_PI / dSpin; /* convert to angular frequency */
	vectorSet(in_spin.spin, 0.0, 0.0, dSpin); /* spin restricted to z direction */
	pstSpinup(msr->pst, &in_spin, sizeof(in_spin), NULL, NULL);
	}
#endif /* !SPINUP_WTIH_AGGS */

#ifdef SPINUP_WITH_AGGS
void msrSpinupAccel(MSR msr, double dStep)
{
	static int bFirstCall = 1;
	static int bFirstInter = 1;
	static int nSuspendedSteps = 0; /*DEBUG this will not work with restarts!*/

	/* define following for convenience */
	static double dStart, dInter, dFinal, dInner, dOuter, dCoheCoe, dAngMomInc, dSpin_old;
	static int bMode, nStart, nInter, nFinal, nHold, nSpinup, nRestNum, nRestTime, nCoheTime;
	static int nRest = 1; /* record the current rest times */
	
	Aggregate agg;
    Vector vOmegadot;
	double dSpin, dSpinstep;
	int nSteps;

	if (bFirstCall == 1) {
		/* convert parameters to pkdgrav units... */
		bMode = msr->param.bUseAngMomSpinup;
		dStart = msr->param.dSpinupPeriodStart * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		dInter = msr->param.dSpinupPeriodInter * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		dFinal = msr->param.dSpinupPeriodFinal * 2.0 * M_PI / (365.25 * 24.0); /* h -> yr/2pi */
		nStart = msr->param.nSpinupStepsStart;
		nInter = msr->param.nSpinupStepsInter;
		nFinal = msr->param.nSpinupStepsFinal;
		nHold = msr->param.nSpinupStepsHold;
		if (nFinal == 0)
			nFinal = msr->param.nSteps;
		nRestNum = msr->param.nSpinupStepsRestNum;
		nSpinup = nFinal/(nRestNum + 1);
		nRestTime = msr->param.nSpinupStepsRestTime;
		dInner = msr->param.dSpinupRadiusInner / 1.495978707e8; /* km -> AU */
		dOuter = msr->param.dSpinupRadiusOuter / 1.495978707e8; /* km -> AU */
#ifdef DEM_COHESION
		dCoheCoe = msr->param.dSpinupCoheCoeff;
		nCoheTime = msr->param.nSpinupCoheTime;
		if (nCoheTime == 0 || dCoheCoe == 0.0)
			nCoheTime = msr->param.nSteps;
#endif
		bFirstCall = 0;
		}

	nSteps = (int) dStep - nSuspendedSteps; /* number of steps so far */

#ifdef DEM_COHESION
	if((int) dStep > 0 && (int) dStep%nCoheTime == 0) {
		msr->param.FO.DP.dCohesiveCoeff += dCoheCoe;
		printf("%i %g\n",(int) dStep,msr->param.FO.DP.dCohesiveCoeff); /*DEBUG!!! YUN*/
		}
#endif

	/* stop spin-up during rest periods */
	if (nSteps >= nStart + nInter + nRest * nSpinup && nRest <= nRestNum) {
		nSuspendedSteps++;
		nRestTime--;
		if ( nRestTime >= 0)
			return;
		else {
			nRest++;
			nRestTime = msr->param.nSpinupStepsRestTime;
			}
		}

	/* stop spin up after rubble pile achieves final spin */
	if (nSteps >= nStart + nInter + nFinal + nHold) 
		return;

	/* identify particles in the main body to make them into an aggregate */

	struct inSpinupGetAgg in_get;
	struct outSpinupGetAgg out_get;

	in_get.dRInner = dInner;
	in_get.dROuter = dOuter;
	pstSpinupGetAgg(msr->pst, &in_get, sizeof(in_get), &out_get, NULL);
	if (out_get.bSuspend == 1) {
		++nSuspendedSteps; /* at least one stray particle ==> no spinup (could add diagnostic message here) */
        pstSpinupResetAggs(msr->pst, NULL, 0, NULL, NULL);
		return; 
		}
	
	/* get center of mass of aggregate (should be close to 0,0,0) */
	
	struct inAggsGetCOM in_com;
	struct outAggsGetCOM out_com;

	in_com.iAggIdx = INT_MAX; /* requires all particles' org idx = -1 - INT_MAX */
	pstAggsGetCOM(msr->pst, &in_com, sizeof(in_com), &out_com, NULL);
	assert(out_com.m > 0.0);

	agg.mass = out_com.m;
	vectorScale(out_com.mr, 1.0 / out_com.m, agg.r_com);
	vectorScale(out_com.mv, 1.0 / out_com.m, agg.v_com);
    
	/* (could check to see how far r_com & v_com are from 0,0,0 here) */

    /* Loop over aggregates to remove any drift in the rubble pile COM */

    if (msr->param.bSpinupAlignAxes != 0) {
        int i;
        Aggregate *trueagg;

        for (i=0;i<msr->nAggs;i++) {
            trueagg = &msr->pAggs[i];
            if (trueagg->bAssigned) {
                vectorSub(trueagg->r_com, agg.r_com, trueagg->r_com);
                vectorSub(trueagg->v_com, agg.v_com, trueagg->v_com);
            }
        }
    } 

	/* now get inertia tensor and angular momentum */
	
	struct inAggsGetAxesAndSpin in_axes;
	struct outAggsGetAxesAndSpin out_axes;

	in_axes.iAggIdx = INT_MAX;
	vectorCopy(agg.r_com, in_axes.r_com); /* could also copy from out_com */
	vectorCopy(agg.v_com, in_axes.v_com);
	pstAggsGetAxesAndSpin(msr->pst, &in_axes, sizeof(in_axes), &out_axes, NULL);
    
	/* fill out matrix (see msrAggsGetAxesAndSpin()) */

	out_axes.I[1][0] = out_axes.I[0][1];
	out_axes.I[2][0] = out_axes.I[0][2];
	out_axes.I[2][1] = out_axes.I[1][2];

    matrixCopy(out_axes.I, agg.lambda);
	matrixInvert(agg.lambda); /* lambda = inertia tensor inverse */
	matrixTransform(agg.lambda, out_axes.L, agg.omega); /* omega current spin */
    
	struct inSpinupAccel in_spinaccel;

	vectorCopy(agg.r_com, in_spinaccel.r_com);
    in_spinaccel.bSpinupAlignAxes = msr->param.bSpinupAlignAxes;

	/* Calculate new bulk spin and omegadot required to reach it */

	if (nSteps <= nStart) 
		dSpin = dStart; /* keep forcing start spin */
	else if (nSteps <= nStart + nInter) 
		dSpin = dStart + (nSteps - nStart) * (dInter - dStart) / nInter; /* interpolate between start & inter */
	else if (nSteps <= nStart + nInter + nFinal) {
		if (bMode == 0)
			dSpin = dInter + (nSteps - (nStart + nInter)) * (dFinal - dInter) / nFinal; /* interpolate between inter & final */
		else { /* using the angular momentum increment */
			if (bFirstInter == 1) {
				dAngMomInc = out_axes.I[2][2] / nFinal * (1.0/dFinal - 1.0/dInter);
				bFirstInter = 0;
				dSpin_old = dInter;
				}
			//dSpin = 2.0 * M_PI / (agg.omega[2] + dAngMomInc * agg.lambda[2][2] ); /* depend on current spin state */
			dSpin = 1.0 / ( 1.0 / dSpin_old + dAngMomInc * agg.lambda[2][2] );
			dSpin_old = dSpin;
			//printf("%.10e %.10e %.10e %.10e %.10e\n",dSpin,dSpin_old,agg.lambda[2][2],out_axes.I[2][2],dAngMomInc); /*DEBUG!!! YUN*/
			}
		}
	else
		dSpin = dFinal; /* keep forcing final spin */
	dSpin = 2.0 * M_PI / dSpin; /* convert to angular frequency */

    if (msr->param.bSpinupAlignAxes != 0) {
        vOmegadot[0] = -agg.omega[0] / msr->param.dDelta; /* If desired, set rate of change of x and y-components of spin needed to zero them. */
        vOmegadot[1] = -agg.omega[1] / msr->param.dDelta;
    } else
        vOmegadot[0] = vOmegadot[1] = 0.0;
   
    vOmegadot[2] = (dSpin - agg.omega[2]) / msr->param.dDelta; /* Set rate of change of z-component of spin to achieve desired speed. */
	vectorSet(in_spinaccel.omegadot, vOmegadot[0], vOmegadot[1], vOmegadot[2]);
	pstSpinupAccel(msr->pst, &in_spinaccel, sizeof(in_spinaccel), NULL, NULL);
    pstSpinupResetAggs(msr->pst, NULL, 0, NULL, NULL);
	}

#endif /* SPINUP_WITH_AGGS */

#endif /* SPINUP */

#endif /* COLLISIONS */

#ifdef AGGS

#ifdef AGGS_IN_PATCH

void msrAggsInPatchGetUnwrapped(MSR msr,int iAggIdx,double dTime)
{
	struct inAggsInPatchGetRef inGetRef;
	struct outAggsInPatchGetRef outGetRef;
	struct inAggsInPatchGetUnwrapped inGetUnwrapped;

	inGetRef.iAggIdx = inGetUnwrapped.iAggIdx = iAggIdx;
	pstAggsInPatchGetRef(msr->pst,&inGetRef,sizeof(inGetRef),&outGetRef,NULL);
	inGetUnwrapped.dTime = dTime;
	inGetUnwrapped.PP = msr->param.PP; /* struct copy */
	vectorCopy(outGetRef.r_max,inGetUnwrapped.r_ref);
	pstAggsInPatchGetUnwrapped(msr->pst,&inGetUnwrapped,sizeof(inGetUnwrapped),NULL,NULL);
	}

void msrAggsInPatchCheckCOM(MSR msr,int iAggIdx,Aggregate *agg,double dTime)
{
	/*
	** The purpose of this function is to ensure the aggregate
	** center-of-mass position is located within the patch boundaries
	** (the center-of-mass velocity may need to be adjusted as well).
	** This is needed in order to properly compute the torque on the
	** aggregate, although strictly speaking it is probably ok if the
	** aggregate wanders a bit out of the patch, so long as a
	** sufficient number of ghost cells are being used.  It is
	** necessary to update the unwrapped constituent particle
	** positions and velocities as well.  This routine need only be
	** called after the center-of-mass position of the aggregate has
	** been computed or adjusted, presently only done in
	** msrAggsUpdate() and msrAggsAdvance().
	*/

	const PATCH_PARAMS *PP = &msr->param.PP;

	struct inAggsInPatchOffset in;
	FLOAT dx=0.0,dy=0.0,dvy=0.0;

	if (agg->r_com[0] > 0.5*PP->dWidth) {
		agg->r_com[0] += (dx = -PP->dWidth);
		agg->r_com[1] += (dy = SHEAR(-1,dTime,PP));
		agg->v_com[1] += (dvy = 1.5*PP->dOrbFreq*PP->dWidth);
		}
	else if (agg->r_com[0] < -0.5*PP->dWidth) {
		agg->r_com[0] += (dx = PP->dWidth);
		agg->r_com[1] += (dy = SHEAR(1,dTime,PP));
		agg->v_com[1] += (dvy = -1.5*PP->dOrbFreq*PP->dWidth);
		}
	if (agg->r_com[1] > 0.5*PP->dLength) {
		agg->r_com[1] -= PP->dLength;
		dy -= PP->dLength;
		}
	else if (agg->r_com[1] < -0.5*PP->dLength) {
		agg->r_com[1] += PP->dLength;
		dy += PP->dLength;
		}

	if (dx != 0.0 || dy != 0.0) {
		assert(!PP->bRandAzWrap); /* can't handle aggs random wrapping! */
		in.iAggIdx = iAggIdx;
		in.dx = dx;
		in.dy = dy;
		in.dvy = dvy;
		in.bDoUnwrapped = 1;
		pstAggsInPatchOffset(msr->pst,&in,sizeof(in),NULL,NULL);
		}

	/* update momentum for symplectic patch integrator, radial wrap only (Cf. Randall's e-mail 7/30/08) */
	agg->dPy_com -= dvy/3.0;
	}

void msrAggsInPatchRotateAxes(MSR msr,Aggregate *agg,double dt)
{
	/* 
	** Rotates aggregate orientation through angle Omega * dt.  Needed
	** for resyncing between space & patch frames.  Used in
	** msrAggsAdvance() and msrAggsBackDrift().
	*/

	Matrix mPatchRot,mTmp;
	double alpha,c,s;

	assert(agg != NULL);
	if (dt == 0.0)
		return;
	alpha = - msr->param.PP.dOrbFreq*dt;
	c = cos(alpha);
	s = sin(alpha);
	matrixIdentity(mPatchRot);
	mPatchRot[0][0] = c;
	mPatchRot[0][1] = -s;
	mPatchRot[1][0] = s;
	mPatchRot[1][1] = c;
	matrixMultiply(mPatchRot,agg->lambda,mTmp);
	matrixCopy(mTmp,agg->lambda);
	}

#endif /* AGGS_IN_PATCH */

void msrAggsSetBodyPos(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* to be called ONLY by msrAggsUpdate() */
	/* requires call to msrAggsGetAxesAndSpin() first */

	struct inAggsSetBodyPos in;
	in.iAggIdx = iAggIdx;
	matrixTranspose(agg->lambda,in.spaceToBody);
	pstAggsSetBodyPos(msr->pst,&in,sizeof(in),NULL,NULL);
	}

#define EPS 1.0e-12 /* ratios below this treated as zero in inertia tensor */
/*#define VERBOSE*/ /* uncomment to display inertia tensor & diagonalization */

void msrAggsGetAxesAndSpin(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* to be called ONLY by msrAggsUpdate() */
	/* requires call to msrAggsGetCOM() first */

	struct inAggsGetAxesAndSpin in;
	struct outAggsGetAxesAndSpin out;
	Matrix spaceToBody;
	Vector h;
	int k;
	/* compute inertia tensor I and angular momentum L (not spin, yet!) */
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);
	vectorCopy(agg->v_com,in.v_com);
#ifdef AGGS_IN_PATCH
	in.PP = msr->param.PP; /* struct copy */
#endif
	pstAggsGetAxesAndSpin(msr->pst,&in,sizeof(in),&out,NULL);

	/*
	** The inertia tensor computed is only valid for the diagonal and
	** upper triangle.  Use symmetry to fill it out for diagonalization.
	*/
	out.I[1][0] = out.I[0][1]; /* row-column notation */
	out.I[2][0] = out.I[0][2];
	out.I[2][1] = out.I[1][2];


	{
	    /* zero small inertia tensor elements */
	    /* (helps with integration stability) */

	    double fMax = 0.0;
	    int i, j;
	    for (i = 0; i < 3; i++)
		for (j = i; j < 3; j++)
		    if (fabs(out.I[i][j]) > fMax)
			fMax = fabs(out.I[i][j]);
	    assert(fMax > 0.0);
	    fMax = 1.0 / fMax;
	    for (i = 0; i < 3; i++)
		for (j = i; j < 3; j++)
		    if (fabs(out.I[i][j]) * fMax < EPS)
			out.I[i][j] = out.I[j][i] = 0.0;
	    }


#ifdef VERBOSE
	printf("Aggregate %i\n", iAggIdx);
	printf(" Inertia tensor:\n");
	printf("  Ixx = %23.16e Ixy = %23.16e Ixz = %23.16e\n",
	       out.I[0][0], out.I[0][1], out.I[0][2]);
	printf("  Iyx = %23.16e Iyy = %23.16e Iyz = %23.16e\n",
	       out.I[1][0], out.I[1][1], out.I[1][2]);
	printf("  Izx = %23.16e Izy = %23.16e Izz = %23.16e\n",
	       out.I[2][0], out.I[2][1], out.I[2][2]);
#endif /* VERBOSE */

#ifdef AGGS_IN_PATCH
	/*
	** Need to add the component due to frame rotation for computing
	** the body spin.  This gives the space-frame spin, which is
	** transformed to the body frame next.
	*/
	out.L[0] += msr->param.PP.dOrbFreq*out.I[0][2];
	out.L[1] += msr->param.PP.dOrbFreq*out.I[1][2];
	out.L[2] += msr->param.PP.dOrbFreq*out.I[2][2];
#endif

	/* compute moments and principal axes (lambda) */
	matrixDiagonalize(out.I,agg->moments,agg->lambda,0/* no sort */);

#ifdef VERBOSE
	printf(" Principal moments:\n");
	printf("  p1m = %.16e p2m = %.16e p3m = %.16e\n",
	       agg->moments[0], agg->moments[1], agg->moments[2]);
	printf(" Principal axes:\n");
	printf("  p1x = %23.16e p1y = %23.16e p1z = %23.16e\n",
	       agg->lambda[0][0], agg->lambda[1][0], agg->lambda[2][0]);
	printf("  p2x = %23.16e p2y = %23.16e p2z = %23.16e\n",
	       agg->lambda[0][1], agg->lambda[1][1], agg->lambda[2][1]);
	printf("  p3x = %23.16e p3y = %23.16e p3z = %23.16e\n",
	       agg->lambda[0][2], agg->lambda[1][2], agg->lambda[2][2]);
	printf(" Principal axis magnitudes (should be 1):\n");
	printf("  p1 = %.16e\n", sqrt(agg->lambda[0][0] * agg->lambda[0][0] +
				      agg->lambda[1][0] * agg->lambda[1][0] +
				      agg->lambda[2][0] * agg->lambda[2][0]));
	printf("  p2 = %.16e\n", sqrt(agg->lambda[0][1] * agg->lambda[0][1] +
				      agg->lambda[1][1] * agg->lambda[1][1] +
				      agg->lambda[2][1] * agg->lambda[2][1]));
	printf("  p3 = %.16e\n", sqrt(agg->lambda[0][2] * agg->lambda[0][2] +
				      agg->lambda[1][2] * agg->lambda[1][2] +
				      agg->lambda[2][2] * agg->lambda[2][2]));
	printf(" Dot products (should be 0):\n");
	printf("  p1 . p2 = %.16e\n",
	       agg->lambda[0][0] * agg->lambda[0][1] +
	       agg->lambda[1][0] * agg->lambda[1][1] +
	       agg->lambda[2][0] * agg->lambda[2][1]);
	printf("  p1 . p3 = %.16e\n",
	       agg->lambda[0][0] * agg->lambda[0][2] +
	       agg->lambda[1][0] * agg->lambda[1][2] +
	       agg->lambda[2][0] * agg->lambda[2][2]);
	printf("  p2 . p3 = %.16e\n",
	       agg->lambda[0][1] * agg->lambda[0][2] +
	       agg->lambda[1][1] * agg->lambda[1][2] +
	       agg->lambda[2][1] * agg->lambda[2][2]);
	{
	    Vector dVec1, dVec2;
	    printf(" Eigensystem check:\n");
	    vectorSet(dVec1, agg->lambda[0][0], agg->lambda[1][0], agg->lambda[2][0]);
	    matrixTransform(out.I, dVec1, dVec2);
	    vectorScale(dVec2, 1.0 / agg->moments[0], dVec1);
	    printf("  P p1 / p1m = %23.16e %23.16e %23.16e (compare p1x,y,z)\n",
		   dVec1[0], dVec1[1], dVec1[2]);
	    vectorSet(dVec1, agg->lambda[0][1], agg->lambda[1][1], agg->lambda[2][1]);
	    matrixTransform(out.I, dVec1, dVec2);
	    vectorScale(dVec2, 1.0 / agg->moments[1], dVec1);
	    printf("  P p2 / p2m = %23.16e %23.16e %23.16e (compare p2x,y,z)\n",
		   dVec1[0], dVec1[1], dVec1[2]);
	    vectorSet(dVec1, agg->lambda[0][2], agg->lambda[1][2], agg->lambda[2][2]);
	    matrixTransform(out.I, dVec1, dVec2);
	    vectorScale(dVec2, 1.0 / agg->moments[2], dVec1);
	    printf("  P p3 / p3m = %23.16e %23.16e %23.16e (compare p3x,y,z)\n",
		   dVec1[0], dVec1[1], dVec1[2]);
	    }
#endif /* VERBOSE */

	/* transform L to body frame (not stored) */
	matrixTranspose(agg->lambda,spaceToBody);
	matrixTransform(spaceToBody,out.L,h);

	/* now compute spin vector omega = I^{-1} h */
	for (k=0;k<3;k++) {/* shortcut for diagonal matrix in body frame */
		assert(agg->moments[k] > 0.0);
		agg->omega[k] = h[k]/agg->moments[k];
		}
	/* maximum "radius" based on maximum particle distance from center */
	assert(out.rad_max_sq > 0.0);
	if (msr->param.bPeriodic) {	/* No aggs w/ radius > 0.5 * boundary length allowed */
		assert(out.rad_max_sq < msr->param.dxPeriod * msr->param.dxPeriod);
		assert(out.rad_max_sq < msr->param.dyPeriod * msr->param.dyPeriod);
		assert(out.rad_max_sq < msr->param.dzPeriod * msr->param.dzPeriod);
		}
	agg->rad_max = sqrt(out.rad_max_sq); /* stored as square during calc */

	/* effective radius based on volume of constituent particles */
	assert(out.volume > 0.0); /* note: volume does not contain 4/3 PI factor */
	agg->rad_eff = pow(2.0*out.volume,1.0/3); /* 2.0 = packing fudge factor */

#ifdef AGGS_IN_PATCH
	/* warn if agg is getting comparable to patch size! */
	if (agg->rad_max > 0.2*msr->param.PP.dWidth ||
		agg->rad_max > 0.2*msr->param.PP.dLength) {
#if (INTERNAL_WARNINGS != 0)
		static int nWarn = 1;
		if (nWarn == 1 || nWarn%INTERNAL_WARNINGS == 0)
			(void) fprintf(stderr,"AGGS_IN_PATCH WARNING #%i: iAggIdx=%i, rad_max=%g, patch dimens: (%g, %g)\n",
				       nWarn,iAggIdx,agg->rad_max,
				       msr->param.PP.dWidth,
				       msr->param.PP.dLength);
		++nWarn;
		//		assert(agg->rad_max < 1.0*msr->param.PP.dWidth && agg->rad_max < 1.0*msr->param.PP.dLength); /* RP-DEBUG 7-13-09 */
#endif /* INTERNAL_WARNINGS */
				}
#endif /* AGGS_IN_PATCH */
	}

#undef VERBOSE
#undef EPS

#ifdef FAST_AGGS
void msrAggsInitUnwrapped(MSR msr,int iAggIdx,Aggregate *agg)	// JVD
{
	/* Description */

	struct inAggsGetRef in_ref;
	struct outAggsGetRef out_ref;
	struct inAggsGetPeriodicCOM in_COM;
	struct outAggsGetPeriodicCOM out_COM;
	struct inAggsSetUnwrapped in_unwrapped;
	FLOAT fPeriod[3];
	int i;	// Iterable

	in_ref.iAggIdx = in_COM.iAggIdx = in_unwrapped.iAggIdx = iAggIdx;

	/* Initialize periodic boundary lengths */
	fPeriod[0] = msr->param.dxPeriod;
	fPeriod[1] = msr->param.dyPeriod;
	fPeriod[2] = msr->param.dzPeriod;


	/* Get reference vector -- Position of any particle in agg */
	pstAggsGetRef(msr->pst,&in_ref,sizeof(in_ref),&out_ref,NULL);
	vectorCopy(out_ref.r_ref,in_COM.r_ref);

	/* Calculate COM using the reference particle as the agg system origin */
	pstAggsGetPeriodicCOM(msr->pst,&in_COM,sizeof(in_COM),&out_COM,NULL);

	assert(out_COM.m > 0.0);
	agg->mass = out_COM.m;
	vectorScale(out_COM.mr,1.0/agg->mass,agg->r_com);
	vectorScale(out_COM.mv,1.0/agg->mass,agg->v_com);

	/* Fix agg->r_com in case it lies outside the periodic boundaries */
	for (i=0;i<3;i++) {
		while (agg->r_com[i] >= msr->fCenter[i] + 0.5 * fPeriod[i])
			agg->r_com[i] -= fPeriod[i];
		while (agg->r_com[i] < msr->fCenter[i] - 0.5 * fPeriod[i])
			agg->r_com[i] += fPeriod[i];
		}

	/* Set the appropriate r_unwrapped for every particle in the simulation */
	vectorCopy(agg->r_com,in_unwrapped.r_com);

	pstAggsSetUnwrapped(msr->pst,&in_unwrapped,sizeof(in_unwrapped),NULL,NULL);

	}
#endif /*FAST_AGGS*/

void msrAggsGetCOM(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* to be called ONLY by msrAggsUpdate() */

	struct inAggsGetCOM in;
	struct outAggsGetCOM out;
	in.iAggIdx = iAggIdx;
#ifdef AGGS_IN_PATCH
	in.PP = msr->param.PP; /* struct copy */
#endif
	pstAggsGetCOM(msr->pst,&in,sizeof(in),&out,NULL);
	assert(out.m > 0.0);

	agg->mass = out.m;
	vectorScale(out.mr,1.0/out.m,agg->r_com);
	vectorScale(out.mv,1.0/out.m,agg->v_com);

	/* 
	   RP-DEBUG: why don't we check the COM and 
	   adjust if outside the patch..? 
	   We'd also need to adjust the unwrapped positions
	   and velocities...
	*/

	}
/*
void msrAggsGetPeriodicCOM(MSR msr,int iAggIdx,Aggregate *agg)	// JVD
{

	struct inAggsGetPeriodicCOM in;
	struct outAggsGetPeriodicCOM out;
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);

	pstAggsGetPeriodicCOM(msr->pst,&in,sizeof(in),&out,NULL);

	vectorScale(out.mr,1.0/agg->mass,agg->r_com);

	}
*/
void msrAggsCountPart(MSR msr,int iAggIdx,Aggregate *agg,int *nPart)
{
	/* counts particles (if any) in agg */
	struct inAggsCountPart in;
	struct outAggsCountPart out;

	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	in.iAggIdx = iAggIdx;
	/* note: agg may not be assigned -- that's ok */
	pstAggsCountPart(msr->pst,&in,sizeof(in),&out,NULL);
	*nPart = out.nPart;
	}

void msrAggsDelete(MSR msr,int iAggIdx,Aggregate *agg,double dEventTime)
{
	struct inAggsDelete in;
	struct outAggsDelete out;

	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	in.iAggIdx = iAggIdx;
#ifdef AGGS_IN_PATCH
	in.dStepTime = msr->dTime;
	in.dEventTime = dEventTime;
	in.PP = msr->param.PP; /* struct copy */
#endif
	pstAggsDelete(msr->pst,&in,sizeof(in),&out,NULL);
	assert(out.bFound == 1);
	agg->bAssigned = 0;
	if (msr->param.bVDetails)
		(void) printf("Aggregate %i deleted.\n",iAggIdx);
	}

#ifndef DEM
void msrAggsGetAccel(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* particle accelerations must be up to date */
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	vectorZero(agg->a_com);
	if (msr->param.bDoGravity)
	{
		struct inAggsGetAccel in;
		struct outAggsGetAccel out;
		/* get acceleration of center of mass */
		in.iAggIdx = iAggIdx;
		pstAggsGetAccel(msr->pst,&in,sizeof(in),&out,NULL);
		assert(out.m > 0.0);
		if (fabs(out.m - agg->mass) > CONSERVE_FRAC*agg->mass) {
			(void) printf("WARNING: Aggregate %i mass not conserved: %.16e != %.16e!\n",iAggIdx,out.m,agg->mass);
			}
		vectorScale(out.ma,1.0/out.m,agg->a_com);
		}
	}


void msrAggsGetTorque(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* aggregate c-o-m pos, acceleration, & body orientation assumed up to date */

	Matrix spaceToBody;
	Vector vTorque;
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	vectorZero(vTorque);
	if (msr->param.bDoSelfGravity)
	{
		/* get torque around center of mass */
		struct inAggsGetTorque in;
		struct outAggsGetTorque out;
		in.iAggIdx = iAggIdx;
		vectorCopy(agg->r_com,in.r_com);
		vectorCopy(agg->a_com,in.a_com);
		pstAggsGetTorque(msr->pst,&in,sizeof(in),&out,NULL);
		vectorCopy(out.torque,vTorque);
		}
	matrixTranspose(agg->lambda,spaceToBody); /* this is also the inverse, since lambda is orthogonal */
#ifdef AGGS_IN_PATCH
	{
		/*
		** We have to un-diagonalize the inertia tensor in order to
		** compute tides.  For that, we need the inverse of lambda
		** (spaceToBody).
		*/
		double preFac,y_component,z_component;
		/*
		** The following isolates the two terms in the patch-frame
		** inertia tensor that we need to compute the tides.
		*/
		preFac = 3.0*msr->param.PP.dOrbFreq*msr->param.PP.dOrbFreq;
		y_component = - preFac*(agg->lambda[2][0]*spaceToBody[0][0]*agg->moments[0] + agg->lambda[2][1]*spaceToBody[1][0]*agg->moments[1] + agg->lambda[2][2]*spaceToBody[2][0]*agg->moments[2]);
		z_component = preFac*(agg->lambda[1][0]*spaceToBody[0][0]*agg->moments[0] + agg->lambda[1][1]*spaceToBody[1][0]*agg->moments[1] + agg->lambda[1][2]*spaceToBody[2][0]*agg->moments[2]);
		/* actually add the tides... */
		vTorque[1] += y_component;
		vTorque[2] += z_component;
		}
#endif /* AGGS_IN_PATCH */
	/* transform torque vector to body frame */
	matrixTransform(spaceToBody,vTorque,agg->torque);
	}
#endif /* !DEM */

#ifdef DEM	/* Don't need accelerations as inputs to torque calc w/ DEM */
void msrAggsGetAccelAndTorqueDEM(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* Particle & agg accels, agg COM pos, & agg body orientation must be up to date */
	Matrix spaceToBody;
	Vector vTorque;
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	vectorZero(agg->a_com);
	vectorZero(vTorque);

	/* Fill input structs */
	struct inAggsGetAccelAndTorqueDEM in;
	struct outAggsGetAccelAndTorqueDEM out;
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);
	pstAggsGetAccelAndTorqueDEM(msr->pst,&in,sizeof(in),&out,NULL);

	/* Check outputs */
	assert(out.m > 0.0);
	if (fabs(out.m - agg->mass) > CONSERVE_FRAC*agg->mass) {
		(void) printf("WARNING: Aggregate %i mass not conserved: %.16e != %.16e!\n",iAggIdx,out.m,agg->mass);
		}

	/* Fill output structs */
	vectorScale(out.ma,1.0/out.m,agg->a_com);
	vectorCopy(out.torque,vTorque);

	/* Transform torque vector to body frame (after AGGS_IN_PATCH) */
	matrixTranspose(agg->lambda,spaceToBody);
#ifdef AGGS_IN_PATCH
	{
		/*
		** We have to un-diagonalize the inertia tensor in order to
		** compute tides.  For that, we need the inverse of lambda
		** (spaceToBody).
		*/
		double preFac,y_component,z_component;
		/*
		** The following isolates the two terms in the patch-frame
		** inertia tensor that we need to compute the tides.
		*/
		preFac = 3.0*msr->param.PP.dOrbFreq*msr->param.PP.dOrbFreq;
		y_component = - preFac*(agg->lambda[2][0]*spaceToBody[0][0]*agg->moments[0] + agg->lambda[2][1]*spaceToBody[1][0]*agg->moments[1] + agg->lambda[2][2]*spaceToBody[2][0]*agg->moments[2]);
		z_component = preFac*(agg->lambda[1][0]*spaceToBody[0][0]*agg->moments[0] + agg->lambda[1][1]*spaceToBody[1][0]*agg->moments[1] + agg->lambda[1][2]*spaceToBody[2][0]*agg->moments[2]);
		/* actually add the tides... */
		vTorque[1] += y_component;
		vTorque[2] += z_component;
		}
#endif /* AGGS_IN_PATCH */
	matrixTransform(spaceToBody,vTorque,agg->torque);
	}
#endif /* DEM */

void msrAggsUpdate(MSR msr,int iAggIdx,Aggregate *agg,int bDoAccel,double dUpdateTime)
{
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
#ifdef AGGS_IN_PATCH
	/* here dUpdateTime is measured from start of step; msr->dTime is set in msrGravity() at the start of the step */
	msrAggsInPatchGetUnwrapped(msr,iAggIdx,msr->dTime + dUpdateTime);
#endif
#ifdef FAST_AGGS
	/* get center of mass info, needed before axes */
	if (msr->param.bPeriodic)	// To properly initalize COM w/ Aggs + PBCs
		msrAggsInitUnwrapped(msr,iAggIdx,agg);
	else
		msrAggsGetCOM(msr,iAggIdx,agg);
#else /*!FAST_AGGS*/
	msrAggsGetCOM(msr,iAggIdx,agg);
#endif /*!FAST_AGGS*/

#ifdef AGGS_IN_PATCH
	/* recalculate Py using new velocity and position info (can be
	   outside the patch; ...CheckCOM() fixes this */
	agg->dPy_com = agg->v_com[1] + 2.0*msr->param.PP.dOrbFreq*
	  (agg->r_com[0] + agg->v_com[0]*((msr->param.dDelta/2.0) - dUpdateTime)); /*RP-DEBUG-dPy*/
	/* don't do this in the middle of a step; 
	   don't want particle positions adjusted without cause */
	if (dUpdateTime == 0.0)
		msrAggsInPatchCheckCOM(msr,iAggIdx,agg,msr->dTime);
#endif
#ifdef DEM
	{
		struct inAggsSetMassDEM in;
		in.iAggIdx = iAggIdx;
		in.fMass = agg->mass;
		pstAggsSetMassDEM(msr->pst,&in,sizeof(in),NULL,NULL);
		}
#endif
	/* get moments, principal axes, and spin vector */
	msrAggsGetAxesAndSpin(msr,iAggIdx,agg);
	/* transform particle positions to body frame */
	msrAggsSetBodyPos(msr,iAggIdx,agg);
	/*
	 ** NOTE: particle velocities are not updated at this point since
	 ** this is done during the kick step (or back drift after a
	 ** collision)---this works only for the KDK scheme!  Similarly,
	 ** particle spins are not updated since they're currently only
	 ** needed just before a collision, or at the end of the step (for
	 ** output purposes).  Hence the spins are updated in
	 ** msrAggsAdvance().
	 */
	if (bDoAccel) {
		/* assumes particle acclerations are up to date! */
#ifdef DEM
		msrAggsGetAccelAndTorqueDEM(msr,iAggIdx,agg);
#else /* !DEM */
		msrAggsGetAccel(msr,iAggIdx,agg);
		msrAggsGetTorque(msr,iAggIdx,agg);
#endif /* !DEM */
#ifdef AGGS_IN_PATCH
		/* add Hill's terms -- see msrAggsGravity() for details */
		agg->a_com[0] -= msr->param.PP.dOrbFreq*msr->param.PP.dOrbFreq*agg->r_com[0];
		agg->a_com[2] -= msr->param.PP.dOrbFreqZ2*agg->r_com[2];
#endif /* AGGS_IN_PATCH */
		}
	/* set update clock */
	agg->dLastUpdate = dUpdateTime;
	}

#ifndef DEM
void msrAggsCheckStress(MSR msr,int iAggIdx,Aggregate *agg,int *bLostMass)
{
	struct inAggsCheckStress in;
	struct outAggsCheckStress out;

	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	assert(bLostMass != NULL);
	*bLostMass = 0;
	if (msr->param.CP.SP.dTensileCoef == DBL_MAX && msr->param.CP.SP.dShearCoef == DBL_MAX)
		return; /* no point in doing this check if aggs are perfectly rigid */
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);
	vectorCopy(agg->a_com,in.a_com);
	matrixTransform(agg->lambda,agg->omega,in.omega); /* space frame */
	in.SP = msr->param.CP.SP; /* strength parameters; struct copy */
	in.mass = agg->mass;
	in.rad_max = agg->rad_max;
	in.rad_eff = agg->rad_eff;
#ifdef AGGS_IN_PATCH
	in.dTime = msr->dTime;
	in.PP = msr->param.PP; /* struct copy */
#endif
	pstAggsCheckStress(msr->pst,&in,sizeof(in),&out,NULL);
	assert(out.nLost >= 0 && out.nLeft >= 0);
	if (out.nLost > 0) {
		int nPart;

		msrAggsCountPart(msr,iAggIdx,agg,&nPart); /* safety check */
		assert(nPart == out.nLeft);

		if (msr->param.bVDetails)
			(void) printf("Aggregate %i lost %i particle(s), %i left.\n",iAggIdx,out.nLost,out.nLeft);
		if (nPart > 1)
			msrAggsUpdate(msr,iAggIdx,agg,0/*accel already up to date*/,msr->param.dDelta/*RP: END of step*/);
		else if (nPart == 1)
			msrAggsDelete(msr,iAggIdx,agg,msr->param.dDelta/*RP: END of step*/);
		else { /* (nPart == 0) */
			agg->bAssigned = 0;
			if (msr->param.bVDetails)
				(void) printf("Aggregate %i deleted.\n",iAggIdx);
			}
		*bLostMass = 1;
		}
	}
#endif /* !DEM */

void msrAggsGravity(MSR msr)
{
	/* to be called after particle accelerations have been computed */

	Aggregate *agg;
	int i,bLostMass;

#ifdef FAST_AGGS
	msrAggsReorder(msr);	/* Reorder for FAST_AGGS */
#endif /* FAST_AGGS */
	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned)
#ifdef DEM
			msrAggsGetAccelAndTorqueDEM(msr,i,agg);
#else /* !DEM */
			do {
				msrAggsGetAccel(msr,i,agg);
				msrAggsCheckStress(msr,i,agg,&bLostMass);
				if (!agg->bAssigned)
					break; /* this means aggregate fragmented entirely */
				} while (bLostMass == 1);
		if (agg->bAssigned)
			msrAggsGetTorque(msr,i,agg);
#endif /* !DEM */
#ifdef AGGS_IN_PATCH
		/* finally, add the Hill's terms, if applicable */
		if (agg->bAssigned && msr->param.PP.bPatch) {
			/*
			** Following is from Dan Scheere's derivation of rotating
			** frame terms:
			**
			** 3 n^2 (R-hat dot r) R-hat - n^2 (z dot r) z - 2 n z cross r-dot
			** 
			** where n = mean motion (= sqrt(G M_planet/R^3)), R = orbital
			** radius of patch, R-hat = (1,0,0) in our system, r =
			** center-of-mass position of aggregate in patch frame, and z
			** = (0,0,1) in our system.  BUT, to do this symplectically,
			** we use Tom's method of separating the velocity-dependent
			** terms.
			*/
			agg->a_com[0] -= msr->param.PP.dOrbFreq*msr->param.PP.dOrbFreq*agg->r_com[0];
			agg->a_com[2] -= msr->param.PP.dOrbFreqZ2*agg->r_com[2];
			}
#endif /* AGGS_IN_PATCH */
#ifdef DEM
		if (msr->param.FO.DP.dAccCrit > 0.0) {
			/* zero small motions for agg -- cf. pkdDEMZeroSmallMotions() */
			double dAccCritSq = msr->param.FO.DP.dAccCrit*msr->param.FO.DP.dAccCrit;
			double dDeltaSq = msr->param.dDelta*msr->param.dDelta; /* multistepping...? */
			double r2,dp,a2,v2;
			r2 = agg->rad_max*agg->rad_max;
			dp = vectorDot(agg->a_com,agg->v_com);
			if (dp < 0.0) {
				a2 = vectorMagSq(agg->a_com);
				if (a2 > 0.0 && a2 < dAccCritSq) {
					v2 = vectorMagSq(agg->v_com);
					if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq) {
						vectorZero(agg->v_com);
						vectorZero(agg->a_com);
						}
					}
				}
			Vector wDot;
			vectorSet(wDot, /* we can do this in the body frame -- easiest! */
					  agg->torque[0]/agg->moments[0],
					  agg->torque[1]/agg->moments[1],
					  agg->torque[2]/agg->moments[2]);
			dp = vectorDot(wDot,agg->omega);
			if (dp < 0.0) {
				a2 = vectorMagSq(wDot);
				if (a2 > 0.0 && a2 < dAccCritSq*r2) {
					v2 = vectorMagSq(agg->omega);
					if (v2 > 0.0 && v2 < dAccCritSq*dDeltaSq/r2) {
						vectorZero(agg->omega);
						vectorZero(agg->torque);
						}
					}
				}
			}
#endif
		/*DEBUG!if (agg->bAssigned) printf("AGG |ACC| = %g |TRQ| = %g\n",vectorMag(agg->a_com),vectorMag(agg->torque));*/
		}
	}

void msrAggsGetNewIdx(MSR msr)
{
	int i;

	/* find next unassigned aggregate; grow buffer if needed */

	for (i=0;i<msr->nAggs;i++)
		if (!msr->pAggs[i].bAssigned) {
			msr->iAggNewIdx = i;
			return;
			}

	msr->iAggNewIdx = msr->nAggs;
	msr->nAggs <<= 1; /* double buffer size */
	if (msr->param.bVDetails)
		(void) printf("msrAggsGetNewIdx(): Doubled aggregate buffer to %i.\n",msr->nAggs);
	msr->pAggs = (Aggregate *) realloc(msr->pAggs,msr->nAggs*sizeof(Aggregate));
	assert(msr->pAggs != NULL);
	for (i=msr->iAggNewIdx;i<msr->nAggs;i++) /* just to be sure */
		msr->pAggs[i].bAssigned = 0;
	}

void msrAggsFind(MSR msr)
{
	/*
	 ** Finds and initializes aggregates.  Must be called immediately
	 ** after initial conditions are loaded.
	 */

	struct outAggsFind outFind;
	Aggregate *agg;
	int i,nPart,nAssigned = 0;

#ifdef FAST_AGGS
	/*
	** If using FAST_AGGS, initialize PKD-level variables to enable cache
	** line approach for agg constituent searching.
	*/
	msrAggsCacheInit(msr);
#endif

	if (msr->param.bVDetails)
		(void) printf("msrAggsFind(): Allocating space for aggregates...\n");

	/*
	 ** We don't know in advance how many aggregates there will be.
	 ** Worse, the aggregate indices could be in any order, and may
	 ** not even be continuous.  Strategy: loop through all particles
	 ** to find the largest aggregate index (stored as -1 - iOrgIdx)
	 ** and assume storage for that many aggregates (plus one, since
	 ** aggregate indices start at 0) is required.  Some slots may be
	 ** left empty though.
	 */

	pstAggsFind(msr->pst,NULL,0,&outFind,NULL);
	msr->nOrg = outFind.iMaxOrg;
	if (outFind.iMaxIdx == -1) {
		/* no aggregates found -- make some space anyway */
		msr->nAggs = AGGS_INIT_BUFF_SIZE;
		}
	else {
		/* round up buffer size to next power of 2 */
		msr->nAggs = 1 << (int) (log((outFind.iMaxIdx + 1) << 1)/M_LN2);
		if (msr->nAggs < AGGS_INIT_BUFF_SIZE)
			msr->nAggs = AGGS_INIT_BUFF_SIZE;
		}
	assert(msr->nAggs > 1);

	msr->pAggs = (Aggregate *) malloc(msr->nAggs*sizeof(Aggregate));
	assert(msr->pAggs != NULL);

	if (msr->param.bVDetails)
		(void) printf("msrAggsFind(): Space for %i aggregates allocated.\n",msr->nAggs);


	/* initialize */
#ifdef FAST_AGGS
	msrAggsReorder(msr);	/* Reorder for FAST_AGGS */
#endif
	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		agg->bAssigned = 0;
		/* determine which aggregates are actually occupied */
		if (i <= outFind.iMaxIdx) {
			msrAggsCountPart(msr,i,agg,&nPart);
			if (nPart > 1) {
				agg->bAssigned = 1;
				++nAssigned;
				msrAggsUpdate(msr,i,agg,0/*accel computed later*/,0.0/*start of step*/);
				}
			else if (nPart == 1) {
				fprintf(stderr,"ERROR: aggregate %i has only one particle -- aborting.\n",i);
				_msrExit(msr,1);
				}
			}
		}

	if (msr->param.bVDetails)
		(void) printf("msrAggsFind(): %i aggregate%s initialized.\n",nAssigned,
					  nAssigned == 1 ? "" : "s");
	/* determine what next aggregate index should be, when needed */

	msrAggsGetNewIdx(msr);	//This is a hard-sphere thing. DeMartini debug
	}

#ifdef AGGS_IN_PATCH
void msrAggsSetSpacePos(MSR msr,int iAggIdx,Aggregate *agg,double dTime)
#else
void msrAggsSetSpacePos(MSR msr,int iAggIdx,Aggregate *agg)
#endif
{
	/* aggregate c-o-m pos and body orientation must be up to date */

	struct inAggsSetSpacePos in;
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	in.iAggIdx = iAggIdx;
	vectorCopy(agg->r_com,in.r_com);
	matrixCopy(agg->lambda,in.lambda);
#ifdef AGGS_IN_PATCH
	in.dTime = dTime;
	in.PP = msr->param.PP;
#endif
	pstAggsSetSpacePos(msr->pst,&in,sizeof(in),NULL,NULL);
}

void msrAggsSetSpaceVel(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* aggregate c-o-m vel, spin, & body orientation must be up to date */

	struct inAggsSetSpaceVel in;
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	in.iAggIdx = iAggIdx;
    in.dDelta = msr->param.dDelta;
	vectorCopy(agg->v_com,in.v_com);
    vectorCopy(agg->a_com,in.a_com);
	vectorCopy(agg->omega,in.omega);
    /* Calculate omegadot */
    in.omegadot[0] = (agg->torque[0] + agg->omega[1]*agg->omega[2]*(agg->moments[1] - agg->moments[2]))/agg->moments[0];
    in.omegadot[1] = (agg->torque[1] + agg->omega[2]*agg->omega[0]*(agg->moments[2] - agg->moments[0]))/agg->moments[1];
    in.omegadot[2] = (agg->torque[2] + agg->omega[0]*agg->omega[1]*(agg->moments[0] - agg->moments[1]))/agg->moments[2];
	matrixCopy(agg->lambda,in.lambda);
#ifdef AGGS_IN_PATCH
	in.PP = msr->param.PP;
#endif
	pstAggsSetSpaceVel(msr->pst,&in,sizeof(in),NULL,NULL);
}

void msrAggsSetSpaceSpins(MSR msr,int iAggIdx,Aggregate *agg)
{
	/* aggregate spin and body orientation must be up to date */

	struct inAggsSetSpaceSpins in;
	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);
	in.iAggIdx = iAggIdx;
	matrixTransform(agg->lambda,agg->omega,in.omega); /* space frame */ /* RP-DEBUG: TAKE NOTE! Particle spins are in the space frame - or at least, they are trying to be.  But with aggs_in_patch, this transforms from body to patch frame? */
	pstAggsSetSpaceSpins(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrAggsKick(MSR msr,double dt)
{
	/* dt is one half current drift interval */

	Aggregate *agg;
	int i,k;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			for (k=0;k<3;k++)
				agg->v_com[k] += agg->a_com[k]*dt;
#ifndef DEM
			msrAggsSetSpaceVel(msr,i,agg);
#endif
		}
	}
}

#ifdef AGGS_IN_PATCH

void msrAggsInPatchKick(MSR msr,double dt,int bOpen)
{
	/* dt is one half current drift interval */

	double dOrbFreq = msr->param.PP.dOrbFreq;

	Aggregate *agg;
	int i,k;

	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			if (!bOpen) {
				/* perform Cross Hamiltonian */
				agg->v_com[0] += 2.0*dt*dOrbFreq*agg->dPy_com;
				agg->v_com[1] = agg->dPy_com - 2.0*dOrbFreq*agg->r_com[0];
				}
			for (k=0;k<3;++k)
				agg->v_com[k] += agg->a_com[k]*dt;
			if (bOpen) {
				agg->dPy_com = agg->v_com[1] + 2.0*dOrbFreq*agg->r_com[0];
				/* perform Cross Hamiltonian */
				agg->v_com[0] += 2.0*dt*dOrbFreq*agg->dPy_com;
				agg->v_com[1] = agg->dPy_com - dOrbFreq*agg->r_com[0]
				  - dOrbFreq*(agg->r_com[0] + 2.0*dt*agg->v_com[0]);
				}
			msrAggsSetSpaceVel(msr,i,agg);
			}
		}
	}

#endif /* AGGS_IN_PATCH */

void msrAggsAdvanceOpen(MSR msr)
{
	Aggregate *agg;
	int i,nAssigned;

	nAssigned = 0;
	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned) {
			agg->dLastUpdate = 0.0;
			++nAssigned;
			}
		}
	if (msr->param.bVStep)
		(void) printf("%i aggregate%s currently active.\n",nAssigned,nAssigned==1?"":"s");
	}

#define EPS 1.0e-6 /* for sanity check */

void msrAggsAdvance(MSR msr,int iAggIdx,Aggregate *agg,double dToTime)
{
	/*
	** Integrate to get new lambda, omega, and COM position.
	*/

	double dt;
	int k;

	FLOAT fPeriod[3]; /* for periodic boundary conditions only */

	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);

	dt = dToTime - agg->dLastUpdate;

	if (dt == 0.0)
		return;

	if (msr->param.CP.iOverlapOption != OverlapBackstep) {
		assert(dt >= 0.0);
		assert(dt <= msrDelta(msr));
		}

	/*
	** Advance aggregate spin & orientation using adaptive
	** Runge-Kutta.
	*/

	aggsRungeAdvance(agg,msrDelta(msr),dt); /* smaller than dDelta for multistepping? */

	/*
	** Sanity check: lambda should be orthonormal -- easiest to just
	** ensure no element > 1 (with a fudge factor for roundoff).
	*/

	{
		static int bFirstWarning = 1;
		int j,bOK=1;
		for (k=0;k<3;k++)
			for (j=0;j<3;j++)
			        if (fabs(agg->lambda[k][j]) > 1.0 + EPS) {
					bOK = 0;
					if (bFirstWarning && msr->param.bVWarnings) {
						bFirstWarning = 0;
						fprintf(stderr,"WARNING: msrAggsAdvance(): unphysical orientation matrix, agg %i, elem [%i][%i] = %.16e\n",
								iAggIdx,k,j,agg->lambda[k][j]);
						}
					}
		if (!bOK) { /* enforce orthonormality */
			matrixOrthonormalize(agg->lambda,0); /*Should just need this, instead of block below */

/*			Vector vTmp;
			double dDot;
			dDot = vectorDot(agg->lambda[0],agg->lambda[1]);
			vectorScale(agg->lambda[0],dDot,vTmp);
			vectorSub(agg->lambda[1],vTmp,agg->lambda[1]);*/ /* remove component in 1st direction from 2nd vector */
/*			dDot = vectorDot(agg->lambda[0],agg->lambda[2]);
			vectorScale(agg->lambda[0],dDot,vTmp);
			vectorSub(agg->lambda[2],vTmp,agg->lambda[2]);*/ /* ditto for 3rd vector */
/*			dDot = vectorDot(agg->lambda[1],agg->lambda[2]);
			vectorScale(agg->lambda[1],dDot,vTmp);
			vectorSub(agg->lambda[2],vTmp,agg->lambda[2]);*/ /* finally, remove component in 2nd direction from 3rd vector */

			/* vectors should now be orthogonal (unless 2 or more are the same!) -- renormalize */
/*			vectorNorm(agg->lambda[0]);
			vectorNorm(agg->lambda[1]);
			vectorNorm(agg->lambda[2]);*/

			/* sanity check */
			assert(fabs(vectorDot(agg->lambda[0],agg->lambda[1])) < EPS);
			assert(fabs(vectorDot(agg->lambda[0],agg->lambda[2])) < EPS);
			assert(fabs(vectorDot(agg->lambda[1],agg->lambda[2])) < EPS);
			}
		}

	/* COM "drift" + periodic boundary conditions */
	if (msr->param.bPeriodic) {
		fPeriod[0] = msr->param.dxPeriod;
		fPeriod[1] = msr->param.dyPeriod;
		fPeriod[2] = msr->param.dzPeriod;
		}
	else {	/* Don't do PBCs if bPeriodic not set; override given boundary lengths */
		fPeriod[0] = FLOAT_MAXVAL;
		fPeriod[1] = FLOAT_MAXVAL;
		fPeriod[2] = FLOAT_MAXVAL;
		}

	for (k=0;k<3;k++) {
		agg->r_com[k] += agg->v_com[k]*dt;
		if (msr->param.bPeriodic) {	/* Aggs PBCs */
			while (agg->r_com[k] >= msr->fCenter[k] + 0.5*fPeriod[k]) {
				agg->r_com[k] -= fPeriod[k];
				}
			while (agg->r_com[k] < msr->fCenter[k] - 0.5*fPeriod[k]) {
				agg->r_com[k] += fPeriod[k];
				}
			}
		}

#ifdef AGGS_IN_PATCH
	msrAggsInPatchCheckCOM(msr,iAggIdx,agg,msr->dTime + dToTime);
#endif

#ifdef AGGS_IN_PATCH
	/*
	** Need to adjust principal axes orientations (the lambda matrix)
	** to account for rotation of the patch frame.  Do this now,
	** *before* computing the new "space" (really "patch") positions,
	** velocities, and spins.  The rotation matrix represents a
	** clockwise rotation through an angle Omega * (time interval), as
	** seen looking down on the patch from the positive z axis.
	*/
	msrAggsInPatchRotateAxes(msr,agg,dt);
#endif

	/*
	** Get new particle positions in space frame.  These are needed
	** during the drift step (instead of just at the end) for
	** resolving collisions, i.e. determining the point of impact.
	*/

#ifdef AGGS_IN_PATCH
	msrAggsSetSpacePos(msr,iAggIdx,agg,msr->dTime + dToTime);
#else
	msrAggsSetSpacePos(msr,iAggIdx,agg);
#endif

	/*
	** Get new particle velocities in space frame.  These are also set
	** during kicks and as part of post-collision backdrifting, but it
	** doesn't hurt to do it here as well.  This ensures that
	** velocities are updated prior to stress checks (so released
	** particles have the correct velocity).  Note stress checks are
	** done after the drift as part of the gravity calculation, so we
	** could just update velocities in msrAggsAdvanceClose(), but
	** updating them now has the minor advantage of providing more
	** accurate stick-vs-frag checks in pkdAggsDoCollision().
	*/

	msrAggsSetSpaceVel(msr,iAggIdx,agg);

	/*
	** Get new particle spins in space frame.  These are simply set to
	** the new spin vector of the aggregate itself.  This update is
	** needed in order to compute angular momentum correctly, both for
	** outputs and for msrAggsGetAxesAndSpin(), which uses the total
	** angular momentum to determine the spin vector of a new
	** aggregate formed by merger, etc.
	*/

	msrAggsSetSpaceSpins(msr,iAggIdx,agg);
	/* finally, update the clock */

	agg->dLastUpdate = dToTime;
	}

#undef EPS

void msrAggsAdvanceClose(MSR msr,double dt)
{
	/* dt is current drift interval */

	Aggregate *agg;
	int i;

#ifdef FAST_AGGS
	msrAggsReorder(msr);	/* Reorder for FAST_AGGS, to be safe */
#endif /* FAST_AGGS */
	for (i=0;i<msr->nAggs;i++) {
		agg = &msr->pAggs[i];
		if (agg->bAssigned)
			msrAggsAdvance(msr,i,agg,dt);
		}

#ifdef FAST_AGGS  /* This is necessary to reset the tree */
	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
        msrDomainDecomp(msr,0,1);
#else /* FAST_AGGS */
        if (msr->param.bPeriodic) {
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
                msrDomainDecomp(msr,0,1);
		}
#endif /* !FAST_AGGS */
	}

void msrAggsBackDrift(MSR msr,int iAggIdx,Aggregate *agg,double dt)
{
	/*
	** After a collision, it is necessary to drift aggregate "ghost"
	** particle positions back to start of step so we can continue to
	** do forward collision prediction during the step.  First we get
	** the post-collision space velocities of each particle (including
	** the second-order term), and then we carry out the back drift.
	** We do *not* drift the aggregate center of mass backward, so
	** there is no need for fancy timekeeping between this function
	** and msrAggsAdvance() (which drifts the center of mass forward
	** from the last update time).  NOTE: affected particles also have
	** SMOOTHACTIVE set by this routine to force recomputation of
	** collision circumstances.
	*/

	struct inAggsBackDrift in;

	assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
	assert(agg == &msr->pAggs[iAggIdx]);
	assert(agg->bAssigned);

	/* first update particle velocities to second order */

	msrAggsSetSpaceVel(msr,iAggIdx,agg);

	/* RP-DEBUG-dPy: Positions set below in pkdAggsBackDrift().
	   If that backdrift takes particles outside the
	   patch, the collision search should still match 
	   colliding pairs correctly.  (c.f. Tom Q. email, Oct 2009)  */

	/*
	** Also update particle spins because agg rotation may have
	** changed.  This is needed because this agg may undergo another
	** collision during this step.
	*/

	msrAggsSetSpaceSpins(msr,iAggIdx,agg);
	/* now drift */

	in.iAggIdx = iAggIdx;
	in.dt = dt;
	pstAggsBackDrift(msr->pst,&in,sizeof(in),NULL,NULL);

	/*
	** An alternative to the second-order back drift is to actually do
	** a reverse advance using the Runge-Kutta integrator.  To
	** implement, comment out the call to msrAggsSetSpaceVel() above
	** and also comment out the back drift loop over k in
	** pkdAggsBackDrift() (but still set SMOOTHACTIVE), and uncomment
	** the following code.  Presumably this is more expensive, and
	** seems to be overkill so long as the back step is small
	** (e.g. less than dDelta; it can be larger when position fixes
	** using negative timesteps are allowed---but if it's too large
	** even the reverse advance will not prevent noticeably unphysical
	** results).
	*/
	/*
	aggsRungeAdvance(agg,-dt);
	agg->r_com[0] -= agg->v_com[0]*dt;
	agg->r_com[1] -= agg->v_com[1]*dt;
	agg->r_com[2] -= agg->v_com[2]*dt;
	agg->dLastUpdate -= dt;
	NOTE: if AGGS_IN_PATCH, would need to derotate body axes too...
	*/
}

void msrAggsMerge(MSR msr,const COLLIDER *c1,const COLLIDER *c2,double dImpactTime,COLLIDER *cOut)
{
	/*
	  Note: msr->dTime => total sim time elapsed
	  msr->param.dDelta => timestep
	  dImpactTime => elapsed time since start of step
	*/
	struct inAggsMerge inMerge;
	Aggregate *agg;
	int iAggIdx;

	if (COLLIDER_IS_AGG(c1) && COLLIDER_IS_AGG(c2)) {
		assert(COLLIDER_AGG_IDX(c1) != COLLIDER_AGG_IDX(c2));
		if (COLLIDER_AGG_IDX(c1) < COLLIDER_AGG_IDX(c2)) {
			inMerge.iOldIdx = COLLIDER_AGG_IDX(c2);
			inMerge.iNewIdx = COLLIDER_AGG_IDX(c1);
			*cOut = *c1; /* struct copy */
			}
		else { /* i.e., COLLIDER_AGG_IDX(c2) < COLLIDER_AGG_IDX(c1) */
			inMerge.iOldIdx = COLLIDER_AGG_IDX(c1);
			inMerge.iNewIdx = COLLIDER_AGG_IDX(c2);
			*cOut = *c2;
			}
		pstAggsMerge(msr->pst,&inMerge,sizeof(inMerge),NULL,NULL);
		msr->pAggs[inMerge.iOldIdx].bAssigned = 0;
		iAggIdx = inMerge.iNewIdx;
		if (msr->param.bVDetails)
			(void) printf("Aggregate %i merged with %i.\n",inMerge.iOldIdx,inMerge.iNewIdx);
		}
	else if (COLLIDER_IS_AGG(c1) && !COLLIDER_IS_AGG(c2)) {
		iAggIdx = COLLIDER_AGG_IDX(c1);
		*cOut = *c1;
		}
	else if (COLLIDER_IS_AGG(c2) && !COLLIDER_IS_AGG(c1)) {
		iAggIdx = COLLIDER_AGG_IDX(c2);
		*cOut = *c2;
		}
	else { /* i.e., neither collider is an aggregate */
		msr->pAggs[iAggIdx = msr->iAggNewIdx].bAssigned = 1;
		if (msr->param.bVDetails)
			(void) printf("Aggregate %i created.\n",iAggIdx);
		*cOut = *c1; /* doesn't really matter which */
		cOut->id.iOrgIdx = -1 - iAggIdx; /* so output log knows this is an agg */
		msrAggsGetNewIdx(msr);
		}

	agg = &msr->pAggs[iAggIdx];
	assert(agg->bAssigned);

	/*
	 ** We need to compute the torque on the new merged body so that
	 ** it can be advanced properly.  However, we only have the
	 ** start-of-step particle accelerations, and the particles have
	 ** moved.  There's no easy fix to this, so we'll just use the
	 ** pre-existing accelerations and hope the timestep is small
	 ** enough that any introduced error is tolerable.  Note that we
	 ** do NOT check for excessive stress here, because the
	 ** accelerations are wrong and any liberated particles would need
	 ** to be flagged for future collision checks during this
	 ** interval, and would need to be backdrifted, etc. -- basically
	 ** too messy for now.
	*/
	
	msrAggsUpdate(msr,iAggIdx,agg,1/*do accel*/,dImpactTime);

#ifdef AGGS_IN_PATCH	
	  /* RP-DEBUG-dPy fix 9/20/09 - Py now calculated in msrAggsUpdate() */
#endif /*AGGS_IN_PATCH*/

	msrAggsBackDrift(msr,iAggIdx,agg,dImpactTime);

	cOut->agg = *agg; /* struct copy for output log */
	}

void msrAggsBounce(MSR msr,const COLLIDER *c1,const COLLIDER *c2,double dImpactTime,int bDoUpdate)
{
	const COLLIDER *c;
	Aggregate *agg;
	int i,iAggIdx;

	for (i=1;i<=2;i++) {
		c = (i == 1 ? c1 : c2);
		if (COLLIDER_IS_AGG(c)) {
			iAggIdx = COLLIDER_AGG_IDX(c);
			assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
			agg = &msr->pAggs[iAggIdx];
			assert(agg->bAssigned);
			*agg = c->agg; /* struct copy to get new velocity/spin/Py */
			if (bDoUpdate) { /* needed if overlap adjusted particle position(s) */
				assert(dImpactTime == 0.0);
				msrAggsUpdate(msr,iAggIdx,agg,1/*do accel because positions changed*/,0.0/*overlap errors can only be detected at start of step (SOS)*/); /* RP-DEBUG-dPy: Will this 0.0 argument be a problem? AdjPos moves particles in aggs apart (but NOT COMs) at SOS, thus the COM's Py needs to be computed at SOS. Moral: adjPos and AGGS don't play nice together*/
				}
			msrAggsBackDrift(msr,iAggIdx,agg,dImpactTime);
			}
		}
	}

void msrAggsFrag(MSR msr,const COLLIDER *c1,const COLLIDER *c2,double dImpactTime)
{
	/*
	** Updates any remaining aggs following a frag event (any
	** liberated particles have already been dealt with).  Here a
	** "frag event" consists of a free particle hitting an agg
	** particle (or two particles -- one from each of two colliding
	** aggs -- hitting each other) at sufficient speed to knock the
	** agg particle(s) loose.  The input COLLIDER structs contain the
	** identifiers of the one or two aggs involved; these are used to
	** update and backdrift the aggs that have lost particles.
	*/

	const COLLIDER *c;
	Aggregate *agg;
	int i,iAggIdx,nPart;

	for (i=1;i<=2;i++) {
		c = (i == 1 ? c1 : c2);
		if (COLLIDER_IS_AGG(c)) {
			iAggIdx = COLLIDER_AGG_IDX(c);
			assert(iAggIdx >= 0 && iAggIdx < msr->nAggs);
			agg = &msr->pAggs[iAggIdx];
			assert(agg->bAssigned);
			msrAggsCountPart(msr,iAggIdx,agg,&nPart);
			assert(nPart > 0); /* must be at least 1 particle left */
			if (nPart > 1) {
				msrAggsUpdate(msr,iAggIdx,agg,1/*do accel (torques change because COM changed)*/,dImpactTime);
				msrAggsBackDrift(msr,iAggIdx,agg,dImpactTime);
				}
			else
				msrAggsDelete(msr,iAggIdx,agg,dImpactTime); /* if only 1 particle, already backdrifted in pkdAggsDoCollision() */
			}
		}
	}

#endif /* AGGS */

#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
void msrDustBinsApply(MSR msr)
{
	/*
	 ** Adds mass to planetesimals from dust bins in two steps.
	 ** 1. Compute total mass accreted regardless of mass in bins.
	 ** 2. Adjust mass accreted so total does not exceed bin mass,
	 ** then apply change to each affected planetesimal.
	 */

	DUST_BINS_PARAMS *DB;
	struct inDustBinsGetMass in_get;
	struct outDustBinsGetMass out_get;
	struct outDustBinsInclAvg out_avg;
	struct outDustBinsInclMax out_max;
	struct inDustBinsApply in_do;
	double incl_avg;
	int i;
#ifdef ORIGIN_HISTOGRAM
	int j;
#endif /* ORIGIN_HISTOGRAM */

	/* step 1a: get total mass accreted by planetesimals */

	DB = &msr->param.CP.DB; /* shorthand */
	in_get.DB = *DB; /* struct copy */
	for (i=0;i<DB->nDustBins;i++)
		in_get.aDustBins[i] = in_do.aDustBins[i] = msr->aDustBins[i]; /* struct copy */
	in_get.dTimeInt = in_do.dTimeInt = DB->iDustBinsApplyInt*msr->param.dDelta;
	in_get.dCentMass = in_do.dCentMass = msr->param.dCentMass;
	pstDustBinsGetMass(msr->pst,&in_get,sizeof(in_get),&out_get,NULL);
	/* step 1b: determine debris inclination dispersion */
	/* do parameter test ...*/
	if (DB->iDustBinsVelDispOpt == 0)
	  in_do.dVdisp = 0.0;
	else if (DB->iDustBinsVelDispOpt == 1) {
	  pstDustBinsInclAvg(msr->pst,NULL,0,&out_avg,NULL);
	  incl_avg = out_avg.dDustBinsInclTot/msr->N;
	  in_do.dVdisp = SQ(sin(incl_avg));
	} else if (DB->iDustBinsVelDispOpt == 2) {
	  pstDustBinsInclMax(msr->pst,NULL,0,&out_max,NULL);
	  in_do.dVdisp = out_max.dDustBinsInclMax;
	}

	/* step 2a: adjust values to not exceed bin masses */

	for (i=0;i<DB->nDustBins;i++) {
		if (out_get.aDustBinsMassLoss[i] > msr->aDustBins[i].dMass) {
			printf("iBin = %i, mass = %e, mass loss = %e\n", 
			       i,msr->aDustBins[i].dMass, out_get.aDustBinsMassLoss[i]);	   
			in_do.aMassIncrFrac[i] = msr->aDustBins[i].dMass / out_get.aDustBinsMassLoss[i];
			msr->aDustBins[i].dMass = 0.0;
#ifdef ORIGIN_HISTOGRAM
			for (j=0;j<NUM_ORIGIN_BINS;j++)
				msr->aDustBins[i].origin_bins[j] = 0.0;
#endif /* ORIGIN_HISTOGRAM */
			}
		else {
			in_do.aMassIncrFrac[i] = 1.0;
			msr->aDustBins[i].dMass -= out_get.aDustBinsMassLoss[i];
			}
/*		printf("iBin = %i, MassIncr = %f\n", i, in_do.aMassIncrFrac[i]);*/
		}

	in_do.nBins = DB->nDustBins;

#ifdef ORIGIN_HISTOGRAM
	for (i=0;i<DB->nDustBins;i++)
		for (j=0;j<NUM_ORIGIN_BINS;j++)
			in_do.aOriginBins[i][j] = msr->aDustBins[i].origin_bins[j];
#endif /* ORIGIN_HISTOGRAM */

	/* step 2b: apply changes to planetesimals */
	pstDustBinsApply(msr->pst,&in_do,sizeof(in_do),NULL,NULL);
	}
#endif /* RUBBLE_ZML || COLLMOD_ZML */

#ifdef RUBBLE_ZML

void msrRubbleResetColFlag(MSR msr)
{
	pstRubbleResetColFlag(msr->pst,NULL,0,NULL,NULL);
	}

void msrRubbleStep(MSR msr)
{
	/* this is a bit inefficient: 1) we don't need to recompute
	  dMinStep each time, since it's a fixed value; 2) we really just
	  want to set the particle rungs directly, but it's probably safer
	  to stick with the usual method of setting the timesteps first,
	  then calling msrDtToRung(). */

	struct inRubbleStep in;

	in.dMaxStep = msr->param.dDelta;
	in.dMinStep = msr->param.dDelta/(1<<(msr->param.iMaxRung - 1));
	pstRubbleStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#elif defined(COLLMOD_ZML)

void msrCollModResetColFlag(MSR msr)
{
	pstCollModResetColFlag(msr->pst,NULL,0,NULL,NULL);
	}

void msrCollModStep(MSR msr)
{
	/*DEBUG this is a bit inefficient: 1) we don't need to recompute
	  dMinStep each time, since it's a fixed value; 2) we really just
	  want to set the particle rungs directly, but it's probably safer
	  to stick with the usual method of setting the timesteps first,
	  then calling msrDtToRung().*/

	struct inCollModStep in;

	in.dMaxStep = msr->param.dDelta;
	in.dMinStep = msr->param.dDelta/(1<<(msr->param.iMaxRung - 1));
	pstCollModStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#endif /* RUBBLE_ZML, COLLMOD_ZML */
	
#ifdef RUBBLE_ZML
	
void msrRubCleanup(MSR msr, double dTime)
{
	/*
	 ** Need to check if any rubble pile particles should be converted back
	 ** to full-fledged planetesimals. For expediency, only do this check at
	 ** the start of every full (rung 0) timestep. This means some particles
	 ** may stay as rubble pile particles a little longer than expected (but
	 ** never more than one full timestep).
	 */

	DUST_BINS_PARAMS *DB;
	struct inRubCleanup in;
	struct outRubCleanup out;
	int bCleanup = 0;
	int i,j;

	in.DB = *(DB = &msr->param.CP.DB); /* shorthand */
	in.dCentMass = msr->param.dCentMass;

/*	printf("In RubCleanup msr->re.nEvents = %d\n", msr->re.nEvents);*/

/*	for (i=msr->re.nEvents - 1;i>=0;i--) */
	for (i=0;i<msr->re.nEvents;i++)
		if (dTime >= msr->re.rc[i].dTEndRubblePhase) {
			bCleanup = 1;
			in.iColor = msr->re.rc[i].iColor;
			printf("iColor = %d\n", in.iColor);
			assert(in.iColor != PLANETESIMAL);
			pstRubCleanup(msr->pst,&in,sizeof(in),&out,NULL);
			/*
			** A note about dust and origin histograms: each processor
			** must determine its own contribution to the global dust
			** profile separately, then these contributions must be
			** mass-weighted together to determine the final change to
			** the global dust distribution following the cleanup.
			** This is why we can't simply send the global data down
			** to the pkd level: the processors would not know about
			** each other's changes, and the final distribution would
			** be wrong.
			*/
			for (j=0;j<DB->nDustBins;j++) {
				if (out.aDustBins[j].dMass > 0.0) {
#ifdef ORIGIN_HISTOGRAM
					MergeHistograms(msr->aDustBins[j].dMass, msr->aDustBins[j].origin_bins, out.aDustBins[j].dMass, out.aDustBins[j].origin_bins);
#endif /* ORIGIN_HISTOGRAM */
					msr->aDustBins[j].dMass += out.aDustBins[j].dMass;
					}
				}
			if (out.DustBinsTrash.dMass > 0.0) {
#ifdef ORIGIN_HISTOGRAM
				MergeHistograms(msr->DustBinsTrash.dMass, msr->DustBinsTrash.origin_bins, out.DustBinsTrash.dMass, out.DustBinsTrash.origin_bins);
#endif /* ORIGIN_HISTOGRAM */
				msr->DustBinsTrash.dMass += out.DustBinsTrash.dMass; /* dust outside bin range goes into trash */
				}
			for (j=i--;j<msr->re.nEvents - 1;j++)
				msr->re.rc[j] = msr->re.rc[j + 1];
			--msr->re.nEvents; /* decrease the number of active rubble events */
			}
	if (bCleanup) {
		msrAddDelParticles(msr); /* clean up any deletions */
		if (msr->param.nSmooth > msr->N) {
			msr->param.nSmooth = msr->N;
			if (msr->param.bVWarnings)
				printf("WARNING: msrDoCollisions() RUBBLE_ZML: nSmooth reduced to %i\n",msr->param.nSmooth);
			}
		}
	}

#endif /* RUBBLE_ZML */

#ifdef ORIGIN_HISTOGRAM

void msrInitializeOriginBins(MSR msr) {
	struct inInitializeOriginBins in;

	printf("Initializing origin bins with original positions.\n");
	in.DB = msr->param.CP.DB; /* struct copy */
	pstInitializeOriginBins(msr->pst, &in, sizeof(in), 0, 0);
	/* dust bin origin histograms initialized in msrInitialize() */
	}

#endif /* ORIGIN_HISTOGRAM */ 

#ifdef WALLS

void msrWallsGetData(MSR msr,const char achFilenameWithQuotes[])
{
	const int bVerbose = (msr->param.bVWarnings || msr->param.bVDetails);

	WALL_PARAMS *WP = &msr->param.CP.WP;
	WALL_DATA *wd = NULL; /* store locally first */
	WALL *w;
	FILE *fp;
	char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
	double dTime;
	int i;

	WP->nWalls = 0;
	assert(achFilenameWithQuotes != NULL);
	_msrStripQuotes(achFilenameWithQuotes,achFilename);
	if (!strlen(achFilename))
		return;
	_msrMakePath(msr->param.achDataSubPath,achFilename,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
	if (bVerbose)
		printf("Reading wall data from \"%s\"...\n",achFilename);
	if (!(fp = fopen(achFilename,"r"))) {
		fprintf(stderr,"Unable to open \"%s\".\n",achFilename);
		goto abort;
		}
	if (wallsParseWallsFile(fp,&WP->nWalls,&wd,&dTime,bVerbose) != 0) {
		fprintf(stderr,"Error occurred with parsing wall data.\n");
		goto abort;
		}
	if (WP->nWalls > MAX_NUM_WALLS) {
		fprintf(stderr,"Number of walls (%i) exceeds maximum (%i) set in walls.h.\n",WP->nWalls,MAX_NUM_WALLS);
		goto abort;
		}
	if (bVerbose && dTime != 0.0)
		fprintf(stderr,"WARNING: time stamp (%g) in walls file ignored.\n",dTime); /*DEBUG fix this?*/
	for (i=0;i<WP->nWalls;i++) {
		w = &WP->pWalls[i];
		w->iWallID = i;
		w->wd = wd[i]; /* struct copy */
		vectorZero(w->vOscVel);
		vectorZero(w->vTravel);
		vectorZero(w->vTotVel);
		}
	fclose(fp);
	free((void *) wd);
	if (bVerbose)
		printf("Number of walls read = %i.\n",WP->nWalls);
#ifdef CHARGE
	/*
	** This is a HACK to record the area of the cylinder end caps
	** (needed in pkdChargeZApplyMoments()).
	*/
	msr->param.CP.CP.dArea = 0.0;
	for (i=0;i<WP->nWalls;i++) {
		w = &WP->pWalls[i];
		if (w->wd.iType == WallDisk) { /* assuming first disk found is the one we want! */
			msr->param.CP.CP.dArea = M_PI*w->wd.dRadius*w->wd.dRadius;
			printf("CHARGE: wall %i is a disk of area %g\n",i,msr->param.CP.CP.dArea);
			break;
			}
		}
	assert(msr->param.CP.CP.dArea > 0.0);
#endif /* CHARGE */
	return;
 abort:
	if (fp != NULL)
		fclose(fp);
	if (wd != NULL)
		free((void *) wd);
	_msrExit(msr,1);
	}

void msrWallsInitPosAndVel(MSR msr,double dTime)
{	
	/*
	** Initializes wall positions and velocities while accounting for
	** a possible non-zero start time (i.e., after a restart, or if
	** iStartStep > 0), which would mean moving walls likely need to
	** be offset from the positions read from the walls data file.
	** This function also (re)initializes velocities and spins of any
	** particles stuck to moving walls (including rotating cylinders).
	**
	** Compare this function with msrWallsMove() below.
	*/

	struct inWallsUpdateStuckParticles in;
	WALL_PARAMS *WP = &msr->param.CP.WP;
	WALL *w;
	WALL_DATA *wd;
	Vector vTmp;
	int i;

	for (i=0;i<WP->nWalls;i++) {
		w = &WP->pWalls[i];
		vectorZero(w->vTravel);
		vectorZero(w->vTotVel);
		wd = &w->wd;
#ifdef DEM
		/* default behavior for SSDEM is to step at dDelta interval */
		/*DEBUG if dDelta ever changes during a run, this will break*/
		if (wd->dStep == 0.0)
			wd->dStep = msr->param.dDelta;
#endif /* DEM */
		if (wd->dOscAmp != 0.0) {
			/*
			** If wall stepping is off, treat wall motion as smooth
			** (even though it's only updated every timestep).
			** Otherwise, treat the wall as stationary between steps.
			*/
			double t,dt = msr->param.dDelta;
			if (wd->dStep == 0.0) {
				/*
				** Because in msrWallsMove() we update oscillating wall
				** positions in linear increments based on the walls'
				** velocities in the previous step, we can't just set the
				** offsets to be A*sin(w*t); this would generally
				** introduce a slight error.  Instead, we loop from time 0
				** to the current time in steps of dDelta, updating the
				** walls' positions as we go.  Nonetheless, roundoff error
				** will creep in here.
				*/
				/* account for starting phase... */
				vectorScale(wd->vOscVec,wd->dOscAmp*sin(wd->dOscPhase),w->vTravel); /* assumes travel is zero */
        if (wd->dMass ==0) {
          vectorAdd(wd->vOrigin,w->vTravel,wd->vOrigin);
          vectorAdd(wd->vCoG,w->vTravel,wd->vCoG);
          }
				vectorScale(wd->vOscVec,wd->dOscAmp*wd->dOscFreq*cos(wd->dOscPhase),w->vOscVel);
				for (t=0.0;t<dTime;t+=dt) {
					vectorScale(w->vOscVel,dt,vTmp);
					vectorAdd(w->vTravel,vTmp,w->vTravel);
          if (wd->dMass ==0) {
            vectorAdd(wd->vOrigin,vTmp,wd->vOrigin);
            vectorAdd(wd->vCoG,vTmp,wd->vCoG);
            }
					vectorScale(wd->vOscVec,wd->dOscAmp*wd->dOscFreq*cos(wd->dOscFreq*t + wd->dOscPhase),w->vOscVel);
					}
				/*
				** Allow for the possibility that dDelta has changed.
				** Note this will introduce a slight glitch, but it's
				** probably better than not attempting the correction
				** at all.  (Even better is to not change dDelta when
				** using oscillating walls!)  This is not an issue for
				** walls with stepped motion.
				*/
				vectorScale(w->vOscVel,dTime - ((int) (dTime/dt))*dt,vTmp);
				vectorAdd(w->vTravel,vTmp,w->vTravel);
        if (wd->dMass==0) {
          vectorAdd(wd->vOrigin,vTmp,wd->vOrigin);
          vectorAdd(wd->vCoG,vTmp,wd->vCoG);
          }
        vectorScale(wd->vOscVec,wd->dOscAmp*wd->dOscFreq*cos(wd->dOscFreq*dTime + wd->dOscPhase),w->vOscVel);
				}
			else {
				/* in this case we CAN use A*sin(w*t)... */
				t = ((int) (dTime/wd->dStep))*wd->dStep;
				vectorScale(wd->vOscVec,wd->dOscAmp*sin(wd->dOscFreq*t + wd->dOscPhase),w->vTravel); /* assumes travel is zero */
				if (wd->dMass==0) {
          vectorAdd(wd->vOrigin,w->vTravel,wd->vOrigin);
          vectorAdd(wd->vCoG,w->vTravel,wd->vCoG);
          }
				vectorScale(wd->vOscVec,wd->dOscAmp*wd->dOscFreq*cos(wd->dOscFreq*t + wd->dOscPhase),w->vOscVel);
				}
			vectorCopy(w->vOscVel,w->vTotVel); /* assumes total velocity is zero */
			}
		/* steady motion */
		if (wd->dStep == 0.0)
			vectorScale(wd->vVel,dTime,vTmp); /* smooth */
		else 
			vectorScale(wd->vVel,((int) (dTime/wd->dStep))*wd->dStep,vTmp); /* stepped */
		vectorAdd(w->vTravel,vTmp,w->vTravel);
    if (wd->dMass==0) {
      vectorAdd(wd->vOrigin,vTmp,wd->vOrigin);
      vectorAdd(wd->vCoG,vTmp,wd->vCoG);
      }
		vectorAdd(w->vTotVel,wd->vVel,w->vTotVel);
		}
	/* update particle velocities and spins only */
	in.WP = *WP; /* struct copy */
	in.bUpdatePos = 0;
	in.dDelta = 0.0; /* even for rotating wall (particle positions assumed correct) */
	pstWallsUpdateStuckParticles(msr->pst,&in,sizeof(in),NULL,NULL);
	}

#ifdef DEM_WALLS_REACT
#define EPS 1.0e-6 /* desired fractional accuracy */
#define MIN_STEP 0.001 /* times dDelta */
#endif /* DEM_WALLS_REACT */

void msrWallsMove(MSR msr,double dTime,double dDelta)
{
    /*
    ** Given wall positions and velocities at dTime, advances wall
    ** positions and velocities to dTime + dDelta.  Also updates
    ** positions and velocities of any particles stuck to the walls.
    ** NOTE: hard-sphere needs vOscVel to have the correct value at
    ** the start of the NEXT drift (for collision detection) whereas
    ** soft-sphere needs it at the end of the CURRENT drift (after
    ** walls have moved, for applying friction forces).  Fortunately,
    ** these are equivalent!
    **
    ** This routine also updates position and orientation of
    ** reactive walls/assemblies.
    */

    struct inWallsUpdateStuckParticles in;
    WALL_PARAMS *WP = &msr->param.CP.WP;
    WALL *w;
    WALL_DATA *wd;
    Vector vTmp;
    double dNewTime = dTime + dDelta;
    int i;
#ifdef DEM_WALLS_REACT
    FLOAT hmin = MIN_STEP * dDelta;
    FLOAT fVars[12];
    int k;
    double dTmp;
#endif

    if (dDelta == 0.0)
	return; /* ignore initial drift probing boundary conditions */

#ifdef SEISMIC_WAVE
    double dStime;
    Vector vSVel;
    vectorZero(vSVel);
    FILE *seismic_file;
    seismic_file=fopen("seismic_file.txt","r");
#endif

    for (i=0;i<WP->nWalls;i++) {
	w = &WP->pWalls[i];
	vectorZero(w->vTravel);
	vectorZero(w->vTotVel);
	wd = &w->wd;
	if (wd->dOscAmp != 0.0) {
	    /*
	    ** Handle oscillation first, if applicable, because to
	    ** simplify collision prediction (in the hard-sphere
	    ** smooth travel case) we treat the motion as constant
	    ** during the drift step, so the new wall position is the
	    ** old position plus the step times the old oscillation
	    ** velocity.
	    */
	    if (wd->dStep == 0.0)
		vectorScale(w->vOscVel,dDelta,w->vTravel); /* smooth; assumes travel is zero */
	    else { /* stepped */
		Vector vTmp2;
		vectorScale(wd->vOscVec,wd->dOscAmp*sin(wd->dOscFreq*(((int) (dNewTime/wd->dStep))*wd->dStep) + wd->dOscPhase),vTmp);
		vectorScale(wd->vOscVec,wd->dOscAmp*sin(wd->dOscFreq*(((int) (dTime/wd->dStep))*wd->dStep) + wd->dOscPhase),vTmp2);
		vectorSub(vTmp,vTmp2,w->vTravel); /* assumes travel is zero */
		}
	    vectorScale(wd->vOscVec,wd->dOscAmp*wd->dOscFreq*cos(wd->dOscFreq*dNewTime + wd->dOscPhase),w->vOscVel);
	    vectorCopy(w->vOscVel,w->vTotVel);
	    }
#ifdef SEISMIC_WAVE
	{
	    double dDt2=0.5*dDelta;

	    while (!feof(seismic_file)) {
    
		fscanf(seismic_file,"%lf %lf %lf %lf",&dStime,&vSVel[0],&vSVel[1],&vSVel[2]);
		//printf("%.8g %.8g\n",dTime,dStime);
      
		if (dTime >= (dStime-dDt2) && dTime <= (dStime+dDt2)) {
		    vectorCopy(vSVel,wd->vVel);
		    //printf("%.8g %.8g %.8g %.8g\n",dStime, wd->vVel[0],wd->vVel[1],wd->vVel[2]);
		    break;
		    }

		}

	    fclose(seismic_file);
	    }
#endif /* SEISMIC WAVE */
	/* steady motion */
	if (wd->dStep == 0.0)
	    vectorScale(wd->vVel,dDelta,vTmp); /* smooth */
	else { /* stepped */
	    Vector vTmp2;
	    vectorScale(wd->vVel,((int) (dNewTime/wd->dStep))*wd->dStep,vTmp);
	    vectorScale(wd->vVel,((int) (dTime/wd->dStep))*wd->dStep,vTmp2);
	    vectorSub(vTmp,vTmp2,vTmp);
	    }
	vectorAdd(w->vTravel,vTmp,w->vTravel);
	vectorAdd(wd->vOrigin,w->vTravel,wd->vOrigin);
	vectorAdd(wd->vCoG,w->vTravel,wd->vCoG);
	vectorAdd(w->vTotVel,wd->vVel,w->vTotVel); /*DEBUG! wd->vVel here should really be wd->vVelPred for consistency*/

#ifdef DEM_WALLS_REACT
        if (wd->dMass != 0.0) {
            for (k=0;k<3;k++) {
                fVars[k] = wd->vRot[k];
                fVars[k + 3] = wd->vPxj[k]; /* q1 */
                fVars[k + 6] = wd->vPyj[k]; /* q2 */
                fVars[k + 9] = wd->vPzj[k]; /* q3 */
		}

            /* wall data structure provides torque and moments */

            diffeqIntegrate(fVars, 12, 0.0, dDelta, EPS, dDelta, hmin, wallsEulerDerivs, (void *) wd);

            /* Unpack fVars */
            for (k=0;k<3;k++) {
                wd->vRot[k] = fVars[k];
                wd->vPxj[k] = fVars[k + 3];
                wd->vPyj[k] = fVars[k + 6];
                wd->vPzj[k] = fVars[k + 9];
		}

	    /*The Following Block Rotates Coordinates based on wd->vRot*/

	    /* Assert that Vector1(dot)Vector2 = 0?*/
	    if ((wd->iType==WallTriangle) || (wd->iType==WallRectangle)) {
#ifdef QUATERNIONS /*DEBUG! needs testing*/
		Matrix a, b;
		Vector v;
		FLOAT fTheta = vectorMag(wd->vRot) * dDelta;
		vectorCopy(wd->vPxj, a[0]);
		vectorCopy(wd->vPyj, a[1]);
		vectorCopy(wd->vPzj, a[2]);
		matrixTranspose(a, b);
		matrixTransform(b, wd->vRot, vA);
		vectorNorm(vAxis);
		matrixQuatRot(v, fTheta, a); /* rotation matrix */
		vectorSub(wd->vOrigin, wd->vCoG, wd->vOrigin); //converting vectors relative to space to vector relative to body frame
		vectorSub(wd->vVertex1, wd->vCoG, wd->vVertex1); 
		vectorSub(wd->vVertex2, wd->vCoG, wd->vVertex2);
		vectorCopy(wd->vOrigin, v);
		matrixTransform(a, v, wd->vOrigin);
		vectorCopy(wd->vVertex1, v);
		matrixTransform(a, v, wd->vVertex1);
		vectorCopy(wd->vVertex2, v);
		matrixTransform(a, v, wd->vVertex2);
		vectorAdd(wd->vOrigin, wd->vCoG, wd->vOrigin); //converting vectors relative to space to vector relative to body frame
		vectorAdd(wd->vVertex1, wd->vCoG, wd->vVertex1); 
		vectorAdd(wd->vVertex2, wd->vCoG, wd->vVertex2);
#else /* QUATERNIONS */
		vectorAdd(wd->vVertex1,wd->vOrigin,wd->vVertex1); //converting vectors relative to origin to vectors relative to space frame
		vectorAdd(wd->vVertex2,wd->vOrigin,wd->vVertex2);

		vectorSub(wd->vOrigin, wd->vCoG, wd->vOrigin); //converting vectors relative to space to vector relative to body frame
		vectorRotate(wd->vOrigin, wd->vPxj, wd->vRot[0]*dDelta);
		vectorRotate(wd->vOrigin, wd->vPyj, wd->vRot[1]*dDelta);
		vectorRotate(wd->vOrigin, wd->vPzj, wd->vRot[2]*dDelta);
		vectorAdd(wd->vOrigin, wd->vCoG, wd->vOrigin);

		vectorSub(wd->vVertex1, wd->vCoG, wd->vVertex1);
		vectorRotate(wd->vVertex1, wd->vPxj, wd->vRot[0]*dDelta);
		vectorRotate(wd->vVertex1, wd->vPyj, wd->vRot[1]*dDelta);
		vectorRotate(wd->vVertex1, wd->vPzj, wd->vRot[2]*dDelta);
		vectorAdd(wd->vVertex1, wd->vCoG, wd->vVertex1);

		vectorSub(wd->vVertex2, wd->vCoG, wd->vVertex2);
		vectorRotate(wd->vVertex2, wd->vPxj, wd->vRot[0]*dDelta);
		vectorRotate(wd->vVertex2, wd->vPyj, wd->vRot[1]*dDelta);
		vectorRotate(wd->vVertex2, wd->vPzj, wd->vRot[2]*dDelta);
		vectorAdd(wd->vVertex2, wd->vCoG, wd->vVertex2);

		vectorSub(wd->vVertex1,wd->vOrigin,wd->vVertex1);
		vectorSub(wd->vVertex2,wd->vOrigin,wd->vVertex2);
#endif /* !QUATERNIONS */
		vectorCross(wd->vVertex1,wd->vVertex2,wd->vOrient);
		dTmp = vectorMag(wd->vOrient);
		assert(dTmp > 0.0); /* error condition already checked for */
		vectorScale(wd->vOrient,1.0/dTmp, wd->vOrient);
		}
	    if ((wd->iType==WallCylinderFinite) || (wd->iType==WallDisk) || (wd->iType==WallShell)) {
		/* Rotate the Orientation */
		vectorSub(wd->vOrigin, wd->vCoG, wd->vOrigin); //converting vectors relative to space to vector relative to body frame
		vectorRotate(wd->vOrigin, wd->vPxj, wd->vRot[0]*dDelta);
		vectorRotate(wd->vOrigin, wd->vPyj, wd->vRot[1]*dDelta);
		vectorRotate(wd->vOrigin, wd->vPzj, wd->vRot[2]*dDelta);
		vectorAdd(wd->vOrigin, wd->vCoG, wd->vOrigin);

		vectorRotate(wd->vOrient, wd->vPxj, wd->vRot[0]*dDelta);
		vectorRotate(wd->vOrient, wd->vPyj, wd->vRot[1]*dDelta);
		vectorRotate(wd->vOrient, wd->vPzj, wd->vRot[2]*dDelta);
		}
	    }
#endif /* DEM_WALLS_REACT */
	}
    in.WP = *WP; /* struct copy */
    in.bUpdatePos = 1;
    in.dDelta = dDelta;
    pstWallsUpdateStuckParticles(msr->pst,&in,sizeof(in),NULL,NULL);
    }

#ifdef DEM_WALLS_REACT
#undef MIN_STEP
#undef EPS
#endif /* DEM_WALLS_REACT */

#endif /* WALLS */

#ifdef SPECIAL_PARTICLES

void msrSpecialGetData(MSR msr,const char achFilenameWithQuotes[])
{
	FILE *fp;
	struct parameters *p;
	SPECIAL_PARTICLE_DATA *s;
	char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
	int i,n=0;

	assert(msr != NULL && achFilenameWithQuotes != NULL);
	p = &msr->param;
	_msrStripQuotes(achFilenameWithQuotes,achFilename);
	if (!strlen(achFilename)) {
		p->nSpecial = 0;
		return;
		}
	_msrMakePath(p->achDataSubPath,achFilename,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
	if (!(fp = fopen(achFilename,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",achFilename);
		goto abort;
		}
	if (fscanf(fp,"%i",&p->nSpecial) != 1) {
		(void) fprintf(stderr,"Expected no. special particles in \"%s\"\n",achFilename);
		goto abort;
		}
	if (p->nSpecial <= 0) {
		(void) fprintf(stderr,"Invalid no. special particles in \"%s\"\n",achFilename);
		goto abort;
		}
	if (p->nSpecial > MAX_NUM_SPECIAL_PARTICLES) {
		(void) fprintf(stderr,"Number of special particles (%i) exceeds maximum (%i)\n",
					   p->nSpecial,MAX_NUM_SPECIAL_PARTICLES);
		goto abort;
		}
	s = p->sSpecialData;
	for (i=0;i<p->nSpecial;i++) {
		if (fscanf(fp,"%i%i",&p->iSpecialId[i],&s[i].iType) != 2) {
			(void) fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
						   achFilename,i);
			goto abort;
			}
		if (p->iSpecialId[i] < -1) {
			(void) fprintf(stderr,"Invalid original index (%i) in \"%s\", entry %i\n",
						   p->iSpecialId[i],achFilename,i);
			goto abort;
			}
		if (p->iSpecialId[i] == -1) {
			if (!msr->param.bHeliocentric) {
				puts("ERROR: Negative special particle index requires non-inertial frame");
				goto abort;
				}
			if (++n > 1) {
				puts("ERROR: Only one special particle can have index -1");
				goto abort;
				}
			}
		if (!((s[i].iType & SPECIAL_OBLATE) || (s[i].iType & SPECIAL_GR) ||
			 (s[i].iType & SPECIAL_FORCE)|| (s[i].iType & SPECIAL_NOGHOSTPERT))) {
			(void) fprintf(stderr,"Invalid data type (%i) in \"%s\", entry %i\n",
						   s[i].iType,achFilename,i);
			goto abort;
			}
		if (s[i].iType & SPECIAL_OBLATE) {
			if (fscanf(fp,"%lf%lf%lf%lf%lf%lf",&s[i].oblate.dRadEq,&s[i].oblate.J2,
					   &s[i].oblate.J4,&s[i].oblate.p[0],&s[i].oblate.p[1],
					   &s[i].oblate.p[2]) != 6) {
				(void) fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
							   achFilename,i);
				goto abort;
				}
			}
		if (s[i].iType & SPECIAL_GR) {
			(void) fprintf(stderr,"GR not supported yet.\n");
			goto abort;
			}
		if (s[i].iType & SPECIAL_FORCE) {
			if (fscanf(fp,"%lf",&s[i].force.dMag) != 1) {
				(void) fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
							   achFilename,i);
				goto abort;
				}
			}
		if (s[i].iType & SPECIAL_NOGHOSTPERT) {
		        /* Do Nothing */
		        }
		}
	(void) fclose(fp);
	return;
 abort:
	if (fp) (void) fclose(fp);
	_msrExit(msr,1);
	}

#endif /* SPECIAL_PARTICLES */

#ifdef DEM_COHESION

void msrCohesionGetData(MSR msr,const char achFilenameWithQuotes[])
{
    FILE *fp;
    DEM_PARAMS *DP;
    DEM_COHESION_VAR_PARAMS *cvp;
    char achFilename[MAXPATHLEN],achTmp[MAXPATHLEN];
    int i;

    assert(msr != NULL && achFilenameWithQuotes != NULL);
    DP = &msr->param.FO.DP;
    cvp = &DP->sCohesionVarParams;
    _msrStripQuotes(achFilenameWithQuotes,achFilename);
    if (!strlen(achFilename)) {
        cvp->nCohesiveCoeff = 0;
        return;
        }
    _msrMakePath(msr->param.achDataSubPath,achFilename,achTmp);
    _msrMakePath(msr->lcl.pszDataPath,achTmp,achFilename);
    if (!(fp = fopen(achFilename,"r"))) {
        fprintf(stderr,"Unable to open \"%s\"\n",achFilename);
        goto abort;
        }
    if (fscanf(fp,"%i",&cvp->nCohesiveCoeff) != 1) {
        fprintf(stderr,"Expected number of cohesive coefficients in \"%s\"\n",achFilename);
        goto abort;
        }
    if (cvp->nCohesiveCoeff <= 0) {
        fprintf(stderr,"Invalid number of cohesive coefficients in \"%s\"\n",achFilename);
        goto abort;
        }
    if (fscanf(fp,"%i",&cvp->nConnetCoeff) != 1) {
        fprintf(stderr,"Rule to calculate cohesive coefficient between different colors \"%s\"\n",achFilename);
        goto abort;
        }
    if (cvp->nConnetCoeff < 0 && cvp->nConnetCoeff > 3) {
        fprintf(stderr,"Invalid rule in \"%s\"\n",achFilename);
        goto abort;
        }
    cvp->iColor = (int *) malloc(cvp->nCohesiveCoeff*sizeof(int));
    assert(cvp->iColor != NULL);
    cvp->dCohesiveCoeff = (double *) malloc(cvp->nCohesiveCoeff*sizeof(double));
    assert(cvp->dCohesiveCoeff != NULL);
    for (i=0;i<cvp->nCohesiveCoeff;i++) {
        if (fscanf(fp,"%i%lf",&cvp->iColor[i],&cvp->dCohesiveCoeff[i]) != 2) {
            fprintf(stderr,"Invalid/missing data in \"%s\" (entry %i)\n",
                    achFilename,i);
            fprintf(stderr,"Make sure the first line gives the correct number of vectors to use.\n");
            goto abort;
            }
        if (cvp->dCohesiveCoeff[i] < 0.0) {
            fprintf(stderr,"Invalid cohesive coefficient in \"%s\", entry %i: %g\n"
                    "(DEM cohesive coefficient must not be negative.)\n",
                    achFilename,i,cvp->dCohesiveCoeff[i]);
            goto abort;
            }
        }
    fclose(fp);
    return;
 abort:
    if (fp) fclose(fp);
    _msrExit(msr,1);
    }

#endif /* DEM_COHESION */

#ifdef DEM_TIDAL_SPACE

void msrDEMTidalGetData(MSR msr,DEM_TIDAL *d)
{
	struct inDEMTidalFindPlanet inFindPlanet;
	struct outDEMTidalFindPlanet outFindPlanet;
	struct inDEMTidalFindMarker inFindMarker;
	struct outDEMTidalFindMarker outFindMarker;

	Vector vSpinDot;
	int i;

	/* find the planet and store its info */

	inFindPlanet.iColor = PLANET_COLOR;
	pstDEMTidalFindPlanet(msr->pst,&inFindPlanet,sizeof(inFindPlanet),&outFindPlanet,NULL);
	assert(outFindPlanet.bFound);
	vectorCopy(outFindPlanet.dPos,d->vPlanetPos);
	vectorCopy(outFindPlanet.dVel,d->vPlanetVel);

	/* store aggregate center-of-mass data */

	assert(msr->pAggs[0].bAssigned); /* the asteroid should be agg 0 */
	vectorCopy(msr->pAggs[0].r_com,d->vAggPos);
	vectorCopy(msr->pAggs[0].v_com,d->vAggVel);
	vectorCopy(msr->pAggs[0].a_com,d->vAggAcc);
	/* following converts spin & torque per unit mass from body frame to space frame */
	matrixTransform(msr->pAggs[0].lambda,msr->pAggs[0].omega,d->vAggSpin);
	vectorSet(vSpinDot,
			  msr->pAggs[0].torque[0]/msr->pAggs[0].moments[0],
			  msr->pAggs[0].torque[1]/msr->pAggs[0].moments[1],
			  msr->pAggs[0].torque[2]/msr->pAggs[0].moments[2]);
	matrixTransform(msr->pAggs[0].lambda,vSpinDot,d->vAggSpinDot);

	/* find the markers and store their info */

	for (i=0;i<3;i++) {
		switch (i) {
		case 0:
			inFindMarker.iOrder = IORDER_MARKER1;
			break;
		case 1:
			inFindMarker.iOrder = IORDER_MARKER2;
			break;
		case 2:
			inFindMarker.iOrder = IORDER_MARKER3;
			break;
		default:
			assert(0); /* shouldn't get here! */
			}
		pstDEMTidalFindMarker(msr->pst,&inFindMarker,sizeof(inFindMarker),&outFindMarker,NULL);
		assert(outFindMarker.bFound);
		switch (i) {
		case 0:
			vectorCopy(outFindMarker.dPos,d->vMarker1Pos);
			vectorCopy(outFindMarker.dVel,d->vMarker1Vel);
			break;
		case 1:
			vectorCopy(outFindMarker.dPos,d->vMarker2Pos);
			vectorCopy(outFindMarker.dVel,d->vMarker2Vel);
			break;
		case 2:
			vectorCopy(outFindMarker.dPos,d->vMarker3Pos);
			vectorCopy(outFindMarker.dVel,d->vMarker3Vel);
			break;
		default:
			assert(0);
			}
		}
	}

#endif /* DEM_TIDAL_SPACE */

#ifdef DEM_TIDAL_LOCAL

#define DEM_TIDAL_INIT_BUFF_SIZE 1000

void msrDEMTidalReadData(MSR msr)
{
	/* read in acceleration data */

	DEM_TIDAL **D = &msr->DEMTidalData;

	FILE *fp;
	int nBuffSize = DEM_TIDAL_INIT_BUFF_SIZE,i;

	fp = fopen(DEM_TIDAL_FILENAME,"r");
	assert(fp != NULL);
	*D = (DEM_TIDAL *) malloc(DEM_TIDAL_INIT_BUFF_SIZE*sizeof(DEM_TIDAL));
	assert(*D != NULL);
	i = 0;
	while (!feof(fp)) {
		fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			   &((*D)[i]).dTime,
			   &((*D)[i]).vPlanetPos[0],&((*D)[i]).vPlanetPos[1],&((*D)[i]).vPlanetPos[2],
			   &((*D)[i]).vPlanetVel[0],&((*D)[i]).vPlanetVel[1],&((*D)[i]).vPlanetVel[2],
			   &((*D)[i]).vAggPos[0],&((*D)[i]).vAggPos[1],&((*D)[i]).vAggPos[2],
			   &((*D)[i]).vAggVel[0],&((*D)[i]).vAggVel[1],&((*D)[i]).vAggVel[2],
			   &((*D)[i]).vAggSpin[0],&((*D)[i]).vAggSpin[1],&((*D)[i]).vAggSpin[2],
			   &((*D)[i]).vAggAcc[0],&((*D)[i]).vAggAcc[1],&((*D)[i]).vAggAcc[2],
			   &((*D)[i]).vAggSpinDot[0],&((*D)[i]).vAggSpinDot[1],&((*D)[i]).vAggSpinDot[2],
			   &((*D)[i]).vMarker1Pos[0],&((*D)[i]).vMarker1Pos[1],&((*D)[i]).vMarker1Pos[2],
			   &((*D)[i]).vMarker1Vel[0],&((*D)[i]).vMarker1Vel[1],&((*D)[i]).vMarker1Vel[2],
			   &((*D)[i]).vMarker2Pos[0],&((*D)[i]).vMarker2Pos[1],&((*D)[i]).vMarker2Pos[2],
			   &((*D)[i]).vMarker2Vel[0],&((*D)[i]).vMarker2Vel[1],&((*D)[i]).vMarker2Vel[2],
			   &((*D)[i]).vMarker3Pos[0],&((*D)[i]).vMarker3Pos[1],&((*D)[i]).vMarker3Pos[2],
			   &((*D)[i]).vMarker3Vel[0],&((*D)[i]).vMarker3Vel[1],&((*D)[i]).vMarker3Vel[2]);
		if (++i == nBuffSize) {
			nBuffSize <<= 1; /* double buffer size */
			*D = (DEM_TIDAL *) realloc(*D,nBuffSize*sizeof(DEM_TIDAL));
			assert(*D != NULL);
			}
		}
	fclose(fp);
	printf("DEMTidal: # entries = %i\n",msr->nDEMTidalDataEntries = i);
	assert(msr->nDEMTidalDataEntries > 0);
	}

#undef DEM_TIDAL_INIT_BUFF_SIZE

void msrDEMTidal(MSR msr,double dTime)
{
	static int i = 0;

	struct inDEMTidal in;

	while (i < msr->nDEMTidalDataEntries && msr->DEMTidalData[i].dTime < dTime)
		++i;

	in.d0 = msr->DEMTidalData[0]; /* struct copy */
	in.d = msr->DEMTidalData[i]; /* ditto */ /* REALLY NEED TO INTERPOLATE HERE! */
	pstDEMTidal(msr->pst,&in,sizeof(in),NULL,NULL);
	}

#endif /* DEM_TIDAL_LOCAL */

/* Note if dDataOut is NULL it just counts the number of valid input lines */
int msrReadASCII(MSR msr, char *extension, int nDataPerLine, double *dDataOut)
{
	char achFile[PST_FILENAME_SIZE];
	char ach[PST_FILENAME_SIZE];
	LCL *plcl = &msr->lcl;
	FILE *fp;
	int i,ret;
	char achIn[160];
	double *dData;

	if (dDataOut == NULL) 
		dData = malloc(sizeof(double)*nDataPerLine);
	else
		dData = dDataOut;
	
	assert(nDataPerLine > 0 && nDataPerLine <= 10);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	achFile[0] = '\0';
	snprintf(achFile,PST_FILENAME_SIZE,"%s/%s.%s",msr->param.achDataSubPath,
			 msr->param.achOutName, extension);
	/*
	 ** Add local Data Path.
	 */
	if (plcl->pszDataPath) {
		strcpy(ach,achFile);
		snprintf(achFile,PST_FILENAME_SIZE,"%s/%s",plcl->pszDataPath,ach);
		}
	fp = fopen(achFile,"r");
	if (!fp) {
		if (msr->param.bVWarnings)
			printf("WARNING: Could not open .%s input file:%s\n",extension,achFile);
		return 0;
		}

	i = 0;
	while (1) {
		if (!fgets(achIn,160,fp)) goto Done;
		switch (nDataPerLine) {
		case 1:
			ret = sscanf(achIn,"%lf",dData); 
			break;
		case 2:
			ret = sscanf(achIn,"%lf %lf",dData,dData+1); 
			break;
		case 3:
			ret = sscanf(achIn,"%lf %lf %lf",dData,dData+1,dData+2); 
			break;
		case 4:
			ret = sscanf(achIn,"%lf %lf %lf %lf",dData,dData+1,dData+2,dData+3); 
			break;
		case 5:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4); 
			break;
		case 6:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4,dData+5); 
			break;
		case 7:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6); 
			break;
		case 8:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7); 
			break;
		case 9:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8); 
			break;
		case 10:
			ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
						 dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8,dData+9); 
			break;
		default:
			ret = EOF;
			assert(0);
			}
		if (ret != nDataPerLine) goto Done;
		++i;
		if (dDataOut != NULL) dData += nDataPerLine;
		}
 Done:
	fclose(fp);
	if (dDataOut != NULL && msr->param.bVDetails) printf("Read %i lines from %s\n",i,achFile);
	if (dDataOut == NULL) free(dData);
	return i;
	}

int msrSetTypeFromFile(MSR msr, char *file, int iSetMask)
{
	FILE *fp;
	struct inSetTypeFromFile in;
	struct outSetTypeFromFile out;

	in.iSetMask = iSetMask;
	in.biGasOrder = 1; /* Set from Parent Gas iOrder for stars? */
	assert(strlen(file) < PST_SETTYPEFROMFILEMAXLEN);
	strcpy( &(in.file[0]), file );

	fp = fopen( file, "r" );
	if (!fp) {
	  fprintf(stderr,"ERROR: Could not open iOrder list file:%s\n",file);
	  assert(0);
	  }
	fclose(fp);

	pstSetTypeFromFile(msr->pst,&in,sizeof(in),&out,NULL);
	
	if (msr->param.bVDetails) printf("%d iOrder numbers read.  %d direct and %d Gas Parent iOrder photogenic particles selected.",out.niOrder,out.nSet,out.nSetiGasOrder);
	
	return out.nSet+out.nSetiGasOrder;
	}

FILE *LogTimingInit( MSR msr, const char *fileflag )
    {
    char achFile[MAXPATHLEN];
    FILE *fpLogTiming;
    
    if (!msr->param.bLogTiming) return NULL;

    snprintf(achFile,MAXPATHLEN,"%s.timing",msrOutName(msr));
    fpLogTiming = fopen(achFile,fileflag);
    assert(fpLogTiming != NULL);
    setbuf(fpLogTiming,(char *) NULL); /* no buffering */

    msr->iRungStat = 0;
    msr->RungStat = (struct RungData *) malloc(sizeof(struct RungData)*msr->param.iMaxRung);
    assert( msr->RungStat != NULL );

    return fpLogTiming;
    }

void LogTimingZeroCounters( MSR msr ) 
    {
    int i,j;
    struct RungData *r;

    if (!msr->param.bLogTiming) return;
    for (i=0;i<msr->param.iMaxRung;i++) {
	r = &(msr->RungStat[i]);
	r->nPart = 0;
	r->nPartMin = 2*msr->N;
	r->nPartMax = 0;
	r->nUses = 0;
	r->nPartTot = 0;
	r->nPartMinTot = 2*msr->N;
	r->nPartMaxTot = 0;
	r->nUsesTot = 0;
	for (j=0;j<TIMING_N;j++) {
	    r->nCall[j] = 0;
	    r->t[j] = 0;
	    r->tTot[j] = 0;
	    }
	}
    }

void LogTimingSetRung ( MSR msr, int iRung )
    {
    msr->iRungStat = iRung;
    if (msr->param.bLogTiming)
		printf("Timing: rung: %d set\n",msr->iRungStat);
    }

void LogTimingSetN( MSR msr, int n ) 
    {
    struct RungData *r;

    if (!msr->param.bLogTiming) return;
    r = &(msr->RungStat[msr->iRungStat]);
    r->nPart += n;
    r->nUses ++;
    if (n < r->nPartMin) r->nPartMin = n;
    if (n > r->nPartMax) r->nPartMax = n;
    printf("Timing: rung: %d set n %d\n",msr->iRungStat,n);
    }


void LogTimingOutput( MSR msr, FILE *fpLogTiming, double dTime, int bAll )
    {
    int i,j,nStep;
    double f,Stept[TIMING_N],SteptTot[TIMING_N];
    struct RungData *r;
    
    if (!msr->param.bLogTiming) return;

    /* Note: some entries (e.g. Total for now have zero for nCall! ) */
#ifdef TIMINGDEBUG
    printf("Timing: Output and Zero step timers\n");
#endif

    for (j=0;j<TIMING_N;j++) {
	Stept[j] = 0;
	SteptTot[j] = 0;
	}

    for (i=0;i<msr->param.iMaxRung;i++) {
	r = &(msr->RungStat[i]);
	r->t[0] = 0;
	for (j=1;j<TIMING_N;j++) {
	    r->t[0] += r->t[j];
	    Stept[j] += r->t[j];
	    }
	Stept[0] += r->t[0];
	}

    for (i=0;i<msr->param.iMaxRung;i++) {
	r = &(msr->RungStat[i]);
	r->nPartTot += r->nPart;
	if (r->nPartMin < r->nPartMinTot)  r->nPartMinTot = r->nPartMin;
	if (r->nPartMax > r->nPartMaxTot)  r->nPartMaxTot = r->nPartMax;
	r->nUsesTot += r->nUses;
	for (j=0;j<TIMING_N;j++) {
	    r->tTot[j] += r->t[j];
	    SteptTot[j] += r->tTot[j];
	    }
	}

    nStep = msr->RungStat[0].nUsesTot;
    fprintf( fpLogTiming,"%d %e %e, %f %f\n",
	     nStep,dTime,1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,Stept[0],SteptTot[0]/nStep);

    fprintf( fpLogTiming,"STEP TOT    : " );
    for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",Stept[j] );
    fprintf( fpLogTiming,"\n" );
    if (bAll) {
	fprintf( fpLogTiming,"STEP TOT AVG: " );
	f = 1.0/nStep;
	for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",SteptTot[j]*f );
	fprintf( fpLogTiming,"\n" );
	}
    
    if (msr->param.bLogTimingSubStep) {
	for (i=0;i<msr->param.iMaxRung;i++) {
	    r = &(msr->RungStat[i]);
	    if (r->nUses) {
		f = 1.0/r->nUses;
		fprintf( fpLogTiming,"%d: %lld, %f %lld %lld,   ", i,r->nUses, r->nPart*f, r->nPartMin, r->nPartMax);
		for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",r->t[j]*f );
		fprintf( fpLogTiming,"\n");
		}
	    }
	}
    
    if (msr->param.bLogTimingStep) {
	for (i=0;i<msr->param.iMaxRung;i++) {
	    r = &(msr->RungStat[i]);
	    if (r->nUses) {
		f = 1.0/r->nUses;
		fprintf( fpLogTiming,"%d: %lld, %f %lld %lld,  S ", i,r->nUses, r->nPart*f, r->nPartMin, r->nPartMax);
		for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",r->t[j] );
		fprintf( fpLogTiming,"\n");
		}
	    }
	}

    if (msr->param.bLogTimingSubStepTot || bAll) {
	for (i=0;i<msr->param.iMaxRung;i++) {
	    r = &(msr->RungStat[i]);
	    if (r->nUsesTot) {
		f = 1.0/r->nUsesTot;
		fprintf( fpLogTiming,"%d: %lld, %f %lld %lld,  T ", i,r->nUsesTot, r->nPartTot*f, r->nPartMinTot, r->nPartMaxTot);
		for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",r->tTot[j]*f );
		fprintf( fpLogTiming,"\n");
		}
	    }
	}
    
    if (msr->param.bLogTimingStepTot || bAll) {
	for (i=0;i<msr->param.iMaxRung;i++) {
	    r = &(msr->RungStat[i]);
	    if (r->nUsesTot) {
		f = 1.0/r->nUsesTot;
		fprintf( fpLogTiming,"%d: %lld, %f %lld %lld, ST ", i,r->nUsesTot, r->nPartTot*f, r->nPartMinTot, r->nPartMaxTot);
		f = 1.0/nStep;
		for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",r->tTot[j]*f );
		fprintf( fpLogTiming,"\n");
		}
	    }
	}
    

    for (i=0;i<msr->param.iMaxRung;i++) {
	r = &(msr->RungStat[i]);
	r->nPart = 0;
	r->nPartMin = 2*msr->N;
	r->nPartMax = 0;
	r->nUses = 0;
	for (j=0;j<TIMING_N;j++) {
	    r->t[j] = 0;
	    }
	}
    }

void LogTimingFinish( MSR msr, FILE *fpLogTiming, double dTime )
    {
    double f;
    struct RungData *r;
    int i,j;

    if (!msr->param.bLogTiming) return;

    fprintf( fpLogTiming,"Log Complete\n");
    for (i=0;i<msr->param.iMaxRung;i++) {
	r = &(msr->RungStat[i]);
	if (r->nUsesTot) {
	    f = 1.0/r->nUsesTot;
	    fprintf( fpLogTiming,"%d: Calls/Use ", i );
	    for (j=0;j<TIMING_N;j++) fprintf( fpLogTiming,"%f ",r->nCall[j]*f );
	    fprintf( fpLogTiming,"\n");
	    }
	}
    

    /* Note this corrupts the Data for a final output but we are deallocating it */
    LogTimingOutput( msr, fpLogTiming, dTime, 1 );

    fclose(fpLogTiming);
    
    free( msr->RungStat );
    }
