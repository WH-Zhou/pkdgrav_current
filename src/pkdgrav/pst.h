#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#include "floattype.h"
#include "dumpframe.h"

#ifdef COLLISIONS
#include "collision.h"
#endif

#ifdef AGGS
#include "aggs.h"
#endif

#ifdef RUBBLE_ZML
#include "rubble.h"
#elif defined(COLLMOD_ZML)
#include "collmod.h"
#endif

#ifdef STARFORM
#include "starform.h"
#include "feedback.h"
#include "supernova.h"
#endif

#ifdef DEM_WALLS_REACT
#include "dem.h"
#endif

typedef struct lclBlock {
	char *pszDataPath;
	PKD	pkd;
	int nPstLvl;
	int iWtFrom;
	int iWtTo;
	int iPart;
	int iOrdSplit;
	FLOAT fSplit;
	FLOAT fWtLow;
	FLOAT fWtHigh;
	FLOAT fLow;
	FLOAT fHigh;
	int nWriteStart;
	int nDarkWriteStart;
	int nGasWriteStart;
	int nStarWriteStart;
	} LCL;

typedef struct pstContext {
	struct pstContext *pstLower;
	MDL mdl;
	LCL *plcl;
	int idSelf;
	int idUpper;
	int nLeaves;
	int nLower;
	int nUpper;
	int iLvl;
	BND bnd;
	BND bndActive;
	BND bndTreeActive;
	int iSplitDim;
	int iOrdSplit;
	FLOAT fSplit;
	FLOAT fSplitInactive;
	int nTotal;
	int nDark;
	int nGas;
	int nStar;
	/*
	 ** The PST node is also a valid cell for the tree.
	 */
	KDN kdn;
	} * PST;


#define PST_FILENAME_SIZE	512

enum pst_service {
  PST_SRV_STOP,
  PST_SETADD,
  PST_LEVELIZE,
  PST_READTIPSY,
  PST_DOMAINDECOMP,
  PST_CALCBOUND,
  PST_GASWEIGHT,
  PST_RUNGDDWEIGHT,
  PST_WEIGHT,
  PST_WEIGHTWRAP,
  PST_FREESTORE,
  PST_COLREJECTS,
  PST_SWAPREJECTS,
  PST_DOMAINCOLOR,
  PST_COLORDREJECTS,
  PST_DOMAINORDER,
  PST_LOCALORDER,
  PST_OUTARRAY,
  PST_OUTVECTOR,
  PST_OUTNCVECTOR,
  PST_WRITETIPSY,
  PST_BUILDTREE,
  PST_SMOOTH,
  PST_GRDRAGGETSUNACCEL,
  PST_GRINTEGRATECLOSEPARTICLES,
  PST_GRCORRECTASUN,
  PST_GRAVITY,
  PST_GRAVEXTERNAL,
  PST_CALCEANDL,
  PST_CALCEANDLEXT,
  PST_DRIFT,
  PST_KICK,
  PST_READCHECK,
  PST_WRITECHECK,
  PST_SETSOFT,
  PST_PHYSICALSOFT,
  PST_PREVARIABLESOFT,
  PST_POSTVARIABLESOFT,
  PST_SETTOTAL,
  PST_SETTOTALS,
  PST_CALCCELL,
  PST_COLCELLS,
  PST_DISTRIBCELLS,
  PST_CALCROOT,
  PST_DISTRIBROOT,
  PST_ONENODEREADINIT,
  PST_SWAPALL,
  PST_MASSCHECK,
  PST_GETMASSCHANGE,
  PST_MASSMETALSENERGYCHECK,
  PST_ACTIVETYPEORDER,
  PST_ACTIVEORDER,
  PST_SETRUNG,
  PST_ACTIVERUNG,
  PST_CURRRUNG,
  PST_GRAVSTEP,
  PST_ACCELSTEP,
  PST_DENSITYSTEP,
  PST_RUNGSTATS,
  PST_GETMAP,
  PST_RANDSEEDGENERATOR,
  PST_COOLVELOCITY,
  PST_RESETTOUCHRUNG,
  PST_ACTIVEEXACTTYPE,
  PST_ACTIVETYPE,
  PST_SETTYPE,
  PST_RESETTYPE,
  PST_COUNTTYPE,
  PST_ACTIVEMASKRUNG,
  PST_ACTIVETYPERUNG,
  PST_SETPARTICLETYPES,
  PST_SOUGHTPARTICLELIST,
  PST_COOLUSINGPARTICLELIST,
  PST_SETTYPEFROMFILE,
  PST_MARKSMOOTH,
  PST_RESMOOTH,
  PST_INITACCEL,
  PST_DTTORUNG,
  PST_INITDT,
  PST_ORDWEIGHT,
  PST_SETWRITESTART,
  PST_SETNCWRITESTART,
  PST_COLNPARTS,
  PST_NEWORDER,
  PST_GETNPARTS,
  PST_SETNPARTS,
  PST_UPDATEUDOT,
  PST_UPDATESHOCKTRACKER,
  PST_GETGASPRESSURE,
  PST_LOWERSOUNDSPEED,
  PST_INITENERGY,
  PST_KICKVPRED,
  PST_KICKRHOPRED,
  PST_SPHCURRRUNG,
  PST_RANDOMVELOCITIES,
  PST_DENSCHECK,
#ifdef SPECIAL_PARTICLES
  PST_GETSPECIALPARTICLES,
  PST_DOSPECIALPARTICLES,
#endif
#ifdef COLLISIONS
  PST_NUMREJECTS,
  PST_READSS,
  PST_WRITESS,
  PST_ADDUNIFGRAV,
  PST_NEXTCOLLISION,
  PST_GETCOLLIDERINFO,
  PST_DOCOLLISION,
  PST_RESETCOLLIDERS,
#ifdef RORY_EXTENSION
  PST_SETBALL,
  PST_FINDTIGHESTBINARY,
  PST_MERGEBINARY,
  PST_FINDLM,
  PST_GETNEIGHBORS,
  PST_MOVEPART,
#endif /* RORY_EXTENSION */
#endif /* COLLISIONS */
#ifdef JOHNNY
  PST_FINDESCAPERS,
  PST_CHECKFORKEPLER,
#endif /* JOHNNY */
#ifdef AGGS
  PST_AGGSCACHEINIT,
  PST_DOMAINAGGORDER,
  PST_LOCALAGGORDER,
  PST_AGGSFIND,
  PST_AGGSCOUNTPART,
  PST_AGGSDELETE,
  PST_AGGSMERGE,
  PST_AGGSBACKDRIFT,
#ifdef FAST_AGGS
  PST_AGGSGETREF,
  PST_AGGSGETPERIODICCOM,
  PST_AGGSSETUNWRAPPED,
#endif /* FAST_AGGS */
  PST_AGGSGETCOM,
  PST_AGGSGETAXESANDSPIN,
  PST_AGGSSETBODYPOS,
  PST_AGGSSETSPACEPOS,
  PST_AGGSSETSPACEVEL,
  PST_AGGSSETSPACESPINS,
#ifndef DEM
  PST_AGGSGETACCEL,
  PST_AGGSCHECKSTRESS,
  PST_AGGSGETTORQUE,
#else /* DEM */
  PST_AGGSGETACCELANDTORQUEDEM,
  PST_AGGSSETMASSDEM,
  PST_AGGWEIGHT,
  PST_COLAGGREJECTS,
#endif /* DEM */
#endif /* AGGS */
#ifdef WALLS
  PST_WALLSUPDATESTUCKPARTICLES,
#endif /* WALLS */
#ifdef SPRINGS
  PST_BREAKSPRINGS,
#endif /* SPRINGS */
#ifdef DEM
  PST_DEMZEROSMALLMOTIONS,
  PST_DEMSTATS,
  PST_DEMSETBALL,
  PST_DEMGETBALLRADIUS,
  PST_DEMFIXEDBALLNEXTLARGEST,
#endif /* DEM */
#ifdef SPINUP
  PST_SPINUPGETAGG,
  PST_SPINUPALIGN,
  PST_SPINUP,
  /* following for SPINUP_WITH_AGGS */
  PST_SPINUPACCEL,
  PST_SPINUPRESETAGGS,
  /* following for DEM, WALLS, DEM_WALLS_REACT */
#endif /* SPINUP */
#ifdef DEM_WALLS_REACT
  PST_DEMWALLSREACT,
#endif /* DEM_WALLS_REACT */
#ifdef CHARGE
  PST_CHARGEZGETMOMENTS,
  PST_CHARGEZAPPLYMOMENTS,
#endif /* CHARGE */
#ifdef DEM_TIDAL_SPACE
  PST_DEMTIDALFINDPLANET,
  PST_DEMTIDALFINDMARKER,
#endif /* DEM_TIDAL_SPACE */
#ifdef DEM_TIDAL_LOCAL
  PST_DEMTIDAL,
#endif
#ifdef SLIDING_PATCH
  PST_KICKPATCH,
  PST_RANDAZWRAP,
#endif /* SLIDING_PATCH */
#ifdef AGGS_IN_PATCH
  PST_AGGSINPATCHGETREF,
  PST_AGGSINPATCHGETUNWRAPPED,
  PST_AGGSINPATCHOFFSET,
#endif /* AGGS_IN_PATCH */
#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)
  PST_DUSTBINSINCLMAX,
  PST_DUSTBINSINCLAVG,
  PST_DUSTBINSGETMASS,
  PST_DUSTBINSAPPLY,
#endif /* RUBBLE_ZML || COLLMOD_ZML */
#ifdef RUBBLE_ZML
  PST_RUBBLERESETCOLFLAG,
  PST_RUBBLECHECKFORKDKRESTART,
  PST_RUBBLESTEP,
  PST_RUBCLEANUP,
  PST_RUBINTERPCLEANUP,
#endif /* RUBBLE_ZML */
#ifdef COLLMOD_ZML
  PST_COLLMODRESETCOLFLAG,
  PST_COLLMODCHECKFORKDKRESTART,
  PST_COLLMODSTEP,
#endif /* COLLMOD_ZML */
#ifdef ORIGIN_HISTOGRAM
  PST_INITIALIZEORIGINBINS,
#endif /* ORIGIN_HISTOGRAM */
  /* following for SPH, etc. */
  PST_SPHSTEP,
  PST_SPHVISCOSITYLIMITER,
  PST_INITCOOLING,
  PST_COOLTABLEREAD,
  PST_GROWMASS,
  PST_CLEARTIMER,
  PST_MASSINR,
  PST_ROTBARINIT,
  PST_BALLMAX,
  PST_FORMSTARS,
  PST_FEEDBACK,
  PST_SIMPLESTARFORM,
  PST_SSFCREATESTARS,
  PST_DUMPFRAME,
  PST_DUMPVOXEL,
  PST_COM,
  PST_COMBYTYPE,
  PST_OLDESTSTAR,
  PST_SETSINK,
};

void pstAddServices(PST,MDL);
void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

/* PST_SETADD */
struct inSetAdd {
	int id;
	};
void pstSetAdd(PST,void *,int,void *,int *);

/* PST_LEVELIZE */
struct inLevelize {
	int iLvl;
	};
void pstLevelize(PST,void *,int,void *,int *);

/* PST_READTIPSY */
struct inReadTipsy {
	int nFileStart;
	int nFileEnd;
	int nDark;	
	int nGas;
	int nStar;
	int iOrder;
	float fExtraStore;
	FLOAT fPeriod[3];
	int bStandard;
	int iReadIOrder;
	double dvFac;
	double dTuFac;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadTipsy(PST,void *,int,void *,int *);

/* PST_DOMAINDECOMP */
struct inDomainDecomp {
    int bDoRootFind;
    int bDoSplitDimFind;
    };

void pstDomainDecomp(PST,void *,int,void *,int *);

/* PST_CALCBOUND */
struct outCalcBound {
	BND bnd;
	BND bndActive;
	BND bndTreeActive;
	BND bndBall;
	};
void pstCalcBound(PST,void *,int,void *,int *);

void pstGasWeight(PST,void *,int,void *,int *);

struct inRungDDWeight {
	int iMaxRung;
        double dWeight; 
        };

void pstRungDDWeight(PST,void *,int,void *,int *);

/* PST_WEIGHT */
struct inWeight {
	int iSplitDim;
	FLOAT fSplit;
	int iSplitSide;
	int ittr;
	int pFlag;
	};
struct outWeight {
	int nLow;
	int nHigh;
	FLOAT fLow;
	FLOAT fHigh;
	};
void pstWeight(PST,void *,int,void *,int *);
struct inWeightWrap {
	int iSplitDim;
	FLOAT fSplit;
	FLOAT fSplit2;
	int iSplitSide;
	int ittr;
	int pFlag;
	};
void pstWeightWrap(PST,void *,int,void *,int *);

/* PST_FREESTORE */
struct outFreeStore {
	int nFreeStore;
	};
void pstFreeStore(PST,void *,int,void *,int *);

/*
 ** This structure is used by reject collectors and SwapRejects
 */
typedef struct outReject {
	int id;
	int nRejects;
	int nSpace;
	int nLocal;
	} OREJ;

/* PST_COLREJECTS */
struct inColRejects {
	int iSplitDim;
	FLOAT fSplit;
	FLOAT fSplitInactive;
	int iSplitSide;
	};
void pstColRejects(PST,void *,int,void *,int *);

/* PST_SWAPREJECTS */
void pstSwapRejects(PST,void *,int,void *,int *);

/* PST_DOMAINCOLOR */
void pstDomainColor(PST,void *,int,void *,int *);

/* PST_COLORDREJECTS */
struct inColOrdRejects {
	int iOrdSplit;
	int iSplitSide;
	};
void pstColOrdRejects(PST,void *,int,void *,int *);

/* PST_DOMAINORDER */
struct inDomainOrder {
	int iMaxOrder;
	};
void pstDomainOrder(PST,void *,int,void *,int *);

/* PST_LOCALORDER */
void pstLocalOrder(PST,void *,int,void *,int *);

#ifdef FAST_AGGS
/* PST_LOCALAGGORDER */
void pstLocalAggOrder(PST,void *,int,void *,int *);	//JVD
#endif

/* PST_OUTARRAY */
struct inOutArray {
	char achOutFile[PST_FILENAME_SIZE];
        int nStart;
	int iType;
	int iBinaryOutput;
	int bStandard;
	};
void pstOutArray(PST,void *,int,void *,int *);

/*PST_OUTNCVECTOR*/
struct outNC {
    float min[3][3];
    float max[3][3];
    };
void pstOutNCVector(PST,void *,int,void *,int *);

/* PST_OUTVECTOR */
struct inOutput {
	char achOutFile[PST_FILENAME_SIZE];
	int iDim;
	int iType;
	int iBinaryOutput;
	int N;
	int bStandard;
	double duTFac;
	};
void pstOutVector(PST,void *,int,void *,int *);

/* PST_WRITETIPSY */
struct inWriteTipsy {
	int bStandard;
	double dvFac;
	double duTFac;
	int iGasModel;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteTipsy(PST,void *,int,void *,int *);

/* PST_BUILDTREE */
struct inBuildTree {
	int nBucket;
	int iOpenType;
	int iOrder;
	double dCrit;
	int bBinary;
        int bActiveOnly;
	int bTreeActiveOnly;
	int bGravity;
	};
struct outBuildTree {
	KDN kdn;
	};
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_SMOOTH */
struct inSmooth {
  int nSmooth;
  int bPeriodic;
  int bSymmetric;
  int iSmoothType;
  double dfBall2OverSoft2;
  SMF smf;
  };
struct outSmooth {
    int iSmoothFlags;  /* Warning Flags for need to smooth again, etc... */
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstSmooth(PST,void *,int,void *,int *);

#ifdef GR_DRAG
struct inGRDragGetSunAccel {
  double dSunMass;
  };
struct outGRDragGetSunAccel {
  double aSun[3];
  };
void pstGRDragGetSunAccel(PST,void *,int,void *,int *);

struct inGRIntegrateCloseParticles {
  double dSunMass,dDelta,dTime;
  };
struct outGRIntegrateCloseParticles {
  int nMerged;
  double dMergerMassLost;
  double aSunCorrection[3];
  };
void pstGRIntegrateCloseParticles(PST,void *,int,void *,int *);

struct inGRCorrectaSun {
  double aSunCorrection[3];
  };
void pstGRCorrectaSun(PST,void *, int,void *, int *);

#endif /*GR_DRAG*/

/* PST_GRAVITY */
struct inGravity {
	int nReps;
	int bPeriodic;
	int iOrder;
	int bEwald;
	int iEwOrder;
    int bDoSun;
	double dSunSoft;
	double dEwCut;
	double dEwhCut;
#ifdef COLLISIONS
  int bRepel;
  double dRepelFac;
#endif
#ifdef SLIDING_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif
	};
struct outGravity {
	int nActive;
	int nTreeActive;
    double aSun[3];
	double dPartSum;
	double dCellSum;
	double dSoftSum;
	double dFlop;
	/*	
	 ** Collected CPU time stats.
	 */
	double dWSum;
	double dWMax;
	double dWMin;
	double dISum;
	double dIMax;
	double dIMin;
	double dESum;
	double dEMax;
	double dEMin;
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstGravity(PST,void *,int,void *,int *);

/* PST_GRAVEXTERNAL */
struct inGravExternal {
    int bIndirect;
    int bDoSun;
	double dTime;
    double dSunMass;
	double dSunSoft;
    double aSun[3];
	/*
	 ** For external galaxy potential stuff
	 */
	int bLogHalo;
	int bHernquistSpheroid;
    int bNFWSpheroid;
        double dNFWm200;
        double dNFWr200;
        double dNFWconc;
        double dNFWsoft;
    int bElliptical;
    int bEllipticalDarkNFW;
	int bHomogSpheroid;
	int bBodyForce;
	int bMiyamotoDisk;
	int bTimeVarying;
	int bRotatingBar;
	double dRotBarAmp;
	double dRotBarPosAng;
	double dRotBarB5;
    FLOAT aCom[3];
    
#ifdef ROT_FRAME
	int bRotFrame;
	double dOmega;
	double dOmegaDot;
#endif
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
#ifdef SIMPLE_GAS_DRAG
	int bSimpleGasDrag;
	int iFlowOpt;
	int bEpstein;
	double dGamma;
#endif
	};
struct outGravExternal {
    double dAcc[3];
    double dTorque[3];
    };

void pstGravExternal(PST,void *,int,void *,int *);

/* PST_CALCEANDL */
struct outCalcEandL {
	double T;
	double U;
	double Eth;
	double L[3];
	};
void pstCalcEandL(PST,void *,int,void *,int *);

/* PST_CALCEANDLEXT */
struct inCalcEandLExt {
    int bHeliocentric;
	};
struct outCalcEandLExt {
	double dMass;
    double dSumMR[3];
    double dSumMV[3];
	double dPot;
	};
void pstCalcEandLExt(PST,void *,int,void *,int *);

/* PST_DRIFT */
struct inDrift {
	double dDelta;
	FLOAT fCenter[3];
	int bPeriodic;
	int bFandG;
	FLOAT fCentMass;
#ifdef SLIDING_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif
	};
void pstDrift(PST,void *,int,void *,int *);

/* PST_UPDATEUDOT */
struct inUpdateuDot {
	double duDelta;
	double dTime;	
	double z;
	int iGasModel;
	int bUpdateY;
	};
struct outUpdateuDot {
	double Time;
	double MaxTime;
	double SumTime;
	int nSum;
	};

void pstUpdateuDot(PST,void *,int,void *,int *);

/* PST_KICK */
struct inKick {
	double dvFacOne;
	double dvFacTwo;
	double dvPredFacOne;
	double dvPredFacTwo;
	double duDelta;
	double duPredDelta;
	double duDotLimit;
	int iGasModel;
	double z;
	};
struct outKick {
	double Time;
	double MaxTime;
	double SumTime;
	int nSum;
	};

void pstKick(PST,void *,int,void *,int *);

/* PST_READCHECK */
struct inReadCheck {
	int iVersion;
	int iOffset;
	int nFileStart;
	int nFileEnd;
	int nDark;
	int nGas;
	int nStar;
	int iOrder;
	float fExtraStore;
	FLOAT fPeriod[3];
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadCheck(PST,void *,int,void *,int *);

/* PST_WRITECHECK */
struct inWriteCheck {
	int iOffset;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteCheck(PST,void *,int,void *,int *);

/* PST_SETSOFT */
struct inSetSoft {
	double dSoft;
	};
void pstSetSoft(PST,void *,int,void *,int *);

#ifdef CHANGESOFT
/* PST_PHYSICALSOFT */
struct inPhysicalSoft {
        double dSoftMax;
        double dFac;
        int bSoftMaxMul;
        };
void pstPhysicalSoft(PST,void *,int,void *,int *);

/* PST_PREVARIABLESOFT */
struct inPreVariableSoft {
        int iVariableSoftType;
        };
void pstPreVariableSoft(PST,void *,int,void *,int *);

/* PST_POSTVARIABLESOFT */
struct inPostVariableSoft {
        double dSoftMax;
        int bSoftMaxMul;
        int iVariableSoftType;
        };
void pstPostVariableSoft(PST,void *,int,void *,int *);
#endif

/* PST_SETTOTAL */
struct outSetTotal {
	int nTotal;
	};
void pstSetTotal(PST,void *,int,void *,int *);

/* PST_SETTOTALS */
struct outSetTotals {
	int nGas;
	int nStar;
	int nDark;
	};
void pstSetTotals(PST,void *,int,void *,int *);

/* PST_CALCCELL */
struct inCalcCell {
	int iOrder;
	FLOAT rcm[3];
	};
struct outCalcCell {
	struct pkdCalcCellStruct mom;
	};
void pstCalcCell(PST,void *,int,void *,int *);

/* PST_COLCELLS */
struct inColCells {
	int iCell;
	int nCell;
	};
void pstColCells(PST,void *,int,void *,int *);

/* PST_DISTRIBCELLS */
void pstDistribCells(PST,void *,int,void *,int *);

/* PST_CALCROOT */
struct ioCalcRoot {
	struct ilCellNewt ilcn;
	};
void pstCalcRoot(PST,void *,int,void *,int *);

/* PST_DISTRIBROOT */
void pstDistribRoot(PST,void *,int,void *,int *);

/* PST_ONENODEREADINIT */
void pstOneNodeReadInit(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_SWAPALL */
void pstSwapAll(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_MASSCHECK */
struct outMassCheck {
	double dMass;
	};
void pstMassCheck(PST,void *,int,void *,int *);

/* PST_GETMASSCHANGE */
struct outGetMassChange {
	double dMassChange;
	};
void pstGetMassChange(PST,void *,int,void *,int *);

struct outMassMetalsEnergyCheck {
	double dTotMass;
        double dTotMetals;
        double dTotOx;
        double dTotFe;
        double dTotEnergy;
	};
void pstMassMetalsEnergyCheck(PST,void *,int,void *,int *);

struct inActiveTypeOrder {
        unsigned int iTestMask;
        };

void pstActiveTypeOrder(PST,void *,int,void *,int *);

/* PST_ACTIVEORDER */
void pstActiveOrder(PST,void *,int,void *,int *);

/* PST_SETRUNG */
struct inSetRung {
    int iRung;
    };
void pstSetRung(PST,void *,int,void *,int *);

struct inDensCheck {
    int iRung;
    int bGreater;
    int iMeasure;
    };
struct outDensCheck {
    double dMaxDensError;
    double dAvgDensError;
    int nError;
    int nTotal;
    };

void pstDensCheck(PST,void *,int,void *,int *);

struct inBallMax {
    int iRung;
    int bGreater;
    double dhFac;
    };

void pstBallMax(PST,void *,int,void *,int *);

/* PST_ACTIVERUNG */
struct inActiveRung {
    int iRung;
    int bGreater;
    };

void pstActiveRung(PST,void *,int,void *,int *);

/* PST_CURRRUNG */
struct inCurrRung {
    int iRung;
    };
struct outCurrRung {
    int iCurrent;
    };
void pstCurrRung(PST,void *,int,void *,int *);

/* PST_GRAVSTEP */
struct inGravStep {
    double dEta;
    };
void pstGravStep(PST,void *,int,void *,int *);

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    int    bDoGravity;
    int    bEpsAcc;
    int    bSqrtPhi;
    double dhMinOverSoft;
    };
void pstAccelStep(PST,void *,int,void *,int *);

/* PST_DENSITYSTEP */
struct inDensityStep {
    double dEta;
    double dRhoFac;
    };
void pstDensityStep(PST,void *,int,void *,int *);

/* PST_RUNGSTATS */
struct inRungStats {
	int iRung;
	};
struct outRungStats {
	int nParticles;
	};
void pstRungStats(PST,void *,int,void *,int *);

/* PST_GETMAP */
struct inGetMap {
	int nStart;
	};
void pstGetMap(PST,void *,int,void *,int *);

/* PST_RANDSEEDGENERATOR */
void pstRandSeedGenerator(PST,void *,int,void *,int *);

/* PST_COOLVELOCITY */
struct inCoolVelocity {
	int nSuperCool;
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
	};
void pstCoolVelocity(PST,void *,int,void *,int *);

struct inActiveType {
	unsigned int iFilterMask;
	unsigned int iTestMask;
	unsigned int iSetMask;
	int iRung;
	int bGreater;
	};

void pstActiveExactType(PST,void *,int,void *,int *);
void pstActiveType(PST,void *,int,void *,int *);
void pstSetType(PST,void *,int,void *,int *);
void pstResetType(PST,void *,int,void *,int *);
void pstCountType(PST,void *,int,void *,int *);
void pstActiveMaskRung(PST,void *,int,void *,int *);
void pstActiveTypeRung(PST,void *,int,void *,int *);

#define PST_SETTYPEFROMFILEMAXLEN 160
struct inSetTypeFromFile {
  int iSetMask;
  int biGasOrder;
  char file[PST_SETTYPEFROMFILEMAXLEN];
};

struct outSetTypeFromFile {
  int niOrder;
  int nSet;
  int nSetiGasOrder;
};

void pstSetTypeFromFile(PST,void *,int,void *,int *);

struct inSetParticleTypes {
	int nSuperCool;
	};

void pstSetParticleTypes(PST,void *,int,void *,int *);

#define MAXSOUGHTPARTICLELIST 10

struct inSoughtParticleList {
	int iTypeSought;
    int nMax;
	};

struct inoutParticleList {
    int n;
    struct SoughtParticle p[MAXSOUGHTPARTICLELIST];
	};

void pstSoughtParticleList(PST,void *,int,void *,int *);
void pstCoolUsingParticleList(PST,void *,int,void *,int *);

/* PST_RESMOOTH */
struct inMarkSmooth {
	int nSmooth;
	int bPeriodic;
	int bSymmetric;
	int iSmoothType;
	int iMarkType;
	SMF smf;
	};
void pstMarkSmooth(PST,void *,int,void *,int *);

struct inReSmooth {
	int nSmooth;
	int bPeriodic;
	int bSymmetric;
	int iSmoothType;
        double dfBall2OverSoft2;
	SMF smf;
	};
struct outReSmooth {
    int iSmoothFlags;  /* Warning Flags for need to smooth again, etc... */
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstReSmooth(PST,void *,int,void *,int *);

/* PST_INITACCEL */
void pstInitAccel(PST,void *,int,void *,int *);

/* PST_DTTORUNG */
struct inDtToRung {
    int iRung;
    double dDelta;
    int iMaxRung;
    int bAll;
    };
struct outDtToRung {
    int iMaxRung;
    int nMaxRung;
    int iMaxRungIdeal;
    };
void pstDtToRung(PST,void *,int,void *,int *);

/* PST_INITDT */
struct inInitDt {
    double dDelta;
    };
void pstInitDt(PST,void *,int,void *,int *);

/* PST_ORDWEIGHT */
struct inOrdWeight {
	int iOrdSplit;
	int iSplitSide;
	int ittr;
	};
struct outOrdWeight {
	int nLow;
	int nHigh;
	};
void pstOrdWeight(PST,void *,int,void *,int *);

#ifdef FAST_AGGS
/* PST_AGGSCACHEINIT */
void pstAggsCacheInit(PST,void *,int,void *,int *);	//JVD
#endif /* FAST_AGGS */

/* PST_SETWRITESTART */
struct inSetWriteStart {
	int nWriteStart;
	};
void pstSetWriteStart(PST,void *,int,void *,int *);

/* PST_SETNCWRITESTART */
struct inSetNCWriteStart {
	int nDarkWriteStart;
	int nGasWriteStart;
	int nStarWriteStart;
	};
void pstSetNCWriteStart(PST,void *,int,void *,int *);

/* PST_COLNPARTS */
struct outColNParts {
    int nNew;
    int nDeltaGas;
    int nDeltaDark;
    int nDeltaStar;
    };
void pstColNParts(PST, void *, int, void *, int *);

/* PST_NEWORDER */
void pstNewOrder(PST, void *, int, void *, int *);

/* PST_GETNPARTS */
/* see pkd.h
 struct outGetNParts { 
	int n;
    int nGas;
    int nDark;
    int nStar;
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    };
*/
void pstGetNParts(PST, void *, int, void *, int *);

/* PST_SETNPARTS */
struct inSetNParts {
    int nGas;
    int nDark;
    int nStar;
	int nMaxOrder;
    int nMaxOrderGas;
    int nMaxOrderDark;
    };
void pstSetNParts(PST, void *, int, void *, int *);

#ifdef GASOLINE

struct inGetGasPressure {
	enum GasModel iGasModel; 
  /* Adiabatic */
	double gamma;
	double gammam1;
  /* Isothermal */

  /* Ion evolving */

#ifdef GLASS
  /* Glass */
	double dGlassPoverRhoL;
	double dGlassPoverRhoR;
	double dGlassxL;
	double dGlassxR;
	double dxBoundL;
	double dxBoundR;
#endif
	};



/* PST_GETGASPRESSURE */
void pstGetGasPressure(PST, void *,int,void *,int *);

struct inLowerSoundSpeed {
	double dhMinOverSoft;
	};

void pstLowerSoundSpeed(PST, void *,int,void *,int *);

/* See cooling_*.h for details on input structs */
void pstInitEnergy(PST, void *,int,void *,int *);
void pstInitCooling(PST,void *,int,void *,int *);
void pstCoolTableRead(PST,void *,int,void *,int *);

/* PST_KICKVPRED */
struct inKickRhopred {
	double dHubbFac;
	double dDelta;
	};
void pstKickRhopred(PST,void *,int,void *,int *);

/* PST_SPHCURRRUNG */
struct inSphCurrRung {
    int iRung;
    int bGreater;
    };
struct outSphCurrRung {
    int iCurrent;
    };
void pstSphCurrRung(PST,void *,int,void *,int *);

#endif

#ifdef GLASS
/* PST_RANDOMVELOCITIES */
struct inRandomVelocities {
    double dMaxVelocityL;
    double dMaxVelocityR;
    };
void pstRandomVelocities(PST,void *,int,void *,int *);
#endif

#ifdef SPECIAL_PARTICLES

#include "special.h"

/* PST_GETSPECIALPARTICLES */
struct inGetSpecial {
	int nSpecial;
	int iId[MAX_NUM_SPECIAL_PARTICLES];
	SPECIAL_MASTER_INFO mInfo;
	};
struct outGetSpecial {
	SPECIAL_PARTICLE_INFO sInfo[MAX_NUM_SPECIAL_PARTICLES];
	};
void pstGetSpecialParticles(PST,void *,int,void *,int *);

/* PST_DOSPECIALPARTICLES */
struct inDoSpecial {
	int nSpecial;
	SPECIAL_MASTER_INFO mInfo;
	SPECIAL_PARTICLE_DATA sData[MAX_NUM_SPECIAL_PARTICLES];
	SPECIAL_PARTICLE_INFO sInfo[MAX_NUM_SPECIAL_PARTICLES];
	};
struct outDoSpecial {
	FLOAT aFrame[3];
	};
void pstDoSpecialParticles(PST,void *,int,void *,int *);

#endif

#ifdef COLLISIONS

/* PST_NUMREJECTS */
struct outNumRejects {
	int nRej;
	};
void pstNumRejects(PST,void *,int,void *,int *);

/* PST_READSS */
struct inReadSS {
	int nFileStart;
	int nFileEnd;
	int nDark;
	int nGas;			/* always zero */
	int nStar;			/* always zero */
	int iOrder;
	float fExtraStore;
	FLOAT fPeriod[3];	/* for compatability */
    char achInFile[PST_FILENAME_SIZE];
#ifdef SPRINGS
  int bReadSpringsData;
#endif
#ifdef DEM
  int bReadDEMData;
#endif
#ifdef DEM_TRACK_ORIENT
  int bReadOrient;
#endif
	};
void pstReadSS(PST,void *,int,void *,int *);

/* PST_WRITESS */
struct inWriteSS {
  char achOutFile[PST_FILENAME_SIZE];
  int bReduced;
#ifdef SPRINGS
  int bNoSpringsOutput;
#endif
#ifdef DEM
  int nDEMOutputs;
#endif
	};
void pstWriteSS(PST,void *,int,void *,int *);

/* PST_ADDUNIFGRAV */
struct inAddUnifGrav {
	double dgx,dgy,dgz;
	};
void pstAddUnifGrav(PST,void *,int,void *,int *);

/* PST_NEXTCOLLISION */
struct outNextCollision {
	double dt;
	int iOrder1,iOrder2;
	};
void pstNextCollision(PST,void *,int,void *,int *);

/* PST_GETCOLLIDERINFO */
struct inGetColliderInfo {
	int iOrder;
	};
struct outGetColliderInfo {
	COLLIDER Collider;
	};
void pstGetColliderInfo(PST,void *,int,void *,int *);

/* PST_DOCOLLISION */
struct inDoCollision {
  double dt;
  COLLIDER Collider1,Collider2;
  int bPeriodic;
  COLLISION_PARAMS CP;
  double dTime;
  double dDelta;
  int iStartStep;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
#ifdef AGGS
  int iAggNewIdx;
#endif
#ifdef COLLMOD_ZML
  DUST_BINS_PARAMS DB;
#endif
  };
struct outDoCollision {
  COLLIDER Out[MAX_NUM_FRAG];
  double dT;
  int iOutcome,nOut;
#ifdef AGGS_IN_PATCH
  FLOAT dx,dy,dvy;
#endif
#ifdef COLLMOD_ZML
  int iDustBin;
  DustBins DustBin;
#endif
  };
void pstDoCollision(PST,void *,int,void *,int *);

/* PST_RESETCOLLIDERS */
struct inResetColliders {
	int iOrder1,iOrder2;
	};
void pstResetColliders(PST,void *,int,void *,int *);

#ifdef JOHNNY

struct inFindEscapers {
  double dCentMass,dStarMass;
  };
struct outFindEscapers {
  int nEsc;
  };
void pstFindEscapers(PST,void *,int,void *,int *);

struct inCheckForKepler {
  double dCentMass,dDelta;
  };
struct outCheckForKepler {
  int nKep;
  };
void pstCheckForKepler(PST,void *,int,void *,int *);

#endif /* JOHNNY */

#endif /* COLLISIONS */

#ifdef RORY_EXTENSION

/* PST_SETBALL */
struct inSetBall {
  double dDelta;
  double dBallVelFact;
};
void pstSetBall(PST,void *,int,void *,int *);

/* PST_FINDTIGHTESTBINARY */
struct inBinary {
        int n;
        };
struct outBinary {
        int n;
        double dBindEn;
        int iOrder1;
        int iOrder2;
        };
void pstFindTightestBinary(PST,void *,int,void *,int *);

/* PST_MERGEBINARY */
struct inMrgBnry {
       COLLIDER c1;
       COLLIDER c2;
       int bPeriodic;
       double dDelta;
       double dTime;
       double dDensity;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
       int iStartStep;
       int iPad; /* Cheat to try to make the code owrk in Linux-ia64 */    
       };
struct outMrgBnry {
       COLLIDER cOut;
       int n;
       };
void pstMergeBinary(PST,void *,int,void *,int *);

/* PST_MOVEPART */
struct inMoveParticle {
    PARTICLE p;
    FLOAT dOrigin[3]; /* Origin of new frame */
    FLOAT dRelx[3]; /* Displacement from newx */
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
    };
void pstMoveParticle(PST,void *,int,void *,int *);

#ifdef SLIDING_PATCH

struct inLargeMass {
    double dMass;
    FLOAT fNumHillSphere;
    double dCentMass;
    double dOrbRad;    
    };
struct outLargeMass
{
       PARTICLE p[MAXLARGEMASS];           /* large mass particles */
       double dRadius[MAXLARGEMASS];        /* radius of sphere to be exchanged */    
       int n;                 /* # of large masses */
    };
void pstFindLargeMasses(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inGetNeighbors
{
    FLOAT x[3];
    double dDist;
  PATCH_PARAMS PP;
    double dTime;
    int id;
    };
struct outGetNeighbors 
{
    PARTICLE p[MAXNEIGHBORS];
    double dSep2[MAXNEIGHBORS];    
    int n;
    };
void pstGetNeighborParticles(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* SLIDING_PATCH */

#endif /* RORY_EXTENSION */

#ifdef AGGS

/* PST_AGGSFIND */
struct outAggsFind {
	int iMaxIdx;
	int iMaxOrg;
	};
void pstAggsFind(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSCOUNTPART */
struct inAggsCountPart {
	int iAggIdx;
	};
struct outAggsCountPart {
	int nPart;
	};
void pstAggsCountPart(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSDELETE */
struct inAggsDelete {
  int iAggIdx;
#ifdef AGGS_IN_PATCH
  double dStepTime;
  double dEventTime;
  double dOmega;
  PATCH_PARAMS PP;
#endif
	};
struct outAggsDelete {
	int bFound;
	};
void pstAggsDelete(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSMERGE */
struct inAggsMerge {
	int iOldIdx,iNewIdx;
	};
void pstAggsMerge(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSBACKDRIFT */
struct inAggsBackDrift {
	int iAggIdx;
	double dt;
	};
void pstAggsBackDrift(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSGETREF */
struct inAggsGetRef {
  int iAggIdx;
  };
struct outAggsGetRef {
  Vector r_ref;
  };
void pstAggsGetRef(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSGETPERIODICCOM */
struct inAggsGetPeriodicCOM {
  int iAggIdx;
  Vector r_ref;
  };
struct outAggsGetPeriodicCOM {
  Scalar m;
  Vector mr,mv;
  };
void pstAggsGetPeriodicCOM(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSSETUNWRAPPED */
struct inAggsSetUnwrapped {
  int iAggIdx;
  Vector r_com;
  };
void pstAggsSetUnwrapped(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSGETCOM */
struct inAggsGetCOM {
  int iAggIdx;
#ifdef AGGS_IN_PATCH
  PATCH_PARAMS PP;
#endif
  };
struct outAggsGetCOM {
  Scalar m;
  Vector mr,mv; /* position & velocity moments */
  };
void pstAggsGetCOM(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSGETAXESANDSPIN */
struct inAggsGetAxesAndSpin {
  int iAggIdx;
  Vector r_com,v_com;
#ifdef AGGS_IN_PATCH
  PATCH_PARAMS PP;
#endif
  };
struct outAggsGetAxesAndSpin {
  Matrix I;
  Vector L;
  Scalar rad_max_sq,volume;
  };
void pstAggsGetAxesAndSpin(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);
 
/* PST_AGGSSETBODYPOS */
struct inAggsSetBodyPos {
	int iAggIdx;
	Matrix spaceToBody; /* transpose of lambda */
	};
void pstAggsSetBodyPos(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSSETSPACEPOS */
struct inAggsSetSpacePos {
  int iAggIdx;
  Vector r_com;
  Matrix lambda;
#ifdef AGGS_IN_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif
	};
void pstAggsSetSpacePos(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSSETSPACEVEL */
struct inAggsSetSpaceVel {
  int iAggIdx;
  double dDelta; // JCM 12.16.20
  Vector v_com,a_com,omega,omegadot;
  Matrix lambda;
#ifdef AGGS_IN_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif
  };
void pstAggsSetSpaceVel(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSSETSPACESPINS */
struct inAggsSetSpaceSpins {
  int iAggIdx;
  Vector omega;
  };
void pstAggsSetSpaceSpins(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#ifndef DEM
/* PST_AGGSGETACCEL */
struct inAggsGetAccel {
	int iAggIdx;
	};
struct outAggsGetAccel {
	Scalar m;
	Vector ma; /* acceleration moment */
	};
void pstAggsGetAccel(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSCHECKSTRESS */
struct inAggsCheckStress {
  int iAggIdx;
  Scalar mass,rad_max,rad_eff;
  Vector r_com,a_com,omega;
  STRENGTH_PARAMS SP;
#ifdef AGGS_IN_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif /* AGGS_IN_PATCH */
  };
struct outAggsCheckStress {
  int nLost,nLeft;
  };
void pstAggsCheckStress(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSGETTORQUE */
struct inAggsGetTorque {
	int iAggIdx;
	Vector r_com,a_com;
	};
struct outAggsGetTorque {
	Vector torque;
	};
void pstAggsGetTorque(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#else /* DEM */

/* PST_AGGSGETACCELANDTORQUEDEM */
struct inAggsGetAccelAndTorqueDEM {
  int iAggIdx;
  Vector r_com;
  };
struct outAggsGetAccelAndTorqueDEM {
  Scalar m;
  Vector ma,torque;
  };
void pstAggsGetAccelAndTorqueDEM(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSSETMASSDEM */
struct inAggsSetMassDEM {
  int iAggIdx;
  FLOAT fMass;
  };
void pstAggsSetMassDEM(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* DEM */

#endif /* AGGS */

#ifdef WALLS

/* PST_WALLSUPDATESTUCKPARTICLES */
struct inWallsUpdateStuckParticles {
  WALL_PARAMS WP;
  int bUpdatePos;
  double dDelta;
  };
void pstWallsUpdateStuckParticles(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* WALLS */

#ifdef SPRINGS

/* PST_BREAKSPRINGS */
struct inBreakSprings {
        int iOrder1,iOrder2,bSpringsChanged;
        };
void pstBreakSprings(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* SPRINGS */

#ifdef DEM

struct inDEMZeroSmallMotions {
  double dAccCritSq,dDeltaSq;
  };
void pstDEMZeroSmallMotions(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct outDEMStats {
  DEM_STATS DEMStats;
};
void pstDEMStats(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#ifdef DEM_FIXED_BALL

struct inDEMSetBall {
  double dFixedRadius;
  int iFixedBallOption;
};
void pstDEMSetBall(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inDEMGetBallRadius {
  int iFixedBallOption;
  };
struct outDEMGetBallRadius {
  double dFixedRadius;
};
void pstDEMGetBallRadius(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

void pstDEMFixedBallNextLargest(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* DEM_FIXED_BALL */

#ifdef DEM_WALLS_REACT

struct outDEMWallsReact {
  double dTotalXForceFromParticles[MAX_NUM_WALL_ASSEMBLIES];
  double dTotalYForceFromParticles[MAX_NUM_WALL_ASSEMBLIES];
  double dTotalZForceFromParticles[MAX_NUM_WALL_ASSEMBLIES];
  double dTotalXTorqueFromParticles[MAX_NUM_WALL_ASSEMBLIES];
  double dTotalYTorqueFromParticles[MAX_NUM_WALL_ASSEMBLIES];
  double dTotalZTorqueFromParticles[MAX_NUM_WALL_ASSEMBLIES];
};
void pstDEMWallsReact(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* DEM_WALLS_REACT */

#endif /* DEM */

#ifdef CHARGE

struct ioChargeZMoments {
  CHARGE_PARAMS CP;
  };
void pstChargeZGetMoments(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);
void pstChargeZApplyMoments(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* CHARGE */

#ifdef DEM_TIDAL_SPACE

struct inDEMTidalFindPlanet {
  int iColor;
  };
struct outDEMTidalFindPlanet {
  int bFound;
  double dPos[3],dVel[3];
  };
void pstDEMTidalFindPlanet(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inDEMTidalFindMarker {
  int iOrder;
  };
struct outDEMTidalFindMarker {
  int bFound;
  double dPos[3],dVel[3];
  };
void pstDEMTidalFindMarker(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* DEM_TIDAL_SPACE */

#ifdef DEM_TIDAL_LOCAL

struct inDEMTidal {
  DEM_TIDAL d0,d;
  };
void pstDEMTidal(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* DEM_TIDAL_LOCAL */

#ifdef SLIDING_PATCH

/* PST_KICKPATCH */
struct inKickPatch {
    int bOpen;
    double dOrbFreq;
    double dvFacOne;
    double dvFacTwo;
#ifdef NEED_VPRED
    double dvPredFacOne;
    double dvPredFacTwo;
#endif /* NEED_VPRED */
    };
void pstKickPatch(PST,void *,int,void *,int *);

/* PST_RANDAZWRAP */
struct inRandAzWrap {
  PATCH_PARAMS PP;
};
struct outRandAzWrap {
	int nRandomized;
};
void pstRandAzWrap(PST,void *,int,void *,int *);

#endif /* SLIDING_PATCH */

#ifdef AGGS_IN_PATCH

/* PST_AGGSINPATCHGETREF */
struct inAggsInPatchGetRef {
  int iAggIdx;
  };
struct outAggsInPatchGetRef {
  double m_max;
  Vector r_max;
  };
void pstAggsInPatchGetRef(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSINPATCHGETUNWRAPPED */
struct inAggsInPatchGetUnwrapped {
  int iAggIdx;
  Vector r_ref;
  double dTime;
  PATCH_PARAMS PP;
  };
void pstAggsInPatchGetUnwrapped(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

/* PST_AGGSINPATCHOFFSET */
struct inAggsInPatchOffset {
  int iAggIdx;
  FLOAT dx,dy,dvy;
  int bDoUnwrapped;
  };
void pstAggsInPatchOffset(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* AGGS_IN_PATCH */

#ifdef SPINUP

struct inSpinupGetAgg {
  double dRInner, dROuter;
  };
struct outSpinupGetAgg {
  int bSuspend;
  };
void pstSpinupGetAgg(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inSpinupAlign {
  Vector r_com, vAxis;
  double dAngle;
  };
void pstSpinupAlign(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inSpinup {
  Vector r_com, spin;
  };
void pstSpinup(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#ifdef SPINUP_WITH_AGGS

struct inSpinupAccel {
    Vector r_com, omegadot;
    int bSpinupAlignAxes;
    };
void pstSpinupAccel(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

void pstSpinupResetAggs(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* SPINUP_WITH_AGGS */
#endif /* SPINUP */

#if defined(RUBBLE_ZML) || defined(COLLMOD_ZML)

struct outDustBinsInclAvg {
  double dDustBinsInclTot;
};
void pstDustBinsInclAvg(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct outDustBinsInclMax {
  double dDustBinsInclMax;
};
void pstDustBinsInclMax(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct inDustBinsGetMass {
  DUST_BINS_PARAMS DB;
  DustBins aDustBins[DUST_BINS_MAX_NUM];
  double dTimeInt;
  double dCentMass;
  };
struct outDustBinsGetMass {
  double aDustBinsMassLoss[DUST_BINS_MAX_NUM];
  };
void pstDustBinsGetMass(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inDustBinsApply {
  double dTimeInt;
  double dCentMass;
  DustBins aDustBins[DUST_BINS_MAX_NUM];
  double aMassIncrFrac[DUST_BINS_MAX_NUM];
  int nBins;
  int dVdisp;
#ifdef ORIGIN_HISTOGRAM
  FLOAT aOriginBins[DUST_BINS_MAX_NUM][NUM_ORIGIN_BINS];
#endif /* ORIGIN_HISTOGRAM */
  };
void pstDustBinsApply(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* RUBBLE_ZML || COLLMOD_ZML */

#ifdef RUBBLE_ZML

void pstRubbleResetColFlag(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct outRubbleCheckForKDKRestart {
	int bRestart;
	};
void pstRubbleCheckForKDKRestart(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inRubbleStep {
	double dMaxStep,dMinStep;
	};
void pstRubbleStep(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inRubCleanup {
	int iColor;
	DUST_BINS_PARAMS DB;
	double dCentMass;
	};
struct outRubCleanup {
	DustBins aDustBins[DUST_BINS_MAX_NUM];
	DustBins DustBinsTrash;
	};
void pstRubCleanup(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inRubInterpCleanup {
	DUST_BINS_PARAMS DB;
	double dCentMass;
	int iOrder;
	};
struct outRubInterpCleanup {
	DustBins DustBinsInterp;
	int iBin;
	};
void pstRubInterpCleanup(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#elif defined(COLLMOD_ZML)

void pstCollModResetColFlag(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct outCollModCheckForKDKRestart {
  int bRestart;
	};
void pstCollModCheckForKDKRestart(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

struct inCollModStep {
  double dMaxStep,dMinStep;
  };
void pstCollModStep(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* RUBBLE_ZML, COLLMOD_ZML */

#ifdef ORIGIN_HISTOGRAM

struct inInitializeOriginBins {
	DUST_BINS_PARAMS DB;
};
void pstInitializeOriginBins(PST pst,void *vIn,int nIn,void *vOut,int *pnOut);

#endif /* ORIGIN_HISTOGRAM */ 

#ifdef GASOLINE

/* PST_SPHSTEP */
struct inSphStep {
    double dCosmoFac;
    double dEtaCourant;
    double dEtauDot;
    int bViscosityLimitdt;
    };
void pstSphStep(PST,void *,int,void *,int *);

struct inSphViscosityLimiter {
    int bOn;
    int bShockTracker;
    };
void pstSphViscosityLimiter(PST,void *,int,void *,int *);

struct inUpdateShockTracker {
    double dDelta;
    double dShockTrackerA;
    double dShockTrackerB;
};

void pstUpdateShockTracker(PST,void *,int,void *,int *);

#ifdef STARFORM
/* PST_FORMSTARS */
struct inFormStars
{
    struct stfmContext stfm;
    double dTime;
    };

struct outFormStars 
{
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

void pstFormStars(PST,void *,int,void *,int *);

struct inFeedback
{
    struct fbContext fb;
    struct snContext sn;
    double dTime;
    double dDelta;
    };

struct outFeedback
{
    FBEffects fbTotals[FB_NFEEDBACKS];
    };

void pstFeedback(PST,void *,int,void *,int *);

#endif

#ifdef SIMPLESF
struct inSimpleStarForm
{
    double dRateCoeff;
    double dTMax;
    double dDenMin;
    double dDelta;

	double dTime;
	double dInitStarMass;
	double dESNPerStarMass;
	double dtCoolingShutoff;
    int bdivv;
    };

struct outSimpleStarForm 
{
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

void pstSimpleStarForm(PST,void *,int,void *,int *);

#endif

#endif

/* Return is pixmap */
void pstDumpFrame(PST,void *,int,void *,int *);

void pstCOM(PST,void *,int,void *,int *);

void pstCOMByType(PST,void *,int,void *,int *);

void pstOldestStar(PST,void *,int,void *,int *);

struct inSetSink 
{
  double dSinkMassMin;
  };

struct outSetSink 
{
  int nSink;
  };

void pstSetSink(PST,void *,int,void *,int *);

/* Return is pixmap */
void pstDumpVoxel(PST,void *,int,void *,int *);

/* PST_GROWMASS */
struct inGrowMass 
{
    int nGrowMass;
    double dDeltaM;
    };

void pstGrowMass(PST,void *,int,void *,int *);

/* PST_CLEARTIMER */
struct inClearTimer 
{
    int iTimer;
    };

void pstClearTimer(PST,void *,int,void *,int *);

/* PST_MASSINR */
struct inMassInR 
{
    double R;
    };

struct outMassInR
{
    double dMass;
    FLOAT com[3];
    };
void pstMassInR(PST,void *,int,void *,int *);

struct inRotBar
{
    struct rotbarContext rotbar;
    };
void pstInitRotBar(PST,void *,int,void *,int *);

/* PST_KICKVPRED */
#ifdef NEED_VPRED
struct inKickVpred {
	double dvFacOne;
	double dvFacTwo;
	double duDelta;
	double duDotLimit;
	int iGasModel;
	double z;
	};
void pstKickVpred(PST,void *,int,void *,int *);
#endif
#endif
