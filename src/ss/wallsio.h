#ifndef WALLSIO_HINCLUDED
#define WALLSIO_HINCLUDED

#define WALLS_MAX_STR_LEN 4096 /* must be at least 1 */

#ifndef MAX_NUM_WALL_ASSEMBLIES
#define MAX_NUM_WALL_ASSEMBLIES 2
#endif

/* wall types */

enum {
  WallPlane=0,
  WallTriangle,
  WallRectangle,
  WallDisk,
  WallCylinderInfinite,
  WallCylinderFinite,
  WallShell
};

/* verbosity flags */

#define WALLS_SILENT  0
#define WALLS_VERBOSE 1

/* wall data struct */

typedef struct {
  int iType;
  double vOrigin[3]; /* can treat these as 3-vectors */
  double vOrient[3];
  double vVertex1[3];
  double vVertex2[3];
  double vVel[3];
  double dStep;
  double dOscAmp;
  double dOscFreq;
  double dOscPhase;
  double vOscVec[3];
  double dRadius;
  double dHoleRadius;
  double dLength;
  double dTaper;
  double dOpenAngle;
  double dAngSpeed;
  double dEpsN;
  double dEpsT;
  int iColor;
  double dTrans;
  double dKn;
  double dKt;
  double dMuS;
  double dMuR;
  double dMuT;
  double dBetaS;
  int iCohesionModel;
  double dCohesiveCoeff;
  double dKnOuter;
  double dKtOuter;
  double dKnInner;
  double dKtInner;
  double dInnerOverlapBndry;
  double dMass;
  double vCoG[3];
  double vRot[3];
  double vIxj[3];
  double vIyj[3];
  double vIzj[3];
  double vPxj[3];
  double vPyj[3];
  double vPzj[3];
  int iAssembly;
  double accel[3];
  double torque[3];
  double moments[3];
  } WALL_DATA;

/* following values indicate particle property takes precedence */

#define PARTICLE_OVERRIDE_EPSN (99.0)
#define PARTICLE_OVERRIDE_EPST (99.0)
#define PARTICLE_OVERRIDE_KN (-1.0)
#define PARTICLE_OVERRIDE_KT (-1.0)
#define PARTICLE_OVERRIDE_MUS (-1.0)
#define PARTICLE_OVERRIDE_MUR (-1.0)
#define PARTICLE_OVERRIDE_MUT (-1.0)
#define PARTICLE_OVERRIDE_BETAS (-1.0)
#define PARTICLE_OVERRIDE_COMO (-1) /* cohesion model */
#define PARTICLE_OVERRIDE_COCO (-1.0) /* cohesive coefficient */
#define PARTICLE_OVERRIDE_IOB (-1.0) /* inner overlap boundary */
#define PARTICLE_OVERRIDE_KNOUTER (-1.0)
#define PARTICLE_OVERRIDE_KTOUTER (-1.0)
#define PARTICLE_OVERRIDE_KNINNER (-1.0)
#define PARTICLE_OVERRIDE_KTINNER (-1.0)

/* function prototype */

int wallsParseWallsFile(FILE *fp,int *nWalls,WALL_DATA **pWallsData,double *dTime,int bVerbose);

#endif /* !WALLSIO_HINCLUDED */
