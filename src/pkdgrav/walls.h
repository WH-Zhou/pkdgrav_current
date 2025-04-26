#ifndef WALLS_HINCLUDED
#define WALLS_HINCLUDED

#ifdef WALLS

#include "linalg.h"
#include "wallsio.h"

/* #define MAX_NUM_WALLS 5000 *//* drop this value when not using this many walls */
#define MAX_NUM_WALLS 20 /* upper limit because may need to pass wall data to remote processors */

typedef struct {
  int iWallID;
  WALL_DATA wd; 
  Vector vOscVel; /* needed for oscillating walls (hard sphere only) */
  Vector vTravel,vTotVel; /* needed for updating stuck particles */
  } WALL;

typedef struct {
  int nWalls;
  WALL pWalls[MAX_NUM_WALLS];
  int bWallsEdgeDetect,bWallsSolveQuartic;
  } WALL_PARAMS;

/* following for hard sphere only */

extern const int nWallFaces[]; /* defined in walls.c */

#define WALL_MAX_NUM_FACES 10 /* MUST be max of entries in nWallFaces[] */

enum {PlaneUpper=0,PlaneLower}; /* 2 faces */
enum {TriangleUpper=0,TriangleLower,TriangleSide1,TriangleSide2,TriangleSide3,TriangleVertex1,TriangleVertex2,TriangleVertex3}; /* 8 "faces" */
enum {RectangleUpper=0,RectangleLower,RectangleSide1,RectangleSide2,RectangleSide3,RectangleSide4,RectangleCorner1,RectangleCorner2,RectangleCorner3,RectangleCorner4}; /* 10 "faces" */
enum {DiskUpper=PlaneUpper,DiskLower=PlaneLower,DiskRing,DiskHole}; /* 4 "faces" */
enum {CylinderInfiniteInner=0,CylinderInfiniteOuter}; /* 2 faces */
enum {CylinderFiniteInner=CylinderInfiniteInner,CylinderFiniteOuter=CylinderInfiniteOuter,CylinderFiniteRingUpper,CylinderFiniteRingLower}; /* 4 "faces" */
enum {ShellInner=0,ShellOuter=0,ShellOpening}; /* 3 "faces" */

/* handy macros */

#include "pkd.h" /* for PARTICLE */

int PARTICLE_STUCK(const PARTICLE *p);
#define PARTICLE_STUCK(p) ((p)->iColor < 0)

int PARTICLE_WALL_ID(const PARTICLE *p);
#define PARTICLE_WALL_ID(p) (-1 - (p)->iColor)

int PARTICLE_STUCK_ON_ROTATING_WALL(const PARTICLE *p,const WALL_PARAMS *WP);
#define PARTICLE_STUCK_ON_ROTATING_WALL(p,WP) (PARTICLE_STUCK(p) && (WP)->pWalls[PARTICLE_WALL_ID(p)].wd.dAngSpeed != 0.0)

#include "collision.h" /* for COLLIDER */

int COLLIDER_STUCK(const COLLIDER *c);
#define COLLIDER_STUCK(c) ((c)->iColor < 0)

int COLLIDER_WALL_ID(const COLLIDER *c);
#define COLLIDER_WALL_ID(c) (-1 - (c)->iColor)

int COLLIDER_STUCK_ON_ROTATING_WALL(const COLLIDER *c,const WALL_PARAMS *WP);
#define COLLIDER_STUCK_ON_ROTATING_WALL(c,WP) (COLLIDER_STUCK(c) && (WP)->pWalls[COLLIDER_WALL_ID(c)].wd.dAngSpeed != 0.0)

/* function prototypes */

void pkdWallsUpdateStuckParticles(PKD pkd,const WALL_PARAMS *WP,int bUpdatePos,double dDelta);
void wallsRotate(const WALL_PARAMS *WP,int iWallID,Vector vPos,Vector vVel,Vector vOmegaV,double dt);
void wallsGetTimesToIntersect(const WALL_PARAMS *WP,const WALL *w,PARTICLE *p,double dStart,double dTmax,double dtCol[WALL_MAX_NUM_FACES],int *bOverlap);

#ifdef DEM_WALLS_REACT
int wallsEulerDerivs(FLOAT t,const FLOAT vars[],FLOAT derivs[],void *wd_as_void);
#endif /* DEM_WALLS_REACT */

#endif /* WALLS */

#endif /* !WALLS_HINCLUDED */
