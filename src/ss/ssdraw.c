/*
** ssdraw.c -- DCR 01-08-08
** ========
** Generates movie frames from ss data.
**
** 3/16/07: added tree (stolen from rpa) to locate potential minimum.
** 6/24/17: added support for sphere orientations
*/

#include <ss.h>
#include <string.h> /* for strlen() */
#include <unistd.h> /* for getpid() */
#include <libgen.h> /* for basename() */
#include <math.h>
#include <assert.h>
#include <boolean.h>
#include <vector.h>
#include <mem_util.h>
#include <raster.h>
#include <rdpar.h>
#include <wallsio.h>

#ifndef INT_MAX
# define INT_MAX 2147483647
#endif

#define PARFILE_DFLT "ssdraw.par"
#define DRAW_EXT     ".ras"
#define POV_EXT      ".pov"
#define ORIENT_EXT   ".ori"

#define MAX_NUM_WARNINGS 1000

/* light/camera/view options */

enum {LightFIX=0,LightLRG,LightPOT,LightCOL,LightPID,LightOID,LightCAM};
enum {CamFIX=0,CamLRG,CamPOT,CamCOL,CamPID,CamOID};
enum {ViewFIX=0,ViewLRG,ViewPOT,ViewCOL,ViewPID,ViewOID,ViewDIR};

/* shape options */

typedef enum {Point,SolidSphere,POV} SHAPE;

/* gravitational potential sampling order */

enum {PotRan=0,PotOrd,PotLrg};

/* handy limits */

#define MIN_NUM_PIX 2 /* for drawing a solid shape */
#define MAX_NUM_PIX (0.1*SQ(p->iFrameSize)) /* arbitrarily < max image pixels */

#define MIN_POV_SIZE 3.0e-4 /* minimum size (any dimension) in units of camera distance for POV-Ray drawing of particles and/or degenerate objects */

/* for the patch */

#define SHEAR(ix,lx,ly,w,t)\
  ((ix) < 0 ? fmod(0.5*(ly) - 1.5*(ix)*(w)*(lx)*(t),(ly)) - 0.5*(ly):\
   (ix) > 0 ? 0.5*(ly) - fmod(0.5*(ly) + 1.5*(ix)*(w)*(lx)*(t),(ly)): 0.0)

/* handy structs */

typedef struct {
  int nFramesActive;
  int iID;
  BOOLEAN bIsCOM;
  } Inertia;

typedef struct {
  int iOpt;
  int iLock;
  VECTOR vOrig;
  VECTOR vPos;
  VECTOR vPosNorm; /* needed for POV-Ray renormalizations */
  Inertia Target;
  } Control;

typedef struct {
  BOOLEAN b24bit; /* true if input colors are to be interpreted as 24-bit RGB tuples (POV only) */
  BOOLEAN bDoOrient; /* true to draw sphere orientations (POV only) */
  BOOLEAN bWallsOnly,bNoGZip;
  /* light, camera, view control (combines supplied & derived params) */
  Control Light,Camera,LookAt;
  /* inertial control for view size */
  Inertia ViewSizeTarget;
  /* other supplied parameters */
  int iFrameSize;
  double dLightInten;
  double dFocusDist;
  VECTOR vSkyVec;
  double dViewSize;
  double dViewScaleCrit;
  double dViewScale;
  BOOLEAN bScaleLengthByRadius;
  double dViewTargetZoomFactor;
  SHAPE Shape;
  double dRScale;
  int iNewColor;
  BOOLEAN bZSort;
  int iFirstColor;
  /* view motion controls */
  BOOLEAN bStartAtCOM;
  double dInertia;
  double dMinTargetMassFraction;
  /* POV-Ray stuff */
  char achShapeFile[MAXPATHLEN];
  BOOLEAN bRenorm;
  double dViewSizeNorm;
  char achAspectRatio[MAXPATHLEN];
  double dCamInten;
  double dBlob;
  int iHighlightIndex;
  /* patch stuff */
  BOOLEAN bUsePatch;
  double dLx;
  double dLy;
  double dOmega;
  int nRepl;
  int iGhostColor;
  BOOLEAN bDrawPatchLines;
  /* wall stuff */
  char achWallFile[MAXPATHLEN];
  int bDoWalls;
  double dWallTimeOffset;
  int nWalls;
  WALL_DATA *pWalls;
  double dTime; /* read from input data files */
  /* gravitational potential stuff */
  int nSample;
  int iSampleOrder;
  double dTheta;
  /* units */
  double dLenUnit;
  double dTimeUnit;
  /* derived parameters */
  VECTOR vPotPos;
  double dOldViewSize;
  double dViewRatio;
  double dPixelSize;
  VECTOR vXHat;
  VECTOR vYHat;
  VECTOR vZHat;
  } PARAMS;

typedef struct {
  double vdotz;
  int zindex;
  } ZSort;

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  VECTOR vCenter,vComPos;
  double dSize,dHalfSize,dMass;
  struct node_s *pCell[CELLS_PER_NODE];
  const SSDATA *pLeaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

struct mass_sort_s {
  int iIndex;
  double dMass;
  };

enum {POV_DRAW_REGULAR,POV_DRAW_HIGHLIGHT,POV_DRAW_TRANSPARENT};

/* functions follow */

static const char *pov_color(const PARAMS *p,int iColor,BOOLEAN bConv8to24)
{
	/*
	** NOTE: The bConv8to24 flag is so that walls can still use 8-bit
	** color specifications while particles use 24 bit.
	*/

	static char sColor[64]; /* enough space for red, green, blue values */

	if (p->b24bit && !bConv8to24) {
		float fRed,fGreen,fBlue;
		assert(iColor >= 0 && iColor <= 0xffffff);
		fRed = (float) (iColor & 0xff)/255;
		fGreen = (float) (iColor & 0xff00)/(256*255);
		fBlue = (float) (iColor & 0xff0000)/(256*256*255);
		snprintf(sColor,64,"red %g green %g blue %g",fRed,fGreen,fBlue);
		return sColor;
		}

	/* 8 bit... */

	assert(iColor >= 0);

	if (iColor >= NUM_COLORS) {
		fprintf(stderr,"pov_color(): invalid color (%i). Try using '-c' option.\n",iColor);
		exit(1);
		}

	switch (iColor) { /* color names require #include "colors.inc" */
	case BLACK:
		return "Black";
	case WHITE:
		return "White";
	case RED:
		return "Red";
	case GREEN:
		return "Green";
	case BLUE:
		return "Blue";
	case YELLOW:
		return "Yellow";
	case MAGENTA:
		return "Magenta";
	case CYAN:
		return "Cyan";
	case GOLD:
		return "Gold";
	case PINK:
		return "Pink";
	case ORANGE:
		return "Orange";
	case KHAKI:
		return "Khaki";
	case VIOLET:
		return "Violet";
	case MAROON:
		return "Maroon";
	case AQUA:
		return "Aquamarine";
	case NAVY:
		return "Navy";
	default:
	  {
		  float fColor = (float) (iColor - 16)/(255 - 16);
		  snprintf(sColor,64,"red %g green %g blue %g",fColor,fColor,fColor);
		  return sColor;
		  }
		}
	}

static void pov_draw_particle(const PARAMS *p,FILE *fp,const VECTOR vPos,double dRadius,int iColor,int iTransColor,int iStatus,MATRIX o)
{
	double dTransmit = 0.0;

	if (iColor < 0 && p->bDoWalls) { /* particle stuck to wall! (inherits wall color & transparency) */
		int iWallID = -1 - iColor;
		if (iWallID >= p->nWalls) {
			fprintf(stderr,"pov_draw_particle(): color %i incompatible with wall data.\n",iColor);
			exit(1);
			}
		iColor = p->pWalls[iWallID].iColor;
		dTransmit = p->pWalls[iWallID].dTrans;
		}

	/*
	** Blob strategy: exaggerate size of sphere (by factor of 4) to
	** ensure smooth overlaps. This means object sizes will appear
	** larger than they would without blobbing (by ~ a factor of 2).
	*/

	fprintf(fp,"sphere {<0,0,0>,%s texture {StyleParticle pigment {color ",
			p->dBlob ? "4 1" : "1");

	switch(iStatus) {
	case POV_DRAW_REGULAR:
		fprintf(fp,"%s transmit %g}%s}",pov_color(p,iColor,FALSE),
				iColor == iTransColor ? 0.95 : dTransmit,
				iColor == 0 ? " finish {phong 0.05}" : "");
		break;
	case POV_DRAW_HIGHLIGHT:
		fprintf(fp,"Red} finish {ambient 1.0}}");
		break;
	case POV_DRAW_TRANSPARENT:
		fprintf(fp,"White filter 0.9}} interior {ior 1.01}");
		break;
	default:
		assert(0); /* shouldn't be here */
		}
	if (p->bDoOrient)
		fprintf(fp, " texture {StyleOrient rotate <90,0,0> matrix <%g,%g,%g, %g,%g,%g, %g,%g,%g, 0,0,0>}",o[0][0],o[0][1],o[0][2],o[1][0],o[1][1],o[1][2],o[2][0],o[2][1],o[2][2]);
	fprintf(fp," scale %g translate <%g,%g,%g>}\n",
			dRadius*p->dRScale,vPos[X],vPos[Y],vPos[Z]);
}

static void pov_blob_open(const PARAMS *p,FILE *fp)
{
	fprintf(fp,"blob {\nthreshold %g\n",p->dBlob);
}

static void pov_blob_close(FILE *fp)
{
	fprintf(fp,"finish {phong 0.5}\n}\n");
}

/*DEBUG not implemented
static void pov_draw_square(const PARAMS *p,FILE *fp,double x1,double y1,
							double x2,double y2,int color)
{
}
*/

static void pov_warn_min_size(const PARAMS *p,const char achType[])
{
	static int nWarnings = 0;

	if (nWarnings < 10) {
		fprintf(stderr,"WARNING: %s -- scaled to minimum.\n",achType);
		if (!p->bRenorm)
			fprintf(stderr,"WARNING: without renormalization, object may not be visible.\n");
		++nWarnings;
		}
	else if (nWarnings == 10) {
		fprintf(stderr,"(Further scaling warnings suppressed.)\n");
		++nWarnings;
		}
	}

static void pov_draw_wall(const PARAMS *p,FILE *fp,const WALL_DATA *w)
{
	VECTOR vOrigin,vTravel;
	double dPhi,dTheta,dTime,dScale;

	if (w->iColor == BLACK)
		return;

	COPY_VEC(w->vOrigin,vOrigin); /* so we can modify it in the case of moving walls... */

	/* rotation angles (not needed in all cases) */

	dPhi = atan2(w->vOrient[Y],w->vOrient[X])*180/M_PI;
	dTheta = atan2(sqrt(w->vOrient[X]*w->vOrient[X] + w->vOrient[Y]*w->vOrient[Y]),w->vOrient[Z])*180.0/M_PI;

	/* account for any wall motion */

	dTime = (w->dStep == 0.0 ? p->dTime : ((int) (p->dTime/w->dStep))*w->dStep);
	
	COPY_VEC(w->vVel,vTravel);
	SCALE_VEC(vTravel,dTime);
	ADD_VEC(vOrigin,vTravel,vOrigin);

	/* and oscillation */

	COPY_VEC(w->vOscVec,vTravel);
	SCALE_VEC(vTravel,w->dOscAmp*sin(w->dOscFreq*dTime + w->dOscPhase));
	ADD_VEC(vOrigin,vTravel,vOrigin);

	/* useful scaling factor */

	dScale = (p->bRenorm ? p->dViewSizeNorm : p->dViewSize)*p->dViewRatio;

	/*
	** Strategy: POV-Ray patterns, used in many textures, are
	** generally defined relative to the origin (0,0,0) and on scales
	** of order unity.  So, we construct our basic objects at the
	** origin and scaled to order unity first, then we rescale,
	** texturize, rotate, and translate to the desired configuration.
	*/

	switch (w->iType) {
	case WallPlane:
		fprintf(fp,"box {<-0.5,-0.5,0>,<0.5,0.5,0> scale %g*<1,1,1> texture {StylePlane} rotate <0,%g,%g>",100.0*dScale,dTheta,dPhi);
		break;
	case WallTriangle: { /* treated a lot differently from the rest */
		VECTOR v1,v2;
		ADD_VEC(vOrigin,w->vVertex1,v1);
		ADD_VEC(vOrigin,w->vVertex2,v2);
		fprintf(fp,"triangle {<%g,%g,%g>,<%g,%g,%g>,<%g,%g,%g> texture {StyleTriangle pigment {color %s transmit %g}}}\n",vOrigin[X],vOrigin[Y],vOrigin[Z],v1[X],v1[Y],v1[Z],v2[X],v2[Y],v2[Z],pov_color(p,w->iColor,p->b24bit),w->dTrans);
		return;
		}
	case WallRectangle: {/* patterns in this case won't be centered */
		VECTOR n,v1,v2;
		COPY_VEC(w->vOrient,n);
		COPY_VEC(w->vVertex1,v1);
		COPY_VEC(w->vVertex2,v2);
		fprintf(fp,"box {<0,0,0>,<0,1,1> matrix <%g,%g,%g,%g,%g,%g,%g,%g,%g,0.,0.,0.> texture {StyleRectangle}",n[X],n[Y],n[Z],v1[X],v1[Y],v1[Z],v2[X],v2[Y],v2[Z]);
		/* description of translation matrix: http://bobobobo.wordpress.com/2009/03/15/rotating-a-vector-to-match-up-with-another-vector/ */
		break;
		}
	case WallDisk:
		if (w->dRadius == 0.0) { /* degenerate case: singularity! */
			fprintf(fp,"sphere {<0,0,0>,1 scale %g texture {Bright}",MIN_POV_SIZE*dScale);
			pov_warn_min_size(p,"degenerate disk (point)");
			}
		else
			fprintf(fp,"disc {<0,0,0>,<0,0,1>,1,%g scale %g texture {StyleDisk} rotate <0,%g,%g>",w->dHoleRadius/w->dRadius,w->dRadius,dTheta,dPhi);
		break;
	case WallCylinderInfinite: {
		double dRadius = w->dRadius;
		if (dRadius == 0.0)
			dRadius = MIN_POV_SIZE*dScale; /* degenerate case: line */
		fprintf(fp,"cylinder {<0,0,-0.5>,<0,0,0.5>,1 open scale <%g,%g,%g> texture {StyleCylinderInfinite} rotate <0,%g,%g>",dRadius,dRadius,10.0*dScale,dTheta,dPhi);
		}
		break;
	case WallCylinderFinite: {
		double dRadius,dLength;
		/* handle degenerate cases */
		if (w->dRadius == 0.0) {
			dRadius = MIN_POV_SIZE*dScale;
			if (w->dLength == 0.0) {
				fprintf(fp,"sphere {<0,0,0>,1 scale %g texture {Bright}",dRadius);
				pov_warn_min_size(p,"degenerate cylinder (point)");
				break;
				}
			pov_warn_min_size(p,"degenerate cylinder (line)");
			}
		else
			dRadius = w->dRadius;
		if (w->dLength == 0.0) {
			dLength = MIN_POV_SIZE*dScale;
			pov_warn_min_size(p,"degenerate cylinder (ring)");
			}
		else
			dLength = w->dLength;
		fprintf(fp,"cone {<0,0,-0.5>,1,<0,0,0.5>,%g open scale <%g,%g,%g> texture {StyleCylinderFinite} rotate <0,%g,%g>",1.0 - w->dTaper,dRadius,dRadius,dLength,dTheta,dPhi);
		}
		break;
	case WallShell:
		if (w->dRadius == 0.0) { /* handle degenerate case */
			fprintf(fp,"sphere {<0,0,0>,1 scale %g texture {Bright}",MIN_POV_SIZE*dScale);
			pov_warn_min_size(p,"degenerate shell (point)");
			}
		else
			fprintf(fp,"difference {difference {sphere {<0,0,0>,1} sphere {<0,0,0>,0.99 hollow}} box {<-1,-1,%g>,<1,1,1>} scale %g texture {StyleShell} rotate <0,%g,%g>",cos(w->dOpenAngle*M_PI/180.0),w->dRadius,dTheta,dPhi);
		break;
	default:
		assert(0);
		}
	/* the end portion is always the same... */
	fprintf(fp," translate <%g,%g,%g> texture {pigment {color %s transmit %g}}}\n",vOrigin[X],vOrigin[Y],vOrigin[Z],pov_color(p,w->iColor,p->b24bit),w->dTrans);
	}

static void pov_set_scene(PARAMS *p,FILE *fp)
{
	VECTOR vPos;

	/* outputs scene data suitable for use with POV-Ray 3.x */

	fprintf(fp,"#include \"%s\"\n",p->achShapeFile);

	/* light source info */

	if (p->bRenorm) {
		COPY_VEC(p->Light.vPosNorm,vPos);
		}
	else {
		COPY_VEC(p->Light.vPos,vPos);
		}

	fprintf(fp,"light_source {<%g,%g,%g> color %g}\n",
			vPos[X],vPos[Y],vPos[Z],p->dLightInten);

	if (p->bRenorm) {
		COPY_VEC(p->Camera.vPosNorm,vPos);
		}
	else {
		COPY_VEC(p->Camera.vPos,vPos);
		}

	/* extra light on camera (optional) */

	if (p->dCamInten > 0.0)
		fprintf(fp,"light_source {<%g,%g,%g> color %g}\n",
				vPos[X],vPos[Y],vPos[Z],p->dCamInten);

	/* camera info */

	/*DEBUG could use "angle" (FOV) in place of "direction" here...*/
	fprintf(fp,"camera {location <%g,%g,%g> sky <%g,%g,%g> direction <0,0,%g> right -x*%s ",
			vPos[X],vPos[Y],vPos[Z],p->vSkyVec[X],p->vSkyVec[Y],p->vSkyVec[Z],p->dViewRatio,p->achAspectRatio);

	if (p->bRenorm) {
		ZERO_VEC(vPos);
		}
	else {
		COPY_VEC(p->LookAt.vPos,vPos);
		}

	fprintf(fp,"look_at <%g,%g,%g>}\n",vPos[X],vPos[Y],vPos[Z]);
	}

static BOOLEAN project(const PARAMS *p,VECTOR pos,int *ix,int *iy)
{
	/*
	** Projects "pos" onto coordinates "*ix" and "*iy" of projection
	** screen.  Returns FALSE if object behind camera or displaced
	** more than INT_MAX from the screen origin, otherwise TRUE.
	*/

	VECTOR vProjVec;
	double rdotz,q;

	SUB_VEC(pos,p->Camera.vPos,vProjVec);
	rdotz = DOT(vProjVec,p->vZHat);
	if (rdotz <= 0)
		return FALSE;
	NORM_VEC(vProjVec,rdotz);
	SUB_VEC(vProjVec,p->vZHat,vProjVec);
	SCALE_VEC(vProjVec,p->dViewRatio);
	q = (0.5 + DOT(vProjVec,p->vXHat))*p->iFrameSize;
	if (fabs(q) > INT_MAX)
		return FALSE;
	*ix = (int) q;
	q = (0.5 - DOT(vProjVec,p->vYHat))*p->iFrameSize;
	if (fabs(q) > INT_MAX)
		return FALSE;
	*iy = (int) q;
	return TRUE;
	}

static void draw_point(const PARAMS *p,IMAGE_T *image,VECTOR pos,int color)
{
	int ix,iy;

	if (project(p,pos,&ix,&iy) == FALSE)
		return;
	DrawPoint(image,ix,iy,(COLOR) color);
}

static void draw_line(const PARAMS *p,IMAGE_T *image,VECTOR v1,VECTOR v2,int color)
{
	int x1,y1,x2,y2;

	if (project(p,v1,&x1,&y1) == FALSE)
		return;
	if (project(p,v2,&x2,&y2) == FALSE)
		return;
	DrawVector(image,x1,y1,x2,y2,(COLOR) color);
}

static void draw_square(const PARAMS *p,IMAGE_T *image,double x1,double y1,
						double x2,double y2,int color)
{
	VECTOR v1,v2;

	SET_VEC(v1,x1,y1,0);
	SET_VEC(v2,x2,y1,0);
	draw_line(p,image,v1,v2,(COLOR) color);
	SET_VEC(v1,x2,y2,0);
	draw_line(p,image,v2,v1,(COLOR) color);
	SET_VEC(v2,x1,y2,0);
	draw_line(p,image,v1,v2,(COLOR) color);
	SET_VEC(v1,x1,y1,0);
	draw_line(p,image,v2,v1,(COLOR) color);
}

static void draw_solid(const PARAMS *p,IMAGE_T *image,VECTOR pos,double r,int color)
{
	/*
	** Draws a shaded sphere (based on algorithm "spheres.c" 1.4 88/02/05
	** Copyright 1986 Sun Microsystems).
	*/

	MATRIX xform;
	VECTOR a,b,v,spos,vv;
	double dNumPix,dLimit,m,s,rr,dstheta,stheta,ctheta,dphi,phi;
	int ix,iy;

	/* Adjust pixel size to compensate for particle distance -- may be SLOW */

	SUB_VEC(p->Camera.vPos,pos,v);
	SUB_VEC(p->LookAt.vPos,p->Camera.vPos,vv);
	s = p->dPixelSize*sqrt(MAG_SQ(v)/MAG_SQ(vv));

	dNumPix = PI*SQ(r/s); /* approximate */

	if (dNumPix < MIN_NUM_PIX) {
		draw_point(p,image,pos,(COLOR) color);
		return;
	}

	if (dNumPix > MAX_NUM_PIX) {
		BOOLEAN bNotDrawn = (dNumPix/MAX_NUM_PIX > 10); /* should be param? */
		fprintf(stderr,"draw_solid(): Object %s",
				(bNotDrawn?"too close/big to draw":"only partially drawn"));
		fprintf(stderr,"\n");
		if (bNotDrawn)
			return;
		else
			s = sqrt(PI/MAX_NUM_PIX)*r; /* this will cause transparency */
	}

	dLimit = r*MAG(p->Light.vPos)*(1 - p->dLightInten);

	/* Get unit vector from sphere centre to camera */

	SUB_VEC(p->Camera.vPos,pos,v);
	if (!(m = MAG(v)))
		return; /* camera and particle position coincide */
	NORM_VEC(v,m);

	/*
	** Load matrix for rotational transformation.
	**
	** Strategy: only top hemisphere is drawn, hence rotate (0,0,1) to
	** point directly at camera. Then construct basis such that (1,0,0)
	** and (0,1,0) rotate to (arbitrary) mutually-orthonormal points.
	*/

	SET_VEC(a,v[Y],-v[X],0); /* orthogonal to (v[X],v[Y],v[Z]) */
	if (!(m = MAG(a))) {
		UNIT_MAT(xform); /* exactly head-on case */
	}
	else {
		NORM_VEC(a,m);
		SET_VEC(b,-v[X]*v[Z],-v[Y]*v[Z],SQ(v[X])+SQ(v[Y])); /* a x v */
		NORM_VEC(b,MAG(b));

		SET_VEC(xform[X],a[X],b[X],v[X]);
		SET_VEC(xform[Y],a[Y],b[Y],v[Y]);
		SET_VEC(xform[Z],a[Z],b[Z],v[Z]);
	}

	dstheta = 0.65*s/r; /* factor of 0.65 just to be safe */

	for (stheta=0.0;stheta<=1;stheta+=dstheta) { /* only 1 hemisphere */
		ctheta = sqrt(1.0 - SQ(stheta));
		rr = r*stheta;
		dphi = 0.65*(rr == 0.0 ? TWO_PI : s/rr); /* factor of 0.65 to be safe */
		for (phi=0;phi<TWO_PI;phi+=dphi) {
			/* get location on sphere */
			SET_VEC(spos,stheta*cos(phi),stheta*sin(phi),ctheta);
			SCALE_VEC(spos,r);
			/* transform to hemisphere visible by observer */
			Transform(xform,spos,v);
			COPY_VEC(v,spos);
			ADD_VEC(pos,spos,v);
			/* get projection */
			if (project(p,v,&ix,&iy) == FALSE)
				continue;
			/* Determine shading */
			SUB_VEC(p->Light.vPos,pos,v);
			if (dLimit*((double) random())/((double) RAND_MAX) <= DOT(v,spos))
				DrawPoint(image,ix,iy,(COLOR) color);
			else
				DrawPoint(image,ix,iy,BLACK);
			}
		}
	}

static int draw_color(const PARAMS *p,int iColor)
{
	if ((iColor < 0 && !p->bDoWalls) || iColor >= NUM_COLORS) {
		if (iColor >= NUM_COLORS && p->b24bit)
			fprintf(stderr,"draw_color(): 24-bit color (%i) only supported for POV-Ray output.\n",iColor);
		else {
			fprintf(stderr,"draw_color(): invalid particle color (%i).",iColor);
			if (iColor >= NUM_COLORS)
				fprintf(stderr," Try using '-c' option.");
			fprintf(stderr,"\n");
			}
		exit(1);
		}
	if (iColor < 0) { /* particle stuck to wall! (inherits color) */
		int iWallID = -1 - iColor;
		if (iWallID >= p->nWalls) {
			fprintf(stderr,"draw_color(): particle color %i incompatible with wall data.\n",iColor);
			exit(1);
			}
		iColor = p->pWalls[iWallID].iColor;
		}
	return iColor;
	}

static void draw_object(const PARAMS *p,IMAGE_T *image,VECTOR pos,
						double r,SHAPE Shape,int color)
{
	switch (Shape) {
	case Point:
		draw_point(p,image,pos,draw_color(p,color));
		break;
	case SolidSphere:
		draw_solid(p,image,pos,r*p->dRScale,draw_color(p,color));
		break;
	case POV: /* handled in pov_draw() */
	default:
		assert(0);
		}
	}

static int zcmpr(const void *e1,const void *e2)
{
	return ((const ZSort *)e1)->vdotz < ((const ZSort *)e2)->vdotz ? -1 :
		   ((const ZSort *)e1)->vdotz > ((const ZSort *)e2)->vdotz ?  1 : 0;
}

static void z_sort(const PARAMS *p,const SSDATA *d,int n,int **zindex)
{
	int i;

	*zindex = (int *) Array((size_t) n,sizeof(int));
	for (i=0;i<n;i++)
		(*zindex)[i] = i;
	if (p->bZSort) {
		ZSort *tbl;
		VECTOR v;

		tbl = (ZSort *) malloc((size_t) n*sizeof(ZSort));
		assert(tbl != NULL);
		for (i=0;i<n;i++) {
			tbl[i].zindex = (*zindex)[i];
			SUB_VEC(p->Camera.vPos,d[i].pos,v);
			tbl[i].vdotz = DOT(v,p->vZHat);
			}
		qsort((void *) tbl,(size_t) n,sizeof(ZSort),zcmpr);
		for (i=0;i<n;i++)
			(*zindex)[i] = tbl[i].zindex;
		free((void *) tbl);
		}
	}

static void make_node(const VECTOR vCenter,double dSize,NODE **pNode)
{
	int i;

	assert(dSize > 0.0);
	*pNode = (NODE *) malloc(sizeof(NODE));
	assert(*pNode != NULL);

	COPY_VEC(vCenter,(*pNode)->vCenter);
	(*pNode)->dSize = dSize;
	(*pNode)->dHalfSize = 0.5*dSize;

	for (i=0;i<CELLS_PER_NODE;i++) {
		(*pNode)->pCell[i] = NULL;
		(*pNode)->pLeaf[i] = NULL;
		}
	}

static void add_to_tree(NODE *pNode,const SSDATA *d)
{
	/* note: assumes particle inside node! */

	int i,idx,idy,idz;

	assert(pNode != NULL && d != NULL);

	idx = (d->pos[X] < pNode->vCenter[X] ? -1 : 1);
	idy = (d->pos[Y] < pNode->vCenter[Y] ? -1 : 1);
	idz = (d->pos[Z] < pNode->vCenter[Z] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1));

	if (pNode->pCell[i] != NULL)
		add_to_tree(pNode->pCell[i],d);
	else if (pNode->pLeaf[i] != NULL) {
		VECTOR v;
		SET_VEC(v,idx,idy,idz);
		SCALE_VEC(v,0.5*pNode->dHalfSize);
		ADD_VEC(v,pNode->vCenter,v);
		make_node(v,pNode->dHalfSize,&pNode->pCell[i]);
		add_to_tree(pNode->pCell[i],pNode->pLeaf[i]);
		add_to_tree(pNode->pCell[i],d);
		pNode->pLeaf[i] = NULL; /* for completeness */
		}
	else
		pNode->pLeaf[i] = d;
	}

static double get_pot(const NODE *pNode,const SSDATA *d,double dTheta)
{
	VECTOR v;
	double dPot,r;
	int i;

	assert(pNode != NULL && d != NULL);

	dPot = 0.0;
	for (i=0;i<CELLS_PER_NODE;i++)
		if (pNode->pCell[i] != NULL) {
			SUB_VEC(pNode->pCell[i]->vComPos,d->pos,v);
			r = MAG(v);
			if (r < d->radius)
				r = d->radius; /* pseudo-softening */
			if (pNode->pCell[i]->dSize/r > dTheta)
				dPot += get_pot(pNode->pCell[i],d,dTheta);
			else
				dPot -= pNode->pCell[i]->dMass/r;
			}
		else if (pNode->pLeaf[i] != NULL && pNode->pLeaf[i] != d) {
			SUB_VEC(pNode->pLeaf[i]->pos,d->pos,v);
			r = MAG(v);
			if (r < pNode->pLeaf[i]->radius + d->radius)
				r = pNode->pLeaf[i]->radius + d->radius;
			dPot -= pNode->pLeaf[i]->mass/r;
		}

	return dPot;
	}

static void calc_moments(NODE *pNode)
{
	VECTOR v;
	int i;

	assert(pNode != NULL);

	pNode->dMass = 0.0;
	ZERO_VEC(pNode->vComPos);
	for (i=0;i<CELLS_PER_NODE;i++)
		if (pNode->pCell[i] != NULL) {
			calc_moments(pNode->pCell[i]);
			pNode->dMass += pNode->pCell[i]->dMass;
			COPY_VEC(pNode->pCell[i]->vComPos,v);
			SCALE_VEC(v,pNode->pCell[i]->dMass);
			ADD_VEC(pNode->vComPos,v,pNode->vComPos);
			}
		else if (pNode->pLeaf[i] != NULL) {
			pNode->dMass += pNode->pLeaf[i]->mass;
			COPY_VEC(pNode->pLeaf[i]->pos,v);
			SCALE_VEC(v,pNode->pLeaf[i]->mass);
			ADD_VEC(pNode->vComPos,v,pNode->vComPos);
			}

	assert(pNode->dMass > 0.0);
	NORM_VEC(pNode->vComPos,pNode->dMass);
	}

static void kill_node(NODE *pNode)
{
	int i;

	assert(pNode != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (pNode->pCell[i] != NULL)
			kill_node(pNode->pCell[i]);

	free((void *) pNode);
	}

static void get_max_cell(const SSDATA *d,int n,VECTOR vPos,double *dMax,int bAdjForRad)
{
	VECTOR vMin,vMax,dv;

	double c;
	int i,k;

	SET_VEC(vMin,HUGE_VAL,HUGE_VAL,HUGE_VAL);
	SET_VEC(vMax,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
	for (i=0;i<n;i++)
		for (k=0;k<N_DIM;k++) {
			c = d[i].pos[k];
			if (bAdjForRad)
				c -= d[i].radius;
			if (c < vMin[k])
				vMin[k] = c;
			c = d[i].pos[k];
			if (bAdjForRad)
				c += d[i].radius;
			if (c > vMax[k])
				vMax[k] = c;
			}
	SUB_VEC(vMax,vMin,dv);
	assert(dv[X] >= 0.0 && dv[Y] >= 0.0 && dv[Z] >= 0.0);
	*dMax = (dv[X] > dv[Y] ? (dv[X] > dv[Z] ? dv[X] : dv[Z]) :
			 (dv[Y] > dv[Z] ? dv[Y] : dv[Z]));
	SCALE_VEC(dv,0.5);
	ADD_VEC(vMin,dv,vPos);
	}

static int mass_sort_compar(const void *vArg1,const void *vArg2)
{
	const struct mass_sort_s *sArg1 = (const struct mass_sort_s *) vArg1;
	const struct mass_sort_s *sArg2 = (const struct mass_sort_s *) vArg2;

	if (sArg1->dMass < sArg2->dMass)
		return -1;

	if (sArg1->dMass > sArg2->dMass)
		return 1;

	return 0;
	}

static void get_pot_pos(PARAMS *p,const SSDATA *d,int n,int nSample)
{
	/*
	** Uses tree to find gravitational potential minimum of system
	** (for light/camera/view positioning).
	*/

	NODE *pRoot = NULL;
	VECTOR v,vRootCenter;
	double dRootSize,dPot,dPotTot;
	int i,*iSampleIndex;

	printf("Computing potential minimum...\n");

	/*
	** First compute root cell location and size that contains all
	** particles.
	*/

	get_max_cell(d,n,vRootCenter,&dRootSize,FALSE);

	printf("Root cell: center = (%.1g,%.1g,%.1g), size = %.1g\n",
		   vRootCenter[X],vRootCenter[Y],vRootCenter[Z],dRootSize);

	/* make the root node */

	make_node(vRootCenter,dRootSize,&pRoot);

	/* fill the tree with particles */

	for (i=0;i<n;i++)
		add_to_tree(pRoot,&d[i]);

	/* compute tree moments */

	calc_moments(pRoot);

	printf("System: center of mass = (%.1g,%.1g,%.1g), mass = %.1g\n",
		   pRoot->vComPos[X],pRoot->vComPos[Y],pRoot->vComPos[Z],pRoot->dMass);

	/* construct list of particles to use as sample points */

	iSampleIndex = (int *) malloc((size_t) nSample*sizeof(int));
	assert(iSampleIndex != NULL);
	switch (p->iSampleOrder) {
	case PotRan:
	  {
		  int *bChosen;
		  bChosen = (int *) malloc((size_t) n*sizeof(int));
		  assert(bChosen != NULL);
		  for (i=0;i<n;i++)
			  bChosen[i] = FALSE; /* not needed with malloc() */
		  for (i=0;i<nSample;i++) {
			  do {
				  iSampleIndex[i] = rand()%n;
				  } while (bChosen[iSampleIndex[i]]);
			  bChosen[iSampleIndex[i]] = TRUE;
			  }
		  free((void *) bChosen);
		  break;
		  }
	case PotOrd:
		for (i=0;i<nSample;i++)
			iSampleIndex[i] = i;
		break;
	case PotLrg:
	  {
		  struct mass_sort_s *sMassSort;
		  sMassSort = (struct mass_sort_s *) malloc((size_t) n*sizeof(struct mass_sort_s));
		  assert(sMassSort != NULL);
		  for (i=0;i<n;i++) {
			  sMassSort[i].iIndex = i;
			  sMassSort[i].dMass = d[i].mass;
			  }
		  qsort((void *) sMassSort,(size_t) n,sizeof(struct mass_sort_s),mass_sort_compar);
		  for (i=0;i<nSample;i++)
			  iSampleIndex[i] = sMassSort[i].iIndex;
		  free((void *) sMassSort);
		  break;
		  }
	default:
		assert(0);
		}

	/*
	** Now compute gravitational potential "moment" based on potential
	** at particle positions.
	*/

	dPotTot = 0.0;
	ZERO_VEC(p->vPotPos);
	for (i=0;i<nSample;i++) {
		dPot = get_pot(pRoot,&d[iSampleIndex[i]],p->dTheta);
		dPotTot += dPot;
		COPY_VEC(d[iSampleIndex[i]].pos,v);
		SCALE_VEC(v,dPot);
		ADD_VEC(p->vPotPos,v,p->vPotPos);
		}

	assert(dPotTot < 0.0);

	NORM_VEC(p->vPotPos,dPotTot);

	printf("Grav. pot.: wgt center = (%.1g,%.1g,%.1g), sum potential = %.1g\n",
		   p->vPotPos[X],p->vPotPos[Y],p->vPotPos[Z],dPotTot);

	/* free storage */

	free((void *) iSampleIndex);

	kill_node(pRoot); /* deallocate tree */
	}

static void get_lrg(const SSDATA *d,int n,int *iID,double *dMass,double *dRadius,VECTOR vPos)
{
	/* returns position and "radius" of largest (most massive) particle/aggregate */

	struct agg_s {
	  int iOrgIdx;
	  int nPart;
	  double dMass;
	  VECTOR vPos;
	  } *agg;

	VECTOR v;
	int i,j,nAggs,bAggFound;

	/* examine individual particles first */

	assert(d != NULL && n > 0);
	*iID = 0;
	*dMass = d->mass;
	for (i=1;i<n;i++)
		if (d[i].mass > *dMass) {
			*iID = i;
			*dMass = d[i].mass;
			}
	*dRadius = d[*iID].radius;
	COPY_VEC(d[*iID].pos,vPos);
	printf("Largest particle: idx=%i, mass=%g, radius=%g, pos=(%g,%g,%g)\n",
		   *iID,*dMass,*dRadius,vPos[X],vPos[Y],vPos[Z]);

	/* now consider aggregates, if any */

	agg = malloc((size_t) n*sizeof(struct agg_s));
	assert(agg != NULL);

	for (i=nAggs=0;i<n;i++)
		if (d[i].org_idx < 0) {
			bAggFound = 0;
			for (j=0;j<nAggs;j++)
				if (d[i].org_idx == agg[j].iOrgIdx) {
					++agg[j].nPart;
					agg[j].dMass += d[i].mass;
					COPY_VEC(d[i].pos,v);
					SCALE_VEC(v,d[i].mass);
					ADD_VEC(agg[j].vPos,v,agg[j].vPos);
					bAggFound = 1;
					break;
					}
			if (!bAggFound) {
				agg[nAggs].iOrgIdx = d[i].org_idx;
				agg[nAggs].nPart = 1;
				agg[nAggs].dMass = d[i].mass;
				COPY_VEC(d[i].pos,agg[nAggs].vPos);
				SCALE_VEC(agg[nAggs].vPos,agg[nAggs].dMass);
				++nAggs;
				}
			}

	if (nAggs > 0)
		printf("%i aggregate%s found.\n",nAggs,nAggs==1?"":"s");

	/* see if any of the aggregates have larger mass than the largest particle */

	for (i=0;i<nAggs;i++)
		if (agg[i].dMass > *dMass) {
			*iID = -1 - i;
			*dMass = agg[i].dMass;
			}

	/* if so, return its position and estimated (maximum) "radius" */

	if (*iID < 0) {
		double r;
		COPY_VEC(agg[-1 - *iID].vPos,vPos);
		NORM_VEC(vPos,*dMass); /* center-of-mass position */
		*dRadius = 0.0;
		for (i=0;i<n;i++)
			if (d[i].org_idx == agg[-1 - *iID].iOrgIdx) {
				SUB_VEC(d[i].pos,vPos,v);
				r = MAG(v) + d[i].radius;
				if (r > *dRadius)
					*dRadius = r;
				}
		printf("Largest aggregate: idx=%i, N=%i, mass=%g, radius=%g, pos=(%g,%g,%g)\n",
			   agg[-1 - *iID].iOrgIdx,agg[-1 - *iID].nPart,*dMass,*dRadius,
			   vPos[X],vPos[Y],vPos[Z]);
		*iID = agg[-1 - *iID].iOrgIdx;
		}

	free((void *) agg);
	}

static void get_com_by_COL(const PARAMS *p,const SSDATA *d,int n,int iColor,VECTOR vComPos)
{
	VECTOR v;
	double m;
	int i;

	ZERO_VEC(vComPos);
	m = 0.0;
	for (i=0;i<n;i++) {
		if (iColor != 0 && draw_color(p,d[i].color) != iColor)
			continue;
		m += d[i].mass;
		COPY_VEC(d[i].pos,v);
		SCALE_VEC(v,d[i].mass);
		ADD_VEC(vComPos,v,vComPos);
		}
	assert(m >= 0.0);
	if (m == 0.0) {
		fprintf(stderr,"No particles of color %i found...reverting to origin.\n",iColor);
		return;
		}
	NORM_VEC(vComPos,m);
	}

static void get_com_by_OID(const PARAMS *p,const SSDATA *d,int n,int iOrgIdx,VECTOR vComPos,double *dRadius)
{
	VECTOR v;
	double m,r;
	int i;

	ZERO_VEC(vComPos);
	m = 0.0;
	for (i=0;i<n;i++)
		if (d[i].org_idx == iOrgIdx) {
			m += d[i].mass;
			COPY_VEC(d[i].pos,v);
			SCALE_VEC(v,d[i].mass);
			ADD_VEC(vComPos,v,vComPos);
			}
	assert(m >= 0.0);
	if (m == 0.0) {
		fprintf(stderr,"No particles of original index %i found...\n"
				"Defaulting to system center of mass.\n",iOrgIdx);
		get_com_by_COL(p,d,n,0,vComPos);
		*dRadius = 0.0;
		return;
		}

	NORM_VEC(vComPos,m);

	/* get "radius" estimate */

	*dRadius = 0.0;
	for (i=0;i<n;i++)
		if (d[i].org_idx == iOrgIdx) {
			SUB_VEC(d[i].pos,vComPos,v);
			r = MAG(v) + d[i].radius;
			if (r > *dRadius)
				*dRadius = r;
			}

	assert(*dRadius > 0.0);
	}

static void set_scene(PARAMS *p,const SSDATA *d,int n,int nMoveFrames,double dMinMass)
{
	static BOOLEAN bFirstCall = TRUE;
	static int nWarnings = 0;

	VECTOR v,vPos,vGoTo;
	BOOLEAN bTargetChanged;
	double dMass,dRadius,dNewViewSize,dv;
	int iID;

	if (bFirstCall) {
		p->Light.Target.nFramesActive = p->Camera.Target.nFramesActive = p->LookAt.Target.nFramesActive = p->ViewSizeTarget.nFramesActive = nMoveFrames;
		p->Light.Target.iID = p->Camera.Target.iID = p->LookAt.Target.iID = p->ViewSizeTarget.iID = INT_MAX;
		p->Light.Target.bIsCOM = p->Camera.Target.bIsCOM = p->LookAt.Target.bIsCOM = p->ViewSizeTarget.bIsCOM = p->bStartAtCOM;
		}

	/* do camera first because it may be needed for light source position */

	dMass = dRadius = 0.0; /* reset for LRG,PID,OID options */

	bTargetChanged = FALSE;

	switch (p->Camera.iOpt) {
	case CamFIX:
		COPY_VEC(p->Camera.vOrig,vGoTo);
		if (!bFirstCall && p->Camera.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Camera.Target.bIsCOM = FALSE;
			}
		break;
	case CamLRG:
		assert(dMass == 0.0); /* should not be assigned yet */
		get_lrg(d,n,&iID,&dMass,&dRadius,vPos);
		assert(dMass > 0.0 && dRadius > 0.0);
		if (dMass >= dMinMass) {
			COPY_VEC(p->Camera.vOrig,v);
			if (p->bScaleLengthByRadius)
				SCALE_VEC(v,dRadius);
			ADD_VEC(vPos,v,vGoTo); /* does not account for particle/aggregate rotation */
			if (iID != p->Camera.Target.iID) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->Camera.Target.iID = iID;
				}
			}
		else { /* default to system center of mass if target not big enough */
			if (p->Camera.Target.iID != INT_MAX || !p->Camera.Target.bIsCOM) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->Camera.Target.iID = INT_MAX;
				p->Camera.Target.bIsCOM = TRUE;
				}
			}
		break;
	case CamPOT:
		ADD_VEC(p->vPotPos,p->Camera.vOrig,vGoTo);
		if (!bFirstCall && p->Camera.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Camera.Target.bIsCOM = FALSE;
			}
		break;
	case CamCOL:
		break; /* handled below in camera motion section */
	case CamPID:
		COPY_VEC(p->Camera.vOrig,v);
		assert(p->Camera.iLock < n);
		dRadius = d[p->Camera.iLock].radius;
		if (p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(d[p->Camera.iLock].pos,v,vGoTo);
		if (d[p->Camera.iLock].org_idx != p->Camera.Target.iID) {
			if (!bFirstCall)
				bTargetChanged = TRUE;
			p->Camera.Target.iID = d[p->Camera.iLock].org_idx;
			}
		break;
	case CamOID:
		get_com_by_OID(p,d,n,p->Camera.iLock,vGoTo,&dRadius);
		assert(dRadius >= 0.0);
		COPY_VEC(p->Camera.vOrig,v);
		if (dRadius > 0.0 && p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(vGoTo,v,vGoTo);
		if (!bFirstCall && p->Camera.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Camera.Target.bIsCOM = FALSE;
			}
		break;
	default:
		assert(0);
		}

	/* camera motion */

	if ((p->Camera.Target.bIsCOM && p->Camera.Target.iID == INT_MAX) || p->Camera.iOpt == CamCOL) {
		get_com_by_COL(p,d,n,p->Camera.iLock,vGoTo);
		ADD_VEC(vGoTo,p->Camera.vOrig,vGoTo);
		}
	if (bTargetChanged)
		p->Camera.Target.nFramesActive = 0;
	SUB_VEC(vGoTo,p->Camera.vPos,v);
	if (p->Camera.Target.nFramesActive < nMoveFrames) {
		NORM_VEC(v,nMoveFrames - p->Camera.Target.nFramesActive);
		++p->Camera.Target.nFramesActive;
		}
	ADD_VEC(p->Camera.vPos,v,p->Camera.vPos);

	bTargetChanged = FALSE;

	switch (p->Light.iOpt) {
	case LightFIX:
		COPY_VEC(p->Light.vOrig,vGoTo);
		if (!bFirstCall && p->Light.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Light.Target.bIsCOM = FALSE;
			}
		break;
	case LightLRG:
		if (dMass == 0.0) {
			get_lrg(d,n,&iID,&dMass,&dRadius,vPos);
			assert(dMass > 0.0 && dRadius > 0.0);
			}
		if (dMass >= dMinMass) {
			COPY_VEC(p->Light.vOrig,v);
			if (p->bScaleLengthByRadius)
				SCALE_VEC(v,dRadius);
			ADD_VEC(vPos,v,vGoTo);
			if (iID != p->Light.Target.iID) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->Light.Target.iID = iID;
				}
			}
		else {
			if (p->Light.Target.iID != INT_MAX || !p->Light.Target.bIsCOM) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->Light.Target.iID = INT_MAX;
				p->Light.Target.bIsCOM = TRUE;
				}
			}
		break;
	case LightPOT:
		ADD_VEC(p->vPotPos,p->Light.vOrig,vGoTo);
		if (!bFirstCall && p->Light.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Light.Target.bIsCOM = FALSE;
			}
		break;
	case LightCOL:
		break;
	case LightPID:
		COPY_VEC(p->Light.vOrig,v);
		assert(p->Light.iLock < n);
		dRadius = d[p->Light.iLock].radius;
		if (p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(d[p->Light.iLock].pos,v,vGoTo);
		if (d[p->Light.iLock].org_idx != p->Light.Target.iID) {
			if (!bFirstCall)
				bTargetChanged = TRUE;
			p->Light.Target.iID = d[p->Light.iLock].org_idx;
			}
		break;
	case LightOID:
		get_com_by_OID(p,d,n,p->Light.iLock,vGoTo,&dRadius);
		assert(dRadius >= 0.0);
		COPY_VEC(p->Light.vOrig,v);
		if (dRadius > 0.0 && p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(vGoTo,v,vGoTo);
		if (!bFirstCall && p->Light.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->Light.Target.bIsCOM = FALSE;
			}
		break;
	case LightCAM:
		ADD_VEC(p->Camera.vPos,p->Light.vOrig,vGoTo);
		COPY_VEC(vGoTo,p->Light.vPos);
		if (p->Light.Target.bIsCOM)
			p->Light.Target.bIsCOM = FALSE; /* no inertia relative to camera */
		break;
	default:
		assert(0);
		}

	if ((p->Light.Target.bIsCOM && p->Light.Target.iID == INT_MAX) || p->Light.iOpt == LightCOL) {
		get_com_by_COL(p,d,n,p->Light.iLock,vGoTo);
		ADD_VEC(vGoTo,p->Light.vOrig,vGoTo);
		}
	if (bTargetChanged)
		p->Light.Target.nFramesActive = 0;
	SUB_VEC(vGoTo,p->Light.vPos,v);
	if (p->Light.Target.nFramesActive < nMoveFrames) {
		NORM_VEC(v,nMoveFrames - p->Light.Target.nFramesActive);
		++p->Light.Target.nFramesActive;
		}
	ADD_VEC(p->Light.vPos,v,p->Light.vPos);

	bTargetChanged = FALSE;

	dNewViewSize = p->dViewSize;

	switch (p->LookAt.iOpt) {
	case ViewFIX:
		COPY_VEC(p->LookAt.vOrig,vGoTo);
		if (!bFirstCall && p->LookAt.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->LookAt.Target.bIsCOM = FALSE;
			}
		break;
	case ViewLRG:
		if (dMass == 0.0) {
			get_lrg(d,n,&iID,&dMass,&dRadius,vPos);
			assert(dMass > 0.0 && dRadius > 0.0);
			}
		if (dMass >= dMinMass) {
			COPY_VEC(p->LookAt.vOrig,v);
			if (p->bScaleLengthByRadius)
				SCALE_VEC(v,dRadius);
			ADD_VEC(vPos,v,vGoTo);
			if (iID != p->LookAt.Target.iID) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->LookAt.Target.iID = iID;
				}
			if (p->dViewTargetZoomFactor > 0.0)
				dNewViewSize = 2.0*dRadius/p->dViewTargetZoomFactor;
			}
		else {
			if (p->LookAt.Target.iID != INT_MAX || !p->LookAt.Target.bIsCOM) {
				if (!bFirstCall)
					bTargetChanged = TRUE;
				p->LookAt.Target.iID = INT_MAX;
				p->LookAt.Target.bIsCOM = TRUE;
				}
			}
		break;
	case ViewPOT:
		ADD_VEC(p->vPotPos,p->LookAt.vOrig,vGoTo);
		if (!bFirstCall && p->LookAt.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->LookAt.Target.bIsCOM = FALSE;
			}
		break;
	case ViewCOL:
		break;
	case ViewPID:
		COPY_VEC(p->LookAt.vOrig,v);
		assert(p->LookAt.iLock < n);
		dRadius = d[p->LookAt.iLock].radius;
		if (p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(d[p->LookAt.iLock].pos,v,vGoTo);
		if (d[p->LookAt.iLock].org_idx != p->LookAt.Target.iID) {
			if (!bFirstCall)
				bTargetChanged = TRUE;
			p->LookAt.Target.iID = d[p->LookAt.iLock].org_idx;
			}
		if (p->dViewTargetZoomFactor > 0.0)
			dNewViewSize = 2.0*dRadius/p->dViewTargetZoomFactor;
		break;
	case ViewOID:
		get_com_by_OID(p,d,n,p->LookAt.iLock,vGoTo,&dRadius);
		assert(dRadius >= 0.0);
		COPY_VEC(p->LookAt.vOrig,v);
		if (dRadius > 0.0 && p->bScaleLengthByRadius)
			SCALE_VEC(v,dRadius);
		ADD_VEC(vGoTo,v,vGoTo);
		if (!bFirstCall && p->LookAt.Target.bIsCOM) {
			bTargetChanged = TRUE;
			p->LookAt.Target.bIsCOM = FALSE;
			}
		if (dRadius > 0.0 && p->dViewTargetZoomFactor > 0.0)
			dNewViewSize = 2.0*dRadius/p->dViewTargetZoomFactor;
		break;
	case ViewDIR:
		COPY_VEC(p->LookAt.vOrig,p->vZHat);
		SCALE_VEC(p->vZHat,p->dFocusDist);
		ADD_VEC(p->Camera.vPos,p->vZHat,vGoTo);
		if (p->LookAt.Target.bIsCOM)
			p->LookAt.Target.bIsCOM = FALSE; /* no inertia */
		break;
	default:
		assert(0);
		}

	if ((p->LookAt.Target.bIsCOM && p->LookAt.Target.iID == INT_MAX) || p->LookAt.iOpt == ViewCOL) {
		get_com_by_COL(p,d,n,p->LookAt.iLock,vGoTo);
		ADD_VEC(vGoTo,p->LookAt.vOrig,vGoTo);
		}
	if (bTargetChanged)
		p->LookAt.Target.nFramesActive = 0;
	SUB_VEC(vGoTo,p->LookAt.vPos,v);
	if (p->LookAt.Target.nFramesActive < nMoveFrames) {
		NORM_VEC(v,nMoveFrames - p->LookAt.Target.nFramesActive);
		++p->LookAt.Target.nFramesActive;
		}
	ADD_VEC(p->LookAt.vPos,v,p->LookAt.vPos);

	/* handle any zoom effects (note: camera/light positions not changed) */

	if (bTargetChanged)
		p->ViewSizeTarget.nFramesActive = 0;
	dv = dNewViewSize - p->dOldViewSize;
	if (p->ViewSizeTarget.nFramesActive < nMoveFrames) {
		dv /= (nMoveFrames - p->ViewSizeTarget.nFramesActive);
		++p->ViewSizeTarget.nFramesActive;
		}
	p->dViewSize = (p->dOldViewSize += dv);

	/* set remaining view parameters for scene */

	SUB_VEC(p->LookAt.vPos,p->Camera.vPos,p->vZHat);
	assert(MAG(p->vZHat) > 0.0); /* camera can't look at itself!... */
	if (p->dFocusDist > 0.0) { /* shorten or lengthen focus as appropriate */
		double dScale = p->dFocusDist/MAG(p->vZHat);
		SCALE_VEC(p->vZHat,dScale);
		}
	p->dViewRatio = MAG(p->vZHat)/p->dViewSize;
	if (p->dViewRatio < 1.0 && nWarnings < MAX_NUM_WARNINGS) {
		fprintf(stderr,"WARNING: View may be distorted.\n");
		if (++nWarnings == MAX_NUM_WARNINGS)
			fprintf(stderr,"(Further warnings suppressed.)\n");
		}
	p->dPixelSize = p->dViewSize/p->iFrameSize;
	COPY_VEC(p->vSkyVec,p->vYHat);
	CROSS(p->vYHat,p->vZHat,p->vXHat);
	if (MAG_SQ(p->vXHat) == 0.0) {
		fprintf(stderr,"ERROR: View direction and sky vector are parallel.\n");
		exit(1);
		}
	SCALE_VEC(p->vXHat,-1); /* left-hand to right-hand transformation */
	MakeBasis(p->vZHat,p->vXHat,p->vYHat);

	bFirstCall = FALSE;
	}

static void check_view(PARAMS *p,SSDATA *d,int n,int w,int h,int nMoveFrames,double dMinMass)
{
	BOOLEAN reset;
	int i;

	do {
		set_scene(p,d,n,nMoveFrames,dMinMass);
		reset = FALSE;
		if (p->dViewScaleCrit > 0.0) {
			int n_out,ix,iy;
			/* quick check to see if view should be zoomed out */
			n_out = 0;
			for (i=0;i<n;i++)
				if (project(p,d[i].pos,&ix,&iy) == FALSE ||
					ix < 0 || ix >= w || iy < 0 || iy >= h)
					++n_out;
			if ((double) n_out/n > p->dViewScaleCrit) {
				VECTOR v1,v2;
				COPY_VEC(p->LookAt.vPos,v1);
				if (p->dFocusDist > 0.0) {
					SCALE_VEC(v1,p->dFocusDist);
					ADD_VEC(p->Camera.vPos,v1,v1);
					}
				SUB_VEC(p->Light.vPos,v1,v2);
				SCALE_VEC(v2,p->dViewScale);
				ADD_VEC(v1,v2,p->Light.vPos);
				SUB_VEC(p->Camera.vPos,v1,v2);
				SCALE_VEC(v2,p->dViewScale);
				ADD_VEC(v1,v2,p->Camera.vPos);
				p->dViewSize *= p->dViewScale;
				reset = TRUE;
				}
			}
		} while (reset);
	}

static void pov_draw(PARAMS *p,SSHEAD *h,SSDATA *d,MATRIX *o,const char *outfile,int nMoveFrames,double dMinMass,int iTransColor)
{
	SSDATA *dd;
	WALL_DATA *w;
	FILE *fp;
	VECTOR vParPos;
	double dRadius;
	double dScale=0.0;/*DEBUG! doesn't work with patch drawing below??*/
	int *zindex,i,c,n;

	dRadius = 0.0; /* to keep gcc -Wall happy */
	ZERO_VEC(vParPos);
	n = (h == NULL ? 0 : h->n_data);
	check_view(p,d,n,p->iFrameSize,p->iFrameSize,nMoveFrames,dMinMass);
	if (n > 0)
		z_sort(p,d,n,&zindex);
	if (p->bRenorm) {
		/*
		** Strategy: Since POV-Ray only renders in single precision,
		** scale everything to the distance between the camera and the
		** view point.
		*/
		VECTOR vOffset;
		COPY_VEC(p->LookAt.vPos,vOffset);
		if (p->dFocusDist > 0.0) {
			SCALE_VEC(vOffset,p->dFocusDist);
			ADD_VEC(p->Camera.vPos,vOffset,vOffset); /* the new "look-at" position */
			dScale = p->dFocusDist;
			}
		else {
			VECTOR v;
			SUB_VEC(p->Camera.vPos,vOffset,v);
			dScale = MAG(v);
			}
		assert(dScale > 0.0);
		dScale = 1.0/dScale;
		for (i=0;i<n;i++) {
			SUB_VEC(d[i].pos,vOffset,d[i].pos);
			SCALE_VEC(d[i].pos,dScale);
			d[i].radius *= dScale;
			if (d[i].radius < MIN_POV_SIZE) {
				d[i].radius = MIN_POV_SIZE;
				pov_warn_min_size(p,"tiny particle");
				}
			}
		/*
		** In case view motion is enabled, save rescaled light and
		** camera positions in separate storage.
		*/
		SUB_VEC(p->Light.vPos,vOffset,p->Light.vPosNorm);
		SCALE_VEC(p->Light.vPosNorm,dScale);
		SUB_VEC(p->Camera.vPos,vOffset,p->Camera.vPosNorm);
		SCALE_VEC(p->Camera.vPosNorm,dScale);
		/* also scale view size */
		p->dViewSizeNorm = p->dViewSize*dScale;
		/* scale walls too (need to reread each time!) */
		for (i=0;i<p->nWalls;i++) {
			w = &p->pWalls[i];
			SUB_VEC(w->vOrigin,vOffset,w->vOrigin);
			SCALE_VEC(w->vOrigin,dScale);
			SCALE_VEC(w->vVertex1,dScale);
			SCALE_VEC(w->vVertex2,dScale);
			SCALE_VEC(w->vVel,dScale);
			w->dOscAmp *= dScale;
			w->dRadius *= dScale;
			w->dLength *= dScale;
			w->dHoleRadius *= dScale;
			}
		}

	/* now write POV info */

	if ((fp = fopen(outfile,"w")) == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",outfile);
		return;
		}
	pov_set_scene(p,fp);
	if (p->bUsePatch) {
		VECTOR vPos;
		double dLx=p->dLx,dLy=p->dLy;
		int ix,iy;
		if (p->bRenorm) {
			dLx *= dScale;
			dLy *= dScale;
			}
		if (p->bDrawPatchLines) {
			/*DEBUG not implemented
			double dy;
			for (ix=-p->nRepl;ix<=p->nRepl;ix++)
				for (iy=-p->nRepl;iy<=p->nRepl;iy++) {
					dy = SHEAR(ix,dLx,dLy,p->dOmega,h->time);
					pov_draw_square(p,fp,
									ix*dLx - 0.5*dLx,
									iy*dLy + dy - 0.5*dLy,
									ix*dLx + 0.5*dLx,
									iy*dLy + dy + 0.5*dLy,
									YELLOW);
					}
			*/
			}
		if (p->dBlob)
			pov_blob_open(p,fp);
		for (i=0;i<n;i++) {
			vPos[Z] = d[i].pos[Z];
			for (ix=-p->nRepl;ix<=p->nRepl;ix++) {
				vPos[X] = d[i].pos[X] + ix*dLx;
				for (iy=-p->nRepl;iy<=p->nRepl;iy++) {
					vPos[Y] = d[i].pos[Y] + iy*dLy +
						SHEAR(ix,dLx,dLy,p->dOmega,h->time);
					if (p->iGhostColor && (ix || iy))
						c = p->iGhostColor;
					else
						c = d[i].color;
					pov_draw_particle(p,fp,vPos,d[i].radius,c,iTransColor,POV_DRAW_REGULAR/*DEBUG*/,o[i]);
					}
				}
			}
		if (p->dBlob) pov_blob_close(fp);
		}
	else {
		int bFoundParticle = FALSE;
		for (i=0;i<p->nWalls;i++)
			pov_draw_wall(p,fp,&p->pWalls[i]);
		if (p->dBlob)
			pov_blob_open(p,fp);
		for (i=0;i<n;i++) {
			dd = &d[zindex[i]];
			if (p->iHighlightIndex >= 0 && dd->org_idx == p->iHighlightIndex) {
				COPY_VEC(dd->pos,vParPos);
				dRadius = dd->radius;
				pov_draw_particle(p,fp,dd->pos,dd->radius,dd->color,iTransColor,POV_DRAW_HIGHLIGHT,o[zindex[i]]);
				bFoundParticle = TRUE;
				}
			else if (bFoundParticle) {
				VECTOR vRelPos,z,vPerp;
				double dVdotZ;
				SUB_VEC(vParPos,dd->pos,vRelPos);
				COPY_VEC(p->vZHat,z);
				dVdotZ = DOT(vRelPos,z);
				SCALE_VEC(z,dVdotZ);
				SUB_VEC(vRelPos,z,vPerp);
				if (MAG(vPerp) < dRadius + dd->radius)					
					pov_draw_particle(p,fp,dd->pos,dd->radius,dd->color,iTransColor,POV_DRAW_TRANSPARENT,o[zindex[i]]);
				else
					pov_draw_particle(p,fp,dd->pos,dd->radius,dd->color,iTransColor,POV_DRAW_REGULAR,o[zindex[i]]);
				}
			else
				pov_draw_particle(p,fp,dd->pos,dd->radius,dd->color,iTransColor,POV_DRAW_REGULAR,o[zindex[i]]);
			}
		if (p->dBlob)
			pov_blob_close(fp);
		}
	fclose(fp);
	if (!p->bNoGZip) {
		/* attempt to gzip the POV file to save space... */
		char achSysStr[MAXPATHLEN];
		snprintf(achSysStr,MAXPATHLEN,"gzip -f %s",outfile);
		system(achSysStr);
		}
	if (n > 0)
		FreeArray(zindex);
	}

#define WALL_SCALING_FACTOR 0.225

static void draw_wall(const PARAMS *p,IMAGE_T *image,const WALL_DATA *w)
{
	VECTOR vn,vx,vy,vOrigin,v1,v2;
	int ix,iy;

	if (w->iType != WallPlane) return; /* for now */

	if (w->iColor == BLACK)
		return;

	/* connect corners of a square around wall origin, plus 8 surrounding squares */

	COPY_VEC(w->vOrient,vn);
	GetBasis(vn,vx,vy);
	SCALE_VEC(vx,WALL_SCALING_FACTOR*p->dViewSize);
	SCALE_VEC(vy,WALL_SCALING_FACTOR*p->dViewSize);

	for (ix=-1;ix<=1;ix++)
		for (iy=-1;iy<=1;iy++) {
			/* offset origin */
			COPY_VEC(w->vOrigin,vOrigin);
			COPY_VEC(vx,v1); /* temporary storage */
			SCALE_VEC(v1,2*ix);
			ADD_VEC(vOrigin,v1,vOrigin);
			COPY_VEC(vy,v1);
			SCALE_VEC(v1,2*iy);
			ADD_VEC(vOrigin,v1,vOrigin);
			/* now connect vertices */
			ADD_VEC(vOrigin,vx,v1);
			ADD_VEC(v1,vy,v1);
			ADD_VEC(vOrigin,vx,v2);
			SUB_VEC(v2,vy,v2);
			draw_line(p,image,v1,v2,w->iColor);
			SUB_VEC(vOrigin,vx,v1);
			SUB_VEC(v1,vy,v1);
			draw_line(p,image,v1,v2,w->iColor);
			SUB_VEC(vOrigin,vx,v2);
			ADD_VEC(v2,vy,v2);
			draw_line(p,image,v1,v2,w->iColor);
			ADD_VEC(vOrigin,vx,v1);
			ADD_VEC(v1,vy,v1);
			draw_line(p,image,v1,v2,w->iColor);
			/* connect diagonals */
			SUB_VEC(vOrigin,vx,v2);
			SUB_VEC(v2,vy,v2);
			draw_line(p,image,v1,v2,w->iColor);
			SUB_VEC(vOrigin,vx,v1);
			ADD_VEC(v1,vy,v1);
			ADD_VEC(vOrigin,vx,v2);
			SUB_VEC(v2,vy,v2);
			draw_line(p,image,v1,v2,w->iColor);
			}
	}

#undef WALL_SCALING_FACTOR

static void draw(PARAMS *p,SSHEAD *h,SSDATA *d,IMAGE_T *image,int nMoveFrames,double dMinMass,int iTransparentColor)
{
	SSDATA *dd;
	int *zindex,i,c,n;

	n = (h == NULL ? 0 : h->n_data);
	check_view(p,d,n,image->w,image->h,nMoveFrames,dMinMass);
	if (n > 0)
		z_sort(p,d,n,&zindex);
	if (p->bUsePatch) {
		VECTOR pos;
		int ix,iy;
		if (p->bDrawPatchLines && iTransparentColor >= 0) { /* negative transparent color indicates second pass through */
			double dy;
			for (ix=-p->nRepl;ix<=p->nRepl;ix++)
				for (iy=-p->nRepl;iy<=p->nRepl;iy++) {
					dy = SHEAR(ix,p->dLx,p->dLy,p->dOmega,h->time);
					draw_square(p,image,
								ix*p->dLx - 0.5*p->dLx,iy*p->dLy + dy - 0.5*p->dLy,
								ix*p->dLx + 0.5*p->dLx,iy*p->dLy + dy + 0.5*p->dLy,
								YELLOW);
					}
			}
		for (i=0;i<n;i++) {
			dd = &d[zindex[i]];
			pos[Z] = dd->pos[Z];
			for (ix=-p->nRepl;ix<=p->nRepl;ix++) {
				pos[X] = dd->pos[X] + ix*p->dLx;
				for (iy=-p->nRepl;iy<=p->nRepl;iy++) {
					pos[Y] = dd->pos[Y] + iy*p->dLy +
						SHEAR(ix,p->dLx,p->dLy,p->dOmega,h->time);
					if (p->iGhostColor && (ix || iy))
						c = p->iGhostColor;
					else
						c = dd->color;
					if (iTransparentColor == 0 || (iTransparentColor > 0 && draw_color(p,dd->color) == iTransparentColor) || (iTransparentColor < 0 && draw_color(p,dd->color) != -iTransparentColor))
						draw_object(p,image,pos,dd->radius,p->Shape,c);
					}
				}
			}
		}
	else {
		if (iTransparentColor >= 0)
			for (i=0;i<p->nWalls;i++)
				draw_wall(p,image,&p->pWalls[i]);
		for (i=0;i<n;i++) {
			dd = &d[zindex[i]];
			if (iTransparentColor == 0 || (iTransparentColor > 0 && draw_color(p,dd->color) == iTransparentColor) || (iTransparentColor < 0 && draw_color(p,dd->color) != -iTransparentColor))
				draw_object(p,image,dd->pos,dd->radius,p->Shape,dd->color);
			}
		}
	if (n > 0)
		FreeArray(zindex);
	}

static void get_view_size(PARAMS *p,int n,SSDATA *d)
{
	/*
	** Sets view size to maximum x/y/z particle separation and points
	** the camera at the middle. The camera position is scaled to give
	** a favorable view ratio.  The light position is similarly scaled.
	** This might do strange things if the camera and/or light positions
	** and/or view option are other than CamFIX, LightFIX, and ViewFIX.
	*/

	VECTOR v;

	if (p->LookAt.iOpt != ViewFIX)
		fprintf(stderr,"WARNING: View option may be incompatible with auto-sizing.\n");
	get_max_cell(d,n,p->LookAt.vOrig,&p->dViewSize,TRUE);
	printf("View size = %g\n",p->dViewSize);
	printf("Look at = %g %g %g\n",p->LookAt.vOrig[X],p->LookAt.vOrig[Y],p->LookAt.vOrig[Z]);
	if (p->Camera.iOpt != CamFIX)
		fprintf(stderr,"WARNING: Camera option may be incompatible with auto-sizing.\n");
	SUB_VEC(p->Camera.vOrig,p->LookAt.vOrig,v);
	if (MAG(v) == 0.0)
		fprintf(stderr,"WARNING: Camera is looking at itself!\n");
	else
		SCALE_VEC(p->Camera.vOrig,3.0*p->dViewSize/MAG(v));
	/* camera (& light) position should be relative to look-at point */
	ADD_VEC(p->Camera.vOrig,p->LookAt.vOrig,p->Camera.vOrig);
	printf("Camera position = %g %g %g\n",p->Camera.vOrig[X],p->Camera.vOrig[Y],p->Camera.vOrig[Z]);
	if (p->Light.iOpt != LightFIX)
		fprintf(stderr,"WARNING: Light source option may be incompatible with auto-sizing.\n");
	SUB_VEC(p->Light.vOrig,p->LookAt.vOrig,v);
	if (MAG(v) > 0.0)
		SCALE_VEC(p->Light.vOrig,3.0*p->dViewSize/MAG(v));
	ADD_VEC(p->Light.vOrig,p->LookAt.vOrig,p->Light.vOrig);
	printf("Light position = %g %g %g\n",p->Light.vOrig[X],p->Light.vOrig[Y],p->Light.vOrig[Z]);
	}

static int read_data(char *infile,SSHEAD *head,SSDATA **data)
{
	SSIO ssio;
	int i;

	if (ssioOpen(infile,&ssio,SSIO_READ)) {
		fprintf(stderr,"Unable to open \"%s\".\n",infile);
		return 1;
		}
	if (ssioHead(&ssio,head) || head->n_data < 0) {
		fprintf(stderr,"Corrupt header.\n");
		ssioClose(&ssio);
		return 1;
		}
	if (head->n_data == 0) {
		fprintf(stderr,"No particles found!\n");
		ssioClose(&ssio);
		return 1;
		}
	*data = (SSDATA *) malloc((size_t) head->n_data*sizeof(SSDATA));
	assert(*data != NULL);
	switch(head->iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		for (i=0;i<head->n_data;i++)
			if (ssioData(&ssio,&(*data)[i])) {
				fprintf(stderr,"Corrupt data.\n");
				free((void *) *data);
				ssioClose(&ssio);
				return 1;
				}
		break;
	case SSIO_MAGIC_REDUCED:
	    {
			SSRDATA dr;
			for (i=0;i<head->n_data;i++) {
				if (ssioDataReduced(&ssio,&dr)) {
					fprintf(stderr,"Corrupt data.\n");
					free((void *) *data);
					ssioClose(&ssio);
					return 1;
					}
				ssioReducedToStandard(&dr,&(*data)[i]);
				}
			}
		break;
	default:
		assert(0);
		}
	ssioClose(&ssio);
	return 0;
}

static int read_orient(char *infile,const SSHEAD *h,MATRIX **o)
{
	SSIO ssio;
	ORIHEAD head;
	union {
	  ORIDATA ori;
	  ORIRDATA orir;
	  } oridata;
	char orientfile[MAXPATHLEN];
	BOOLEAN bReduced;
	int i,k;

	snprintf(orientfile,MAXPATHLEN,"%s%s",infile,ORIENT_EXT);
	if (ssioOpen(orientfile,&ssio,SSIO_READ) != 0) {
		/* this may be the initial conditions, so allow it for now */
		fprintf(stderr,"Unable to open \"%s\"...assuming default orientations.\n",orientfile);
		if (*o != NULL) {
			free((void *) *o);
			*o = NULL;
			}
		return 0;
		}
	if (ssioHead(&ssio,&head) != 0) {
		fprintf(stderr,"Unable to read header info from \"%s\".\n",orientfile);
		goto abort;
		}
	if (head.n_data != h->n_data) {
		fprintf(stderr,"Number of particles in .ori file (%i) and ss file (%i) do not match.\n",head.n_data,h->n_data);
		goto abort;
		}
	switch (head.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		bReduced = FALSE;
		break;
	case SSIO_MAGIC_REDUCED:
		bReduced = TRUE;
		break;
	default:
		fprintf(stderr,"Unrecognized magic number (%i) in \"%s\".\n",head.iMagicNumber,orientfile);
		goto abort;
		}
	if ((*o = (MATRIX *) realloc(*o,(size_t) h->n_data*sizeof(MATRIX))) == NULL) {
		fprintf(stderr,"Unable to allocate memory for orientation matrices.\n");
		goto abort;
		}
	for (i=0;i<h->n_data;i++) {
		if (bReduced) {
			if (ssioOriDataReduced(&ssio,&oridata.orir) != 0) {
				fprintf(stderr,"Problem reading data from \"%s\" (particle %i).\n",orientfile,i);
				goto abort;
				}
			for (k=0;k<N_DIM;k++) {
				(*o)[i][0][k] = oridata.orir.p1[k];
				(*o)[i][1][k] = oridata.orir.p2[k];
				}
			}
		else {
			if (ssioOriData(&ssio,&oridata.ori) != 0) {
				fprintf(stderr,"Problem reading data from \"%s\" (particle %i).\n",orientfile,i);
				goto abort;
				}
			for (k=0;k<N_DIM;k++) {
				(*o)[i][0][k] = oridata.ori.p1[k];
				(*o)[i][1][k] = oridata.ori.p2[k];
				}
			}
		/* construct 3rd vector */
		CROSS((*o)[i][X], (*o)[i][Y], (*o)[i][Z]);
		}
	ssioClose(&ssio);
	return 0;
 abort:
	if (*o != NULL)
		free((void *) *o);
	ssioClose(&ssio);
	return 1;
	}

static void get_wall_data(PARAMS *p,int bVerbose)
{
	FILE *fp = NULL;
	double dTime;

	assert(p != NULL);
	p->nWalls = 0;
	p->pWalls = NULL;
	if (!p->bDoWalls)
		return;
	if (bVerbose)
		printf("Reading wall data from \"%s\"...\n",p->achWallFile);
	if ((fp = fopen(p->achWallFile,"r")) == NULL) {
		fprintf(stderr,"Unable to open \"%s\".\n",p->achWallFile);
		goto abort;
		}
	if (wallsParseWallsFile(fp,&p->nWalls,&p->pWalls,&dTime,bVerbose) != 0) {
		fprintf(stderr,"Error occurred while parsing wall data.\n");
		goto abort;
		}
	fclose(fp);
	if (bVerbose)
		printf("Number of walls read = %i.\n",p->nWalls);
	assert(p->pWalls != NULL);
	if (p->dWallTimeOffset == 0.0)
		p->dWallTimeOffset = dTime;
	return;
 abort:
	if (fp != NULL)
		fclose(fp);
	if (p->pWalls != NULL)
		free((void *) p->pWalls);
	exit(1);
	}

static int get_params(const char *parfile,PARAMS *p)
{
	int di;

	OpenPar(parfile);
	/* Units need to be read first */
	ReadDbl("Length unit",&p->dLenUnit);
	if (p->dLenUnit <= 0.0) {
		fprintf(stderr,"ERROR: Length unit must be positive.\n");
		return 1;
		}
	ReadDbl("Time unit",&p->dTimeUnit);
	if (p->dLenUnit <= 0.0) {
		fprintf(stderr,"ERROR: Time unit must be positive.\n");
		return 1;
		}
	/* Drawing parameters */
	ReadInt("Frame size",&p->iFrameSize);
	if (p->iFrameSize <= 0) {
		fprintf(stderr,"ERROR: Frame size must be positive.\n");
		return 1;
		}
	ReadInt("Light source option",&p->Light.iOpt);
	p->Light.iLock = 0; /* so get_com_by_COL() will work if called */
	switch (p->Light.iOpt) {
	case LightFIX:
	case LightLRG:
	case LightPOT:
	case LightCAM:
		break;
	case LightCOL:
	case LightPID:
	case LightOID:
		ReadInt("Light source lock",&p->Light.iLock);
		if (p->Light.iOpt != LightOID && p->Light.iLock < 0) {
			fprintf(stderr,"ERROR: Light source lock must be >=0 for selected light option.\n");
			return 1;
			}
		break;
	default:
		fprintf(stderr,"ERROR: Unrecognized light source option (%i).\n",p->Light.iOpt);
		return 1;
		}
	ReadNDbl("Light source position",p->Light.vOrig,N_DIM);
	ReadDbl("Light source intensity",&p->dLightInten);
	if (p->dLightInten < 0.0) {
		fprintf(stderr,"ERROR: Light source intensity cannot be negative.\n");
		return 1;
		}
	ReadInt("Camera option",&p->Camera.iOpt);
	p->Camera.iLock = 0;
	switch (p->Camera.iOpt) {
	case CamFIX:
	case CamLRG:
	case CamPOT:
		break;
	case CamCOL:
	case CamPID:
	case CamOID:
		ReadInt("Camera lock",&p->Camera.iLock);
		if (p->Camera.iOpt != CamOID && p->Camera.iLock < 0) {
			fprintf(stderr,"ERROR: Camera lock must be >=0 for selected camera option.\n");
			return 1;
			}
		break;
	default:
		fprintf(stderr,"ERROR: Unrecognized camera option (%i).\n",p->Camera.iOpt);
		return 1;
		}
	ReadNDbl("Camera position",p->Camera.vOrig,N_DIM);
	ReadInt("View option",&p->LookAt.iOpt);
	p->LookAt.iLock = 0;
	switch (p->LookAt.iOpt) {
	case ViewFIX:
	case ViewLRG:
	case ViewPOT:
	case ViewDIR:
		break;
	case ViewCOL:
	case ViewPID:
	case ViewOID:
		ReadInt("View lock",&p->LookAt.iLock);
		if (p->LookAt.iOpt != ViewOID && p->LookAt.iLock < 0) {
			fprintf(stderr,"ERROR: View lock must be >=0 for selected view option.\n");
			return 1;
			}
		break;
	default:
		fprintf(stderr,"ERROR: Unrecognized view option (%i).\n",p->LookAt.iOpt);
		return 1;
		}
	ReadNDbl("Look at/view direction",p->LookAt.vOrig,N_DIM);
	if (p->LookAt.iOpt == ViewDIR) {
		if (MAG(p->LookAt.vOrig) == 0.0) {
			fprintf(stderr,"ERROR: View direction must have non-zero magnitude.\n");
			return 1;
			}
		NORM_VEC(p->LookAt.vOrig,MAG(p->LookAt.vOrig)); /* a different approach would be to interpret the magnitude as the focus distance */
		}
	ReadDbl("Focus distance",&p->dFocusDist);
	if (p->dFocusDist < 0.0)
		p->dFocusDist /= -p->dLenUnit;
	if (p->LookAt.iOpt == ViewDIR && p->dFocusDist == 0.0) {
		fprintf(stderr,"ERROR: Focus distance must be non-zero.\n");
		return 1;
		}
	ReadNDbl("Sky vector",p->vSkyVec,N_DIM);
	if (MAG(p->vSkyVec) == 0.0) {
		fprintf(stderr,"ERROR: Sky vector must have non-zero magnitude.\n");
		return 1;
		}
	NORM_VEC(p->vSkyVec,MAG(p->vSkyVec));
	ReadDbl("Starting view size",&p->dViewSize);
	if (p->dViewSize < 0.0)
		p->dViewSize /= -p->dLenUnit;
	if (p->dViewSize == 0.0) {
		printf("Auto-zoom invoked: View direction fixed.\n");
		p->LookAt.iOpt = ViewFIX; /* LookAt.vPos set in get_view_size() */
		}
	ReadDbl("View scaling criterion",&p->dViewScaleCrit);
	if (p->dViewScaleCrit < 0.0 || p->dViewScaleCrit > 1.0) {
		fprintf(stderr,"ERROR: View scaling criterion must be between 0 and 1 (inclusive).\n");
		return 1;
		}
	p->dViewScale = 1.0;
	if (p->dViewScaleCrit > 0.0) {
		ReadDbl("View scaling factor",&p->dViewScale);
		if (p->dViewScale <= 1.0) {
			fprintf(stderr,"ERROR: View scaling factor must be greater than 1.\n");
			return 1;
			}
		}
	ReadInt("Scale length by radius?",&di);
	p->bScaleLengthByRadius = di;
	ReadDbl("View target zoom factor",&p->dViewTargetZoomFactor);
	if (p->dViewTargetZoomFactor < 0.0) {
		fprintf(stderr,"ERROR: View target zoom factor cannot be negative.\n");
		return 1;
		}
	if (p->dViewTargetZoomFactor > 0.0 && p->LookAt.iOpt != ViewLRG &&
		p->LookAt.iOpt != ViewPID && p->LookAt.iOpt != ViewOID) {
		printf("WARNING: View target zoom factor ignored.\n");
		p->dViewTargetZoomFactor = 0.0;
		}
	ReadInt("Particle shape",&di);
	p->Shape = (SHAPE) di;
	switch (p->Shape) {
	case Point:
	case SolidSphere:
	case POV:
		break;
	default:
		fprintf(stderr,"ERROR: Unrecognized particle shape (%i).\n",p->Shape);
		return 1;
		}
	ReadDbl("Particle radius scaling",&p->dRScale);
	if (p->Shape != Point && p->dRScale <= 0.0) {
		fprintf(stderr,"ERROR: Particle radius scaling must be positive.\n");
		return 1;
		}
	ReadInt("Color override",&p->iNewColor);
	if (p->iNewColor < 0 || p->iNewColor >= NUM_COLORS) {
		fprintf(stderr,"ERROR: Invalid color override value (%i).\n",p->iNewColor);
		return 1;
		}
	p->bZSort = FALSE;
	if (p->Shape != POV) {
		ReadInt("Hide blocked objects?",&di);
		p->bZSort = (BOOLEAN) di;
		}
	ReadInt("Draw color first",&p->iFirstColor);
	if (p->iFirstColor < 0 || p->iFirstColor >= NUM_COLORS) {
		fprintf(stderr,"ERROR: Invalid first color value (%i).\n",p->iFirstColor);
		return 1;
		}
	/*
	** View motion controls...
	*/
	ReadInt("Start at COM?",&di);
	p->bStartAtCOM = (BOOLEAN) di;
	ReadDbl("Inertia control",&p->dInertia);
	if (p->dInertia > 1.0) {
		fprintf(stderr,"ERROR: Inertia control must be <= 1.\n");
		return 1;
		}
	if (p->Light.iOpt == LightLRG || p->Camera.iOpt == CamLRG || p->LookAt.iOpt == ViewLRG) {
		ReadDbl("Min. target mass frac.",&p->dMinTargetMassFraction);
		if (p->dMinTargetMassFraction < 0.0 || p->dMinTargetMassFraction > 1.0) {
			fprintf(stderr,"ERROR: Min. target mass frac. must be between 0 & 1 (inclusive).\n");
			return 1;
			}
		}
	/*
	** POV-Ray stuff...
	*/
	if (p->Shape == POV) {
		ReadStr("Shape file",p->achShapeFile,MAXPATHLEN);
		assert(strlen(p->achShapeFile) > 0);
		ReadInt("Renormalize?",&di);
		p->bRenorm = (BOOLEAN) di;
		ReadStr("Aspect ratio",p->achAspectRatio,MAXPATHLEN);
		assert(strlen(p->achAspectRatio) > 0);
		ReadDbl("Camera light intensity",&p->dCamInten);
		assert(p->dCamInten >= 0.0);
		ReadDbl("Blob threshold",&p->dBlob);
		assert(p->dBlob >= 0.0);
		ReadInt("Highlight index",&p->iHighlightIndex);
		if (p->iHighlightIndex >= 0)
			p->bZSort = TRUE;
		}
	else if (p->dLightInten > 1.0) {
		p->dLightInten = 1.0;
		fprintf(stderr,"WARNING: light intensity reset to 1 (max).\n");
		}
	/*
	** Sliding patch_lines stuff...
	*/
	ReadInt("Use sliding patches?",&di);
	p->bUsePatch = (BOOLEAN) di;
	if (p->bUsePatch) {
		ReadDbl("Patch width",&p->dLx);
		if (p->dLx < 0.0) p->dLx /= -p->dLenUnit;
		assert(p->dLx > 0.0);
		ReadDbl("Patch length",&p->dLy);
		if (p->dLy < 0.0) p->dLy /= -p->dLenUnit;
		assert(p->dLy > 0.0);
		ReadDbl("Orbital frequency",&p->dOmega);
		if (p->dOmega < 0.0) p->dOmega *= -p->dTimeUnit;
		ReadInt("Number of replicas",&p->nRepl);
		assert(p->nRepl > 0);
		ReadInt("Ghost color",&p->iGhostColor);
		ReadInt("Draw patch lines?",&di);
		p->bDrawPatchLines = (BOOLEAN) di;
		}
	ReadStr("Wall data file",p->achWallFile,MAXPATHLEN);
	p->bDoWalls = (strlen(p->achWallFile) > 0);
   	if (p->bDoWalls)
		ReadDbl("Wall time offset",&p->dWallTimeOffset);
	/*
	** Minimum gravitational potential stuff...
	*/
	if (p->Light.iOpt == LightPOT || p->Camera.iOpt == CamPOT || p->LookAt.iOpt == ViewPOT) {
		ReadInt("Potentials to sample",&p->nSample);
		/* (check for valid negative value in main() */
		if (p->nSample == 0 || p->nSample > 100) {
			fprintf(stderr,"ERROR: %% potentials to sample must be > 0 and <= 100.\n");
			return 1;
			}
		ReadInt("Potential sample order",&p->iSampleOrder);
		switch (p->iSampleOrder) {
		case PotRan:
			if (p->nSample == 100) {
				fprintf(stderr,"ERROR: %% potentials must be < 100 for this sample order.\n");
				return 1;
				}
		case PotOrd:
		case PotLrg:
			break;
		default:
			fprintf(stderr,"ERROR: Unrecognized potential sample order option.\n");
			return 1;
			}
		ReadDbl("Tree opening angle",&p->dTheta);
		if (p->dTheta < 0.0) {
			fprintf(stderr,"ERROR: Tree opening angle cannot be negative.\n");
			return 1;
			}
		if (p->dTheta >= 1.0)
			fprintf(stderr,"WARNING: Large tree opening angle (%g rad).\n",p->dTheta);
		}
	ClosePar();
	/* initialize inertia control */
	COPY_VEC(p->Light.vOrig,p->Light.vPos);
	COPY_VEC(p->Camera.vOrig,p->Camera.vPos);
	COPY_VEC(p->LookAt.vOrig,p->LookAt.vPos);
	return 0;
	}

int main(int argc,char *argv[])
{
	PARAMS params;
	COLORMAP_T colormap;
	IMAGE_T *image;
	SSHEAD head;
	SSDATA *data=NULL;
	MATRIX *orient=NULL;
	double dMinMass;
	const char *parfile=PARFILE_DFLT;
	char *scenefile=NULL,outfile[MAXPATHLEN];
	BOOLEAN usage=FALSE,bNeedPotPos=FALSE,auto_zoom=FALSE;
	int c,i,j,nSample,nMoveFrames;

	params.b24bit = params.bDoOrient = params.bWallsOnly = params.bNoGZip = FALSE;

	while (!usage && (c = getopt(argc,argv,"cop:s:wz")) != EOF)
		switch (c) {
		case 'c':
			params.b24bit = TRUE;
			break;
		case 'o':
			params.bDoOrient = TRUE;
			break;
		case 'p':
			parfile = optarg;
			break;
		case 's':
			scenefile = optarg;
			break;
		case 'w':
			params.bWallsOnly = TRUE;
			break;
		case 'z':
			params.bNoGZip = TRUE;
			break;
		default:
			usage = TRUE;
			}

	if (optind == argc && !params.bWallsOnly)
		usage = TRUE;

	if (usage) {
		fprintf(stderr,"Usage: %s [ -c ] [ -o ] [ -p par-file ] [ -s scene-file ] [ -w ] [ -z ] ss-file [ ss-file ... ]\n",argv[0]);
		fprintf(stderr,"Use -c if colors in ss file(s) are in 24-bit format (POV-Ray only).\n");
		fprintf(stderr,"Use -o to draw sphere orientations (POV-Ray only, needs .ori files).\n");
		fprintf(stderr,"Use -w to only draw walls (no particles) -- ss-file(s) optional.\n");
		fprintf(stderr,"Use -z to suppress gzipping of pov files.\n");
		fprintf(stderr,"NOTE: scene-file only used for autozooming and LRG option.\n");
		return 1;
		}

	setbuf(stdout,(char *)NULL);
	srandom((unsigned int) getpid());

	if (get_params(parfile,&params) != 0)
		return 1;

	if (params.bDoOrient && params.Shape != POV) {
		fprintf(stderr,"ERROR: particle shape must be %i to draw sphere orientations.\n",POV);
		return 1;
	}

	if (params.bWallsOnly && !params.bDoWalls) {
		fprintf(stderr,"WARNING: -w option ignored (no walls file specified).\n");
		params.bWallsOnly = FALSE;
		}

	if (params.Light.iOpt == LightPOT || params.Camera.iOpt == CamPOT || params.LookAt.iOpt == ViewPOT)
		bNeedPotPos = TRUE;
	if (params.dViewSize == 0.0)
		auto_zoom = TRUE;
	if (params.bWallsOnly && auto_zoom && scenefile == NULL && optind == argc) {
		fprintf(stderr,"Need scene file or ss file for walls-only auto zooming.\n");
		return 1;
		}
	if (params.Shape != POV)
		MakeColormap(&colormap);
	get_wall_data(&params,WALLS_VERBOSE);
	if (params.dInertia < 0.0)
		nMoveFrames = - (int) params.dInertia;
	else
		nMoveFrames = (int) (params.dInertia*(argc - optind));

	/* if auto-zooming, and a scenefile is specified, fix the view size from it */

	if (auto_zoom && scenefile != NULL) {
		if (read_data(scenefile,&head,&data) != 0)
			return 1;
		printf("Obtaining view info from scene file...\n");
		get_view_size(&params,head.n_data,data);
		free((void *) data);
		}

	/* if LRG option is set, get largest mass either from scenefile, or last frame */

	dMinMass = 0.0;
	if (params.Light.iOpt == LightLRG || params.Camera.iOpt == CamLRG || params.LookAt.iOpt == ViewLRG) {
		int rv;
		if (params.bWallsOnly && scenefile == NULL && optind == argc) {
			fprintf(stderr,"Need scene file or ss file to find largest mass.\n");
			return 1;
			}
		printf("Obtaining info for largest mass ");
		if (scenefile != NULL) {
			printf("from scene file %s...\n",scenefile);
			rv = read_data(scenefile,&head,&data);
			}
		else {
			printf("from last frame %s...\n",argv[argc - 1]);
			rv = read_data(argv[argc - 1],&head,&data);
			}
		if (rv == 0) {
			VECTOR vDum;
			double dDum;
			int iDum;
			get_lrg(data,head.n_data,&iDum,&dMinMass,&dDum,vDum);
			free((void *) data);
			dMinMass = params.dMinTargetMassFraction*dMinMass;
			}
		else
			return 1;
		}

	/* final parameter checking for walls-only option */

	if (params.bWallsOnly && optind == argc) {
		const PARAMS *p = &params;
		if (p->Light.iOpt == LightPOT || p->Light.iOpt == LightCOL || p->Light.iOpt == LightPID || p->Light.iOpt == LightOID || p->Camera.iOpt == CamPOT || p->Camera.iOpt == CamCOL || p->Camera.iOpt == CamPID || p->Camera.iOpt == CamOID || p->LookAt.iOpt == ViewPOT || p->LookAt.iOpt == ViewCOL || p->LookAt.iOpt == ViewPID || p->LookAt.iOpt == ViewOID) {
			fprintf(stderr,"Invalid light/camera/view option for walls-only drawing.\n");
			return 1;
			}
		}

	/* special case: just drawing walls, with no ss files */

	if (params.bWallsOnly && optind == argc) {
		/* no ss files */
		printf("Drawing wall file...\n");
		params.dTime = params.dWallTimeOffset;
		params.dOldViewSize = params.dViewSize;
		if (ssioNewExt("walls",SS_EXT,outfile,
					   params.Shape == POV ? POV_EXT : DRAW_EXT)) {
			fprintf(stderr,"Unable to generate output filename for walls.\n");
			goto done;
			}
		if (params.Shape == POV)
			pov_draw(&params,NULL,NULL,NULL,outfile,0,0.0,0);
		else {
			AllocImage(&image,params.iFrameSize,params.iFrameSize);
			draw(&params,NULL,NULL,image,0,0.0,0);
			DumpRaster(&colormap,image,outfile);
			FreeImage(image);
			}
		printf("Done!\n");
		goto done;
		}

	/* loop through frames */

	for (i=optind;i<argc;i++) {
		printf("%s: ",argv[i]);
		if (read_data(argv[i],&head,&data) != 0)
			return 1; /* was "continue" but that can break scripts */
		if (params.bDoOrient) {
			if (read_orient(argv[i],&head,&orient) != 0)
				return 1;
			if (orient == NULL) {
				/* if undefined, set to identity matrices (a bit of a hack) */
				assert(head.n_data > 0);
				orient = (MATRIX *) malloc((size_t) head.n_data*sizeof(MATRIX));
				assert(orient != NULL);
				for (j=0;j<head.n_data;j++) {
					UNIT_MAT(orient[j]);
					}
				}
			}
		params.dTime = head.time; /* currently needed only for moving walls */
		if (params.bDoWalls) {
			if (params.bRenorm && i > optind)
				get_wall_data(&params,WALLS_SILENT); /* renormalization overwrites wall data each time */
			params.dTime += params.dWallTimeOffset;
			}
		assert(data != NULL);
		if (auto_zoom && scenefile == NULL)
			get_view_size(&params,head.n_data,data);
		if (i == optind) /* first frame */
			params.dOldViewSize = params.dViewSize;
		if (bNeedPotPos) {
			/* needed to know number of particles before checking these params */
			assert(params.nSample != 0);
			if (params.nSample < 0)
				nSample = -params.nSample;
			else {
				nSample = (int) (0.01*params.nSample*head.n_data);
				if (nSample == 0) {
					fprintf(stderr,"WARNING: Num. potentials to sample too small -- forced to 1.\n");
					nSample = 1;
					}
				}
			if (nSample > head.n_data) {
				fprintf(stderr,"ERROR: Num. potentials to sample (%i) > num. particles.\n",nSample);
				return 1;
				}
			if (nSample == head.n_data && params.iSampleOrder == PotRan) {
				fprintf(stderr,"ERROR: Num. potentials must be < num. particles for this sample order.\n");
				return 1;
				}
			get_pot_pos(&params,data,head.n_data,nSample);
			}
		if (params.iNewColor)
			for (j=0;j<head.n_data;j++)
				data[j].color = params.iNewColor;
		if (ssioNewExt(argv[i],SS_EXT,outfile,
					   params.Shape == POV ? POV_EXT : DRAW_EXT)) {
			fprintf(stderr,"Unable to generate output filename for %s.\n",argv[i]);
			free((void *) data);
			continue;
			}
		if (params.Shape == POV)
			pov_draw(&params,&head,data,orient,outfile,nMoveFrames,dMinMass,params.iFirstColor);
		else {
			AllocImage(&image,params.iFrameSize,params.iFrameSize);
			draw(&params,&head,data,image,nMoveFrames,dMinMass,params.iFirstColor);
			if (params.iFirstColor > 0)
				draw(&params,&head,data,image,nMoveFrames,dMinMass,-params.iFirstColor); /* draws everything except first color */
			DumpRaster(&colormap,image,outfile);
			FreeImage(image);
			}
		free((void *) data);
		printf("Drawing done.\n");
		}

 done:

	if (params.pWalls != NULL)
		free((void *) params.pWalls);

	if (orient != NULL)
		free((void *) orient);

	return 0;
	}

/* ssdraw.c */
