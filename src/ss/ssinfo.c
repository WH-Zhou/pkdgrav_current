/*
** ssinfo.c -- DCR 3/21/07
** ========
** Outputs simple information about the contents of ss files.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <vector.h>
#include <boolean.h>

typedef struct {
  int bReduced;
  double dTime;
  int nData;
  double dTotMass,*dMass,*dRadius,*dDensity;
  VECTOR *vPos,*vVel,*vSpin;
  int *iColor,*iOrgIdx;
  } AllData;

typedef struct {
  double dMin,dMax;
  int iMinIdx,iMaxIdx,iMinOrg,iMaxOrg;
  } sLimits;

typedef struct {
  VECTOR vMin,vMax,vCom;
  sLimits mag,mag_com;
  } vLimits;

typedef struct {
  int iMin,iMax,iMinIdx,iMaxIdx,iMinOrg,iMaxOrg,nMin,nMax;
  } iLimits;

static int read_ssfile(const char *achFile,AllData *a)
{
	SSIO ssio;
	SSHEAD h;
	union {SSDATA std; SSRDATA red;} data;
	int i;

	if (ssioOpen(achFile,&ssio,SSIO_READ)) {
		fprintf(stderr,"Unable to open ss file for reading.\n");
		return 1;
		}

	if (ssioHead(&ssio,&h) || h.n_data <= 0) {
		fprintf(stderr,"Corrupt header.\n");
		ssioClose(&ssio);
		return 1;
		}

	switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
	case SSIO_MAGIC_REDUCED:
		break;
	default:
		fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}

	a->bReduced = (h.iMagicNumber == SSIO_MAGIC_REDUCED);
	a->dTime = h.time;
	a->nData = h.n_data;

	a->dMass = (double *) realloc(a->dMass,(size_t) a->nData*sizeof(double));
	assert(a->dMass != NULL);
	a->dRadius = (double *) realloc(a->dRadius,(size_t) a->nData*sizeof(double));
	assert(a->dRadius != NULL);
	a->dDensity = (double *) realloc(a->dDensity,(size_t) a->nData*sizeof(double));
	assert(a->dDensity != NULL);
	a->vPos = (VECTOR *) realloc(a->vPos,(size_t) a->nData*sizeof(VECTOR));
	assert(a->vPos != NULL);
	if (a->bReduced)
		a->vVel = a->vSpin = NULL;
	else {
		a->vVel = (VECTOR *) realloc(a->vVel,(size_t) a->nData*sizeof(VECTOR));
		assert(a->vVel != NULL);
		a->vSpin = (VECTOR *) realloc(a->vSpin,(size_t) a->nData*sizeof(VECTOR));
		assert(a->vSpin != NULL);
		}
	a->iColor = (int *) realloc(a->iColor,(size_t) a->nData*sizeof(int));
	assert(a->iColor != NULL);
	a->iOrgIdx = (int *) realloc(a->iOrgIdx,(size_t) a->nData*sizeof(int));
	assert(a->iOrgIdx != NULL);

	a->dTotMass = 0.0;

	for (i=0;i<a->nData;i++) {
		if (a->bReduced) {
			if (ssioDataReduced(&ssio,&data.red)) {
				fprintf(stderr,"Corrupt data (particle %i).\n",i);
				ssioClose(&ssio);
				return 1;
				}
			a->dTotMass += (a->dMass[i] = data.red.fMass);
			a->dRadius[i] = data.red.fRadius;
			COPY_VEC(data.red.vPos,a->vPos[i]);
			a->iColor[i] = data.red.iColor;
			a->iOrgIdx[i] = data.red.iOrgIdx;
			}
		else {
			if (ssioData(&ssio,&data.std)) {
				fprintf(stderr,"Corrupt data (particle %i).\n",i);
				ssioClose(&ssio);
				return 1;
				}
			a->dTotMass += (a->dMass[i] = data.std.mass);
			a->dRadius[i] = data.std.radius;
			COPY_VEC(data.std.pos,a->vPos[i]);
			COPY_VEC(data.std.vel,a->vVel[i]);
			COPY_VEC(data.std.spin,a->vSpin[i]);
			a->iColor[i] = data.std.color;
			a->iOrgIdx[i] = data.std.org_idx;
			}
		a->dDensity[i] = a->dMass[i] / (4.0/3.0*M_PI*CUBE(a->dRadius[i]));
		}

	return 0;
	}

static void get_scalar_limits(const AllData *a,const double d[],sLimits *lim)
{
	int i;

	lim->dMin = HUGE_VAL;
	lim->dMax = -HUGE_VAL;
	lim->iMinIdx = lim->iMaxIdx = INT_MAX;
	lim->iMinOrg = lim->iMaxOrg = INT_MAX;

	for (i=0;i<a->nData;i++) {
		if (d[i] < lim->dMin) {
			lim->dMin = d[i];
			lim->iMinIdx = i;
			lim->iMinOrg = a->iOrgIdx[i];
			}
		if (d[i] > lim->dMax) {
			lim->dMax = d[i];
			lim->iMaxIdx = i;
			lim->iMaxOrg = a->iOrgIdx[i];
			}
		}
	}

static void get_vector_limits(const AllData *a,VECTOR v[],vLimits *lim,BOOLEAN bDoCom)
{
	VECTOR vTmp;
	double *dMag;
	int i;

	SET_VEC(lim->vMin,HUGE_VAL,HUGE_VAL,HUGE_VAL);
	SET_VEC(lim->vMax,-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
	if (bDoCom) {
		ZERO_VEC(lim->vCom);
		}

	for (i=0;i<a->nData;i++) {
		if (v[i][X] < lim->vMin[X])
			lim->vMin[X] = v[i][X];
		if (v[i][Y] < lim->vMin[Y])
			lim->vMin[Y] = v[i][Y];
		if (v[i][Z] < lim->vMin[Z])
			lim->vMin[Z] = v[i][Z];
		if (v[i][X] > lim->vMax[X])
			lim->vMax[X] = v[i][X];
		if (v[i][Y] > lim->vMax[Y])
			lim->vMax[Y] = v[i][Y];
		if (v[i][Z] > lim->vMax[Z])
			lim->vMax[Z] = v[i][Z];
		if (bDoCom) {
			COPY_VEC(v[i],vTmp);
			SCALE_VEC(vTmp,a->dMass[i]);
			ADD_VEC(lim->vCom,vTmp,lim->vCom);
			}
		}

	if (bDoCom) {
		if (a->dTotMass > 0.0) {
			NORM_VEC(lim->vCom,a->dTotMass);
			}
		else {
			NORM_VEC(lim->vCom,a->nData);
			}
		}

	dMag = (double *) malloc((size_t) a->nData*sizeof(double));
	assert(dMag != NULL);

	for (i=0;i<a->nData;i++) {
		dMag[i] = MAG(v[i]);
		}

	get_scalar_limits(a,dMag,&lim->mag);

	if (bDoCom) {
		for (i=0;i<a->nData;i++) {
			SUB_VEC(v[i],lim->vCom,vTmp);
			dMag[i] = MAG(vTmp);
			}
		get_scalar_limits(a,dMag,&lim->mag_com);
		}

	free((void *) dMag);
	}

static void get_int_limits(const AllData *a,const int iVal[],iLimits *lim)
{
	int i;

	lim->iMin = INT_MAX;
	lim->iMax = INT_MIN;
	lim->iMinIdx = lim->iMaxIdx = INT_MAX;
	lim->iMinOrg = lim->iMaxOrg = INT_MAX;

	for (i=0;i<a->nData;i++) {
		if (iVal[i] < lim->iMin) {
			lim->iMin = iVal[i];
			lim->iMinIdx = i;
			lim->iMinOrg = a->iOrgIdx[i];
			}
		if (iVal[i] > lim->iMax) {
			lim->iMax = iVal[i];
			lim->iMaxIdx = i;
			lim->iMaxOrg = a->iOrgIdx[i];
			}
		}

	lim->nMin = lim->nMax = 0;

	for (i=0;i<a->nData;i++) {
		if (iVal[i] == lim->iMin)
			++lim->nMin;
		if (iVal[i] == lim->iMax)
			++lim->nMax;
		}
	}

//This function is modified by Shoucun 07/11/19
static void show_scalar_limits(const sLimits *lim,const char *achType,double dScale,const char *achUnits)
{
	if (lim->dMin <= lim->dMax){
		printf("%s range: %g (%g %s) [%i,%i] to %g (%g %s) [%i,%i]\n",achType,
			   lim->dMin,lim->dMin*dScale,achUnits,lim->iMinIdx,lim->iMinOrg,
			   lim->dMax,lim->dMax*dScale,achUnits,lim->iMaxIdx,lim->iMaxOrg);
	}
	else{
		printf("%s range: %g (%g %s) [%i,%i] to %g (%g %s) [%i,%i]\n",achType,
			   lim->dMax,lim->dMax*dScale,achUnits,lim->iMaxIdx,lim->iMaxOrg,
			   lim->dMin,lim->dMin*dScale,achUnits,lim->iMinIdx,lim->iMinOrg);
	}
}

static void show_vector_limits(const vLimits *lim,const char *achType,double dScale,const char *achUnits,BOOLEAN bDoCom)
{
	printf("%s X range: %g (%g %s) to %g (%g %s)\n",achType,
		   lim->vMin[X],lim->vMin[X]*dScale,achUnits,
		   lim->vMax[X],lim->vMax[X]*dScale,achUnits);
	printf("%s Y range: %g (%g %s) to %g (%g %s)\n",achType,
		   lim->vMin[Y],lim->vMin[Y]*dScale,achUnits,
		   lim->vMax[Y],lim->vMax[Y]*dScale,achUnits);
	printf("%s Z range: %g (%g %s) to %g (%g %s)\n",achType,
		   lim->vMin[Z],lim->vMin[Z]*dScale,achUnits,
		   lim->vMax[Z],lim->vMax[Z]*dScale,achUnits);

	if (bDoCom)
		printf("%s com: %g,%g,%g (%g,%g,%g %s)\n",achType,
			   lim->vCom[X],lim->vCom[Y],lim->vCom[Z],
			   lim->vCom[X]*dScale,lim->vCom[Y]*dScale,lim->vCom[Z]*dScale,achUnits);
	}

static void show_int_limits(const iLimits *lim,const char *achType,int n)
{
	assert(n > 0);
	printf("%s range: %i [%i,%i] (%g%%) to %i [%i,%i] (%g%%)\n",achType,
		   lim->iMin,lim->iMinIdx,lim->iMinOrg,100.0*lim->nMin/n,
		   lim->iMax,lim->iMaxIdx,lim->iMaxOrg,100.0*lim->nMax/n);
	}				  

int main(int argc,char *argv[])
{
	AllData all;
	sLimits sLim;
	vLimits vLim;
	iLimits iLim;
	int i;

	if (argc < 2) {
		fprintf(stderr,"Usage: %s ss-file [ ss-file ... ]\n",argv[0]);
		fprintf(stderr,"Output best viewed in a wide terminal!\n");
		return 1;
		}

	/* mandatory initializations */
	all.dMass = all.dRadius = all.dDensity = NULL;
	all.vPos = all.vVel = all.vSpin = NULL;
	all.iColor = all.iOrgIdx = NULL;
	
	for (i=1;i<argc;i++) {
		printf("%s:\n",argv[i]);
		if (read_ssfile(argv[i],&all) == 0) {
			if (all.bReduced)
				printf("REDUCED FORMAT: some floats are single precision.\n");
			printf("Time = %g (%g yr)\n",all.dTime,all.dTime*T_SCALE/SID_YR);
			printf("Number of particles = %i\n",all.nData);
			printf("Total mass = %g (%g kg)\n",all.dTotMass,all.dTotMass*M_SCALE);
			get_scalar_limits(&all,all.dMass,&sLim);
			show_scalar_limits(&sLim,"Mass",M_SCALE,"kg");
			get_scalar_limits(&all,all.dRadius,&sLim);
			show_scalar_limits(&sLim,"Radius",0.001*L_SCALE,"km");
			get_scalar_limits(&all,all.dDensity,&sLim);
			show_scalar_limits(&sLim,"Density",0.001*D_SCALE,"g/cc");
			get_vector_limits(&all,all.vPos,&vLim,TRUE);
			show_vector_limits(&vLim,"Position",0.001*L_SCALE,"km",TRUE);
			show_scalar_limits(&vLim.mag,"Pos abs mag",0.001*L_SCALE,"km");
			show_scalar_limits(&vLim.mag_com,"Pos com mag",0.001*L_SCALE,"km");
			if (!all.bReduced) {
				get_vector_limits(&all,all.vVel,&vLim,TRUE);
				show_vector_limits(&vLim,"Velocity",V_SCALE,"m/s",TRUE);
				show_scalar_limits(&vLim.mag,"Vel abs mag",V_SCALE,"m/s");
				show_scalar_limits(&vLim.mag_com,"Vel com mag",V_SCALE,"m/s");
				get_vector_limits(&all,all.vSpin,&vLim,FALSE);
				show_vector_limits(&vLim,"Spin",1.0/T_SCALE,"rad/s",FALSE);
				show_scalar_limits(&vLim.mag,"Spin mag",1.0/T_SCALE,"rad/s");
				//The three lines below are added by Shoucun 07/11/19
				vLim.mag.dMin = 6.283185307179586/vLim.mag.dMin;
				vLim.mag.dMax = 6.283185307179586/vLim.mag.dMax;
				show_scalar_limits(&vLim.mag,"Spin period",T_SCALE/3600,"hours");
				}
			get_int_limits(&all,all.iColor,&iLim);
			show_int_limits(&iLim,"Color",all.nData);
			get_int_limits(&all,all.iOrgIdx,&iLim);
			show_int_limits(&iLim,"Original index",all.nData);
			}
		}

	if (all.iOrgIdx == NULL)
		return 1;
	
	free((void *) all.iOrgIdx);
	free((void *) all.iColor);
	if (!all.bReduced) {
		free((void *) all.vSpin);
		free((void *) all.vVel);
		}
	free((void *) all.vPos);
	free((void *) all.dDensity);
	free((void *) all.dRadius);
	free((void *) all.dMass);

	return 0;
	}

/* ssinfo.c */
