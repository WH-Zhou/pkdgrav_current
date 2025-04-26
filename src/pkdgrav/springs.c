/* suppress pedantic "empty translation unit" warning if !SPRINGS */
typedef int make_iso_compilers_happy;

#ifdef SPRINGS

#include <math.h>

int sign(double A);
#define sign(A) ((A)<0.?-1:(A)>0.?1:0)
#include "smoothfcn.h" /* includes pkd.h, which includes springs.h */

#define iParticles 147000 /*DEBUG!*/

#ifdef SPRINGS_STRESS_MOVE
/*
void springsCalcPressure(PARTICLE *p,int iParticles,float fZeroStrainLength,double dCubeStrain,unsigned char start_of_step,double dConvert)
{
	static double x_force,y_force,z_force,dPressure;  // this will be the total restoring force in the -x(-y,-z) direction on the plane of blue particles pulling in the +x(+y,+z) direction (so 2x this quantity for total force)
        if (start_of_step) {
		x_force = 0.;
		y_force = 0.;
		z_force = 0.;
	}

        // iTugLevel = iTugLevelStart + iTugCounter; alternate way of calculating area (would need to pass variables into this function)
	if (p->iColor == 4 && p->r[0] > 0) { // x
	        x_force += p->fMass*p->a[0]; // this is summing ACTUAL force, not acc
	}
	if ((p->iColor == 6 || p->iColor == 4) && p->r[1] > 0) { // y 
		y_force += p->fMass*p->a[1];
	}
	if ((p->iColor == 1 || p->iColor == 6 || p->iColor == 4) && p->r[2] > 0) { // z
		z_force += p->fMass*p->a[2];
	}
	if (!(start_of_step)) {
	        // NOTE: assumes equal-sized particles, equal initial seperation!
	        int iParticles_on_side = pow(iParticles,(1./3.)) + 1.e-8;
		double dTotalArea = (fZeroStrainLength*dCubeStrain*(iParticles_on_side-1)+(1.98*p->fSoft)); // assumes fZeroStrainLength and dstrain are the same for all springs
		//double dTotalArea = (dNeighbor_dx[p->iOrder]*(1.+dStretchStep*iTugLevel))*(iParticles_on_side-1)+1.98*RADIUS(p);
		//printf("%e %e\n",dNeighbor_dx[p->iOrder],dTotalArea);
		dTotalArea *= dTotalArea;
		//printf("%e\n",dTotalArea);
		dTotalArea *= 6.;
		//printf("%e\n",dTotalArea);
		dPressure = -(x_force+x_force+y_force+y_force+z_force+z_force)/(dTotalArea*dConvert);
		//printf("Step = %li; Fx; = %e; Fy = %e; Fz = %e; side = %d; Area = %e; Total Pressure estimate = %e\n",iStep,x_force,y_force,z_force,iParticles_on_side,dTotalArea,dPressure);
	}
}
*/

void springsMoveSpecials(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
        /*
	** MOVE_TENSILE - move outermost x-y(y-z,z-x) planes in x(y,z)-dimension away from x(y,z)=0 plane (tensile stress)
	** Special Particle Colors(3): 4,6,1
        */

	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	double dMove = 5e-11; /* *-1.*sin(iStep*2.*M_PI*0.0001); hack for oscillation */
	/* if ((0.0001 * iStep) >= 1) dMove = 0; */ /* hack for oscillation */
	static const int iTugFreq = iParticles * 1; /* increase distance of special particle by dStretchStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	//	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
		//      start_of_step = 1;
		if (iStep == 0) printf("SPRINGS_STRESS_TEST %d\n",SPRINGS_STRESS_TEST);
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

	static const double dStretchStep = 0.005; /* increase in strain per (iTugFreq/iParticles) steps */

	/*
	** the use of the following arrays (size of iParticles - could
	** optimize by making array only as large as number of special
	** particles) allow the special particles to keep the specified
	** strain with their neighbors from step to step.  This allows
	** the entire object to approach the specified strain rapidly.
	*/

	static double dPosition_x[iParticles];                                                    /* initial position in x-dimension */
	static double dPosition_y[iParticles];                                                    /* initial position in y-dimension */
	static double dPosition_z[iParticles];                                                    /* initial position in z-dimension */
	static int iNeighbors_colorx[iParticles];                                                 /* color assigned to x-neighbor */
	static int iNeighbors_colory[iParticles];                                                 /* color assigned to y-neighbor */
	static int iNeighbors_colorz[iParticles];                                                 /* color assigned to z-neighbor */
	static int iNeighbors_number[iParticles];                                                 /* number of attached neighbors (i.e. number of attached springs) */

	/* zero-out neighbor positions */
	if (iStep == 1)
	        iNeighbors_colorx[p->iOrder] = iNeighbors_colory[p->iOrder] = iNeighbors_colorz[p->iOrder] = iNeighbors_number[p->iOrder] = 0.;

	const PARTICLE *pn;
	SPRING *ps,pns;
	int i,j,bFoundNbr;

	if (p->iColor != 4 && p->iColor != 6 && p->iColor != 1)
		return; /* save time by ignoring non-special particles */

	/* loop over this particle's neighbors */
	if (iStep == 1) {
	        dPosition_x[p->iOrder] = p->r[0];
		dPosition_y[p->iOrder] = p->r[1];
		dPosition_z[p->iOrder] = p->r[2];
	        for (i=0;i<nSmooth;i++) {
			pn = nnList[i].pPart;
			if (p == pn) continue;
			if (pn->iOrder < p->iOrder) {
			        for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
				        pns = pn->springs[j];
					if (pns.iOrder == p->iOrder && pns.fZeroStrainLength != 0.0) {
					        iNeighbors_number[p->iOrder]++;
						if (pn->iColor == 4) iNeighbors_colorx[p->iOrder]++;
						if (pn->iColor == 6) iNeighbors_colory[p->iOrder]++;
						if (pn->iColor == 1) iNeighbors_colorz[p->iOrder]++;
					}
				}
			}
			else {
			        for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
				        ps = &p->springs[j];
					if (ps->iOrder == pn->iOrder && ps->fZeroStrainLength != 0.0) {
					        iNeighbors_number[p->iOrder]++;
						if (pn->iColor == 4) iNeighbors_colorx[p->iOrder]++;
						if (pn->iColor == 6) iNeighbors_colory[p->iOrder]++;
						if (pn->iColor == 1) iNeighbors_colorz[p->iOrder]++;
					}
				}
			}
		}
	}

#if (SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE)
        iTugLevel = iTugLevelStart + iTugCounter; /* to hack for oscillation, comment out "+ iTugCounter" */
	if (p->iColor == 4) { /* x */
	        p->r[0] = dPosition_x[p->iOrder] + sign(p->r[0])*dMove*(dStretchStep*iTugLevel);
	        p->r[1] = dPosition_y[p->iOrder]; /* hack: need to fix for more than one special color */
	        p->r[2] = dPosition_z[p->iOrder]; /* hack: need to fix for more than one special color */
		p->a[0] = 0.;         
		p->v[0] = 0.;
		/* if (p->r[0] < 0) p->r[0] = dPosition_x[p->iOrder]; include to hack for oscillation */
	}
	if (p->iColor == 6) { /* y */
	        p->r[1] = dPosition_y[p->iOrder] + sign(p->r[1])*dMove*(dStretchStep*iTugLevel);
		p->a[1] = 0.;
		p->v[1] = 0.;
	}
	if (p->iColor == 10) { /* DEBUG:HACK!! */
	        p->a[0] = 0.0005*(p->r[0] - dPosition_x[p->iOrder]);
		p->a[1] = 0.0005*(p->r[1] - dPosition_y[p->iOrder]);
	}
	if (p->iColor == 1) { /* z */
	        p->r[2] = dPosition_z[p->iOrder] + sign(p->r[2])*dMove*(dStretchStep*iTugLevel);
		p->a[2] = 0.;
		p->v[2] = 0.;
	}
#endif /* !(SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_TENSILE) */

#if (SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_SHEAR)
        iTugLevel = iTugLevelStart + iTugCounter;
	if (p->iColor == 4) { /* x */
	        p->r[0] = dPosition_x[p->iOrder];
	        p->r[1] = dPosition_y[p->iOrder];
	        p->r[2] = dPosition_z[p->iOrder] - sign(p->r[0])*dMove*(dStretchStep*iTugLevel);
		p->a[0] = 0.;         
		p->v[0] = 0.;
		p->a[1] = 0.;         
		p->v[1] = 0.;
		p->a[2] = 0.;         
		p->v[2] = 0.;
	}
#endif /* !(SPRINGS_STRESS_TEST == SPRINGS_STRESS_MOVE_SHEAR) */

}

#endif /* SPRINGS_STRESS_MOVE */

void springsForceSpecialsTensile(PARTICLE *p)
{
        /*
	** SPRINGS_STRESS_FORCE_TENSILE - force outermost x-y(y-z,z-x) planes in x(y,z)-dimension away from x(y,z)=0 plane (tensile stress)
	** Special Particle Colors(3): 4,6,1
        */

	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	static const int iTugFreq = iParticles * 2500; /* increase acceleration of special particle by dTugStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	//	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
		//      start_of_step = 1;
		if (iStep == 0) printf("SPRINGS_STRESS_TEST %d\n",SPRINGS_STRESS_TEST);
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

	static const double dTugStep = 0.005; /* increase in force per (iTugFreq/iParticles) steps */

        iTugLevel = iTugLevelStart + iTugCounter;
	if (p->iColor == 4) { /* x */
		/* if (p->r[0] > 0) x_force += p->fMass*p->a[0]; */ /* this sums force, not acc */
	        p->a[0] += sign(p->r[0])*dTugStep*iTugLevel; /* for pulling face off */
#ifdef YORP
#include <time.h>
#include <unistd.h>
	        double dTheta,dAlpha,dForceMag;
		unsigned int uSeed;
		static double dRadiative[iParticles][3];
		if (iCounter == 1) {
		        uSeed = time(NULL)%getpid() + getppid();
			(void) printf("YORP random seed = %u\n",uSeed);
			srandom(uSeed);
			for (i=0;i<iParticles;i++) {
				double dTheta = M_PI * (randUniform() - 0.5); /* a double between -pi/2 & pi/2 */
				double dAlpha = M_PI * (randUniform() - 0.5); /* a double between 0 & pi/2 */
				double dForceMag = 1e-3 * randUniform();
				dRadiative[i][0] = -dForceMag * cos(dTheta) * cos(dAlpha); /* x */
				dRadiative[i][1] = dForceMag * sin(dTheta) * cos(dAlpha); /* y */
				dRadiative[i][2] = dForceMag * cos(dTheta) * sin(dAlpha); /* z */
				printf("dTheta = %f  dAlpha = %f  dForceMag = %e  dRadiative[i][0] = %f  dRadiative[i][1] = %f  dRadiative[i][2] = %f\n",dTheta,dAlpha,dForceMag,dRadiative[i][0],dRadiative[i][1],dRadiative[i][2]);
			}
		}
	}
#else /* YORP */
        }
	if (p->iColor == 6) p->a[1] += sign(p->r[1])*dTugStep*iTugLevel; /* y */
	if (p->iColor == 1) p->a[2] += sign(p->r[2])*dTugStep*iTugLevel; /* z */
#endif /* !YORP */
}

void springsForceSpecialsShear(PARTICLE *p)
{
        /*
	** SPRINGS_STRESS_FORCE_SHEAR - force outermost x-y planes into z-dimension (opposite directions)
	** Special Particle Colors(1): 4
        */

	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	static const int iTugFreq = iParticles * 2500; /* increase acceleration of special particle by dTugStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	//	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
		//      start_of_step = 1;
		if (iStep == 0) printf("SPRINGS_STRESS_TEST %d\n",SPRINGS_STRESS_TEST);
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

	static const double dTugStep = 0.01; /* increase in force per (iTugFreq/iParticles) steps */

        iTugLevel = iTugLevelStart + iTugCounter;
	if (p->iColor == 4) { /* x */

		static double d_pos_x[iParticles]; /* initial position (could optimize by only making arrays as large as number of special particles) */
		static double d_pos_y[iParticles];
		if (iStep == 1) {
		        d_pos_x[p->iOrder] = p->r[0];
			d_pos_y[p->iOrder] = p->r[1];
		}
		p->r[0] = d_pos_x[p->iOrder];
		p->r[1] = d_pos_y[p->iOrder];
		p->v[0] = 0.;
		p->v[1] = 0.;
		p->a[0] = 0.;
		p->a[1] = 0.;
		p->a[2] -= sign(p->r[0])*dTugStep*iTugLevel;
	}
}

void springsForceSpecialsTwist(PARTICLE *p)
{
        /*
	** SPRINGS_STRESS_FORCE_TWIST - tangential force applied to outermost x-y planes around x-axis
	** Special Particle Colors(1): 4
        */

	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	static const int iTugFreq = iParticles * 2500; /* increase acceleration of special particle by dTugStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	//	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
		//      start_of_step = 1;
		if (iStep == 0) printf("SPRINGS_STRESS_TEST %d\n",SPRINGS_STRESS_TEST);
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

	static const double dTugStep = 0.01; /* increase in force per (iTugFreq/iParticles) steps */

        iTugLevel = iTugLevelStart + iTugCounter;
	if (p->iColor == 4) { /* x */

		static double d_pos_x[iParticles]; /* initial position (could optimize by only making arrays as large as number of special particles) */
		static double d_pos_y[iParticles];
		if (iStep == 1) {
		        d_pos_x[p->iOrder] = p->r[0];
			d_pos_y[p->iOrder] = p->r[1];
		}
		static double d_pos_z[iParticles];
		if (iStep == 1) d_pos_z[p->iOrder] = p->r[2];
		if (d_pos_y[p->iOrder] == 0. && d_pos_z[p->iOrder] == 0.) {
		        p->r[1] = d_pos_y[p->iOrder];
			p->r[2] = d_pos_z[p->iOrder];
		}
		p->r[0] = d_pos_x[p->iOrder];
		p->a[0] = 0.;
		p->a[1] += sign(p->r[0])*dTugStep*iTugLevel*1.e9*p->r[2];
		p->a[2] -= sign(p->r[0])*dTugStep*iTugLevel*1.e9*p->r[1];		
	}
}

void springsForceSpecialsRadial(PARTICLE *p)
{
        /*
	** SPRINGS_STRESS_FORCE_RADIAL - applied force directed out radially from origin
	** Special Particle Colors(1): 4
        */

	static long int iCounter = 0;
	static long int iTugCounter = 0;
	static long int iStep = 0;
	++iCounter;

	static const int iTugFreq = iParticles * 2500; /* increase acceleration of special particle by dTugStep after every iTugFreq counts (tracked by iCounter) */
	static const int iTugLevelStart = 1;
	int iTugLevel;
	//	unsigned char start_of_step = 0;
	if (iCounter % iParticles == 1) {
		//      start_of_step = 1;
		if (iStep == 0) printf("SPRINGS_STRESS_TEST %d\n",SPRINGS_STRESS_TEST);
		iStep ++;
	}
	if (iCounter > iTugFreq) {
		iCounter -= iTugFreq;
		++iTugCounter;
	}

	static const double dTugStep = 0.01; /* increase in force per (iTugFreq/iParticles) steps */

        iTugLevel = iTugLevelStart + iTugCounter;
	if (p->iColor == 4) { /* x */

		double dist = sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]+p->r[2]*p->r[2]);
	        p->a[0] += p->r[0]/dist*dTugStep*iTugLevel;
	        p->a[1] += p->r[1]/dist*dTugStep*iTugLevel;
	        p->a[2] += p->r[2]/dist*dTugStep*iTugLevel;
	}
}

void springsAssignColor(PARTICLE *p,int nSmooth,NN *nnList, int ALLOW_COLOR_CHANGE)
{
        SPRING *ps,pns;
	const PARTICLE *pn;
	int i,j,nLinks=0;

#ifndef REFORM_SPRINGS
	if (p->iColor == 2 || !ALLOW_COLOR_CHANGE) /* skip if unlinked or if color changes are not allowed */
#else /* REFORM_SPRINGS */
	if (!ALLOW_COLOR_CHANGE)
#endif /* REFORM_SPRINGS */
	        return;

	/* check for springs attached to particle located in particle's springs array (if they have a higher iOrder number) */
        for (i=0;i<MAX_NUM_SPRINGS_PER_PARTICLE;i++) {
	        ps = &p->springs[i];
		if (ps->iOrder >= 0 && ps->iOrder > p->iOrder && ps->fZeroStrainLength != 0.0)
		        ++nLinks;
	}

	/* check for springs attached to particle located in each of particle's nSmooth neighbors' springs arrays (if they have a lower iOrder number) */
	for (i=0;i<nSmooth;i++) {
	        pn = nnList[i].pPart;
		if (pn->iOrder < p->iOrder)
		        for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
			        pns = pn->springs[j];
				if (pns.iOrder == p->iOrder && pns.fZeroStrainLength != 0.0)
				        ++nLinks;
			}
	}
	/* if (nLinks == 0 && p->iColor != 12) */ /* unattached violet particles will remain violet */

#ifdef SPRINGS_STRESS_MOVE
	if (p->iColor == 4 || p->iColor == 1 || p->iColor == 6)
	        return; /* leave color unchanged */
#endif /* SPRINGS_STRESS_MOVE_TENSILE */
	if (p->iColor < 0)
	        return; /* leave color unchanged */
	else if (nLinks == 0 && p->iColor != 12)
	        p->iColor = 2; /* red: visual indication of unlinked particle */
	else if (nLinks == 1)
	        p->iColor = 10; /* orange */
	else if (nLinks == 2)
	        p->iColor = 5; /* yellow */
#ifdef SPRINGS_COLOR_GREYSCALE
	else {
	        p->iColor = 30 + 15*nLinks;
		if (p->iColor > 255) p->iColor = 255;
	}
#else /* SPRINGS_COLOR_GREYSCALE */ 
	else
		p->iColor = p->iColor;
	        /* p->iColor = 3; */
#endif /* SPRINGS_COLOR_GREYSCALE */ 
}

int pkdBreakSprings(PKD pkd,int iOrder1,int iOrder2)
{
	PARTICLE *p;
	int i,j,nLocal = pkdLocal(pkd),bSpringsChanged = 0;

	for (i=0;i<nLocal;i++) {
		p = &pkd->pStore[i];
		/*
		** Either remove all springs from this particle (if it's one
		** of the designated particles), or remove any springs
		** connected to designated particles.
		 */
		for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++)
			if (p->iOrder == iOrder1 || p->iOrder == iOrder2 ||
				p->springs[j].iOrder == iOrder1 ||
				p->springs[j].iOrder == iOrder2) {
				p->springs[j].fZeroStrainLength = 0.0; /* remove spring */
				bSpringsChanged = 1;
			}
		}
	return bSpringsChanged;
}

#endif /* SPRINGS */
