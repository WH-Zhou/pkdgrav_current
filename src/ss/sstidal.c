
/*
 ** sstidal.c -- DCR 7/11/12
 ** ========
 ** Helper for generating tidal encounter simulations.
 ** See etc/sstidal.pdf for geometry.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ss.h>

#define OUTFILENAME "sstidal.ss"
#define LOGFILENAME "sstidal.log"

#define R_EARTH 6371100.0 /* mean Earth radius in m */ 

double SQ(double x);
#define SQ(x) ((x)*(x))

double CUBE(double x);
#define CUBE(x) ((x)*SQ(x))

typedef struct {
  double dPeriapse,dSpeed,dStartDist,dStartSpeed,dSma,dEcc;
  } PARAMS;

static void write_log(const PARAMS *p,const SSDATA *d)
{
	FILE *fp;
	double dDelta,dTimeToPeriapse;
	int nSteps;

	fp = fopen(LOGFILENAME,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",LOGFILENAME);
		exit(1);
		}

	fprintf(fp,"Intruder mass = %g M_Sun (%g kg or %g M_Earth)\n",
			d->mass,d->mass*M_SCALE,d->mass*M_SCALE/M_EARTH);
	fprintf(fp,"Intruder radius = %g AU (%g km or %g R_Earth)\n",
			d->radius,0.001*d->radius*L_SCALE,d->radius*L_SCALE/R_EARTH);
	fprintf(fp,"Close-approach distance = %g AU (%g km or %g R_Earth)\n",
			p->dPeriapse,0.001*p->dPeriapse*L_SCALE,p->dPeriapse*L_SCALE/R_EARTH);
	fprintf(fp,"Close-approach speed = %g x 30 km/s (%g km/s)\n",
			p->dSpeed,0.001*p->dSpeed*V_SCALE);
	fprintf(fp,"Starting distance = %g AU (%g km or %g R_Earth)\n",
			p->dStartDist,0.001*p->dStartDist*L_SCALE,p->dStartDist*L_SCALE/R_EARTH);
	fprintf(fp,"Starting speed = %g x 30 km/s (%g km/s)\n",
			p->dStartSpeed,0.001*p->dStartSpeed*V_SCALE);
	fprintf(fp,"Orbit semimajor axis = %g AU (%g km or %g R_Earth)\n",
			p->dSma,0.001*p->dSma*L_SCALE,p->dSma*L_SCALE/R_EARTH);		
	fprintf(fp,"Orbit eccentricity = %g (%s)\n",p->dEcc,
			p->dEcc == 0.0 ? "circle" : p->dEcc < 1.0 ? "ellipse" : p->dEcc == 1.0 ? "parabola" : "hyperbola");
	fprintf(fp,"Starting position = %g %g %g AU\n",d->pos[0],d->pos[1],d->pos[2]);
	fprintf(fp,"Starting velocity = %g %g %g x 30 km/s\n",d->vel[0],d->vel[1],d->vel[2]);

	/* recommended time step: 1/30th of circular orbit period at periapse */

	dDelta = TWO_PI*sqrt(CUBE(p->dPeriapse)/d->mass)/30.0;

	/* time to encounter: see e.g. Appendix A of Richardson et al. (1998) */

	if (p->dEcc == 0.0)
		dTimeToPeriapse = 0.5*PI*sqrt(CUBE(p->dSma)/d->mass); /* quarter circle */
	else if (p->dEcc < 1.0) {
		double dEccAnom = acos((1.0 - p->dStartDist/p->dSma)/p->dEcc);
		dTimeToPeriapse = sqrt(CUBE(p->dSma)/d->mass)*(p->dEcc*sin(dEccAnom) - dEccAnom);
		}
	else if (p->dEcc == 1.0) {
		double dParAnom = sqrt(p->dStartDist/p->dPeriapse - 1.0);
		dTimeToPeriapse = sqrt(2.0*CUBE(p->dPeriapse)/d->mass)*dParAnom*(1.0 + SQ(dParAnom)/3.0);
		}
	else { /* note: p->dSma is NEGATIVE here... */
		double dHypAnom = acosh((1.0 - p->dStartDist/p->dSma)/p->dEcc);
		dTimeToPeriapse = sqrt(-CUBE(p->dSma)/d->mass)*(p->dEcc*sinh(dHypAnom) - dHypAnom);
		}

	nSteps = (int) (2.0*dTimeToPeriapse/dDelta);

	fprintf(fp,"Recommended timestep = %g yr/2pi (%g s)\n",
			dDelta,dDelta*T_SCALE);
	fprintf(fp,"Time to periapse = %g yr/2pi (%g s)\n",
			dTimeToPeriapse,dTimeToPeriapse*T_SCALE);
	fprintf(fp,"Number of timesteps for full encounter = %i\n",nSteps);

	fclose(fp);
}

static void write_data(SSDATA *d)
{
	SSIO ssio;
	SSHEAD head;

	if (ssioOpen(OUTFILENAME,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",OUTFILENAME);
		exit(1);
		}

	head.time = 0.0;
	head.n_data = 1;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Unable to write header.\n");
		ssioClose(&ssio);
		exit(1);
		}

	d->spin[0] = d->spin[1] = d->spin[2] = 0.0;
	d->color = TEST;
	d->org_idx = 0;

	if (ssioData(&ssio,d)) {
		fprintf(stderr,"Error writing data.\n");
		ssioClose(&ssio);
		exit(1);
		}

	ssioClose(&ssio);
	}

static void set_orbit(PARAMS *p,SSDATA *d)
{
	double dInvSma,dAngMom,dEccSq,dTheta,dPhi,dAlpha;

	dInvSma = 2.0/p->dPeriapse - SQ(p->dSpeed)/d->mass;

	assert(dInvSma != 0.0);

	p->dSma = 1.0/dInvSma;

	dAngMom = p->dPeriapse*p->dSpeed;

	dEccSq = 1.0 - SQ(dAngMom)*dInvSma/d->mass;

	assert(dEccSq >= 0.0);

	p->dEcc = sqrt(dEccSq);

	/* get starting distance or speed from conservation of energy */

	if (p->dStartDist == 0.0) {
		if (p->dStartSpeed == p->dSpeed)
			p->dStartDist = p->dPeriapse;
		else {
			double dRecip =  0.5*(SQ(p->dStartSpeed) + SQ(p->dSpeed))/d->mass + 1.0/p->dPeriapse;
			assert(dRecip > 0.0);
			p->dStartDist = 1.0/dRecip;
			}
		}
	else {
		if (p->dStartDist == p->dPeriapse)
			p->dStartSpeed = p->dSpeed;
		else {
			double dArg = SQ(p->dSpeed) + 2.0*d->mass*(1.0/p->dStartDist - 1.0/p->dPeriapse);
			assert(dArg >= 0.0);
			p->dStartSpeed = sqrt(dArg);
			}
		}

	/* get starting position and velocity from orbit equation and CoAM */

	if (p->dEcc == 0.0) {
		assert(p->dStartDist == p->dPeriapse && p->dStartSpeed == p->dSpeed);
		dTheta = 0.5*PI; /* special case: circular orbit */
		}
	else if (p->dEcc == 1.0)
		dTheta = acos(2.0*p->dPeriapse/p->dStartDist - 1.0); /* parabola: r(1 + cos(theta)) = 2q */
	else {
		if (p->dEcc < 1.0)
			assert(p->dStartDist <= p->dSma*(1.0 + p->dEcc));
		/* general equation: r = a(1 - e^2)/(1 + e cos(theta)) */
		dTheta = acos((p->dSma*(1.0 - dEccSq) - p->dStartDist)/(p->dEcc*p->dStartDist));
		}

	dPhi = asin(dAngMom/(p->dStartDist*p->dStartSpeed)); /* CoAM */

	d->pos[0] = -p->dStartDist*sin(dTheta);
	d->pos[1] = p->dStartDist*cos(dTheta);
	d->pos[2] = 0.0;

	dAlpha = dTheta + dPhi - 0.5*PI;

	d->vel[0] = p->dStartSpeed*cos(dAlpha);
	d->vel[1] = p->dStartSpeed*sin(dAlpha);
	d->vel[2] = 0.0;
	}

static void get_params(PARAMS *p,SSDATA *d)
{
	do {
		printf("Enter the disrupting object (intruder) mass, in Earth masses:\n");
		scanf("%lf",&d->mass);
		if (d->mass <= 0.0)
			fprintf(stderr,"Object mass must be positive.\n");
		} while (d->mass <= 0.0);

	d->mass *= M_EARTH/M_SCALE; /* Earth masses -> solar masses */

	do {
		printf("Enter the disrupting object (intruder) radius, in Earth radii:\n");
		scanf("%lf",&d->radius);
		if (d->radius <= 0.0)
			fprintf(stderr,"Object radius must be positive.\n");
		} while (d->radius <= 0.0);

	d->radius *= R_EARTH/L_SCALE; /* Earth radii -> AU */

	do {
		printf("Enter the close-approach distance (periapse), in km (or < 0 for Earth radii):\n");
		scanf("%lf",&p->dPeriapse);
		if (p->dPeriapse == 0.0)
			fprintf(stderr,"Periapse distance cannot be zero.\n");
		else if (p->dPeriapse < 0.0)
			p->dPeriapse *= -R_EARTH/L_SCALE; /* Earth radii -> AU */
		else
			p->dPeriapse *= 1000.0/L_SCALE; /* km -> AU */
		} while (p->dPeriapse == 0.0);

	do {
		printf("Enter the close-approach speed, in km/s:\n");
		scanf("%lf",&p->dSpeed);
		if (p->dSpeed <= 0.0)
			fprintf(stderr,"Close-approach speed must be positive.\n");
		} while (p->dSpeed <= 0.0);

	p->dSpeed *= 1000.0/V_SCALE; /* km/s -> 2(pi)AU/yr */

	do {
		printf("Enter the starting distance, in km (or < 0 for Earth radii, or = 0 to use starting speed instead):\n");
		scanf("%lf",&p->dStartDist);
		if (p->dStartDist < 0.0)
			p->dStartDist *= -R_EARTH/L_SCALE; /* Earth radii -> AU */
		else if (p->dStartDist > 0.0)
			p->dStartDist *= 1000.0/L_SCALE; /* km -> AU */
		if (p->dStartDist > 0.0 && p->dStartDist < p->dPeriapse)
			fprintf(stderr,"Starting distance cannot be less than periapse distance.\n");
		} while (p->dStartDist > 0.0 && p->dStartDist < p->dPeriapse);

	if (p->dStartDist == 0.0)
		do {
			printf("Enter the starting speed, in km/s:\n");
			scanf("%lf",&p->dStartSpeed);
			if (p->dStartSpeed > p->dSpeed)
				fprintf(stderr,"Starting speed cannot be greater than periapse speed.\n");
			} while (p->dStartSpeed > p->dSpeed);
	else
		p->dStartSpeed = 0.0;
	}

int main(int argc,char *argv[])
{
	PARAMS params;
	SSDATA planet;

	if (argc > 1) {
		fprintf(stderr,"%s takes no arguments.\n",argv[0]);
		return 1;
		}

	printf("This code helps with tidal encounter initial conditions.\n"
		   "The encounter occurs in the body center-of-mass frame.\n"
		   "In this sense, the disrupting object is the \"intruder.\"\n\n");

	get_params(&params,&planet);

	set_orbit(&params,&planet);

	printf("Writing data...\n");
	write_data(&planet);

	write_log(&params,&planet);
	printf("Log written to %s.\n",LOGFILENAME);

	printf("Now run rpg or ssgen to generate the small body initial conditions.\n"
		   "Then run rpx to combine with \"%s\".\n",OUTFILENAME);

	return 0;
}

/* sstidal.c */
