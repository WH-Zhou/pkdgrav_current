/*
 ** demparams.c -- DCR 6/27/12
 ** ===========
 ** Code to assist with setting DEM parameters.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h> /* for tolower() */
#include <math.h>
#include <assert.h>
#include <signal.h>

#define LOGFILE "demparams.log"

double SQ(double);
#define SQ(x) ((x)*(x))

double CUBE(double);
#define CUBE(x) ((x)*SQ(x))

#define L_SCALE 1.496e11 /* AU in m */
#define M_SCALE 1.989e30 /* solar mass in kg */
#define T_SCALE 5.0225e6 /* yr/2pi in s */
#define V_SCALE ((L_SCALE)/(T_SCALE))
#define D_SCALE ((M_SCALE)/CUBE(L_SCALE))
#define COS_45 0.707106781186548
#define CUBEROOT3 1.442249570307408

#define G 6.672e-11 /* gravitation constant in mks */

enum {FALSE=0,TRUE};
typedef int BOOLEAN;

typedef struct {
  int bPKDUnits,nStepsPerOverlap,bCohemode;
  double dMass,dRadius,dVmax,dXmaxOverR,dGravAcc,dHeight,dPackEff,dCoheCoe;
  } PARAMS;

typedef struct {
  double dDensity,dFreeFallSpeed,dFreeFallTime,dBulkDensity,dKnKinetic,dKnStatic,dKnCohesion,dKn,dTauOverlap,dSoundSpeed,dDelta,dDynTimePart,dDynTimeBulk,dCoheStrength;
  } DERIVED;

enum {Quit,Mass,Radius,Density,Speed,Overlap,Gravity,Height,Packing,BulkDensity,Steps,Cohesion,CohesionC};

static BOOLEAN get_yn(const char *str,const char *dflt_str)
{
	enum {none,yes,no} dflt = none;
	
	int c;

	if (dflt_str && strlen(dflt_str)) {
		if (tolower(*dflt_str) == 'y') dflt = yes;
		else if (tolower(*dflt_str) == 'n') dflt = no;
		}

	do {
		printf("%s [%s]? ",str,
			   (dflt == none ? "y/n" : dflt == yes ? "Y/n" : "y/N"));
		c = tolower(getchar());
		if (c == '\n' && dflt != none)
			return dflt == yes;
		else if (c == 'y' || c == 'n') {
			BOOLEAN is_yes = (c == 'y');
			do c = getchar(); /* eat any leftover characters */
			while (c != '\n');
			return is_yes;
			}
		while (c != '\n') c = getchar();
		} while (/*CONSTCOND*/1);
	}

static void analyze(const PARAMS *p,DERIVED *d)
{
        d->dDensity = p->dMass/(4.0/3.0*M_PI*CUBE(p->dRadius));
	d->dFreeFallSpeed = sqrt(2.0*p->dGravAcc*p->dHeight);
	if (p->dGravAcc > 0.0)
		d->dFreeFallTime = d->dFreeFallSpeed/p->dGravAcc;
	else
		d->dFreeFallTime = 0.0;
	d->dBulkDensity = d->dDensity*p->dPackEff;
	d->dKnKinetic = p->dMass*SQ(p->dVmax/(p->dXmaxOverR*p->dRadius));
	d->dKnStatic = p->dGravAcc*d->dBulkDensity*p->dHeight*p->dRadius/p->dXmaxOverR;
        if (p->bCohemode == 0)
                d->dCoheStrength = p->dCoheCoe*4.0*SQ(0.5*p->dRadius); /* assumes beta_s = 0.5 */
        else
                d->dCoheStrength = p->dCoheCoe*M_PI*SQ(p->dRadius);
        d->dKnCohesion = d->dCoheStrength/(p->dXmaxOverR*p->dRadius);
        d->dKn = d->dKnKinetic > d->dKnStatic ? d->dKnKinetic : d->dKnStatic;
        d->dKn = d->dKnCohesion > d->dKn ? d->dKnCohesion : d->dKn;
	d->dTauOverlap = M_PI*sqrt(0.5*p->dMass/d->dKn); /* assumes equal masses and zero damping coefficient */
	d->dSoundSpeed = 2.0*p->dRadius/d->dTauOverlap;
	d->dDelta = d->dTauOverlap/p->nStepsPerOverlap;
	d->dDynTimePart = 1.0/sqrt(G*d->dDensity);
        d->dDynTimeBulk = 1.0/sqrt(G*d->dBulkDensity);
	}

static double dDelta(const PARAMS *p, double dKn)
{
    /*
    ** Returns recommended timestep for given normal spring constant,
    ** given the typical partical mass and the desired number of steps
    ** per overlap.  If dKn == 0, zero is returned for dDelta.
    */

    if (dKn == 0.0)
	return 0.0;

    return M_PI * sqrt(0.5 * p->dMass / dKn) / p->nStepsPerOverlap;
    }

static void show_menu(const PARAMS *p,const DERIVED *d,FILE *fp)
{
	fprintf(fp,"%2i. Typical (energetic) particle mass = ",Mass);
	if (p->bPKDUnits) fprintf(fp,"%g M_Sun",p->dMass/M_SCALE);
	else fprintf(fp,"%g g",p->dMass*1000.0);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. Typical (energetic) particle radius = ",Radius);
	if (p->bPKDUnits) fprintf(fp,"%g AU",p->dRadius/L_SCALE);
	else fprintf(fp,"%g cm",p->dRadius*100.0);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. [DERIVED] Particle density = ",Density);
	if (p->bPKDUnits) fprintf(fp,"%g M_Sun/AU^3",d->dDensity/D_SCALE);
	else fprintf(fp,"%g g/cc",d->dDensity/1000.0);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. Maximum expected particle speed = ",Speed);
	if (p->bPKDUnits) fprintf(fp,"%g x 30 km/s",p->dVmax/V_SCALE);
	else fprintf(fp,"%g cm/s",p->dVmax*100.0);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. Maximum fractional overlap = %g (times typical radius)\n",Overlap,p->dXmaxOverR);

	fprintf(fp,"%2i. Acceleration due to gravity (uniform or surface) = ",Gravity);
	if (p->bPKDUnits) fprintf(fp,"%g AU/(yr/2pi)^2",p->dGravAcc/L_SCALE*SQ(T_SCALE));
	else fprintf(fp,"%g m/s^2",p->dGravAcc);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. Maximum system height (or bulk radius) = ",Height);
	if (p->bPKDUnits) fprintf(fp,"%g AU",p->dHeight/L_SCALE);
	else fprintf(fp,"%g cm",p->dHeight*100.0);
	fprintf(fp,"\n");

	fprintf(fp,"    [DERIVED: free-fall speed = ");
	if (p->bPKDUnits) fprintf(fp,"%g x 30 km/s",d->dFreeFallSpeed/V_SCALE);
	else fprintf(fp,"%g cm/s",d->dFreeFallSpeed*100.0);
	fprintf(fp,"; free-fall time = ");
	if (p->bPKDUnits) fprintf(fp,"%g yr/2pi",d->dFreeFallTime/T_SCALE);
	else fprintf(fp,"%g s",d->dFreeFallTime);
	fprintf(fp,"]\n");

	fprintf(fp,"%2i. Packing efficiency = %g\n",Packing,p->dPackEff);

	fprintf(fp,"%2i. [DERIVED] Bulk density = ",BulkDensity);
	if (p->bPKDUnits) fprintf(fp,"%g M_Sun/AU^3",d->dBulkDensity/D_SCALE);
	else fprintf(fp,"%g g/cc",d->dBulkDensity/1000.0);
	fprintf(fp,"\n");

	fprintf(fp,"%2i. Steps per overlap = %i\n",Steps,p->nStepsPerOverlap);

        fprintf(fp,"%2i. iCohesionModel = %i (0:contact area depends on beta_s, 1:constant area)\n",Cohesion,p->bCohemode);
        fprintf(fp,"%2i. Cohesive strength between particles in contact = %g Pa\n",CohesionC,p->dCoheCoe);

	fprintf(fp,"\nCurrent best estimates for DEM parameters:\n\n");

	fprintf(fp,"Kinetic dKn = ");
	if (!p->bPKDUnits) fprintf(fp,"%e kg/s^2 = ",d->dKnKinetic);
	fprintf(fp,"%e M_Sun/(yr/2pi)^2",d->dKnKinetic/M_SCALE*SQ(T_SCALE));
        if (d->dKnKinetic > d->dKnStatic && d->dKnKinetic > d->dKnCohesion) fprintf(fp," (recommended)");
	fprintf(fp,"\n");

	fprintf(fp,"Static dKn  = ");
	if (!p->bPKDUnits) fprintf(fp,"%e kg/s^2 = ",d->dKnStatic);
	fprintf(fp,"%e M_Sun/(yr/2pi)^2",d->dKnStatic/M_SCALE*SQ(T_SCALE));
        if (d->dKnKinetic <= d->dKnStatic && d->dKnCohesion <= d->dKnStatic) fprintf(fp," (recommended)");
	fprintf(fp,"\n");

        fprintf(fp,"Cohesion dKn  = ");
        if (!p->bPKDUnits) fprintf(fp,"%e kg/s^2 = ",d->dKnCohesion);
        fprintf(fp,"%e M_Sun/(yr/2pi)^2",d->dKnCohesion/M_SCALE*SQ(T_SCALE));
        if (d->dKnKinetic < d->dKnCohesion && d->dKnStatic < d->dKnCohesion) fprintf(fp," (recommended)");
        fprintf(fp,"\n");

	fprintf(fp, "\n");

	fprintf(fp,"Kinetic dDelta = ");
	if (!p->bPKDUnits) fprintf(fp, "%e s = ", dDelta(p, d->dKnKinetic));
	fprintf(fp, "%e yr/2pi", dDelta(p, d->dKnKinetic) / T_SCALE);
        if (d->dKnKinetic > d->dKnStatic && d->dKnKinetic > d->dKnCohesion) fprintf(fp," (recommended)");
	fprintf(fp, "\n");

	fprintf(fp,"Static dDelta = ");
	if (!p->bPKDUnits) fprintf(fp, "%e s = ", dDelta(p, d->dKnStatic));
	fprintf(fp,"%e yr/2pi", dDelta(p, d->dKnStatic) / T_SCALE);
        if (d->dKnKinetic <= d->dKnStatic && d->dKnCohesion <= d->dKnStatic) fprintf(fp," (recommended)");
	fprintf(fp,"\n");

	fprintf(fp,"Cohesion dDelta = ");
	if (!p->bPKDUnits) fprintf(fp, "%e s = ", dDelta(p, d->dKnCohesion));
	fprintf(fp, "%e yr/2pi", dDelta(p, d->dKnCohesion)/T_SCALE);
        if (d->dKnKinetic < d->dKnCohesion && d->dKnStatic < d->dKnCohesion) fprintf(fp," (recommended)");
	fprintf(fp, "\n");

	/* paranoia check to make sure recommended dDelta matched dKn */
	assert(fabs(dDelta(p, d->dKn) - d->dDelta) / d->dDelta < 1.0e-0);

	fprintf(fp, "\n");

        fprintf(fp,"dCohesiveCoeff = ");
        if (!p->bPKDUnits) fprintf(fp,"%e Pa = ",p->dCoheCoe);
        fprintf(fp,"%e M_Sun/AU/(yr/2pi)^2",p->dCoheCoe/M_SCALE*L_SCALE*T_SCALE*T_SCALE);
	fprintf(fp,"\n");

	fprintf(fp,"\n");

	fprintf(fp,"NOTE: Typical overlap duration  = ");
	if (p->bPKDUnits) fprintf(fp,"%g yr/2pi",d->dTauOverlap/T_SCALE);
	else fprintf(fp,"%g s",d->dTauOverlap);
	fprintf(fp,"\n");

	fprintf(fp,"      Typical propagation speed = ");
	if (p->bPKDUnits) fprintf(fp,"%g 2pi AU/yr",d->dSoundSpeed/V_SCALE);
	else fprintf(fp,"%g m/s",d->dSoundSpeed);
	fprintf(fp,"\n");

	fprintf(fp,"      Particle dynamical time = ");
	if (p->bPKDUnits) fprintf(fp,"%g yr/2pi",d->dDynTimePart/V_SCALE);
	else fprintf(fp,"%g s",d->dDynTimePart);
	fprintf(fp," (number of timesteps = %lli)\n",(long long int) (d->dDynTimePart/d->dDelta));

	fprintf(fp,"      Bulk dynamical time = ");
	if (p->bPKDUnits) fprintf(fp,"%g yr/2pi",d->dDynTimeBulk/V_SCALE);
	else fprintf(fp,"%g s",d->dDynTimeBulk);
	fprintf(fp," (number of timesteps = %lli)\n",(long long int) (d->dDynTimeBulk/d->dDelta));

	fprintf(fp,"      Cohesive force = ");
	if (p->bPKDUnits) fprintf(fp,"%g M_Sun*AU/(yr/2pi)^2",d->dCoheStrength/M_SCALE/L_SCALE*T_SCALE*T_SCALE);
	else fprintf(fp,"%g N",d->dCoheStrength);
	if (p->dGravAcc > 0) fprintf(fp," = %g particle's weight",d->dCoheStrength/(p->dMass*p->dGravAcc));
	if (!p->bCohemode) fprintf(fp," (for beta_s = 0.5)");
	fprintf(fp,"\n");

	}

static void sighandler(int signumber)
{
	(void) signumber; /* to suppress unused parameter warning */
	printf("\nEnter number to change (or 0 when done): ");
	fflush(stdout);
	}

static int process(void)
{
	PARAMS params;
	DERIVED derived;
	FILE *fp = NULL;
	BOOLEAN bTest;
	double d;
	int iChoice,i;

	printf("\n*** Invoke with any argument to get a helpful usage message! ***\n\n");

	params.bPKDUnits = get_yn("Use PKDGRAV units (AU, M_sun, etc.)","n");

	params.dMass = 4.0/3.0*M_PI/1.0e3; /* kg (1 cm radius @ 1 g/cc) */
	params.dRadius = 0.01; /* m (1 cm) */
	params.dVmax = 0.01; /* m/s (1 cm/s) */
	params.dXmaxOverR = 0.01;
	params.dGravAcc = 9.8; /* m/s^2 */
	params.dHeight = 1.0; /* m */
	params.dPackEff = 0.65;
	params.nStepsPerOverlap = 30;
        params.bCohemode = 0;
        params.dCoheCoe = 0;
	
	while (/*CONSTCOND*/1) {

		analyze(&params,&derived);

		printf("\n");

		show_menu(&params,&derived,stdout);

		printf("\n");

		printf("Enter number to change (or 0 when done): ");
		signal(SIGINT,sighandler); /* intercept CTRL-C */
		scanf("%i",&iChoice);
		signal(SIGINT,SIG_DFL);

		getchar(); /* chomp \n since getchar() may be needed later */

		switch (iChoice) {
		case Quit:
			fp = fopen(LOGFILE,"w");
			assert(fp != NULL);
			show_menu(&params,&derived,fp);
			fclose(fp);
			return 0;
			}

		switch (iChoice) {
		case Mass:
			bTest = get_yn("Keep particle density constant","n");
			do {
				printf("Enter new mass (-ve ==> scale mass instead): ");
				scanf("%lf",&d);
				if (d == 0.0) printf("Mass cannot be zero.\n");
				} while (d == 0.0);
			getchar();
			if (d < 0.0)
				d *= -params.dMass;
			else {
				if (params.bPKDUnits) d *= M_SCALE; /* M_Sun -> kg */
				else d /= 1000.0; /* g -> kg */
				}
			if (bTest) params.dRadius *= pow(d/params.dMass,1.0/3.0);
			params.dMass = d;
			break;
		case Radius:
			bTest = get_yn("Keep particle density constant","n");
			do {
				printf("Enter new radius (-ve ==> scale radius instead): ");
				scanf("%lf",&d);
				if (d == 0.0) printf("Radius cannot be zero.\n");
				} while (d == 0.0);
			getchar();
			if (d < 0.0)
				d *= -params.dRadius;
			else {
				if (params.bPKDUnits) d *= L_SCALE; /* AU -> m */
				else d /= 100.0; /* cm -> m */
				}
			if (bTest) params.dMass *= CUBE(d/params.dRadius);
			params.dRadius = d;
			break;
		case Density:
			bTest = get_yn("Keep radius constant","y");
			do {
				printf("Enter new density (-ve ==> scale density instead): ");
				scanf("%lf",&d);
				if (d == 0.0) printf("Density cannot be zero.\n");
				} while (d == 0);
			getchar();
			if (d < 0.0)
				d *= -derived.dDensity;
			else {
				if (params.bPKDUnits) d *= D_SCALE; /* M_Sun/AU^3 -> kg/m^3 */
				else d *= 1000.0; /* g/cc -> kg/m^3 */
				}
			if (bTest) params.dMass *= d/derived.dDensity;
			else params.dRadius *= pow(derived.dDensity/d,1.0/3.0);
			break;
		case Speed:
			printf("Enter new speed (-ve ==> scale speed instead): ");
			scanf("%lf",&d);
			getchar();
			if (d < 0.0)
				d *= -params.dVmax;
			else {
				if (params.bPKDUnits) d *= V_SCALE; /* AU/(yr/2pi) -> m/s */
				else d /= 100.0; /* cm/s -> m/s */
				}
			params.dVmax = d;
			break;
		case Overlap:
			do {
				printf("Enter new fractional overlap: ");
				scanf("%lf",&d);
				if (d <= 0.0) printf("Fraction overlap must be positive.\n");
				} while (d <= 0.0);
			getchar();
			params.dXmaxOverR = d;
			break;
		case Gravity:
			printf("Enter new gravitational acceleration (-ve ==> scale gravity instead): ");
			scanf("%lf",&d);
			getchar();
			if (d < 0.0)
				d *= -params.dGravAcc;
			else
				if (params.bPKDUnits) d *= L_SCALE/SQ(T_SCALE); /* AU/(yr/2pi)^2 --> m/s^2 */
			params.dGravAcc = d;
			break;
		case Height:
			printf("Enter new system height/bulk radius (-ve ==> scale instead): ");
			scanf("%lf",&d);
			getchar();
			if (d < 0.0)
				d *= -params.dHeight;
			else {
				if (params.bPKDUnits) d *= L_SCALE; /* AU --> m */
				else d /= 100.0; /* cm -> m */
				}
			params.dHeight = d;
			break;
		case Packing:
			do {
				printf("Enter new packing efficiency: ");
				scanf("%lf",&d);
				if (d < 0.0 || d > 1.0) printf("Packing efficiency must be between 0 and 1.\n");
				} while (d < 0.0 || d > 1.0);
			getchar();
			params.dPackEff = d;
			break;
		case BulkDensity:
			bTest = get_yn("OK to change packing efficiency","y");
			if (!bTest) break;
			do {
				printf("Enter new bulk density (-ve ==> scale instead): ");
				scanf("%lf",&d);
				if (d == 0.0) printf("Bulk density cannot be zero.\n");
				} while (d == 0);
			getchar();
			if (d < 0.0)
				d *= -derived.dBulkDensity;
			else {
				if (params.bPKDUnits) d *= D_SCALE; /* M_Sun/AU^3 -> kg/m^3 */
				else d *= 1000.0; /* g/cc -> kg/m^3 */
				}
			params.dPackEff = d/derived.dDensity;
			break;
		case Steps:
			do {
				printf("Enter new number of steps per overlap: ");
				scanf("%i",&i);
				if (i < 1) printf("Number of steps must be positive.\n");
				} while (i < 1);
			getchar();
			params.nStepsPerOverlap = i;
			break;
		case Cohesion:
			do {
				printf("Enter the cohesion model (0:contact area depends on beta_s, 1: constant contact area): ");
				scanf("%i",&i);
				if ( i != 1 && i !=0 ) printf("The cohesion model must be 0 or 1.\n");
				} while ( i != 1 && i !=0 );
			getchar();
			params.bCohemode = i;
			break;
		case CohesionC:
			printf("Enter the cohesive coefficient (-ve ==> scale cohesive coefficient instead): ");
			scanf("%lf",&d);
			getchar();
			if (d < 0.0)
				d *= -params.dCoheCoe;
			else if (params.bPKDUnits) d *= M_SCALE/L_SCALE/T_SCALE/T_SCALE;
			params.dCoheCoe = d;
			break;
		default:
			printf("Invalid menu choice.\n");
			}
		}
	}

int main(int argc,char *argv[])
{
	if (argc > 1) {
		fprintf(stderr,"\n%s takes no arguments\n\n",argv[0]);
		fprintf(stderr,"Helpful hints:\n\n");
		fprintf(stderr,"1. Kinetic dKn is proportional to M, V^2, 1/R^2, 1/f^2.\n");
		fprintf(stderr,"2. Static dKn is proportional to g, rho, h, R, 1/f.\n");
		fprintf(stderr,"3. Cohesion dKn is proportional to c, 1/R, 1/f.\n");
		fprintf(stderr,"4. dDelta is proportional to sqrt(M), 1/sqrt(dKn), 1/n.\n");
		fprintf(stderr,"\nHere M = particle mass, V = particle speed, R = particle radius,\nf = fractional overlap, g = acceleration, rho = bulk density,\nh = height, c = cohesive force, and n = no. steps per overlap.\nNote therefore dDelta does NOT depend on M in the kinetic case.\n");
		fprintf(stderr,"\nAlso, any friction will tend to lower dKn and increase dDelta,\nbut excessive friction may ultimately require smaller dDelta.\n");
		fprintf(stderr,"\nFor polydisperse particles (including aggregates and reactive\nwalls), take the most massive particle (or aggregate or wall)\nfor M and the smallest particle radius for R (worst case).\n");
		fprintf(stderr,"\nFor a single particle on a horizontal wall, use dKn = M * g / x,\nwhere x = f * R; then dDelta = pi * sqrt(0.5 * M / dKn) / n.\n");
		fprintf(stderr,"\nFinally, remember wall motions for maximum relative speeds.\nWant (v_max + v_wall_max) * step << smallest particle radius.\n");
		fprintf(stderr,"\nCalculations use certain assumptions -- see code for details.\n");
		fprintf(stderr,"\n");
		return 1;
		}

	process();

	return 0;
	}

/* demparams.c */
