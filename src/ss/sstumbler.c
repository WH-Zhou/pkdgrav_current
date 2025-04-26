/*
 ** sstumbler.c -- DCR 02/07/12
 ** ===========
 ** Generates initial conditions in a uniform grid on the inner
 ** surface of a cylinder.  NOTE: cylinder origin at (0,0,0) by
 ** default.
 **
 ** UPDATE 2/19/12: now optionally also fills cylinder on uniform
 ** (square) grid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <time.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <random.h>
#include <colors.h> /* for NUM_COLORS */

#define OUTFILENAME "sstumbler.ss"
#define LOGFILENAME "sstumbler.log"

typedef struct {
  int bOffset,nStuckTarget,nStuckR,nStuckZ,nStuck;
  int nFreeTarget,nFreeXY,nFreeZ,nFree,nTarget,nPart,iColor,iWallID;
  double dCylRad,dCylLen,dStuckAvg,dStuckDev,dDepthAvg,dDepthDev;
  double dFreeAvg,dFreeDev,dCutoff,dDensity,dStuckR,dStuckZ,dFreeXY,dFreeZ;
  } PARAMS;

typedef struct {
  int color;
  double pos[N_DIM],radius;
  } DATA;

double SQ(double);
#define SQ(x) ((x)*(x))

static void write_data(const PARAMS *p,const DATA *d)
{
	SSIO ssio;
	SSHEAD head;
	SSDATA data;
	int i;

	if (ssioOpen(OUTFILENAME,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",OUTFILENAME);
		exit(1);
		}

	head.time = 0.0;
	head.n_data = p->nPart;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Unable to write header.\n");
		ssioClose(&ssio);
		exit(1);
		}

	for (i=0;i<p->nPart;i++) {
		data.radius = d[i].radius;
		data.mass = 4.0/3.0*PI*data.radius*data.radius*data.radius*p->dDensity;
		data.pos[0] = d[i].pos[0];
		data.pos[1] = d[i].pos[1];
		data.pos[2] = d[i].pos[2] + (p->bOffset ? 0.5*p->dCylLen : 0.0);
		data.vel[0] = data.vel[1] = data.vel[2] = 0.0;
		data.spin[0] = data.spin[1] = data.spin[2] = 0.0;
		if (i < p->nStuck)
			data.color = -1 - p->iWallID;
		else
			data.color = p->iColor;
		data.org_idx = i;
		if (ssioData(&ssio,&data)) {
			fprintf(stderr,"Error writing data (particle %i).\n",i);
			ssioClose(&ssio);
			exit(1);
			}
		}

	ssioClose(&ssio);
	}

static int compar(const void *v1,const void *v2)
{
	if (((const DATA *)v1)->radius < ((const DATA *)v2)->radius)
		return -1;
	else if (((const DATA *)v1)->radius > ((const DATA *)v2)->radius)
		return 1;
	else
		return 0;
	}

static void write_log(const PARAMS *p,DATA *d)
{
	FILE *fp;
	double x,s;
	int i;

	fp = fopen(LOGFILENAME,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",LOGFILENAME);
		exit(1);
		}

	fprintf(fp,"Cylinder radius = %g cm\n",100.0*p->dCylRad*AU);
	fprintf(fp,"Cylinder length = %g cm\n",100.0*p->dCylLen*AU);
	fprintf(fp,"Cylinder offset = %g cm\n",p->bOffset*50.0*p->dCylRad*AU);
	fprintf(fp,"Target average stuck-particle radius = %g cm\n",100.0*p->dStuckAvg*AU);
	fprintf(fp,"Target stuck-particle radius uncertainty = %g cm\n",100.0*p->dStuckDev*AU);
	fprintf(fp,"Target average depth fraction = %g\n",p->dDepthAvg);
	fprintf(fp,"Target depth fraction uncertainty = %g\n",p->dDepthDev);
	fprintf(fp,"Target average free-particle radius = %g cm\n",100.0*p->dFreeAvg*AU);
	fprintf(fp,"Target free-particle radius uncertainty = %g cm\n",100.0*p->dFreeDev*AU);
	fprintf(fp,"Gaussian distribution cutoff = %g sigma\n",p->dCutoff);
	fprintf(fp,"Particle density = %g g/cc\n",1.0e-3*p->dDensity*M_SUN/(AU*AU*AU));
	fprintf(fp,"Particle color = %i\n",p->iColor);
	fprintf(fp,"Wall ID = %i\n",p->iWallID);
	fprintf(fp,"Cylinder goes from %g to %g cm in z\n",p->bOffset ? 0.0 : -50.0*p->dCylLen*AU,p->bOffset ? 100.0*p->dCylLen*AU : 50.0*p->dCylLen*AU);
	fprintf(fp,"Total number of stuck particles = %i",p->nStuck);
	if (p->nStuckTarget == -1)
		fprintf(fp," (maximum)");
	else if (p->nStuckTarget != 0)
		fprintf(fp," (aimed for %i)",p->nStuckTarget);
	fprintf(fp,"\n");
	if (p->nStuck > 0) {
		fprintf(fp,"Dimensions of stuck-particle grid = %ix(%i/%i)\n",p->nStuckR,p->nStuckZ,p->nStuckZ > 1 ? p->nStuckZ - 1 : 1);
		fprintf(fp,"Dimensions of stuck-particle grid cell = %gx%g cm (fill factor %g%% of maximum)\n",100.0*p->dStuckR*AU,100.0*p->dStuckR*AU,400.0*SQ(p->dStuckAvg)/SQ(p->dStuckR));
		}
	fprintf(fp,"Total number of free particles = %i",p->nFree);
	if (p->nFreeTarget == -1)
		fprintf(fp," (maximum)");
	else if (p->nFreeTarget != 0)
		fprintf(fp," (aimed for %i)",p->nFreeTarget);
	fprintf(fp,"\n");
	if (p->nFree > 0) {
		fprintf(fp,"Dimensions of free-particle grid = %ix%ix%i\n",p->nFreeXY,p->nFreeXY,p->nFreeZ);
		fprintf(fp,"Dimensions of free-particle grid cell = %gx%gx%g cm (fill factor %g%% of maximum)\n",100.0*p->dFreeXY*AU,100.0*p->dFreeXY*AU,100*p->dFreeZ*AU,800.0*p->dFreeAvg*SQ(p->dFreeAvg)/(SQ(p->dFreeXY)*p->dFreeZ));
		}
	fprintf(fp,"Total number of particles = %i\n",p->nPart);

	if (p->nStuck > 1) {
		x = 0.0;
		for (i=0;i<p->nStuck;i++)
			x += d[i].radius;
		x /= p->nStuck;
		s = 0.0;
		for (i=0;i<p->nStuck;i++)
			s += SQ(d[i].radius - x);
		s = sqrt(s/(p->nStuck - 1));
		fprintf(fp,"Average stuck-particle radius = %g +/- %g cm\n",100.0*x*AU,100.0*s*AU);
		}

	if (p->nFree > 1) {
		x = 0.0;
		for (i=p->nStuck;i<p->nPart;i++)
			x += d[i].radius;
		x /= p->nFree;
		s = 0.0;
		for (i=p->nStuck;i<p->nPart;i++)
			s += SQ(d[i].radius - x);
		s = sqrt(s/(p->nFree - 1));
		fprintf(fp,"Average free-particle radius = %g +/- %g cm\n",100.0*x*AU,100.0*s*AU);
		}

	qsort(d,(size_t) p->nPart,sizeof(DATA),compar); /* put particles in ascending-radius order */

	fprintf(fp,"Minimum particle radius = %g cm\n",100.0*d[0].radius*AU);
	fprintf(fp,"Maximum particle radius = %g cm\n",100.0*d[p->nPart - 1].radius*AU);
	fprintf(fp,"Suggested minimum nSmooth = %i\n",4*((int) (SQ(d[p->nPart - 1].radius/d[0].radius) + 0.5)));

	fclose(fp);
	}

static double gaussian(double dSig,double dLim)
{
	/* truncated Gaussian */

	double dVal;

	do {
		dVal = dSig*randGaussian();
		} while (dLim > 0.0 && fabs(dVal) > dLim*dSig);

	return dVal;
	}

static void generate(const PARAMS *p,DATA *d)
{
	double dTheta,dSpaceTheta,dSpaceZ,dDepth,dSpaceXY;
	int i,iR,iZ,iX,iY;

	randSeedGenerator((Ullong) (time(NULL) % getpid() + getppid()));

	i = 0;

	if (p->nStuck > 0) {
		printf("Generating %i stuck particle%s in %ix(%i/%i) grid.\n",
			   p->nStuck,p->nStuck == 1 ? "" : "s",
			   p->nStuckR,p->nStuckZ,p->nStuckZ > 1 ? p->nStuckZ - 1 : 1);

		dTheta = TWO_PI/p->nStuckR;

		printf("Stuck-particle cell dimensions: %g cm x %g cm (angular increment = %g deg).\n",
			   100.0*p->dStuckR*AU,100.0*p->dStuckZ*AU,dTheta*180.0/PI);

		for (iZ=0;iZ<p->nStuckZ;iZ++)
			for (iR=0;iR<p->nStuckR;iR++) {
				if (p->nStuckZ > 1 && iR % 2 == 1 && iZ == 0)
					continue; /* to stagger the grid */
				do {
					d[i].radius = p->dStuckAvg + gaussian(p->dStuckDev,p->dCutoff);
					} while (d[i].radius <= 0.0); /* reject non-physical values */
				dSpaceTheta = (p->dStuckR - 2.0*d[i].radius)/p->dCylRad; /* available longitudinal space in cell */
				if (dSpaceTheta < 0.0)
					dSpaceTheta = 0.0;
				dSpaceZ = p->dStuckZ - 2.0*d[i].radius; /* available vertical space in cell */
				if (dSpaceZ < 0.0)
					dSpaceZ = 0.0;
				dDepth = 2.0*(0.5 - (p->dDepthAvg + gaussian(p->dDepthDev,p->dCutoff)))*d[i].radius;
				d[i].pos[0] = (p->dCylRad - dDepth)*cos(dTheta*iR + (randUniform() - 0.5)*dSpaceTheta);
				d[i].pos[1] = (p->dCylRad - dDepth)*sin(dTheta*iR + (randUniform() - 0.5)*dSpaceTheta);
				d[i].pos[2] = -0.5*p->dCylLen + ((double) iZ + 0.5)*p->dStuckZ + (randUniform() - 0.5)*dSpaceZ - (p->nStuckZ > 1 ? (iR % 2)*0.5*p->dStuckZ : 0.0);
				++i;
				}
		}

	assert(i == p->nStuck);

	if (p->nFree > 0) {
		printf("Generating %i free particle%s in %ix%ix%i grid.\n",
			   p->nFree,p->nFree == 1 ? "" : "s",
			   p->nFreeXY,p->nFreeXY,p->nFreeZ);

		printf("Free-particle cell dimensions: %g cm x %g cm x %g cm.\n",
			   100.0*p->dFreeXY*AU,100.0*p->dFreeXY*AU,100.0*p->dFreeZ*AU);

		for (iZ=0;iZ<p->nFreeZ;iZ++)
			for (iX=0;iX<p->nFreeXY;iX++)
				for (iY=0;iY<p->nFreeXY;iY++) {
					do {
						do {
							d[i].radius = p->dFreeAvg + gaussian(p->dFreeDev,p->dCutoff);
							} while (d[i].radius <= 0.0);
						dSpaceXY = (p->dFreeXY - 2.0*d[i].radius); /* available xy-space in cell */
						} while (dSpaceXY <= 0.0);
					dSpaceZ = p->dFreeZ - 2.0*d[i].radius; /* available z-space in cell */
					if (dSpaceZ < 0.0)
						dSpaceZ = 0.0;
					d[i].pos[0] = -sqrt(0.5)*p->dCylRad + ((double) iX + 0.5)*p->dFreeXY + (randUniform() - 0.5)*dSpaceXY;
					d[i].pos[1] = -sqrt(0.5)*p->dCylRad + ((double) iY + 0.5)*p->dFreeXY + (randUniform() - 0.5)*dSpaceXY;
					d[i].pos[2] = -0.5*p->dCylLen + ((double) iZ + 0.5)*p->dFreeZ + (randUniform() - 0.5)*dSpaceZ;
					++i;
					}
		}

	assert(i == p->nPart);
	}

static void usage(const char *achProgName)
{
	fprintf(stderr,
			"Usage: %s -r cylrad -l cyllen [ -v ] -m nstuck [ -s stuckrad [ -a stuckdev ] -h height [ -z hgtdev ] ] -n nfree [ -f freerad [ -b freedev ] ] [ -q cutoff ] -d density [ -c color ] [ -w wall ]\n", achProgName);
	fprintf(stderr, "where: cylrad is the cylinder radius (cm);\n");
	fprintf(stderr, "       cyllen is the cylinder length (cm, along z-axis);\n");
	fprintf(stderr, "       -v indicates cylinder bottom is to be at z = 0;\n");
	fprintf(stderr, "       nstuck is the target number of stuck particles (default -1)\n");
	fprintf(stderr, "          (0 means no stuck particles, -1 means maximum);\n");
	fprintf(stderr, "       stuckrad is the mean stuck-particle radius (cm);\n");
	fprintf(stderr, "       stuckdev is the uncertainty in stuck-particle radius (cm);\n");
	fprintf(stderr, "       height is the mean particle depth on cylinder surface\n");
	fprintf(stderr, "          (0 means just touching, 1 means completely buried);\n");
	fprintf(stderr, "       hgtdev is the uncertainty in particle depth;\n");
	fprintf(stderr, "       nfree is the target number of free particles (default -1)\n");
	fprintf(stderr, "          (0 means no free particles, -1 means maximum);\n");
	fprintf(stderr, "       freerad is the mean free-particle radius (cm);\n");
	fprintf(stderr, "       freedev is the uncertainty in free-particle radius (cm);\n");
	fprintf(stderr, "       cutoff is max. no. sigma for dispersions (default 1; 0 means no limit);\n");
	fprintf(stderr, "       density is the particle density (g/cc);\n");
	fprintf(stderr, "       color is the free-particle color (default 3 = green);\n");
	fprintf(stderr, "       and wall is the wall ID to use (default 0).\n");
	fprintf(stderr, "NOTE: either nstuck or nfree must be non-zero.\n");
	exit(0);
}

int main(int argc,char *argv[])
{
	PARAMS params;
	DATA *data = NULL;
	double d;
	int c,i;

	/* defaults */

	params.dCylRad = 0.0;
	params.dCylLen = 0.0;
	params.bOffset = 0;
	params.nStuckTarget = -1;
	params.dStuckAvg = 0.0;
	params.dStuckDev = 0.0;
	params.dDepthAvg = 0.0;
	params.dDepthDev = 0.0;
	params.nFreeTarget = -1;
	params.dFreeAvg = 0.0;
	params.dFreeDev = 0.0;
	params.dCutoff = 1.0;
	params.dDensity = 0.0;
	params.iColor = PLANETESIMAL; /* 3 = green */
	params.iWallID = 0;

	params.nStuck = params.nStuckZ /* for compiler */ = params.nFree = 0;

	/* parse command-line arguments */

	while ((c = getopt(argc,argv,"r:l:m:s:a:h:z:n:f:b:q:d:w:c:v")) != EOF)
		switch (c) {
		case 'r':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dCylRad = 0.01*d/AU; /* cm -> AU */
			break;
		case 'l':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dCylLen = 0.01*d/AU; /* cm -> AU */
			break;
		case 'v':
			params.bOffset = 1;
			break;
		case 'm':
			i = atoi(optarg);
			if (i < -1)
				usage(argv[0]);
			params.nStuckTarget = i;
			break;
		case 's':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dStuckAvg = 0.01*d/AU; /* cm -> AU */
			break;
		case 'a':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dStuckDev = 0.01*d/AU; /* cm -> AU */
			break;
		case 'h':
			d = atof(optarg);
			if (d < 0.0 || d > 1.0)
				usage(argv[0]);
			params.dDepthAvg = d;
			break;
		case 'z':
			d = atof(optarg);
			if (d < 0.0 || d > 1.0)
				usage(argv[0]);
			params.dDepthDev = d;
			break;
		case 'n':
			i = atoi(optarg);
			if (i < -1)
				usage(argv[0]);
			params.nFreeTarget = i;
			break;
		case 'f':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dFreeAvg = 0.01*d/AU; /* cm -> AU */
			break;
		case 'b':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dFreeDev = 0.01*d/AU; /* cm -> AU */
			break;
		case 'q':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dCutoff = d;
			break;
		case 'd':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dDensity = 1.0e3*d/M_SUN*AU*AU*AU; /* g/cc -> M_SUN/AU^3 */
			break;
		case 'c':
			i = atoi(optarg);
			if (i < 0 || i >= NUM_COLORS) {
				fprintf(stderr,"%s: Invalid color specified (%i).\n",argv[0],i);
				usage(argv[0]);
				}
			params.iColor = i;
			break;
		case 'w':
			i = atoi(optarg);
			if (i < 0)
				usage(argv[0]);
			params.iWallID = i;
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind > argc || params.dCylRad == 0.0 || params.dDensity == 0.0)
		usage(argv[0]);

	if (TWO_PI*params.dCylRad <= params.dCylLen)
		fprintf(stderr,"WARNING: Algorithm assumes cylinder perimeter exceeds cylinder length.\n");

	if (params.nStuckTarget != 0) {
		if (params.dStuckAvg == 0.0)
			usage(argv[0]);
		if (params.dCylRad < params.dStuckAvg)
			fprintf(stderr,"WARNING: cylinder radius < average stuck-particle radius -- very small!\n");
		else if (params.dCylRad < 10.0*params.dStuckAvg)
			fprintf(stderr,"WARNING: cylinder radius < 10x average stuck-particle radius.\n");
		if (params.dCylLen < params.dStuckAvg)
			fprintf(stderr,"WARNING: stuck-particle radius exceeds cylinder length.\n");
		if (params.dStuckDev >= params.dStuckAvg)
			fprintf(stderr,"WARNING: large stuck-particle radius uncertainty---infinite loop possible.\n");
		}

	if (params.nFreeTarget != 0) {
		if (params.dFreeAvg == 0.0)
			usage(argv[0]);
		if (params.dCylRad < params.dFreeAvg)
			fprintf(stderr,"WARNING: cylinder radius < average free-particle radius -- very small!\n");
		else if (params.dCylRad < 10.0*params.dFreeAvg)
			fprintf(stderr,"WARNING: cylinder radius < 10x average free-particle radius.\n");
		if (params.dCylLen < params.dFreeAvg)
			fprintf(stderr,"WARNING: free-particle radius exceeds cylinder length.\n");
		if (params.dFreeDev >= params.dFreeAvg)
			fprintf(stderr,"WARNING: large free-particle radius uncertainty---infinite loop possible.\n");
		}

	if (params.nStuckTarget == 0 && params.nFreeTarget == 0)
		usage(argv[0]);

	params.nStuckR = params.nStuckZ = 0; /* to satisfy compiler */
	params.dStuckR = params.dStuckZ = 0.0; /* ditto */
	if (params.nStuckTarget != 0) {
		/* currently allow stuck particles to overlap */
		params.nStuckR = (int) (PI*params.dCylRad/params.dStuckAvg);
		assert(params.nStuckR > 0);
		params.nStuckZ = (int) (0.5*params.dCylLen/params.dStuckAvg) + 1;
		params.nStuck = params.nStuckR*params.nStuckZ;
		if (params.nStuckZ > 1)
			params.nStuck = (int) ((double) params.nStuck - (double) params.nStuckR * 0.5 + 0.5);
		if (params.nStuckTarget != -1 && params.dCylLen > 0.0) {
			if (params.nStuckTarget > params.nStuck)
				fprintf(stderr,"WARNING: Cannot achieve target number of stuck particles.\n");
			else if (params.nStuckTarget < params.nStuck) {
				double s = sqrt(TWO_PI*params.dCylRad*params.dCylLen/params.nStuckTarget);
				params.nStuckR = (int) (TWO_PI*params.dCylRad/s);
				params.nStuckZ = (int) (params.dCylLen/s) + 1;
				params.nStuck = params.nStuckR*params.nStuckZ;
				if (params.nStuckZ > 1)
					params.nStuck = (int) ((double) params.nStuck - (double) params.nStuckR/2 + 0.5);
				if (params.nStuck > params.nStuckTarget)
					fprintf(stderr,"WARNING: Number of stuck particles exceeds target.\n");
				}
			}
		params.dStuckR = TWO_PI*params.dCylRad/params.nStuckR;
		params.dStuckZ = params.dCylLen/params.nStuckZ;
		data = (DATA *) malloc((size_t) params.nStuck*sizeof(DATA));
		assert(data != NULL);
		}

	if (params.nFreeTarget != 0) {
		params.nFreeXY = (int) (sqrt(0.5)*params.dCylRad/(params.dFreeAvg + params.dCutoff*params.dFreeDev));
		assert(params.nFreeXY > 0);
		params.nFreeZ = (int) (0.5*params.dCylLen/(params.dFreeAvg + params.dCutoff*params.dFreeDev));
		assert(params.nFreeZ > 0);
		params.nFree = SQ(params.nFreeXY)*params.nFreeZ;
		if (params.nFreeTarget != -1 && params.dCylLen > 0.0) {
			if (params.nFreeTarget > params.nFree)
				fprintf(stderr,"WARNING: Cannot achieve target number of free particles.\n");
			else if (params.nFreeTarget < params.nFree) {
				double s = pow(2.0*SQ(params.dCylRad)*params.dCylLen/params.nFreeTarget,1.0/3.0);
				params.nFreeXY = (int) (sqrt(2.0)*params.dCylRad/s);
				params.nFreeZ = (int) (params.dCylLen/s) + 1;
				params.nFree = SQ(params.nFreeXY)*params.nFreeZ;
				}
			}
		params.dFreeXY = sqrt(2.0)*params.dCylRad/params.nFreeXY;
		params.dFreeZ = params.dCylLen/params.nFreeZ;
		if (params.dFreeXY <= 2.0*(params.dFreeAvg + params.dCutoff*params.dFreeDev) ||
			(params.nFreeZ > 1 && params.dFreeZ <= 2.0*(params.dFreeAvg + params.dCutoff*params.dFreeDev))) {
			fprintf(stderr,"Cells for free particles too small.\n");
			return 1;
			}

		data = (DATA *) realloc(data,(size_t) (params.nStuck + params.nFree)*sizeof(DATA));
		assert(data != NULL);
		}

	params.nPart = params.nStuck + params.nFree;

	generate(&params,data);

	write_log(&params,data); /* NOTE: reorders particles! */
	printf("Log written to %s.\n",LOGFILENAME);

	write_data(&params,data);

	free(data);

	return 0;
	}

/* sstumbler.c */
