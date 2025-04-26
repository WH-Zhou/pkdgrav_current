/*
 ** sscube.c -- DCR 6/8/12
 ** ========
 ** Generates initial conditions in cubic close pack.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <time.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <random.h>

#define OUTFILENAME "sscube.ss"
#define LOGFILENAME "sscube.log"

typedef struct {                
  int nPartX,nPartY,nPartZ,nPart,bRandPos;
  double dRadius,dDensity,dCellFac,dVelDev;
  } PARAMS;

static void write_log(const PARAMS *p)
{
  FILE *fp;

  fp = fopen(LOGFILENAME,"w");
  if (fp == NULL) {
    fprintf(stderr,"Unable to open \"%s\" for writing.\n",LOGFILENAME);
    exit(1);
  }

  fprintf(fp,"Number of particles in x = %i\n",p->nPartX);
  fprintf(fp,"Number of particles in y = %i\n",p->nPartY);
  fprintf(fp,"Number of particles in z = %i\n",p->nPartZ);
  fprintf(fp,"Number of particles = %i\n",p->nPart);
  fprintf(fp,"Randomize position = %s\n",p->bRandPos ? "yes" : "no");
  fprintf(fp,"Particle radius = %g cm\n",100.0*p->dRadius*AU);
  fprintf(fp,"Particle density = %g g/cc\n",1.0e-3*p->dDensity*M_SUN/(AU*AU*AU));
  fprintf(fp,"Cell multiplier factor = %g\n",p->dCellFac);
  fprintf(fp,"Velocity component deviate = %g cm/s\n",p->dVelDev*(JUL_YR/(TWO_PI*AU)));
  fclose(fp);
}

static void write_data(const PARAMS *p,SSDATA *d)
{
  SSIO ssio;	
  SSHEAD head;
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
    if (ssioData(&ssio,&d[i])) {
      fprintf(stderr,"Error writing data (particle %i).\n",i);
      ssioClose(&ssio);
      exit(1);
    }
  }

  ssioClose(&ssio);
}

/* not used
static int compar(const void *v1,const void *v2)
{
  if (((const SSDATA *)v1)->radius < ((const SSDATA *)v2)->radius)
    return -1;
  else if (((const SSDATA *)v1)->radius > ((const SSDATA *)v2)->radius)
    return 1;
  else
    return 0;
}
*/

static void generate(const PARAMS *p,SSDATA *d)
{
  double dLengthX,dLengthY,dLengthZ,dGap;
  int i,j,k,n;

  randSeedGenerator((Ullong) (time(NULL) % getpid() + getppid()));

  dLengthX = 2.0*p->nPartX*p->dRadius*p->dCellFac;
  dLengthY = 2.0*p->nPartY*p->dRadius*p->dCellFac;
  dLengthZ = 2.0*p->nPartZ*p->dRadius*p->dCellFac;

  dGap = (p->dCellFac - 1.0)*p->dRadius;

  if (!p->bRandPos || dGap < 0.0)
    dGap = 0.0;

  for (i=0;i<p->nPartX;i++)
    for (j=0;j<p->nPartY;j++)
      for (k=0;k<p->nPartZ;k++) {
	n = (i*p->nPartY + j)*p->nPartZ + k;
	d[n].mass = 4.0/3.0*M_PI*p->dRadius*p->dRadius*p->dRadius*p->dDensity;
	d[n].radius = p->dRadius;
	d[n].pos[0] = -0.5*dLengthX + (2.0*i + 1.0)*p->dRadius*p->dCellFac + (2.0*randUniform() - 1.0)*dGap;
	d[n].pos[1] = -0.5*dLengthY + (2.0*j + 1.0)*p->dRadius*p->dCellFac + (2.0*randUniform() - 1.0)*dGap;
	d[n].pos[2] = -0.5*dLengthZ + (2.0*k + 1.0)*p->dRadius*p->dCellFac + (2.0*randUniform() - 1.0)*dGap;
	d[n].vel[0] = p->dVelDev*randGaussian();
	d[n].vel[1] = p->dVelDev*randGaussian();
	d[n].vel[2] = p->dVelDev*randGaussian();
	d[n].spin[0] = d[n].spin[1] = d[n].spin[2] = 0.0;
	d[n].color = PLANETESIMAL;
	d[n].org_idx = n;
      }

  /* sort particles by increasing radius */
  /*qsort(d,p->nPart,sizeof(SSDATA),compar);*/
}

static void usage(const char *achProgName)
{
  fprintf(stderr,
	  "Usage: %s -x nX -y nY -z nZ -r radius -d density [ -c cell-scaling-factor [ -p ] ] [ -v velocity-component-deviate ]\n"
	  "where: nX, nY, nZ are the number of particles in each axis direction\n"
	  "       radius is the particle radius in cm\n"
	  "       density is the particle density in g/cc\n"
	  "       cell-scaling-factor is 1 by default\n"
	  "       -p indicates randomizing position\n"
	  "       and velocity-component-deviate is 0 by default (cm/s).\n",
	  achProgName);

  exit(1);
}

int main(int argc,char *argv[])
{
  PARAMS params;
  SSDATA *data;
  double d;
  int c,n;

  /* defaults */

  params.nPartX = params.nPartY = params.nPartZ = 0;
  params.bRandPos = FALSE;
  params.dRadius = 0.0;
  params.dDensity = 0.0;
  params.dCellFac = 1.0;
  params.dVelDev = 0.0;

  /* parse command-line arguments */

  while ((c = getopt(argc,argv,"x:y:z:r:d:c:pv:")) != EOF)
    switch (c) {
    case 'x':
      n = atoi(optarg);
      if (n <= 0)
	usage(argv[0]);
      params.nPartX = n;
      break;
    case 'y':
      n = atoi(optarg);
      if (n <= 0)
	usage(argv[0]);
      params.nPartY = n;
      break;
    case 'z':
      n = atoi(optarg);
      if (n <= 0)
	usage(argv[0]);
      params.nPartZ = n;
      break;
    case 'r':
      d = atof(optarg);
      if (d <= 0.0)
	usage(argv[0]);
      params.dRadius = 0.01*d/AU; /* cm -> AU */
      break;
    case 'd':
      d = atof(optarg);
      if (d <= 0.0)
	usage(argv[0]);
      params.dDensity = 1.0e3*d/M_SUN*AU*AU*AU; /* g/cc -> M_SUN/AU^3 */
      break;
    case 'c':
      d = atof(optarg);
      if (d <= 0.0)
	usage(argv[0]);
      params.dCellFac = d;
      break;
    case 'p':
      params.bRandPos = TRUE;
      break;
    case 'v':
      d = atof(optarg);
      if (d < 0.0)
	usage(argv[0]);
      params.dVelDev = 0.01*d/(TWO_PI*AU/JUL_YR); /* cm/s -> 2pi AU/YR */
      break;
    case '?':
    default:
      usage(argv[0]);
    }

  params.nPart = params.nPartX*params.nPartY*params.nPartZ;

  if (argc > optind + 1 || params.nPart == 0 || params.dRadius == 0.0 || params.dDensity == 0.0)
    usage(argv[0]);

  data = (SSDATA *) malloc((size_t) params.nPart*sizeof(SSDATA));
  assert(data != NULL);

  printf("Generating data...\n");
  generate(&params,data);

  write_log(&params);
  printf("Log written to %s.\n",LOGFILENAME);

  printf("Writing data...\n");
  write_data(&params,data);

  free(data);

  printf("Now run rpx on %s to reset center of mass, etc.\n",OUTFILENAME);

  return 0;
}

/* sscube.c */
