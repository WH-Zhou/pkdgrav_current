/*
 ** ranics.c -- DCR 1/13/04
 ** ========
 ** ss utility to generate uniform random ICs.
 ** (based on old box.c).
 **
 ** UNITS: standard ss units (G=1,etc.).
 */

#include <ss.h>
#include <stdlib.h>
#include <unistd.h> /* for getpid() */
#include <math.h>
#include <assert.h>
#include <rdpar.h>
#include <vector.h>

#define PARAM_FILE	"ranics.par"
#define OUTPUT_FILE	"ranics.ss"

#define RAN ((double) rand()/RAND_MAX)

typedef struct {
	int N;
	double m,R,Rd,KEf,radf;
	} PARAMS;

static void gen_body(PARAMS *p,SSDATA *d,int iOrgIdx)
{
	double phi,theta,v;
	int k;

	d->mass = p->m;
	d->radius = p->R;
	for (k=0;k<N_DIM;k++) {
		do {
			SET_VEC(d->pos,(2*RAN - 1)*p->Rd,(2*RAN - 1)*p->Rd,(2*RAN - 1)*p->Rd);
			} while (MAG(d->pos) + p->R > p->Rd);
/* use following for more centrally condensed distribution...
		r = RAN*(p->Rd - p->R);
		phi = TWO_PI*RAN;
		theta = PI*RAN;
		SET_VEC(d->pos,sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
		SCALE_VEC(d->pos,r);
*/
		v = RAN; /* will be rescaled later */
		phi = TWO_PI*RAN;
		theta = PI*RAN;
		SET_VEC(d->vel,sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta));
		SCALE_VEC(d->vel,v); /* isotropic initially */
		}
	ZERO_VEC(d->spin);
	d->color = PLANETESIMAL;
	d->org_idx = iOrgIdx;
	}

static double pot(PARAMS *p,SSDATA *d,int i)
{
	/* computes potential at particle i position */

	double pot=0.0,r;
	int j;

	for (j=0;j<p->N;j++) {
		if (j == i) continue;
		r = DIST(d[i].pos,d[j].pos);
		assert(r > 0.0); /* particles should not be at same spot */
		pot -= d[j].mass/r;
		}

	return pot;
	}

static void reset_vel(PARAMS *p,double f,SSDATA *d)
{
	VECTOR r_hat,vr_vec,t_hat,vt_vec;
	double v,vr,vt;

	/* set speed scaled by desired virial fraction */

	SCALE_VEC(d->vel,sqrt(f));

	/* adjust radial and tangential components as desired */

	if (p->radf < 0) return;

	COPY_VEC(d->pos,r_hat);
	NORM_VEC(r_hat,MAG(r_hat)); /* radial unit vector */

	vr = DOT(d->vel,r_hat);
	COPY_VEC(r_hat,vr_vec);
	SCALE_VEC(vr_vec,vr); /* current radial component of vel */

	SUB_VEC(d->vel,vr_vec,t_hat);
	NORM_VEC(t_hat,MAG(t_hat)); /* tangential unit vector */

	v = MAG(d->vel);
	vr = sqrt(p->radf)*v;
	vt = sqrt(1 - p->radf)*v;

	COPY_VEC(r_hat,vr_vec);
	SCALE_VEC(vr_vec,vr); /* new radial component */

	COPY_VEC(t_hat,vt_vec);
	SCALE_VEC(vt_vec,vt); /* new tangential component */

	ADD_VEC(vr_vec,vt_vec,d->vel);
	}

static void generate(PARAMS *p,SSDATA *d)
{
	double GPE,KE,f;
	int i,j;

	for (i=0;i<p->N;i++) {
		gen_body(p,&d[i],i);
		for (j=0;j<i;j++)
			if (OVERLAP(d[i].pos,d[i].radius,d[j].pos,d[j].radius)) {
				(void) printf("Particle %i overlaps particle %i\n"
							  "-- regenerating particle %i\n",i,j,i);
				--i;
				break;
				}
		}

	GPE = KE = 0.0;
	for (i=0;i<p->N;i++) {
		GPE += d[i].mass*pot(p,d,i);
		KE += 0.5*d[i].mass*MAG_SQ(d[i].vel);
		}

	(void) printf("Starting KE = %g = %g times GPE\n",KE,KE/GPE);

	assert(KE > 0.0);

	f = -0.5*p->KEf*GPE/KE;

	assert(f >= 0.0);

	for (i=0;i<p->N;i++)
		reset_vel(p,f,&d[i]);

	GPE = KE = 0.0;
	for (i=0;i<p->N;i++) {
		GPE += d[i].mass*pot(p,d,i);
		KE += 0.5*d[i].mass*MAG_SQ(d[i].vel);
		}

	(void) printf("Adjusted KE = %g = %g times GPE (= %g times virial)\n",KE,KE/GPE,-2*KE/GPE);

	{
	double t_dyn,t_step;

	t_dyn = 2*sqrt(CUBE(p->Rd)/(p->N*p->m));
	(void) printf("Dynamical time ~ %g\n",t_dyn);
	t_step = 2*sqrt(CUBE(p->R)/p->m)/33;
	(void) printf("Recommended time step < %g\n",t_step);
	(void) printf("Estimated number of steps needed per t_dyn = %i\n",(int) (t_dyn/t_step + 0.5));
	}
	}

static void output_data(PARAMS *p,SSDATA *d)
{
	SSIO ssio;
	SSHEAD head;
	int i;

	if (ssioOpen(OUTPUT_FILE,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",OUTPUT_FILE);
		exit(1);
		}

	head.time = 0.0;
	head.n_data = p->N;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	(void) ssioHead(&ssio,&head);

	for (i=0;i<p->N;i++)
		(void) ssioData(&ssio,&d[i]);

	(void) ssioClose(&ssio);
	}

static void fix_rejects(PARAMS *p)
{
	FILE *fp;
	SSIO ssio;
	SSDATA data;
	int sh,sd;
	int rej,nRej=0,i;

	if (!(fp = fopen(REJECTS_FILE,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",REJECTS_FILE);
		exit(1);
		}

	(void) fscanf(fp,"%i",&i);
	if (i != p->N) {
		(void) fprintf(stderr,"Incompatible rejects file.\n");
		exit(1);
		}

	if (ssioOpen(OUTPUT_FILE,&ssio,SSIO_UPDATE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",OUTPUT_FILE);
		exit(1);
		}

	sh = SSHEAD_SIZE; /* NOT sizeof(SSHEAD)! XDR has special data sizes */
	sd = SSDATA_SIZE; /* NOT sizeof(SSDATA)! Ditto */

	(void) ssioSetPos(&ssio,(size_t) sh);

	i = 0;

	while (fscanf(fp,"%i",&rej) != EOF) {
		if (rej >= 0) {
			if (i != rej) {
				(void) fprintf(stderr,"Corrupted rejects data.\n");
				exit(1);
				}
			gen_body(p,&data,i);
			(void) ssioSetPos(&ssio,(size_t) (sh + i*sd));
			(void) ssioData(&ssio,&data);
			(void) printf("Particle %i regenerated.\n",rej);
			++nRej;
			}
		++i;
		}

	(void) ssioClose(&ssio);
	(void) fclose(fp);

	(void) printf("%i particle%s regenerated.\n",nRej,nRej==1?"":"s");
	}

static void get_params(PARAMS *params)
{
	OpenPar(PARAM_FILE);

	ReadInt("Number of particles",&params->N);
	ReadDbl("Particle mass",&params->m);
	ReadDbl("Particle radius",&params->R);
	ReadDbl("Domain radius",&params->Rd);
	ReadDbl("KE fraction of virial",&params->KEf);
	ReadDbl("Radial fraction of speed",&params->radf);

	ClosePar();

	assert(params->N > 0);
	assert(params->m > 0.0);
	assert(params->R > 0.0);
	assert(params->Rd >= params->R);
	assert(params->KEf >= 0.0);
	assert(params->radf <= 1.0);
	}

int main(int argc,char *argv[])
{
	PARAMS params;
	SSDATA *data = NULL;
	long int pid;
	int force=0,rejects=0,usage=0,c;

	/* Disable stdout buffering */

	setbuf(stdout,(char *) NULL);

	/* Check arguments */

	while ((c = getopt(argc,argv,"fr")) != EOF)
		switch (c) {
		case 'f':
			force = 1;
			break;
		case 'r':
			rejects = 1;
			break;
		default:
			usage = 1;
			}

	if (usage || optind < argc) {
		(void) fprintf(stderr,"Usage: %s [ -f ] [ -r ]\n",argv[0]);
		exit(1);
		}

	/* Get model parameters */

	get_params(&params);

	if (!force && !rejects) {
		FILE *fp = fopen(OUTPUT_FILE,"r");
		if (fp) {
			(void) fprintf(stderr,"%s NOT OVERWRITTEN.\n",OUTPUT_FILE);
			(void) fclose(fp);
			return 1;
			}
		}

	/* Generate initial conditions */

	pid = getpid();
	(void) printf("Random number seed = %li\n",pid);
	srand((unsigned int) pid);

	if (rejects) {
		assert(0); /* not supported currently: do it the N^2 way for now */
		fix_rejects(&params);
		return 0;
		}

	if (!(data = (SSDATA *) malloc((size_t) params.N*sizeof(SSDATA)))) {
		(void) fprintf(stderr,"Unable to allocate data memory\n");
		return 1;
		}

	generate(&params,data);

	/* Save data */

	output_data(&params,data);

	/* All done */

	free((void *) data);

	return 0;
	}

/* ranics.c */
