/*
 ** patcha.c -- DCR 99-01-22
 ** ========
 ** Analyzes orbiting patch model output.
 */

#include <ss.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() & getpid() */
#include <string.h>
#include <math.h>
#include <assert.h>
#include <boolean.h>
#include <rdpar.h>
#include <vector.h>

#define LOG_FILE "patchic.log"
#define OUT_FILE "patcha.out"

#define N_SKEWERS 10000 /* for physical optical depth calculation */
#define N_SAMPLES 10 /* total no. skewers = N_SAMPLES X N_SKEWERS */

/* Tree code definitions */

/* (possible improvement: use buckets) */

#define CELLS_PER_NODE 4 /* Barnes & Hut quad-tree */

struct node_s {
	VECTOR pos;
	double half_size,eff_half_size;
	struct node_s *cell[CELLS_PER_NODE];
	const SSDATA *leaf[CELLS_PER_NODE];
	};

typedef struct node_s NODE;

static int
skewered(const NODE *node,const VECTOR pos)
{
	/*
	 ** Returns 1 if a vertical skewer at pos intercepts a particle in node; 0 otherwise.
	 */

	/* note: checks EACH node cell in case particle not wholly contained in cell */

	int i;

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i]) {
			if (fabs(pos[X] - node->cell[i]->pos[X]) < node->eff_half_size &&
				fabs(pos[Y] - node->cell[i]->pos[Y]) < node->eff_half_size)
				if (skewered(node->cell[i],pos))
					return 1;
			}
		else if (node->leaf[i]) { /* project onto x-y plane */
			if (sqrt(SQ(pos[X] - node->leaf[i]->pos[X]) + SQ(pos[Y] - node->leaf[i]->pos[Y])) <= node->leaf[i]->radius)
				return 1;
			}
	return 0;
	}

static void
make_node(const VECTOR pos,double size,NODE **node)
{
	int i;

	assert(size > 0.0);
	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	COPY_VEC(pos,(*node)->pos);
	(*node)->half_size = (*node)->eff_half_size = 0.5*size;

	for (i=0;i<CELLS_PER_NODE;i++) {
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}
	}

static void
check_eff_size(NODE *node,const SSDATA *d)
{
	double r;

	/* use Manhattan metric because cells are checked on rectangular grid in skewered() */
	r = fabs(d->pos[X] - node->pos[X]) + fabs(d->pos[Y] - node->pos[Y]) + d->radius;
	if (r > node->eff_half_size)
		node->eff_half_size = r;
	}

static void
add_to_tree(NODE *node,const SSDATA *d)
{
	/* note: assumes particle inside node! */

	int i,idx,idy;

	idx = (d->pos[X] < node->pos[X] ? -1 : 1);
	idy = (d->pos[Y] < node->pos[Y] ? -1 : 1);

	i = (idx + 1)/2 + idy + 1; /* 2-D tree */

	if (node->cell[i])
		add_to_tree(node->cell[i],d);
	else if (node->leaf[i]) {
		VECTOR v;
		assert(node->leaf[i]->pos[X] != d->pos[X] || node->leaf[i]->pos[Y] != d->pos[Y]); /* otherwise would get infinite loop */
		SET_VEC(v,idx,idy,0.0);
		SCALE_VEC(v,0.5*node->half_size);
		ADD_VEC(v,node->pos,v);
		make_node(v,node->half_size,&node->cell[i]);
		add_to_tree(node->cell[i],node->leaf[i]);
		add_to_tree(node->cell[i],d);
		node->leaf[i] = NULL; /* for completeness */
		}
	else
		node->leaf[i] = d;

	/* update the current cell's effective size in all cases */
	
	check_eff_size(node,d);
	}

static void
kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}

static void
analyze(FILE *summary,double omega,double lx,double ly,double r,
		double time,int n_data,const SSDATA *data)
{
	const SSDATA *d;
	NODE *root = NULL;
	VECTOR v,com_pos,com_vel,vel_disp;
	double total_mass,ff,zh,r2,z2,x,y,tau_dyn,tau_phys[N_SAMPLES],tau_phys_avg,tau_phys_std,xsig;
	int i,j,k,nOut;

	total_mass = ff = 0;

	ZERO_VEC(com_pos);
	ZERO_VEC(com_vel);

	for (i=0;i<n_data;i++) {
		d = &data[i];
		total_mass += d->mass;
		COPY_VEC(d->pos,v);
		SCALE_VEC(v,d->mass);
		ADD_VEC(com_pos,v,com_pos);
		COPY_VEC(d->vel,v);
		v[Y] += 1.5*omega*d->pos[X]; /* remove shear component */
		SCALE_VEC(v,d->mass);
		ADD_VEC(com_vel,v,com_vel);
		r2 = SQ(d->radius);
		z2 = SQ(d->pos[Z]);
		if (z2 <= r2) ff += r2 - z2;
		}

	NORM_VEC(com_pos,total_mass);
	NORM_VEC(com_vel,total_mass);
	ff *= PI/(lx*ly);

	ZERO_VEC(vel_disp);
	zh = xsig = 0;

	for (i=0;i<n_data;i++) {
		d = &data[i];
		COPY_VEC(d->vel,v);
		v[Y] += 1.5*omega*d->pos[X];
		SUB_VEC(v,com_vel,v);
		for (k=0;k<N_DIM;k++)
			vel_disp[k] += d->mass*SQ(v[k]);
		zh += d->mass*SQ(d->pos[Z] - com_pos[Z]);
		xsig += SQ(d->pos[X] - com_pos[X]); /* only useful if not periodic in x (radial direction) */
		}

	NORM_VEC(vel_disp,total_mass);

	for (k=0;k<N_DIM;k++)
		vel_disp[k] = sqrt(vel_disp[k]);

	zh = sqrt(zh/total_mass);
	xsig = sqrt(xsig/n_data);

	/* Optical depth computations: start by building quad tree */

	/*assert(ly >= lx);*/ /* for now, make ly*ly tree and compute with central lx*lx square */
	(void) printf("Building tree...\n");
	ZERO_VEC(v); /* tree center coincides with patch center */
	make_node(v,(lx > ly ? lx : ly),&root);
	tau_dyn = 0.0;
	nOut = 0;
	for (i=0;i<n_data;i++) {
		if (fabs(data[i].pos[X]) > root->half_size ||
			fabs(data[i].pos[Y]) > root->half_size) {
			++nOut;
			continue;
			}
		add_to_tree(root,&data[i]);
		/*
		 ** Dynamical optical depth calculation assumes periodic boundary conditions,
		 ** i.e., any portion of a particle that extends beyond the patch boundary is
		 ** wrapped around to the other side.
		 */
		tau_dyn += SQ(data[i].radius);
		}
	if (nOut > 0) {
		(void) fprintf(stderr,"WARNING: %i particle%s (%.2f%%) outside patch!\n",nOut,nOut==1?"":"s",100*(double)nOut/n_data);
		if (nOut == n_data) {
			(void) fprintf(stderr,"All particles outside patch!  Skipping optical depth computation.\n");
			tau_phys_avg = tau_phys_std = 0.0;
			goto finish;
			}
		}
	tau_dyn *= M_PI/(lx*ly);
	(void) printf("Estimating physical optical depth...\n");
	srand48(getpid()); /* seed random number generator */
	/* shoot random skewers into patch */
	assert(N_SAMPLES > 1);
	tau_phys_avg = 0.0;
	for (i=0;i<N_SAMPLES;i++) {
		for (k=j=0;j<N_SKEWERS;j++) {
			x = (drand48() - 0.5)*lx;
			y = (drand48() - 0.5)*ly;
			SET_VEC(v,x,y,0.0);
			if (!skewered(root,v))
				++k;
			}
		tau_phys[i] = -log((double)k/N_SKEWERS); /* k = 0 ==> tau_phys = INF */
		tau_phys_avg += tau_phys[i];
		}
	tau_phys_avg /= N_SAMPLES;
	tau_phys_std = 0.0;
	for (i=0;i<N_SAMPLES;i++)
		tau_phys_std += SQ(tau_phys[i] - tau_phys_avg);
	tau_phys_std = sqrt(tau_phys_std/(N_SAMPLES - 1));
	(void) printf("Done...tau_phys = %g +/- %g (tau_dyn = %g)\n",
				  tau_phys_avg,tau_phys_std,tau_dyn);
	kill_node(root);

 finish:

	/* Scale units */

	time *= JUL_YR*omega/(TWO_PI*T_SCALE); /* in orbital periods */
	NORM_VEC(com_pos,sqrt(lx*ly));
	NORM_VEC(com_vel,omega*sqrt(lx*ly));
	NORM_VEC(vel_disp,omega*r);
	zh /= r;
	xsig *= L_SCALE/1000; /* in km */

	/* Output */

	(void) fprintf(summary,
				   "%e "
				   "%e %e %e "
				   "%e %e %e "
				   "%e %e %e "
				   "%e %e %e %e %e %e\n",
				   time,
				   com_pos[X],com_pos[Y],com_pos[Z],
				   com_vel[X],com_vel[Y],com_vel[Z],
				   vel_disp[X],vel_disp[Y],vel_disp[Z],
				   tau_dyn,tau_phys_avg,tau_phys_std,ff,zh,xsig);
	}

static int
read_data(const char *filename,double *time,int *n_data,SSDATA **data)
{
	SSIO ssio;
	SSHEAD h;
	int i;

	*time = 0;
	*n_data = 0;
	*data = NULL;

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"Corrupt header.\n");
		return 1;
		}

	*time = h.time;
	*n_data = h.n_data;

	*time *= (T_SCALE/JUL_YR); /* get time in Julian years */

	(void) printf("time %e Julian yr, n_data = %i\n",*time,*n_data);

	if (*n_data <= 0) {
		(void) fprintf(stderr,"Missing or corrupt data.\n");
		return 1;
		}

    switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		ssioClose(&ssio);
		return 1;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}

	*data = (SSDATA *)malloc((size_t) *n_data*sizeof(SSDATA));
	assert(*data != NULL);

	for (i=0;i<*n_data;i++)
		if (ssioData(&ssio,&(*data)[i])) {
			(void) fprintf(stderr,"Corrupt data.\n");
			return 1;
			}

	(void) ssioClose(&ssio);

	return 0;
	}

static void
get_params(double *omega,double *lx,double *ly,double *r)
{
	OpenPar(LOG_FILE);
	ReadDbl("Orbital frequency",omega);
	ReadDbl("Patch width",lx);
	ReadDbl("Patch length",ly);
	ReadDbl("Mean particle radius",r);
	ClosePar();
	}

int
main(int argc,char *argv[])
{
	SSDATA *data;
	FILE *fp = NULL;
	double omega,lx,ly,r,time;
	BOOLEAN usage=FALSE,append=FALSE,force=FALSE;
	int n_data,c,i;

	setbuf(stdout,(char *)NULL);

	while (!usage && ((c = getopt(argc,argv,"af")) != EOF))
		switch (c) {
		case 'a':
			if (force || append) usage = TRUE;
			else if (!(fp = fopen(OUT_FILE,"a"))) {
				(void) fprintf(stderr,"Unable to open %s for append.\n",
							   OUT_FILE);
				return 1;
				}
			else append = TRUE;
			break;
		case 'f':
			if (force || append) usage = TRUE;
			else force = TRUE;
			break;
		default:
			usage = TRUE;
			}

	if (optind == argc) usage = TRUE;

	if (usage) {
		(void) fprintf(stderr,"Usage: %s [ -a | -f ] ss-file "
					   "[ ss-file ... ]\n",argv[0]);
		return 1;
		}

	if (!append && !force && (fp = fopen(OUT_FILE,"r"))) {
		(void) fprintf(stderr,"Output file %s already exists.\n",OUT_FILE);
		(void) fclose(fp);
		return 1;
		}

	if (!append && !(fp = fopen(OUT_FILE,"w"))) {
		(void) fprintf(stderr,"Unable to open %s for writing.\n",OUT_FILE);
		return 1;
		}

	get_params(&omega,&lx,&ly,&r);

	for (i=optind;i<argc;i++) {
		(void) printf("%s:\n",argv[i]);
		if (!read_data(argv[i],&time,&n_data,&data))
			analyze(fp,omega,lx,ly,r,time,n_data,data);
		if (data) free((void *)data);
		}

	(void) fclose(fp);

	return 0;
	}

/* patcha.c */
