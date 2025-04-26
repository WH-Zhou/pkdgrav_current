/*
 ** ssa.c -- DCR 97-08-14
 ** =====
 ** Analyzes Solar System particle data output.
 */

#include <ss.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() */
#include <math.h>
#include <assert.h>
#include <boolean.h>
#include <vector.h>
#include <coords.h>

#ifdef sparc
#ifdef sun
#undef sun /* yikes! */
#endif
#endif

/*#define SSTEST*/

/*#define AVOID_RESONANCES*/

#define OUTFILE "ssa.out"

typedef struct {
	double a,e,i,w,W,vdr2,vd2;
	} ELEM_T;

int no_mass_wgt = 0; /*DEBUG hack for TEST particles*/

#define MASS(p) (no_mass_wgt ? 1.0 : (p)->mass)

static int read_data(char *filename,double *time,int *n_data,SSDATA **data,SSDATA **sun,
		  SSDATA **jup,SSDATA **sat,SSDATA **ura,SSDATA **nep)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA *d;
	int i,n,n_planets,nt;

	*n_data = n = n_planets = nt = 0;
	*data = *sun = *jup = *sat = *ura = *nep = NULL;

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",filename);
		return 1;
		}

	(void) ssioHead(&ssio,&h);

	*time = h.time;
	*n_data = h.n_data;

	*time *= (T_SCALE/JUL_YR); /* get time in Julian years */

	(void) printf("time %e Julian yr, n_data = %i -- ",*time,*n_data);

	if (*n_data <= 0) {
		(void) fprintf(stderr,"Missing or corrupted data -- skipped\n");
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

	*data = (SSDATA *) malloc((size_t) *n_data*sizeof(SSDATA));
	assert(*data != NULL);

	for (i=0;i<*n_data;++i) {
		d = &(*data)[i];
		(void) ssioData(&ssio,d);
		switch (d->color) {
		case SUN:
			assert(*sun == NULL);
			*sun = d;
			break;
		case JUPITER:
			assert(*jup == NULL);
			*jup = d;
			++n_planets;
			break;
		case SATURN:
			assert(*sat == NULL);
			*sat = d;
			++n_planets;
			break;
		case URANUS:
			assert(*ura == NULL);
			*ura = d;
			++n_planets;
			break;
		case NEPTUNE:
			assert(*nep == NULL);
			*nep = d;
			++n_planets;
			break;
		case PLANETESIMAL:
			++n;
			break;
		case FRAGMENT:
			++n;
			break;
		case BINARIES:
			++n;
			break;
		case TEST:
			++nt;
			if (nt == 1) {
				(void) fprintf(stderr,"WARNING: test particle detected -- "
							   "mass weighting disabled.\n");
				no_mass_wgt = 1;
				}
			d->mass = 0;
			break;
		default:
			(void) fprintf(stderr,"Unknown object (item %i, color %i) in "
						   "data file! Aborting...\n",i,d->color);
			return 1;
			}
		}

	(void) printf("%s%i planet%s, %i planetesimal%s\n",(*sun?"Sun,":""),
				  n_planets,(n_planets==1?"":"s"),n,(n==1?"":"s"));

	if (nt) (void) printf("%i test particle%s\n",nt,(nt==1?"":"s"));

	(void) ssioClose(&ssio);

	if (!*sun)
		(void) fprintf(stderr,"Sun not found...assuming heliocentric coords...\n");

	return 0;
	}

typedef struct {
	int idx,color;
	double mass;
	} SORT_ELEM;

static int
sort_compar(const void *ev1,const void *ev2)
{
	const SORT_ELEM *e1 = (const SORT_ELEM *) ev1;
	const SORT_ELEM *e2 = (const SORT_ELEM *) ev2;
	
	double m1,m2;
	int c1,c2;

	c1 = e1->color;
	c2 = e2->color;
	m1 = e1->mass;
	m2 = e2->mass;

	if ((c1 == PLANETESIMAL && c2 == PLANETESIMAL) ||
		(c1 != PLANETESIMAL && c2 != PLANETESIMAL))
		return (m1 < m2 ? 1 : (m1 == m2 ? 0 : -1));

	return (c1 != PLANETESIMAL ? 1 : -1);
	}

static void
do_snap(int n_data,SSDATA *data,ELEM_T *elem)
{
	SORT_ELEM *base;
	SSDATA *d;
	ELEM_T *e;
	FILE *fp = NULL;
	int i;

	base = (SORT_ELEM *) malloc((size_t) n_data*sizeof(SORT_ELEM));
	assert(base != NULL);

	for (i=0;i<n_data;i++) {
		base[i].idx = i;
		base[i].color = data[i].color;
		base[i].mass = data[i].mass;
		}

	/* Sort in descending order of mass, planetesimals first */

	qsort(base,(size_t) n_data,sizeof(SORT_ELEM),sort_compar);

	if (!(fp = fopen("ssa.sss","w"))) {
		(void) fprintf(stderr,"do_snap(): Unable to open sss file\n");
		return;
		}

	for (i=0;i<n_data;i++) {
		d = &data[base[i].idx];
		e = &elem[base[i].idx];
		(void) fprintf(fp, "%i %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e\n",
					   d->color,d->mass,d->radius,e->a,e->e,e->i,
					   d->spin[X],d->spin[Y],d->spin[Z]);
		}

	(void) fclose(fp);

	free((void *)base);
	}

static void
get_vel_disp(double central_mass,double a,VECTOR pos,VECTOR vel,
			 double *vd2,double *vdr2)
{
	/* for vel disp, subtract from equivalent circular orbit */

	/*DEBUG assumes central_mass >> object mass!*/

	double omega,vdr;
	VECTOR v,dv;

	omega = sqrt(central_mass/CUBE(a));
	SET_VEC(v,-pos[Y],pos[X],0);
	SCALE_VEC(v,omega);
	SUB_VEC(vel,v,dv);

	*vd2 = MAG_SQ(dv);

	/* get radial component in plane */

	SET_VEC(v,pos[X],pos[Y],0);
	vdr = DOT(dv,v)/MAG(v);
	*vdr2 = SQ(vdr);
	}

#define N_ABINS 100
#define MAX_N_MBINS 100

static void
do_bins(int n_data,SSDATA *data,ELEM_T *elem,double a_inner,double a_outer,
		double m_min,double m_max)
{
	SSDATA *d;
	ELEM_T *e;
	FILE *fp;
	double aw,mw,mb[N_ABINS],vdb[N_ABINS],vdrb[N_ABINS],*rmseb,*rmsib;
	int i,j,n_mbins,nab[N_ABINS],*nmb;

	if (!(fp = fopen("ssa.abin","w"))) {
		(void) fprintf(stderr,"Unable to open abin file\n");
		return;
		}

	aw = (a_outer - a_inner)/(N_ABINS - 1);

	for (i=0;i<N_ABINS;i++) {
		nab[i] = 0;
		mb[i] = vdb[i] = vdrb[i] = 0.0;
		}

	for (i=0;i<n_data;i++) {
		d = &data[i];
		if (d->color != PLANETESIMAL && d->color != TEST) continue;
		e = &elem[i];
		if (e->a < a_inner || e->a > a_outer) continue;
		j = (aw == 0 ? 0 : (int) (0.5 + (e->a - a_inner)/aw));
		++nab[j];
		mb[j] += MASS(d);
		vdb[j] += MASS(d)*e->vd2;
		vdrb[j] += MASS(d)*e->vdr2;
		}

	for (i=0;i<N_ABINS;i++) {
		vdb[i]  = (mb[i] == 0 ? 0 : sqrt( vdb[i]/mb[i]));
		vdrb[i] = (mb[i] == 0 ? 0 : sqrt(vdrb[i]/mb[i]));
		}

	for (i=0;i<N_ABINS;i++) {
		(void) fprintf(fp,"%e %e %e %e\n",a_inner + i*aw,mb[i],vdb[i],vdrb[i]);
		if (aw == 0) break;
		}

	(void) fclose(fp);

	if (!(fp = fopen("ssa.mbin","w"))) {
		(void) fprintf(stderr,"Unable to open mbin file\n");
		return;
		}

	n_mbins = (int) (m_max/m_min + 0.5); /* 0.5 to fix truncation */
	mw = m_min;
	if (n_mbins > MAX_N_MBINS) {
		n_mbins = MAX_N_MBINS;
		mw = (m_max - m_min)/(n_mbins - 1);
		}

	rmseb = (double *) malloc((size_t) n_mbins*sizeof(double));
	rmsib = (double *) malloc((size_t) n_mbins*sizeof(double));
	nmb = (int *) malloc((size_t) n_mbins*sizeof(int));

	assert(rmseb != NULL && rmsib != NULL && nmb != NULL);

	for (i=0;i<n_mbins;i++) {
		nmb[i] = 0;
		rmseb[i] = rmsib[i] = 0.0;
		}

	for (i=0;i<n_data;i++) {
		d = &data[i];
		if (d->color != PLANETESIMAL && d->color != TEST) continue;
		e = &elem[i];
		j = (int) ((d->mass - m_min)/mw);
		++nmb[j];
		rmseb[j] += SQ(e->e);
		rmsib[j] += SQ(e->i);
		}

	for (i=0;i<n_mbins;i++) {
		rmseb[i] = (nmb[i] == 0 ? 0 : sqrt(rmseb[i]/nmb[i]));
		rmsib[i] = (nmb[i] == 0 ? 0 : sqrt(rmsib[i]/nmb[i]));
		}

	for (i=0;i<n_mbins;i++)
		(void) fprintf(fp,"%e %e %e\n",m_min + i*mw,rmseb[i],rmsib[i]);

	(void) fclose(fp);

	free((void *) nmb);
	free((void *) rmsib);
	free((void *) rmseb);
	}

#undef MAX_N_MBINS
#undef N_ABINS

static void
analyze(FILE *summary,double time,int n_data,SSDATA *data,SSDATA *sun,
		SSDATA *jup,SSDATA *sat,SSDATA *ura,SSDATA *nep,double central_mass,
		BOOLEAN last_file)
{
	struct helio hel;
	struct delaunay del;

	SSDATA *d;
	ELEM_T *elem,*e;
	VECTOR com_pos,com_vel,v;
	double total_mass,a_inner,m_min,a_outer,m_max,a_max,e_max,i_max,
	       spin_max,spinz_max,vd,rmse,rmsi;
	int i,n;

#ifdef AVOID_RESONANCES
	double zone_mass = 0;
#endif

	elem = (ELEM_T *) malloc((size_t) n_data*sizeof(ELEM_T));
	assert(elem != NULL);

	/* Get barycentre position & velocity */

	total_mass = 0;

    ZERO_VEC(com_pos);
    ZERO_VEC(com_vel);

    for (i=0;i<n_data;++i) {
		d = &data[i];
		total_mass += d->mass;
		COPY_VEC(d->pos,v);
		SCALE_VEC(v,d->mass);
		ADD_VEC(com_pos,v,com_pos);
		COPY_VEC(d->vel,v);
		SCALE_VEC(v,d->mass);
		ADD_VEC(com_vel,v,com_vel);
		}

	NORM_VEC(com_pos,total_mass);
	NORM_VEC(com_vel,total_mass);

	a_inner = m_min = HUGE_VAL;

	total_mass = a_outer = m_max = a_max = e_max = i_max = spin_max =
		spinz_max = vd = rmse = rmsi = 0;

	for (i=0;i<n_data;++i) {
		d = &data[i];
		e = &elem[i];
		if (d == sun)
			e->a = e->e = e->i = e->w = 0;
		else {
			if (sun) {
				SUB_VEC(d->pos,sun->pos,d->pos); /* refer to heliocentre */
				SUB_VEC(d->vel,sun->vel,d->vel);
				}
			hel.x = d->pos[X];
			hel.y = d->pos[Y];
			hel.z = d->pos[Z];
			hel.vx = d->vel[X];
			hel.vy = d->vel[Y];
			hel.vz = d->vel[Z];
			heltodel(central_mass + d->mass,&hel,&del);
			e->a = del.sma;
			e->e = del.ecc;
			e->i = del.inc;
			e->w = del.lop;
			e->W = del.lan;
#ifdef SSTEST
			{
			VECTOR h;
			double m,r,v2,h2,sma,ecc,inc;

			m = central_mass + d->mass;
			r = MAG(d->pos);
			v2 = MAG_SQ(d->vel);
			sma = 1/(2/r - v2/m);
			CROSS(d->pos,d->vel,h);
			h2 = MAG_SQ(h);
			ecc = sqrt(fabs(1 - h2/(m*sma))); /* fabs(): no tiny -ve args */
			inc = acos(h[Z]/sqrt(h2));
			(void) printf("Expected %f %f %f got %f %f %f\n",sma,ecc,inc,
						  e->a,e->e,e->i);
			}
#endif /*SSTEST*/
			}
		NORM_VEC(d->spin,365.25); /* convert spin to inverse days */
		if (d->color == PLANETESIMAL || d->color == TEST) {
			total_mass += MASS(d);
			a_inner = (e->a < a_inner ? e->a : a_inner);
			a_outer = (e->a > a_outer ? e->a : a_outer);
			m_min = (d->mass < m_min ? d->mass : m_min);
			if (d->mass > m_max) {
				m_max = d->mass;
				a_max = e->a;
				e_max = e->e;
				i_max = e->i;
				spin_max = MAG(d->spin);
				spinz_max = d->spin[Z];
				}
			get_vel_disp(central_mass,e->a,d->pos,d->vel,&e->vd2,&e->vdr2);
			vd   += MASS(d)*e->vd2;
#ifdef AVOID_RESONANCES
			if (e->a < 2) {
				zone_mass += MASS(d);
#endif
			rmse += MASS(d)*SQ(e->e);
			rmsi += MASS(d)*SQ(e->i);
#ifdef AVOID_RESONANCES
				}
#endif
			}
		}

	if (last_file)
		do_snap(n_data,data,elem);

	/* Binning */

	if (last_file && total_mass > 0)
		do_bins(n_data,data,elem,a_inner,a_outer,m_min,m_max);

	/* Summary */

	if (total_mass > 0) {
		vd = sqrt(vd/total_mass);
#ifdef AVOID_RESONANCES
		rmse = sqrt(rmse/zone_mass);
		rmsi = sqrt(rmsi/zone_mass);
#else
		rmse = sqrt(rmse/total_mass);
		rmsi = sqrt(rmsi/total_mass);
#endif
		(void) printf("rms e = %f, rms i = %f rad\n",rmse,rmsi);
		}

	n = n_data - (sun != NULL) - (jup != NULL) - (sat != NULL) -
		(ura != NULL) - (nep != NULL); /* no. of planetesimals */

	(void) fprintf(summary,"%e %i %e %e %e %e %e %e %e %e %e %e %e %e %e %e "
				   "%e %e",time,n,com_pos[X],com_pos[Y],com_pos[Z],com_vel[X],
				   com_vel[Y],com_vel[Z],(n==0?0.0:total_mass/n),m_max,vd,rmse,
				   rmsi,a_max,e_max,i_max,spin_max,spinz_max);

	for (i=0;i<n_data;++i)
		if (data[i].color == JUPITER) {
			(void) fprintf(summary," %e %e %e %e %e",elem[i].a,elem[i].e,
						   elem[i].i,elem[i].w,elem[i].W);
			break;
			}

	if (i == n_data)
		(void) fprintf(summary," 0 0 0 0 0");

	for (i=0;i<n_data;++i)
		if (data[i].color == SATURN) {
			(void) fprintf(summary," %e %e %e %e %e",elem[i].a,elem[i].e,
						   elem[i].i,elem[i].w,elem[i].W);
			break;
			}

	if (i == n_data)
		(void) fprintf(summary," 0 0 0 0 0");

	for (i=0;i<n_data;++i)
		if (data[i].color == URANUS) {
			(void) fprintf(summary," %e %e %e %e %e",elem[i].a,elem[i].e,
						   elem[i].i,elem[i].w,elem[i].W);
			break;
			}

	if (i == n_data)
		(void) fprintf(summary," 0 0 0 0 0");

	for (i=0;i<n_data;++i)
		if (data[i].color == NEPTUNE) {
			(void) fprintf(summary," %e %e %e %e %e",elem[i].a,elem[i].e,
						   elem[i].i,elem[i].w,elem[i].W);
			break;
			}

	if (i == n_data)
		(void) fprintf(summary," 0 0 0 0 0");

	(void) fprintf(summary,"\n");

	/* Free resources */

	free((void *) elem);
	}

int
main(int argc, char *argv[])
{
	FILE *fp = NULL;
	SSDATA *data,*sun,*jup,*sat,*ura,*nep;
	double time,cm = 0;
	int i,c,usage,append,force,n_data;

#ifdef AVOID_RESONANCES
	(void) fprintf(stderr,"WARNING: rms e & i calculated away from resonances.\n");
#endif

	/* Disable stdout buffering */

	setbuf(stdout,(char *)NULL);

	/* Check arguments */

	usage = append = force = FALSE;

	while (!usage && ((c = getopt(argc,argv,"afm:")) != EOF))
		switch (c) {
		case 'a':
			if (force || append) usage = TRUE;
			else if (!(fp = fopen(OUTFILE,"a"))) {
				(void) fprintf(stderr,"Unable to open %s for append\n",OUTFILE);
				exit(1);
				}
			append = TRUE;
			break;
		case 'f':
			if (force || append) usage = TRUE;
			force = TRUE;
			break;
		case 'm':
			cm = atof(optarg);
			if (cm < 0) usage = TRUE;
			break;
		default:
			usage = TRUE;
			}

	if (optind == argc) usage = TRUE;

	if (usage) {
		(void) fprintf(stderr,"Usage: %s [ -a | -f ] [ -m central_mass ] datafile [ datafile ... ]\n",argv[0]);
		return 1;
		}

	/* Better to abort than accidentally overwrite file... */

	if (!append && !force && (fp = fopen(OUTFILE,"r"))) {
		(void) fprintf(stderr,"Output file %s already exists\n",OUTFILE);
		(void) fclose(fp);
		return 1;
		}

	/* Check for default */

	if (!cm) {
		(void) fprintf(stderr,"No -m flag...assuming central mass == 1...\n");
		cm = 1;
		}

	/* Open summary file for output */

	if (!append && !(fp = fopen(OUTFILE,"w"))) {
		(void) fprintf(stderr,"Unable to open %s for writing\n",OUTFILE);
		return 1;
		}

	/* Loop through supplied files */

	for (i=optind;i<argc;++i) {
		(void) printf("%s:\n",argv[i]);
		if (!read_data(argv[i],&time,&n_data,&data,&sun,&jup,&sat,&ura,&nep))
			analyze(fp,time,n_data,data,sun,jup,sat,ura,nep,cm,i==argc-1);
		if (data)
			free((void *) data);
		}

	(void) fclose(fp);

	/* All done */

	return 0;
	}

/* ssa.c */
