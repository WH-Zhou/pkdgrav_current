/*
 ** ssic.c -- DCR 97-08-06
 ** ======
 ** Generates Solar System initial conditions.
 **
 ** 11/8/13: This version includes binary star support, written by Stefan Lines
 */

#include <ss.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> /* for getopt(), on some systems */
#include <boolean.h>
#include <rdpar.h>
#include <vector.h>
#include <coords.h>
#include <random.h>

#ifdef sparc
#ifdef sun
#undef sun /* sheesh */
#endif
#endif

/*#define SSTEST*/

/*#define TILT_RING*/ /*DEBUG used to tilt a portion of the disk at the start*/

#ifdef TILT_RING
#define TILT_R 1.275
#define TILT_DR 0.025
#define TILT_I 0.0125
#endif

/* Definitions */

#define RAYLEIGH 0
#define GAUSSIAN 1

typedef struct {

	/* Following data read in from parameter file... */

	double central_mass;
	BOOLEAN heliocentric;
	
	int binaryop;
	double binm1, binm2, separation, eccentricity;

	int n,dst_fnc;
	double total_mass,density,radius,scaling,r_inner,r_outer,surf_den_exp;
	double ecc_dsp,inc_dsp,ecc_max,inc_max;

	double seed_mass,seed_density,seed_radius,seed_scaling;
	double seed_sma,seed_ecc,seed_inc,seed_gap_scale;

	int n_planets;
	BOOLEAN tilt;
	char planet_data[MAXPATHLEN];

	double time;
	BOOLEAN adjust_com,softening;
	char output_file[MAXPATHLEN];

	/* Following derived from supplied parameters... */

	int n_data;
	double mass,red_hill,esc_vel,seed_half_gap;

    /* used in gen_planetesimal */
    double q, a, b;

} PARAMS;

static void
output_data(PARAMS *p,SSDATA *data)
{
	SSIO ssio;
	SSHEAD head;
	int i;

	if (ssioOpen(p->output_file,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",p->output_file);
		exit(1);
		}

	head.time = p->time;
	head.n_data = p->n_data;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	(void) ssioHead(&ssio,&head);

	for (i=0;i<p->n_data;i++) {
		data[i].org_idx = i; /* easiest time to store this */
		(void) ssioData(&ssio,&data[i]);
		}

	(void) ssioClose(&ssio);

#ifdef SSTEST
	{
	FILE *fp;
	long nb;

	fp = fopen(p->output_file,"r");
	(void) fseek(fp,0L,SEEK_END);
	nb = ftell(fp);
	/*DEBUG invalid: XDR uses SSHEAD_SIZE and SSDATA_SIZE*/
	(void) printf("I/O: nb=%li ==> n_data=%i compare with %i ",nb,
				  (int)(nb-sizeof(SSHEAD))/sizeof(SSDATA),p->n_data);
	(void) printf("(ss_head=%i, ss_data=%i)\n",sizeof(SSHEAD),sizeof(SSDATA));
	(void) fclose(fp);
	}
#endif
	}

static void
add_mom(SSDATA *data,VECTOR mom_pos,VECTOR mom_vel)
{
	VECTOR v;

	COPY_VEC(data->pos,v);
	SCALE_VEC(v,data->mass);
	ADD_VEC(mom_pos,v,mom_pos);
	COPY_VEC(data->vel,v);
	SCALE_VEC(v,data->mass);
	ADD_VEC(mom_vel,v,mom_vel);
	}

static void
get_com(SSDATA *data,int n,VECTOR com_pos,VECTOR com_vel,double *mass)
{
	int i;

	ZERO_VEC(com_pos);
	ZERO_VEC(com_vel);

	*mass = 0;

	for (i=0;i<n;i++) {
		add_mom(&data[i],com_pos,com_vel);
		*mass += data[i].mass;
		}

	NORM_VEC(com_pos,*mass);
	NORM_VEC(com_vel,*mass);
	}

static void
sub_com(SSDATA *data,int n,VECTOR com_pos,VECTOR com_vel)
{
	int i;

	for (i=0;i<n;i++) {
		SUB_VEC(data[i].pos,com_pos,data[i].pos);
		SUB_VEC(data[i].vel,com_vel,data[i].vel);
		}
	}

static void
calc_ang_mom(SSDATA *data, int n, VECTOR ang_mom)
{
	VECTOR v;
	int i;

	ZERO_VEC(ang_mom);

	for (i=0;i<n;i++) {
		CROSS(data[i].pos,data[i].vel,v);
		SCALE_VEC(v,data[i].mass);
		ADD_VEC(ang_mom,v,ang_mom);
		}
	}

static void
tilt(PARAMS *p,SSDATA *data)
{
	SSDATA *ptr;
	MATRIX r,r0;
	VECTOR zp,z,x,u;
	double cth,sth;
	int i;

	calc_ang_mom(data,p->n_planets,zp);
	NORM_VEC(zp,MAG(zp));
	SET_VEC(z,0,0,1);

	(void) printf("Inclination of invariable plane to ecliptic = %f deg\n",
				  acos(zp[Z])*(double)RAD_TO_DEG);

	/* Construct rotation matrix */
				  
	CROSS(z,zp,x);
	COPY_VEC(x,u);
	NORM_VEC(u,MAG(u));

	cth = DOT(z,zp);
	sth = MAG(x);

	UNIT_MAT(r);
	SCALE_MAT(r,cth);

	VecToMat(u,u,r0);
	SCALE_MAT(r0,(1 - cth));

	ADD_MAT(r,r0,r);

	SKEW_SYM_MAT(u,r0);
	SCALE_MAT(r0,sth);

	ADD_MAT(r,r0,r);

	/* Rotate planet positions and velocities */

	for (i=0;i<p->n_planets;i++) {
		ptr = &data[i];
		Transform(r,ptr->pos,u);
		COPY_VEC(u,ptr->pos);
		Transform(r,ptr->vel,u);
		COPY_VEC(u,ptr->vel);
		}
	}

static void
conv_orb_elem(double mu,double a,double e,double i,double lop,
			  double lan,double mean_anom,VECTOR pos,VECTOR vel)
{
	/* Wrapper for Stadel's Delaunay to heliocentric coordinates converter */

	struct delaunay input;
	struct helio output;

	input.mea = mean_anom;
	input.lop = lop;
	input.lan = lan;
	input.sma = a;
	input.ecc = e;
	input.inc = i;

	deltohel(mu,&input,&output);

	SET_VEC(pos,output.x,output.y,output.z);
	SET_VEC(vel,output.vx,output.vy,output.vz);

#ifdef SSTEST
	{
    double r,v2,h2;
    VECTOR h;

    (void) printf("Wanted %f %f %f got ",a,e,i);
    r = MAG(pos);
    v2 = MAG_SQ(vel);
    a = 1.0/(2/r - v2/mu);
    CROSS(pos,vel,h);
    h2 = MAG_SQ(h);
    e = (h2/(a*mu) > 1 ? 0 : sqrt(1 - h2/(a*mu)));
    i = acos(h[Z]/sqrt(h2));
    (void) printf("%f %f %f\n",a,e,i);
	}
#endif
	}

static double
radius(double m,double p)
{
	return pow(m/(4*PI*p/3),1.0/3);
	}

static void
get_planets(PARAMS *p,SSDATA *data)
{
	/* Reads "trq" style data file for planets */

	static const double v_conv = AU/(86400*V_SCALE); /* AU/d --> V_Earth */

	FILE *fp;
	SSDATA planets[9];
	double d;
	int i,j;

	if (!(fp = fopen(p->planet_data,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",p->planet_data);
		exit(1);
		}

	(void) printf("Reading planet data...");

	/* inverse masses, extra Earth stuff, and time */

	for (i=0;i<9;i++) {
		(void) fscanf(fp,"%lf",&d);
		planets[i].mass = 1/d;
		if (i == 2) (void) fscanf(fp,"%lf%lf%lf",&d,&d,&d);
		}

	(void) fscanf(fp,"%lf",&d);
			
	(void) printf("JD = %f...",d);

	/* positions and velocities */

	for (i=0;i<9;i++) {
		for (j=0;j<N_DIM;j++)
			(void) fscanf(fp,"%lf",&planets[i].pos[j]);
		for (j=0;j<N_DIM;j++)
			(void) fscanf(fp,"%lf",&planets[i].vel[j]);
		SCALE_VEC(planets[i].vel,v_conv);
		}

	/* save desired planet data */

	for (i=0;i<p->n_planets;i++) {
		data[i] = planets[4 + i]; /* start at Jupiter */
		ZERO_VEC(data[i].spin);
		}

	if (p->n_planets > 0) {
		data[0].radius = 4.7727e-4;
		data[0].color = JUPITER;
		if (p->n_planets > 1) {
			data[1].radius = 1.5e-4;
			data[1].color = SATURN;
			if (p->n_planets > 2) {
				data[2].radius = 1.7e-4;
				data[2].color = URANUS;
				if (p->n_planets > 3) {
					data[3].radius = 1.62e-4;
					data[3].color = NEPTUNE;
					}
				}
			}
		}

	(void) fclose(fp);
	(void) printf("done!\n");
	}

static void
gen_planets(PARAMS *p,SSDATA *data)
{
	SSDATA *ptr;
	double t,mean_anom;
	int i;

	/*
	 * Heliocentric orbital elements of the planets, epoch 1980.0
	 * [Table 7 of Practical Astronomy with your Calculator, 2nd ed.].
	 *
	 */

	typedef struct {
		double	mass;	/* mass in M_SUN */
		double	radius;	/* radius in AU */
		double	period;	/* period in TROP_YR */
		double	a;		/* semi-major axis in AU */
		double	e;		/* eccentricity */
		double	i;		/* inclination in deg */
		double	lop;	/* longitude of perihelion in deg */
		double	lan;	/* longitude of ascending node in deg */
		double	lng;	/* longitude at epoch in deg */
		int		color;
		} ELEM_T;

	ELEM_T elem[MAX_NUM_PLANETS] = {
		{9.5480e-4, 4.7727e-4,  11.86224,  5.202561, 0.0484658, 1.3041819,
			  14.0095493, 100.2520175, 146.966365 , JUPITER},
		{2.8588e-4, 1.5e-4   ,  29.45771,  9.554747, 0.0556155, 2.4893741,
			  92.6653974, 113.4888341, 165.322242 , SATURN},
		{4.355e-5 , 1.7e-4   ,  84.01247, 19.21814 , 0.0463232, 0.7729895,
			 172.7363288,  73.8768642, 228.0708551, URANUS},
		{5.1777e-5, 1.62e-4  , 164.79558, 30.10957 , 0.0090021, 1.7716017,
			  47.8672148, 131.5606494, 260.3578998, NEPTUNE}
		};

	assert(p->central_mass == 1.0);

	/* Giant planets */

	t = 100*randUniform(); /*DEBUG between 0 and 100 tropical years since 1980 */

	for (i=0;i<p->n_planets;i++) {
		ptr = &data[i];
		ptr->mass = elem[i].mass;
		ptr->radius = elem[i].radius;
		ptr->color = elem[i].color;

		ZERO_VEC(ptr->spin);

		/* Convert from degrees to radians */

		elem[i].i	*= DEG_TO_RAD;
		elem[i].lop	*= DEG_TO_RAD;
		elem[i].lan	*= DEG_TO_RAD;
		elem[i].lng	*= DEG_TO_RAD;

		/* Use random time for mean anomaly to get initial pos & vel */

		mean_anom = TWO_PI*t/elem[i].period + elem[i].lng - elem[i].lop;

		conv_orb_elem(p->central_mass + ptr->mass,elem[i].a,elem[i].e,
					  elem[i].i,elem[i].lop,elem[i].lan,mean_anom,
					  ptr->pos,ptr->vel);
		}
	}

static void
gen_planetesimal(const PARAMS * const p, SSDATA * const data)
{
	double f,sma,ecc=0,inc=0;

	data->mass = p->mass;
	data->radius = p->radius;
	data->color = PLANETESIMAL;
	ZERO_VEC(data->spin);

	/*
	 * Semi-major axis (strategy: the surface density is assumed to be a
	 * kind of "average" formed by the planetesimal semi-major axes,
	 * rather than the exact value formed by the instantaneous positions).
	 *
	 */

	do {
		f = randUniform();
		sma = pow((1.0 - f)*p->a + f*p->b,1.0/p->q);
		} while (p->seed_mass && fabs(sma - p->seed_sma) <= p->seed_half_gap);
	/*
	 * Eccentricity & inclination (strategy: assume uncorrelated and
	 * use user-supplied distribution function and dispersions).
	 *
	 */

	do {
		switch (p->dst_fnc) {
		case RAYLEIGH:
			ecc = randRayleigh()*p->ecc_dsp;
			assert(ecc >= 0.0);
			break;
		case GAUSSIAN:
			do
				ecc = randGaussian()*p->ecc_dsp;
			while (ecc < 0);
			break;
		default:
			assert(0);
			}
		}
	while (ecc > p->ecc_max);

	do {
		switch (p->dst_fnc) {
		case RAYLEIGH:
			inc = randRayleigh()*p->inc_dsp;
			assert(inc >= 0.0);
			break;
		case GAUSSIAN:
			do
				inc = randGaussian()*p->inc_dsp;
			while (inc < 0);
			break;
		default:
			assert(0);
			}
		}
	while (inc > p->inc_max);

	/*
	 * Convert to instantaneous position and velocity, using uniform
	 * random longitude of perihelion, longitude of ascending node, and
	 * instantaneous longitude (valid if ecc & inc are small).
	 *
	 */

#ifdef TILT_RING
	if (fabs(sma - TILT_R) < TILT_DR) {
		conv_orb_elem(p->central_mass + data->mass,sma,ecc,TILT_I + inc,
					  TWO_PI*randUniform(),0,TWO_PI*randUniform(),data->pos,data->vel);
		}
	else
		conv_orb_elem(p->central_mass + data->mass,sma,ecc,inc,TWO_PI*randUniform(),
					  TWO_PI*randUniform(),TWO_PI*randUniform(),data->pos,data->vel);
#else
	conv_orb_elem(p->central_mass + data->mass,sma,ecc,inc,TWO_PI*randUniform(),
				  TWO_PI*randUniform(),TWO_PI*randUniform(),data->pos,data->vel);
#endif

	}

static void
gen_planetesimals(PARAMS *p,SSDATA *data)
{
	SSDATA *d;
	int i;

	if (p->binaryop==1) {
	  int x;
	  double apoa, apob, va, vb, binma, binmb, dsep;
	  for (x=0;x<2;x++) {
	    data[x].radius = 6e8/AU;
	    data[x].pos[1] = 0.0;
        data[x].pos[2] = 0.0;
        data[x].vel[0] = 0.0;
        data[x].vel[2] = 0.0;
        data[x].spin[0] = 0.0;
        data[x].spin[1] = 0.0;
        data[x].spin[2] = 0.0;
        data[x].color = 8;
	  }
	  
	  binma=p->binm1*1.98855e30;
	  binmb=p->binm2*1.98855e30;
	  dsep=p->separation*AU;
	  
	  printf("Binary Masses - A: %f B: %f\n",p->binm1,p->binm2);
	  
	  /* Apoapses in SI */

	  apoa = (dsep*binmb*(1+p->eccentricity))/(binma+binmb);
	  apob = (dsep*binma*(1+p->eccentricity))/(binma+binmb);
	  
	  /* Calculate the y-component velocity at apoapses */
	  
	  va = sqrt((G*(1-p->eccentricity)*binmb*binmb)/((binma+binmb)*(1+p->eccentricity)*dsep));
	  printf("Va is %e km/s\n",(va/1000));
	  
	  /* Convert back into system units */
	  
	  va = va/29785.4;
	  
	  /* Ditto */
	  
	  vb = sqrt((G*(1-p->eccentricity)*binma*binma) /
                ((binma+binmb)*(1+p->eccentricity)*dsep));
	  printf("Vb is %e km/s\n",(vb/1000));
	  vb = vb/29785.4;
	  
	  printf("Apoapsis Velocity A: %f\n",va);
	  printf("Apoapsis Velocity B: %f\n",vb);
	  
	  /* Convert back into AU */
	  
	  apoa = apoa / AU;
	  apob = apob / AU;
	  
	  printf("Apo A: %e\n",apoa);
	  printf("Apo B: %e\n",apob);
	  
	  /* Assign particle 0 and 1 as the binary system */
	  /* With all calculated or read-in values */
	  
	  data[0].mass = p->binm1;
      data[0].pos[0] = apoa*(-1.0);
      data[0].vel[1] = va*(-1.0);
      data[1].mass = p->binm2;
      data[1].pos[0] = apob;
      data[1].vel[1] = vb;
	  
	  for (i=2;i<p->n;i++) {
		d = &data[i];
		gen_planetesimal(p,d);
	  }
	}
	else {
	  for (i=0;i<p->n;i++) {
		d = &data[i];
		gen_planetesimal(p,d);
		}
	}
}

static void
gen_seed(PARAMS *p,SSDATA *d)
{
	d->mass = p->seed_mass;
	d->radius = p->seed_radius;
	conv_orb_elem(p->central_mass + p->seed_mass,p->seed_sma,p->seed_ecc,
				  p->seed_inc,0,0,0,d->pos,d->vel);
	ZERO_VEC(d->spin);
	d->color = JUPITER; /*DEBUG*/
	}

static void
generate(PARAMS *p,SSDATA *data)
{
	SSDATA *sun,*planets,*seed,*planetesimals,*ptr;
	VECTOR com_pos,com_vel;
	double total_mass;
	int i;

	/* Get handy pointers */

	if (p->heliocentric) {
		sun = NULL;
		planets = data;
		}
	else {
		sun = data;
		planets = sun + 1;
		}

	seed = planets + p->n_planets;

	if (p->seed_mass)
		planetesimals = seed + 1;
	else
		planetesimals = seed;

	/* Start with the Sun in the heliocentric frame */

	if (!p->heliocentric) {
		ZERO_VEC(sun->pos);
		ZERO_VEC(sun->vel);
		ZERO_VEC(sun->spin);
		sun->mass = 1;
		sun->radius = R_SUN/L_SCALE;
		sun->color = SUN;
    }

	/* Add planets, seed, and planetesimals as appropriate */

	if (p->n_planets) {
		if (strlen(p->planet_data))
            get_planets(p,planets);
		else
            gen_planets(p,planets);
    }

	if (p->seed_mass)
        gen_seed(p,seed);

	if (p->n)
        gen_planetesimals(p,planetesimals);

	/*
	 * Must now tilt the planets so the planetesimal disk lies in the
	 * invariable plane defined by the net angular momenta of the planets.
	 *
	 */

	if (p->tilt) {
		tilt(p,planets);
#ifdef SSTEST
		{
		VECTOR z,zp;

		assert(p->seed_mass == 0.0);
		calc_ang_mom(planetesimals,p->n,z);
		calc_ang_mom(planets,p->n_planets,zp);

		NORM_VEC(z,MAG(z));
		NORM_VEC(zp,MAG(zp));

		(void) printf("net planetesimal ang mom direction = (%e,%e,%e)\n",
					  z[X],z[Y],z[Z]);

		(void) printf("-- will assume (0,0,1)\n");

		SET_VEC(z,0,0,1);

		SUB_VEC(z,zp,z);

		(void) printf("invar plane: z - zp = (%e,%e,%e)\n",z[X],z[Y],z[Z]);
		}
#endif
		}

	if (p->adjust_com) {
		/* get planetesimal moments and estimate of mass centering error */
		if (p->n) {
			get_com(planetesimals,p->n,com_pos,com_vel,&total_mass);
			(void) printf("Planetesimal c-o-m (pos,vel) offset mag = (%e,%e)\n",
						  MAG(com_pos),MAG(com_vel));
			SCALE_VEC(com_pos,total_mass);
			SCALE_VEC(com_vel,total_mass);
			}
		else {
			ZERO_VEC(com_pos);
			ZERO_VEC(com_vel);
			total_mass = 0;
			}

		/* now shift everything so centre-of-mass is at origin */

		assert(p->seed_mass == 0.0);

		if (!p->heliocentric)
			total_mass += sun->mass;

		for (i=0;i<p->n_planets;i++) {
			ptr = &planets[i];
			add_mom(ptr,com_pos,com_vel);
			total_mass += ptr->mass;
			}

		NORM_VEC(com_pos,total_mass);
		NORM_VEC(com_vel,total_mass);

		sub_com(data,p->n_data,com_pos,com_vel);
		}
	}

static void
fix_rejects(PARAMS *p)
{
	FILE *fp;
	SSIO ssio;
	SSDATA data;
#ifdef SSIO_USE_MPI
	MPI_Offset sh,sd;
#else
    unsigned int sh,sd;
#endif
    int rej,nRej=0,i;

	if (!(fp = fopen(REJECTS_FILE,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",REJECTS_FILE);
		exit(1);
		}

	(void) fscanf(fp,"%i",&i);
	if (i != p->n_data) {
		(void) fprintf(stderr,"Incompatible rejects file.\n");
		exit(1);
		}

	if (ssioOpen(p->output_file,&ssio,SSIO_UPDATE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",p->output_file);
		exit(1);
		}

	sh = SSHEAD_SIZE; /* NOT sizeof(SSHEAD)! XDR has special data sizes */
	sd = SSDATA_SIZE; /* NOT sizeof(SSDATA)! Ditto */

	(void) ssioSetPos(&ssio,sh);

	i = 0;

	while (fscanf(fp,"%i",&rej) != EOF) {
		if (rej >= 0) {
			if (i != rej) {
				(void) fprintf(stderr,"Corrupted rejects data.\n");
				exit(1);
				}
			assert(p->seed_mass == 0.0 || i > 0);
			gen_planetesimal(p,&data);
			data.org_idx = i;
			(void) ssioSetPos(&ssio,sh + (unsigned int) i*sd);
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

static void
get_params(PARAMS *p)
{
	/* Read in supplied parameters */

	OpenPar("ssic.par");

	/* central object */

	ReadDbl("Central mass",&p->central_mass);
	assert(p->central_mass >= 0.0);
	ReadInt("Use heliocentric frame?",&p->heliocentric);

	/* binary */
	
	ReadInt("Insert binary?",&p->binaryop);
	assert(p->binaryop == 0 || p->binaryop == 1);
	ReadDbl("Star A Mass",&p->binm1);
	assert(p->binm1 >= 0.0);
	ReadDbl("Star B Mass",&p->binm2);
	assert(p->binm2 >= 0.0);
	ReadDbl("Separation",&p->separation);
	assert(p->separation >= 0.0);
	ReadDbl("Eccentricity",&p->eccentricity);
	assert(p->eccentricity >= 0.0 && p->eccentricity < 1.0);

	/* planetesimals */

	ReadInt("Number of planetesimals",&p->n);
	assert(p->n >= 0);

	if (p->n) {
		ReadDbl("Total mass",						&p->total_mass);
		assert(p->total_mass > 0.0);
		ReadDbl("Planetesimal density",				&p->density);
		assert(p->density >= 0.0);
		ReadDbl("Planetesimal radius",				&p->radius);
		assert(p->density || p->radius);
		assert(p->density == 0.0 || p->radius == 0.0);
		ReadDbl("Planetesimal radius scaling",		&p->scaling);
		assert(p->scaling > 0.0);
		ReadDbl("Inner orbital radius",				&p->r_inner);
		assert(p->r_inner > 0.0);
		ReadDbl("Outer orbital radius",				&p->r_outer);
		assert(p->r_outer >= p->r_inner);
		ReadDbl("Projected surface density exp.",	&p->surf_den_exp);
		assert(p->surf_den_exp != -2.0);
		ReadDbl("Eccentricity dispersion",			&p->ecc_dsp);
		assert(p->ecc_dsp >= 0.0);
		ReadDbl("Inclination dispersion",			&p->inc_dsp);
		assert(p->inc_dsp >= 0.0);
		ReadDbl("Maximum eccentricity",				&p->ecc_max);
		assert(p->ecc_max >= 0.0 && p->ecc_max <= 1.0);
		ReadDbl("Maximum inclination",				&p->inc_max);
		assert(p->inc_max >= 0.0 && p->inc_max <= 90.0);
		ReadInt("Distribution functions",			&p->dst_fnc);
		assert(p->dst_fnc == RAYLEIGH || p->dst_fnc == GAUSSIAN);
		}

	/* seed */

	ReadDbl("Seed mass",&p->seed_mass);
	if (p->seed_mass) {
		ReadDbl("Seed density",&p->seed_density);
		assert(p->seed_density >= 0.0);
		ReadDbl("Seed radius",&p->seed_radius);
		assert(p->seed_density || p->seed_radius);
		assert(p->seed_density == 0.0 || p->seed_radius == 0.0);
		ReadDbl("Seed radius scaling",&p->seed_scaling);
		assert(p->seed_scaling > 0.0);
		ReadDbl("Seed semi-major axis",&p->seed_sma);
		assert(p->seed_sma > 0.0);
		ReadDbl("Seed eccentricity",&p->seed_ecc);
		assert(p->seed_ecc >= 0.0 && p->seed_ecc < 1.0);
		ReadDbl("Seed inclination",&p->seed_inc);
		assert(p->seed_inc >= 0.0 && p->seed_inc <= 180.0);
		ReadDbl("Seed gap scaling",&p->seed_gap_scale);
		assert(p->seed_gap_scale >= 0.0);
		}

	/* planets */

	ReadInt("Number of planets",&p->n_planets);
	assert(p->n_planets >= 0);
	assert(p->n_planets == 0 || p->central_mass == 1.0);
	p->tilt = FALSE;
	if (p->n_planets) {
		if (p->n)
            ReadInt("Align invariable plane?",&p->tilt);

		ReadStr("Planet data file name",p->planet_data,MAXPATHLEN);
		}

	/* miscellaneous */

	ReadDbl("Start time",				&p->time);
	ReadInt("Adjust c-o-m to origin?",	&p->adjust_com);
	ReadInt("Radius is softening?",		&p->softening);
	ReadStr("Output file name",			p->output_file,MAXPATHLEN);

	ClosePar();

	/* Calculate remaining parameters */

	p->n_data = p->n_planets + p->n; /* planets + planetsimals */

	if (!p->heliocentric)
        ++p->n_data; /* add Sun */

	if (p->seed_mass)
        ++p->n_data; /* add seed */

	if (p->n) {

		p->mass = (p->total_mass*M_EARTH/p->n); /* equal mass planetesimals */
		(void) printf("Planetesimal mass = %e g = %g M_Sun\n",1000*p->mass,
					  p->mass/M_SCALE);

		if (p->density)
			p->radius = radius(p->mass,1000*p->density);
		else if (p->radius < 0)
			p->radius *= -1000;
		else
			p->radius *= L_SCALE; /* for now */

		if (p->softening)
            p->scaling *= 2; /* for pkdgrav COLLISIONS */

		p->radius *= p->scaling;
		(void) printf("Planetesimal radius = %g km",0.001*p->radius);

		if (p->scaling != 1)
            (void) printf(" (scaled)");

		(void) printf(" = %g AU\n",p->radius/L_SCALE);
		/*
		 * IMPORTANT NOTE: Kokubo & Ida (1998) define the reduced Hill radius
		 * in terms of the SUM of the two masses involved, in this case twice
		 * the planetesimal mass ASSUMING EQUAL-MASS PLANETESIMALS. It's
		 * important to be consistent for the initial eccentricity and
		 * inclination dispersion scaling.
		 */
		p->red_hill = pow(2*p->mass/(3*p->central_mass*M_SUN),1.0/3);
		(void) printf("Reduced Hill radius = %g\n",p->red_hill);
		(void) printf("Hill radius at 1 AU = %g km\n",0.001*p->red_hill*AU);
		p->esc_vel = sqrt(2*G*p->mass/p->radius);
		(void) printf("Planetesimal surface escape speed = %g km/s\n",
					  0.001*p->esc_vel);
		if (p->r_outer > p->r_inner)
			(void) printf("Mean surface density = %g g/cm^2\n",
						  1000*p->mass*p->n/
						  (PI*(SQ(p->r_outer) - SQ(p->r_inner))*SQ(AU)*10000));
		p->mass /= M_SCALE;
		p->radius /= L_SCALE;
		p->ecc_dsp *= p->red_hill;
		p->inc_dsp *= p->red_hill;

		if (p->ecc_max > 0.1)
			(void) fprintf(stderr,"WARNING: Large maximum eccentricity\n");

		if (p->inc_max > 10)
			(void) fprintf(stderr,"WARNING: Large maximum inclination\n");

		p->inc_max *= DEG_TO_RAD;
		(void) printf("Predicted velocity dispersion = %g\n",
					  sqrt(SQ(p->ecc_dsp) + SQ(p->inc_dsp)));

		p->q = p->surf_den_exp + 2.0;
		p->a = pow(p->r_inner, p->q);
		p->b = pow(p->r_outer, p->q);
		}

	if (p->seed_mass) {

		if (p->seed_density)
			p->seed_radius = radius(p->seed_mass*M_SCALE,p->seed_density*1000);
		else if (p->seed_radius < 0)
			p->seed_radius *= -1000;

		if (p->softening)
            p->seed_scaling *= 2;

		p->seed_radius *= p->seed_scaling;
		(void) printf("Seed radius = %g km ",p->seed_radius/1000);

		if (p->seed_scaling)
            (void) printf("(scaled) ");

		p->seed_radius /= L_SCALE;
		(void) printf("= %g AU\n",p->seed_radius);
		p->seed_half_gap = 1.5*p->seed_sma*p->seed_gap_scale*
			pow(p->seed_mass/p->central_mass,2.0/7); /* chaotic zone */
		assert(p->r_outer - p->r_inner > 2*p->seed_half_gap);
		(void) printf("Seed gap width = %g AU\n",2*p->seed_half_gap);
		p->seed_inc *= DEG_TO_RAD;
    }

	p->time *= (JUL_YR/T_SCALE);

}

void
generate_chunk
(PARAMS orig_params,
 SSDATA *orig_data);

int
main(int argc,char *argv[])
{
	PARAMS params;
	SSDATA *data = NULL;
	int seed,force=0,rejects=0,usage=0,c;

	/* Disable stdout buffering */

	setbuf(stdout,(char *) NULL);

	/* Generate a random number seed */

	seed = getpid();

	/* Check arguments */

	while ((c = getopt(argc,argv,"frs:")) != EOF)
		switch (c) {
		case 'f':
			force = 1;
			break;
		case 'r':
			rejects = 1;
			break;
		case 's':
			seed = atoi(optarg);
			if (seed <= 0) usage = 1;
			break;
		default:
			usage = 1;
			}

#ifdef SSIO_USE_MPI
    /* need to make sure the seed is different in each process! */
    ssioInitialise();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    seed = rank + seed;
#endif

	if (usage || optind < argc) {
		(void) fprintf(stderr,"Usage: %s [ -f ] [ -r ] [ -s seed ]\nwhere seed > 0 is random number generator seed\n",argv[0]);
		exit(1);
		}

	/* Get model parameters */

	get_params(&params);

	if (rejects && params.adjust_com)
		(void) fprintf(stderr,"WARNING: adj_com ignored for rejects.\n");

	if (!force && !rejects) {
		FILE *fp = fopen(params.output_file,"r");
		if (fp) {
			(void) fprintf(stderr,"%s NOT OVERWRITTEN.\n",params.output_file);
			(void) fclose(fp);
			exit(1);
			}
		}

	/* Generate initial conditions */

    randSeedGenerator(randReadUrandom());
    
	if (rejects) {
		fix_rejects(&params);
		return 0;
		}

#ifdef SSIO_USE_MPI
    ssioInitialise();
    generate_chunk(params, data);
    return 0;
#endif

	if (!(data = (SSDATA *) malloc((size_t) params.n_data*sizeof(SSDATA)))) {
		(void) fprintf(stderr,"Unable to allocate data memory\n");
		exit(1);
		}

	generate(&params,data);

	/* Save data */

	output_data(&params,data);

	/* All done */

	free((void *) data);

	return 0;
	}

/* ssic.c */
