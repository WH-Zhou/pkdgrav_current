/*
 ** rpg.c -- DCR 98-09-18
 ** =====
 ** Rubble pile generator. Currently assumes equal-size particles.
 */

#include <rpu.h>
#include <math.h>
#include <assert.h>
#include <rdpar.h>
#include <random.h>

/* 
** FPE trapping for Linux/glibc 2.2 and later.  Some compilers (e.g.,
** gcc, icc) recognize feenableexcept() without fully conforming to
** the C99 standard.  Since FPE trapping is expensive, we keep it off
** by default; set TRAP_FPE to 1 during compilation to turn it on.
** WARNING! The Intel compiler (icc) will still ignore some FPEs
** unless the "-mp" compile flag is also used!
*/
#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
/* feenableexcept() is part of the C99 standard */
# if !__APPLE__ /* FPE handling not supported under macOS */
#  define _GNU_SOURCE
#  include <fenv.h>
int feenableexcept(int);
# endif /* !__APPLE__ */
#endif

/*#define EFFICIENCY*/

#define DFLT_PARAMS_FILE "rpg.par"

typedef struct {
	char output_file[MAXPATHLEN];
	double mass,radius,scaling,speed_max,density,spin; /* particle properties */
	} PARAMS;

static void
write_data(PARAMS *p,RUBBLE_PILE *rp)
{
	SSIO ssio;
	SSHEAD head;
	int i;

	if (ssioOpen(p->output_file,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",p->output_file);
		exit(1);
		}

	head.time = 0;
	head.n_data = rp->n_particles;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&head)) {
		(void) fprintf(stderr,"Unable to write header.\n");
		exit(1);
		}

	for (i=0;i<rp->n_particles;i++)
		if (ssioData(&ssio,&rp->data[i])) {
			(void) fprintf(stderr,"Error writing data (particle %i).\n",i);
			exit(1);
			}

	(void) ssioClose(&ssio);
	}

static void
ran_vel(double speed_max,VECTOR v)
{
	double v_max = speed_max/sqrt(N_DIM);

	do {
		SET_VEC(v,randUniform(),randUniform(),randUniform());
		SCALE_VEC(v,v_max);
		} while (MAG_SQ(v) > SQ(speed_max));
	}

static void
calc_data(PARAMS *p,RUBBLE_PILE *rp)
{
	VECTOR pos;
	double dx,dy,dz,rmax,speed_esc;
	int nx,ny,nz,np,ix,iy,iz,i;

	rp->n_particles = (int) (1.1 * rp->n_particles); /* safety factor */
	rpuMalloc(rp);

	/* hexagonal closest packing */

	dx = 2*p->radius;
	dy = sqrt(3.0)*p->radius;
	dz = (2.0/3)*sqrt(6.0)*p->radius;

	rmax = rp->radius - p->radius;

	nx = (int) (rmax/dx);
	ny = (int) (rmax/dy);
	nz = (int) (rmax/dz);

	np = 0;

	for (ix=-nx;ix<=nx;ix++)
		for (iy=-ny;iy<=ny;iy++)
			for (iz=-nz;iz<=nz;iz++) {
				pos[X] = ix*dx + ((iy + iz)%2)*0.5*dx;
				pos[Y] = iy*dy - ((iz%4 - 2*SGN(iz))%2)*(dy/3);
				pos[Z] = iz*dz;
				if (rpuInEllipsoid(pos,p->radius,rp->axis_len)) {
					if (np == rp->n_particles) {
						(void) fprintf(stderr,"hcp(): Too many particles.\n");
						exit(1);
						}
					COPY_VEC(pos,rp->data[np].pos);
					++np;
					}
				}

	rp->n_particles = np;
	rpuRealloc(rp); /* remove any wasted storage */

	/* Fill in remaining particle data fields */

	p->radius *= p->scaling;

	for (i=0;i<rp->n_particles;i++) {
		rp->data[i].mass = p->mass; /* readjusted below */
		rp->data[i].radius = p->radius;
		ZERO_VEC(rp->data[i].vel); /* recalculated below */
		ZERO_VEC(rp->data[i].spin); /* ditto */
		}

	rpuApplyColor(rp); /* sets rp->data[i].color */
	rpuApplyAggID(rp); /* sets rp->data[i].org_idx */

	/*
	 ** Adjust bulk properties. Generally the reported packing
	 ** efficiency will be a bit too small, and the particle density
	 ** slightly too large, owing to the difficulty in fitting
	 ** ellipsoids to the data.
	 */

	rpuAnalyze(rp);

	rp->radius = rp->axis_len[X]; /* first axis is longest */
	rp->packing = rp->n_particles*rpuVolSph(p->radius)/rpuVolEll(rp->axis_len);
	p->density = rp->density/rp->packing;
	p->mass = rpuVolSph(p->radius)*p->density;
	rp->mass = rp->n_particles*p->mass;

	/* Set new particle masses */

	for (i=0;i<rp->n_particles;i++)
		rp->data[i].mass = p->mass;

	/* Escape speed from particle surface */

	speed_esc = sqrt(2*rp->data->mass/rp->data->radius);

	/* Set new particle velocities */

	for (i=0;i<rp->n_particles;i++)
		ran_vel(p->speed_max*speed_esc,rp->data[i].vel);

	/* Adjust c-o-m position & velocity to zero */

	rpuAnalyze(rp);
	SCALE_VEC(rp->pos,-1.0);
	rpuApplyPos(rp);
	SCALE_VEC(rp->vel,-1.0);
	rpuApplyVel(rp);

	/* Set spin */

	SCALE_VEC(rp->spin,-1.0);
	rpuAddSpin(rp,FALSE);
	ZERO_VEC(rp->spin);
	rp->spin[Z] = p->spin;
	rpuAddSpin(rp,FALSE);

	/* Output new parameters */

#ifdef EFFICIENCY
	(void) printf("actual R/r %f N %i eff %f\n",rp->radius/p->radius,
				  rp->n_particles,rp->packing);
#else
	(void) printf("Adjusted parameters:\n");
	(void) printf("(a,b,c)=(%g,%g,%g) N=%i r=%g rho_p=%g\n",
		      rp->axis_len[X]*L_SCALE,rp->axis_len[Y]*L_SCALE,
		      rp->axis_len[Z]*L_SCALE,rp->n_particles,
		      p->radius*L_SCALE,p->density*D_SCALE);
	(void) printf("bulk mass=%g bulk density=%g\n",rp->mass*M_SCALE,rp->density*D_SCALE);
	(void) printf("Actual packing efficiency = %f\n",rp->packing);
	(void) printf("Additional characteristics:\n");
	(void) printf("Confining radius = %g (surface gravity = %g)\n",rp->radius*L_SCALE,rp->mass/SQ(rp->radius)*V_SCALE/T_SCALE);
	{
		double dRadEq = pow(rp->axis_len[X]*rp->axis_len[Y]*rp->axis_len[Z],1.0/3.0);
		(void) printf("Radius of equivalent sphere = %g (surface gravity = %g)\n",dRadEq*L_SCALE,rp->mass/SQ(dRadEq)*V_SCALE/T_SCALE);
		}
	(void) printf("Escape speed from surface of particle = %g\n",
				  speed_esc*V_SCALE);
	(void) printf("Dynamical time for particle = %g\n",
				  sqrt(1/p->density)*T_SCALE);
	(void) printf("Escape speed from surface of rubble pile = %g\n",
				  sqrt(2*rp->mass/rp->radius)*V_SCALE);
	(void) printf("Surface speed at equator due to rotation = %g\n",
				  MAG(rp->spin)*rp->radius*V_SCALE);
	(void) printf("Dynamical time for rubble pile = %g\n",
				  sqrt(1/rp->density)*T_SCALE);
#endif /* !EFFICIENCY */
	}

#define PEMAX 0.74 /* theoretical maximum packing efficiency */

/*
 ** Following packing-efficiency formulae derived empirically...
 ** (cf. "EFFICIENCY" macro and sm script "rpgeff.sm")
 */

#define PEN(n) ((n) == 1 ? 1.0 : 2*PEMAX/PI*atan(0.496*pow((n),0.264)))
#define PER(r) ((r) == 1 ? 1.0 : 2*PEMAX/PI*atan(0.368*pow((r),0.843)))

static void
get_params(const char *parfile,PARAMS *p,RUBBLE_PILE *rp)
{
	double vf;
	int color,agg_id;

	OpenPar(parfile);

	ReadNDbl("Bulk semi-axes",rp->axis_len,N_DIM);
	ReadDbl("Bulk density",&rp->density);
	ReadInt("Number of particles",&rp->n_particles);
	ReadInt("Color",&color);
	ReadDbl("Particle radius",&p->radius);
	ReadDbl("Particle radius scaling",&p->scaling);
	ReadDbl("Particle max speed",&p->speed_max);
	ReadDbl("Particle density",&p->density);
	ReadDbl("Initial spin period",&p->spin);
	ReadInt("Aggregate ID",&agg_id);
	ReadStr("Output file",p->output_file,MAXPATHLEN);

	ClosePar();

	/* Sanity checks */

	assert(SGN(rp->axis_len[X]) == SGN(rp->axis_len[Y]) &&
		   SGN(rp->axis_len[Y]) == SGN(rp->axis_len[Z]));
	assert(rp->density >= 0.0);
	assert(rp->n_particles >= 0);
	assert(color <= 255); /* allow negative colors <--> particles stuck on walls */
	assert(p->radius >= 0.0);
	assert(p->scaling > 0.0);
	assert(p->density >= 0.0);
	assert(agg_id >= -1);

	assert(rp->axis_len[X] > 0.0 || (rp->n_particles > 0 && p->radius > 0.0));
	assert(rp->density > 0.0 || p->density > 0.0);
	assert(rp->n_particles > 0 || p->radius > 0.0);

	if (p->scaling > 1) {
		(void) fprintf(stderr,"WARNING: particle radius scaling > 1...\n");
		(void) fprintf(stderr,"...particles may overlap.\n");
		}

	rp->color = color;
	rp->agg_id = agg_id;

	/* Fill in any missing parameters */

	if (rp->n_particles)
		rp->packing = PEN(rp->n_particles);
	else if (rp->radius && p->radius)
		rp->packing = PER(rp->radius/p->radius);
	else
		rp->packing = 0.0; /*DEBUG this breaks! for now */

	if (rp->axis_len[X] <= 0) {
		double x;
		if (rp->axis_len[X] == 0) { /* assume sphere */
			x = p->radius*pow(rp->n_particles/rp->packing,1.0/3);
			SET_VEC(rp->axis_len,x,x,x);
			}
		else { /* axes are actually ratios in this case */
			SCALE_VEC(rp->axis_len,-1);
			x = p->radius*pow(rp->n_particles/
							  (rp->packing*rp->axis_len[X]*rp->axis_len[Y]*
							   rp->axis_len[Z]),1.0/3);
			SCALE_VEC(rp->axis_len,x);
			}
		}

	/* bulk radius */
	rp->radius = MAX(MAX(rp->axis_len[X],rp->axis_len[Y]),rp->axis_len[Z]);

	vf = rpuVolEll(rp->axis_len)/rpuVolSph(rp->radius);

	if (rp->density == 0)
		rp->density = rp->packing*p->density;
	if (!rp->n_particles)
		rp->n_particles =
		  (int) (vf*rp->packing/CUBE(p->radius/rp->radius));
	if (p->radius == 0)
		p->radius =
		  rp->radius*pow(vf*PEN(rp->n_particles/vf)/rp->n_particles,1.0/3);
	if (p->density == 0)
		p->density = rp->density/rp->packing;
	
	/* Compute remaining derived parameters */

	p->mass = rpuVolSph(p->radius)*p->density; /* particle mass */
	rp->mass = rp->n_particles*p->mass; /* bulk mass */

	(void) printf("Nominal rubble pile parameters (all units MKS):\n");
	(void) printf("(a,b,c)=(%g,%g,%g) rho_b=%g N=%i r=%g rho_p=%g\n",
				  rp->axis_len[X],rp->axis_len[Y],rp->axis_len[Z],rp->density,
				  rp->n_particles,p->radius,p->density);
	(void) printf("bulk mass=%g particle mass=%g color=%i\n",
				  rp->mass,p->mass,rp->color);
	(void) printf("Estimated packing efficiency: e_N=%f e_R/r=%f\n",
				  PEN(rp->n_particles),PER(rp->radius/p->radius));

	/* Convert units */

	NORM_VEC(rp->axis_len,L_SCALE);
	rp->density /= D_SCALE;
	rp->radius /= L_SCALE;
	rp->mass /= M_SCALE;
	if (p->spin != 0.0)
		p->spin = TWO_PI/(p->spin*3600.0)*T_SCALE;
	p->density /= D_SCALE;
	p->radius /= L_SCALE;
	p->mass /= M_SCALE;

	if (rp->n_particles == 1) /* special case for single particle: avoid roundoff */
		rp->axis_len[X] = rp->axis_len[Y] = rp->axis_len[Z] = rp->radius = p->radius;
	}

#undef PER
#undef PEN
#undef PEMAX

#ifdef EFFICIENCY
int
main(int argc,char *argv[])
{
	PARAMS params;
	RUBBLE_PILE rp;
	double x;

	setbuf(stdout,(char *)NULL);

	params.radius = params.scaling = params.mass = 1.0;
	params.speed_max = 0.0;

	for (x=0;x<=2;x+=0.02) {
		rp.n_particles = 1000000; /* safe value */
		rp.radius = pow(10.0,x);
		SET_VEC(rp.axis_len,rp.radius,rp.radius,rp.radius);
		(void) printf("target R/r %f ",rp.radius/params.radius);
		rp.color = 3; /* dummy */
		calc_data(&params,&rp);
		rpuFree(&rp);
		}

	return 0;
	}
#else
int
main(int argc,char *argv[])
{
	PARAMS params;
	RUBBLE_PILE rp;

	setbuf(stdout,(char *)NULL);

#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
	/* enable explicit FPE trapping -- see #include above */
# if __APPLE__
	#error FPE handling not supported under macOS
# else
	feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
# endif
#else
	(void) fprintf(stderr,"WARNING: FPE trapping not enabled.\n");
#endif

	if (argc < 1 || argc > 2) {
		(void) fprintf(stderr,"Usage: %s [ rpg-par-file ]\n",argv[0]);
		return 1;
		}

    randSeedGenerator(randReadUrandom());

	get_params(argc == 2 ? argv[1] : DFLT_PARAMS_FILE,&params,&rp);
	calc_data(&params,&rp);
	write_data(&params,&rp);

	rpuFree(&rp);

	return 0;
	}
#endif /* !EFFICIENCY */

/* rpg.c */
