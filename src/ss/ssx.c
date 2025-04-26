/*
 ** ssx.c -- DCR 00-11-29
 ** =====
 ** Solar System data file transformer.
 */

#include <ss.h>
#include <unistd.h> /* for getpid() */
#include <ctype.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <boolean.h>
#include <vector.h>
#include <colors.h>
#include <matrix.h>

typedef struct {
	double total_mass;
	VECTOR bnd_min,bnd_max,com_pos,com_vel,ang_mom,vel_dsp;
	MATRIX inertia;
	int color;
	} PROPERTIES;

enum {NegativeOK,PositiveOnly}; /* for get_scaling() */

enum {FirstFile,MergeFile}; /* for read_data() */

static BOOLEAN
get_yn(const char *str,const char *dflt_str)
{
	enum {none,yes,no} dflt = none;

	int c;

	if (dflt_str && strlen(dflt_str)) {
		if (tolower(*dflt_str) == 'y') dflt = yes;
		else if (tolower(*dflt_str) == 'n') dflt = no;
		}

	do {
		(void) printf("%s [%s]? ",str,
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

static void
write_data(char *filename,SSDATA *d,int n,double t)
{
	SSIO ssio;
	SSHEAD h;
	int i;

	if (ssioOpen(filename,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",filename);
		return;
		}

	(void) printf("Number of particles to output = %i, time %g\n",n,t);

	h.time = t;
	h.n_data = n;
	h.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"Error writing header\n");
		return;
		}

	for (i=0;i<n;i++)
		if (ssioData(&ssio,&d[i])) {
			(void) fprintf(stderr,"Error writing data element %i\n",i);
			return;
			}

	(void) ssioClose(&ssio);
	}

static const char *
color_str(int color)
{
	switch (color) {
	case FIRST_GRAY:
	case BLACK:			return "black,reserved";
	case LAST_GRAY:
	case WHITE:			return "white";
	case RED:			return "red,reserved";
	case GREEN:			return "green";
	case BLUE:			return "blue";
	case YELLOW:		return "yellow,reserved";
	case MAGENTA:		return "magenta,reserved";
	case CYAN:			return "cyan,reserved";
	case GOLD:			return "gold";
	case PINK:			return "pink";
	case ORANGE:		return "orange";
	case KHAKI:			return "khaki,reserved";
	case VIOLET:		return "violet";
	case MAROON:		return "maroon";
	case AQUA:			return "aqua";
	case NAVY:			return "navy";
	default:			return "gray value";
		}
	}

static void
calc_inertia(SSDATA *d,int n,MATRIX inertia)
{
	int i,j,k;

	if (n == 1) {
		UNIT_MAT(inertia);
		SCALE_MAT(inertia,0.4*d->mass*SQ(d->radius));
		return;
		}

	ZERO_MAT(inertia);
	for (i=0;i<n;i++)
		for (j=0;j<N_DIM;j++)
			for (k=0;k<N_DIM;k++)
				inertia[j][k] +=
					d[i].mass*((j == k ? MAG_SQ(d[i].pos) : 0) -
							   d[i].pos[j]*d[i].pos[k]);
	}

static void
ss_analyze(SSDATA *data,int n_data,PROPERTIES *p)
{
	SSDATA *d;
	VECTOR r,v,l;
	int i,k,nc[NUM_COLORS],ncmax;

	assert(n_data > 0);

	p->total_mass = 0;
	ZERO_VEC(p->bnd_min);
	ZERO_VEC(p->bnd_max);
	ZERO_VEC(p->com_pos);
	ZERO_VEC(p->com_vel);
	ZERO_VEC(p->ang_mom);
	ZERO_VEC(p->vel_dsp);
	p->color = 0;

	for (i=0;i<NUM_COLORS;i++) nc[i] = 0;

	for (i=0;i<n_data;i++) {
		d = &data[i];
		p->total_mass += d->mass;
		COPY_VEC(d->pos,r);
		COPY_VEC(d->vel,v);
		if (i==0) {
			COPY_VEC(r,p->bnd_min);
			COPY_VEC(r,p->bnd_max);
			}
		else for (k=0;k<N_DIM;k++) {
			if (r[k] < p->bnd_min[k]) p->bnd_min[k] = r[k];
			if (r[k] > p->bnd_max[k]) p->bnd_max[k] = r[k];
			}
		SCALE_VEC(r,d->mass);
		ADD_VEC(p->com_pos,r,p->com_pos);
		SCALE_VEC(v,d->mass);
		ADD_VEC(p->com_vel,v,p->com_vel);
		CROSS(d->pos,d->vel,l);
		SCALE_VEC(l,d->mass);
		ADD_VEC(p->ang_mom,l,p->ang_mom);
		++nc[(int) d->color];
		}

	if (p->total_mass == 0) {
		(void) printf("WARNING: total mass = 0...will divide by N\n");
		p->total_mass = n_data;
		}

	NORM_VEC(p->com_pos,p->total_mass);
	NORM_VEC(p->com_vel,p->total_mass);
	NORM_VEC(p->ang_mom,p->total_mass);

	for (i=0;i<n_data;i++) {
		d = &data[i];
		SUB_VEC(d->vel,p->com_vel,v);
		for (k=0;k<N_DIM;k++)
			p->vel_dsp[k] += d->mass*SQ(v[k]);
		}

	for (k=0;k<N_DIM;k++)
		p->vel_dsp[k] = sqrt(p->vel_dsp[k]/p->total_mass);

	ncmax = 0;

	for (i=0;i<NUM_COLORS;i++)
		if (nc[i] > ncmax) {
			p->color = i;
			ncmax = nc[i];
			}

	calc_inertia(data,n_data,p->inertia);
	}

static void
scale_mass(SSDATA *d,int n,double f)
{
	int i;

	for (i=0;i<n;i++)
		d[i].mass *= f;
	}

static void
adj_com_pos(SSDATA *d,int n,PROPERTIES *p,VECTOR v)
{
	int i;

	for (i=0;i<n;i++) {
		SUB_VEC(d[i].pos,p->com_pos,d[i].pos);
		ADD_VEC(d[i].pos,v,d[i].pos);
		}
	}

static void
adj_com_vel(SSDATA *d,int n,PROPERTIES *p,VECTOR v)
{
	int i;

	for (i=0;i<n;i++) {
		SUB_VEC(d[i].vel,p->com_vel,d[i].vel);
		ADD_VEC(d[i].vel,v,d[i].vel);
		}
	}

static void
adj_ang_mom(SSDATA *d,int n,PROPERTIES *p,VECTOR v)
{
	VECTOR u,w;
	int i;

	matrixInvert(p->inertia);
	SUB_VEC(v,p->ang_mom,v);
	Transform(p->inertia,v,u);
	SCALE_VEC(u,p->total_mass);
	for (i=0;i<n;i++) {
		CROSS(u,d[i].pos,w);
		ADD_VEC(d[i].vel,w,d[i].vel);
		}
	}

static void
scale_vel_dsp(SSDATA *d,int n,PROPERTIES *p,VECTOR v)
{
	int i,k;

	for (i=0;i<n;i++) {
		SUB_VEC(d[i].vel,p->com_vel,d[i].vel);
		for (k=0;k<N_DIM;k++)
			d[i].vel[k] *= v[k];
		}

	for (i=0;i<n;i++) {
		ADD_VEC(d[i].vel,p->com_vel,d[i].vel);
		}
	}

static void
change_color(SSDATA *d,int n,int c)
{
	int i;

	for (i=0;i<n;i++)
		d[i].color = c;
	}

static void scale_masses(SSDATA *d,int n,double f)
{
	int i;

	if (f < 0.0) {
		f = -f;
		for (i=0;i<n;i++)
		  d[i].mass = f;
	}
	else
	  for (i=0;i<n;i++)
		d[i].mass *= f;
	}

static void scale_radii(SSDATA *d,int n,double f)
{
	int i;

	if (f < 0.0) {
		f = -f;
		for (i=0;i<n;i++)
		  d[i].radius = f;
	}
	else
	  for (i=0;i<n;i++)
		d[i].radius *= f;
	}

static void
get_scaling(double *f,int option)
{
	(void) printf("Enter scaling factor");
	if (option == NegativeOK) printf(" (-ve for abs val)");
	printf(": ");
	do {
		(void) scanf("%lf",f);
		(void) getchar();
		} while (option == PositiveOnly && *f <= 0);
	}

static void
get_components(VECTOR v)
{
	(void) printf("Enter vector components (x y z): ");
	(void) scanf("%lf%lf%lf",&v[X],&v[Y],&v[Z]);
	(void) getchar();
	}

static void
get_component_scaling(VECTOR v)
{
	(void) printf("Enter component scaling factors (x y z): ");
	(void) scanf("%lf%lf%lf",&v[X],&v[Y],&v[Z]);
	(void) getchar();
}	

static void
process(SSDATA *d,int n,double *t)
{
	PROPERTIES p;
	int choice;

	enum {Next,Time,Mass,Bounds,ComPos,ComVel,AngMom,VelDsp,Color,Units,
			  Offsets,Masses,Radii,End};

	while (/*CONSTCOND*/1) {

		ss_analyze(d,n,&p);

		(void) printf("%2i. Time = %g\n",Time,*t);
		(void) printf("%2i. Total mass = %g\n",Mass,p.total_mass);
		(void) printf("%2i. Bounds: x=[%g,%g]\n"
					  "            y=[%g,%g]\n"
					  "            z=[%g,%g]\n",Bounds,
					  p.bnd_min[X],p.bnd_max[X],
					  p.bnd_min[Y],p.bnd_max[Y],
					  p.bnd_min[Z],p.bnd_max[Z]);
		(void) printf("%2i. Centre-of-mass position = %g %g %g\n",ComPos,
					  p.com_pos[X],p.com_pos[Y],p.com_pos[Z]);
		(void) printf("%2i. Centre-of-mass velocity = %g %g %g\n",ComVel,
					  p.com_vel[X],p.com_vel[Y],p.com_vel[Z]);
		(void) printf("%2i. Specific angular momentum = %g %g %g\n",AngMom,
					  p.ang_mom[X],p.ang_mom[Y],p.ang_mom[Z]);
		(void) printf("%2i. Velocity dispersion = %g %g %g\n",VelDsp,
					  p.vel_dsp[X],p.vel_dsp[Y],p.vel_dsp[Z]);
		(void) printf("%2i. Dominant color = %i (%s)\n",Color,p.color,
					  color_str(p.color));
		(void) printf("%2i. Units\n",Units);
		(void) printf("%2i. Offsets\n",Offsets);
		(void) printf("%2i. Particle masses\n",Masses);
		(void) printf("%2i. Particle radii\n",Radii);

		do {
			(void) printf("Enter number to change (or 0 to continue): ");
			(void) scanf("%i",&choice);
			} while (choice < Next || choice >= End);
		
		(void) getchar();

		if (choice == Next) return;

		switch(choice) {
		case Time:
			do {
				(void) printf("Enter new time: ");
				(void) scanf("%lf",t);
				(void) getchar();
				} while (*t < 0);
			break;
		case Mass:
			{
			double f;
			do get_scaling(&f,NegativeOK);
			while (f == 0);
			if (f < 0) f = -f/p.total_mass;
			scale_mass(d,n,f);
			break;
			}
		case Bounds:
			{
			double f,min,max;
			int i,choice2;
			do {
				do {
					(void) printf("%i. Change x bounds (now [%g,%g])\n",X + 1,
								  p.bnd_min[X],p.bnd_max[X]);
					(void) printf("%i. Change y bounds (now [%g,%g])\n",Y + 1,
								  p.bnd_min[Y],p.bnd_max[Y]);
					(void) printf("%i. Change z bounds (now [%g,%g])\n",Z + 1,
								  p.bnd_min[Z],p.bnd_max[Z]);
					(void) printf("Your choice (or 0 when done): ");
					(void) scanf("%i",&choice2);
					(void) getchar();
					} while (choice2 < 0 || choice2 > N_DIM);
				if (choice2 == 0) break;
				--choice2; /* put back in range [X,Z] */
				if (p.bnd_min[choice2] == p.bnd_max[choice2]) {
					(void) printf("Chosen dimension is degenerate\n");
					continue;
					}
				do {
					(void) printf("Enter new bounds (min max): ");
					(void) scanf("%lf%lf",&min,&max);
					(void) getchar();
					} while (min > max);
				if (min == max &&
					get_yn("Zero velocities for this component","y"))
					for (i=0;i<n;i++)
						d[i].vel[choice2] = 0;
				f = (max - min)/(p.bnd_max[choice2] - p.bnd_min[choice2]);
				for (i=0;i<n;i++)
					d[i].pos[choice2] =
						(d[i].pos[choice2] - p.bnd_min[choice2])*f + min;
				p.bnd_min[choice2] = min;
				p.bnd_max[choice2] = max;
				} while (/*CONSTCOND*/1);
			break;
			}
		case ComPos:
			{
			VECTOR v;
			if (MAG(p.com_pos) && get_yn("Scale the magnitude","y")) {
				double f;
				get_scaling(&f,NegativeOK);
				COPY_VEC(p.com_pos,v);
				if (f < 0) f = -f/MAG(v);
				SCALE_VEC(v,f);
				}
			else get_components(v);
			adj_com_pos(d,n,&p,v);
			break;
			}
		case ComVel:
			{
			VECTOR v;
			if (MAG(p.com_vel) && get_yn("Scale the magnitude","y")) {
				double f;
				get_scaling(&f,NegativeOK);
				COPY_VEC(p.com_vel,v);
				if (f < 0) f = -f/MAG(v);
				SCALE_VEC(v,f);
				}
			else get_components(v);
			adj_com_vel(d,n,&p,v);
			break;
			}
		case AngMom:
			{
			VECTOR v;
			(void) printf("NOTE: specific angular momentum is measured with\n"
						  "respect to fixed space frame centred at (0,0,0)\n"
						  "and does not take particle spins into account\n");
			if (MAG(p.ang_mom) && get_yn("Scale the magnitude","y")) {
				double f;
				get_scaling(&f,NegativeOK);
				COPY_VEC(p.ang_mom,v);
				if (f < 0) f = -f/MAG(p.ang_mom);
				SCALE_VEC(v,f);
				}
			else if (get_yn("Scale the components","y")) {
				VECTOR u;
				int k;
				get_component_scaling(u);
				for (k=0;k<N_DIM;k++)
					v[k] = u[k]*p.ang_mom[k];
				}
			else get_components(v);
			adj_ang_mom(d,n,&p,v);
			break;
			}
		case VelDsp:
			{
			VECTOR v;
			(void) printf("NOTE: velocity dispersion is context dependent,\n"
						  "for now relative ONLY to center-of-mass velocity,\n"
						  "i.e. without considering bulk rotation or shear\n");
			if (!MAG(p.vel_dsp)) {
				(void) printf("Zero velocity dispersion -- cannot adjust\n");
				break;
				}
			if (get_yn("Scale the magnitude","y")) {
				double f;
				get_scaling(&f,NegativeOK);
				if (f < 0) f = -f/MAG(p.vel_dsp);
				SET_VEC(v,f,f,f);
				}
			else get_component_scaling(v);
			scale_vel_dsp(d,n,&p,v);
			break;
			}
		case Color:
			{
			int c;
			(void) printf("Color scheme:\n");
			for (c=BLACK;c<FIRST_GRAY;c++)
				(void) printf("%2i. %s\n",c,color_str(c));
			(void) printf("[values from %i to %i are levels of gray]\n",
						  FIRST_GRAY,LAST_GRAY);
			do {
				(void) printf("Enter new color: ");
				(void) scanf("%i",&c);
				(void) getchar();
				} while (c < 0 || c >= NUM_COLORS);
			change_color(d,n,c);
			break;
			}
		case Units:
			{
			enum {N,M,L,T,V,E};
			double f;
			int i,choice2;
			(void) printf("NOTE: It is up to you to ensure dimensions are\n"
						  "internally consistent. pkdgrav assumes G == 1.\n");
			do {
				do {
					(void) printf("%i. Mass (particle masses)\n",M);
					(void) printf("%i. Length (particle radii, pos'ns)\n",L);
					(void) printf("%i. Time (time,particle spins)\n",T);
					(void) printf("%i. Velocity (particle velocities)\n",V);
					(void) printf("Select dimension to scale "
								  "(or 0 when done): ");
					(void) scanf("%i",&choice2);
					(void) getchar();
					} while (choice2 < N || choice2 >= E);
				if (choice2 == N) break;
				switch (choice2) {
				case M:
					(void) printf("M_Sun     = 1.9891e30 kg\n"
								  "M_Earth   = 5.9742e24 kg\n"
								  "M_Jupiter = 1.8992e27 kg\n"
								  "M_Saturn  = 5.6864e26 kg\n");
					get_scaling(&f,PositiveOnly);
					for (i=0;i<n;i++)
						d[i].mass *= f;
					break;
				case L:
					(void) printf("1 AU    = 1.49597892e11 m\n"
								  "R_Earth = 6.37814e6 m\n");
					get_scaling(&f,PositiveOnly);
					for (i=0;i<n;i++) {
						d[i].radius *= f;
						SCALE_VEC(d[i].pos,f);
						}
					break;
				case T:
					(void) printf("1 yr        = 3.15576e7 s\n"
								  "1 yr / 2 pi = 5.02255e6 s\n");
					get_scaling(&f,PositiveOnly);
					*t *= f;
					for (i=0;i<n;i++)
						NORM_VEC(d[i].spin,f);
					break;
				case V:
					(void) printf("V_Earth = 2.97852586e4 m/s\n");
					get_scaling(&f,PositiveOnly);
					for (i=0;i<n;i++)
						SCALE_VEC(d[i].vel,f);
					break;
				default:
					assert(0);
					}
				} while (/*CONSTCOND*/1);
			break;
			}
		case Offsets:
			{
			VECTOR v;
			int i;
			(void) printf("POSITION OFFSET (0 0 0 for none)...\n");
			get_components(v);
			for (i=0;i<n;i++)
				ADD_VEC(d[i].pos,v,d[i].pos);
			(void) printf("VELOCITY OFFSET (0 0 0 for none)...\n");
			get_components(v);
			for (i=0;i<n;i++)
				ADD_VEC(d[i].vel,v,d[i].vel);
			break;
			}
		case Masses:
			{
			double f;
			do get_scaling(&f,NegativeOK);
			while (f == 0);
			scale_masses(d,n,f);
			break;
			}
		case Radii:
			{
			double f;
			do get_scaling(&f,NegativeOK);
			while (f == 0);
			scale_radii(d,n,f);
			break;
			}
		default:
			assert(0);
			}
		}
	}

static int
read_data(char *filename,SSDATA **d,int *n,double *t,int mode)
{
	SSIO ssio;
	SSHEAD h;
	int i;

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h) || h.n_data < 0) {
		(void) fprintf(stderr,"Corrupt header\n");
		(void) ssioClose(&ssio);
		return 1;
		}

	if (h.n_data == 0) {
		(void) fprintf(stderr,"No data found!");
		(void) ssioClose(&ssio);
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

	if (mode == FirstFile) {
		*d = NULL;
		*n = 0;
		*t = 0.0;
		}

	(void) printf("Number of particles = %i, time %g\n",h.n_data,h.time);

	if (mode == MergeFile && h.time != *t && !get_yn("Keep old time","y")) {
		if (get_yn("Use time of merged file","y"))
			*t = h.time;
		else {
			(void) printf("Enter new time (old time was %g): ",*t);
			(void) scanf("%lf",t);
			}
		}

	*d = (SSDATA *) realloc(*d,(size_t) (*n + h.n_data)*sizeof(SSDATA));

	assert(*d != NULL);

	for (i=*n;i<*n + h.n_data;i++)
		if (ssioData(&ssio,&((*d)[i]))) {
			(void) fprintf(stderr,"Corrupt data\n");
			(void) ssioClose(&ssio);
			return 1;
			}

	(void) ssioClose(&ssio);

	*n += h.n_data;

	return 0;
	}

int
main(int argc,char *argv[])
{
	SSDATA *d;
	double t;
	int n;
	FILE *fp;
	char file[MAXPATHLEN],first_file[MAXPATHLEN];

	setbuf(stdout,(char *)NULL);

	srand((unsigned int) getpid());

	if (argc > 1) {
		(void) fprintf(stderr,"%s takes no arguments\n",argv[0]);
		return 1;
		}

	do {
		file[0] = '\0';
		(void) printf("Enter file to process: ");
		(void) fgets(file,MAXPATHLEN,stdin);
		assert(strlen(file));
		file[strlen(file) - 1] = '\0'; /* get rid of newline at end */
		} while (!strlen(file) || read_data(file,&d,&n,&t,FirstFile));

	process(d,n,&t);

	(void) strcpy(first_file,file);

	do {
		file[0] = '\0';
		(void) printf("Enter file to merge [or ENTER if done]: ");
		(void) fgets(file,MAXPATHLEN,stdin);
		assert(strlen(file));
		file[strlen(file) - 1] = '\0';
		if (!strlen(file)) break;
		if (!read_data(file,&d,&n,&t,MergeFile))
			(void) printf("Total number of particles read so far = %i\n",n);
		process(d,n,&t);
		} while (/*CONSTCOND*/1);

	do {
		(void) printf("Output file [default %s]: ",first_file);
		(void) fgets(file,MAXPATHLEN,stdin);
		assert(strlen(file));
		file[strlen(file) - 1] = '\0';
		if (!strlen(file)) (void) strcpy(file,first_file);
		if (!(fp = fopen(file,"r"))) break;
		(void) fclose(fp);
		if (get_yn("Output file already exists...overwrite","n")) break;
		} while (/*CONSTCOND*/1);

	write_data(file,d,n,t);

	if (d) free((void *) d);

	(void) printf("Done!\n");

	return 0;
	}

/* ssx.c */
