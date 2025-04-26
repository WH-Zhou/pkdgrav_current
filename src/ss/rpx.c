/*
 ** rpx.c -- DCR 98-09-22
 ** =====
 ** Rubble pile transformer.
 **
 ** DEBUG 4/1/13: note that rpx does not do anything with the original
 ** indices of the particles in each ss file that is read, so the
 ** merged file will likely have duplicated original indices.  One way
 ** to fix this would be to add (to each original index of each
 ** particle in each new rubble pile read in) the number of particles
 ** read from the previous rubble pile.  This could be the default
 ** behavior, with a flag to disable it.
 */

#include <rpu.h>
#include <stdlib.h>    /* for rand() & srand() */
#include <unistd.h>    /* for getopt() & getpid() */
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#define LOG_FILE "rpx.log"

#define MAX_NUM_FILES 128

enum {Descending=-1,None=0,Ascending=+1}; /* for sort_order */

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

#pragma GCC diagnostic push /* require GCC 4.6 */
#pragma GCC diagnostic ignored "-Wcast-qual"
static int
ascend(const void *a,const void *b)
{
	/*
	** a & b are elements to be compared; in our case, they are
	** generic pointers to SSDATA pointers to SSDATA elements, hence
	** the awkward unpacking below...
	 */

	double Ra = (*((SSDATA **) a))->radius;
	double Rb = (*((SSDATA **) b))->radius;

	return (Ra < Rb ? -1 : Ra > Rb ? 1 : 0);
	}

static int
descend(const void *a,const void *b)
{
	double Ra = (*((SSDATA **) a))->radius;
	double Rb = (*((SSDATA **) b))->radius;
	
	return (Ra > Rb ? -1 : Ra < Rb ? 1 : 0);
	}
#pragma GCC diagnostic pop

static int
write_data(char *filename,RUBBLE_PILE *rp,int n,double time,BOOLEAN reindex,int sort_order)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA **data;
	int n_data,i,j,k;

	if (ssioOpen(filename,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\".\n",filename);
		return 1;
		}

	n_data = 0;
	for (i=0;i<n;i++)
		n_data += rp[i].n_particles;

	printf("Number of particles to output = %i.\n",n_data);

	h.time = time;
	h.n_data = n_data;
	h.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&h)) {
		fprintf(stderr,"Error writing header.\n");
		return 1;
		}

	/* to enable sorting, make reference data structure to all data */

	data = (SSDATA **) malloc((size_t) n_data*sizeof(SSDATA *));
	if (data == NULL) {
		fprintf(stderr,"Unable to allocate memory for data referencing.\n");
		return 1;
		}
	for (i=k=0;i<n;i++)
		for (j=0;j<rp[i].n_particles;j++,k++)
			data[k] = &(rp[i].data[j]); /* save memory: store data pointers */
	assert(k == n_data); /* paranoia check */
	
	switch (sort_order) {
	case Ascending:
		qsort((void *) data,(size_t) n_data,sizeof(SSDATA *),ascend);
		break;
	case Descending:
		qsort((void *) data,(size_t) n_data,sizeof(SSDATA *),descend);
		break;
	case None:
		break;
	default:
		assert(0); /* should not get here */
		}

	if (reindex)
		for (i=0;i<n_data;i++)
			data[i]->org_idx = i;
	
	for (i=0;i<n_data;i++)
		if (ssioData(&ssio,data[i])) {
			fprintf(stderr,"Error writing data.\n");
			return 1;
			}

	/* old
	for (i=0;i<n;i++)
		for (j=0;j<rp[i].n_particles;j++)
			if (ssioData(&ssio,&rp[i].data[j])) {
				fprintf(stderr,"Error writing data.\n");
				return 1;
				}
	*/

	free((void *) data);
	ssioClose(&ssio);

	return 0;
	}

static void
show_encounter(RUBBLE_PILE *rp,int n,BOOLEAN bUseSimUnits)
{
	VECTOR rel_pos,rel_vel,ang_mom;
	double m,sr,r,v,eb,a,e,q,Q,vq,vQ;

	assert(n == 2);

	m = rp[0].mass + rp[1].mass;
	sr = rp[0].radius + rp[1].radius;

	SUB_VEC(rp[1].pos,rp[0].pos,rel_pos);
	SUB_VEC(rp[1].vel,rp[0].vel,rel_vel);
	CROSS(rel_pos,rel_vel,ang_mom);

	r = MAG(rel_pos);
	v = MAG(rel_vel);

	if (r == 0) eb = 0;
	else eb = 0.5*SQ(v) - m/r;
	if (eb == 0) {
		a = 0;
		e = 1;
		}
	else {
		a = -0.5*m/eb;
		e = sqrt(1 - MAG_SQ(ang_mom)/(a*m));
		}
	q = (1 - e)*a;
	Q = (1 + e)*a;

	(void) printf("2-body encounter parameters:\n");
	(void) printf("Binding energy per unit reduced mass = %g sim units\n",eb);
	(void) printf("Semi-major axis a = ");
	if (bUseSimUnits) (void) printf("%g AU",a);
	else (void) printf("%g km",a*0.001*L_SCALE);
	(void) printf("\n");
	(void) printf("Eccentricity e = %g\n",e);
	(void) printf("Close-approach distance q = ");
	if (bUseSimUnits) (void) printf("%g AU",q);
	else (void) printf("%g km",q*0.001*L_SCALE);
	if (q <= sr) (void) printf(" (impact possible)");
	(void) printf("\n");
	(void) printf("Initial relative speed v = ");
	if (bUseSimUnits) (void) printf("%g x 30 km/s",v);
	else (void) printf("%g m/s",v*V_SCALE);
	(void) printf("\n");
	(void) printf("Relative speed at ");
	if (eb < 0) {
		vQ = sqrt(2*(eb + m/Q));
		(void) printf("apoapse vQ = ");
		}
	else {
		vQ = sqrt(2*eb);
		(void) printf("infinity v_inf = ");
		}
	if (bUseSimUnits) (void) printf("%g x 30 km/s",vQ);
	else (void) printf("%g m/s",vQ*V_SCALE);
	(void) printf("\n");
	(void) printf("Relative speed at ");
	if (q <= sr) {
		if (r > sr) q = sr;
		(void) printf("impact (approx)");
		}
	else (void) printf("periapse");
	(void) printf(" vq = ");
	if (q < sr) vq = v;
	else vq = sqrt(2*(eb + m/q));	
	if (bUseSimUnits) (void) printf("%g x 30 km/s",vq);
	else (void) printf("%g m/s",vq*V_SCALE);
	(void) printf("\n");
	(void) printf("Relative orb. ang. mom. per unit reduced mass h = "
				  "(%g,%g,%g) sim units\n",ang_mom[X],ang_mom[Y],ang_mom[Z]);
	}

static const char *
color_str(int color)
{
	switch (color) {
	case FIRST_GRAY:
	case BLACK:		return "black";
	case LAST_GRAY:
	case WHITE:		return "white";
	case RED:		return "red";
	case GREEN:		return "green";
	case BLUE:		return "blue";
	case YELLOW:	return "yellow";
	case MAGENTA:	return "magenta";
	case CYAN:		return "cyan";
	case GOLD:		return "gold";
	case PINK:		return "pink";
	case ORANGE:	return "orange";
	case KHAKI:		return "khaki";
	case VIOLET:	return "violet";
	case MAROON:	return "maroon";
	case AQUA:		return "aqua";
	case NAVY:		return "navy";
	default: 		return "gray value";
		}
	}

static int
process(char *filename,RUBBLE_PILE *rp,BOOLEAN sim_units,double *time,int *reindex)
{
	enum {next,mass,radius,density,pos,vel,orient,spin,color,agg_id,par_id};

	SSIO ssio;
	SSHEAD h;
	int i,choice;

	assert(rp != NULL);

	*time = 0.0;

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

	rp->n_particles = h.n_data;
	*time = h.time;
	(void) printf("Number of particles = %i (time %g)\n",rp->n_particles,*time);
	rpuMalloc(rp);

	for (i=0;i<rp->n_particles;i++)
		if (ssioData(&ssio,&rp->data[i])) {
			(void) fprintf(stderr,"Corrupt data\n");
			(void) ssioClose(&ssio);
			return 1;
			}

	(void) ssioClose(&ssio);

	while (/*CONSTCOND*/1) {

		rpuAnalyze(rp);

		(void) printf("%i. Total mass = ",mass);
		if (sim_units) (void) printf("%g M_sun",rp->mass);
		else (void) printf("%g kg",rp->mass*M_SCALE);
		(void) printf("\n");

		(void) printf("%i. Bulk radius = ",radius);
		if (sim_units) (void) printf("%g AU",rp->radius);
		else (void) printf("%g km",rp->radius*0.001*L_SCALE);
		(void) printf("\n");

		(void) printf("   [Bulk semi-axes: ");
		if (sim_units) (void) printf("%g %g %g AU",
									 rp->axis_len[rp->axis_ord[X]],
									 rp->axis_len[rp->axis_ord[Y]],
									 rp->axis_len[rp->axis_ord[Z]]);
		else (void) printf("%g %g %g km",
						   rp->axis_len[rp->axis_ord[X]]*0.001*L_SCALE,
						   rp->axis_len[rp->axis_ord[Y]]*0.001*L_SCALE,
						   rp->axis_len[rp->axis_ord[Z]]*0.001*L_SCALE);
		(void) printf("]\n");

		printf("   [DEEVE semi-axes: ");
		if (sim_units)
		    printf("%g %g %g AU",
			   rp->DEEVE[0],
			   rp->DEEVE[1],
			   rp->DEEVE[2]);
		else
		    printf("%g %g %g km",
			   rp->DEEVE[0]*0.001*L_SCALE,
			   rp->DEEVE[1]*0.001*L_SCALE,
			   rp->DEEVE[2]*0.001*L_SCALE);
		printf("]\n");

		printf("   [DEEVE equivalent radius = ");
		if (sim_units)
		    printf("%g AU",
			   pow(rp->DEEVE[0] * rp->DEEVE[1] * rp->DEEVE[2], 1.0 / 3.0));
		else
		    printf("%g km",
			   pow(rp->DEEVE[0] * rp->DEEVE[1] * rp->DEEVE[2], 1.0 / 3.0) * 0.001 * L_SCALE);
		printf("]\n");

		(void) printf("%i. Bulk density = ",density);
		if (sim_units) (void) printf("%g M_sun/AU^3",rp->density);
		else (void) printf("%g g/cc",rp->density*0.001*D_SCALE);
		(void) printf("\n");

		(void) printf("   [DEEVE bulk density = ");
		if (sim_units)
		    printf("%g M_sun/AU^3",
			   rp->mass / (4.0 / 3.0 * M_PI * rp->DEEVE[0] * rp->DEEVE[1] * rp->DEEVE[2]));
		else
		    printf("%g g/cc",
			   rp->mass / (4.0 / 3.0 * M_PI * rp->DEEVE[0] * rp->DEEVE[1] * rp->DEEVE[2]) * 0.001 * D_SCALE);
		printf("]\n");

		(void) printf("%i. Centre-of-mass position = ",pos);
		if (sim_units) (void) printf("%g %g %g AU",rp->pos[X],rp->pos[Y],
									 rp->pos[Z]);
		else (void) printf("%.2f %.2f %.2f km",rp->pos[X]*0.001*L_SCALE,
						   rp->pos[Y]*0.001*L_SCALE,rp->pos[Z]*0.001*L_SCALE);
		(void) printf("\n");

		(void) printf("%i. Centre-of-mass velocity = ",vel);
		if (sim_units) (void) printf("%g %g %g x 30 km/s",rp->vel[X],rp->vel[Y],
									 rp->vel[Z]);
		else (void) printf("%.2f %.2f %.2f m/s",rp->vel[X]*V_SCALE,
						   rp->vel[Y]*V_SCALE,rp->vel[Z]*V_SCALE);
		(void) printf("\n");

		(void) printf("%i. Orientation: a1 = %6.6f %6.6f %6.6f\n",orient,
					  rp->axes[rp->axis_ord[X]][X],
					  rp->axes[rp->axis_ord[X]][Y],
					  rp->axes[rp->axis_ord[X]][Z]);
		(void) printf("                a2 = %6.6f %6.6f %6.6f\n",
					  rp->axes[rp->axis_ord[Y]][X],
					  rp->axes[rp->axis_ord[Y]][Y],
					  rp->axes[rp->axis_ord[Y]][Z]);
		(void) printf("                a3 = %6.6f %6.6f %6.6f\n",
					  rp->axes[rp->axis_ord[Z]][X],
					  rp->axes[rp->axis_ord[Z]][Y],
					  rp->axes[rp->axis_ord[Z]][Z]);

		(void) printf("%i. Spin = ",spin);
		if (sim_units) (void) printf("%g %g %g x 2pi rad/yr",
									 rp->spin[X],rp->spin[Y],rp->spin[Z]);
		else {
			double scale = 3600/(TWO_PI*T_SCALE),w;
			(void) printf("%.2f %.2f %.2f 1/h (",rp->spin[X]*scale,
						  rp->spin[Y]*scale,rp->spin[Z]*scale);
			w = MAG(rp->spin);
			if (w) (void) printf("period %g h",1/(w*scale));
			else (void) printf("no spin");
			(void) printf(")");
			}
		(void) printf("\n");

		(void) printf("   [Ang mom = ");
		if (sim_units) (void) printf("%g %g %g (sys units)",rp->ang_mom[X],
									 rp->ang_mom[Y],rp->ang_mom[Z]);
		else {
			double scale = M_SCALE*SQ(L_SCALE)/T_SCALE;
			(void) printf("%.5e %.5e %.5e N m/s",rp->ang_mom[X]*scale,
						  rp->ang_mom[Y]*scale,rp->ang_mom[Z]*scale);
			}
		(void) printf("]\n");

		(void) printf("   [Effective spin = ");
		if (sim_units) (void) printf("%g x 2pi rad/yr",rp->eff_spin);
		else {
			double scale = 3600/(TWO_PI*T_SCALE);
			(void) printf("%.2f 1/h (",rp->eff_spin*scale);
			if (rp->eff_spin) printf("period %g h",1/(rp->eff_spin*scale));
			else (void) printf("no spin");
			(void) printf(")");
			}
		(void) printf("]\n");

		(void) printf("   [Rotation index = %.2f (",rp->rot_idx);
		if (rp->eff_spin == 0.0) (void) printf("undefined");
		else if (rp->rot_idx == 1.0) (void) printf("unif rot about max moment");
		else if (rp->rot_idx < 1.0 && rp->rot_idx > 0.0) (void) printf("SAM");
		else if (rp->rot_idx == 0.0) (void) printf("unif rot about mid moment");
		else if (rp->rot_idx < 0.0 && rp->rot_idx > -1.0) (void) printf("LAM");
		else if (rp->rot_idx == -1.0) (void) printf("unif rot about min moment");
		else assert(0);
		(void) printf(")]\n");

		(void) printf("%i. Color = %i (%s)\n",color,(int) rp->color,
					  color_str(rp->color));

		(void) printf("%i. Aggregate ID = ",agg_id);
		if (rp->agg_id < 0)
			(void) printf("N/A");
		else
			(void) printf("%i",(int) rp->agg_id);
		(void) printf("\n");

		(void) printf("%i. Particle ID\n",par_id);

		do {
			(void) printf("Enter number to change (or 0 to continue): ");
			(void) scanf("%i",&choice);
			(void) getchar();
			} while (choice < next || choice > par_id);
		
		if (choice == next) return 0;

		switch(choice) {
		case mass:
			{
			double f;
			BOOLEAN const_den = get_yn("Keep bulk density constant","n"),proceed=TRUE;
			do {
				(void) printf("Enter mass scaling factor (-ve ==> abs val): ");
				(void) scanf("%lf",&f);
				if (f == 0.0) {
					getchar();
					proceed = get_yn("WARNING: This will set all particle masses to zero...continue","n");
					}
				} while (!proceed);
			if (f < 0) {
				if (!sim_units) f /= M_SCALE;
				f = -f/rp->mass;
				}
			rpuScaleMass(rp,f);
			if (const_den) rpuScaleRadius(rp,pow(f,1.0/3),FALSE);
			break;
			}
		case radius:
			{
			double f;
			BOOLEAN just_particles = get_yn("Just scale particles","n"),const_den = FALSE;
			if (!just_particles)
				const_den = get_yn("Keep bulk density constant","n");
			do {
				(void) printf("Enter radius scaling factor "
							  "(-ve ==> abs val): ");
				(void) scanf("%lf",&f);
				} while (f == 0);
			getchar();
			if (f < 0) {
				if (!sim_units) f *= 1000/L_SCALE;
				if (!just_particles)
					f = -f/rp->radius;
				}
			rpuScaleRadius(rp,f,just_particles);
			if (just_particles)
				rpuCalcRadius(rp); /* because outer edge may have changed */
			if (const_den) rpuScaleMass(rp,CUBE(f));
			break;
			}
		case density:
			{
			double f;
			BOOLEAN const_radius = get_yn("Keep bulk radius constant","y");
			do {
				(void) printf("Enter density scaling factor (-ve ==> abs val): ");
				(void) scanf("%lf",&f);
				} while (f == 0);
			getchar();
			if (f < 0) {
				if (!sim_units) f *= 1000/D_SCALE;
				f = -f/rp->density;
				}
			if (const_radius) rpuScaleMass(rp,f);
			else rpuScaleRadius(rp,pow(f,-1.0/3),FALSE);
			break;
			}
		case pos:
		    {
			BOOLEAN absolute = get_yn("Specify absolute position","y");
			if (absolute) {
				SCALE_VEC(rp->pos,-1.0);
				rpuApplyPos(rp); /* reset COM to (0,0,0) first */
				(void) printf("Enter new position [x y z in ");
				}
			else
				(void) printf("Enter position offset [x y z in ");
			if (sim_units) (void) printf("AU");
			else (void) printf("km");
			(void) printf("]: ");
			(void) scanf("%lf%lf%lf",&rp->pos[X],&rp->pos[Y],&rp->pos[Z]);
			(void) getchar();
			if (!sim_units) NORM_VEC(rp->pos,0.001*L_SCALE);
			rpuApplyPos(rp);
			break;
			}
		case vel:
		    {
			BOOLEAN absolute = get_yn("Specify absolute velocity","y");
			if (absolute) {
				SCALE_VEC(rp->vel,-1.0);
				rpuApplyVel(rp);
				(void) printf("Enter new velocity [vx vy vz in ");
				}
			else
				(void) printf("Enter velocity offset [vx vy vz in ");
			if (sim_units) (void) printf("units of 30 km/s");
			else (void) printf("m/s");
			(void) printf("]: ");
			(void) scanf("%lf%lf%lf",&rp->vel[X],&rp->vel[Y],&rp->vel[Z]);
			(void) getchar();
			if (!sim_units) NORM_VEC(rp->vel,V_SCALE);
			rpuApplyVel(rp);
			break;
			}
		case orient:
			{
			enum {x=1,y,z};
			MATRIX rot;
			double angle;
			BOOLEAN align,body,rndm;
			int choice2;
			align = get_yn("Align body axes with coordinate axes","n");
			if (align) {
				MATRIX m;
				COPY_VEC(rp->axes[MAJOR(rp)],m[X]);
				COPY_VEC(rp->axes[INTER(rp)],m[Y]);
				COPY_VEC(rp->axes[MINOR(rp)],m[Z]);
				rpuRotate(rp,m,FALSE);
				break;
				}
			body = get_yn("Use body axes","n");
			(void) printf("%i. Rotate about %s axis\n",x,body?"major":"x");
			(void) printf("%i. Rotate about %s axis\n",y,body?"intermediate":"y");
			(void) printf("%i. Rotate about %s axis\n",z,body?"minor":"z");
			do {
				(void) printf("Enter choice: ");
				(void) scanf("%i",&choice2);
				} while (choice2 < x || choice2 > z);
			--choice2; /* to conform with X,Y,Z macros */
			(void) getchar();
			rndm = get_yn("Use random angle","n");
			if (rndm) angle = TWO_PI*rand()/RAND_MAX;
			else {
				(void) printf("Enter rotation angle in ");
				if (sim_units) (void) printf("radians");
				else (void) printf("degrees");
				(void) printf(" (-ve=clockwise): ");
				(void) scanf("%lf",&angle);
				(void) getchar();
				if (!sim_units) angle *= DEG_TO_RAD;
				}
			UNIT_MAT(rot);
			if (body) choice2 = rp->axis_ord[choice2];
			switch (choice2) {
			case X:
				rot[Y][Y] = cos(angle); rot[Y][Z] = -sin(angle);
				rot[Z][Y] = -rot[Y][Z]; rot[Z][Z] = rot[Y][Y];
				break;
			case Y:
				rot[X][X] = cos(angle); rot[X][Z] = sin(angle);
				rot[Z][X] = -rot[X][Z]; rot[Z][Z] = rot[X][X];
				break;
			case Z:
				rot[X][X] = cos(angle); rot[X][Y] = -sin(angle);
				rot[Y][X] = -rot[X][Y]; rot[Y][Y] = rot[X][X];
				break;
			default:
				assert(0);
				}
			rpuRotate(rp,rot,body);
			break;
			}
		case spin:
		    {
			VECTOR old_spin,d;
			double w_max = 2*PI/sqrt(3*PI/rp->density);
			BOOLEAN incr_only,ang_mom,body;
			COPY_VEC(rp->spin,old_spin); /* needed for spin increment */
			/*DEBUG following only removes net spin---it does not
			  ensure that every particle has wxr_i = 0 and w_i = 0*/
			SCALE_VEC(rp->spin,-1.0);
			rpuAddSpin(rp,FALSE); /* remove current spin */
			(void) printf("Classical breakup limit for measured density = ");
			if (sim_units) (void) printf("%g x 2pi rad/yr",w_max);
			else (void) printf("%g 1/h (P_min = %g h)",
							   3600*w_max/(TWO_PI*T_SCALE),
							   TWO_PI*T_SCALE/(3600*w_max));
			(void) printf("\n");
			incr_only = get_yn("Specify increment only, instead of absolute value","n");
			if (incr_only)
				ang_mom = get_yn("Specify angular momentum increment","n");
			else
				ang_mom = get_yn("Specify angular momentum","n");
			if (ang_mom) {
				if (incr_only)
					(void) printf("Enter angular momentum increment [dlx dly dlz in ");
				else
					(void) printf("Enter new angular momentum [lx ly lz in ");
				if (sim_units) (void) printf("sys units");
				else (void) printf("N m/s");
				(void) printf("]: ");
				(void) scanf("%lf%lf%lf",&d[X],&d[Y],&d[Z]);				
				(void) getchar();
				if (!sim_units)
					SCALE_VEC(d,T_SCALE/(SQ(L_SCALE)*M_SCALE));
				if (incr_only) {
					ADD_VEC(rp->ang_mom,d,rp->ang_mom);
					}
				else {
					COPY_VEC(d,rp->ang_mom);
					}
				rpuAddAngMom(rp);
				if (MAG(rp->spin) > w_max)
					(void) printf("WARNING: exceeds classical breakup limit\n");
				break;
				}
			body = get_yn("Use body axes","n");
			if (incr_only)
				(void) printf("Enter spin increment [dwx dwy dwz in ");
			else
				(void) printf("Enter new spin [wx wy wz in ");
			if (sim_units) (void) printf("units of 2pi rad/yr");
			else (void) printf("1/h");
			(void) printf("]: ");
			(void) scanf("%lf%lf%lf",&d[X],&d[Y],&d[Z]);
			(void) getchar();
			if (!sim_units) SCALE_VEC(d,TWO_PI*T_SCALE/3600);
			if (incr_only) {
				ADD_VEC(old_spin,d,rp->spin);
				}
			else {
				COPY_VEC(d,rp->spin);
				}
			if (MAG(rp->spin) > w_max)
				(void) printf("WARNING: exceeds classical breakup limit\n");
			rpuAddSpin(rp,body);
			break;
			}
		case color:
			{
			BOOLEAN invalid_color;
			int c;
			(void) printf("Color scheme:\n");
			for (i=BLACK;i<FIRST_GRAY;i++)
				(void) printf("%2i. %s\n",i,color_str(i));
			do {
				invalid_color = FALSE;
				(void) printf("Enter new color: ");
				(void) scanf("%i",&c);
				(void) getchar();
				if (c >= NUM_COLORS) { /* allow negative colors */
					(void) printf("Invalid color\n");
					invalid_color = TRUE;
					continue;
					}
				if (c == BLACK || c == FIRST_GRAY)
					(void) printf("WARNING: Particles may be invisible!\n");
				} while (invalid_color);
			rp->color = c;
			rpuApplyColor(rp);
			break;
			}
		case agg_id:
			{
			do {
				(void) printf("Enter new aggregate ID (or -1 to reset): ");
				(void) scanf("%i",&rp->agg_id);
				(void) getchar();
				if (rp->agg_id < -1)
					(void) printf("Invalid ID\n");
				else {
					rpuApplyAggID(rp);
					if (*reindex) {
						printf("\n*** Particle reindexing overridden ***\n\n");
						*reindex = FALSE;
						}
					}
				} while(rp->agg_id < -1);
			break;
			}
		case par_id:
			{
			SSDATA *p,*pmax;
			VECTOR r;
			double lng,lat,dmax,d;
			int c;
			(void) printf("Enter angular coordinates [lng lat in deg]: ");
			(void) scanf("%lf%lf",&lng,&lat);
			(void) getchar();
			lng *= DEG_TO_RAD;
			lat *= DEG_TO_RAD;
			SET_VEC(r,cos(lng)*cos(lat),sin(lng)*cos(lat),sin(lat));
			dmax = 0.0;
			pmax = NULL;
			for (i=0;i<rp->n_particles;i++) {
				p = &rp->data[i];
				/* projected distance along vector minus distance perpendicular
				   to vector favors particles closest to vector */
				d = DOT(p->pos,r) - sqrt(MAG_SQ(p->pos) - SQ(DOT(p->pos,r)));
				if (d > dmax) {
					dmax = d;
					pmax = p;
				}
			}
			if (!pmax) {
				(void) printf("No particle found.\n");
				continue;
			}
			(void) printf("Found particle (original index %i, color %i [%s])\n",
						  pmax->org_idx,pmax->color,color_str(pmax->color));
			while (/*CONSTCOND*/1) {
				(void) printf("Enter new color: ");
				(void) scanf("%i",&c);
				(void) getchar();
				if (c >= NUM_COLORS) { /* allow negative colors */
					(void) printf("Invalid color\n");
					goto invalid;
					}
				if (c == BLACK || c == RED || c == YELLOW || c == MAGENTA ||
					c == CYAN || c == KHAKI || c == FIRST_GRAY) {
					(void) printf("Color %i is reserved\n",c);
					goto invalid;
					}
				break;
			invalid:
				(void) printf("Color scheme:\n");
				for (i=BLACK;i<FIRST_GRAY;i++)
					(void) printf("%2i. %s\n",i,color_str(i));
				}; /* while */
			pmax->color = c;
			break;
			}
		default:
			assert(0);
			}
		}
	}

static void
usage(const char *progname)
{
	printf("Usage: %s [-f] [-h] [-q] [-r] [-s -1|+1] [-u]\n"
		   "where -f forces overwriting of %s if it exists,\n"
		   "      -h prints this help message,\n"
		   "      -q turns on option queries,\n"
		   "      -r renumbers original indices,\n"
		   "      -s sorts by radius: descending (-1) or ascending (+1),\n"
		   "      -u uses simulation units.\n"
		   "By default all options are off and sorting is set to none.\n"
		   "Note -q overrides other options.\n",
		   progname,LOG_FILE);
	exit(1);
	}

int
main(int argc,char *argv[])
{
	int c;

	FILE *fp;
	BOOLEAN force=FALSE,query=FALSE,reindex=FALSE,sim_units=FALSE;
	int sort_order=None;
	RUBBLE_PILE rp[MAX_NUM_FILES],*p;
	char infile[MAXPATHLEN],last_infile[MAXPATHLEN],outfile[MAXPATHLEN];
	double time = 0.0,old_time = 0.0;
	int n_files;

	setbuf(stdout,(char *)NULL);

	srand((unsigned int) getpid());

	/* parse command-line arguments */

	while ((c = getopt(argc,argv,"fhqrs:u")) != EOF) {
		switch (c) {
		case 'f':
			force = TRUE;
			break;
		case 'h':
			usage(argv[0]);
			break; /* to prevent fall-through warning */
		case 'q':
			query = TRUE;
			if (force)
				printf("-q option overrides -f option.\n");
			if (reindex)
				printf("-q option overrides -r option.\n");
			if (sim_units)
				printf("-q option overrides -u option.\n");
			if (sort_order != None)
				printf("-q option overrides -s option.\n");
			force = reindex = sim_units = FALSE;
			sort_order = None;
			break;
		case 'r':
			reindex = TRUE;
			printf("Particle reindexing requested.\n");
			break;
		case 's':
			sort_order = atoi(optarg);
			if (sort_order != Descending && sort_order != Ascending) {
				fprintf(stderr,"Invalid sort order.\n");
				usage(argv[0]);
				}
			printf("Particle size sorting requested (%s).\n",
			       sort_order == Descending ? "descending" :
			       "ascending");
			break;
		case 'u':
			sim_units = TRUE;
			printf("System units selected.\n");
			break;
		case '?':
		default:
			usage(argv[0]);
			}
		}

	if (optind != argc) {
	    fprintf(stderr,"Invalid argument.\n");
		usage(argv[0]);
		}

	/*
	** For backwards compatibility, overwriting and units questions
	** are asked even if query (-q) is not specified unless explicitly
	** overridden with -f and/or -u.  Other questions are not asked
	** unless query is active.
	*/
	
	if (query || !force) {
		fp = fopen(LOG_FILE,"r");
		if (fp != NULL) {
			fclose(fp);
			if (!get_yn("Overwrite log file","y"))
				return 0;
			}
		}

	fp = fopen(LOG_FILE,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open %s for writing.\n",LOG_FILE);
		return 1;
		}

	if (query && !reindex)
		reindex = get_yn("Renumber particle original indices","y");

	if (query && sort_order == None) {
		sort_order = get_yn("Sort in ascending particle size order","y");
		if (sort_order)
			sort_order = Ascending;
		else {
			sort_order = get_yn("Sort in descending particle size order","y");
			if (sort_order)
				sort_order = Descending;
			}
		}

	if (query || !sim_units)
		sim_units = get_yn("Use simulation units (AU, M_sun, etc.)","n");

	n_files = 0;

	while (n_files < MAX_NUM_FILES) {
		infile[0] = '\0';
		printf("File %i [or RETURN to quit]: ",n_files + 1);
		fgets(infile,MAXPATHLEN,stdin);
		assert(strlen(infile) > 0);
		infile[strlen(infile) - 1] = '\0'; /* get rid of newline at end */
		if (strlen(infile) == 0)
			break;
		p = &rp[n_files];
		if (process(infile,p,sim_units,&time,&reindex))
			continue;
		if (n_files == 0)
			old_time = time;
		else if (time != old_time)
			time = 0.0; /* unless all times are the same, set to zero */
		fprintf(fp,
				"File number\t\t%i\n"
				"Filename\t\t%s\n"
				"Mass\t\t\t%e\n"
				"Bulk radius\t\t%e\n"
				"Bulk density\t\t%e\n"
				"Position\t\t%+e\t%+e\t%+e\n"
				"Velocity\t\t%+e\t%+e\t%+e\n"
				"Spin\t\t\t%+e\t%+e\t%+e\n"
				"Major axis\t\t%+f\t%+f\t%+f\n"
				"Inter axis\t\t%+f\t%+f\t%+f\n"
				"Minor axis\t\t%+f\t%+f\t%+f\n"
				"Color\t\t\t%i\n"
				"Aggregate ID\t\t%i\n\n",
				n_files,infile,p->mass,p->radius,p->density,
				p->pos[X],p->pos[Y],p->pos[Z],
				p->vel[X],p->vel[Y],p->vel[Z],
				p->spin[X],p->spin[Y],p->spin[Z],
				p->axes[rp->axis_ord[X]][X],
				p->axes[rp->axis_ord[X]][Y],
				p->axes[rp->axis_ord[X]][Z],
				p->axes[rp->axis_ord[Y]][X],
				p->axes[rp->axis_ord[Y]][Y],
				p->axes[rp->axis_ord[Y]][Z],
				p->axes[rp->axis_ord[Z]][X],
				p->axes[rp->axis_ord[Z]][Y],
				p->axes[rp->axis_ord[Z]][Z],
				p->color,p->agg_id);
		strcpy(last_infile,infile);
		if (++n_files == MAX_NUM_FILES)
			printf("File limit reached.\n");
		}

	if (n_files == 0) {
		fprintf(stderr,"No files specified!\n");
		return 1;
		}

	if (get_yn("Recenter barycentric position and velocity","y")) {
		VECTOR pos,vel,r,v;
		double total_mass;
		int i;

		total_mass = 0;
		ZERO_VEC(pos);
		ZERO_VEC(vel);
		for (i=0;i<n_files;i++) {
			COPY_VEC(rp[i].pos,r);
			SCALE_VEC(rp[i].pos,-1);
			rpuApplyPos(&rp[i]);
			SCALE_VEC(r,rp[i].mass);
			ADD_VEC(pos,r,pos);
			COPY_VEC(rp[i].vel,v);
			SCALE_VEC(rp[i].vel,-1);
			rpuApplyVel(&rp[i]);
			SCALE_VEC(v,rp[i].mass);
			ADD_VEC(vel,v,vel);
			total_mass += rp[i].mass;
			}
		NORM_VEC(pos,total_mass);
		NORM_VEC(vel,total_mass);
		for (i=0;i<n_files;i++) {
			SCALE_VEC(rp[i].pos,-1);
			SUB_VEC(rp[i].pos,pos,rp[i].pos);
			rpuApplyPos(&rp[i]);
			SCALE_VEC(rp[i].vel,-1);
			SUB_VEC(rp[i].vel,vel,rp[i].vel);
			rpuApplyVel(&rp[i]);
			}
		fprintf(fp,"BARYCENTRIC CORRECTION APPLIED\n");
		}

	fclose(fp);

	if (n_files == 2)
		show_encounter(rp,n_files,sim_units);

	do {
		printf("Output file [default %s]: ",last_infile);
		fgets(outfile,MAXPATHLEN,stdin);
		assert(strlen(outfile) > 0);
		outfile[strlen(outfile) - 1] = '\0';
		if (strlen(outfile) == 0)
			strcpy(outfile,last_infile);
		outfile[MAXPATHLEN - 1] = '\0';
		if ((fp = fopen(outfile,"r"))) {
			fclose(fp);
			if (!get_yn("Output file already exists...overwrite","n"))
				continue;
			}
		break;
		} while (/*CONSTCOND*/1);

	if (write_data(outfile,rp,n_files,time,reindex,sort_order) != 0)
		return 1;

	for (--n_files;n_files>=0;n_files--)
		rpuFree(&rp[n_files]);

	printf("Done!\n");

	return 0;
	}

/* rpx.c */
