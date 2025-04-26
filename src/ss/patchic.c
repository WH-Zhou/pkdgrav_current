/*
 ** patchic.c -- DCR 99-01-11
 ** =========
 ** Generates initial conditions for orbiting patch model, in ss format.
 */

#include <ss.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <unistd.h> /* for getopt(), on some systems */
#include <boolean.h>
#include <rdpar.h>
#include <vector.h>
#include <random.h>

#ifndef INT_MAX
#define INT_MAX 2147483647
#endif

#define PAR_FILE "patchic.par"
#define LOG_FILE "patchic.log"

#define GEN_INIT (-1)

enum {LoVerb,MidVerb,HiVerb};          /* for verbosity level */
enum {ThickZ,DispZ,WgtDispZ,Rayleigh}; /* for patch height option */
enum {MaxV,DispV,WgtDispV};            /* for velocity limits option */

typedef struct {

	/* Following read in from parameter file... */

	int iVerbosity,iHgtOpt,iVelOpt;
	double dCentralMass,dOrbDist,dTime,dTau,dSurfDen,dDensity,dRmin,dRmax;
	double dSDExp,dMScaling,dRScaling,dStartTime;
	BOOLEAN bSmooth;
	VECTOR vPatchDim,vVelLim;
	char achOutputFile[MAXPATHLEN];

	/* Following derived from supplied parameters... */

	double dRavg,dR2avg,dR3avg,dMavg,dRHavg,dVesc,dOmega,dLambdaCrit,dt,dVFF;
	int nData;

	/* Run-time flags */

	BOOLEAN bRejects,bAdjCom;

	} PARAMS;

static void
output_data(const PARAMS *p,SSDATA *d)
{
	SSIO ssio;
	SSHEAD h;
	int i;

	if (ssioOpen(p->achOutputFile,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",p->achOutputFile);
		return;
		}

	h.time = p->dStartTime;
	h.n_data = p->nData;
	h.iMagicNumber = SSIO_MAGIC_STANDARD;
	(void) ssioHead(&ssio,&h);

	for (i=0;i<p->nData;i++)
		(void) ssioData(&ssio,&d[i]);

	(void) ssioClose(&ssio);

	if (p->iVerbosity > LoVerb) {
		FILE *fp;
		long nb;

		fp = fopen(p->achOutputFile,"r");
		(void) fseek(fp,0L,SEEK_END);
		nb = ftell(fp);
		(void) printf("I/O: nb=%li ==> nData=%i compare with %i ",nb,
					  (int) (((size_t) nb-sizeof(SSHEAD))/sizeof(SSDATA)),p->nData);
		(void) printf("(ss_head=%i, ss_data=%i)\n",
					  (int)sizeof(SSHEAD),(int)sizeof(SSDATA));
		(void) fclose(fp);
		}
	}

static void
write_log(const PARAMS *p)
{
	FILE *fp;

	fp = fopen(LOG_FILE,"w");
	assert(fp != NULL);
	(void) fprintf(fp,"Number of particles\t%i\n",p->nData);
	(void) fprintf(fp,"Number of steps\t\t%i\n",
				   (int)(TWO_PI*p->dTime/(p->dOmega*p->dt) + 0.5));
	(void) fprintf(fp,"Timestep\t\t%.16e\t# %g s\n",p->dt,p->dt*T_SCALE);
	(void) fprintf(fp,"Central mass\t\t%.16e\t# %g kg\n",p->dCentralMass,
				   p->dCentralMass*M_SCALE);
	(void) fprintf(fp,"Orbital distance\t%.16e\t# %g km\n",p->dOrbDist,p->dOrbDist*L_SCALE/1000);
	(void) fprintf(fp,"Orbital frequency\t%.16e\t# period = %g h\n",p->dOmega,
				   TWO_PI*T_SCALE/(3600*p->dOmega));
	(void) fprintf(fp,"Patch width\t\t%.16e\t# %g m\n",p->vPatchDim[X],
				   p->vPatchDim[X]*L_SCALE);
	(void) fprintf(fp,"Patch length\t\t%.16e\t# %g m\n",p->vPatchDim[Y],
				   p->vPatchDim[Y]*L_SCALE);
	if (p->iHgtOpt == ThickZ)
		(void) fprintf(fp,"Patch thickness\t\t%.16e\t# %g m\n",p->vPatchDim[Z],
					   p->vPatchDim[Z]*L_SCALE);
	else if (p->iHgtOpt == Rayleigh)
		(void) fprintf(fp,"Mean oscillation amp.\t%.16e\t# %g m\n",
					   p->vPatchDim[Z],p->vPatchDim[Z]*L_SCALE);
	else
		(void) fprintf(fp,"Scale height\t\t%.16e\t# %g m\n",p->vPatchDim[Z],
					   p->vPatchDim[Z]*L_SCALE);
	(void) fprintf(fp,"Mean particle radius\t%.16e\t# %g m\n",p->dRavg,
				   p->dRavg*L_SCALE);
	(void) fprintf(fp,"Mean particle mass\t%.16e\t# %g kg\n",p->dMavg,
				   p->dMavg*M_SCALE);
	(void) fprintf(fp,"Mean Hill radius\t%.16e\t# %g <R>\n",p->dRHavg,
				   p->dRHavg/p->dRavg);
	(void) fprintf(fp,"Mean escape speed\t%.16e\t# %g Omega <R>\n",p->dVesc,
				   p->dVesc/(p->dOmega*p->dRavg));
	(void) fprintf(fp,"Optical depth\t\t%.16e\t# %g\n",p->dTau,p->dTau);
	(void) fprintf(fp,"Surface density\t\t%.16e\t# %g kg/m^2\n",p->dSurfDen,
				   p->dSurfDen*M_SCALE/SQ(L_SCALE));
	(void) fprintf(fp,"Particle density\t%.16e\t# %g g/cc\n",p->dDensity,
				   p->dDensity*D_SCALE/1000);
	(void) fprintf(fp,"Critical wavelength\t%.16e\t# %g m\n",p->dLambdaCrit,
				   p->dLambdaCrit*L_SCALE);
	(void) fprintf(fp,"Volume filling factor\t%.16e\t# %g\n",p->dVFF,p->dVFF);
	(void) fprintf(fp,"Mass scaling\t\t%.16e\t# %g\n",p->dMScaling,
				   p->dMScaling);
	(void) fprintf(fp,"Radius scaling\t\t%.16e\t# %g\n",p->dRScaling,
				   p->dRScaling);
	(void) fclose(fp);
	}

static void
gen_particle(const PARAMS *p,SSDATA *d,int i)
{
	/*
	 ** It is important to check whether this particle has been initialized
	 ** before. If so, any initial conditions generated from NON-uniform
	 ** deviates should be preserved in some fashion to prevent biasing.
	 */

	if (d[i].color == GEN_INIT) {
		static double q,a,b;
		static BOOLEAN bFirstCall = TRUE;
		double f;

		if (bFirstCall) {
			q = p->dSDExp + 1;
			if (q == 0)
				a = b = 0;
			else {
				a = pow(p->dRmin,q);
				b = pow(p->dRmax,q);
				}
			bFirstCall = FALSE;
			}
		if (p->bSmooth)
			f = 1 - (double)i/(p->nData - 1); /* largest first */
		else
			f = randUniform();
		if (q == 0)
			d[i].radius = p->dRmin*pow(p->dRmax/p->dRmin,f);
		else
			d[i].radius = pow((1 - f)*a + f*b,1/q);
		d[i].mass = 4.0/3*PI*p->dDensity*CUBE(d[i].radius);
		d[i].color = PLANETESIMAL;
		ZERO_VEC(d[i].spin);
		d[i].mass *= p->dMScaling;
		d[i].radius *= p->dRScaling;
		/*
		 ** Initially pos[X,Y,Z] contains the location of the *guiding center*
		 ** of the particle. Later once the velocity components have been
		 ** determined, this will be turned into the actual particle position.
		 ** This means it's possible for a guiding center to be outside the
		 ** patch boundary, even though the particle is inside. Conversely, the
		 ** particle may end up outside the boundary, in which case we need to
		 ** wrap it around. See the e-mail exchange of 8/21/03 and 8/22/03 with
		 ** Rory Barnes.
		 */
		d[i].pos[X] = (randUniform() - 0.5)*p->vPatchDim[X];
		d[i].pos[Y] = (randUniform() - 0.5)*p->vPatchDim[Y];
		switch (p->iHgtOpt) { /* various options in vertical direction */
		case ThickZ:
			d[i].pos[Z] = (randUniform() - 0.5)*p->vPatchDim[Z];
			break;
		case DispZ:
			d[i].pos[Z] = randGaussian()*p->vPatchDim[Z];
			break;
		case WgtDispZ:
			d[i].pos[Z] = randGaussian()*p->vPatchDim[Z]*sqrt(p->dMavg/d[i].mass);
			break;
		case Rayleigh:
			break; /* case handled below */
		default:
			assert(0);
			}
		switch (p->iVelOpt) { /* various velocity options */
		case MaxV:
			SET_VEC(d[i].vel,
					(2*randUniform() - 1)*p->vVelLim[X],
					(2*randUniform() - 1)*p->vVelLim[Y],
					(2*randUniform() - 1)*p->vVelLim[Z]);
			break;
		case DispV:
			SET_VEC(d[i].vel,
					randGaussian()*p->vVelLim[X],
					randGaussian()*p->vVelLim[Y],
					randGaussian()*p->vVelLim[Z]);
			break;
		case WgtDispV:
			SET_VEC(d[i].vel,
					randGaussian()*p->vVelLim[X],
					randGaussian()*p->vVelLim[Y],
					randGaussian()*p->vVelLim[Z]);
			SCALE_VEC(d[i].vel,sqrt(p->dMavg/d[i].mass));
			break;
		default:
			assert(0);
			}
		/* handle special case of Rayleigh disrtibution */
		if (p->iHgtOpt == Rayleigh) {
			const double Sqrt2OverPi = M_2_SQRTPI/M_SQRT2; /* sqrt(2/PI) */
			double amp = p->vPatchDim[Z]*Sqrt2OverPi*sqrt(-2.0*log(randUniform()))*sqrt(p->dMavg/d[i].mass); /* mass-weighted Rayleigh deviate */
			double phase = TWO_PI*randUniform();
			d[i].pos[Z] = amp*cos(phase);
			d[i].vel[Z] = - amp*p->dOmega*sin(phase);
			}
		/* now determine actual particle position (z values don't change) */
		d[i].pos[X] -= 0.5*d[i].vel[Y]/p->dOmega;
		d[i].pos[Y] += 2*d[i].vel[X]/p->dOmega;
		/* check if particle now outside patch; if so, fix it */
		if (d[i].pos[X] < -0.5*p->vPatchDim[X])
			d[i].pos[X] += p->vPatchDim[X]; /* no y-offset at t=0 */
		else if (d[i].pos[X] > 0.5*p->vPatchDim[X])
			d[i].pos[X] -= p->vPatchDim[X];
		if (d[i].pos[Y] < -0.5*p->vPatchDim[Y])
			d[i].pos[Y] += p->vPatchDim[Y];
		else if (d[i].pos[Y] > 0.5*p->vPatchDim[Y])
			d[i].pos[Y] -= p->vPatchDim[Y];
		assert(d[i].pos[X] >= -0.5*p->vPatchDim[X] && d[i].pos[X] <= 0.5*p->vPatchDim[X]);
		assert(d[i].pos[Y] >= -0.5*p->vPatchDim[Y] && d[i].pos[Y] <= 0.5*p->vPatchDim[Y]);
		d[i].vel[Y] -= 1.5*p->dOmega*d[i].pos[X]; /* shear flow */
		d[i].org_idx = i;
		}
	else {
		d[i].vel[Y] += 1.5*p->dOmega*d[i].pos[X];
		d[i].pos[X] = (randUniform() - 0.5)*p->vPatchDim[X];
		d[i].pos[Y] = (randUniform() - 0.5)*p->vPatchDim[Y];
		d[i].pos[X] -= 0.5*d[i].vel[Y]/p->dOmega;
		d[i].pos[Y] += 2*d[i].vel[X]/p->dOmega;
		if (d[i].pos[X] < -0.5*p->vPatchDim[X])
			d[i].pos[X] += p->vPatchDim[X]; /* no y-offset at t=0 */
		else if (d[i].pos[X] > 0.5*p->vPatchDim[X])
			d[i].pos[X] -= p->vPatchDim[X];
		if (d[i].pos[Y] < -0.5*p->vPatchDim[Y])
			d[i].pos[Y] += p->vPatchDim[Y];
		else if (d[i].pos[Y] > 0.5*p->vPatchDim[Y])
			d[i].pos[Y] -= p->vPatchDim[Y];
		assert(d[i].pos[X] >= -0.5*p->vPatchDim[X] && d[i].pos[X] <= 0.5*p->vPatchDim[X]);
		assert(d[i].pos[Y] >= -0.5*p->vPatchDim[Y] && d[i].pos[Y] <= 0.5*p->vPatchDim[Y]);
		if (p->iHgtOpt == ThickZ) /* other options would introduce bias */
			/*
			** NOTE: if the thickness isn't much bigger than the
			** particle radius, you'll get a slight bimodal
			** distribution in z values if a lot of particles are
			** regenerated.  This is because in the patch model there
			** are no particles above or below the vertical fill
			** limits, so there will always be slightly more room near
			** those limits to place particles. The solution is to
			** choose a larger thickness, or ignore the issue since it
			** will resolve itself quickly once the system starts to
			** evolve.
			*/
			d[i].pos[Z] = (randUniform() - 0.5)*p->vPatchDim[Z];
		d[i].vel[Y] -= 1.5*p->dOmega*d[i].pos[X];
		}
	}

#ifdef BAD_CODE
static int
apply_bcs(const PARAMS *p,SSDATA *d,const VECTOR offset)
{
	static double hx,hy;
	static BOOLEAN bFirstCall = TRUE,bc,redo;

	VECTOR r;
	double xmin,xmax,ymin,ymax,dSafeDist,dx,dy;
	BOOLEAN bNearEdgeX,bNearEdgeY;
	int i,j,ix,iy,n = 0;

	if (bFirstCall) {
		hx = 0.5*p->vPatchDim[X];
		hy = 0.5*p->vPatchDim[Y];
		assert(hx > p->dRmax && hy > p->dRmax);
		bFirstCall = FALSE;
		}

	xmin = -hx; xmax = hx; ymin = -hy; ymax = hy;

	bc = FALSE;

	/*
	 ** Strategy: compensate for depletion of area shifted out of simulation
	 ** region by repopulating equal area at opposite boundary. Unfortunately,
	 ** we must check for overlaps manually with these particles, again to
	 ** eliminate repopulation bias (pkdgrav won't distinguish regions).
	 */

	for (i=0;i<p->nData;i++) {
		if (d[i].pos[X] < xmin) {bc = TRUE; xmin = xmax - offset[X];}
		if (d[i].pos[X] > xmax) {bc = TRUE; xmax = xmin - offset[X];}
		if (d[i].pos[Y] < ymin) {bc = TRUE; ymin = ymax - offset[Y];}
		if (d[i].pos[Y] > ymax) {bc = TRUE; ymax = ymin - offset[Y];}
		if (bc) {
			dSafeDist = d[i].radius + p->dRmax;
			d[i].vel[Y] += 1.5*p->dOmega*d[i].pos[X]; /* remove shear */
			do {
				d[i].pos[X] = xmin + randUniform()*(xmax - xmin);
				d[i].pos[Y] = ymin + randUniform()*(ymax - ymin);
				if (!p->bRejects) break;
				bNearEdgeX = bNearEdgeY = TRUE; /* initialization */
				if ((dx = hx + d[i].pos[X]) > dSafeDist)
					if ((dx = hx - d[i].pos[X]) > dSafeDist)
						bNearEdgeX = FALSE;
				if (bNearEdgeX) dx -= d[i].radius; /* dist. to nearest edge */
				if ((dy = hy + d[i].pos[Y]) > dSafeDist)
					if ((dy = hy - d[i].pos[Y]) > dSafeDist)
						bNearEdgeY = FALSE;
				if (bNearEdgeY) dy -= d[i].radius;
				redo = FALSE;
				for (j=0;j<p->nData;j++) {
					if (j == i) continue;
					if ((d[j].pos[X] >= xmin && d[j].pos[X] <= xmax) ||
						(d[j].pos[Y] >= ymin && d[j].pos[Y] <= ymax)) {
						if (OVERLAP(d[i].pos,d[i].radius,
									d[j].pos,d[j].radius)) {
							redo = TRUE;
							break;
							}
						}
					if ((bNearEdgeX || bNearEdgeY) &&
						(hx + d[j].pos[X] <= d[j].radius - dx ||
						 hx - d[j].pos[X] <= d[j].radius - dx ||
						 hy + d[j].pos[Y] <= d[j].radius - dy ||
						 hy - d[j].pos[Y] <= d[j].radius - dy)) {
						r[Z] = d[j].pos[Z]; /* no boundary in z */
						for (ix=-1;ix<=1;ix++) {
							r[X] = d[j].pos[X];
							if (ix) r[X] += ix*p->vPatchDim[X];
							for (iy=-1;iy<=1;iy++) {
								if (ix == 0 && iy == 0) continue;
								r[Y] = d[j].pos[Y];
								if (iy) r[Y] += iy*p->vPatchDim[Y];
								if (OVERLAP(d[i].pos,d[i].radius,r,d[j].radius)) {
									redo = TRUE;
									break;
									}
								}
							if (redo) break;
							}
						}
					if (redo) break;
					}
				} while(redo);
			d[i].vel[Y] -= 1.5*p->dOmega*d[i].pos[X]; /* put shear back */
			xmin = -hx; xmax = hx; ymin = -hy; ymax = hy;
			bc = FALSE;
			++n;
			}
		}

	return n;
	}
#endif

static void
get_com_pos(const PARAMS *p,const SSDATA *d,VECTOR vComPos)
{
	VECTOR v;
	double dTotalMass;
	int i;

	dTotalMass = 0;
	ZERO_VEC(vComPos);
	for (i=0;i<p->nData;i++) {
		dTotalMass += d[i].mass;
		COPY_VEC(d[i].pos,v);
		SCALE_VEC(v,d[i].mass);
		ADD_VEC(vComPos,v,vComPos);
		}
	assert(dTotalMass > 0.0);
	NORM_VEC(vComPos,dTotalMass);
	}

static void
sub_pos(const PARAMS *p,SSDATA *d,const VECTOR vPos)
{
	int i;

	for (i=0;i<p->nData;i++) {
		SUB_VEC(d[i].pos,vPos,d[i].pos);
		d[i].vel[Y] += 1.5*p->dOmega*vPos[X]; /* shear correction */
		}
	}

static void
get_com_vel(const PARAMS *p,const SSDATA *d,VECTOR vComVel)
{
	VECTOR v;
	double dTotalMass;
	int i;

	dTotalMass = 0;
	ZERO_VEC(vComVel);
	for (i=0;i<p->nData;i++) {
		dTotalMass += d[i].mass;
		COPY_VEC(d[i].vel,v);
		v[Y] += 1.5*p->dOmega*d[i].pos[X]; /* shear correction */
		SCALE_VEC(v,d[i].mass);
		ADD_VEC(vComVel,v,vComVel);
		}
	assert(dTotalMass > 0.0);
	NORM_VEC(vComVel,dTotalMass);
	}

static void
sub_vel(const PARAMS *p,SSDATA *d,const VECTOR vVel)
{
	int i;

	for (i=0;i<p->nData;i++)
		SUB_VEC(d[i].vel,vVel,d[i].vel);
	}

static void
adj_com(const PARAMS *p,SSDATA *d)
{
	VECTOR vComPos,vComVel;

	if (p->iVerbosity > LoVerb)
		(void) printf("Adjusting centre-of-mass position...\n");
	get_com_pos(p,d,vComPos);
	if (p->iVerbosity > LoVerb)
		(void) printf("dx = %g Lx = %g <R>\ndy = %g Ly = %g <R>\n"
					  "dz = %g Lz = %g <R>\n",
					  vComPos[X]/p->vPatchDim[X],
					  vComPos[X]/p->dRavg,
					  vComPos[Y]/p->vPatchDim[Y],
					  vComPos[Y]/p->dRavg,
					  vComPos[Z]/p->vPatchDim[Z],
					  vComPos[Z]/p->dRavg);
	sub_pos(p,d,vComPos);
	if (p->iVerbosity > LoVerb)
		(void) printf("Adjusting centre-of-mass velocity...\n");
	get_com_vel(p,d,vComVel);
	if (p->iVerbosity > LoVerb)
		(void) printf("|dv| = %g Omega <R>\n",
					  MAG(vComVel)/(p->dOmega*p->dRavg));
	sub_vel(p,d,vComVel);
/*
   if (p->iVerbosity == HiVerb)
   (void) printf("Checking boundary conditions...\n");
   n = apply_bcs(p,d,vComPos);
   if (p->iVerbosity == HiVerb)
   (void) printf("%i correction%s applied.\n",n,n==1?"":"s");
   } while (n);
*/
	}

static void
generate(const PARAMS *p,SSDATA *d)
{
	int i;

	if (p->iVerbosity == HiVerb)
		(void) printf("Generating %i particles...\n",p->nData);
	for (i=0;i<p->nData;i++)
		gen_particle(p,d,i);
	if (p->bAdjCom) adj_com(p,d);
	}

static int
fix_rejects(const PARAMS *p,SSDATA *d)
{
	SSIO ssio;
	SSHEAD h;
	FILE *fp;
	int nRej=0,iRej,i;

	if (ssioOpen(p->achOutputFile,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",p->achOutputFile);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"Corrupt header.\n");
		(void) ssioClose(&ssio);
		return 1;
		}

	if (h.n_data != p->nData) {
		(void) fprintf(stderr,"Incompatible data file.\n");
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

	for (i=0;i<p->nData;i++)
		if (ssioData(&ssio,&d[i]) ||
			(d[i].color != PLANETESIMAL && d[i].color != TEST)) {
			(void) fprintf(stderr,"Corrupt or invalid data.\n");
			(void) ssioClose(&ssio);
			return 1;
			}

	(void) ssioClose(&ssio);

	if (!(fp = fopen(REJECTS_FILE,"r"))) {
		(void) fprintf(stderr,"WARNING: Unable to open \"%s\".\n"
					   "Assuming no rejects found.\n",REJECTS_FILE);
		goto finish;
		}

	(void) fscanf(fp,"%i",&i);
	if (i != p->nData) {
		(void) fprintf(stderr,"Incompatible rejects file.\n");
		return 1;
		}

	i = 0;

	while (fscanf(fp,"%i",&iRej) != EOF) {
		if (iRej >= 0) {
			if (i != iRej) {
				(void) fprintf(stderr,"Corrupted reject data.\n");
				(void) fclose(fp);
				return 1;
				}
			gen_particle(p,d,iRej);
			if (p->iVerbosity == HiVerb)
				(void) printf("Particle %i regenerated.\n",iRej);
			++nRej;
			}
		++i;
		}

	(void) fclose(fp);

 finish:

	if (p->iVerbosity > LoVerb)
		(void) printf("%i particle%s regenerated.\n",nRej,nRej==1?"":"s");

	if (p->bAdjCom) adj_com(p,d);

	return 0;
	}

static double
sd_mom(double dRmin,double dRmax,double dAlpha,int n)
{
	/* for calculating moments of size distribution */

	double x;

	assert(dRmax > dRmin);
	assert(n > 0);
	if (dAlpha == -1) {
		x = pow(dRmax/dRmin,n);
		return pow(dRmin,n)*(x - 1)/log(x);
		}
	else {
		double dBeta = dAlpha + 1;
		x = pow(dRmax,dBeta) - pow(dRmin,dBeta);
		if (dBeta + n == 0)
			return dBeta*log(dRmax/dRmin)/x;
		else
			return (dBeta/(dBeta + n))*
				(pow(dRmax,dBeta + n) - pow(dRmin,dBeta + n))/x;
		}
	}

static void
get_params(PARAMS *p)
{
	/* Read in supplied parameters */

	int k;

	OpenPar(PAR_FILE);
	ReadInt("Verbosity level",&k);
	assert (k == LoVerb || k == MidVerb || k == HiVerb);
	p->iVerbosity = k;
	ReadDbl("Central mass",&p->dCentralMass);
	assert(p->dCentralMass > 0.0);
	p->dCentralMass /= M_SCALE;
	ReadDbl("Orbital distance",&p->dOrbDist);
	assert(p->dOrbDist > 0.0);
	p->dOrbDist *= 1000/L_SCALE;
	ReadDbl("Number of orbits",&p->dTime);
	assert(p->dTime >= 0.0);
	ReadDbl("Dynamical optical depth",&p->dTau);
	assert(p->dTau >= 0.0);
	ReadDbl("Surface density",&p->dSurfDen);
	assert(p->dSurfDen >= 0.0);
	p->dSurfDen *= SQ(L_SCALE)/M_SCALE;
	ReadDbl("Particle density",&p->dDensity);
	assert(p->dDensity >= 0.0);
	p->dDensity *= 1000/D_SCALE;
	assert((p->dTau == 0.0 && p->dSurfDen > 0.0 && p->dDensity > 0.0) ||
		   (p->dSurfDen == 0.0 && p->dTau > 0.0 && p->dDensity > 0.0) ||
		   (p->dDensity == 0.0 && p->dTau > 0.0 && p->dSurfDen > 0.0));
	ReadDbl("Minimum particle radius",&p->dRmin);
	assert(p->dRmin > 0.0);
	ReadDbl("Maximum particle radius",&p->dRmax);
	assert(p->dRmax > 0.0);
	assert(p->dRmax >= p->dRmin);
	p->dRmin /= L_SCALE;
	p->dRmax /= L_SCALE;
	ReadDbl("Size distribution exponent",&p->dSDExp);
	if (p->iVerbosity > LoVerb && p->dSDExp > 0)
		(void) fprintf(stderr,"WARNING: Size distribution exponent > 0...\n"
					   "...adjust & reject routines may perform poorly.\n");
	ReadInt("Use smooth distribution?",&k);
	assert(k == 0 || k == 1);
	p->bSmooth = k;
	ReadNDbl("Patch dimensions",p->vPatchDim,N_DIM);
	assert(p->vPatchDim[X] != 0.0 && p->vPatchDim[Y] != 0.0);
	ReadInt("Patch height option",&k);
	assert(k == ThickZ || k == DispZ || k == WgtDispZ || k == Rayleigh);
	p->iHgtOpt = k;
	ReadNDbl("Velocity limits",p->vVelLim,N_DIM);
	ReadInt("Velocity limits option",&k);
	assert(k == MaxV || k == DispV || k == WgtDispV);
	p->iVelOpt = k;
	ReadDbl("Mass scaling",&p->dMScaling);
	assert(p->dMScaling > 0.0);
	ReadDbl("Radius scaling",&p->dRScaling);
	assert(p->dRScaling > 0.0);
	ReadDbl("Start time",&p->dStartTime);
	ReadStr("Output file name",p->achOutputFile,MAXPATHLEN);
	ClosePar();

	/* Calculate remaining parameters */

	if (p->dRmax > p->dRmin) {
		p->dRavg  = sd_mom(p->dRmin,p->dRmax,p->dSDExp,1);
		p->dR2avg = sd_mom(p->dRmin,p->dRmax,p->dSDExp,2);
		p->dR3avg = sd_mom(p->dRmin,p->dRmax,p->dSDExp,3);
		}
	else {
		p->dRavg  = p->dRmax;
		p->dR2avg = SQ(p->dRavg);
		p->dR3avg = CUBE(p->dRavg);
		}
	if (p->iVerbosity > LoVerb)
		(void) printf("Mean %sparticle radius = %g m\n",
					  p->dRScaling != 1 ? "(unscaled) " : "",p->dRavg*L_SCALE);
	if (p->dTau == 0) {
		p->dTau = 0.75*p->dSurfDen/p->dDensity*(p->dR2avg/p->dR3avg);
		if (p->iVerbosity > LoVerb)
			(void) printf("Dynamical optical depth = %g\n",p->dTau);
		}
	else if (p->dSurfDen == 0) {
		p->dSurfDen = 4.0/3*p->dDensity*p->dTau*(p->dR3avg/p->dR2avg);
		if (p->iVerbosity > LoVerb)
			(void) printf("Surface density = %g kg/m^2\n",
						  p->dSurfDen*M_SCALE/SQ(L_SCALE));
		}
	else {
		p->dDensity = 0.75*p->dSurfDen/p->dTau*(p->dR2avg/p->dR3avg);
		if (p->iVerbosity > LoVerb)
			(void) printf("Particle density = %g g/cc\n",
						  p->dDensity*D_SCALE/1000);
		}
	p->dMavg = 4.0/3*PI*p->dDensity*p->dR3avg;
	if (p->iVerbosity > LoVerb)
		(void) printf("Mean %sparticle mass = %g kg\n",
					  p->dMScaling != 1 ? "(unscaled) " : "",p->dMavg*M_SCALE);
	p->dRHavg =
		p->dOrbDist*pow(8*PI*p->dDensity/(9*p->dCentralMass),1.0/3)*p->dRavg;
	if (p->iVerbosity > LoVerb)
		(void) printf("Mean Hill radius = %g <R>\n",p->dRHavg/p->dRavg);
	p->dOmega = sqrt(p->dCentralMass/CUBE(p->dOrbDist));
	assert(p->dOmega > 0.0);
	if (p->iVerbosity > LoVerb)
		(void) printf("Orbital frequency = %g rad/s (orbital period = %g h)\n",
					  p->dOmega/T_SCALE,TWO_PI*T_SCALE/(p->dOmega*3600));
	p->dLambdaCrit = 4*SQ(PI)*p->dSurfDen/SQ(p->dOmega);
	if (p->iVerbosity > LoVerb)
		(void) printf("Critical wavelength = %g m (= %g <R>)\n",
					  p->dLambdaCrit*L_SCALE,p->dLambdaCrit/p->dRavg);
	for (k=0;k<N_DIM;k++)
		if (p->vPatchDim[k] < 0)
			p->vPatchDim[k] /= -L_SCALE;
		else if (k < Z)
			p->vPatchDim[k] *= p->dLambdaCrit;
		else
			p->vPatchDim[k] *= p->dRavg;
	if (p->vPatchDim[X] > p->vPatchDim[Y])
		(void) fprintf(stderr,"WARNING: Lx > Ly\n");
	p->nData = (int) (p->dTau*p->vPatchDim[X]*p->vPatchDim[Y]/(PI*p->dR2avg));
	if (p->iVerbosity > LoVerb)
		(void) printf("N = %i\n",p->nData);
	assert(p->nData > 1 && p->nData < INT_MAX);
	p->dt = 0.01/sqrt(p->dDensity); /*DEBUG! was 0.03*/
	if (p->dt > 0.01*TWO_PI/p->dOmega) /* 1/100th of an orbit (conservative) */
	  p->dt = 0.01*TWO_PI/p->dOmega;
	if (p->iVerbosity > LoVerb)
		(void) printf("Recommended timestep = %g (1 orbit = %g steps)\n",
					  p->dt,TWO_PI/(p->dOmega*p->dt));
	if (p->iVerbosity > LoVerb && p->dMScaling != 1)
		(void) printf("Particle masses will be scaled by %g\n",p->dMScaling);
	if (p->iVerbosity > LoVerb && p->dRScaling != 1)
		(void) printf("Particle radii will be scaled by %g\n",p->dRScaling);
	p->dRavg *= p->dRScaling;
	p->dR2avg *= SQ(p->dRScaling);
	p->dR3avg *= CUBE(p->dRScaling);
	p->dMavg *= p->dMScaling;
	p->dVesc = p->dRavg*sqrt(8.0/3*PI*p->dDensity*
							 p->dMScaling/CUBE(p->dRScaling));
	if (p->iVerbosity > LoVerb)
		(void) printf("%s escape speed = %g m/s = %g Omega <R>\n",
					  p->dMScaling != 1 || p->dRScaling != 1 ? "Scaled mean" :
					  "Mean",p->dVesc*V_SCALE,p->dVesc/(p->dOmega*p->dRavg));
	for (k=0;k<N_DIM;k++)
		if (p->vVelLim[k] < 0)
			p->vVelLim[k] /= -V_SCALE;
		else
			p->vVelLim[k] *= p->dOmega*p->dRavg;
	p->dVFF = (p->iHgtOpt == ThickZ ? 1 : 0.68)*p->nData*4.0/3*PI*p->dR3avg/
		(p->vPatchDim[X]*p->vPatchDim[Y]*p->vPatchDim[Z]*
		 (p->iHgtOpt == ThickZ ? 1 : 2));
	if (p->iVerbosity > LoVerb)
		(void) printf("Predicted %svolume filling factor = %g\n",
					  p->dRScaling != 1 ? "scaled " : "",p->dVFF);
	}

int
main(int argc,char *argv[])
{
	PARAMS params;
	SSDATA *data = NULL;
	BOOLEAN bAdjCom=FALSE,bForce=FALSE,bRejects=FALSE,bUsage=FALSE;
	int c;

	/* Disable stdout buffering */

	setbuf(stdout,(char *)NULL);

	/* Check arguments */

	while ((c = getopt(argc,argv,"afr")) != EOF)
		switch (c) {
		case 'a':
			bAdjCom = TRUE;
			break;
		case 'f':
			bForce = TRUE;
			break;
		case 'r':
			bRejects = TRUE;
			break;
		default:
			bUsage = TRUE;
			}

	if (bUsage || optind < argc) {
		(void) fprintf(stderr,"Usage: %s [ -a ] [ -f ] [ -r ]\n",argv[0]);
		return 1;
		}

	/* Get model parameters (needed now for output file name) */

	get_params(&params);
	params.bAdjCom = bAdjCom;
	params.bRejects = bRejects;

	if (!bForce && !bRejects) {
		FILE *fp = fopen(params.achOutputFile,"r");
		if (fp) {
			(void) fprintf(stderr,"%s exists -- use \"-f\" to overwrite.\n",
						   params.achOutputFile);
			(void) fclose(fp);
			return 1;
			}
		}

	/* Generate initial conditions */

    randSeedGenerator(randReadUrandom());

	if (!(data = (SSDATA *) malloc((size_t) params.nData*sizeof(SSDATA)))) {
		(void) fprintf(stderr,"Unable to allocate data memory.\n");
		return 1;
		}

	if (bRejects) {
		if (fix_rejects(&params,data)) {
			free((void *)data);
			return 1;
			}
		}
	else {
		int i;
		for (i=0;i<params.nData;i++)
			data[i].color = GEN_INIT; /* initialize */
		generate(&params,data);
		write_log(&params);
		}

	if (params.iVerbosity > LoVerb) {
		VECTOR vComPos,vComVel,vVelDisp,v;
		double h,lz,dTotalMass,dArea,dVol,xy,xyz,ls,vs;
		int i,k;
		h = params.vPatchDim[Z];
		if (params.iHgtOpt == ThickZ) h *= 0.5; /* half thickness */
		lz = 2*h;
		dTotalMass = dArea = dVol = 0;
		for (i=0;i<params.nData;i++) {
			dTotalMass += data[i].mass;
			dArea += SQ(data[i].radius);
			if (fabs(data[i].pos[Z]) < h) dVol += CUBE(data[i].radius);
			}
		get_com_pos(&params,data,vComPos);
		get_com_vel(&params,data,vComVel);
		ZERO_VEC(vVelDisp);
		for (i=0;i<params.nData;i++) {
			COPY_VEC(data[i].vel,v);
			v[Y] += 1.5*params.dOmega*data[i].pos[X];
			SUB_VEC(data[i].vel,vComVel,data[i].vel);
			for (k=0;k<N_DIM;k++)
				vVelDisp[k] += data[i].mass*SQ(v[k]);
			}
		NORM_VEC(vVelDisp,dTotalMass);
		for (k=0;k<N_DIM;k++)
			vVelDisp[k] = sqrt(vVelDisp[k]);
		NORM_VEC(vVelDisp,params.dOmega*params.dRavg);
		xy = params.vPatchDim[X]*params.vPatchDim[Y];
		xyz = xy*lz;
		ls = sqrt(xy);
		vs = params.dOmega*ls;
		NORM_VEC(vComPos,ls);
		NORM_VEC(vComVel,vs);
		(void) printf("Centre-of-mass position = (%g,%g,%g)\n",
					  vComPos[X],vComPos[Y],vComPos[Z]);
		(void) printf("Centre-of-mass velocity = (%g,%g,%g)\n",
					  vComVel[X],vComVel[Y],vComVel[Z]);
		(void) printf("Actual dynamical optical depth = %g\n",PI*dArea/xy);
		(void) printf("Actual surface density = %g kg/m^2\n",
					  dTotalMass/xy*M_SCALE/SQ(L_SCALE));
		(void) printf("Actual volume filling factor = %g\n",
					  4.0/3*PI*dVol/xyz);
		(void) printf("Actual velocity dispersion = (%g,%g,%g) Omega <R>\n",
					  vVelDisp[X],vVelDisp[Y],vVelDisp[Z]);
		(void) printf("                 magnitude = %g Omega <R> = %g m/s\n",
					  MAG(vVelDisp),MAG(vVelDisp)*params.dOmega*params.dRavg*
					  V_SCALE);
		(void) printf("***IMPORTANT NOTE***\n"
				"Spins are taken to be in the rotating frame.\n"
				"This means particles initially with zero spin\n"
				"actually have space z-spins of Omega, such that\n"
				"they always show the same face to the planet.\n");
		}

	/* Save data */

	output_data(&params,data);

	/* All done */

	free((void *)data);

	return 0;
	}

/* patchic.c */
