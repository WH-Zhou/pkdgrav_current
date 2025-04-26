/*
 ** ssgen_patch.c -- Ronald Ballouz 02/06/17
 ** =======
 ** Generates initial conditions for SLIDING_PATCH runs.
 ** Written to replace dopatch script which used pkdgrav to find overlaps.
*/

#include <ss.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>	/* for time() & getpid() */
#include <time.h>		/* for time() */
#include <unistd.h>		/* for getpid(), and getopt() if needed */
#include <boolean.h>
#include <rdpar.h>
#include <vector.h>
#include <random.h>

#define PAR_FILE "patchic.par"
#define LOG_FILE "patchic.log"

enum {LoVerb,MidVerb,HiVerb};          /* for verbosity level */
enum {ThickZ,DispZ,WgtDispZ,Rayleigh}; /* for patch height option */
enum {MaxV,DispV,WgtDispV};            /* for velocity limits option */

double DISTSQ(double v1[3],double v2[3]);
#define DISTSQ(v1,v2) ((SQ(v1[0] - v2[0]) + SQ(v1[1] - v2[1]) + SQ(v1[2] - v2[2])))

#define MAX_NUM_PASS_DFLT 500

typedef struct {                
  int nPassMax;
  /*double dRadAvg,dBRad,dBDensity,dRadDev,dCutoff,dPIndex,dPRmax,dPFlag;*/
  
  /* Following read in from parameter file... */
  int iVerbosity,iHgtOpt,iVelOpt,iParticleColor; 
  double dCentralMass,dOrbDist,dTime,dTau,dSurfDen,dDensity,dRmin,dRmax;
  double dSDExp,dMScaling,dRScaling,dStartTime;
  BOOLEAN bSmooth;
  VECTOR vPatchDim,vVelLim;
  char achOutputFile[MAXPATHLEN];
  
  /*New Features from SSGEN*/
  int bDensDist;
  double dDensityMin,dDensityMinFrac,dDensityAvg;
  
  /* Following derived from supplied parameters... */
  double dRavg,dR2avg,dR3avg,dMavg,dRHavg,dVesc,dOmega,dLambdaCrit,dt,dVFF;
  int nData;
  
  /* SSDEM Parameters */
  double dDelta, dKn;
  
  /* Run-time flags */
  BOOLEAN bReverseSort, bAdjCom;
  
  } PARAMS;

typedef struct {
  int keep;
  double mass;
  double radius;
  double pos[N_DIM];
  double vel[N_DIM];
  double spin[N_DIM];
  int color, org_idx;
  } DATA;

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  double cmin[N_DIM],cmin_eff[N_DIM],cmax[N_DIM],cmax_eff[N_DIM],pos[N_DIM]; /*NDIM=3 defined in ssio.h*/
  struct node_s *cell[CELLS_PER_NODE];
  DATA *leaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

/*** START TREE CODE ROUTINES ***/
static void kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}

/*function that squezzes particles into box by checking for overlaps. Uses node_s structure*/
static int reject(const NODE *node,const DATA *d)
{
	NODE *n;
	double r2;
	int i,k;

	/* intersect test -- check each subcell in turn */

	/*
	** Algorithm adapted from...
	** "A Simple Method for Box-Sphere Intersection Testing",
	** by Jim Arvo, in "Graphics Gems", Academic Press, 1990.
	** (http://www.ics.uci.edu/~arvo/code/BoxSphereIntersect.c)
	*/

	for (i=0;i<CELLS_PER_NODE;i++) {
		n = node->cell[i];
		if (n != NULL) {
			r2 = 0.0;
			for (k=0;k<N_DIM;k++)
				if (d->pos[k] < n->cmin_eff[k])
					r2 += SQ(d->pos[k] - n->cmin_eff[k]);
				else if (d->pos[k] > n->cmax_eff[k])
					r2 += SQ(d->pos[k] - n->cmax_eff[k]);
			if (r2 <= SQ(d->radius) && reject(n,d)) /*reject called through root than loops through all the nodes till it finds a NULL*/
				return 1;
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != d &&
				 node->leaf[i]->keep &&
				 DISTSQ(d->pos,node->leaf[i]->pos) <=
				 SQ(d->radius + node->leaf[i]->radius))
			return 1;
		}

	return 0;
	}

/* makes node given xyz min's and maxes from user input, to be used by rejection algorithm */
static void make_node(double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,NODE **node)
{
	int i;

	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	(*node)->cmin[0] = (*node)->cmin_eff[0] = xmin; /* define boundaries of rectangular cell */
	(*node)->cmin[1] = (*node)->cmin_eff[1] = ymin;
	(*node)->cmin[2] = (*node)->cmin_eff[2] = zmin;
	(*node)->cmax[0] = (*node)->cmax_eff[0] = xmax;
	(*node)->cmax[1] = (*node)->cmax_eff[1] = ymax;
	(*node)->cmax[2] = (*node)->cmax_eff[2] = zmax;

	(*node)->pos[0] = 0.5*((*node)->cmin[0] + (*node)->cmax[0]); /* define center of cell */
	(*node)->pos[1] = 0.5*((*node)->cmin[1] + (*node)->cmax[1]);
	(*node)->pos[2] = 0.5*((*node)->cmin[2] + (*node)->cmax[2]);

	for (i=0;i<CELLS_PER_NODE;i++) { /* initalize the components of the struct, also required to end recursion */
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}
	}

/* part of rejection algorithm */
static void add_to_tree(NODE *node,DATA *d)
{
	int i,idx,idy,idz;

	idx = (d->pos[0] < node->pos[0] ? -1 : 1); /*tri-conditional*/
	idy = (d->pos[1] < node->pos[1] ? -1 : 1);
	idz = (d->pos[2] < node->pos[2] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1)); /* determine which octant the node/particle belongs in */

	if (node->cell[i] != NULL)
		add_to_tree(node->cell[i],d);
	else if (node->leaf[i] != NULL) {
		make_node(idx < 0 ? node->cmin[0] : node->pos[0],
				  idx < 0 ? node->pos[0] : node->cmax[0],
				  idy < 0 ? node->cmin[1] : node->pos[1],
				  idy < 0 ? node->pos[1] : node->cmax[1],
				  idz < 0 ? node->cmin[2] : node->pos[2],
				  idz < 0 ? node->pos[2] : node->cmax[2],
				  &node->cell[i]);
		add_to_tree(node->cell[i],node->leaf[i]);
		add_to_tree(node->cell[i],d);
		node->leaf[i] = NULL;
		}
	else
		node->leaf[i] = d;
	
	if (node->cmin[0] - d->radius < node->cmin_eff[0])
		node->cmin_eff[0] = node->cmin[0] - d->radius;
	if (node->cmax[0] + d->radius > node->cmax_eff[0])			
		node->cmax_eff[0] = node->cmax[0] + d->radius;
	if (node->cmin[1] - d->radius < node->cmin_eff[1])
		node->cmin_eff[1] = node->cmin[1] - d->radius;
	if (node->cmax[1] + d->radius > node->cmax_eff[1])
		node->cmax_eff[1] = node->cmax[1] + d->radius;
	if (node->cmin[2] - d->radius < node->cmin_eff[2])
		node->cmin_eff[2] = node->cmin[2] - d->radius;
	if (node->cmax[2] + d->radius > node->cmax_eff[2])
		node->cmax_eff[2] = node->cmax[2] + d->radius;
	}

/**** END TREE CODE ROUTINES ****/

/*** START SLIDING_PATCH ROUTINES ***/
static void
get_com_pos(const PARAMS *p,const DATA *d,VECTOR vComPos)
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
sub_pos(const PARAMS *p,DATA *d,const VECTOR vPos)
{
  int i;
  
  for (i=0;i<p->nData;i++) {
    SUB_VEC(d[i].pos,vPos,d[i].pos);
    d[i].vel[Y] += 1.5*p->dOmega*vPos[X]; /* shear correction */
		}
}

static void
get_com_vel(const PARAMS *p,const DATA *d,VECTOR vComVel)
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
sub_vel(const PARAMS *p,DATA *d,const VECTOR vVel)
{
  int i;
  
  for (i=0;i<p->nData;i++)
    SUB_VEC(d[i].vel,vVel,d[i].vel);
}

static void
adj_com(const PARAMS *p,DATA *d)
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
/*** END SLIDING PATCH ROUTINES ***/


/**** START RADIUS SORTING ROUTINES ***/
/*making sure that radius matchup in data. (ordering by increasing radius)*/
static int compar(const void *v1,const void *v2)
{
	if (((const DATA *)v1)->radius < ((const DATA *)v2)->radius)
		return -1;
	else if (((const DATA *)v1)->radius > ((const DATA *)v2)->radius)
		return 1;
	else
		return 0;
	}

static int rev_compar(const void *v1,const void *v2)
{
  if (((const DATA *)v1)->radius > ((const DATA *)v2)->radius)
    return -1;
  else if (((const DATA *)v1)->radius < ((const DATA *)v2)->radius)
    return 1;
  else
    return 0;
}
/*** END RADIUS SORTING ROUTINES***/

/*** START PARTICLE GENERATING ROUTINES****/
#ifdef UNUSED
static double gaussian(double dSig,double dLim)
{
	/* truncated Gaussian */

	double dVal;

	do {
		dVal = dSig*randGaussian();
		} while (dLim > 0.0 && fabs(dVal) > dLim*dSig);

	return dVal;
	}
#endif /* UNUSED */

static void get_pos(const PARAMS *p,DATA *d)
{
	/*
	** NOTE: for the ellipsoidal or cylindrical region, this routine
	** guarantees that the particle center -- but not necessarily the
	** entire particle -- is contained within the ellipsoid (both are
	** guaranteed for the rectangular region).
	*/

	int k;
	double a;

  for (k=0;k<N_DIM;k++) {
    a = p->vPatchDim[k] - 2.0*d->radius;
    if (a < 0.0) {
      fprintf(stderr,"get_pos(): Region too small for particle.\n");
      exit(1);
      }
    if (k!=Z) {
      d->pos[k] = (randUniform() - 0.5)*a;
      }
    else if (k==Z) {
      switch (p->iHgtOpt) { /* various options in vertical direction */
        case ThickZ:
          d->pos[Z] = (randUniform() - 0.5)*p->vPatchDim[Z];
          break;
        case DispZ:
          d->pos[Z] = randGaussian()*p->vPatchDim[Z];
          break;
        case WgtDispZ:
          d->pos[Z] = randGaussian()*p->vPatchDim[Z]*sqrt(p->dMavg/d->mass);
          break;
        case Rayleigh:
          break; /* case handled below */
        default:
          assert(0);
        }
      }
    }
  /* handle special case of Rayleigh disrtibution */
  if (p->iHgtOpt == Rayleigh) {
    const double Sqrt2OverPi = M_2_SQRTPI/M_SQRT2; /* sqrt(2/PI) */
    double amp = p->vPatchDim[Z]*Sqrt2OverPi*sqrt(-2.0*log(randUniform()))*sqrt(p->dMavg/d->mass); /* mass-weighted Rayleigh deviate */
    double phase = TWO_PI*randUniform();
    d->pos[Z] = amp*cos(phase);
    d->vel[Z] = - amp*p->dOmega*sin(phase);
    }
  /* now determine actual particle position (z values don't change)
  d->pos[X] -= 0.5*d->vel[Y]/p->dOmega;
  d->pos[Y] += 2*d->vel[X]/p->dOmega;
  if (d->pos[X] < -0.5*p->vPatchDim[X])
    d->pos[X] += p->vPatchDim[X];
  else if (d->pos[X] > 0.5*p->vPatchDim[X])
    d->pos[X] -= p->vPatchDim[X];
  if (d->pos[Y] < -0.5*p->vPatchDim[Y])
    d->pos[Y] += p->vPatchDim[Y];
  else if (d->pos[Y] > 0.5*p->vPatchDim[Y])
    d->pos[Y] -= p->vPatchDim[Y];
  assert(d->pos[X] >= -0.5*p->vPatchDim[X] && d->pos[X] <= 0.5*p->vPatchDim[X]);
  assert(d->pos[Y] >= -0.5*p->vPatchDim[Y] && d->pos[Y] <= 0.5*p->vPatchDim[Y]);*/
}

/* generating particles */
static void generate(const PARAMS *p,DATA *d)
{
	NODE *root;
	long i,j,k=0,res;
	int nPass,nRej;
	double Nb,Ns,unidev;
	randSeedGenerator((Ullong) (time(NULL) % getpid() + getppid()));
  
	/* pregenerate particle radii */
	res=(long)p->nData; /* set resolution to total number of particles (i.e. each particle has a unique radius)
			      Seems to work best for smooth distribution. If this is optimal, then remove second for loop
			      Currently keeping this here if a resolution change is desired in future (possibly as user-input)
			      However, may be best to remove resolution parameter completely so as to remove possibility of resolution > N_particles*/
  
  /*****BEGIN POWER LAW SIZE DISTRIBUTION *********/
  if (p->dRmin != p->dRmax) {
    double dPRmax = p->dRmax/p->dRmin;
    
    if (p->bSmooth == 1.0){
      if (p->dSDExp != -1){
        for (i=0;i<res;i++){
          for (j=k;j<p->nData*(i+1)/res;j++){
            d[j].keep=0;
            do{
			  unidev=(float) (i+1) / (float) res;
              d[j].radius=pow(((1-unidev)*pow(p->dRmin,p->dSDExp +1))+(unidev*pow(dPRmax*p->dRmin,p->dSDExp+1)),1/(p->dSDExp+1));
            } while (d[j].radius <= 0.0);
          }
          k=j;
        }}
      
      if (p->dSDExp == -1){
        for (i=0;i<res;i++){
          for (j=k;j<p->nData*(i+1)/res;j++){
            d[j].keep=0;
            do{
			  unidev=(float) (i+1) / (float) res;
              d[j].radius=pow(dPRmax*p->dRmin,unidev)/pow(p->dRmin,unidev-1);
            } while (d[j].radius <= 0.0);
          }
          k=j;
        }}
    }
    
    if (p->bSmooth == 0.0){
      if (p->dSDExp != -1){
        for (i=0;i<p->nData;i++){
          d[i].keep = 0;
          do{
            unidev=randUniform();
            d[i].radius=pow(((1-unidev)*pow(p->dRmin,p->dSDExp +1))+(unidev*pow(dPRmax*p->dRmin,p->dSDExp+1)),1/(p->dSDExp+1));
          } while (d[i].radius <= 0.0); /* reject unphysical values */
        }}
      
      if (p->dSDExp == -1){
        for (i=0;i<p->nData;i++){
          d[i].keep = 0;
          do{
            unidev=randUniform();
            d[i].radius=pow(dPRmax*p->dRmin,unidev)/pow(p->dRmin,unidev-1);
          } while (d[i].radius <= 0.0); /* reject unphysical values */
        }}
    }
  }
  else {
    for (i=0;i<p->nData;i++) {
      d[i].keep = 0;
      do {
		d[i].radius = p->dRavg; /*+ gaussian(p->dRadDev,p->dCutoff);*/
      } while (d[i].radius <= 0.0); /* reject unphysical values */
    }
  }
  /*****END POWER LAW SIZE DISTRIBUTION *********/
  
  
  /*****BEGIN SETTING MASS BASED ON DENSITY DISTRIBUTION *********/
  if (p->bDensDist == 1) {
    Ns = p->dDensityMinFrac*p->nData;
    Nb = p->nData-Ns;
    }
  
  /* Initializing Radii and Masses based on Density Distributions */
  if (p->bDensDist == 0) {
    for (i=0;i<p->nData;i++) {
      d[i].keep = 0;
      do {
		/*d[i].radius = p->dRavg;*//* + gaussian(p->dRadDev,p->dCutoff);*/
        d[i].mass = (4.0/3)*PI*p->dDensityAvg*CUBE(d[i].radius);
        } while (d[i].mass <= 0.0); /* reject unphysical values */
      }
    }
  else if (p->bDensDist == 1) {
    for (i=0;i<Ns;i++) {
      d[i].keep = 0;
      do {
		/*d[i].radius = p->dRavg;*//* + gaussian(p->dRadDev,p->dCutoff);*/
        d[i].mass = (4.0/3)*PI*p->dDensityMin*CUBE(d[i].radius);
        } while (d[i].mass <= 0.0); /* reject unphysical values */
      }
    for (j=i;j<Ns+Nb;j++) {
      d[j].keep = 0;
      do {
		/*d[j].radius = p->dRavg;*//* + gaussian(p->dRadDev,p->dCutoff);*/
        d[j].mass = (4.0/3)*PI*p->dDensityAvg*CUBE(d[j].radius);
        } while (d[j].mass <= 0.0); /* reject unphysical values */
      }
    }
  /*****END MASS DISTRIBUTION *********/
  
  
	/*Speeds, Spin, and org_idx data initialization*/
  for (i=0;i<p->nData;i++) {
    ZERO_VEC(d[i].spin);
    d[i].org_idx = (int) i;
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
    }

  /*sorting in size order speeds up the rejection algorithm*/
  qsort(d,(size_t) p->nData,sizeof(DATA),compar);
  
    /* rejection algorithm */
    nPass = 0;
	do {
		if (++nPass > MAX_NUM_PASS_DFLT) {
			fprintf(stderr,"Unable to converge after %i pass%s.\n",
					nPass - 1,nPass == 2 ? "" : "es");
			exit(1);
			}
    printf("Pass %i ",nPass);
    if (nPass == 1)
      printf("(nData = %i)\n",p->nData);
    else
      printf("(nRej = %i)\n",nRej);
    /* make root cell */
    make_node(-0.5*p->vPatchDim[0],0.5*p->vPatchDim[0],
          -0.5*p->vPatchDim[1],0.5*p->vPatchDim[1],
          -0.5*p->vPatchDim[2],0.5*p->vPatchDim[2],&root);
    /* generate particle positions (as needed) and add to tree */
    for (i=0;i<p->nData;i++) {
      if (!d[i].keep) {
        get_pos(p,&d[i]);
        d[i].keep = 1;
        }
      add_to_tree(root,&d[i]);
      }
    /* check for overlaps, rejecting smallest particles first */
    nRej = 0;
    for (i=0;i<p->nData;i++)
      if (reject(root,&d[i])) {
        d[i].keep = 0;
        ++nRej;
        }
    kill_node(root);
    } while (nRej > 0);
  
  /*Apply Shear Flow to final particle position */
  for (i=0;i<p->nData;i++) {
    d[i].vel[Y] -= 1.5*p->dOmega*d[i].pos[X];
    }
  if (p->bReverseSort) /* sort particles by decreasing radius */
	  qsort(d,(size_t) p->nData,sizeof(DATA),rev_compar);
}
/***** END PARTICLE GENERATING ROUTINES ****/


/***** START IO ROUTINES ***/
static void write_data(const PARAMS *p,const DATA *d)
{
	/* these data-structures are defined in ssio.h */
  SSIO ssio;
  SSHEAD head;
  SSDATA data;
  int i;
  /*checks that ssgen.ss is open for writing?*/
  if (ssioOpen(p->achOutputFile,&ssio,SSIO_WRITE)) {
    fprintf(stderr,"Unable to open \"%s\" for writing.\n",p->achOutputFile);
    exit(1);
		}
  
  head.time = 0.0;
  head.n_data = p->nData;
  head.iMagicNumber = SSIO_MAGIC_STANDARD;
  /*makes sure header is properly written (ss standard header?) ?*/
  if (ssioHead(&ssio,&head)) {
    fprintf(stderr,"Unable to write header.\n");
    ssioClose(&ssio);
    exit(1);
		}
  /*writing data... to ssdata data structure*/
  /*using top defined DATA structure to define SSDATA members (named data.), defined within this function*/
  for (i=0;i<p->nData;i++) {
	data.mass = d[i].mass; /*mass given by particle radius and density (given by user)*/
    data.radius = d[i].radius;
    data.pos[0] = d[i].pos[0];
    data.pos[1] = d[i].pos[1];
    data.pos[2] = d[i].pos[2];
    data.vel[0] = d[i].vel[0];
    data.vel[1] = d[i].vel[1];
    data.vel[2] = d[i].vel[2];
    data.spin[0] = d[i].spin[0];
    data.spin[1] = d[i].spin[1];
    data.spin[2] = d[i].spin[2];

    data.color = p->iParticleColor;
    data.org_idx = i;
    if (ssioData(&ssio,&data)) {
      fprintf(stderr,"Error writing data (particle %i).\n",i);
      ssioClose(&ssio);
      exit(1);
    }
		}
  
  ssioClose(&ssio);
}

static double sd_mom(double dRmin,double dRmax,double dAlpha,int n)
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

static void get_params(PARAMS *p)
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
  ReadInt("Particle color", &p->iParticleColor);
  assert(p->iParticleColor >= 0 && p->iParticleColor <= 255);       
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
  ReadInt("Particle density distribution",&p->bDensDist);
  assert(p->bDensDist == 0 || p->bDensDist == 1);
  ReadDbl("Minimum particle density",&p->dDensityMin);
  p->dDensityMin *= 1000/D_SCALE;
  assert(p->dDensity != 0 || (p->dDensityMin > 0 && p->dDensityMin < p->dDensity) || (p->dDensity == 0 && p->bDensDist == 0));
  ReadDbl("Minimum density fraction",&p->dDensityMinFrac);
  assert((p->dDensityMinFrac > 0.0 && p->dDensityMinFrac < 1.0) || (p->bDensDist==0));
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
  
  if (p->bDensDist == 0 && p->dDensity != 0) {
    p->dDensityAvg = p->dDensity;
    }
  else if (p->bDensDist == 1) {
    p->dDensityAvg = (p->dDensityMin*p->dDensityMinFrac)+((1-p->dDensityMinFrac)*p->dDensity);
    }
  
  if (p->iVerbosity > LoVerb)
    (void) printf("Mean %sparticle radius = %g m\n",
                  p->dRScaling != 1 ? "(unscaled) " : "",p->dRavg*L_SCALE);
  if (p->dTau == 0) {
    p->dTau = 0.75*p->dSurfDen/p->dDensityAvg*(p->dR2avg/p->dR3avg);
    if (p->iVerbosity > LoVerb)
      (void) printf("Dynamical optical depth = %g\n",p->dTau);
		}
  else if (p->dSurfDen == 0) {
    p->dSurfDen = 4.0/3*p->dDensityAvg*p->dTau*(p->dR3avg/p->dR2avg);
    if (p->iVerbosity > LoVerb)
      (void) printf("Surface density = %g kg/m^2\n",
                    p->dSurfDen*M_SCALE/SQ(L_SCALE));
		}
  else if (p->dDensity == 0) {
    p->dDensityAvg = 0.75*p->dSurfDen/p->dTau*(p->dR2avg/p->dR3avg);
    if (p->iVerbosity > LoVerb)
      (void) printf("Particle density = %g g/cc\n",
                    p->dDensity*D_SCALE/1000);
		}
  p->dMavg = 4.0/3*PI*p->dDensityAvg*p->dR3avg;
  if (p->iVerbosity > LoVerb)
    (void) printf("Mean %sparticle mass = %g kg\n",
                  p->dMScaling != 1 ? "(unscaled) " : "",p->dMavg*M_SCALE);
  p->dRHavg =
		p->dOrbDist*pow(8*PI*p->dDensityAvg/(9*p->dCentralMass),1.0/3)*p->dRavg;
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
  /*TimeStep Determined for HardSphere*/
  p->dt = 0.01/sqrt(p->dDensityAvg); /*DEBUG! was 0.03*/
  if (p->dt > 0.01*TWO_PI/p->dOmega) /* 1/100th of an orbit (conservative) */
    p->dt = 0.01*TWO_PI/p->dOmega;	
  /*SSDEM Timestep - - - and Kn and Etc.*/
  p->dKn = p->dMavg * SQ(4*p->dOmega*p->dRavg/(0.01*p->dRavg)); /*0.01 is the maximum fractional overlap*/
  p->dDelta = M_PI*sqrt(0.5*p->dMavg/p->dKn)/30; /*30 is the number of steps per overlap*/
  /** END SSDEM parameters **/
  if (p->iVerbosity > LoVerb)
    (void) printf("Recommended SSDEM timestep = %g (1 orbit = %g steps)\n",
                  p->dDelta,TWO_PI/(p->dOmega*p->dDelta));
  if (p->iVerbosity > LoVerb && p->dMScaling != 1)
    (void) printf("Particle masses will be scaled by %g\n",p->dMScaling);
  if (p->iVerbosity > LoVerb && p->dRScaling != 1)
    (void) printf("Particle radii will be scaled by %g\n",p->dRScaling);
  p->dRavg *= p->dRScaling;
  p->dR2avg *= SQ(p->dRScaling);
  p->dR3avg *= CUBE(p->dRScaling);
  p->dMavg *= p->dMScaling;
  p->dVesc = p->dRavg*sqrt(8.0/3*PI*p->dDensityAvg*
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

static void
write_log(const PARAMS *p)
{
  FILE *fp;
  
  fp = fopen(LOG_FILE,"w");
  assert(fp != NULL);
  (void) fprintf(fp,"Number of particles\t%i\n",p->nData);
  (void) fprintf(fp,"Number of SSDEM steps\t%i\n",
                 (int)(TWO_PI*p->dTime/(p->dOmega*p->dDelta) + 0.5));
  (void) fprintf(fp,"SSDEM Timestep\t\t%.16e\t# %g s\n",p->dDelta,p->dDelta*T_SCALE);
  (void) fprintf(fp,"Recommended Kn\t\t%.16e\t# %g kg/s^2\n",p->dKn,p->dKn*M_SCALE/(T_SCALE*T_SCALE));
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
  (void) fprintf(fp,"Mean particle density\t%.16e\t# %g g/cc\n",p->dDensityAvg,
                 p->dDensityAvg*D_SCALE/1000);
  (void) fprintf(fp,"Critical wavelength\t%.16e\t# %g m\n",p->dLambdaCrit,
                 p->dLambdaCrit*L_SCALE);
  (void) fprintf(fp,"Volume filling factor\t%.16e\t# %g\n",p->dVFF,p->dVFF);
  (void) fprintf(fp,"Mass scaling\t\t%.16e\t# %g\n",p->dMScaling,
                 p->dMScaling);
  (void) fprintf(fp,"Radius scaling\t\t%.16e\t# %g\n",p->dRScaling,
                 p->dRScaling);
  (void) fclose(fp);
}
/***** END IO ROUTINES ****/

int main(int argc,char *argv[])
{
  PARAMS params;
  DATA *data = NULL;
  
  BOOLEAN bAdjCom=FALSE,bForce=FALSE,bUsage=FALSE,bReverseSort=FALSE;
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
        bReverseSort = TRUE;
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
  params.bReverseSort = bReverseSort;
  
  if (!bForce) {
    FILE *fp = fopen(params.achOutputFile,"r");
    if (fp) {
      (void) fprintf(stderr,"%s exists -- use \"-f\" to overwrite.\n",
                     params.achOutputFile);
      (void) fclose(fp);
      return 1;
      }
		}
  
  /* Generate initial conditions */
  
  if (!(data = (DATA *) malloc((size_t) params.nData*sizeof(DATA)))) {
    (void) fprintf(stderr,"Unable to allocate data memory.\n");
    return 1;
		}
  
  generate(&params,data);
  if (params.bAdjCom) adj_com(&params,data);
  /*Writing PATCHIC LOGFILE _FIX*/
  write_log(&params);
  
  /*Screen Output*/
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
  write_data(&params,data);
  
  /* All done */
  
  free((void *)data);
  
  return 0;
}
