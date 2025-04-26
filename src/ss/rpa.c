/*
 ** rpa.c -- DCR 09/18/01
 ** =====
 ** Rubble pile analyzer Mark II (tree version).
 **
 ** RP: Added capability for local (patch) geometry -- Aug 2009 
 ** RP: Added option to calculate rp radius as a volume-equivalent
 **  sphere (R = third root of (abc)), not a minimally-enclosing 
 **  sphere (R = a). -- Sept 2009
 */

#include <rpu.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <math.h>
#include <assert.h>
#include <rdpar.h> /* for parameter-reading functions (OpenPar... etc) */

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

/*#define COLOR_BY_DIST*/

#define DFLT_LINKING_SCALE            1.1
#define DFLT_THETA                    0.5
#define DFLT_NUM_BIGGEST              1
#define DFLT_DO_END_STATS             FALSE
#define DFLT_SAVE_GROUPS_AS_AGGS      FALSE
#define DFLT_DO_MIXING_STATS          FALSE
#define DFLT_SAVE_BIGGEST_IN_ONE_FILE FALSE
#define DFLT_LINK_BY_PARTICLE         FALSE
#define DFLT_MERGE_AGGS               FALSE
#define DFLT_STOP_AT_AGGS             FALSE
#define DFLT_SAVE_SINGLES             TRUE
#define DFLT_SAVE_DOUBLES             TRUE
#define DFLT_DO_NOT_RECOLOR           FALSE
#define DFLT_USE_RESERVED_COLORS      TRUE
#define DFLT_VERBOSE                  FALSE
#define DFLT_DOPATCH                  FALSE
#define DFLT_ALTERNATER               FALSE
#define DFLT_ERROR                    FALSE

#define NUM_MIXING_CALCS 100

#define SUMMARY    "rpa.out"
#define GROUPS_SS  "rpa_groups.ss"
#define CENTERS_BT "rpa_centers.bt"
#define CENTERS_DD "rpa_centers.dd"
#define CENTERS_OE "rpa_centers.oe"

#define BIGGEST_IN_ONE_FILE_SS	"rpa_biggest.ss"
#define BIGGEST_ONE_PER_FILE_SS	"rpa_biggest%i.ss"

#define PATCH_LOG "patchic.log"

typedef struct
{
  double dLength;
  double dWidth;
  double dOrbFreq;
}  PATCH_PARAMS;

typedef struct {
  double dLinkingScale;
  double dThetaSq;
  double dTime; /* from data file */
  int nBiggest;
  BOOLEAN bDoEndStats;
  BOOLEAN bSaveGroupsAsAggs;
  BOOLEAN bDoMixingStats;
  BOOLEAN bSaveBiggestInOneFile;
  BOOLEAN bLinkByParticle;
  BOOLEAN bMergeAggs;
  BOOLEAN bSaveSingles;
  BOOLEAN bSaveDoubles;
  BOOLEAN bDoNotRecolor;
  BOOLEAN bUseReservedColors;
  BOOLEAN bVerbose;
  BOOLEAN bLastFile;
  BOOLEAN bDoPatch;
  BOOLEAN bAlternateR;
  BOOLEAN bError;
  } PARAMS;

typedef struct {
  double a,e,q;
  } OSC_ELEM;

typedef struct {
  int i;
#ifdef COLOR_BY_DIST
  double r;
#endif
  } DSS;

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  VECTOR pos;
  double size,half_size,eff_half_size,eff_size_sq;
  struct node_s *cell[CELLS_PER_NODE];
  RUBBLE_PILE *leaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

static void write_ss(const PARAMS *p,const char *filename,RUBBLE_PILE *rp,int skip,int n)
{
	SSIO ssio;
	SSHEAD h;
	int i,j,n_data;

	if (p->bVerbose)
		(void) printf("\nWriting to \"%s\"...\n",filename);

	if (ssioOpen(filename,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\" for writing\n",filename);
		return;
		}

	for (n_data=0,i=skip;i<skip+n;i++) {
		if (!p->bSaveSingles && rp[i].n_particles == 1)
			continue;
		if (!p->bSaveDoubles && rp[i].n_particles == 2)
			continue;
		n_data += rp[i].n_particles;
		}
	h.time = 0;
	h.n_data = n_data;
	h.iMagicNumber = SSIO_MAGIC_STANDARD;

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"Error writing header\n");
		return;
		}

	for (i=skip;i<skip+n;i++) {
		assert(rp[i].data != NULL);
		if (!p->bSaveSingles && rp[i].n_particles == 1)
			continue;
		if (!p->bSaveDoubles && rp[i].n_particles == 2)
			continue;
		for (j=0;j<rp[i].n_particles;j++) {
			if (p->bSaveGroupsAsAggs)
				AGG_SET_IDX(&rp[i].data[j],i);
			if (ssioData(&ssio,&rp[i].data[j])) {
				(void) fprintf(stderr,"Error writing data\n");
				return;
				}
			}
		}

	(void) ssioClose(&ssio);
	
	if (p->bVerbose)
	  (void) printf("%i particle%s written.\n",n_data,n_data==1?"":"s");
}

/* Obtains the necessary patch geometry parameters from a file */
static int get_params(double* dOrbFreq, double* dLx, double* dLy)
{
  OpenPar(PATCH_LOG);
   
  ReadDbl("Orbital frequency",dOrbFreq); /* read in radians*2pi/yr */
  assert(*dOrbFreq > 0.0);

  ReadDbl("Patch width",dLx); /* read in AU */
  assert(*dLx > 0.0);

  ReadDbl("Patch length",dLy); /* read in AU */
  assert(*dLy > 0.0);

  ClosePar();

  return 0;
}

/* Stolen from pkd.h */
#define SHEAR(ix,t,pp)\
	((ix) < 0 ? fmod(0.5*(pp)->dLength - 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength) - 0.5*(pp)->dLength:\
	 (ix) > 0 ? 0.5*(pp)->dLength - fmod(0.5*(pp)->dLength + 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength): 0.0)

/* Edits particle data in light of patch geometry:
     -Convert spins from 'space' frame into the patch frame (to be consistent with particle velocities). 
     -"Unwrap" particle positions (and y-velocities) so that aggregates
     are grouped in space, rather than within patch boundaries.
*/
static int convert_to_patch_coords(PARAMS *p,RUBBLE_PILE *wrp)
{
  int i,j; /* loop variables */
  int thisAgg; /* stores an aggregate's index */
  PATCH_PARAMS PP; /* Stores patch geometry params - for 'SHEAR' */
  double cosOmegaTime,sinOmegaTime,wx,wy;

  get_params(&PP.dOrbFreq,&PP.dWidth,&PP.dLength);

  if (p->bVerbose)
    (void) printf("Patch parameters: omega %g, Lx %g, Ly %g \n",
		  PP.dOrbFreq,PP.dWidth,PP.dLength);

  /* Prepare to rotate all the particle spins by the same angle, -Omega*time.... */
  /* Stolen from PKDGRAV: pkdDoCollision(). */
  cosOmegaTime = cos(-PP.dOrbFreq*p->dTime);
  sinOmegaTime = sin(-PP.dOrbFreq*p->dTime);
  
  for(i=0; i<wrp->n_particles; i++)
    {
      /* Rotate spin axes */
      wx =  cosOmegaTime*wrp->data[i].spin[0] - sinOmegaTime*wrp->data[i].spin[1];
      wy =  sinOmegaTime*wrp->data[i].spin[0] + cosOmegaTime*wrp->data[i].spin[1];
      wrp->data[i].spin[0] = wx;
      wrp->data[i].spin[1] = wy;
      
      /* Compensate for rotating coordinate system */
      wrp->data[i].spin[2] -= PP.dOrbFreq;
      
      /* Loop through particles and unwrap aggregates
	 NOTE: This assumes that "aggregates" share a common index!! I.e., cannot do "rubble piles" */
      thisAgg = AGG_IDX(&wrp->data[i]);
      if(thisAgg >= 0)
	for(j=0; j<wrp->n_particles; j++) {
	  if(thisAgg == AGG_IDX(&wrp->data[j]) && i != j)
	    {
	      /* Stolen from pkdgrav: pkdAggsInPatchGetUnwrapped(),
		 with modifications to use particle i as the reference
		 point. */
	      /* we're assuming a 2-D patch here, with origin at (0,0,0)... */
	      if (wrp->data[j].pos[0] - wrp->data[i].pos[0] > 0.5*PP.dWidth){
		wrp->data[j].pos[0] -= PP.dWidth;
		wrp->data[j].pos[1] += SHEAR(-1,p->dTime,&PP);
		wrp->data[j].vel[1] += 1.5*PP.dOrbFreq*PP.dWidth;
	      }
	      else if (wrp->data[i].pos[0] - wrp->data[j].pos[0] > 0.5*PP.dWidth){
		wrp->data[j].pos[0] += PP.dWidth;
		wrp->data[j].pos[1] += SHEAR(1,p->dTime,&PP);
		wrp->data[j].vel[1] -= 1.5*PP.dOrbFreq*PP.dWidth;
	      }
	      /* finally, check (newly unwrapped) y-position wrt
		 reference. */
	      if (wrp->data[j].pos[1] - wrp->data[i].pos[1] > 0.5*PP.dLength){
		wrp->data[j].pos[1] -= PP.dLength;
	      }
	      else if (wrp->data[i].pos[1] - wrp->data[j].pos[1] > 0.5*PP.dLength){
		wrp->data[j].pos[1] += PP.dLength;
	      }
	    }
	} /* loop over j */
    } /* loop over i */
  return 0;
}


/* Edits particle data in light of patch geometry - this time, converting from patch frame to ss-format:
     -Convert spins from 'patch' frame into the space frame (to be consistent with ss format choice). 
     -"Wrap" particle positions (and y-velocities) so that all particles are located within patch boundaries.
*/
static int convert_to_space_coords(PARAMS *p,RUBBLE_PILE *wrp)
{
  int i; /* loop variable */
  PATCH_PARAMS PP; /* Stores patch geometry params - for 'SHEAR' */
  double cosOmegaTime,sinOmegaTime,wx,wy;
  
  get_params(&PP.dOrbFreq,&PP.dWidth,&PP.dLength);

  /* Prepare to de-rotate all the particle spins by the same angle, Omega*time.... */
  /* Stolen from PKDGRAV: pkdDoCollision(). */
  cosOmegaTime = cos(PP.dOrbFreq*p->dTime);
  sinOmegaTime = sin(PP.dOrbFreq*p->dTime);
  
  for(i=0; i<wrp->n_particles; i++)
    {
      /* Rotate spin axes */
      wx =  cosOmegaTime*wrp->data[i].spin[0] - sinOmegaTime*wrp->data[i].spin[1];
      wy =  sinOmegaTime*wrp->data[i].spin[0] + cosOmegaTime*wrp->data[i].spin[1];
      wrp->data[i].spin[0] = wx;
      wrp->data[i].spin[1] = wy;
      
      /* Compensate for rotating coordinate system */
      wrp->data[i].spin[2] += PP.dOrbFreq;
      
      /* Loop through particles and wrap particles */

      /* we're assuming a 2-D patch here, with origin at (0,0,0)... */
      if (wrp->data[i].pos[0] > 0.5*PP.dWidth){
	wrp->data[i].pos[0] -= PP.dWidth;
	wrp->data[i].pos[1] += SHEAR(-1,p->dTime,&PP);
	wrp->data[i].vel[1] += 1.5*PP.dOrbFreq*PP.dWidth;
      }
      else if (wrp->data[i].pos[0] < 0.5*PP.dWidth){
	wrp->data[i].pos[0] += PP.dWidth;
	wrp->data[i].pos[1] += SHEAR(1,p->dTime,&PP);
	wrp->data[i].vel[1] -= 1.5*PP.dOrbFreq*PP.dWidth;
      }
      if (wrp->data[i].pos[1] > 0.5*PP.dLength){
	wrp->data[i].pos[1] -= PP.dLength;
      }
      else if (wrp->data[i].pos[1] < 0.5*PP.dLength){
	wrp->data[i].pos[1] += PP.dLength;
      }
    } /* loop over i */
  return 0;
}

#ifdef COLOR_BY_DIST

static int dist_sort(const void *v1,const void *v2)
{
	DSS *d1 = (DSS *) v1;
	DSS *d2 = (DSS *) v2;

	return ((d1->r < d2->r) ? -1 : (d1->r > d2->r) ? 1 : 0);
	}

#endif

static void recolor(const PARAMS *p,const RUBBLE_PILE *wrp,RUBBLE_PILE *rp,int n)
{
	/*DEBUG used to be that dominant color would be used: ie. color
	  of rubble pile = most abundant particle color in rubble pile,
	  but now either color by mass or distance from center of mass*/

	DSS *dss;
	int i;

	dss = (DSS *) malloc((size_t) n*sizeof(DSS));
	assert(dss != NULL);

	for (i=0;i<n;i++)
		dss[i].i = i;

#ifdef COLOR_BY_DIST
	{
	VECTOR v;

	for (i=0;i<n;i++) {
		SUB_VEC(rp[i].pos,wrp->pos,v);
		dss[i].r = MAG_SQ(v);
		}

	qsort((void *) dss,n,sizeof(DSS),dist_sort);
	}
#else
	(void) wrp; /* to suppress unused parameter warning */
#endif

	for (i=0;i<n;i++)
		if (rp[dss[i].i].n_particles == 1)
			rp[i].color = WHITE;
		else if (rp[dss[i].i].n_particles == 2)
			rp[i].color = GREEN;
		else if (p->bUseReservedColors) {
			switch (i%13) {
			case 0:
				rp[dss[i].i].color = RED;
				break;
			case 1:
				rp[dss[i].i].color = BLUE;
				break;
			case 2:
				rp[dss[i].i].color = YELLOW;
				break;
			case 3:
				rp[dss[i].i].color = MAGENTA;
				break;
			case 4:
				rp[dss[i].i].color = CYAN;
				break;
			case 5:
				rp[dss[i].i].color = GOLD;
				break;
			case 6:
				rp[dss[i].i].color = PINK;
				break;
			case 7:
				rp[dss[i].i].color = ORANGE;
				break;
			case 8:
				rp[dss[i].i].color = KHAKI;
				break;
			case 9:
				rp[dss[i].i].color = VIOLET;
				break;
			case 10:
				rp[dss[i].i].color = MAROON;
				break;
			case 11:
				rp[dss[i].i].color = AQUA;
				break;
			case 12:
				rp[dss[i].i].color = NAVY;
				break;
			default:
				assert(0);
				}
			}
		else {
			switch (i%11) { /* RED & YELLOW reserved */
			case 0:
				rp[dss[i].i].color = BLUE;
				break;
			case 1:
				rp[dss[i].i].color = MAGENTA;
				break;
			case 2:
				rp[dss[i].i].color = CYAN;
				break;
			case 3:
				rp[dss[i].i].color = GOLD;
				break;
			case 4:
				rp[dss[i].i].color = PINK;
				break;
			case 5:
				rp[dss[i].i].color = ORANGE;
				break;
			case 6:
				rp[dss[i].i].color = KHAKI;
				break;
			case 7:
				rp[dss[i].i].color = VIOLET;
				break;
			case 8:
				rp[dss[i].i].color = MAROON;
				break;
			case 9:
				rp[dss[i].i].color = AQUA;
				break;
			case 10:
				rp[dss[i].i].color = NAVY;
				break;
			default:
				assert(0);
				}
			}
	
	for (i=0;i<n;i++)
		rpuApplyColor(&rp[i]);

	free((void *) dss);
	}

static int get_mixing(const double *wmc,double wm,const RUBBLE_PILE *rp,double *mixing)
{
	/*
	 ** Computes a "mixing" statistic for the distribution of world particle
	 ** colors, weighted by their world mass fraction, in given rubble pile.
	 ** A mixing value of 1 implies the rubble pile contains a perfectly
	 ** homogeneous mixture of the world colors.
	 */

	MATRIX rot;
	VECTOR pos0,pos;
	double v,rs,mp,sum,sum_sq[NUM_COLORS],mc[NUM_COLORS];
	int nmin,ns,np,i,j;

	assert(wmc != NULL && wm > 0.0 && rp != NULL && mixing != NULL);

	*mixing = -1;

	/* Set size of sample region to contain on average N^1/2 particles */

	v = rpuVolEll(rp->axis_len);
	rs = pow(0.75*v/(PI*sqrt(rp->n_particles)),1.0/3);
	nmin = (int) sqrt(sqrt(rp->n_particles));

	/* Make sure every region is sampled once (on average) */

	ns = (int) (v/rpuVolSph(rs)); /* same as sqrt(rp->n_particles) */

	/* Sanity check */

	if (rs > rp->axis_len[MINOR(rp)] || nmin < 2 || ns < 1)
		return 1; /* cannot reliably determine mixing */

	/* Accumulate statistics */

	for (i=0;i<NUM_COLORS;i++)
		sum_sq[i] = 0;

	COPY_MAT(rp->axes,rot);
	Transpose(rot); /* for body-to-Cartesian-axes conversion */

	for (i=0;i<ns;i++) {
		do {
			/*
			 ** Choose a random location within the rectangular prism
			 ** enclosing the rubble pile. At least "nmin" particles
			 ** must be found.
			 */
			for (j=0;j<N_DIM;j++)
				pos0[j] = (2.0*((double) rand()/RAND_MAX) - 1.0)*(rp->axis_len[j] - rs);
			Transform(rot,pos0,pos);
			ADD_VEC(pos,rp->pos,pos);
			/* Sum particle masses by color in the sample region */
			for (j=0;j<NUM_COLORS;j++)
				mc[j] = 0;
			mp = np = 0;
			for (j=0;j<rp->n_particles;j++) {
				SUB_VEC(rp->data[j].pos,pos,pos0);
				if (MAG(pos0) <= rs) {
					mc[rp->data[j].color] += rp->data[j].mass;
					mp += rp->data[j].mass;
					np += 1;
					}
				}
			} while (np < nmin);
		assert(mp > 0.0);
		for (j=0;j<NUM_COLORS;j++) {
			sum = mc[j]/mp - wmc[j]/wm;
			sum_sq[j] += SQ(sum);
			}
		}

	*mixing = 1;
	for (i=0;i<NUM_COLORS;i++)
		*mixing -= sqrt(sum_sq[i]/ns); /* can suffer from round-off */

	return 0;
	}

static int mass_sort(const void *v1,const void *v2)
{
	const RUBBLE_PILE *rp1 = (const RUBBLE_PILE *) v1;
	const RUBBLE_PILE *rp2 = (const RUBBLE_PILE *) v2;

	return ((rp1->mass > rp2->mass) ? -1 : (rp1->mass < rp2->mass) ? 1 : 0);
	}

static BOOLEAN ok_to_merge(RUBBLE_PILE *rp1,RUBBLE_PILE *rp2)
{
	/*
	 ** The following merger strategy is based entirely on geometry and
	 ** does not take the gravitational potential into account. It is
	 ** suitable for searching from the bottom up, that is, for starting
	 ** with individual particles and linking them together into larger
	 ** groups. The search takes the ellipsoidal shapes of the current
	 ** groups into account. In order for rp1 to be merged with rp2,
	 ** spheres drawn with radii equal to the major axes of the bodies
	 ** (times linking-scale) and centred on the bodies must overlap.
	 ** If the scaled minor spheres also overlap, the bodies are merged.
	 ** Failing that, if either body has its centre of mass in the other's
	 ** scaled ellipsoid, the bodies are merged. Otherwise no merge occurs.
	 */

	double r1,r2;

	r1 = rp1->axis_len[MAJOR(rp1)];
	r2 = rp2->axis_len[MAJOR(rp2)];

	/* If scaled major spheres overlap, could have merger */

	if (OVERLAP(rp1->pos,r1,rp2->pos,r2)) {
		VECTOR r,v;

		/* If both "rubble piles" are single particles, have merger */

		if (rp1->n_particles == 1 && rp2->n_particles == 1)
			return TRUE;

		r1 = rp1->axis_len[MINOR(rp1)];
		r2 = rp2->axis_len[MINOR(rp2)];

		/* If scaled minor spheres overlap, definitely have merger */

		if (OVERLAP(rp1->pos,r1,rp2->pos,r2))
			return TRUE;

		/* Otherwise check for one group inside the other */

		SUB_VEC(rp2->pos,rp1->pos,r);
		Transform(rp1->axes,r,v);
		if (rpuInEllipsoid(v,0.0,rp1->axis_len))
			return TRUE;

		SUB_VEC(rp1->pos,rp2->pos,r);
		Transform(rp2->axes,r,v);
		if (rpuInEllipsoid(v,0.0,rp2->axis_len))
			return TRUE;
		}

	return FALSE;
	}

static void find_closest(const NODE *node,double dThetaSq,RUBBLE_PILE *rpi,
						 RUBBLE_PILE **rpc,double *dDistSq)
{
	/*
	 ** Identifies closest rubble pile to rpi, and its distance.
	 */

	VECTOR v;
	double r2;
	int i;

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i]) {
			SUB_VEC(node->cell[i]->pos,rpi->pos,v);
			if (node->cell[i]->eff_size_sq/MAG_SQ(v) > dThetaSq)
				find_closest(node->cell[i],dThetaSq,rpi,rpc,dDistSq);
			}
		else if (node->leaf[i] && node->leaf[i] != rpi &&
				 node->leaf[i]->data) {
			SUB_VEC(node->leaf[i]->pos,rpi->pos,v);
			if ((r2 = MAG_SQ(v)) < *dDistSq) {
				*rpc = node->leaf[i];
				*dDistSq = r2;
				}
			}
	}

static void find_nbrs(const PARAMS *p,const NODE *node,const RUBBLE_PILE *rp,
					  const SSDATA *d,RUBBLE_PILE *rpn[],int *nNbr,int nMax)
{
	/*
	** Returns a list of neighbours rpn (pointers to rubble pile
	** structs) that are within the linking distance of the given
	** particle (rubble pile rp, data d).
	*/

	const double sqrt3 = sqrt(3.0);

	VECTOR v;
	double r2;
	int i;

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i] != NULL) {
			SUB_VEC(node->cell[i]->pos,d->pos,v);
			r2 = MAG_SQ(v);
			/* opening test checks for cell intersection as well */
			if (node->cell[i]->eff_size_sq/r2 > p->dThetaSq ||
				r2 <= SQ(p->dLinkingScale*(sqrt3*node->cell[i]->half_size + d->radius)))
				find_nbrs(p,node->cell[i],rp,d,rpn,nNbr,nMax);
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != rp &&
				 node->leaf[i]->data != NULL && node->leaf[i]->data != d) {
			if (node->leaf[i]->n_particles > 1)
				continue; /* can ignore already-formed links -- they would have found us by now */
			assert(node->leaf[i]->n_particles == 1); /* one particle per "rubble pile" to be tested */
			SUB_VEC(node->leaf[i]->data->pos,d->pos,v);
			r2 = MAG_SQ(v);
			if (r2 <= SQ(p->dLinkingScale*(node->leaf[i]->data->radius + d->radius))) {
				assert(*nNbr < nMax);
				rpn[(*nNbr)++] = node->leaf[i];
				}
			}
	}

static void link_part(const PARAMS *p,const NODE *root,RUBBLE_PILE *rp,int n)
{
	/*
	** Strategy: for each particle, determine if there are any other
	** particles within a center-to-center distance of L*(R_i + R_j),
	** where L is the linking scale, R_i is the current particle's
	** radius, and R_j is the radius of a neighbouring particle found
	** by find_nbrs() above.  If there are, these particles get added
	** to R_i's rubble pile structure and become the centers for
	** neighbour searches later on.  In this way, only 1 complete pass
	** over the particle data is needed.
	*/

	RUBBLE_PILE **rpn;
	int i,j,k,nLink=0,nNbr,nMax,nOld;

	if (p->bVerbose)
		(void) printf("Linking particles...\n");

	assert(n > 0);

	if (n == 1)
		goto link_done;

	nMax = n - 1; /* max num neighbours per particle */

	rpn = (RUBBLE_PILE **) malloc((size_t) nMax*sizeof(RUBBLE_PILE *));
	assert(rpn != NULL);

	for (i=0;i<n;i++) {
		if (rp[i].data == NULL)
			continue;
		for (j=0;j<rp[i].n_particles;j++) {
			nNbr = 0;
			find_nbrs(p,root,&rp[i],&rp[i].data[j],rpn,&nNbr,nMax);
			if (nNbr > 0) {
				nLink += nNbr;
				nOld = rp[i].n_particles;
				rp[i].n_particles += nNbr;
				rpuRealloc(&rp[i]);
				for (k=0;k<nNbr;k++) {
					rpuCopyData(rpn[k],&rp[i],(size_t) (nOld + k));
					rpuFree(rpn[k]); /* frees & NULLs rpn[k]->data */
					}
				}
			}
		}

	free((void *) rpn);

 link_done:

	(void) printf("No. particle links = %i.\n",nLink);
	}

static int merge_rp(const NODE *root,double dThetaSq,RUBBLE_PILE *rp,int n)
{
	/*
	 ** Executes one merging pass over rubble piles.
	 */

	RUBBLE_PILE *rpClosest;
	double r2;
	int nMerge,nOld,i;

	nMerge = 0;
	for (i=0;i<n;i++) {
		if (rp[i].data == NULL)
			continue;
		r2 = HUGE_VAL;
		find_closest(root,dThetaSq,&rp[i],&rpClosest,&r2);
		if (r2 == HUGE_VAL)
			continue; /* no neighbour found */
		assert(rpClosest != NULL);
		if (ok_to_merge(&rp[i],rpClosest)) {
			++nMerge;
			/* move data from second rubble pile to first */
			nOld = rp[i].n_particles;
			rp[i].n_particles += rpClosest->n_particles;
			rpuRealloc(&rp[i]);
			rpuCopyData(rpClosest,&rp[i],(size_t) nOld);
			rpuFree(rpClosest);
			}
		}
	return nMerge;
	}

static void make_node(const VECTOR pos,double size,NODE **node)
{
	int i;

	assert(size > 0.0);
	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	COPY_VEC(pos,(*node)->pos);
	(*node)->size = size;
	(*node)->half_size = (*node)->eff_half_size = 0.5*size;
	(*node)->eff_size_sq = SQ(size);

	for (i=0;i<CELLS_PER_NODE;i++) {
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}
	}

static void add_to_tree(NODE *node,RUBBLE_PILE *rp)
{
	/* note: assumes rubble pile inside node! */

	int i,idx,idy,idz;

	idx = (rp->pos[X] < node->pos[X] ? -1 : 1);
	idy = (rp->pos[Y] < node->pos[Y] ? -1 : 1);
	idz = (rp->pos[Z] < node->pos[Z] ? -1 : 1);

	i = (idx + 1)/2 + (idy + 1 + 2*(idz + 1));

	if (node->cell[i])
		add_to_tree(node->cell[i],rp);
	else if (node->leaf[i]) {
		VECTOR v;
		SET_VEC(v,idx,idy,idz);
		SCALE_VEC(v,0.5*node->half_size);
		ADD_VEC(v,node->pos,v);
		make_node(v,node->half_size,&node->cell[i]);
		add_to_tree(node->cell[i],node->leaf[i]);
		add_to_tree(node->cell[i],rp);
		node->leaf[i] = NULL; /* for completeness */
		}
	else
		node->leaf[i] = rp;
	
	if (rp->radius > node->eff_half_size) {
		node->eff_half_size = rp->radius;
		node->eff_size_sq = 4.0*SQ(rp->radius);
		}
	}

static void kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i=0;i<CELLS_PER_NODE;i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}

static int merge_aggs(const PARAMS *p,const RUBBLE_PILE *wrp,RUBBLE_PILE *rp)
{
	int i,iAggMaxIdx,iAggIdx,nAggMax,nAgg,*iAggIdxArr,*nPartInAgg,iIdx,iRpIdx;

	if (p->bVerbose)
		(void) printf("Linking aggregates...\n");

	/* first determine maximum number of aggregates we need to deal with */

	iAggMaxIdx = -1;
	for (i=0;i<wrp->n_particles;i++) {
		iAggIdx = AGG_IDX(&wrp->data[i]);
		if (iAggIdx > iAggMaxIdx)
			iAggMaxIdx = iAggIdx;
		}
	nAggMax = iAggMaxIdx + 1;

	assert(nAggMax >= 0);

	if (nAggMax == 0) {
		if (p->bVerbose)
			(void) printf("No aggregates found.\n");
		return 0; /* nothing to do */
		}

	/* now allocate & initialize space for aggregate index arrays */

	iAggIdxArr = (int *) malloc((size_t) nAggMax*sizeof(int));
	assert(iAggIdxArr != NULL);
	nPartInAgg = (int *) malloc((size_t) nAggMax*sizeof(int));
	assert(nPartInAgg != NULL);

	for (i=0;i<nAggMax;i++) {
		iAggIdxArr[i] = -1;
		nPartInAgg[i] = 0; /* (should be done automatically by malloc()) */
		}

	/* count number of particles in each aggregate */

	for (i=0;i<wrp->n_particles;i++) {
		iAggIdx = AGG_IDX(&wrp->data[i]);
		if (iAggIdx >= 0) /* i.e., particle belongs to aggregate */
			++nPartInAgg[iAggIdx];
		}

	for (nAgg=i=0;i<nAggMax;i++)
		if (nPartInAgg[i] > 0)
			++nAgg;

	assert(nAgg > 0 && nAgg <= wrp->n_particles);

	/* assign particles or aggregates to initial rubble piles */

	for (iIdx=i=0;i<wrp->n_particles;i++) {
		rp[i].n_particles = 0; /* initialize */
		rp[i].data = NULL;
		iAggIdx = AGG_IDX(&wrp->data[i]);
		if (iAggIdx >= 0) {
			if (iAggIdxArr[iAggIdx] == -1) {
				rp[iIdx].n_particles = nPartInAgg[iAggIdx];
				rpuMalloc(&rp[iIdx]);
				rp[iIdx].n_particles = 1; /* for use as counter */
				*rp[iIdx].data = wrp->data[i]; /* struct copy */
				iAggIdxArr[iAggIdx] = iIdx++;
				}
			else {
				/* the following ugly assignment copies the particle
				   data at wrp->data[i] to the end of the data in the
				   rubble pile corresponding to this aggregate */
				iRpIdx = iAggIdxArr[iAggIdx];
				rp[iRpIdx].data[rp[iRpIdx].n_particles++] = wrp->data[i];
				}
			}
		else {
			rp[iIdx].n_particles = 1;
			rpuMalloc(&rp[iIdx]);
			*rp[iIdx++].data = wrp->data[i];
			}
		}

	/* deallocate temporary arrays */

	free((void *) nPartInAgg);
	free((void *) iAggIdxArr);

	if (p->bVerbose)
		(void) printf("%i aggregate%s linked.\n",nAgg,nAgg==1?"":"s");

	return nAgg;
	}

static void find_rp(const PARAMS *p,RUBBLE_PILE *wrp,RUBBLE_PILE *rp)
{
	NODE *root = NULL;
	int i,nAgg=0,nMerge;

	/*
	** If desired, do a first pass putting any aggregate particles
	** into their own rubble piles.  Otherwise initialize with 1
	** particle per rubble pile.
	*/

	if (p->bMergeAggs)
		nAgg = merge_aggs(p,wrp,rp);

	if (nAgg == 0) {
		for (i=0;i<wrp->n_particles;i++) {
			rp[i].n_particles = 1;
			rpuMalloc(&rp[i]);
			*rp[i].data = wrp->data[i];
			}
		if (p->dLinkingScale == 0.0) {
			(void) fprintf(stderr,"WARNING: Group finding terminated but no aggregates found.\n");
			return;
			}
		}

	if (p->dLinkingScale == 0.0)
		return; /* no linking to be performed */

	if (p->bLinkByParticle) {
		/* only 1 iteration needed for this option */
		assert(!p->bMergeAggs); /* pre-linking aggs not supported */
		make_node(wrp->pos,2.0*wrp->radius,&root);
		for (i=0;i<wrp->n_particles;i++) {
			assert(rp[i].data != NULL);
			/* set position & radius data for insertion into tree */
			COPY_VEC(rp[i].data->pos,rp[i].pos);
			rp[i].radius = rp[i].data->radius;
			add_to_tree(root,&rp[i]);
			}
		link_part(p,root,rp,wrp->n_particles);
		goto find_done;
		}

	/*
	 ** Iterate, first replacing any existing tree, then computing
	 ** each rubble pile's info and inserting it into new tree.
	 **/

	if (p->bVerbose)
		(void) printf("Grouping... nMerge = ");

	do {
		if (root) kill_node(root);
/*		(void) printf("Building tree...\n");*/
		make_node(wrp->pos,2.0*wrp->radius,&root);
		for (i=0;i<wrp->n_particles;i++) {
			if (rp[i].data == NULL)
				continue;
			rpuCalcMass(&rp[i]);
			rpuCalcPos(&rp[i]);
			rpuCalcInertia(&rp[i]);
			rpuCalcAxes(&rp[i]);
			SCALE_VEC(rp[i].axis_len,p->dLinkingScale);
			add_to_tree(root,&rp[i]);
			}
/*		(void) printf("Done!\n");*/
		nMerge = merge_rp(root,p->dThetaSq,rp,wrp->n_particles);
		if (p->bVerbose)
			(void) printf("%10i\b\b\b\b\b\b\b\b\b\b",nMerge);
		}
	while (nMerge > 0);

	if (p->bVerbose)
		(void) printf("\n");

 find_done:

	kill_node(root); /* deallocate the tree */
	}

static void analyze_world(PARAMS *p,RUBBLE_PILE *wrp,FILE *summary_fp)
{
	RUBBLE_PILE *tmp,*rp,*rpMax;
	OSC_ELEM *oe;
	VECTOR r,v,h;
	double mi,m1,m2,m,vd1,vd2,vd,mTotal;
	double mAcc,mOrb,mEsc,sm,sr,eb;
	int nTotal,n1,n2,n,i,j;

	{
	  if(p->bDoPatch) {
	    if(convert_to_patch_coords(p,wrp)) {
	      (void) printf("Error while converting from space to patch frame!\n");
	      return;
	    }
	  }
	}

	if (p->bVerbose)
		(void) printf("Finding groups...\n");

	/* Analyze world file */

	rpuAnalyze(wrp);

	/* Find groupings using temporary storage */

	tmp = (RUBBLE_PILE *) malloc((size_t) wrp->n_particles*sizeof(RUBBLE_PILE));
	assert(tmp != NULL);
	find_rp(p,wrp,tmp);

	/* Shrink to minimum fitting array */

	for (nTotal=i=0;i<wrp->n_particles;i++)
		if (tmp[i].data != NULL)
			++nTotal;

	rp = (RUBBLE_PILE *) malloc((size_t) nTotal*sizeof(RUBBLE_PILE));
	assert(rp != NULL);

	for (i=j=0;i<wrp->n_particles;i++)
		if (tmp[i].data != NULL)
			rp[j++] = tmp[i]; /* struct copy */

	assert(j == nTotal);

	free((void *) tmp); /* (tmp[i].data released later when rp freed) */

	/* Get grouping properties */

	if (p->bVerbose)
		(void) printf("Analyzing...\n");

	for (i=0;i<nTotal;i++) { /* only need some basic quantities */
		rpuAnalyze(&rp[i]);
#ifdef OLD /*DEBUG previously didn't need to calc axes, so that's why we didn't bother with rpuAnalyze() here, but new version of spin calc needs axes (to compute Dan's spin params), which can lead to a seg fault here*/
		rpuCalcMass(&rp[i]);
		rpuCalcPos(&rp[i]);
		rpuCalcVel(&rp[i]);
		rpuCalcRadius(&rp[i]);
		rpuCalcInertia(&rp[i]);
		rpuCalcSpin(&rp[i]);
		rpuCalcAggID(&rp[i]);
		/*rpuCalcColor(&rp[i]);*//*DEBUG not currently needed*/
#endif
		if(p->bAlternateR)
		  /* Overwrite radius with volume-equivalent sphere */
		  rp[i].radius = cbrt(rp[i].axis_len[0]*
				      rp[i].axis_len[1]*
				      rp[i].axis_len[2]);
	}
	
	/* Sort in descending order of mass */

	qsort((void *) rp,(size_t) nTotal,sizeof(RUBBLE_PILE),mass_sort);

	/* Generate grouping stats */

	n1 = n2 = n = 0;
	m1 = m2 = m = vd1 = vd2 = vd = 0.0;
	for (i=0;i<nTotal;i++) {
		mi = rp[i].mass;
		SUB_VEC(rp[i].vel,wrp->vel,v);
		if (rp[i].n_particles == 1) {
			++n1;
			m1 += mi;
			vd1 += mi*MAG_SQ(v);
			}
		else if (rp[i].n_particles == 2) {
			++n2;
			m2 += mi;
			vd2 += mi*MAG_SQ(v);
			}
		else {
			assert(rp[i].n_particles > 2);
			++n;
			m += mi;
			vd += mi*MAG_SQ(v);
			}
		}
	assert(n1 + n2 + n == nTotal);
	mTotal = m1 + m2 + m;
	if (m1 > 0) vd1 = sqrt(vd1/m1);
	if (m2 > 0) vd2 = sqrt(vd2/m2);
	if (m > 0) vd = sqrt(vd/m);

	rpuAnalyze(rpMax = &rp[0]); /* rpMax points to first element in rp */

	/* Get osculating element stats wrt largest clump */

	oe = (OSC_ELEM *) malloc((size_t) (nTotal - 1)*sizeof(OSC_ELEM));
	assert(oe != NULL);
	mAcc = mOrb = mEsc = 0;
	for (i=1;i<nTotal;i++) {
		j = i - 1;
		mi = rp[i].mass;
		sm = mi + rpMax->mass;
		sr = rp[i].radius + rpMax->radius;
		SUB_VEC(rp[i].pos,rpMax->pos,r);
		SUB_VEC(rp[i].vel,rpMax->vel,v);
		eb = 0.5*MAG_SQ(v) - sm/MAG(r);
		CROSS(r,v,h);
		oe[j].a = -0.5*sm/eb;
		oe[j].e = sqrt(1 - MAG_SQ(h)/(oe[j].a*sm));
		oe[j].q = (1 - oe[j].e)*oe[j].a;
		if (eb >= 0) mEsc += mi;
		else if (oe[j].q <= sr) mAcc += mi;
		else mOrb += mi;
		}

	/* Write line to summary file */

	if (p->bVerbose)
		(void) printf("Writing...\n");

	if (summary_fp) {
		(void) fprintf(summary_fp,
					   "%.16e %.16e "
					   "%.16e %.16e %.16e "
					   "%.16e %.16e "
					   "%.16e %.16e %.16e "
					   "%.16e %.16e "
					   "%.16e %.16e "
					   "%.16e %.16e %.16e "
					   "%.16e "
					   "%i %i %i %.16e %.16e %.16e %.16e %.16e %.16e "
					   "%.16e %.16e %.16e\n",
					   p->dTime,wrp->radius,
					   MAG(wrp->pos),MAG(wrp->vel),MAG(wrp->spin),
					   rpMax->mass,MAG(rpMax->vel),
					   rpMax->spin[X],rpMax->spin[Y],rpMax->spin[Z],
					   rpMax->eff_spin,rpMax->rot_idx,
					   rpMax->ang_mom_mag,rpMax->kin_energy,
					   rpMax->axis_len[MAJOR(rpMax)],
					   rpMax->axis_len[INTER(rpMax)],
					   rpMax->axis_len[MINOR(rpMax)],
					   rpMax->density,
					   n1,n2,n,m1,m2,m,vd1,vd2,vd,
					   mAcc,mOrb,mEsc);		
		(void) fflush(summary_fp);
		}

	/* Output grouping stats if last file and end stats desired */

	if (p->bDoEndStats && p->bLastFile) {
		COLOR c;
		double world_m_color[NUM_COLORS],max_m_color[NUM_COLORS];
		int world_n_color[NUM_COLORS],max_n_color[NUM_COLORS];

		(void) printf("\nFinal grouping summary:\n\n");
		(void) printf("                  Number Number%% Mass%%\n");
		(void) printf("Rubble piles      %6i %7.1f %5.1f\n",
					  n,100.0*n/nTotal,100*m/mTotal);
		(void) printf("2-particle groups %6i %7.1f %5.1f\n",
					  n2,100.0*n2/nTotal,100*m2/mTotal);
		(void) printf("Free particles    %6i %7.1f %5.1f\n",
					  n1,100.0*n1/nTotal,100*m1/mTotal);
		(void) printf("                  --------------------\n");
		(void) printf("Total             %6i   100.0 100.0\n",nTotal);

		/* Get world color distribution */

		for (i=0;i<NUM_COLORS;i++)
			world_m_color[i] = world_n_color[i] = 0;
		for (i=0;i<wrp->n_particles;i++) {
			c = (COLOR) wrp->data[i].color;
			world_m_color[c] += wrp->data[i].mass;
			++world_n_color[c];
			}

		/* Output world stats */

		(void) printf("\nWorld color distribution summary:\n\n");
		(void) printf("Color Number Number%% Mass%%\n");
		for (i=0;i<NUM_COLORS;i++)
			if (world_n_color[i])
				(void) printf("%5i %6i %7.1f %5.1f\n",i,world_n_color[i],
							  100.0*world_n_color[i]/wrp->n_particles,
							  100*world_m_color[i]/wrp->mass);
		(void) printf("      --------------------\n");
		(void) printf("Total %6i   100.0 100.0\n",wrp->n_particles);

		/* Generate stats for largest group */

		for (i=0;i<NUM_COLORS;i++)
			max_m_color[i] = max_n_color[i] = 0;
		for (i=0;i<rpMax->n_particles;i++) {
			c = (COLOR) rpMax->data[i].color;
			max_m_color[c] += rpMax->data[i].mass;
			++max_n_color[c];
			}

		/* Output largest group stats */

		(void) printf("\nSummary for most massive group:\n\n");
		(void) printf("Color Number Number%% Mass%%\n");
		for (i=0;i<NUM_COLORS;i++)
			if (max_n_color[i])
				(void) printf("%5i %6i %7.1f %5.1f\n",i,max_n_color[i],
							  100.0*max_n_color[i]/rpMax->n_particles,
							  100*max_m_color[i]/rpMax->mass);
		(void) printf("      --------------------\n");
		(void) printf("Total %6i   100.0 100.0\n",rpMax->n_particles);

		/* Compute mixing stats for largest group */

		if (p->bDoMixingStats && NUM_MIXING_CALCS > 0) {
			double mixing[NUM_MIXING_CALCS],mix_avg = 0;

			for (i=0;i<NUM_MIXING_CALCS;i++)
				if (get_mixing(world_m_color,wrp->mass,rpMax,&mixing[i])) {
					(void) printf("\nRubble pile too small to determine mixing.\n");
					i = 0;
					break;
					}
				else mix_avg += mixing[i];

			if (i) {
				double mix_err = 0;
				assert(i == NUM_MIXING_CALCS);
				mix_avg /= NUM_MIXING_CALCS;
				if (NUM_MIXING_CALCS > 1) {
					for (i=0;i<NUM_MIXING_CALCS;i++)
						mix_err += SQ(mixing[i] - mix_avg);
					mix_err = sqrt(mix_err/(NUM_MIXING_CALCS - 1));
					}
				(void) printf("\nCOLOR MIXING IN LARGEST RUBBLE PILE = "
							  "%.1f +/- %.2f%%\n",mix_avg*100,mix_err*100);
				}

			(void) printf("\n");
			}

		/* Re-color by groupings (cycling outwards from c.o.m. if COLOR_BY_DIST defined) */

		if (!p->bDoNotRecolor)
			recolor(p,wrp,rp,nTotal);

		/* Undo patch manipulations, if needed */
		if(p->bDoPatch)
		  convert_to_space_coords(p,wrp);

		/* Save to groups files */

		{
			FILE *fp = fopen(CENTERS_BT,"w");
			int nn;

			if (p->bVerbose)
				(void) printf("\nWriting to \"%s\"...\n",CENTERS_BT);
			assert(fp != NULL);
			for (nn=i=0;i<nTotal;i++) {
				int iOrgIdx;
				if (!p->bSaveSingles && rp[i].n_particles == 1)
					continue;
				if (!p->bSaveDoubles && rp[i].n_particles == 2)
					continue;
				iOrgIdx = (p->bSaveGroupsAsAggs ? -1 - i :
						   (rp[i].agg_id >= 0 ? -1 - rp[i].agg_id : i));
				(void) fprintf(fp,"%i %i %.16e %.16e %.16e %.16e %.16e "
							   "%.16e %.16e %.16e %.16e %.16e %.16e %i\n",
							   i,iOrgIdx,rp[i].mass,rp[i].radius,
							   rp[i].pos[X],rp[i].pos[Y],rp[i].pos[Z],
							   rp[i].vel[X],rp[i].vel[Y],rp[i].vel[Z],
							   rp[i].spin[X],rp[i].spin[Y],rp[i].spin[Z],
							   rp[i].color);
				++nn;
				}
			(void) fclose(fp);
			if (p->bVerbose)
				(void) printf("%i entr%s written.\n",nn,nn==1?"y":"ies");
			}

		write_ss(p,GROUPS_SS,rp,0,nTotal);

		/* Write largest rubble pile(s) separately as well */
		/*DEBUG NOTE: we do not check if these are true rubble piles (N>2)...*/
		
		if (p->nBiggest > nTotal) p->nBiggest = nTotal;

		if (p->bSaveBiggestInOneFile || p->nBiggest == 1)
			write_ss(p,BIGGEST_IN_ONE_FILE_SS,rp,0,p->nBiggest);
		else {
			char filename[MAXPATHLEN];
			int rv;

			for (i=0;i<p->nBiggest;i++) {
				rv = snprintf(filename,MAXPATHLEN,BIGGEST_ONE_PER_FILE_SS,i);
				assert(rv > 0);
				write_ss(p,filename,rp,i,1);
				}
			}

		/* Write distribution data */

		/*
		** FORMAT: INDEX MASS RADIUS SPEED SPIN ROT_IDX COLOR
		**
		** INDEX   = mass order position (most massive = 0)
		** MASS    = mass in kilograms
		** RADIUS  = effective radius in metres
		** SPEED   = speed wrt com in metres/second
		** SPIN    = spin period in hours
		** ROT_IDX = rotation index
		** COLOR   = rubble pile color
		*/

		{
			FILE *fp = fopen(CENTERS_DD,"w");
			VECTOR vv;
			double w;
			int nn;

			if (p->bVerbose)
				(void) printf("\nWriting to \"%s\"...\n",CENTERS_DD);
			assert(fp != NULL);
			for (nn=i=0;i<nTotal;i++) {
				if (!p->bSaveSingles && rp[i].n_particles == 1)
					continue;
				if (!p->bSaveDoubles && rp[i].n_particles == 2)
					continue;
				rpuAnalyze(&rp[i]);
				SUB_VEC(rp[i].vel,wrp->vel,vv);
				w = MAG(rp[i].spin);
				(void) fprintf(fp,"%i %g %g %g %g %g %i\n",i,
							   rp[i].mass*M_SCALE,
							   pow(rp[i].axis_len[0]*rp[i].axis_len[1]*rp[i].axis_len[2],1.0/3)*L_SCALE,
							   MAG(vv)*V_SCALE,w==0.0?0.0:TWO_PI/(w/T_SCALE)/3600.0,
							   rp[i].rot_idx,rp[i].color);
				++nn;
				}
			(void) fclose(fp);
			if (p->bVerbose)
				(void) printf("%i entr%s written.\n",nn,nn==1?"y":"ies");
			}

		/* And now, osculating elements wrt largest rubble pile */

		/*
		** FORMAT: INDEX MASS SMA ECC PERI COLOR
		**
		** INDEX = mass order position (most massive = 0)
		** MASS  = mass in units of largest rubble pile mass
		** SMA   = semimajor axis in units of largest rubble pile maximum radius
		** ECC   = eccentricity
		** PERI  = periapsis in units of sum of maximum radii
		** COLOR = rubble pile color
		*/

		{
			FILE *fp = fopen(CENTERS_OE,"w");
			int nn;

			if (p->bVerbose)
				(void) printf("\nWriting to \"%s\"...\n",CENTERS_OE);
			assert(fp != NULL);
			for (nn=i=0;i<nTotal - 1;i++) {
				if (!p->bSaveSingles && rp[i+1].n_particles == 1)
					continue;
				if (!p->bSaveDoubles && rp[i+1].n_particles == 2)
					continue;
				(void) fprintf(fp,"%i %g %g %g %g %i\n",i,
							   rp[i+1].mass/rpMax->mass,
							   oe[i].a/rpMax->radius,oe[i].e,
							   oe[i].q/(rpMax->radius + rp[i+1].radius),
							   rp[i+1].color);
				++nn;
				}
			(void) fclose(fp);
			if (p->bVerbose)
				(void) printf("%i entr%s written.\n",nn,nn==1?"y":"ies");
			}
		} /* if end stats & last file */

	/* Free resources */

	free((void *) oe);
	for (i=0;i<nTotal;i++) {
		assert(rp[i].data != NULL);
		rpuFree(&rp[i]);
		}
	free((void *) rp);

	if (p->bVerbose)
		(void) printf("\nDone!\n");
	}

static int read_data(PARAMS *p,const char *filename,RUBBLE_PILE *wrp)
{
	SSIO ssio;
	SSHEAD h;
	int i,n;

	(void) printf("%s: ",filename);

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\" for reading\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"Corrupt header\n");
		(void) ssioClose(&ssio);
		return 1;
		}

	if (h.n_data <= 0) {
		(void) fprintf(stderr,"Invalid file size\n");
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

	p->dTime = h.time;
	(void) printf("time %g yr...",p->dTime*T_SCALE/SID_YR);
	if(p->bDoPatch)
	  if(p->dTime == 0.0)
	    (void) printf("WARNING: Time in header is zero.  Potential problem for spin reorientation!\n");
	
	n = wrp->n_particles = h.n_data;
	rpuMalloc(wrp);

	for (i=0;i<n;i++) {
		if (ssioData(&ssio,&wrp->data[i])) {
			(void) fprintf(stderr,"Corrupt data\n");
			rpuFree(wrp);
			(void) ssioClose(&ssio);
			return 1;
			}
		else if (!p->bUseReservedColors && RESERVED_COLOR(wrp->data[i].color)) {
			--n;
			--i;
			}
	}
	
	wrp->n_particles = n;

	rpuRealloc(wrp);

	(void) printf("read %i particle%s (%i skipped)\n",wrp->n_particles,
				  wrp->n_particles == 1 ? "" : "s",h.n_data - wrp->n_particles);

	(void) ssioClose(&ssio);

	return 0;
	}

static void write_params_header(const PARAMS *p,FILE *fp)
{
	/* write headers for rpa.py parsing... */
	fprintf(fp,"# rpa.out summary file\n");
	fprintf(fp,"# main options:");
	fprintf(fp," -l %g",p->dLinkingScale);
	if (!p->bSaveSingles) fprintf(fp," -1");
	if (!p->bSaveDoubles) fprintf(fp," -2");
	if (p->bMergeAggs) fprintf(fp," -A");
	if (p->bLinkByParticle) fprintf(fp," -P");
	if (p->bDoPatch) fprintf(fp," -p");
	if (p->bAlternateR) fprintf(fp, " -r");
	fprintf(fp," -t %g",sqrt(p->dThetaSq));
	fprintf(fp,"\n");
	fprintf(fp,"# time world_radius world_pos world_vel world_spin max_mass max_vel max_spinx max_spiny max_spinz eff_spin rot_idx ang_mom ke max_major max_inter max_minor max_den n1 n2 n m1 m2 m vd1 vd2 vd mAcc mOrb mEsc\n");
	}

static void usage(const char *progname)
{
	fprintf(stderr, "Usage: %s [-e [-b #] [-g] [-m] [-o]] [-A|-B|-P] [-p] [-1|-2] [-a|-f|-n] [-l #] [-t #] [-d|-x] [-E] [-v] ss-file [ss-file ...]\n", progname);
	fprintf(stderr, "where -e = generate end stats and files\n");
	fprintf(stderr, "      -b = number of biggest groups to save\n");
	fprintf(stderr, "      -g = mark biggest as aggregates\n");
	fprintf(stderr, "      -m = include mixing stats\n");
	fprintf(stderr, "      -o = save biggest in one file instead of individual files\n");
	fprintf(stderr, "      -A = automatically link particles in aggregates\n");
	fprintf(stderr, "      -B = do not group find beyond aggregates (assumes -A and sets -l 0)\n");
	fprintf(stderr, "      -P = link by particles (not by groups)\n");
	fprintf(stderr, "      -p = assume local coordinate system (patch) geometry\n");
	fprintf(stderr, "      -1 = do not save single particles\n");
	fprintf(stderr, "      -2 = do not save singles or doubles\n");
	fprintf(stderr, "      -a = append to summary file\n");
	fprintf(stderr, "      -f = force overwrite of summary file\n");
	fprintf(stderr, "      -n = no output to summary file (use with -e)\n");
	fprintf(stderr, "      -l = linking scale (>= 1, or 0)\n");
	fprintf(stderr, "      -t = opening angle (>= 0)\n");
	fprintf(stderr, "      -r = use alternate rp radius (volume-equivalent sphere)\n");
	fprintf(stderr, "      -d = do not recolor\n");
	fprintf(stderr, "      -x = exclude reserved colors\n");
	fprintf(stderr, "      -E = return value equals number of unreadable input files (max 255)\n");
	fprintf(stderr, "      -v = verbose output\n");
	fprintf(stderr, "NOTE: -e assumed if only one ss-file provided\n");
	fprintf(stderr, "      -o assumed if number of biggest to save = 1\n");
	fprintf(stderr, "      -l 0 means do not group find\n");
	fprintf(stderr, "Defaults: -b %i -l %g -t %g\n", DFLT_NUM_BIGGEST, DFLT_LINKING_SCALE, DFLT_THETA);
	fprintf(stderr, "Summary file: %s\n", SUMMARY);
	exit(1);
}

int main(int argc,char *argv[])
{
	enum {None,Append,Force,NoSummary} mode=None;

	PARAMS params;
	RUBBLE_PILE wrp; /* "world" rubble pile */
	FILE *fp = NULL;
	double d;
	int c,n;
	int returnValue = 0; /* main's return value (see 'bError' parameter) */

	setbuf(stdout,(char *) NULL); /* no buffering */

#if TRAP_FPE/* || __STDC_VERSION__ >= 199901L*/
	/* enable explicit FPE trapping -- see #include above */
# if __APPLE__
	#error FPE handing not supported under macOS
# else
	feenableexcept(FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);
# endif
#else
	(void) fprintf(stderr,"WARNING: FPE trapping not enabled.\n");
#endif

	(void) srand((unsigned int) getpid()); /* seed random number generator */

	/* Defaults */

	params.dLinkingScale = DFLT_LINKING_SCALE;
	params.dThetaSq = SQ(DFLT_THETA);
	params.nBiggest = DFLT_NUM_BIGGEST;
	params.bDoEndStats = DFLT_DO_END_STATS;
	params.bSaveGroupsAsAggs = DFLT_SAVE_GROUPS_AS_AGGS;
	params.bDoMixingStats = DFLT_DO_MIXING_STATS;
	params.bSaveBiggestInOneFile = DFLT_SAVE_BIGGEST_IN_ONE_FILE;
	params.bLinkByParticle = DFLT_LINK_BY_PARTICLE;
	params.bMergeAggs = DFLT_MERGE_AGGS;
	params.bSaveSingles = DFLT_SAVE_SINGLES;
	params.bSaveDoubles = DFLT_SAVE_DOUBLES;
	params.bDoNotRecolor = DFLT_DO_NOT_RECOLOR;
	params.bUseReservedColors = DFLT_USE_RESERVED_COLORS;
	params.bVerbose = DFLT_VERBOSE;
	params.bDoPatch = DFLT_DOPATCH;
	params.bAlternateR = DFLT_ALTERNATER;
	params.bError = DFLT_ERROR;

	/* Parse command-line arguments */

	while ((c = getopt(argc,argv,"12ABEPab:defgl:mn:oprt:vxyz")) != EOF)
		switch (c) {
		case '1':
			params.bSaveSingles = FALSE;
			break;
		case '2':
			params.bSaveSingles = FALSE;
			params.bSaveDoubles = FALSE;
			break;
		case 'A':
			if (params.bLinkByParticle)
				usage(argv[0]);
			params.bMergeAggs = TRUE;
			break;
		case 'B':
			if (params.bLinkByParticle)
				usage(argv[0]);
			params.bMergeAggs = TRUE;
			params.dLinkingScale = 0.0;
			break;
		case 'E':
		        params.bError = TRUE;
		        break;
		case 'P':
			if (params.bMergeAggs)
				usage(argv[0]);
			params.bLinkByParticle = TRUE;
			break;
		case 'a':
			if (mode != None)
				usage(argv[0]);
			mode = Append;
			break;
		case 'b':
			n = atoi(optarg);
			if (n < 1)
				usage(argv[0]);
			params.nBiggest = n;
			break;
		case 'd':
			if (!params.bUseReservedColors)
				usage(argv[0]);
			params.bDoNotRecolor = TRUE;
			break;
		case 'e':
			params.bDoEndStats = TRUE;
			break;
		case 'f':
			if (mode != None) usage(argv[0]);
			mode = Force;
			break;
		case 'g':
			params.bSaveGroupsAsAggs = TRUE;
			break;
		case 'l':
			d = atof(optarg);
			if (d < 1.0 && d != 0.0)
				usage(argv[0]);
			params.dLinkingScale = d;
			break;
		case 'm':
			params.bDoMixingStats = TRUE;
			break;
		case 'n':
			if (mode != None)
				usage(argv[0]);
			mode = NoSummary;
			break;
		case 'o':
			params.bSaveBiggestInOneFile = TRUE;
			break;
		case 'p':
			params.bDoPatch = TRUE;
			break;
		case 'r':
			params.bAlternateR = TRUE;
			break;
		case 't':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			if (d > 1.0)
				(void) fprintf(stderr,"Warning: Large opening angle\n");
			params.dThetaSq = SQ(d);
			break;
		case 'v':
			params.bVerbose = TRUE;
			break;
		case 'x':
			if (params.bDoNotRecolor)
				usage(argv[0]);
			params.bUseReservedColors = FALSE;
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind >= argc)
		usage(argv[0]);

	if (optind == argc - 1)
		params.bDoEndStats = TRUE;

	switch (mode) {
	case None:
		if ((fp = fopen(SUMMARY,"r"))) {
			(void) fprintf(stderr,"%s already exists"
						   " -- use \"-a\", \"-f\", or \"-n\"\n",SUMMARY);
			(void) fclose(fp);
			return 1;
			}
		/* FALLTHRU */
	case Force:
		fp = fopen(SUMMARY,"w");
		if (!fp) {
			(void) fprintf(stderr,"Unable to open %s for writing\n",SUMMARY);
			return 0;
			}
		write_params_header(&params,fp);
		break;
	case Append:
		fp = fopen(SUMMARY,"a");
		if (!fp) {
			(void) fprintf(stderr,"Unable to open %s for appending\n",SUMMARY);
			return 0;
			}
		write_params_header(&params,fp);
		break;
	case NoSummary:
		if (!params.bDoEndStats) {
			(void) fprintf(stderr,"Nothing to do\n");
			return 1;
			}
		optind = argc - 1; /* skip to last file */
		break;
	default:
		assert(0);
		}

	for (;optind<argc;optind++) {
		if (read_data(&params,argv[optind],&wrp))
		  {
		    if(params.bError && returnValue < 255) returnValue++;
		    continue;
		  }
		params.bLastFile = (optind == argc - 1); /* last file? */
		analyze_world(&params,&wrp,fp);
		rpuFree(&wrp);
		}

	if (fp) (void) fclose(fp);

	return returnValue;
	}

/* rpa.c */
