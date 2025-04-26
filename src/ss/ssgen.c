/*
 ** ssgen.c -- DCR 02/03/12
 ** =======
 ** Generates initial conditions in a rectangular or ellipsoidal
 ** space, allowing for uncertainty in particle size, and using a tree
 ** code to eliminate particle overlaps.
 ** Update -- RLB 06/04/12
 ** now includes bi-modal and tri-modal radius and density distribution, and
 ** discontinuous and continuous power-law distribution 
 ** Update -- JVD 07/17/19
 ** Now includes gaussian-distributed particle velocities about a mean
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <string.h> /* for strlen() */
#include <ctype.h> /* for tolower() */
#include <time.h>
#include <math.h>
#include <assert.h>
#include <ss.h>
#include <random.h>

#define OUTFILENAME "ssgen.ss"
#define LOGFILENAME "ssgen.log"

#define MAX_NUM_PASS_DFLT 500

typedef int BOOLEAN;
//DeMartini - Edited here to add dVMean, dVSig
typedef struct {                
  int nPart,bEllipsoidal, bCylindrical, bReverseSort, nPassMax;
  double vAxes[N_DIM],dRadAvg,dBRad,dBDensity,dTRad,dTDensity,dRadDev,dCutoff,dDensity,dPIndex,dPRmax,dPFlag,dVMean,dVSig;
  } PARAMS;

/* position of particle, and particle radius? would need to ammend this for bimodality etc... -rb */
typedef struct {
  int keep;
  double pos[N_DIM],radius,density;
  } DATA;

double SQ(double x);
#define SQ(x) ((x)*(x))

double CUBE(double x);
#define CUBE(x) ((x)*SQ(x))

double DISTSQ(double v1[3],double v2[3]);
#define DISTSQ(v1,v2) ((SQ(v1[0] - v2[0]) + SQ(v1[1] - v2[1]) + SQ(v1[2] - v2[2])))

static BOOLEAN get_yn(const char *str,const char *dflt_str)
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

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  double cmin[N_DIM],cmin_eff[N_DIM],cmax[N_DIM],cmax_eff[N_DIM],pos[N_DIM]; /*NDIM=3 defined in ssio.h*/
  struct node_s *cell[CELLS_PER_NODE];
  DATA *leaf[CELLS_PER_NODE];
  };

typedef struct node_s NODE;

/*writing the logfile - rb*/
static void write_log(const PARAMS *p,const DATA *d)
{
	FILE *fp;

	fp = fopen(LOGFILENAME,"w");
	if (fp == NULL) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",LOGFILENAME);
		exit(1);
		}

	fprintf(fp,"Target particle radius = %g cm\n",100.0*p->dRadAvg*AU);
	/* add if statement here if bi-modality or tri-modality is specified. */
	if (p->dBRad > 1 && p->dTRad == 1.0){	
		fprintf(fp,"Bi-modal Distribution specified with:\n");
		fprintf(fp,"	Radius Ratio = %.2f, Density Ratio = %.2f \n",p->dBRad,p->dBDensity);}
	if (p->dBRad > 1 && p->dTRad > 1){	
		fprintf(fp,"Tri-modal Distribution specified with:\n");
		fprintf(fp,"	Radius Ratio = %.2f:%.2f:1, Density Ratio = %.2f:%.2f:1 \n",p->dTRad,p->dBRad,p->dTDensity,p->dBDensity);}
	if (p->dPFlag == 1){
		fprintf(fp,"Discontinuous Power-Law Distribution specified with:\n");
		fprintf(fp,"	Index = %.2f, Rmax/Rmin = %.2f \n",p->dPIndex,p->dPRmax);}
	if (p->dPFlag == 2){
		fprintf(fp,"Continuous Power-Law Distribution specified with:\n");
		fprintf(fp,"	Index = %.2f, Rmax/Rmin = %.2f \n",p->dPIndex,p->dPRmax);}
	//DeMartini - add if statement here if velocity is specified
	if (p->dVMean > 0 || p->dVSig > 0) {
		fprintf(fp,"Velocity Distribution (m/s) specified with:\n");
		fprintf(fp,"	Mean = %g , Sigma = %g \n",p->dVMean*V_SCALE,p->dVSig*V_SCALE);}
	fprintf(fp,"Index of power law distribution = %g \n",p->dPIndex);
	fprintf(fp,"Target particle radius uncertainty = %g cm\n",100.0*p->dRadDev*AU);
	fprintf(fp,"Gaussian distribution cutoff = %g sigma\n",p->dCutoff);
	fprintf(fp,"Particle density = %g g/cc\n",1.0e-3*p->dDensity*M_SUN/(AU*AU*AU));
	fprintf(fp,"Region dimensions = %g x %g x %g m ",p->vAxes[0]*AU,p->vAxes[1]*AU,p->vAxes[2]*AU);
  if (p->bEllipsoidal)
    fprintf(fp,"(ellipsoidal)\n");
  else if (p->bCylindrical)
    fprintf(fp,"(cylindrical)\n");
  else
    fprintf(fp,"(rectangular)\n");
  if (p->bReverseSort)
    fprintf(fp,"Sorted by decreasing radius\n");
	fprintf(fp,"Number of particles = %i\n",p->nPart);
	/*probably need to change this as this now reflects the average radius of the smallest particle*/
	if (p->nPart > 1) {
		double x,s;
		int i;

		x = 0.0;
		for (i=0;i<p->nPart;i++) /* determining the average particle radius by adding up all the radii and diving by the # of particles */
			x += d[i].radius;
		x /= p->nPart;
		s = 0.0;
		for (i=0;i<p->nPart;i++) /* determining the variance in particle radius */
			s += SQ(d[i].radius - x); /* add up all the differences between each particle radius and average */
		s = sqrt(s/(p->nPart - 1)); /* take sqrt of (value / (# of particles -1)) */
		fprintf(fp,"Average particle radius = %g +/- %g cm\n",100.0*x*AU,100.0*s*AU);
		}

	fprintf(fp,"Minimum particle radius = %g cm\n",100.0*d[0].radius*AU); /*particles indices are arranged by radius, determined by other function?*/
	fprintf(fp,"Maximum particle radius = %g cm\n",100.0*d[p->nPart - 1].radius*AU);
  if (p->bReverseSort ==1)
    fprintf(fp,"Reverse Sorted\n");
  fprintf(fp,"Suggested minimum nSmooth:\n");
  /* for HSDEM, ratio of largest particle surface area to smallest particle cross-section */
  /* for SSDEM, number of small particles in sphere of radius small + large particle radii */
  if (p->bReverseSort ==1){
	fprintf(fp,"   HSDEM = %i\n",4*((int) (SQ(d[0].radius/d[p->nPart - 1].radius) + 0.5)));	/*nsmooth  = 4*(sqrt(largest_radius/smallest_radius + 1/2)) - an integer*/
    fprintf(fp,"   SSDEM = %i\n",(int) (CUBE(d[0].radius/d[p->nPart - 1].radius + 1.0) + 0.5));
    }
  else{
	fprintf(fp,"   HSDEM = %i\n",4*((int) (SQ(d[p->nPart - 1].radius/d[0].radius) + 0.5)));	/*nsmooth  = 4*(sqrt(largest_radius/smallest_radius + 1/2)) - an integer*/
    fprintf(fp,"   SSDEM = %i\n",(int) (CUBE(d[p->nPart - 1].radius/d[0].radius + 1.0) + 0.5));
    }
  fprintf(fp,"   (note: a smaller value for SSDEM is possible depending on filling factor)\n");
	fclose(fp);
	}

//DeMartini - Needs to take a dMean arg - change here and in all calls
static double gaussian(double dMean,double dSig,double dLim)
{
	/* truncated Gaussian */

	double dVal;

	do {
		dVal = dSig*randGaussian() + dMean;	//DeMartini - Added dMean
		} while (dLim > 0.0 && fabs(dVal) - fabs(dMean) > dLim*dSig);

	return dVal;
	}

/*function that is responsible for writing to ssgen.ss*/
static void write_data(const PARAMS *p,const DATA *d)
{
	/* these data-structures are defined in ssio.h */
	SSIO ssio;	
	SSHEAD head;
	SSDATA data;
	int i;
	double theta, mag, prefac;
	/*checks that ssgen.ss is open for writing?*/
	if (ssioOpen(OUTFILENAME,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",OUTFILENAME);
		exit(1);
		}

	head.time = 0.0;
	head.n_data = p->nPart;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	/*makes sure header is properly written (ss standard header?) ?*/
	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Unable to write header.\n");
		ssioClose(&ssio);
		exit(1);
		}
	/*writing data... to ssdata data structure*/
	/*using top defined DATA structure to define SSDATA members (named data.), defined within this function*/
	for (i=0;i<p->nPart;i++) {
		data.mass = 4.0/3.0*M_PI*d[i].radius*d[i].radius*d[i].radius*d[i].density; /*mass given by particle radius and density (given by user)*/
		data.radius = d[i].radius;
		data.pos[0] = d[i].pos[0];
		data.pos[1] = d[i].pos[1];
		data.pos[2] = d[i].pos[2];
		//DeMartini - Calculate Velocities
		data.vel[2] = 2.0*randUniform() - 1.0;		//v_z normalized between -1,1
		theta = TWO_PI*randUniform();			//Angle between 0,2pi
		mag = gaussian(p->dVMean, p->dVSig, 3.0);	// v mag, from specified gaussian
		prefac = mag*sqrt(1-data.vel[2]*data.vel[2]);	//Ensure normal vector
		data.vel[0] = prefac*cos(theta);
		data.vel[1] = prefac*sin(theta);
		data.vel[2] = mag*data.vel[2];
		data.spin[0] = data.spin[1] = data.spin[2] = 0.0;
		data.color = PLANETESIMAL;
		data.org_idx = i;
		if (ssioData(&ssio,&data)) {
			fprintf(stderr,"Error writing data (particle %i).\n",i);
			ssioClose(&ssio);
			exit(1);
			}
		}

	ssioClose(&ssio);
	}
/*node's used for rejection algorithm*/
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

static void get_pos(const PARAMS *p,DATA *d)
{
	/*
	** NOTE: for the ellipsoidal or cylindrical region, this routine
	** guarantees that the particle center -- but not necessarily the
	** entire particle -- is contained within the ellipsoid (both are
	** guaranteed for the rectangular region).
	*/

    int k;
    double a,max_radius;
    
    if (p->dBRad > 1.0){
        if (p->dTRad > 1.0)
            max_radius = p->dRadAvg*p->dTRad;
        else
            max_radius = p->dRadAvg*p->dBRad;
    }
    if (p->dPRmax > 1.0){
        max_radius = p->dRadAvg*p->dPRmax;
    }

    if (!p->bCylindrical){
    do {
      for (k=0;k<N_DIM;k++) {
        a = p->vAxes[k] - 2.0*max_radius;
        if (a < 0.0) {
          fprintf(stderr,"get_pos(): Region too small for particle.\n");
          exit(1);
          }
        d->pos[k] = (randUniform() - 0.5)*a;
        }
      } while (p->bEllipsoidal &&
           SQ(d->pos[0])/SQ(p->vAxes[0]) +
           SQ(d->pos[1])/SQ(p->vAxes[1]) +
           SQ(d->pos[2])/SQ(p->vAxes[2]) > 0.25);
    }
    else if (p->bCylindrical) {
    do {
      for (k=0;k<N_DIM;k++) {
        a = p->vAxes[k] - 2.0*max_radius;
        if (a < 0.0) {
	  printf("a=%g x=%g y=%g z=%g m=%g\n",a,p->vAxes[0],p->vAxes[1],p->vAxes[2],max_radius);
	  printf("R=%g max=%g\n",d->radius,p->dPRmax);
          fprintf(stderr,"get_pos(): Region too small for particle.\n");
          exit(1);
        }
        d->pos[k] = (randUniform() - 0.5)*a;
        }
      } while (SQ(d->pos[0])/SQ(p->vAxes[0]) +
               SQ(d->pos[1])/SQ(p->vAxes[1]) > 0.25);
    }
    }

/* generating particles */
static void generate(const PARAMS *p,DATA *d)
{
	NODE *root;
	BOOLEAN bTest;
	long i,j,k=0,res;
	int nPass,nRej=0/*for first pass*/;
	double Nb,Nm,Ns,unidev;
	randSeedGenerator((Ullong) (time(NULL) % getpid() + getppid()));

	/* pregenerate particle radii */
	res=(long)p->nPart; /* set resolution to total number of particles (i.e. each particle has a unique radius)
			      Seems to work best for smooth distribution. If this is optimal, then remove second for loop
			      Currently keeping this here if a resolution change is desired in future (possibly as user-input)
			      However, may be best to remove resolution parameter completely so as to remove possibility of resolution > N_particles*/
	/*Continuous and Discontinuous Power Law Distribution*/
	if (p->dPFlag == 2.0){
				
		if (p->dPIndex != -1){	
		for (i=0;i<res;i++){
			for (j=k;j<p->nPart*(i+1)/res;j++){
				d[j].keep=0;
			do{
				unidev=(float) (i+1) / (float) res;
				d[j].radius=pow(((1-unidev)*pow(p->dRadAvg,p->dPIndex +1))+(unidev*pow(p->dPRmax*p->dRadAvg,p->dPIndex+1)),1/(p->dPIndex+1));
				d[j].density=p->dDensity;
				} while (d[j].radius <= 0.0);			
			}
			k=j;
		}}

		if (p->dPIndex == -1){
		for (i=0;i<res;i++){
			for (j=k;j<p->nPart*(i+1)/res;j++){
				d[j].keep=0;
			do{
				unidev=(float) (i+1) / (float) res;
				d[j].radius=pow(p->dPRmax*p->dRadAvg,unidev)/pow(p->dRadAvg,unidev-1);
				d[j].density=p->dDensity;
				} while (d[j].radius <= 0.0);			
			}
			k=j;
		}}				
	}
	
	if (p->dPFlag == 1.0){
		if (p->dPIndex != -1){	
		for (i=0;i<p->nPart;i++){
			d[i].keep = 0;
			do{	
				unidev=randUniform();
				d[i].radius=pow(((1-unidev)*pow(p->dRadAvg,p->dPIndex +1))+(unidev*pow(p->dPRmax*p->dRadAvg,p->dPIndex+1)),1/(p->dPIndex+1));
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}}

		if (p->dPIndex == -1){
		for (i=0;i<p->nPart;i++){
			d[i].keep = 0;
			do{	
				unidev=randUniform();
				d[i].radius=pow(p->dPRmax*p->dRadAvg,unidev)/pow(p->dRadAvg,unidev-1);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}}

	}

	/*Whether or not bi-modality is slected this will generate particle radii based on equal Volume*/
	if (p->dTRad == 1 && p->dTDensity == 1){
    if (p->dBRad > 1.0 || p->dBDensity != 1){
      bTest = get_yn("Bimodality: Specify Number of Larger (Denser) Particles?","n");
      if (bTest){	
        do {
          printf("Enter Number of Larger Particles: ");
          (void) scanf("%lf",&Nm);
          if (Nm == 0.0) printf("Number of Larger (Denser) Particles cannot be zero.\n");
          } while (Nm == 0.0);
      }	
      else Nm=p->nPart/(1+pow(p->dBRad,3));
      

      Ns=p->nPart-Nm;
      for (i=0;i<Ns;i++) {
        d[i].keep = 0;
        do {
          d[i].radius = p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff);
          d[i].density=p->dDensity;
          } while (d[i].radius <= 0.0); /* reject unphysical values */
        }
      for (j=i;j<p->nPart;j++) {
        d[j].keep = 0;
        do {
          d[j].radius = (p->dBRad)*(p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff));
          d[j].density=(p->dBDensity)*p->dDensity;
                  } while (d[j].radius <= 0.0); /* reject unphysical values */
        }
      }
    }

	if (p->dTRad > 1.0 || p->dTDensity != 1){
    bTest = get_yn("Trimodality: Specify Number of Larger Particles?","n");
    if (bTest){
      do {
        printf("Enter Number of Largest Particles: ");
        (void) scanf("%lf",&Nb);
        if (Nb == 0.0) printf("Number of Largest Particles cannot be zero.\n");
        if (Nb > p->nPart) printf("Number of Largest Particles cannot exceed total number of particles.\n");
      } while (Nb == 0.0 || Nb > p->nPart);
      do {
        printf("Enter Number of Second Largest Particles: ");
        (void) scanf("%lf",&Nm);
        if (Nm == 0.0) printf("Number of Larger (Denser) Particles cannot be zero.\n");
        if (Nm > p->nPart) printf("Number of Second Largest Particles cannot exceed total number of particles.\n");
        if (Nm+Nb > p->nPart-1) printf("At least 1 particle has to be small, cannot have number of larger particles (1st and 2nd) be greater than number of particles - 1.\n");
      } while (Nm == 0.0 || Nb > p->nPart || Nm+Nb > p->nPart-1);
      Ns=p->nPart-Nm-Nb;
    }
    else {
      Nb=p->nPart/( 1+pow(p->dTRad,3)+pow(p->dTRad/p->dBRad,3));
      Nm=p->nPart/( 1+pow(p->dBRad,3)+pow(p->dBRad/p->dTRad,3));
      Ns=p->nPart-Nm-Nb;
      }
		for (i=0;i<Ns;i++) {
			d[i].keep = 0;
			do {
				d[i].radius = p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}
		for (j=i;j<Ns+Nm;j++) {
			d[j].keep = 0;
			do {
				d[j].radius = (p->dBRad)*(p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff));
				d[j].density=(p->dBDensity)*p->dDensity;
				} while (d[j].radius <= 0.0); /* reject unphysical values */
			}

		for (k=j;k<p->nPart;k++) {
			d[k].keep = 0;
			do {
				d[k].radius = (p->dTRad)*(p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff));
				d[k].density=(p->dTDensity)*p->dDensity;
				} while (d[k].radius <= 0.0); /* reject unphysical values */
			}
    }

	/*default arrangment*/
	if (p->dTRad == 1.0 && p->dBRad == 1.0 && p->dPFlag == 0.0 && p->dBDensity == 1.0 && p->dTDensity == 1.0 ){
		for (i=0;i<p->nPart;i++) {
			d[i].keep = 0;
			do {
				d[i].radius = p->dRadAvg + gaussian(0.0,p->dRadDev,p->dCutoff);
				d[i].density=p->dDensity;
				} while (d[i].radius <= 0.0); /* reject unphysical values */
			}
	}
		
	/*putting this here as a precaution for now*/
  for (i=0;i<p->nPart;i++) d[i].keep=0;
  
  /*sorting in size order speeds up the rejection algorithm*/
  qsort(d,(size_t) p->nPart,sizeof(DATA),compar);

  /* rejection algorithm */
	nPass = 0;
	do {
		if (++nPass > p->nPassMax) {
			fprintf(stderr,"Unable to converge after %i pass%s.\n",
					nPass - 1,nPass == 2 ? "" : "es");
			exit(1);
			}
		printf("Pass %i ",nPass);
		if (nPass == 1)
			printf("(nPart = %i)\n",p->nPart);
		else
			printf("(nRej = %i)\n",nRej);
		/* make root cell */
		make_node(-0.5*p->vAxes[0],0.5*p->vAxes[0],
				  -0.5*p->vAxes[1],0.5*p->vAxes[1],
				  -0.5*p->vAxes[2],0.5*p->vAxes[2],&root);
		/* generate particle positions (as needed) and add to tree */
		for (i=0;i<p->nPart;i++) {
			if (!d[i].keep) {
				get_pos(p,&d[i]);
				d[i].keep = 1;
				}
			add_to_tree(root,&d[i]);
			}
		/* check for overlaps, rejecting smallest particles first */
		nRej = 0;
		for (i=0;i<p->nPart;i++)
			if (reject(root,&d[i])) {
				d[i].keep = 0;
				++nRej;
				}
		kill_node(root);
		} while (nRej > 0);
  
  if (p->bReverseSort) /* sort particles by decreasing radius */
	qsort(d,(size_t) p->nPart,sizeof(DATA),rev_compar);
    
  }

static void usage(const char *achProgName)
{
	//DeMartini - Add usage lines for velocity option "-V"
	fprintf(stderr, "Usage: %s -r radius -d density [-b RadRatio,DenRatio] [-t R3:R2,D3:D2] [-p index,Rmax/Rmin,flag] [-s stddev] [-q cutoff] -x lx -y ly -z lz [-e] [-c] [-v] [-V mean, sigma] [-n max] N\n", achProgName);
	fprintf(stderr, "where: radius is the mean particle radius (cm) [min radius for -b, -t, or -p];\n");
	fprintf(stderr, "       density is the particle density (g/cc);\n");
	fprintf(stderr, "       RadRadio,DenRatio set radius & density of 2nd particle type\n");
	fprintf(stderr, "          for bi-modal distribution, e.g., -b 1.4,0.8 means second\n");
	fprintf(stderr, "          particle has 1.4x radius and 0.8x density of first;\n");
	fprintf(stderr, "       R3:R2,D3:D2 is for tri-modal case, e.g., -t 2.4:1.2,1.1:3.5:\n");
	fprintf(stderr, "          means second particle has 2.4x radius and 1.1x density of\n");
	fprintf(stderr, "          first while third has 1.2x radius and 3.5x density of first;\n");
	fprintf(stderr, "       index,Rmax/Rmin,flag are for differential radius distribution:\n");
	fprintf(stderr, "          power law index (-1 allowed), max to min radius ratio,\n");
	fprintf(stderr, "          flag for discontinuous (=1) or continuous (=2) distribution;\n");
	fprintf(stderr, "       stddev is the uncertainty in particle radius (cm);\n");
	fprintf(stderr, "       cutoff is max # sigma for dispersions (default 1; 0 means no limit);\n");
	fprintf(stderr, "       lx, ly, lz are the simulation region dimensions (m);\n");
	fprintf(stderr, "       -e means ellipsoidal, not rectangular, region;\n");
	fprintf(stderr, "       -c means cylindrical region (lx,ly are cross-section dimensions);\n");
	fprintf(stderr, "       -v means reverse sort order (decreasing radius order);\n");
	fprintf(stderr, "       mean and sigma are the positive mean and standard deviation of the\n");
	fprintf(stderr, "          gaussian-distributed, randomly-oriented particle velocities (m/s)\n");
	fprintf(stderr, "          (defaults to 0 for both - all particles will have 0 velocity);\n");
	fprintf(stderr, "       max is the maximum number of passes allowed (default %i);\n",MAX_NUM_PASS_DFLT);
	fprintf(stderr, "       and N is the number of particles to generate.\n");
	fprintf(stderr, "       NOTE: region dimensions are full axes, not semi axes;\n");
	fprintf(stderr, "             particle radius taken into account where possible.\n");

	exit(1);
}

int main(int argc,char *argv[])
{
	PARAMS params;
	DATA *data;
	double d,dd,ddd,dddd;
	char *rstring, *dstring, *r3, *r2, *d3, *d2, *istring, *fstring;	
	int c;
	
	assert(MAX_NUM_PASS_DFLT > 0);
	//DeMartini - Add lines here for new options
	/* defaults */

	params.dBRad = 1.0; /*new parameter need to define*/
	params.dBDensity = 1.0;
	params.dTRad = 1.0;
	params.dTDensity = 1.0;
	params.dPIndex = 0.0;
	params.dPRmax = 1.0;
	params.dPFlag = 0.0;

	params.dVMean = 0.0;
	params.dVSig = 0.0;

	params.dRadAvg = 0.0;
	params.dRadDev = 0.0;
	params.dCutoff = 1.0;
	params.dDensity = 0.0;
	params.vAxes[0] = params.vAxes[1] = params.vAxes[2] = 0.0;
	params.bEllipsoidal = 0;
	params.bCylindrical = 0;
	params.bReverseSort = 0;
	params.nPart = 0;

	params.nPassMax = MAX_NUM_PASS_DFLT;

	//DeMartini - Add new case ('v' already used) for velocities option -> '-V' (for now)
	/* parse command-line arguments */

	while ((c = getopt(argc,argv,"r:b:t:p:s:q:d:x:y:z:vecV:n:")) != EOF)
		switch (c) {
		case 'r':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dRadAvg = 0.01*d/AU; /* cm -> AU */
			break;
		case 'b':
			rstring = strtok (optarg, ","); /*parsing input string*/
			dstring = strtok (NULL, ",");
			d = atof(rstring);
			dd = atof(dstring);
			if (d < 1.0)
				usage(argv[0]);
			if (dd <= 0.0)
				usage(argv[0]);
			
			params.dBRad = d; 
			params.dBDensity =dd; 
			break;
		case 't':
			r3 = strtok (optarg, ":");
			r2 = strtok (NULL, ",");
			d3 = strtok (NULL, ":");
			d2= strtok (NULL, ":");
			d = atof(r2);
			dd = atof(d2);
			ddd = atof(r3);
			dddd = atof (d3);
			if (d <= 1.0)
				usage(argv[0]);
			if (dd <= 0.0)
				usage(argv[0]);
			if (ddd <= 1.0)
				usage(argv[0]);
			if (dddd <= 0.0)
				usage(argv[0]);
			
			params.dBRad = d; 
			params.dBDensity = dd;
			params.dTRad = ddd;
			params.dTDensity = dddd; 
			break;
		case 'p':
			istring = strtok (optarg, ",");
			rstring = strtok (NULL, ",");
			fstring = strtok (NULL, ",");
			d = atof(istring);
			dd = atof(rstring);
			ddd = atof(fstring); 			
			if (dd <= 1.0)
				usage(argv[0]);
			if (ddd != 1 && ddd != 2)
				usage(argv[0]);
			params.dPIndex = d; 
			params.dPRmax = dd;
			params.dPFlag = ddd; 
			break;				
		case 's':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dRadDev = 0.01*d/AU; /* cm -> AU */
			break;
		case 'q':
			d = atof(optarg);
			if (d < 0.0)
				usage(argv[0]);
			params.dCutoff = d;
			break;
		case 'd':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.dDensity = 1.0e3*d/M_SUN*AU*AU*AU; /* g/cc -> M_SUN/AU^3 */
			break;
		case 'x':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[0] = d/AU; /* m -> AU */
			break;
		case 'y':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[1] = d/AU; /* m -> AU */
			break;


		case 'z':
			d = atof(optarg);
			if (d <= 0.0)
				usage(argv[0]);
			params.vAxes[2] = d/AU; /* m -> AU */
			break;
		case 'v':
			params.bReverseSort = 1;
			break;
		case 'e':
			params.bEllipsoidal = 1;
			break;
		case 'c':
			params.bCylindrical = 1;
			break;
		case 'V':
			rstring = strtok (optarg, ","); //Easier to use existing variables r/dstring
			dstring = strtok (NULL, ",");
			d = atof(rstring);
			dd = atof(dstring);
			if (d < 0 || dd < 0)		//Randomly oriented, so only take positives
				usage(argv[0]);
			params.dVMean = d/V_SCALE;
			params.dVSig = dd/V_SCALE;
			break;
		case 'n':
			params.nPassMax = atoi(optarg);
			if (params.nPassMax < 0) /* = 0 is kind of silly too */
				usage(argv[0]);
			break;
		case '?':
		default:
			usage(argv[0]);
			}

	if (optind >= argc || params.dRadAvg == 0.0 || params.dDensity == 0.0 || params.vAxes[0] == 0.0 || params.vAxes[1] == 0.0 || params.vAxes[2] == 0.0)
		usage(argv[0]);

	if (params.dRadDev >= params.dRadAvg)
		fprintf(stderr,"WARNING: large particle radius uncertainty---infinite loop possible.\n");

	params.nPart = atoi(argv[optind]);

	if (params.nPart <= 0 || argc > optind + 1)
		usage(argv[0]);

	data = (DATA *) malloc((size_t) params.nPart*sizeof(DATA));
	assert(data != NULL);

	generate(&params,data);

	write_log(&params,data);
	printf("Log written to %s.\n",LOGFILENAME);

	write_data(&params,data);

	free(data);

	printf("Now run rpx on %s to reset center of mass, etc.\n",OUTFILENAME);

	return 0;
	}

/* ssgen.c */
