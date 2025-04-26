/*
 ** sscolor.c -- DCR 07/28/15
 ** =========
 ** Utility for changing particle colors in ss files (various options).
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <math.h>
#include <assert.h>
#include <time.h> /* for time() */
#include <boolean.h>
#include <vector.h>
#include <ss.h>
#include <random.h>

#define MAX_STR_LEN 256

typedef struct {
  /* user-supplied parameters, in main() order */
  BOOLEAN bColorAggs;
  int iColorMax;
  int iColorMin;
  int iColor;
  double dCycles;
  int iField;
  int iLogOpt;
  BOOLEAN bUseMKS;
  BOOLEAN bRandom;
  BOOLEAN bSortEach;
  BOOLEAN bUnique;
  double dVMax;
  double dVMin;
  int iScheme;
  } PARAMS_T;

typedef struct { /* for sorting */
  int iIndex;
  double dData;
  } SORT_DATA_T;

enum { /* color schemes */
  SchemeSS = 0,
  SchemeHSV,
  SchemeRGB
};

#define NUM_FIELDS 21

const char achFieldTags[NUM_FIELDS][MAX_STR_LEN] = { /* must match order in enum */
  "none",
  "oi",   /* original index */
  "m",    /* mass */
  "r",    /* radius */
  "x",    /* position */
  "y",
  "z",
  "vx",   /* velocity */
  "vy",
  "vz",
  "wx",   /* spin */
  "wy",
  "wz",
  "c",    /* color */
  "rmag", /* distance from origin */
  "vmag", /* absolute speed */
  "wmag", /* spin magnitude */
  "mden", /* mass density */
  "sden", /* spatial density */
  "ovlp", /* overlap fraction */
  "gpot"  /* gravitational potential */
};

enum { /* see achFieldTags above */
  FieldNone = 0,
  FieldOrgIdx,
  FieldMass,
  FieldRadius,
  FieldPosX,
  FieldPosY,
  FieldPosZ,
  FieldVelX,
  FieldVelY,
  FieldVelZ,
  FieldSpinX,
  FieldSpinY,
  FieldSpinZ,
  FieldColor,
  FieldRMag,
  FieldVMag,
  FieldWMag,
  FieldMDen,
  FieldSDen,
  FieldOvlp,
  FieldGPot
};

#define NUM_LOG_OPTS 3

const char achLogOptTags[NUM_LOG_OPTS][MAX_STR_LEN] = { /* must match order in enum */
  "none",
  "std", /* standard (i.e., no log option supplied) */
  "abs",
};

enum { /* see achLogOptTags above */
  LogNone = 0,
  LogStd,
  LogAbs
};

/* tree code definitions */

#define CELLS_PER_NODE 8 /* Barnes & Hut oct-tree */

struct node_s {
  VECTOR cmin, cmax, pos;
  struct node_s *cell[CELLS_PER_NODE];
  const SSDATA *leaf[CELLS_PER_NODE];
  double max_radius;
  };

typedef struct node_s NODE;

/* functions follow */

static double get_max_overlap(const SSDATA *d, const NODE *node)
{
	NODE *n;
	double r2, max_overlap, overlap;
	int i, k;

	/* intersect test -- check each subcell in turn */

	/*
	** Algorithm adapted from...
	** "A Simple Method for Box-Sphere Intersection Testing",
	** by Jim Arvo, in "Graphics Gems", Academic Press, 1990.
	** (http://www.ics.uci.edu/~arvo/code/BoxSphereIntersect.c)
	*/

	max_overlap = 0.0;
	
	for (i = 0; i < CELLS_PER_NODE; i++) {
		n = node->cell[i];
		if (n != NULL) {
			r2 = 0.0;
			for (k = 0; k < N_DIM; k++)
				if (d->pos[k] < n->cmin[k])
					r2 += SQ(d->pos[k] - n->cmin[k]);
				else if (d->pos[k] > n->cmax[k])
					r2 += SQ(d->pos[k] - n->cmax[k]);
			if (r2 <= SQ(d->radius + n->max_radius)) {
				overlap = get_max_overlap(d, n);
				if (overlap > max_overlap)
					max_overlap = overlap;
				}
			}
		else if (node->leaf[i] != NULL && node->leaf[i] != d &&
				 (r2 = DIST_SQ(d->pos, node->leaf[i]->pos)) <= SQ(d->radius + node->leaf[i]->radius)) {
			overlap = 1.0 - sqrt(r2) / (d->radius + node->leaf[i]->radius);
			if (overlap > max_overlap)
				max_overlap = overlap;
			}
		}
	
	return max_overlap;
	}

static void kill_node(NODE *node)
{
	int i;

	assert(node != NULL);

	for (i = 0; i < CELLS_PER_NODE; i++)
		if (node->cell[i])
			kill_node(node->cell[i]);

	free((void *) node);
	}

static void make_node(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, NODE **node)
{
	int i;

	*node = (NODE *) malloc(sizeof(NODE));
	assert(*node != NULL);

	(*node)->cmin[0] = xmin;
	(*node)->cmax[0] = xmax;
	(*node)->cmin[1] = ymin;
	(*node)->cmax[1] = ymax;
	(*node)->cmin[2] = zmin;
	(*node)->cmax[2] = zmax;

	(*node)->pos[0] = 0.5*((*node)->cmin[0] + (*node)->cmax[0]);
	(*node)->pos[1] = 0.5*((*node)->cmin[1] + (*node)->cmax[1]);
	(*node)->pos[2] = 0.5*((*node)->cmin[2] + (*node)->cmax[2]);

	for (i = 0; i < CELLS_PER_NODE; i++) {
		(*node)->cell[i] = NULL;
		(*node)->leaf[i] = NULL;
		}

	(*node)->max_radius = 0.0;
	}

static void add_to_tree(NODE *node, const SSDATA *d)
{
	int i, idx, idy, idz;

	idx = (d->pos[0] < node->pos[0] ? -1 : 1);
	idy = (d->pos[1] < node->pos[1] ? -1 : 1);
	idz = (d->pos[2] < node->pos[2] ? -1 : 1);

	i = (idx + 1) / 2 + (idy + 1 + 2 * (idz + 1));

	if (node->cell[i] != NULL)
		add_to_tree(node->cell[i], d);
	else if (node->leaf[i] != NULL) {
		make_node(idx < 0 ? node->cmin[0] : node->pos[0],
				  idx < 0 ? node->pos[0] : node->cmax[0],
				  idy < 0 ? node->cmin[1] : node->pos[1],
				  idy < 0 ? node->pos[1] : node->cmax[1],
				  idz < 0 ? node->cmin[2] : node->pos[2],
				  idz < 0 ? node->pos[2] : node->cmax[2],
				  &node->cell[i]);
		add_to_tree(node->cell[i], node->leaf[i]);
		add_to_tree(node->cell[i], d);
		node->leaf[i] = NULL;
		}
	else
		node->leaf[i] = d;

	if (d->radius > node->max_radius)
		node->max_radius = d->radius;
	}

static void build_tree(int n, const SSDATA *d, NODE **root)
{
	int i;
	double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, dmin, dmax;

	xmin = ymin = zmin = HUGE_VAL;
	xmax = ymax = zmax = -HUGE_VAL;
	
	/* get bounds of cube that contains all particles */

	for (i = 0; i < n; i++) {
		if (d[i].pos[0] < xmin)
			xmin = d[i].pos[0];
		if (d[i].pos[0] > xmax)
			xmax = d[i].pos[0];
		if (d[i].pos[1] < ymin)
			ymin = d[i].pos[1];
		if (d[i].pos[1] > ymax)
			ymax = d[i].pos[1];
		if (d[i].pos[2] < zmin)
			zmin = d[i].pos[2];
		if (d[i].pos[2] > zmax)
			zmax = d[i].pos[2];
		}

	dx = xmax - xmin;
	dy = ymax - ymin;
	dz = zmax - zmin;

	assert(dx >= 0.0 && dy >= 0.0 && dz >= 0.0);

	/* get minimum corner */
	
	if (xmin < ymin && xmin < zmin)
		dmin = xmin;
	else if (ymin < zmin && ymin < xmin)
		dmin = ymin;
	else
		dmin = zmin;

	dmax = dmin;

	/* get maximum dimension */
	
	if (dx > dy && dx > dz)
		dmax += dx;
	else if (dy > dz && dy > dx)
		dmax += dy;
	else
		dmax += dz;

	if (dmin == dmax) {
		fprintf(stderr, "Degenerate tree: all particles occupy same point in space.\n");
		return; /* root will be NULL */
		}

	/* it might be beneficial to center the root node, but ignore for now */
	
	make_node(dmin, dmax, dmin, dmax, dmin, dmax, root);

	/* fill tree with particles */

	for (i = 0; i < n; i++)
		add_to_tree(*root, &d[i]);
	}

static int rgb_to_color(int r, int g, int b)
{
	return r + 256 * (g + 256 * b);
	}

static void color_to_rgb(int iColor, int *r, int *g, int *b)
{
	*r = iColor & 0xff;
	*g = (iColor & 0xff00) / 256;
	*b = (iColor & 0xff0000) / (256 * 256);
	}

static int f_to_value(int iMin, int iMax, double c, double f)
{
	double f_cycle;
	
	f_cycle = f * c - (int) (f * c);
	if (f_cycle == 0.0 && f > 0.0)
		f_cycle = 1.0;
	
	return iMin + (int) ((iMax - iMin) * f_cycle);
	}

static int interpolate(const PARAMS_T *p, double f)
{
	switch (p->iScheme) {
	case SchemeSS:
		return f_to_value(p->iColorMin, p->iColorMax, p->dCycles, f);
	case SchemeRGB: {
		int rmax, gmax, bmax, rmin, gmin, bmin;

		color_to_rgb(p->iColorMax, &rmax, &gmax, &bmax);
		color_to_rgb(p->iColorMin, &rmin, &gmin, &bmin);

		return rgb_to_color(f_to_value(rmin, rmax, p->dCycles, f),
							f_to_value(gmin, gmax, p->dCycles, f),
							f_to_value(bmin, bmax, p->dCycles, f));
		}
	default:
		assert(0); /* this shouldn't happen */
		}
	}

static int randomize(const PARAMS_T *p)
{
	switch (p->iScheme) {
	case SchemeSS:
		return f_to_value(p->iColorMin, p->iColorMax, p->dCycles, randUniform());
	case SchemeRGB: {
		int rmax, gmax, bmax, rmin, gmin, bmin;
		double fr, fg, fb;

		color_to_rgb(p->iColorMax, &rmax, &gmax, &bmax);
		color_to_rgb(p->iColorMin, &rmin, &gmin, &bmin);

		fr = fg = fb = randUniform();
		if (p->bUnique == TRUE) {
			fg = randUniform();
			fb = randUniform();
			}

		return rgb_to_color(f_to_value(rmin, rmax, p->dCycles, fr),
							f_to_value(gmin, gmax, p->dCycles, fg),
							f_to_value(bmin, bmax, p->dCycles, fb));
		}
	default:
		assert(0); /* this shouldn't happen */
		}
	}

static void change_color(const PARAMS_T *p, int n, const SORT_DATA_T *sd, SSDATA *d)
{
	int i, iColor;

	/* loop over data */
	
	for (i = 0; i < n; i++) {
		iColor = d[sd[i].iIndex].color;
		if (sd[i].dData >= p->dVMin && sd[i].dData <= p->dVMax) {
			if (p->iField != FieldNone)
				if (p->dVMin == p->dVMax)
					iColor = interpolate(p, 0.5);
				else
					iColor = interpolate(p, (sd[i].dData - p->dVMin) / (p->dVMax - p->dVMin));
			else if (p->bRandom == TRUE)
				iColor = randomize(p);
			else if (p->iColor >= 0)
				iColor = p->iColor;
			}
		else if (p->iColor >= 0)
			iColor = p->iColor;
		d[sd[i].iIndex].color = iColor;
		}
	}

static int data_compare(const void *e1, const void *e2)
{
	return ((const SORT_DATA_T *)e1)->dData < ((const SORT_DATA_T *)e2)->dData ? -1 :
		   ((const SORT_DATA_T *)e1)->dData > ((const SORT_DATA_T *)e2)->dData ?  1 : 0;
	}

static void sort_data(const PARAMS_T *p, int n, SSDATA *d, SORT_DATA_T **sd)
{
	NODE *root = NULL;
	double dData;
	int i;

	if (p->iField == FieldSDen || p->iField == FieldOvlp ||	p->iField == FieldGPot)
		build_tree(n, d, &root); /* allocate tree */
	
	*sd = (SORT_DATA_T *) calloc((size_t) n, sizeof(SORT_DATA_T));
	assert(*sd != NULL);

	for (i = 0; i < n; i++) {
		(*sd)[i].iIndex = i;
		switch (p->iField) {
		case FieldNone:
			dData = i; /* i.e., sort will have no effect */
			break;
		case FieldOrgIdx:
			dData = d[i].org_idx;
			break;
		case FieldMass:
			dData = d[i].mass;
			if (p->bUseMKS)
				dData *= M_SCALE;
			break;
		case FieldRadius:
			dData = d[i].radius;
			if (p->bUseMKS)
				dData *= L_SCALE;
			break;
		case FieldPosX:
			dData = d[i].pos[0];
			if (p->bUseMKS)
				dData *= L_SCALE;
			break;
		case FieldPosY:
			dData = d[i].pos[1];
			if (p->bUseMKS)
				dData *= L_SCALE;
			break;
		case FieldPosZ:
			dData = d[i].pos[2];
			if (p->bUseMKS)
				dData *= L_SCALE;
			break;
		case FieldVelX:
			dData = d[i].vel[0];
			if (p->bUseMKS)
				dData *= V_SCALE;
			break;
		case FieldVelY:
			dData = d[i].vel[1];
			if (p->bUseMKS)
				dData *= V_SCALE;
			break;
		case FieldVelZ:
			dData = d[i].vel[2];
			if (p->bUseMKS)
				dData *= V_SCALE;
			break;
		case FieldSpinX:
			dData = d[i].spin[0];
			if (p->bUseMKS)
				dData *= W_SCALE;
			break;
		case FieldSpinY:
			dData = d[i].spin[1];
			if (p->bUseMKS)
				dData *= W_SCALE;
			break;
		case FieldSpinZ:
			dData = d[i].spin[2];
			if (p->bUseMKS)
				dData *= W_SCALE;
			break;
		case FieldColor:
			dData = d[i].color;
			break;
		case FieldRMag:
			dData = sqrt(d[i].pos[0] * d[i].pos[0] +
						 d[i].pos[1] * d[i].pos[1] +
						 d[i].pos[2] * d[i].pos[2]);
			if (p->bUseMKS)
				dData *= L_SCALE;
			break;
		case FieldVMag:
			dData = sqrt(d[i].vel[0] * d[i].vel[0] +
						 d[i].vel[1] * d[i].vel[1] +
						 d[i].vel[2] * d[i].vel[2]);
			if (p->bUseMKS)
				dData *= V_SCALE;
			break;
		case FieldWMag:
			dData = sqrt(d[i].spin[0] * d[i].spin[0] +
						 d[i].spin[1] * d[i].spin[1] +
						 d[i].spin[2] * d[i].spin[2]);
			if (p->bUseMKS)
				dData *= W_SCALE;
			break;
		case FieldMDen:
			dData = d[i].mass / (4.0 / 3.0 * M_PI * d[i].radius * d[i].radius * d[i].radius);
			if (p->bUseMKS)
				dData *= D_SCALE;
			break;
		case FieldSDen:
			fprintf(stderr, "*** sden not yet implemented ***\n");
			exit(1);
		case FieldOvlp:
			dData = get_max_overlap(&d[i], root);
			break;
		case FieldGPot:
			fprintf(stderr, "*** gpot not yet implemented ***\n");
			exit(1);
		default:
			assert(0); /* this shouldn't happen */
			}
		/* handle log options */
		switch (p->iLogOpt) {
		case LogNone:
			break;
		case LogStd:
			if (dData <= 0.0 && p->dVMin == - HUGE_VAL) {
				fprintf(stderr, "*** cannot take logarithm of particle %i data (%g) ***\n", i, dData);
				exit(1);
				}
			if (dData <= 0.0)
				dData = - HUGE_VAL;
			else
				dData = log10(dData);
			break;
		case LogAbs:
			if (dData == 0.0 && p->dVMin == - HUGE_VAL) {
				fprintf(stderr, "*** cannot take logarithm of particle %i data (%g) ***\n", i, dData);
				exit(1);
				}
			if (dData == 0.0)
				dData = - HUGE_VAL;
			else
				dData = log10(fabs(dData));
			break;
		default:
			assert(0); /* this shouldn't happen */
			}
		(*sd)[i].dData = dData;					
		}

	qsort((void *) (*sd), (size_t) n, sizeof(SORT_DATA_T), data_compare);

	if (root != NULL)
		kill_node(root); /* deallocate tree */
	}

static int get_color_from_str(int iScheme, const char *achStr, int *iColor)
{
	switch (iScheme) {
	case SchemeSS:
		if (strstr(achStr, ",") != NULL)
			fprintf(stderr, "*** ss colormap value: everything after ',' ignored (%s) ***\n", achStr);
		*iColor = atoi(achStr);
		if (*iColor < 0 || *iColor > 255) {
			fprintf(stderr, "*** invalid ss color (%i) ***\n", *iColor);
			return 1;
			}
		break;
	case SchemeRGB: {
		float r, g, b;
		if (sscanf(achStr, "%f,%f,%f", &r, &g, &b) != 3) {
			fprintf(stderr, "*** error parsing RGB string (%s) ***\n", achStr);
			return 1;
			}
		if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			fprintf(stderr, "*** invalid RGB value(s) (%g,%g,%g) ***\n", r, g, b);
			return 1;
			}
		/* following convention means color (1,1,1) is not supported */
		r = (r <= 1.0 ? 255 * r : r);
		g = (g <= 1.0 ? 255 * g : g);
		b = (b <= 1.0 ? 255 * b : b);
		*iColor = rgb_to_color((int) r, (int) g, (int) b);
		break;
		}
	default:
		assert(0); /* this shouldn't happen */
		}

	return 0;
	}
	
static int save_data(const char *achFilename, SSHEAD *h, SSDATA *d)
{
	SSIO ssio;
	int i;
	
	if (ssioOpen(achFilename, &ssio, SSIO_UPDATE) != 0) {
		fprintf(stderr, "Unable to open \"%s\" for writing.\n", achFilename);
		return 1;
		}
	
	if (ssioHead(&ssio, h) != 0) {
		fprintf(stderr, "%s: error writing header.\n", achFilename);
		goto abort;
		}
	
	switch(h->iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		for (i = 0; i < h->n_data; i++) {
			if (ssioData(&ssio, &(d[i])) != 0) {
				fprintf(stderr, "%s: error writing data (particle index %i).\n", achFilename, i);
				goto abort;
				}
			}
		break;
	case SSIO_MAGIC_REDUCED: {
		SSRDATA rd;
		for (i = 0; i < h->n_data; i++) {
			ssioStandardToReduced((const SSDATA *) &(d[i]), &rd);
			if (ssioDataReduced(&ssio, &rd) != 0) {
				fprintf(stderr, "%s: error writing data (particle index %i).\n", achFilename, i);
				goto abort;
				}
			}
		}
		break;
	default:
		assert(0); /* this shouldn't happen */
		}

	ssioClose(&ssio);
	return 0;
	
 abort:
	ssioClose(&ssio);
	return 1;
	}

static int read_data(const char *achFilename, SSHEAD *h, SSDATA **d)
{
	SSIO ssio;
	int i;

	*d = NULL;
	
	if (ssioOpen(achFilename, &ssio, SSIO_READ) != 0) {
		fprintf(stderr, "Unable to open \"%s\".\n", achFilename);
		return 1;
		}
	
	if (ssioHead(&ssio, h) != 0) {
		fprintf(stderr, "%s: corrupt header.\n", achFilename);
		goto abort;
		}

	printf("%s: time = %g, n_data = %i%s\n", achFilename, h->time, h->n_data, h->iMagicNumber == SSIO_MAGIC_REDUCED ? " [reduced format]" : "");

	if (h->n_data < 0) {
		fprintf(stderr, "%s: invalid input data format.\n", achFilename);
		goto abort;
		}

	if (h->n_data == 0) {
		fprintf(stderr, "%s: no particles found!\n", achFilename);
		goto abort;
		}
	
	*d = (SSDATA *) calloc((size_t) h->n_data, sizeof(SSDATA));
	if (*d == NULL) {
		fprintf(stderr, "%s: unable to allocate memory for data.\n", achFilename);
		goto abort;
		}
	
	switch(h->iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		for (i = 0; i < h->n_data; i++) {
			if (ssioData(&ssio, &((*d)[i])) != 0) {
				fprintf(stderr, "%s: corrupt data (particle index %i).\n", achFilename, i);
				goto abort;
				}
			}
		break;
	case SSIO_MAGIC_REDUCED: {
		SSRDATA rd;
		for (i = 0; i < h->n_data ; i++) {
			if (ssioDataReduced(&ssio, &rd) != 0) {
				fprintf(stderr, "%s: corrupt data (particle index %i).\n", achFilename, i);
				goto abort;
				}
			ssioReducedToStandard(&rd, &((*d)[i]));
			}
		break;
		}
	default:
		fprintf(stderr, "%s: unrecognized magic number (%i).\n", achFilename, h->iMagicNumber);
		goto abort;
		}

	ssioClose(&ssio);
	return 0;
	
 abort:
	if (*d != NULL)
		free((void *) (*d));
	ssioClose(&ssio);
	return 1;
	}

static void usage(const char *achProgName)
{
	fprintf(stderr,
			"Usage: %s [ options ] [ toggles ] ssfile [ ssfile ... ]\n",
			achProgName);
	fprintf(stderr,	"Use '--help' for a list of options and toggles, with examples.\n");
	exit(1);
	}

static void help(const char *achProgName)
{
	fprintf(stderr, "%s is a tool for changing particle colors in place in ss files.\n", achProgName);
	fprintf(stderr, "Usage: %s [ options ] [ toggles ] ssfile [ ssfile ... ]\n", achProgName);
	fprintf(stderr, "1. The configurable options are:\n");
	fprintf(stderr, "   --color=value [or -o value]\n");
	fprintf(stderr, "     single integer between 0 and 255 for ss color scheme (default)\n");
	fprintf(stderr, "     triples (no spaces, e.g., 180,0.5,0.5 or 0.5,0.5,0.5) for HSV\n");
	fprintf(stderr, "     triples (no spaces, e.g., 127,255,0 or 0.5,1.0,0.0) for RGB\n");
	fprintf(stderr, "   --field=value [or -f value]\n");
	fprintf(stderr, "     oi: original index\n");
	fprintf(stderr, "     m:  particle mass\n");			
	fprintf(stderr, "     r:  particle radius\n");
	fprintf(stderr, "     x:  x position\n");
	fprintf(stderr, "     y:  y position\n");
	fprintf(stderr, "     z:  z position\n");
	fprintf(stderr, "     vx: x velocity component\n");
	fprintf(stderr, "     vy: y velocity component\n");
	fprintf(stderr, "     vz: z velocity component\n");
	fprintf(stderr, "     wx: x spin component\n");
	fprintf(stderr, "     wy: y spin component\n");
	fprintf(stderr, "     wz: z spin component\n");
	fprintf(stderr, "     c:  particle color\n");
	fprintf(stderr, "     Rmag: particle distance from origin\n");
	fprintf(stderr, "     Vmag: particle absolute speed\n");
	fprintf(stderr, "     Wmag: particle spin magnitude\n");
	fprintf(stderr, "     mden: particle mass density\n");
	fprintf(stderr, "     sden: particle spatial density\n");
	fprintf(stderr, "     ovlp: particle overlap fraction\n");
	fprintf(stderr, "     gpot: particle gravitational potential\n");
	fprintf(stderr, "   --cmax=value [default red in SS and HSV/RGB schemes]\n");
  	fprintf(stderr, "     color to assign particles with maximum field value\n");
	fprintf(stderr, "   --cmin=value [default navy in SS, blue in HSV/RGB schemes]\n");
	fprintf(stderr, "     color to assign particles with minimum field value\n");
	fprintf(stderr, "   --cycles=value [or -c value; default 1]\n");
	fprintf(stderr, "     number of color range cycles to use between data range\n");
	fprintf(stderr, "   --vmax=value [pkdgrav units]\n");
	fprintf(stderr, "     maximum field value to associate with cmax\n");
	fprintf(stderr, "   --vmin=value [pkdgrav units]\n");
	fprintf(stderr, "     minimum field value to associate with cmin\n");
	fprintf(stderr, "   --log[=value]: take log of data, map from log(vmin) to log(vmax)\n");
	fprintf(stderr, "     abs: take log of absolute value of field data\n");
	fprintf(stderr, "2. The toggles are:\n");
	fprintf(stderr, "   --aggs [or -a]: color by aggs (if present), not just particles\n");
	fprintf(stderr, "   --mks [or -m]: --vmin and --vmax in MKS, not pkdgrav units\n");
	fprintf(stderr, "   --random [or -r]: color randomly in --cmin, --cmax range\n");
	fprintf(stderr, "   --sorteach [or -s]: sort each file (default: only sort first)\n");
	fprintf(stderr, "   --unique [or -u]: random HSV/RGB colors truly random in requested range\n");
	fprintf(stderr, "   --hsv: use HSV (hue-saturation-value) color space\n");
	fprintf(stderr, "   --rgb: use RGB (red-green-blue) color space\n");
	fprintf(stderr, "   Without --hsv or --rgb, the standard 8-bit ss colormap is used.\n");
	fprintf(stderr, "   Use ssdraw -c with --hsv or --rgb (24-bit support).\n");
	fprintf(stderr, "   Note --hsv or --rgb must appear BEFORE --color, --cmax, or --cmin.\n");
	fprintf(stderr, "3. Examples:\n");
	fprintf(stderr, "   a) Color by initial z position using standard map:\n");
    fprintf(stderr, "      sscolor -f z myfiles*.ss\n");
	fprintf(stderr, "   b) Color each ss file by speed from blue to yellow logarithmically:\n");
	fprintf(stderr, "      sscolor -f vmag --log --rgb --cmin=0,0,255 --cmax=255,255,0 --sorteach myfiles*.ss\n");
	fprintf(stderr, "   c) Color all particles 0.1 m or larger in radius red, otherwise green:\n");
	fprintf(stderr, "      sscolor -f R -o 3 --mks --vmin=0.1 --cmin=2 --cmax=2 myfile.ss\n");
	
	exit(0);
	}

int main(int argc, char *argv[])
{
	/* options descriptor */
	static struct option longopts[] = {
	  {"aggs", no_argument, NULL, 'a'},
	  {"cmax", required_argument, NULL, 'p'},
	  {"cmin", required_argument, NULL, 'q'},
	  {"color", required_argument, NULL, 'o'},
	  {"cycles", required_argument, NULL, 'c'},
	  {"field", required_argument, NULL, 'f'},
	  {"help", no_argument, NULL, 'h'},
	  {"log", optional_argument, NULL, 'l'},
	  {"mks", no_argument, NULL, 'm'},
	  {"random", no_argument, NULL, 'r'},
	  {"sorteach", no_argument, NULL, 's'},
	  {"unique", no_argument, NULL, 'u'},
	  {"vmax", required_argument, NULL, 'v'},
	  {"vmin", required_argument, NULL, 'w'},
	  {"hsv", no_argument, NULL, 'x'},
	  {"rgb", no_argument, NULL, 'y'},
	  {NULL, 0, NULL, 0}
	};

	PARAMS_T p;
	SSHEAD h;
	SSDATA *d;
	SORT_DATA_T *sd = NULL;
	int i, n_sort = -1;

	/* seed random number generator */

	randSeedGenerator((Ullong) (time(NULL) % getpid() + getppid()));
	
	/* defaults (in longopts order for consistency) */

	p.bColorAggs = FALSE;
	p.iColorMax = -1;
	p.iColorMin = -1;
	p.iColor = -1;
	p.dCycles = 1.0;
	p.iField = FieldNone;
	p.iLogOpt = LogNone;
	p.bUseMKS = FALSE;
	p.bRandom = FALSE;
	p.bSortEach = FALSE;
	p.bUnique = FALSE;
	p.dVMax = HUGE_VAL;
	p.dVMin = - HUGE_VAL;
	p.iScheme = SchemeSS;
	
	/* read parameters (arguments) */
	
	while ((i = getopt_long(argc, argv, "ac:f:mo:rsu", longopts, NULL)) != -1) {
		switch (i) {
		case 'a':
			p.bColorAggs = TRUE;
			fprintf(stderr, "*** aggs not yet implemented ***\n");
			return 1;
		case 'c':
			p.dCycles = atof(optarg);
			if (p.dCycles <= 0.0) {
				fprintf(stderr, "*** cycles (%s) must be a positive real number ***\n", optarg);
				usage(argv[0]);
				}
			break;
		case 'f':
			p.iField = -1;
			for (i = 0; i < NUM_FIELDS; i++)
				if (strcasecmp(optarg, achFieldTags[i]) == 0)
					p.iField = i;
			if (p.iField == -1) {
				fprintf(stderr, "*** unrecognized field (%s) ***\n", optarg);
				usage(argv[0]);
				}
			break;
		case 'h':
			help(argv[0]);
			break; /* to suppress fall-through warning */
		case 'l':
			if (p.iLogOpt != LogNone) {
				fprintf(stderr, "*** only one log option may be selected at a time ***\n");
				usage(argv[0]);
				}
			p.iLogOpt = LogStd; /* default */
			if (optarg != NULL) {
				p.iLogOpt = -1;
				for (i = 0; i < NUM_LOG_OPTS; i++)
					if (strcasecmp(optarg, achLogOptTags[i]) == 0)
						p.iLogOpt = i;
				if (p.iLogOpt == -1) {
					fprintf(stderr, "*** unrecognized log option (%s) ***\n", optarg);
					usage(argv[0]);
					}
				}
			break;
		case 'm':
			p.bUseMKS = TRUE;
			break;
		case 'o':
			if (get_color_from_str(p.iScheme, optarg, &p.iColor) != 0)
				usage(argv[0]);
			break;
		case 'p':
			if (get_color_from_str(p.iScheme, optarg, &p.iColorMax) != 0)
				usage(argv[0]);
			break;
		case 'q':
			if (get_color_from_str(p.iScheme, optarg, &p.iColorMin) != 0)
				usage(argv[0]);
			break;
		case 'r':
			p.bRandom = TRUE;
			break;
		case 's':
			p.bSortEach = TRUE;
			break;
		case 'u':
			p.bUnique = TRUE;
			break;
		case 'v':
			p.dVMax = atof(optarg);
			break;
		case 'w':
			p.dVMin = atof(optarg);
			break;
		case 'x':
			if (p.iScheme != SchemeSS) {
				fprintf(stderr, "*** only one color scheme may be selected at a time ***\n");
				usage(argv[0]);
				}
			p.iScheme = SchemeHSV;
			fprintf(stderr, "*** HSV not yet implemented ***\n");
			return 1;
			break;
		case 'y':
			if (p.iScheme != SchemeSS) {
				fprintf(stderr, "*** only one color scheme may be selected at a time ***\n");
				usage(argv[0]);
				}
			p.iScheme = SchemeRGB;
			break;
		default:
			usage(argv[0]);
			}
		}

	/* more defaults */

	if (p.iField != FieldNone || p.bRandom == TRUE) {
		if (p.iColorMax == -1)
			switch (p.iScheme) {
			case SchemeSS:
				p.iColorMax = 2; /* red */
				break;
			case SchemeRGB:
				p.iColorMax = rgb_to_color(255, 0, 0); /* red */
				break;
			default:
				assert(0);
				}
		if (p.iColorMin == -1)
			switch (p.iScheme) {
			case SchemeSS:
				p.iColorMin = 15; /* navy */
				break;
			case SchemeRGB:
				p.iColorMin = rgb_to_color(0, 0, 255); /* blue */
				break;
			default:
				assert(0);
				}
		}
	
	/* check parameter consistency */

	if (p.dVMin > p.dVMax) {
		fprintf(stderr, "*** min value (%g) cannot be larger than max value (%g) ***\n", p.dVMin, p.dVMax);
		usage(argv[0]);
		}

	if (p.bRandom == TRUE && p.iField != FieldNone) {
		fprintf(stderr, "*** cannot use random with sort field ***\n");
		usage(argv[0]);
		}

	if (p.iLogOpt != LogNone && p.iField == FieldNone)
		fprintf(stderr, "*** warning: log option ignored since only 1 color specified ***\n");
		
	if (optind == argc) {
		fprintf(stderr, "*** ss file(s) not specified ***\n");
		usage(argv[0]);
		}

	/* display option summary */

	printf("Options:");
	if (p.bColorAggs == TRUE)
		printf(" --aggs");
	if (p.iColorMax >= 0)
		printf(" --cmax=%i", p.iColorMax);
	if (p.iColorMin >= 0)
		printf(" --cmin=%i", p.iColorMin);
	if (p.iColor >= 0)
		printf(" --color=%i", p.iColor);
	if (p.iField != FieldNone || p.bRandom == TRUE)
		printf(" --cycles=%g", p.dCycles);
	if (p.iField != FieldNone)
		printf(" --field=%s", achFieldTags[p.iField]);
	if (p.iLogOpt != LogNone) {
		printf(" --log");
		if (p.iLogOpt != LogStd)
			printf("=%s", achLogOptTags[p.iLogOpt]);
		}
	if (p.bUseMKS == TRUE)
		printf(" --mks");
	if (p.bRandom == TRUE)
		printf(" --random");
	if (p.bSortEach == TRUE)
		printf(" --sorteach");
	if (p.bUnique == TRUE)
		printf(" --unique");
	if (p.dVMax < HUGE_VAL)
		printf(" --vmax=%g", p.dVMax);
	if (p.dVMin > - HUGE_VAL)
		printf(" --vmin=%g", p.dVMin);
	switch (p.iScheme) {
	case SchemeSS:
		break;
	case SchemeHSV:
		printf(" --hsv");
		break;
	case SchemeRGB:
		printf(" --rgb");
		break;
	default:
		assert(0); /* this shouldn't happen */
		}
	printf("\n");
	
	/* handle log options */

	if (p.iField != FieldNone && (p.dVMax != HUGE_VAL || p.dVMin != - HUGE_VAL))
		switch (p.iLogOpt) {
		case LogNone:
			break;
		case LogStd:
			if (p.dVMax <= 0.0) {
				fprintf(stderr, "*** cannot take logarithm of vmax (%g) ***\n", p.dVMax);
				return 1;
				}
			if (p.dVMax < HUGE_VAL)
				p.dVMax = log10(p.dVMax);
			if (p.dVMin > - HUGE_VAL && p.dVMin <= 0.0) {
				fprintf(stderr, "*** cannot take logarithm of vmin (%g) ***\n", p.dVMin);
				return 1;
				}
			if (p.dVMin > - HUGE_VAL)
				p.dVMin = log10(p.dVMin);
			break;
		case LogAbs:
			if (p.dVMax == 0.0) {
				fprintf(stderr, "*** cannot take logarithm of vmax (%g) ***\n", p.dVMax);
				return 1;
				}
			if (p.dVMax < HUGE_VAL)
				p.dVMax = log10(fabs(p.dVMax));
			if (p.dVMin == 0.0) {
				fprintf(stderr, "*** cannot take logarithm of vmin (%g) ***\n", p.dVMin);
				return 1;
				}
			if (p.dVMin > - HUGE_VAL)
				p.dVMin = log10(fabs(p.dVMin));
			break;
		default:
			assert(0); /* this shouldn't happen */
			}
	
	/* loop through ss files */

	for (i = optind; i < argc; i++) {
		if (read_data(argv[i], &h, &d) != 0) {
			if (i == optind) {
				fprintf(stderr, "*** Unable to read starting file ***\n");
				return 1;
				}
			else
				continue;
			}
		if (i == optind || p.bSortEach == TRUE) {
			sort_data(&p, h.n_data, d, &sd);
			if (p.iField != FieldNone) {
				/* assign default max/min values if needed */
				if (p.dVMax == HUGE_VAL)
					p.dVMax = sd[h.n_data - 1].dData;
				if (p.dVMin == - HUGE_VAL)
					p.dVMin = sd[0].dData;
				printf("Data: max value = %g, min value = %g\n", sd[h.n_data - 1].dData, sd[0].dData);
				}
			n_sort = h.n_data;
			}
		if (h.n_data != n_sort) {
			fprintf(stderr, "*** %s: number of particles differs from first ss file ***\n", argv[i]);
			return 1;
			}
		change_color(&p, h.n_data, sd, d);
		save_data(argv[i], &h, d);
		free((void *) d);
		}

	return 0;
	}

/* sscolor.c */
