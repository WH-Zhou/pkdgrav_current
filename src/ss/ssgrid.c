/*
 ** ssgrid.c -- DCR 7/17/14
 ** ========
 ** Utility to locate uppermost (largest z) particles in an x-y grid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for getopt() */
#include <assert.h>
#include <ss.h>

#define OUTFILENAME "ssgrid.dat"

typedef int BOOLEAN_T;
enum {False=0, True};

typedef struct {
  /* user-supplied parameters */
  double dXMin, dXMax, dYMin, dYMax, dZMax, dConvert;
  int nX, nY;
  BOOLEAN_T bVolume, bXAvg, bYAvg;
  /* derived parameters */
  double dLenX, dLenY, dX, dY;
  } PARAMS_T;

double GRID_XCOOR(int iX);
#define GRID_XCOOR(iX) ((p->dXMin + (iX) * p->dX + 0.5 * p->dX) * dConvert)

double GRID_YCOOR(int iY);
#define GRID_YCOOR(iY) ((p->dYMin + (iY) * p->dY + 0.5 * p->dY) * dConvert)

double GRID_VALUE(const SSDATA *d);
#define GRID_VALUE(d) ((d) == NULL ? 0.0 : ((d)->pos[2] + (d)->radius) * dConvert)

static int output_grid(const PARAMS_T *p, SSDATA ***g) {
	FILE *fp;
	double dConvert;
	int iX, iY, iStatus;

	dConvert = 1.0 / p->dConvert; /* from system to user units */

	fp = fopen(OUTFILENAME, "w");
	if (fp == NULL) {
		fprintf(stderr, "Unable to open \"%s\" for writing.\n", OUTFILENAME);
		return 1;
		}
	iStatus = 0;
	if (fprintf(fp, "# xmin: %g xmax: %g ymin: %g ymax: %g nx: %i ny: %i zmax: %g convert: %g xavg: %i yavg: %i\n"
				"# Each output line is x, y center of cell and maximum z + R value in bin (else 0).\n"
				"# Note x and/or y omitted if averaging over that dimension.\n",
				p->dXMin * dConvert, p->dXMax * dConvert, p->dYMin * dConvert, p->dYMax * dConvert,
				p->nX, p->nY, p->dZMax * dConvert, p->dConvert, p->bXAvg, p->bYAvg) <= 0) {
		iStatus = 1;
		goto fill_grid_finish;
		}
	if (p->bXAvg && p->bYAvg) {
		double dAvg = 0.0;

		for (iX = 0; iX < p->nX; iX++)
			for (iY = 0; iY < p->nY; iY++)
				dAvg += GRID_VALUE(g[iX][iY]);

		dAvg /= (p->nX * p->nY);

		if (fprintf(fp, "%.16e\n", dAvg) <= 0) {
			iStatus = 1;
			goto fill_grid_finish;
			}
		}
	else if (p->bXAvg) {
		double dAvg;

		for (iY = 0; iY < p->nY; iY++) {
			dAvg = 0.0;
			for (iX = 0; iX < p->nX; iX++)
				dAvg += GRID_VALUE(g[iX][iY]);
			dAvg /= p->nX;
			if (fprintf(fp, "%.16e %.16e\n", GRID_YCOOR(iY), dAvg) <= 0) {
				iStatus = 1;
				goto fill_grid_finish;
				}
			}
		}
	else if (p->bYAvg) {
		double dAvg;

		for (iX = 0; iX < p->nX; iX++) {
			dAvg = 0.0;
			for (iY = 0; iY < p->nY; iY++)
				dAvg += GRID_VALUE(g[iX][iY]);
			dAvg /= p->nY;
			if (fprintf(fp, "%.16e %.16e\n", GRID_XCOOR(iX), dAvg) <= 0) {
				iStatus = 1;
				goto fill_grid_finish;
				}
			}
		}
	else {
		for (iX = 0; iX < p->nX; iX++)
			for (iY = 0; iY < p->nY; iY++)
				if (fprintf(fp, "%.16e %.16e %.16e\n", GRID_XCOOR(iX), GRID_YCOOR(iY), GRID_VALUE(g[iX][iY])) <= 0) {
					iStatus = 1;
					goto fill_grid_finish;
					}
 		}
 fill_grid_finish:
	fclose(fp);
	return iStatus;
	}

#undef GRID_VALUE
#define GRID_VALUE(d) ((d) == NULL ? 0.0 : ((d)->pos[2] + (d)->radius))

static int get_volume(const PARAMS_T *p, SSDATA ***g) {
	double dMultFac, dVolume;
	int iX, iY;


	dMultFac = p->dX * p->dY / (p->dConvert * p->dConvert * p->dConvert);
	dVolume = 0.0;
	for (iX = 0; iX < p->nX; iX++)
		for (iY = 0; iY < p->nY; iY++)
			dVolume += dMultFac * (p->dZMax - GRID_VALUE(g[iX][iY]));

	printf("Volume = %g cubic length units.\n", dVolume);

	return 0;
	}

#undef GRID_VALUE

static void fill_grid(const PARAMS_T *p, SSDATA d[], int n, SSDATA ***g)
{
	int iX, iY, i;

	for (iX = 0; iX < p->nX; iX++)
		for (iY = 0; iY < p->nY; iY++)
			g[iX][iY] = NULL;
	for (i = 0; i < n; i++) {
		/* use '>=' in following so bins work out... */
		if (d[i].pos[0] < p->dXMin || d[i].pos[0] >= p->dXMax || d[i].pos[1] < p->dYMin || d[i].pos[1] >= p->dYMax || d[i].pos[2] + d[i].radius > p->dZMax)
			continue;
		iX = (int) (p->nX * (d[i].pos[0] - p->dXMin) / p->dLenX);
		iY = (int) (p->nY * (d[i].pos[1] - p->dYMin) / p->dLenY);
		if (g[iX][iY] == NULL || d[i].pos[2] + d[i].radius > g[iX][iY]->pos[2] + g[iX][iY]->radius)
			g[iX][iY] = &d[i];
		}
	}

static int read_data(const char *achFilename, SSDATA **d, int *n)
{
	SSIO ssio;
	SSHEAD head;
	int iStatus, i;

	printf("%s\n", achFilename);

	if (ssioOpen(achFilename, &ssio, SSIO_READ)) {
		fprintf(stderr, "Unable to open \"%s\".\n", achFilename);
		return 1;
		}

	iStatus = 0;

	if (ssioHead(&ssio,&head)) {
		fprintf(stderr, "Corrupt header.\n");
		iStatus = 1;
		goto read_data_finish;
		}

	printf("[time = %e, n_data = %i]\n", head.time, head.n_data);

	if (head.iMagicNumber == SSIO_MAGIC_REDUCED)
		printf("[reduced format]\n");

	if (head.n_data <= 0) {
		fprintf(stderr, "Invalid input data format.\n");
		iStatus = 1;
		goto read_data_finish;
		}

	*n = head.n_data;
	*d = calloc((size_t) *n, sizeof(SSDATA));
	if (*d == NULL) {
		fprintf(stderr, "Unable to allocate memory for data.\n");
		iStatus = 1;
		goto read_data_finish;
		}

	switch(head.iMagicNumber) {
	case SSIO_MAGIC_STANDARD: {
		for (i = 0; i < *n; i++) {
			if (ssioData(&ssio, &((*d)[i]))) {
				fprintf(stderr, "Corrupt data (i = %i).\n", i);
				iStatus = 1;
				goto read_data_finish;
				}
			}
		break;
		}
	case SSIO_MAGIC_REDUCED: {
		SSRDATA rd;
		for (i = 0; i < *n ; i++) {
			if (ssioDataReduced(&ssio, &rd)) {
				fprintf(stderr, "Corrupt data (i = %i).\n", i);
				iStatus = 1;
				goto read_data_finish;
				}
			(*d)[i].mass = rd.fMass;
			(*d)[i].radius = rd.fRadius;
			(*d)[i].pos[0] = rd.vPos[0];
			(*d)[i].pos[1] = rd.vPos[1];
			(*d)[i].pos[2] = rd.vPos[2];
			(*d)[i].color = rd.iColor;
			(*d)[i].org_idx = rd.iOrgIdx;
			}
		break;
		}
	default:
		fprintf(stderr, "Unrecognized ss file magic number (%i).\n", head.iMagicNumber);
		iStatus = 1;
		}

 read_data_finish:
	ssioClose(&ssio);
	return iStatus;
	}

static void usage(const char *achProgName)
{
	fprintf(stderr, "Usage: %s --xmin=xmin --xmax=xmax --ymin=ymin --ymax=ymax --nx=nx --ny=ny [--zmax=zmax [--volume]] [--convert=lscale] [--xavg] [--yavg] ssfile\n", achProgName);
	fprintf(stderr, "       where xmin, xmax, ymin, ymax define grid,\n");
	fprintf(stderr, "             nx, ny are number of grid cells in each dimension,\n");
	fprintf(stderr, "             zmax is the maximum z value to consider for particles (dflt 0),\n");
	fprintf(stderr, "             volume indicates report volume between zmax and surface,\n");
	fprintf(stderr, "             lscale multiplies length units (else uses system units),\n");
	fprintf(stderr, "             xavg and/or yavg mean take average in that dimension,\n");
	fprintf(stderr, "             and ssfile is the file to process.\n");
	fprintf(stderr, "             (Note: cannot specify both volume and xavg and/or yavg.)\n");
	fprintf(stderr, "Output data is x, y, z, R for particle with largest z + R in each cell.\n");
	fprintf(stderr, "Output filename = \"%s\".\n", OUTFILENAME);

	exit(1);
}

int main(int argc, char *argv[])
{
	/* options descriptor */
	static struct option longopts[] = {
	  {"xmin", required_argument, NULL, 'd'},
	  {"xmax", required_argument, NULL, 'e'},
	  {"ymin", required_argument, NULL, 'f'},
	  {"ymax", required_argument, NULL, 'g'},
	  {"nx", required_argument, NULL, 'x'},
	  {"ny", required_argument, NULL, 'y'},
	  {"zmax", required_argument, NULL, 'z'},
	  {"volume", no_argument, NULL, 'v'},
	  {"convert", required_argument, NULL, 'c'},
	  {"xavg", no_argument, NULL, 'a'},
	  {"yavg", no_argument, NULL, 'b'},
	  {NULL, 0, NULL, 0}
	};

	PARAMS_T p;
	SSDATA *d = NULL, ***g = NULL /* a 2D array of pointers to particle data */ ;
	int n, i, iStatus;
	
	p.dXMin = p.dXMax = p.dYMin = p.dYMax = p.dZMax = 0.0;
	p.nX = p.nY = 0;
	p.dConvert = 1.0;
	p.bVolume = p.bXAvg = p.bYAvg = False;

	while ((i = getopt_long(argc, argv, "d:e:f:g:x:y:z:c:ab", longopts, NULL)) != -1)
		switch (i) {
		case 'd':
			p.dXMin = atof(optarg);
			break;
		case 'e':
			p.dXMax = atof(optarg);
			break;
		case 'f':
			p.dYMin = atof(optarg);
			break;
		case 'g':
			p.dYMax = atof(optarg);
			break;
		case 'x':
			p.nX = atoi(optarg);
			break;
		case 'y':
			p.nY = atoi(optarg);
			break;
		case 'z':
			p.dZMax = atof(optarg);
			break;
		case 'v':
			p.bVolume = True;
			break;
		case 'c':
			p.dConvert = atof(optarg);
			break;
		case 'a':
			p.bXAvg = True;
			break;
		case 'b':
			p.bYAvg = True;
			break;
		default:
			usage(argv[0]);
			}

	if (p.dXMin >= p.dXMax || p.dYMin >= p.dYMax || p.nX <= 0 || p.nY <= 0 || p.dConvert <= 0.0 || (p.bVolume && (p.bXAvg || p.bYAvg))) {
		if (p.dXMin >= p.dXMax)
			fprintf(stderr,"ERROR: xmin must be less than xmax.\n");
		if (p.dYMin >= p.dYMax)
			fprintf(stderr,"ERROR: ymin must be less than ymax.\n");
		if (p.nX <= 0)
			fprintf(stderr,"ERROR: nx must be positive.\n");
		if (p.nY <= 0)
			fprintf(stderr,"ERROR: ny must be positive.\n");
		if (p.dConvert <= 0.0)
			fprintf(stderr,"ERROR: convert must be positive.\n");
		if (p.bVolume && (p.bXAvg || p.bYAvg))
			fprintf(stderr,"ERROR: cannot specify both volume and xavg and/or yavg.\n");
		return 1;
		}
	
	p.dXMin *= p.dConvert;
	p.dXMax *= p.dConvert;
	p.dYMin *= p.dConvert;
	p.dYMax *= p.dConvert;
	p.dZMax *= p.dConvert;

	p.dLenX = p.dXMax - p.dXMin;
	p.dLenY = p.dYMax - p.dYMin;
	p.dX = p.dLenX / p.nX;
	p.dY = p.dLenY / p.nY;

	if (argc - optind != 1) usage(argv[0]);

	iStatus = 0;

	d = NULL;

	if (read_data(argv[optind], &d, &n) != 0) {
		iStatus = 1;
		goto main_finish;
		}

	{
		/* quality check: bins must be at least as big as particles... */

		double dRadMax = 0.0;

		for (i = 0; i < n; i++)
			if (2.0*d[i].radius > p.dX || 2.0*d[i].radius > p.dY)
				if (d[i].radius > dRadMax)
					dRadMax = d[i].radius;

		if (dRadMax > 0.0) {
			fprintf(stderr, "Bin sizes (dX, dY = %g, %g) smaller than largest particle (diameter %g) (scaled units).\n",
					p.dX / p.dConvert, p.dY / p.dConvert, 2.0 * dRadMax / p.dConvert);
			fprintf(stderr, "Recommended maximum (nX, nY) for given dimensions = (%i, %i).\n",
					(int) (0.5 * p.dLenX / dRadMax), (int) (0.5 * p.dLenY / dRadMax));
			iStatus = 1;
			goto main_finish;
			}
		}

	g = calloc((size_t) p.nX, sizeof(void *));
	if (g == NULL) {
		fprintf(stderr, "Unable to allocate memory for grid.\n");
		iStatus = 1;
		goto main_finish;
		}
	for (i = 0; i < p.nX; i++) {
		g[i] = calloc((size_t) p.nY, sizeof(void *));
		if (g[i] == NULL) {
			fprintf(stderr, "Unable to allocate memory for grid (row %i).\n", i);
			iStatus = 1;
			goto main_finish;
			}
		}

	fill_grid(&p, d, n, g);

	if ((p.bVolume && get_volume(&p, g) != 0) || output_grid(&p, g) != 0) {
		iStatus = 1;
		goto main_finish;
		}

 main_finish:
	if (d != NULL)
		free(d);

	if (g != NULL) {
		for (i = 0; i < p.nX; i++)
			if (g[i] != NULL)
				free(g[i]);
		free(g);
		}

	return iStatus;
	}

/* ssgrid.c */
