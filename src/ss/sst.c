/*
 ** sst.c -- DCR 7/1/03
 ** =====
 ** Traces particle evolution from ssg map file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>     /* for getopt() */
#include <sys/param.h>	/* for MAXPATHLEN */
#include <assert.h>
#include <ss.h>

#define VERBOSE /*DEBUG this should really be a run-time option*/

#define TRACEFILENAME "trace.ss"

typedef struct {
	int iIdx,iOrgIdx;
	double dMass;
	} BIG;

static int trace(FILE *fp,char *ssfile_in,int iOrgIdx)
{
	SSIO ssio_in,ssio_out;
	SSHEAD h;
	SSDATA d;
	int i,rv,nin,nout,idx;


	if (ssioOpen(ssfile_in,&ssio_in,SSIO_READ)) {
		(void) fprintf(stderr,"trace(): Unable to open \"%s\" for reading.\n",ssfile_in);
		return 1;
		}

	if (ssioOpen(TRACEFILENAME,&ssio_out,SSIO_WRITE)) {
		(void) fprintf(stderr,"trace(): Unable to open \"%s\" for writing.\n",TRACEFILENAME);
		return 1;
		}

	if (ssioHead(&ssio_in,&h)) {
		(void) fprintf(stderr,"trace(): Error reading \"%s\" header.\n",ssfile_in);
		return -1;
		}

	switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		ssioClose(&ssio_in);
		return 1;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio_in);
		return 1;
		}

	if (ssioHead(&ssio_out,&h)) { /* dummy header */
		(void) fprintf(stderr,"trace(): Error writing \"%s\" header.\n",TRACEFILENAME);
		return -1;
		}

	nin = h.n_data;

	for (i=nout=0;i<nin;i++) {
		/* map file should have as many entries as particles in first ss file */
		rv = fscanf(fp,"%i",&idx);
		assert(rv == 1);
		if (ssioData(&ssio_in,&d)) {
			(void) fprintf(stderr,"trace(): Error reading \"%s\" data.\n",ssfile_in);
			return -1;
			}
		if (idx == iOrgIdx) {
#ifdef VERBOSE
			(void) printf("%i --> %i\n",i,iOrgIdx);
#endif
			if (ssioData(&ssio_out,&d)) {
				(void) fprintf(stderr,"trace(): Error writing \"%s\" data.\n",TRACEFILENAME);
				return -1;
				}
			++nout;
			}
		}

	(void) ssioRewind(&ssio_out);

	h.n_data = nout;

	if (ssioHead(&ssio_out,&h)) { /* correct header */
		(void) fprintf(stderr,"trace(): Error writing \"%s\" header.\n",TRACEFILENAME);
		return -1;
		}

	(void) ssioClose(&ssio_out);
	(void) ssioClose(&ssio_in);

	return 0;
	}

static int compar(const void *v1,const void *v2)
{
	const BIG *p1 = (const BIG *)v1,*p2 = (const BIG *)v2;

	assert(p1->iIdx != p2->iIdx);

	return (p1->dMass > p2->dMass ? -1 : (p1->dMass < p2->dMass ? 1 :
										  p1->iIdx < p2->iIdx ? -1 : 1));
	}

static int biggest(char *filename,int nBig)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA d;
	BIG *p;
	int i;


	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"biggest(): Unable to open \"%s\" for reading.\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"biggest(): Error reading second file header.\n");
		return -1;
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

	assert(h.n_data > 0);
	if (nBig > h.n_data) {
		(void) fprintf(stderr,"biggest(): Too few particles (%i) to find (n=%i)th biggest.\n",h.n_data,nBig);
		return -1;
		}

	p = (BIG *) malloc((size_t) h.n_data*sizeof(BIG));
	assert(p != NULL);

	for (i=0;i<h.n_data;i++) {
		if (ssioData(&ssio,&d)) {
			(void) fprintf(stderr,"biggest(): Error read second file data.\n");
			return -1;
			}
		p[i].iIdx = i;
		p[i].iOrgIdx = d.org_idx;
		p[i].dMass = d.mass;
		}

	(void) ssioClose(&ssio);

	/* sort in decreasing mass order */

	qsort((void *)p,(size_t) h.n_data,sizeof(BIG),compar);

	i = p[nBig - 1].iOrgIdx;

	free((void *)p);

	return i;
	}

static void usage(char *progname)
{
	(void) fprintf(stderr,"Usage: %s ( -i # | -m # ) ssg-map-file\n"
				   "where -i specifies particle original index\n"
				   "      -m specifies nth most massive particle\n",progname);
	exit(1);
	}

int main(int argc,char *argv[])
{
	FILE *fp;
	char ssfile1[MAXPATHLEN],ssfile2[MAXPATHLEN];
	int c,n,iOrgIdx=-1,nBig=-1;

	/* Disable stdout buffering */

	setbuf(stdout,(char *)NULL);

	/* Check arguments */

	while ((c = getopt(argc,argv,"i:m:")) != EOF)
		switch (c) {
		case 'i':
			if (nBig != -1) {
				(void) fprintf(stderr,"main(): Cannot specify both -i and -m.\n");
				usage(argv[0]);
				}
			n = atoi(optarg);
			if (n < 0) {
				(void) fprintf(stderr,"main(): Argument to -i cannot be negative.\n");
				usage(argv[0]);
				}
			iOrgIdx = n;
			break;
		case 'm':
			if (iOrgIdx != -1) {
				(void) fprintf(stderr,"main(): Cannot specify both -i and -m.\n");
				usage(argv[0]);
				}
			n = atoi(optarg);
			if (n <= 0) {
				(void) fprintf(stderr,"main(): Argument to -m must be positive.\n");
				usage(argv[0]);
				}
			nBig = n;
			break;
		default:
			usage(argv[0]);
			}

	if (iOrgIdx == -1 && nBig == -1) usage(argv[0]);

	if (optind != argc - 1) usage(argv[0]);

	/* Open files */

	fp = fopen(argv[optind],"r");
	assert(fp != NULL);

	assert(MAXPATHLEN >= 256);
	(void) fscanf(fp,"%255s",ssfile1);
	ssfile1[255] = '\0';
	(void) fscanf(fp,"%255s",ssfile2);
	ssfile2[255] = '\0';

	/* Determine original index of nth biggest particle in second ss file */

	if (iOrgIdx == -1)
		iOrgIdx = biggest(ssfile2,nBig);

	assert(iOrgIdx >= 0);

	/* Trace particle */

	(void) trace(fp,ssfile1,iOrgIdx);

	/* Close map file */

	(void) fclose(fp);

	/* All done */

	return 0;
	}

/* sst.c */
