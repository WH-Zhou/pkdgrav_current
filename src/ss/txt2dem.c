/*
 ** txt2dem.c -- SRS 8/13/11
 ** =========
 ** Converts human-readable text to binary DEM data file.
 ** Mods by DCR 8/7/16
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <assert.h>
#include <rpc/rpc.h>

#define DEM_EXT "dem"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void process_scalar(FILE *fpi,XDR *xdrs,int n)
{
	double d;
	int i;

	for (i=0;i<n;i++) {
		fscanf(fpi,"%lf",&d);
		xdr_double(xdrs,&d);
		}
	}

static void process_vector(FILE *fpi,XDR *xdrs,int n)
{
	process_scalar(fpi,xdrs,3*n);
	}

static void convert(char *infile)
{
	XDR xdrs;
	FILE *fpi,*fpo;
	char outfile[MAXPATHLEN];
	double dTime;
	int i,j,nWarn,nErr,nPart,nMaxOvlp,bRotationDashpot,iOrder;
	int nMaxOvlpWalls,iWallID; /* walls only */
#if defined(__APPLE__) && defined (__LP64__)
	int liOverlapCounter;
#else
	long int liOverlapCounter;
#endif

	nWarn = nErr = 0;
	printf("%s\n",infile);
	snprintf(outfile,MAXPATHLEN,"%s.%s",infile,DEM_EXT);
	assert(strcmp(infile,outfile));
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	xdrstdio_create(&xdrs,fpo,XDR_ENCODE);
	fscanf(fpi,"%lf",&dTime);
	if (dTime < 0.0) {
		fprintf(stderr,"WARNING: time < 0\n");
		++nWarn;
		}
	fscanf(fpi,"%i",&nPart);
	printf("Number of particles: %i\n",nPart);
	if (nPart < 1) {
		fprintf(stderr,"ERROR: nData < 1\n");
		++nErr;
		}
	fscanf(fpi,"%i",&nMaxOvlp);
	printf("MAX_NUM_OVERLAPS_PER_PARTICLE: %i\n",nMaxOvlp);
	if (nMaxOvlp < 0) {
		fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE < 0\n");
		++nErr;
		}
	fscanf(fpi,"%i",&nMaxOvlpWalls);
	printf("MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS: %i\n",nMaxOvlpWalls);
	if (nMaxOvlpWalls < -1) {
		fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS < -1\n");
		++nErr;
		}
	fscanf(fpi,"%i",&bRotationDashpot);
	printf("Rotation dashpot flag: %i\n",bRotationDashpot);
	if (bRotationDashpot != 0 && bRotationDashpot != 1) {
		fprintf(stderr,"ERROR: Rotation dashpot flag must be 0 or 1\n");
		++nErr;
		}
	if (nWarn > 0)
		fprintf(stderr,"Number of warnings = %i\n",nWarn);
	if (nErr > 0)
		fprintf(stderr,"Number of errors = %i\n",nErr);
	assert(nErr == 0);
	xdr_double(&xdrs,&dTime);
	xdr_int(&xdrs,&nPart);
	xdr_int(&xdrs,&nMaxOvlp);
	xdr_int(&xdrs,&nMaxOvlpWalls);
	xdr_int(&xdrs,&bRotationDashpot);
	for (i=0;i<nPart;i++) {
		process_vector(fpi,&xdrs,2); /* vPred & wPred */
		for (j=0;j<nMaxOvlp;j++) {
			fscanf(fpi,"%i",&iOrder);
			if (iOrder >= nPart) {
				/* negative is ok ==> no overlap */
				fprintf(stderr,"ERROR: Invalid iOrder value (%i), particle %i, overlap %i\n",iOrder,i,j);
				exit(1);
				}
			xdr_int(&xdrs,&iOrder);
			process_vector(fpi,&xdrs,2); /* vShear & vnOld */
			if (bRotationDashpot) {
				process_vector(fpi,&xdrs,1); /* vRoll */
				process_scalar(fpi,&xdrs,1); /* dTwist */
				}
#if defined(__APPLE__) && defined (__LP64__)
			fscanf(fpi,"%i",&liOverlapCounter);
#else
			fscanf(fpi,"%li",&liOverlapCounter);
#endif
			if (iOrder >= 0 && liOverlapCounter <= 0) {
				fprintf(stderr,"ERROR: Invalid overlap counter (%li), particle %i (ID %i), overlap %i\n",(long) liOverlapCounter,i,iOrder,j);
				exit(1);
				}
			xdr_long(&xdrs,&liOverlapCounter);
			}
		for (j=0;j<nMaxOvlpWalls;j++) {
			fscanf(fpi,"%i",&iWallID);
			if (iWallID >= nPart) {
				/* negative is ok ==> no overlap */
				fprintf(stderr,"ERROR: Invalid wall ID (%i), particle %i (ID %i), overlap %i\n",iWallID,i,iWallID,j);
				exit(1);
				}
			xdr_int(&xdrs,&iWallID);
			process_vector(fpi,&xdrs,2); /* vShear & vnOld */
			if (bRotationDashpot) {
				process_vector(fpi,&xdrs,1); /* vRoll */
				process_scalar(fpi,&xdrs,1); /* dTwist */
				}
#if defined(__APPLE__) && defined (__LP64__)
			fscanf(fpi,"%i",&liOverlapCounter);
#else
			fscanf(fpi,"%li",&liOverlapCounter);
#endif
			if (iWallID >= 0 && liOverlapCounter <= 0) {
				fprintf(stderr,"ERROR: Invalid overlap counter (%li), particle %i, wall %i (ID %i)\n",(long) liOverlapCounter,i,j,iWallID);
				exit(1);
				}
			xdr_long(&xdrs,&liOverlapCounter);
			}
		}
	xdr_destroy(&xdrs);
	fclose(fpo);
	fclose(fpi);
	}

int main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		fprintf(stderr,"Usage: %s file [ file ... ]\n",argv[0]);
		return 1;
		}
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* txt2dem.c */
