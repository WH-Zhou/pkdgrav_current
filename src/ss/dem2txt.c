/*
** dem2txt.c -- SRS 7/11/11
** =========
** Converts binary DEM data file to human-readable text.
** Mods by DCR 8/7/16.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <assert.h>
#include <rpc/rpc.h>

#define TXT_EXT "txt"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void process_scalar(XDR *xdrs,FILE *fpo,int n)
{
	double d;
	int i;

	for (i=0;i<n;i++) {
		xdr_double(xdrs,&d);
		fprintf(fpo,"%.16e ",d);
		}
	}

static void process_vector(XDR *xdrs,FILE *fpo,int n)
{
	process_scalar(xdrs,fpo,3*n);
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
	snprintf(outfile,MAXPATHLEN,"%s.%s",infile,TXT_EXT);
	assert(strcmp(infile,outfile));
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	xdrstdio_create(&xdrs,fpi,XDR_DECODE);
	xdr_double(&xdrs,&dTime);
	printf("Time: %.16e\n",dTime);
	if (dTime < 0.0) {
		fprintf(stderr,"WARNING: time < 0\n");
		++nWarn;
		}
	xdr_int(&xdrs,&nPart);
	printf("Number of particles: %i\n",nPart);
	if (nPart < 1) {
		fprintf(stderr,"ERROR: nData < 1\n");
		++nErr;
		}
	xdr_int(&xdrs,&nMaxOvlp);
	printf("MAX_NUM_OVERLAPS_PER_PARTICLE: %i\n",nMaxOvlp);
	if (nMaxOvlp < 0) {
		fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE < 0\n");
		++nErr;
		}
	xdr_int(&xdrs,&nMaxOvlpWalls);
	printf("MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS: %i\n",nMaxOvlpWalls);
	if (nMaxOvlpWalls < -1) {
		fprintf(stderr,"ERROR: MAX_NUM_OVERLAPS_PER_PARTICLE_FOR_WALLS < -1\n");
		++nErr;
		}
	xdr_int(&xdrs,&bRotationDashpot);
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
	fprintf(fpo,"%.16e %i %i %i %i\n",dTime,nPart,nMaxOvlp,nMaxOvlpWalls,bRotationDashpot);
	for (i=0;i<nPart;i++) {
		process_vector(&xdrs,fpo,2); /* vPred & wPred */
		for (j=0;j<nMaxOvlp;j++) {
			xdr_int(&xdrs,&iOrder);
			if (iOrder >= nPart) {
				/* negative is ok ==> no overlap */
				fprintf(stderr,"ERROR: Invalid iOrder value (%i), particle %i, overlap %i\n",iOrder,i,j);
				exit(1);
				}
			fprintf(fpo,"%i ",iOrder);
			process_vector(&xdrs,fpo,2); /* vShear & vnOld */
			if (bRotationDashpot) {
				process_vector(&xdrs,fpo,1); /* vRoll */
				process_scalar(&xdrs,fpo,1); /* dTwist */
				}
			xdr_long(&xdrs,&liOverlapCounter);
			if (iOrder >= 0 && liOverlapCounter <= 0) {
				fprintf(stderr,"ERROR: Invalid overlap counter (%li), particle %i (ID %i), overlap %i\n",(long) liOverlapCounter,i,iOrder,j);
				exit(1);
				}
			fprintf(fpo,"%li ",(long) liOverlapCounter);
			}
		for (j=0;j<nMaxOvlpWalls;j++) {
			xdr_int(&xdrs,&iWallID);
			if (iWallID >= nPart) {
				/* negative is ok ==> no overlap */
				fprintf(stderr,"ERROR: Invalid wall ID (%i), particle %i (ID %i), overlap %i\n",iWallID,i,iWallID,j);
				exit(1);
				}
			fprintf(fpo,"%i ",iWallID);
			process_vector(&xdrs,fpo,2); /* vShear & vnOld */
			if (bRotationDashpot) {
				process_vector(&xdrs,fpo,1); /* vRoll */
				process_scalar(&xdrs,fpo,1); /* dTwist */
				}
			xdr_long(&xdrs,&liOverlapCounter);
			if (iWallID >= 0 && liOverlapCounter <= 0) {
				fprintf(stderr,"ERROR: Invalid overlap counter (%li), particle %i, wall %i (ID %i)\n",(long) liOverlapCounter,i,j,iWallID);
				exit(1);
				}
			fprintf(fpo,"%li ",(long) liOverlapCounter);
			}
		fprintf(fpo,"\n");
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
		fprintf(stderr,"Usage: %s demfile [ demfile ... ]\n",argv[0]);
		return 1;
		}
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* dem2txt.c */
