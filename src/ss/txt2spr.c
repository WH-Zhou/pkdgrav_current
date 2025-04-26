/*
 ** txt2spr.c -- DCR 6/11/08
 ** =========
 ** Converts human-readable text to binary springs data file.
 ** Revamped for new springs file format -- SRS 8/4/09
 ** Mods by DCR 8/7/16
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
#include <springs.h>

#define SPR_EXT "spr"

#ifndef MAXPATHLEN
#define MAXPATHLEN 256
#endif

static void convert(char *infile)
{
	XDR xdrs;
	FILE *fpi,*fpo;
	char outfile[MAXPATHLEN];
	double dTime,vPred[3];
	float fZeroStrainLength,fYoungsModulus,fStressLimit;
	int i,j,nPart,nMaxSpr,iOrder;

	printf("%s\n",infile);
	snprintf(outfile,MAXPATHLEN,"%s.%s",infile,SPR_EXT);
	assert(strcmp(infile,outfile));
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	fpo = fopen(outfile,"w");
	assert(fpo != NULL);
	xdrstdio_create(&xdrs,fpo,XDR_ENCODE);
	fscanf(fpi,"%lf",&dTime);
	if (dTime < 0.0) {
		fprintf(stderr,"%s contains invalid data (dTime = %g < 0.0).\n",infile,dTime);
		exit(1);
		}
	fscanf(fpi,"%i",&nPart);
	if (nPart < 1) {
		fprintf(stderr,"%s contains invalid data (nData = %i < 1).\n",infile,nPart);
		exit(1);
		}
	fscanf(fpi,"%i",&nMaxSpr);
	if (nMaxSpr < 0) {
		fprintf(stderr,"%s contains invalid data (nMaxSpr = %i < 0).\n",infile,nMaxSpr);
		exit(1);
		}
	xdr_double(&xdrs,&dTime);
	xdr_int(&xdrs,&nPart);
	xdr_int(&xdrs,&nMaxSpr);
	for (i=0;i<nPart;i++) {
		for (j=0;j<3;j++)
			fscanf(fpi,"%lf",&vPred[j]);
		for (j=0;j<3;j++)
			xdr_double(&xdrs,&vPred[j]);
		for (j=0;j<nMaxSpr;j++) {
			fscanf(fpi,"%i",&iOrder);
			if (iOrder >= nPart) {
				fprintf(stderr,"%s contains invalid data (particle %i spring %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
				exit(1);
			}
			xdr_int(&xdrs,&iOrder);
			fscanf(fpi,"%f",&fZeroStrainLength);
			if (fZeroStrainLength < 0.0) {
				fprintf(stderr,"%s contains invalid data (particle %i spring %i: fZeroStrainLength = %g < 0).\n",infile,i,j,fZeroStrainLength);
				exit(1);
			}
			xdr_float(&xdrs,&fZeroStrainLength);
			fscanf(fpi,"%f",&fYoungsModulus);
			if (fYoungsModulus < 0.0) {
				fprintf(stderr,"%s contains invalid data (particle %i spring %i: fYoungsModulus = %g < 0).\n",infile,i,j,fYoungsModulus);
				exit(1);
			}
			xdr_float(&xdrs,&fYoungsModulus);
			fscanf(fpi,"%f",&fStressLimit);
			if (fStressLimit < 0.0) {
				fprintf(stderr,"%s contains invalid data (particle %i spring %i: fStressLimit = %g < 0).\n",infile,i,j,fStressLimit);
				exit(1);
			}
			xdr_float(&xdrs,&fStressLimit);
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

/* txt2spr.c */
