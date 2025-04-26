/*
 ** spra.c -- SRS 3/17/10
 ** =========
 ** Groups contiguous particles using springs data file.
 ** Requires input file that with "ss." and the presence of springs file "inputfile".spr
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <rpc/rpc.h>
#include <springs.h>
#include <ssdefs.h>
#include <boolean.h>

#define SPR_EXT ".spr"
#define GROUP_EXT ".group"

#define MAX_STR_LEN 256

typedef struct {
        int iGroup;
        int iOrder1;
        int iOrder2;
        float fZeroStrainLength;
        float fYoungsModulus;
        float fStressLimit;
} spring;

static long int get_num_particles(char *infile) { /* returns maximum possible number of springs in system */
	XDR xdrs;
	FILE *fpi;
	double dTime;
	int nPart;
	fpi = fopen(infile,"r");
	assert(fpi != NULL);
	xdrstdio_create(&xdrs,fpi,XDR_DECODE);
	xdr_double(&xdrs,&dTime);
	if (dTime < 0.0) {
		(void) fprintf(stderr,"%s contains invalid data (dTime = %g < 0.0).\n",infile,dTime);
		exit(1);
		}
	xdr_int(&xdrs,&nPart);
	if (nPart < 1) {
		(void) fprintf(stderr,"%s contains invalid data (nData = %i < 1).\n",infile,nPart);
		exit(1);
		}
	xdr_destroy(&xdrs);
	(void) fclose(fpi);
	return nPart;
}

static int group(char *infile,spring *s,BOOLEAN *frieze) {
	XDR xdrs;
	FILE *fpi;
	double dTime;
	float fZeroStrainLength,fYoungsModulus,fStressLimit;
	int i,j,nSprings,nPart,iOrder,nMaxSpringsInSystem,iSystemSpring,iThisGroup,bDidSomething,bAllGrouped;
	char spr_infile[MAX_STR_LEN];
	snprintf(spr_infile,strlen(infile)+strlen(SPR_EXT),"%s%s",infile,SPR_EXT);
	fpi = fopen(spr_infile,"r");
	assert(fpi != NULL);
	xdrstdio_create(&xdrs,fpi,XDR_DECODE);
	xdr_double(&xdrs,&dTime);
	if (dTime < 0.0) {
		(void) fprintf(stderr,"%s contains invalid data (dTime = %g < 0.0).\n",infile,dTime);
		exit(1);
		}
	xdr_int(&xdrs,&nPart);
	if (nPart < 1) {
		(void) fprintf(stderr,"%s contains invalid data (nData = %i < 1).\n",infile,nPart);
		exit(1);
		}
	xdr_int(&xdrs,&i); /* dummy read */
	nMaxSpringsInSystem = nPart*MAX_NUM_SPRINGS_PER_PARTICLE;
	for (i=0;i<nMaxSpringsInSystem;i++)
		for (j=0;j<MAX_NUM_SPRINGS_PER_PARTICLE;j++) {
			s[i].iGroup = 1;
			s[i].iOrder1 = -1;
			s[i].iOrder2 = -1;
			s[i].fZeroStrainLength = -1.;
			s[i].fYoungsModulus = -1.;
			s[i].fStressLimit = -1.;
		}
	for (i=0;i<nPart;i++) {
	        xdr_int(&xdrs,&nSprings);
		for (j=0;j<nSprings;j++) {
		        iSystemSpring = i*MAX_NUM_SPRINGS_PER_PARTICLE+j;
			xdr_int(&xdrs,&iOrder);
			if (iOrder >= nPart) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: iOrder = %i >= %i).\n",infile,i,j,iOrder,nPart);
				exit(1);
			}
			xdr_float(&xdrs,&fZeroStrainLength);
			if (fZeroStrainLength < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fZeroStrainLength = %g < 0).\n",infile,i,j,fZeroStrainLength);
				exit(1);
			}
			xdr_float(&xdrs,&fYoungsModulus);
			if (fYoungsModulus < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fYoungsModulus = %g < 0).\n",infile,i,j,fYoungsModulus);
				exit(1);
			}
			xdr_float(&xdrs,&fStressLimit);
			if (fStressLimit < 0.) {
				(void) fprintf(stderr,"%s contains invalid data (particle %i spring %i: fStressLimit = %g < 0).\n",infile,i,j,fStressLimit);
				exit(1);
			}
			s[iSystemSpring].iGroup = 0;
			s[iSystemSpring].iOrder1 = i;
			s[iSystemSpring].iOrder2 = iOrder;
			s[iSystemSpring].fZeroStrainLength = fZeroStrainLength;
			s[iSystemSpring].fYoungsModulus = fYoungsModulus;
			s[iSystemSpring].fStressLimit = fStressLimit;
		}
	}
	xdr_destroy(&xdrs);
	(void) fclose(fpi);

	/* make groups */
	iThisGroup = 0;
	bAllGrouped = 0;
	while (!bAllGrouped) {
	  bAllGrouped = 1;
	  for (i=0;i<nMaxSpringsInSystem;i++)
	    if (s[i].iGroup < iThisGroup)
	      iThisGroup = s[i].iGroup;
	  iThisGroup--;
	  for (i=0;i<nMaxSpringsInSystem;i++)
	    if (s[i].iGroup == 0) {
	      s[i].iGroup = iThisGroup;
	      break;
	    }
	  bDidSomething = 1;
	  while (bDidSomething) {
	    bDidSomething = 0;
	    for (i=0;i<nMaxSpringsInSystem;i++)
	      if (s[i].iGroup == iThisGroup)
		for (j=0;j<nMaxSpringsInSystem;j++) {
		  if ((s[i].iOrder1 == s[j].iOrder1 || s[i].iOrder1 == s[j].iOrder2 || s[i].iOrder2 == s[j].iOrder1 || s[i].iOrder2 == s[j].iOrder2) && (s[j].iGroup != iThisGroup)) {
		    s[j].iGroup = iThisGroup;
		    bDidSomething = 1;
		  }
		}
	    for (i=0;i<nMaxSpringsInSystem;i++)
	      if (s[i].iGroup == 0) {
		bAllGrouped = 0;
		break;
	      }
	    break;
	  }
        }

	/* check to see if groups need to be merged */
	bDidSomething = 1;
	while (bDidSomething) {
	  bDidSomething = 0;	
	  for (i=0;i<nMaxSpringsInSystem;i++)
	    for (j=0;j<nMaxSpringsInSystem;j++)
	      if ((s[i].iOrder1 == s[j].iOrder1 || s[i].iOrder1 == s[j].iOrder2 || s[i].iOrder2 == s[j].iOrder1 || s[i].iOrder2 == s[j].iOrder2) && (s[i].iGroup != s[j].iGroup)) {
		if (s[i].iGroup > s[j].iGroup) s[j].iGroup = s[i].iGroup;
		if (s[j].iGroup > s[i].iGroup) s[i].iGroup = s[j].iGroup;
		bDidSomething = 1;
	      }
	}

	/* fill list of free particles */
	for (i=0;i<nPart;i++)
	  for (j=0;j<nMaxSpringsInSystem;j++)
	    if (i == s[j].iOrder1 || i == s[j].iOrder2) {
	      frieze[i] = FALSE;
	      break;
	    }
	    else
	      frieze[i] = TRUE;

	/* find most negative group number */
     	for (i=0;i<nMaxSpringsInSystem;i++)
	  if (s[i].iGroup < iThisGroup) iThisGroup = s[i].iGroup;
	return(0 - iThisGroup); /* could make this a bitwise op! */
}

static void springs_analyze(char *infile,spring *s,int iGroup,BOOLEAN *frieze) {
        SSIO ssio_in,ssio_out;
	SSHEAD head;
	SSDATA data;
	char outfile[MAXPATHLEN];
	int i,j;

	if (ssioNewExt(infile,SS_EXT,outfile,GROUP_EXT)) {
		fprintf(stderr,"Unable to generate output filename.\n");
		return;
	}
	if (ssioOpen(infile,&ssio_in,SSIO_READ)) {
		fprintf(stderr,"Unable to open \"%s\".\n",infile);
		return;
	}
	if (ssioOpen(outfile,&ssio_out,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\".\n",outfile);
		return;
	}
	if (ssioHead(&ssio_in,&head)) {
		fprintf(stderr,"Corrupt header.\n");
		goto finish;
	}
	if (head.iMagicNumber == SSIO_MAGIC_REDUCED) {
		fprintf(stderr,"reduced file format not supported.\n");
		goto finish;
	}
	if (head.iMagicNumber != SSIO_MAGIC_STANDARD) {
		fprintf(stderr,"Unrecognized ss file magic number (%i).\n",head.iMagicNumber);
		goto finish;
	}
	if (head.n_data <= 0) {
		fprintf(stderr,"Invalid input data format.\n");
		goto finish;
		}
	if (ssioHead(&ssio_out,&head)) {
	        fprintf(stderr,"Error writing header.\n");
		goto finish;
	}
	for (i=0;i<head.n_data;i++) {
	        if (ssioData(&ssio_in,&data)) {
		        fprintf(stderr,"Corrupt input data.\n");
			goto finish;
		}
		if (frieze[i] == TRUE) {
		        iGroup--;
		        data.org_idx = iGroup;
		}
		else
		        for (j=0;j<head.n_data*MAX_NUM_SPRINGS_PER_PARTICLE;j++)
			        if (i == s[j].iOrder1 || i == s[j].iOrder2) {
				        data.org_idx = s[j].iGroup;
					break;
				}
		if (data.org_idx >= 0) {
		        fprintf(stderr,"Error finding group (%d).\n",i);
			goto finish;
		}
		if (ssioData(&ssio_out,&data)) {
		        fprintf(stderr,"Unable to write data for particle %i\n",i);
			goto finish;
		}
	}
finish:
	ssioClose(&ssio_in);
	ssioClose(&ssio_out);
}

int main(int argc,char *argv[]) {
        BOOLEAN *frieze; /* array elements will be deemed TRUE for liberated particles */
        int i,nGroups,nPart;
	spring *s;

	setbuf(stdout,(char *)NULL);
	if (argc != 2) {
	        (void) fprintf(stderr,"%s takes one argument...\nsyntax: %s [file]\n",argv[0],argv[0]);
		return 1;
	}
	nPart = (int) get_num_particles(argv[1]);
	s = (spring *) malloc((size_t) nPart*MAX_NUM_SPRINGS_PER_PARTICLE*sizeof(spring));
	frieze = (BOOLEAN *) malloc((size_t) nPart*sizeof(BOOLEAN));
	nGroups = group(argv[1],s,frieze);
	springs_analyze(argv[1],s,-nGroups,frieze);

	/* output springs data to screen */
	for (i=0;i<(long int)nPart * MAX_NUM_SPRINGS_PER_PARTICLE;i++)
	        if (s[i].iGroup != 1)
		        printf("%d %d %d %.7e %.7e %.7e\n",s[i].iGroup,s[i].iOrder1,s[i].iOrder2,s[i].fZeroStrainLength,s[i].fYoungsModulus,s[i].fStressLimit);
	return 0;
}

/* spra.c */
