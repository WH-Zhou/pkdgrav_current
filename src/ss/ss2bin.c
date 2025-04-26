/*
 ** ss2bin.c -- DCR 98-01-13
 ** ========
 ** Converts Solar System binary data to tipsy format.
 **
 */

#include <ssdefs.h>
#include <tipsydefs.h>
#include <vector.h>

#define BIN_EXT ".bin"

/*
 ** recent (00-01-11) addition: "red" particles --> gas,
 ** "blue" particles --> stars, other --> dark.
 **
 */

/*#define USE_GAS_AND_STARS*/

#ifdef USE_GAS_AND_STARS
#include <colors.h>
#endif

static int
get_particle_types(SSIO *ssio, int n, int *nsph, int *ndark, int *nstar)
{
#ifdef USE_GAS_AND_STARS
	SSDATA data;
	int i;

	if (!ssio || !nsph || !ndark || !nstar || n < 0) {
		(void) fprintf(stderr,"get_particle_types(): Invalid argument(s).\n");
		return 1;
		}

	(void) printf("first pass -- counting particle types\n");

	*nsph = *ndark = *nstar = 0;

	for (i=0;i<n;i++) {
		if (ssioData(ssio,&data)) {
			(void) fprintf(stderr,"Corrupt data.\n");
			return 1;
			}
		switch (data.color) {
		case RED:
			++(*nsph);
			break;
		case BLUE:
			++(*nstar);
			break;
		default:
			++(*ndark);
			}
		}

	ssioSetPos(ssio,SSHEAD_SIZE); /* rewind to start of data */

	(void) printf("nsph = %i ndark = %i nstar = %i total = %i\n",
				  *nsph,*ndark,*nstar,*nsph + *ndark + *nstar);
#else
	(void) ssio; /* to suppress unused parameter warning */
	*nsph = *nstar = 0;
	*ndark = n;
#endif

	return 0;
	}

static void
convert(char *infile)
{
	SSIO ssio;
	SSHEAD hi;
	SSDATA data;
	FILE *fp;
	struct dump ho;
	struct dark_particle dp;
	char outfile[MAXPATHLEN];
	int i,nsph,ndark,nstar;

#ifdef USE_GAS_AND_STARS
	struct gas_particle gp;
	struct star_particle sp;
	int isph=0,idark=0,istar=0;
#endif

	(void) printf("%s: ",infile);
	if (ssioNewExt(infile,SS_EXT,outfile,BIN_EXT)) {
		(void) fprintf(stderr,"Unable to generate output filename.\n");
		return;
		}
	if (ssioOpen(infile,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\" for reading.\n",infile);
		return;
		}
	if (!(fp = fopen(outfile,"w"))) {
		(void) fprintf(stderr,"Unable to open \"%s\" for writing.\n",outfile);
		(void) ssioClose(&ssio);
		return;
		}
	(void) fclose(fp); /* empty file */
	if (!(fp = fopen(outfile,"r+"))) {
		(void) fprintf(stderr,"Unable to open \"%s\" for updating.\n",outfile);
		(void) ssioClose(&ssio);
		return;
		}
	if (ssioHead(&ssio,&hi)) {
		(void) fprintf(stderr,"Corrupt header.\n");
		goto finish;
		}
	(void) printf("time = %e, n_data = %i\n",hi.time,hi.n_data);
	if (hi.n_data <= 0) {
		(void) fprintf(stderr,"Invalid input data format.\n");
		goto finish;
		}
	switch(hi.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		goto finish;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",hi.iMagicNumber);
		goto finish;
		}
	if (get_particle_types(&ssio,hi.n_data,&nsph,&ndark,&nstar))
		goto finish;
	if (nsph + ndark + nstar != hi.n_data) {
		(void) fprintf(stderr,"Inconsistent data.\n");
		goto finish;
		}
	ho.time = hi.time;
	ho.nbodies = hi.n_data;
	ho.ndim = N_DIM;
	ho.nsph = nsph;
	ho.ndark = ndark;
	ho.nstar = nstar;
	if (fwrite(&ho,sizeof(ho),1,fp) != 1) {
		(void) fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	for (i=0;i<hi.n_data;i++) {
		if (ssioData(&ssio,&data)) {
			(void) fprintf(stderr,"Corrupt data.\n");
			goto finish;
			}
#ifdef USE_GAS_AND_STARS
		switch (data.color) {
		case RED:
			if (fseek(fp,sizeof(ho) + (isph++)*sizeof(gp),SEEK_SET)) {
				(void) fprintf(stderr,"Seek error.\n");
				goto finish;
				}
			gp.mass = data.mass;
			COPY_VEC(data.pos,gp.pos);
			COPY_VEC(data.vel,gp.vel);
			gp.rho = 1;
			gp.temp = 0;
			gp.hsmooth = data.radius;
			gp.metals = 0;
			gp.phi = 0;
			if (fwrite(&gp,sizeof(gp),1,fp) != 1) {
				(void) fprintf(stderr,"Error writing data.\n");
				goto finish;
				}
			break;
		case BLUE:
			if (fseek(fp,sizeof(ho) + nsph*sizeof(gp) + ndark*sizeof(dp) +
					  (istar++)*sizeof(sp),SEEK_SET)) {
				(void) fprintf(stderr,"Seek error.\n");
				goto finish;
				}
			sp.mass = data.mass;
			COPY_VEC(data.pos,sp.pos);
			COPY_VEC(data.vel,sp.vel);
			sp.metals = 0;
			sp.tform = 0;
			sp.eps = data.radius;
			sp.phi = 0;
			if (fwrite(&sp,sizeof(sp),1,fp) != 1) {
				(void) fprintf(stderr,"Error writing data.\n");
				goto finish;
				}
			break;
		default:
			if (fseek(fp,sizeof(ho) + nsph*sizeof(gp) + (idark++)*sizeof(dp),
					  SEEK_SET)) {
				(void) fprintf(stderr,"Seek error.\n");
				goto finish;
				}
#endif
			dp.mass = (Real) data.mass;
			dp.pos[X] = (Real) data.pos[X];
			dp.pos[Y] = (Real) data.pos[Y];
			dp.pos[Z] = (Real) data.pos[Z];
			dp.vel[X] = (Real) data.vel[X];
			dp.vel[Y] = (Real) data.vel[Y];
			dp.vel[Z] = (Real) data.vel[Z];
			dp.eps = (Real) data.radius;
			dp.phi = 0;
			if (fwrite(&dp,sizeof(dp),1,fp) != 1) {
				(void) fprintf(stderr,"Error writing data.\n");
				goto finish;
				}
			}
#ifdef USE_GAS_AND_STARS
		}
	if (isph != nsph || idark != ndark || istar != nstar) {
		(void) fprintf(stderr,"Inconsistent data.\n");
		goto finish;
		}
#endif
 finish:
	(void) fclose(fp);
	(void) ssioClose(&ssio);
	}

int
main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		(void) fprintf(stderr,"Usage: %s ss-file [ ss-file ... ]\n",argv[0]);
		exit(1);
		}
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* ss2bin.c */
