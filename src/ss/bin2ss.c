/*
 ** bin2ss.c -- DCR 98-02-11
 ** ========
 ** Converts tipsy binary data to Solar System format.
 */

#include <ssdefs.h>
#include <tipsydefs.h>
#include <vector.h>

#define BIN_EXT ".bin"

static void
convert(char *infile)
{
	FILE *fp;
	struct dump hi;
	struct dark_particle dp;
	SSIO ssio;
	SSHEAD ho;
	SSDATA data;
	char outfile[MAXPATHLEN];
	int i;

	(void) printf("%s: ",infile);
	if (ssioNewExt(infile,BIN_EXT,outfile,SS_EXT)) {
		(void) fprintf(stderr,"Unable to generate output filename.\n");
		return;
		}
	if (!(fp = fopen(infile,"r"))) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",infile);
		return;
		}
	if (ssioOpen(outfile,&ssio,SSIO_WRITE)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",outfile);
		(void) fclose(fp);
		return;
		}
	if (fread(&hi,sizeof(struct dump),1,fp) != 1) {
		(void) fprintf(stderr,"Corrupt header.\n");
		goto finish;
		}
	(void) printf("time = %e, n_data = %i\n",hi.time,hi.nbodies);
	if (hi.ndim != N_DIM) {
		(void) fprintf(stderr,"Unexpected input file dimension.\n");
		goto finish;
		}
	if (hi.ndark != hi.nbodies || hi.nbodies <= 0) {
		(void) fprintf(stderr,"Invalid input file format.\n");
		goto finish;
		}
	ho.time = hi.time;
	ho.n_data = hi.nbodies;
	ho.iMagicNumber = SSIO_MAGIC_STANDARD;
	if (ssioHead(&ssio,&ho)) {
		(void) fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	for (i=0;i<hi.nbodies;i++) {
		if (fread(&dp,sizeof(struct dark_particle),1,fp) != 1) {
			(void) fprintf(stderr,"Corrupt data.\n");
			goto finish;
			}
		data.mass = dp.mass;
		data.radius = dp.eps;
		COPY_VEC(dp.pos,data.pos);
		COPY_VEC(dp.vel,data.vel);
		ZERO_VEC(data.spin);
		/* make a guess at Sun/planets/planetesimal color tags */
		data.color = PLANETESIMAL;
		if (data.mass == 0.0) data.color = TEST;
		if (data.mass > 0.5 && data.mass < 1.5) data.color = SUN;
		if (data.mass > 9.0e-4 && data.mass < 1.0e-3) data.color = JUPITER;
		if (data.mass > 2.0e-4 && data.mass < 3.0e-4) data.color = SATURN;
		if (data.mass > 4.0e-5 && data.mass < 5.0e-5) data.color = URANUS;
		if (data.mass > 5.0e-5 && data.mass < 6.0e-5) data.color = NEPTUNE;
		data.org_idx = i;
		if (ssioData(&ssio,&data)) {
			(void) fprintf(stderr,"Error writing data.\n");
			goto finish;
			}
		}
 finish:
	(void) ssioClose(&ssio);
	(void) fclose(fp);
	}

int
main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		(void) fprintf(stderr,"Usage: %s bin-file [ bin-file ... ]\n",argv[0]);
		return 1;
		}
	for (i=1;i<argc;i++)
		convert(argv[i]);
	return 0;
	}

/* bin2ss.c */
