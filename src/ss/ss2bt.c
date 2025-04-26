/*
 ** ss2bt.c -- DCR 97-03-04
 ** =======
 ** Converts Solar System binary data to box_tree ASCII data.
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <vector.h>

#define BT_EXT ".bt"

static void convert(char *infile, int bQuiet)
{
	SSIO ssio;
	SSHEAD head;
	FILE *fp;
	char outfile[MAXPATHLEN];
	int i;

	if (!bQuiet) printf("%s: ",infile);
	if (ssioNewExt(infile,SS_EXT,outfile,BT_EXT)) {
		fprintf(stderr,"Unable to generate output filename.\n");
		return;
		}
	if (ssioOpen(infile,&ssio,SSIO_READ)) {
		fprintf(stderr,"Unable to open \"%s\".\n",infile);
		return;
		}
	if (!(fp = fopen(outfile,"w"))) {
		fprintf(stderr,"Unable to open \"%s\".\n",outfile);
		ssioClose(&ssio);
		return;
		}
	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Corrupt header.\n");
		goto finish;
		}
	if (!bQuiet) printf("time = %e n_data = %i",head.time,head.n_data);
	if (head.iMagicNumber == SSIO_MAGIC_REDUCED && !bQuiet)
		printf(" (reduced format)");
	if (!bQuiet) printf("\n");
	if (head.n_data <= 0) {
		fprintf(stderr,"Invalid input data format.\n");
		goto finish;
		}
	switch(head.iMagicNumber) {
	case SSIO_MAGIC_STANDARD: {
		SSDATA data;
		for (i=0;i<head.n_data;i++) {
			if (ssioData(&ssio,&data)) {
				fprintf(stderr,"Corrupt data.\n");
				goto finish;
				}
			fprintf(fp,"%i %i %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
					"%.16e %.16e %.16e %.16e %i\n",i,data.org_idx,data.mass,
					data.radius,data.pos[X],data.pos[Y],data.pos[Z],
					data.vel[X],data.vel[Y],data.vel[Z],data.spin[X],
					data.spin[Y],data.spin[Z],data.color);
			}
		break;
		}
	case SSIO_MAGIC_REDUCED: {
		SSRDATA data;
		for (i=0;i<head.n_data;i++) {
			if (ssioDataReduced(&ssio,&data)) {
				fprintf(stderr,"Corrupt data.\n");
				goto finish;
				}
			fprintf(fp,"%i %i %e %e %e %e %e 0 0 0 0 0 0 %i\n",
					i,data.iOrgIdx,data.fMass,data.fRadius,
					data.vPos[X],data.vPos[Y],data.vPos[Z],data.iColor);
			}
		break;
		}
	default:
		fprintf(stderr,"Unrecognized ss file magic number (%i).\n",head.iMagicNumber);
		}
 finish:
	fclose(fp);
	ssioClose(&ssio);
	}

int main(int argc,char *argv[])
{
	int bQuiet = 0; /* default = verbose */

	int i;

	setbuf(stdout,(char *) NULL);
	while ((i = getopt(argc, argv, "q")) != EOF)
		switch (i) {
		case 'q':
			bQuiet = 1;
			}
	if (optind >= argc) {
		fprintf(stderr,"Usage: %s [ -q ] file [ file ... ]\n",argv[0]);
		fprintf(stderr,"Where: -q = activate quiet mode\n");
		exit(1);
		}
	for (i=optind;i<argc;++i)
		convert(argv[i], bQuiet);
	return 0;
	}

/* ss2bt.c */
