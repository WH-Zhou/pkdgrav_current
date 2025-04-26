/*
 ** ori2ort.c -- DCR 6/30/17
 ** =========
 ** Converts particle orientation binary data to ASCII text.
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <vector.h>

#define ORT_EXT ".ort"

static void convert(char *infile, int bQuiet)
{
	SSIO ssio;
	ORIHEAD head;
	FILE *fp;
	char outfile[MAXPATHLEN];
	int i;

	if (!bQuiet) printf("%s: ",infile);
	if (ssioNewExt(infile,ORI_EXT,outfile,ORT_EXT)) {
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
	if (ssioOriHead(&ssio,&head)) {
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
		ORIDATA data;
		for (i=0;i<head.n_data;i++) {
			if (ssioOriData(&ssio,&data)) {
				fprintf(stderr,"Corrupt data.\n");
				goto finish;
				}
			fprintf(fp,"%.16e %.16e %.16e %.16e %.16e %.16e\n",
					data.p1[X],data.p1[Y],data.p1[Z],
					data.p2[X],data.p2[Y],data.p2[Z]);
			}
		break;
		}
	case SSIO_MAGIC_REDUCED: {
		ORIRDATA data;
		for (i=0;i<head.n_data;i++) {
			if (ssioOriDataReduced(&ssio,&data)) {
				fprintf(stderr,"Corrupt data.\n");
				goto finish;
				}
			fprintf(fp,"%e %e %e %e %e %e\n",
					data.p1[X],data.p1[Y],data.p1[Z],
					data.p2[X],data.p2[Y],data.p2[Z]);
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

/* ori2ort.c */
