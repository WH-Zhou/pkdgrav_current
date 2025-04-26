/*
 ** ort2ori.c -- DCR 6/30/17
 ** =========
 ** Converts particle orientation ASCII text to binary data.
 **
 ** NOTE: time field in binary data is set to increment by 1.
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <vector.h>

#define ORT_EXT ".ort"

static void convert(char *infile, int bQuiet)
{
	static double time = 0.0; /* start at zero */

	FILE *fp;
	SSIO ssio;
	ORIHEAD head;
	ORIDATA data;
	char outfile[MAXPATHLEN];
	int rv;

	if (!bQuiet) printf("%s: ",infile);
	if (ssioNewExt(infile,ORT_EXT,outfile,ORI_EXT)) {
		fprintf(stderr,"Unable to generate output filename.\n");
		return;
		}
	if (!(fp = fopen(infile,"r"))) {
		fprintf(stderr,"Unable to open \"%s\".\n",infile);
		return;
		}
	if (ssioOpen(outfile,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\".\n",outfile);
		fclose(fp);
		return;
		}
	head.n_data = 0; /* bogus value -- fixed later */
	head.time = time++;
	head.iMagicNumber = SSIO_MAGIC_STANDARD; /* enforce standard format */
	if (ssioOriHead(&ssio,&head)) {
		fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	while ((rv = fscanf(fp,"%lf%lf%lf%lf%lf%lfi",
						&data.p1[X],&data.p1[Y],&data.p1[Z],
						&data.p2[X],&data.p2[Y],&data.p2[Z])) != EOF) {
		if (rv != 6) {
			fprintf(stderr,"Improper input format.\n");
			goto finish;
			}
		if (ssioOriData(&ssio,&data)) {
			fprintf(stderr,"Error writing data.\n");
			goto finish;
			}	
		++head.n_data;
		}
	if (!bQuiet) printf("n_data = %i\n",head.n_data);
	/* redo header */
	ssioRewind(&ssio);
	ssioOriHead(&ssio,&head);
 finish:
	ssioClose(&ssio);
	fclose(fp);
	}

int main(int argc,char *argv[])
{
	int bQuiet = 0; /* default = verbose */

	int i;

	setbuf(stdout,(char *)NULL);
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

/* ort2ori.c */
