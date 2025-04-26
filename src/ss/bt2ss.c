/*
 ** bt2ss.c -- DCR 97-03-04
 ** =======
 ** Converts box_tree ASCII data to Solar System binary data.
 **
 ** NOTE: time field in Solar System binary data is set to increment by 1.
 */

#include <stdio.h>
#include <unistd.h> /* for getopt() and getpid() */
#include <ssdefs.h>
#include <vector.h>

#define BT_EXT ".bt"

static void convert(char *infile, int bQuiet)
{
	static double time = 0;

	FILE *fp;
	SSIO ssio;
	SSHEAD head;
	SSDATA data;
	char outfile[MAXPATHLEN];
	int idx,rv,dum;

	if (!bQuiet) printf("%s: ",infile);
	if (ssioNewExt(infile,BT_EXT,outfile,SS_EXT)) {
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
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	idx = 0;
	while ((rv = fscanf(fp,"%i%i%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%i",&dum,
						&data.org_idx,&data.mass,&data.radius,
						&data.pos[X],&data.pos[Y],&data.pos[Z],
						&data.vel[X],&data.vel[Y],&data.vel[Z],
						&data.spin[X],&data.spin[Y],&data.spin[Z],
						&data.color)) != EOF) {
		if (rv != 14) {
			fprintf(stderr,"Improper input format.\n");
			goto finish;
			}
		if (idx >= 0 && dum != idx++) {
			fprintf(stderr,"WARNING: Non-standard indices ignored.\n");
			idx = -1;
			}
		if (ssioData(&ssio,&data)) {
			fprintf(stderr,"Error writing data.\n");
			goto finish;
			}	
		++head.n_data;
		}
	if (!bQuiet) printf("n_data = %i\n",head.n_data);
	/* redo header */
	ssioRewind(&ssio);
	ssioHead(&ssio,&head);
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

/* bt2ss.c */
