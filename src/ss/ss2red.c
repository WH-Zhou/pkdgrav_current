/*
 ** ss2red.c -- DCR 6/17/08
 ** =======-
 ** Converts Solar System binary data to reduced form.
 */

#include <ssio.h>
#include <vector.h>

#define RED_EXT ".r"

static void convert(char *infile)
{
	SSIO ssio_i,ssio_o;
	SSHEAD head;
	SSDATA d;
	SSRDATA dr;
	char outfile[MAXPATHLEN];
	int i;

	printf("%s: ",infile);
	snprintf(outfile,MAXPATHLEN,"%s%s",infile,RED_EXT);
	if (ssioOpen(infile,&ssio_i,SSIO_READ)) {
		fprintf(stderr,"Unable to open \"%s\" for reading.\n",infile);
		return;
		}
	if (ssioOpen(outfile,&ssio_o,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",outfile);
		ssioClose(&ssio_i);
		return;
		}
	if (ssioHead(&ssio_i,&head)) {
		fprintf(stderr,"Corrupt header.\n");
		goto finish;
		}
	printf("time = %e, n_data = %i\n",head.time,head.n_data);
	if (head.n_data <= 0) {
		fprintf(stderr,"Invalid input data format.\n");
		goto finish;
		}
	switch(head.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		fprintf(stderr,"File already in reduced format.\n");
		goto finish;
	default:
		fprintf(stderr,"Unrecognized ss file magic number (%i).\n",head.iMagicNumber);
		goto finish;
		}
	head.iMagicNumber = SSIO_MAGIC_REDUCED;
	if (ssioHead(&ssio_o,&head)) {
		fprintf(stderr,"Error writing header.\n");
		goto finish;
		}
	for (i=0;i<head.n_data;i++) {
		if (ssioData(&ssio_i,&d)) {
			fprintf(stderr,"Corrupt data.\n");
			goto finish;
			}
		ssioStandardToReduced(&d,&dr);
		if (ssioDataReduced(&ssio_o,&dr)) {
			fprintf(stderr,"Error writing data.\n");
			goto finish;
			}
		}
 finish:
	ssioClose(&ssio_o);
	ssioClose(&ssio_i);
	}

int main(int argc,char *argv[])
{
	int i;

	setbuf(stdout,(char *)NULL);
	if (argc <= 1) {
		fprintf(stderr,"Usage: %s file [ file ... ]\n",argv[0]);
		exit(1);
		}
	for (i=1;i<argc;++i)
		convert(argv[i]);
	return 0;
	}

/* ss2red.c */
