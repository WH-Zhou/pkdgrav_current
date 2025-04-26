
/* rdpar.h -- Definitions for rdpar routines DCR 93-03-23 */

#include <copyright.h>

#define NEND		-1	/* Flags for terminating multiple reads */
#define NEND_LABEL	"NULL"

#define MAX_RDPAR_SIZE	1000 	/* Maximum number of lines in parameter file */

#ifndef MAX_STR_LEN		/* Maximum string length */
#	define MAX_STR_LEN 256
#endif

extern void			/* External voids */
	OpenPar(const char *fpstr),
	ReadInt(const char *label, int *x),
	ReadLng(const char *label, long int *x),
	ReadFlo(const char *label, float *x),
	ReadDbl(const char *label, double *x),
	ReadStr(const char *label, char *str, int len),
	ClosePar(void);

extern int			/* External ints */
	ReadNInt(const char *label, int *x, int n),
	ReadNLng(const char *label, long int *x, int n),
	ReadNFlo(const char *label, float *x, int n),
	ReadNDbl(const char *label, double *x, int n),
	ReadNStr(const char *label, char **str, int n, int len);
