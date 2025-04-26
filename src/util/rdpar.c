/*
 * rdpar ver 1.1 -- DCR 93-03-23
 * =============================
 * A parameter file parsing utility based on the original RDPAR.
 *
 * Modifications by David Earn:
 *  7 May 1993: added __STDIO_H_SEEN and __STDLIB_H_SEEN lines for Convex cc
 *
 * SYNOPSIS:
 *
 * Call OpenPar(const char *filename) to prepare a parameter file for reading,
 * and ClosePar(void) when you're done. The parameter file should be organized
 * with parameter "labels" or key words at the start of each line, followed
 * by the data, separated with whitespace (no commas). Comments may be
 * added anywhere by prefacing with "!" or "#" (the remainder of the line
 * is ignored). Currently, "!" and "#" cannot be used inside character
 * strings. Single integers, long integers, floats, doubles, and strings can
 * be read by calling:
 *    void ReadInt(char *label, int *myint);
 *    void ReadLng(char *label, long int *mylng);
 *    void ReadFlo(char *label, float *myflo);
 *    void ReadDbl(char *label, double *mydbl);
 *    void ReadStr(char *label, char *mystr, int maxstrlen);
 * where "label" is the keyword(s) as used in the parameter file. Note that
 * strings must be given a maximum length. For ReadStr(), strings need not
 * be bracketed by quotes (") (which are removed if present); apostrophes
 * (') are left alone. To read more than one item of a given type, use:
 *    int ReadNInt(char *label, int *myintarray, int numelem);
 *    int ReadNLng(char *label, long int *mylngarray, int numelem);
 *    int ReadNFlo(char *label, float *myfloarray, int numelem);
 *    int ReadNDbl(char *label, double *mydblarray, int numelem);
 *    int ReadNStr(char *label, char **mystrarray, int numelem, int maxstrlen);
 * where "numelem" is the number of elements of each type to be read. These
 * should be on the same line as the key words (up to 256 characters),
 * separated by whitespace. Strings MUST be delimited by quotes in this case.
 * The ReadN commands return the constant NEND (-1) if the keyword
 * NEND_LABEL (the string "NULL", without quotes) appears in place of a
 * data item. This is useful for reading multiple data sets with the same
 * keywords, e.g.,
 *    x, y, z positions		1 2 0
 *    x, y, z positions		3 4 0
 *    x, y, z positions		NULL NULL NULL
 * At the moment it is primarily the responsibility of the programmer to
 * ensure that valid arguments are passed to the rdpar routines. Further
 * functionality may be added in future versions. Comments welcome!
 *
 * USAGE:
 *
 * Build the library by executing the makefile (type "make"). To use the
 * rdpar routines, add the library (-lrdpar) to the argument list of your
 * program compile command. You will also need to specify the path of the
 * library through the "-L" option to cc. Be sure to #include "rdpar.h" in
 * your source file(s). Compilation creates the executable "rptest" which
 * parses the parameter file "rptest.par" (included in the distribution)
 * as a test.
 *
 */

/* Header files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __STDLIB_H_SEEN  /* Convex */
/*#    include <malloc.h>*/
#endif
#include "rdpar.h"

/* Some #define's */

#define VERBOSE 0				/* Non-zero to switch on */

#define EOS	'\0'				/* String termination marker */
#define TAB	'\011'				/* TAB character (whitespace) */
#define SPACE	' '				/* Space (whitespace) */
#define COM1	'!'				/* Comment marker */
#define COM2	'#'				/* Comment marker */
#define QUOTE	'\"'				/* String delimiter */
#define CR	'\n'				/* Carriage return */

/* Useful macros */

#define MIN(x,y) ((x) < (y) ? (x) : (y))	/* Minimum of x & y */
#define MAX(x,y) ((x) > (y) ? (x) : (y))	/* Maximum of x & y */

/* Local variables */

FILE	*fp;					/* Pointer to par file */
char	*data[MAX_RDPAR_SIZE];			/* Storage for file in memory */
int	line;					/* Counter/total no. of lines */

/* Local routines */

static void
	strip(char **str),
	stripSpace(char **str),
	stripLeadingSpace(char **str),
	stripTrailingSpace(char **str),
	cutToSpace(char **str),
	error(const char *str);

static int
	stripQuotes(char **str),
	my_cindex(const char *src, char *tgt);

static char
	*findLabel(const char *label),
	*my_index(char *tgt, char c),
	*my_sindex(char *src, char *tgt);

/* Non-ANSI C may not declare these... */

#ifndef _STDIO_H           /* Sun */
#ifndef __STDIO_H__        /* Silicon Graphics */
#ifndef _H_STDIO           /* IBM RISC 6000 */
#ifndef __STDIO_H_SEEN     /* Convex */
#ifndef _STDIO_H_          /* Alpha */

extern int fclose(FILE *), fprintf(FILE *, const char *, ...), printf(const char *, ...);

#endif
#endif
#endif
#endif
#endif

void OpenPar(const char *fpstr)
{
	/* Opens par file and reads lines into memory */

	static int first_call = 1;
	char buf[MAX_STR_LEN], *bp;

	/* Make sure only one file is open at a time */

	if (first_call) {
		fp = NULL;
		first_call = 0;
	}

	if (fp != NULL)
		error("OpenPar(): File already open...close it first.");

	/* Attempt to assign file pointer */

	if ((fp = fopen(fpstr, "r")) == NULL) {
		char buf2[256];
		(void) snprintf(buf2,256,"OpenPar(): Unable to open \"%s\".",
				fpstr);
		buf2[255] = '\0';
		error(buf2);
	}

	/* Read in data lines, stripping comments */

	line = 0;

	while (fgets(buf, MAX_STR_LEN, fp) != NULL) {
		bp = buf;
		strip(&bp);
		if (strlen(bp) > 0) {
			if (line == MAX_RDPAR_SIZE)
				error("OpenPar(): Max file size exceeded.");
			data[line] = (char *) malloc((unsigned) strlen(bp) + 1);
			(void) strcpy(data[line++], bp);
		}
	}

	/* Close par file */

	(void) fclose(fp);

#if (VERBOSE != 0)
	(void) printf("OpenPar(): %i line(s) read from \"%s\"\n",
			line, fpstr);
#endif
}

void ClosePar(void)
{
	/* Resets file pointer and deallocates storage */

	fp = NULL;

	for (--line; line >= 0; line--)
		free(data[line]);
}

void ReadInt(const char *label, int *x)
{
	/* Reads integer x associated with keyword(s) label */

	*x = atoi(findLabel(label));
}

void ReadLng(const char *label, long int *x)
{
	/* Reads long integer x associated with keyword(s) label */

	*x = atol(findLabel(label));
}

void ReadFlo(const char *label, float *x)
{
	/* Reads float x associated with keyword(s) label */

	*x = (float) atof(findLabel(label));
}
	
void ReadDbl(const char *label, double *x)
{
	/* Reads double x associated with keyword(s) label */

	*x = atof(findLabel(label));
}

void ReadStr(const char *label, char *str, int len)
{
	/* Reads string str (max length len) associated with keyword(s) label */

	char *substr;

	substr = findLabel(label);
	(void) stripQuotes(&substr); /* Remove quotation marks, if any */
	(void) strncpy(str, substr, (size_t) len); /* NOTE: strncpy() automatically pads
					     with 0's if strlen(substr) < len */
}

int ReadNInt(const char *label, int *x, int n)
{
	/* Reads n integers associated with label into array x */

	int i;
	char *substr;

	substr = findLabel(label);

	/* Loop over integers (assume they're all there!) */

	for (i = 0; i < n; i++) {

		/* Check for end of list read */

		if (my_cindex(NEND_LABEL, substr) == 0)
			return NEND;

		/* Assign value */

		x[i] = atoi(substr);

		/* Skip to next value in string */

		if (i < n - 1) {
			cutToSpace(&substr);
			stripLeadingSpace(&substr);
		}
	}

	return 0;
}

int ReadNLng(const char *label, long int *x, int n)
{
	/* Reads n long integers associated with label into array x */

	int i;
	char *substr;

	substr = findLabel(label);

	for (i = 0; i < n; i++) {
		if (my_cindex(NEND_LABEL, substr) == 0)
			return NEND;
		x[i] = atol(substr);
		if (i < n - 1) {
			cutToSpace(&substr);
			stripLeadingSpace(&substr);
		}
	}

	return 0;
}

int ReadNFlo(const char *label, float *x, int n)
{
	/* Reads n floats associated with label into array x */

	int i;
	char *substr;

	substr = findLabel(label);

	for (i = 0; i < n; i++) {
		if (my_cindex(NEND_LABEL, substr) == 0)
			return NEND;
		x[i] = (float) atof(substr);
		if (i < n - 1) {
			cutToSpace(&substr);
			stripLeadingSpace(&substr);
		}
	}

	return 0;
}

int ReadNDbl(const char *label, double *x, int n)
{
	/* Reads n doubles associated with label into array x */

	int i;
	char *substr;

	substr = findLabel(label);

	for (i = 0; i < n; i++) {
		if (my_cindex(NEND_LABEL, substr) == 0)
			return NEND;
		x[i] = atof(substr);
		if (i < n - 1) {
			cutToSpace(&substr);
			stripLeadingSpace(&substr);
		}
	}

	return 0;
}

int ReadNStr(const char *label, char **str, int n, int len)
{
	/* Reads n strings (max length len) assoc'd with label into array str */

	int i, l;
	char *substr;

	substr = findLabel(label);

	for (i = 0; i < n; i++) {
		if (my_cindex(NEND_LABEL, substr) == 0)
			return NEND;

		/* Strip off quotation marks around current string */

		l = stripQuotes(&substr);

                /* Get length of string and copy to array */

                l = MIN(len, l);
                (void) strncpy(str[i], substr, (size_t) l);

                /* Add EOS marker */

                str[i][l - (l == len ? 1 : 0)] = EOS;

		/* Skip to next string */

		if (i < n - 1) {
			substr += l; /*DEBUG what if l = MAX?*/
			stripLeadingSpace(&substr);
		}
	}

	return 0;
}

static char *findLabel(const char *label)
{
	/*
	 * Returns pointer to first data item after label. Note that
	 * the label is then "deleted" from memory, so subsequent searches
	 * for the same label will find the next occurence in the par file.
	 *
	 */

	int l;
	char dum[MAX_STR_LEN],*workspace,*substr;

	workspace = strncpy(dum, label, MAX_STR_LEN);
	
	/* Strip the label, just in case */

	strip(&workspace);

	/* Scan through the stored file a line at a time */

	for (l = 0; l < line; l++)

		/* Look for a match */

		if ((substr = my_sindex(workspace, data[l])) != NULL) {

			/* Strip space up to the first data item */

			stripLeadingSpace(&substr);

			/* Remove the label from mem with the help of an EOS */

			data[l][0] = EOS;

			/* Return pointer to data */

			return substr;
		}

	/* Fatal error if the label is not found */

	{
		char errstr[MAX_STR_LEN];

		(void) sprintf(errstr, "findLabel(): Label \"%s\" not found.",
				label);
		error(errstr);
	}

	/* The following is to keep lint happy */

	return NULL;
}

static char *my_index(char *tgt, char c)
{
	/*
	 * Returns a pointer to the first occurence of character c in string
	 * tgt, or NULL if c does not occur in the string. Note that my_index()
	 * is not an ANSI C function which is why it is given here explicitly.
	 *
	 */

	size_t len, cc;

	/* Automatic failure if tgt has zero length */

	if ((len = strlen(tgt)) == 0)
		return NULL;

	/* Search for the first occurence */

	for (cc = 0; cc < len; cc++)
		if (c == tgt[cc])
			return (tgt + cc);

	/* No match */

	return NULL;
}

static int my_cindex(const char *src, char *tgt)
{
	/*
	 * Returns the numerical position in string tgt where the
	 * entire string src first occurs. The code -1 is returned
	 * if src does not occur anywhere in tgt.
	 *
	 */

	size_t slen, tlen, c, cc;

	/* Automatic failure if src cannot fit inside tgt */

	if ((slen = strlen(src)) > (tlen = strlen(tgt)) || slen == 0)
		return -1;

	/* Do a character-by-character comparison until a match is found */

	for (c = 0; c < tlen - slen + 1; c++) {
		cc = 0;
		while (src[cc] == tgt[c + cc])
			if (++cc == slen)
				return (int) c;
	}

	/* No match */

	return -1;
}

static char *my_sindex(char *src, char *tgt)
{
	/*
	 * Returns a pointer to the position in tgt at which string src
	 * first differs from string tgt.
	 *
	 */

	size_t len, c;

	/* Automatic failure if src is longer than tgt */

	if ((len = strlen(src)) > strlen(tgt))
		return NULL;

	/* Search for the first difference */

	for (c = 0; c < len; c++)
		if (src[c] != tgt[c])
			return NULL;

	/* Return pointer to first differing character */

	return (tgt + c);
}

static void strip(char **str)
{
	/*
	 * Removes comments and carriage returns as well as leading
	 * and trailing whitespace from string str.
	 *
	 */

	/*DEBUG can't handle ! or # or \n inside quotes...*/

	char *ptr;

	/* Return if string already empty */

	if (strlen(*str) == 0)
		return;

	/* Stick EOS's in place of comment markers or CR */

	if ((ptr = my_index(*str, COM1)) != NULL)
		*ptr = EOS;

	if ((ptr = my_index(*str, COM2)) != NULL)
		*ptr = EOS;

	if ((ptr = my_index(*str, CR)) != NULL)
		*ptr = EOS;

	/* Remove whitespace */

	stripSpace(str);
}

static void stripSpace(char **str)
{
	/* Removes leading and trailing whitespace */

	if (strlen(*str) == 0)
		return;

	stripLeadingSpace(str);

	stripTrailingSpace(str);
}

static void stripLeadingSpace(char **str)
{
	/* Removes leading whitespace, incrementing str pointer as required */

	if (strlen(*str) == 0)
		return;

	while ((*str)[0] == SPACE || (*str)[0] == TAB)
		++(*str);
}

static void stripTrailingSpace(char **str)
{
	/* Removes trailing whitespace */

	size_t l;

	if (strlen(*str) == 0)
		return;

	while ((*str)[l = strlen(*str) - 1] == SPACE || (*str)[l] == TAB)
		(*str)[l] = EOS;
}

static int stripQuotes(char **str)
{
	/*
	 * Strips 2 quotes from str, incrementing pointer str if necessary.
	 * The length of the delimited string (or portion thereof) is returned.
	 *
	 */

	size_t l;
	char *ptr;

	/* Leading quote assumed to be in position 0 */

	if ((*str)[0] == QUOTE)
		++(*str);

	/* Find next occurence of a quote, and shift string down to cover it */

	l = strlen(*str);

	if ((ptr = my_index(*str, QUOTE)) != NULL) {
		l -= strlen(ptr);
		for (; strlen(ptr) > 0; ++ptr)
			*ptr = *(ptr + 1);
	}

	/* Return the length of the string that was delimited by quotes */

	return (int) l;
}

static void cutToSpace(char **str)
{
	/* Increments pointer str until the first character is not whitespace */

	if (strlen(*str) == 0)
		return;

	while ((*str)[0] != SPACE && (*str)[0] != TAB)
		++(*str);
}

static void error(const char *str)
{
	/* Fatal error routine; will dump core if possible */

	(void) fprintf(stderr, "\007Rdpar error in %s\n", str);

	abort();
}

/* rdpar.c */
