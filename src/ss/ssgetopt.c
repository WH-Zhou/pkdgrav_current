/*
 ** ssgetopt.c -- DCR 08/15/01
 ** ==========
 ** Command line parser for the C shell.
 ** Note this code is not specific to ss_core.
 **
 ** Usage: ssgetopt option-list option [ arg [ arg ... ] ]
 **
 ** where "option-list" is the allowed options in getopt() format,
 ** "option" is the option to test against, and the "arg"s are the
 ** actual arguments passed.  If one of the arguments matches the
 **	desired option, and the option itself takes an argument, that
 ** argument is printed as a string to stdout().  If the option
 ** doesn't take an argument (and no argument is given), "1" will
 ** be printed as an integer to stdout().  An invalid option or an
 ** option with an invalid or missing argument generates an error.
 **
 ** If "option" is a colon (":"), all remaining arguments (after
 ** parsing all valid options) will be printed out instead.
 */

#include <stdio.h>
#include <stdlib.h> /* getopt() sometimes declared here */
#include <unistd.h> /* sometimes here */
#include <string.h> /* for strlen() */
#include <assert.h>
#include "ss.h"

int main(int argc,char *argv[])
{
	char *optlist;
	int option,c;

	/*DEBUG
	for (c=0;c<argc;c++)
		(void) printf("argc=%i argv=%s\n",c,argv[c]);
	 */

	if (argc < 3) /* minimum: program name, option list, and target option */
		goto usage;

	assert(argv[0] != NULL && argv[1] != NULL && argv[2] != NULL);

	optlist = argv[1];
	if (strlen(argv[2]) != 1)
		goto usage;
	option = argv[2][0];

	optind = 3; /* start at argument list */

	while ((c = getopt(argc,argv,optlist)) != EOF) {
		if (c == '?')
			return 1; /* error encountered */
		if (c == option) {
			if (optarg != NULL)
				(void) printf("%s\n",optarg);
			else
				(void) printf("1\n"); /* flag on */
			return 0;
			}
		}

	if (option != ':') /* ':' is a special case */
		return 0;

	for (;optind<argc;optind++) { /* output remaining arguments */
		(void) printf("%s",argv[optind]);
		if (optind <= argc - 1)
			(void) printf(" ");
		}

	(void) printf("\n");

	return 0;

 usage:
	(void) fprintf(stderr,"Usage: %s option-list option [ arg [ arg ... ] ]\n",argv[0]);
	(void) fprintf(stderr,"option-list uses getopt() format.\n");
	return 1;
	}

/* ssgetopt.c */
