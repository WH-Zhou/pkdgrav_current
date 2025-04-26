/*
 ** ssg.c -- DCR 97-10-21 (revised 7/1/03)
 ** =====
 ** Constructs Solar System particle genealogy from data output.
 **
 ** STRATEGY: a starting file and ending file are provided by the
 ** user. The program determines the final original index for each
 ** particle in the starting file based on collision history. (The
 ** ending file simply provides the timestamp.) A map file is
 ** output that gives the names of the files used and a list of
 ** the final original indices of each starting particle.
 **
 ** To trace an individual particle's evolution, use sst after
 ** creating a map file with ssg.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>     /* for getopt() */
#include <string.h>
#include <rpc/rpc.h>	/* for XDR routines */
#include <sys/param.h>	/* for MAXPATHLEN */
#include <assert.h>
#include <ss.h>

#define MAPFILENAME "ssg.map"

/*#define VERBOSE*/ /*DEBUG this should really be a run-time option*/

enum {UNKNOWN,BINARY,TEXT};
enum {RESTART,MERGE,BOUNCE,FRAG};

typedef struct {
	FILE *fp;
	XDR xdrs;
	int type;
	} LOG;

typedef struct {
	double time;
	int action,idx1,idx2,new;
	} EVENT;

static void write_map(char *ssfile1,char *ssfile2,int n,int *map)
{
	FILE *fp;
	int i;

	/* Output data */

	if (!(fp = fopen(MAPFILENAME,"w"))) {
		(void) fprintf(stderr,"Unable to open %s for writing.\n",MAPFILENAME);
		exit(1);
		}

	(void) fprintf(fp,"%s\n%s\n",ssfile1,ssfile2);

	for (i=0;i<n;i++)
		(void) fprintf(fp,"%i\n",map[i]);

	(void) fclose(fp);
	}

#define MAXLINELEN 256

static void read_action(LOG *log,int *action)
{
	if (log->type == BINARY) {
		int n;
		(void) xdr_int(&log->xdrs,&n);
		if (n == -1)
			*action = RESTART;
		else if (n == 1)
			*action = MERGE;
		else if (n == 2)
			*action = BOUNCE; /* shouldn't happen: bounces aren't logged in binary file */
		else if (n > 2)
			*action = FRAG;
		else
			assert(0); /* unknown action type */
		}
	else {
		char line[MAXLINELEN],*r;

		r = fgets(line,MAXLINELEN,log->fp);
		if (feof(log->fp)) {
			(void) fprintf(stderr,"read_action(): Unexpected end of file.\n");
			exit(1);
			}
		assert(r != NULL);
		if (strstr(line,"MERGE"))
			*action = MERGE;
		else if (strstr(line,"BOUNCE"))
			*action = BOUNCE;
		else if (strstr(line,"FRAG"))
			*action = FRAG;
		else {
			(void) fprintf(stderr,"read_action(): Unknown action: \"%s\".\n",line);
			exit(1);
			}
		}
	}

static void read_idx(LOG *log,int *idx)
{
	if (log->type == BINARY) {
		(void) xdr_int(&log->xdrs,idx);
		assert(*idx >= -1); /* -1 ==> restart */
		}
	else {
		char line[MAXLINELEN],*r;

		r = fgets(line,MAXLINELEN,log->fp);
		if (feof(log->fp)) {
			(void) fprintf(stderr,"read_idx(): Unexpected end of file.\n");
			exit(1);
			}
		assert(r != NULL);
		r = strstr(line,"o=");
		assert(r != NULL);
		*idx = atoi(r + 2);
		assert(*idx >= 0);
		}
	}

static void read_time(LOG *log,double *time,int *restart)
{
	*restart = 0;
	if (log->type == BINARY)
		(void) xdr_double(&log->xdrs,time);
	else {
		char line[MAXLINELEN],*r;

		while (/*CONSTCOND*/1) {
			r = fgets(line,MAXLINELEN,log->fp);
			if (feof(log->fp)) return;
			assert(r != NULL);
			if (strstr(line,"RESTART")) {
				r = strchr(line,'=');
				assert(r != NULL);
				*time = atof(r + 1);
				*restart = 1;
				return;
				}
			if (strstr(line,"COLLISION")) {
				r = strchr(line,'=');
				assert(r != NULL);
				*time = atof(r + 1);
				return;
				}
			};
		}
	/*DEBUG! set negative times to zero otherwise initial overlap mergers will be missed!*/
	if (*time < 0.0)
		*time = 0.0;
	}

#undef MAXLINELEN

static void next_event(LOG *log,EVENT *event)
{
	int restart;

	if (feof(log->fp)) goto eof;

	/*
	 ** Restarts for text and binary mode are handled differently:
	 ** in text mode, a restart is indicated when reading the event time;
	 ** in binary mode, a restart is indicated when reading the action.
	 */

	read_time(log,&event->time,&restart);
	if (feof(log->fp)) goto eof;
	if (restart) {
		assert(log->type == TEXT);
		event->idx1 = event->idx2 = event->new = -1;
		event->action = RESTART;
		return;
		}
	read_idx(log,&event->idx1);
	if (feof(log->fp)) goto eof;
	read_idx(log,&event->idx2);
	if (feof(log->fp)) goto eof;
	read_action(log,&event->action);
	if (feof(log->fp)) goto eof;
	if (event->action == RESTART) {
		assert(log->type == BINARY);
		event->new = -1;
		event->action = RESTART;
		return;
		}

	/* sanity check */

	assert(event->idx1 != event->idx2);

	switch (event->action) {
	case MERGE:
		read_idx(log,&event->new);
		if (feof(log->fp)) goto eof;
		assert(event->new == event->idx1 || event->new == event->idx2);
		break;
	case BOUNCE:
		event->new = -1;
		break;
	case FRAG:
		assert(0); /* not implemented yet */
	default:
		assert(0); /* unknown/unexpected event */
		}

	return;

 eof:
	event->action = EOF;
	}

static int make_map(char *filename,double time0,int n0,int map[],LOG *log)
{
	SSIO ssio;
	SSHEAD h;
	EVENT event,event0;
	double time;

	/* Get ending time stamp */

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"%s: error reading header.\n",filename);
		return 1;
		}

	switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		ssioClose(&ssio);
		return 1;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}

	(void) ssioClose(&ssio);

	time = h.time;

#ifdef VERBOSE
	(void) printf("%s: time %e yr, n = %i\n",filename,time*T_SCALE/SID_YR,h.n_data);
#endif

	assert(time >= time0);

	/* Scroll through log until start time */

	do {
		next_event(log,&event);
		} while (event.action != EOF && event.time < time0);

	/* Loop through events up to end time */

	do {
#ifdef VERBOSE
		(void) printf("Event at t=%e yr: ",event.time*T_SCALE/SID_YR);
#endif
		switch (event.action) {
		case MERGE: {
			int i,nc;
#ifdef VERBOSE
			(void) printf("Merge %i+%i --> %i\n",
						  event.idx1,event.idx2,event.new);
#endif
			assert(event.idx1 >= 0 && event.idx1 < n0);
			assert(event.idx2 >= 0 && event.idx2 < n0);
			/* update map */
			for (nc=i=0;i<n0;i++)
				if (map[i] == event.idx1 || map[i] == event.idx2) {
					map[i] = event.new;
					++nc;
					}
			assert(nc >= 2);
			break;
			}
		case BOUNCE:
			assert(log->type != BINARY); /* bounces not recorded in terse log */
#ifdef VERBOSE
			(void) printf("Bounce %i & %i\n",event.idx1,event.idx2);
#endif
			break;
		case FRAG:
			(void) fprintf(stderr,"Frag -- NOT SUPPORTED.\n");
			exit(1);
		case EOF:
			(void) fprintf(stderr,"UNEXPECTED EOF.\n");
			exit(1);
		default:
			(void) fprintf(stderr,"UNEXPECTED/UNKNOWN EVENT.\n");
			exit(1);
			}
		event0 = event;
		next_event(log,&event);
		/*
		 ** Handle restarts.
		 ** Do it this way because with the FIX_COLLAPSE collision option it's
		 ** possible to have negative timesteps...
		 */
		if (event.action == RESTART) {
#ifdef VERBOSE
			(void) printf("Restart at t=%e yr; skipping past t=%e yr...\n",
						  event.time*T_SCALE/SID_YR,event0.time*T_SCALE/SID_YR);
#endif
			do {
				next_event(log,&event);
				} while (event.action != EOF && event.time <= event0.time); /* assumes no simultaneous events */
			if (event.action == event0.action && event.new == event0.new &&
				((event.idx1 == event0.idx1 && event.idx2 == event0.idx2) ||
				 (event.idx1 == event0.idx2 && event.idx2 == event0.idx1)))
				(void) printf("WARNING: Possible duplicate event detected...\n"); /* with just mergers, can be ignored */
			}
		} while (event.action != EOF && event.time <= time);

	return 0;
	}

static int init_map(char *filename,double *time0,int *n0,int **map)
{
	SSIO ssio;
	SSHEAD h;
	SSDATA d;
	int i;

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\".\n",filename);
		return 1;
		}

	/* Get starting time stamp and number of particles */

	if (ssioHead(&ssio,&h)) {
		(void) fprintf(stderr,"%s: error reading header.\n",filename);
		return 1;
		}

	switch(h.iMagicNumber) {
	case SSIO_MAGIC_STANDARD:
		break;
	case SSIO_MAGIC_REDUCED:
		(void) fprintf(stderr,"Reduced ss format not supported.\n");
		ssioClose(&ssio);
		return 1;
	default:
		(void) fprintf(stderr,"Unrecognized ss file magic number (%i).\n",h.iMagicNumber);
		ssioClose(&ssio);
		return 1;
		}

	*time0 = h.time;
	*n0 = h.n_data;

#ifdef VERBOSE
	(void) printf("%s: time %e yr, n = %i\n",filename,*time0*T_SCALE/SID_YR,*n0);
#endif

	assert(*time0 >= 0.0);
	assert(*n0 > 0);

	/* Initialize map */

	*map = (int *) malloc((size_t) *n0*sizeof(int));
	assert(*map != NULL);

	for (i=0;i<*n0;i++) {
		if (ssioData(&ssio,&d)) {
			(void) fprintf(stderr,"%s: error reading data.\n",filename);
			return 1;
			}
		(*map)[i] = d.org_idx;
		if ((*map)[i] < 0) {
			(void) fprintf(stderr,"%s: particle %i has negative original index (%i).\n",filename,i,(*map)[i]);
			return 1;
			}
		}

	(void) ssioClose(&ssio);

	return 0;
	}

static void usage(char *progname)
{
	(void) fprintf(stderr,"Usage: %s [ -b | -t ] ss-file1 ss-file2\n"
				   "where -b specifies binary (terse) collision log\n"
				   "      -t specifies text (verbose) collision log\n",progname);
	exit(1);
	}

int main(int argc,char *argv[])
{
	LOG log;
	double time0;
	int c,n0=-1,*map=NULL;

	/* Disable stdout buffering */

	setbuf(stdout,(char *)NULL);

	/* Check arguments */

	log.type = UNKNOWN;

	while ((c = getopt(argc,argv,"bt")) != EOF)
		switch (c) {
		case 'b':
			if (log.type == TEXT) {
				(void) fprintf(stderr,"Cannot specify both -b and -t.\n");
				usage(argv[0]);
				}
			log.type = BINARY;
			break;
		case 't':
			if (log.type == BINARY) {
				(void) fprintf(stderr,"Cannot specify both -b and -t.\n");
				usage(argv[0]);
				}
			log.type = TEXT;
			break;
		default:
			usage(argv[0]);
			}

	if (optind == argc) usage(argv[0]);

	if (optind != argc - 2) {
		(void) fprintf(stderr,"Must specify exactly two ss files for mapping\n"
					   "(e.g. initial conditions and final output).\n");
		usage(argv[0]);
		}

	/* Attempt to open text or binary collision log */

	switch (log.type) {
	case BINARY:
		if (!(log.fp = fopen(COLL_LOG_BIN,"r"))) {
			(void) fprintf(stderr,"Unable to open %s.\n",COLL_LOG_BIN);
			return 1;
			}
		break;
	case TEXT:
		if (!(log.fp = fopen(COLL_LOG_TXT,"r"))) {
			(void) fprintf(stderr,"Unable to open %s.\n",COLL_LOG_TXT);
			return 1;
			}
		break;
	case UNKNOWN:
		if ((log.fp = fopen(COLL_LOG_BIN,"r"))) {
			FILE *fp;
			if ((fp = fopen(COLL_LOG_TXT,"r"))) {
				(void) fprintf(stderr,"Both binary and text log files found:\n"
							   "Use -b or -t to specify which one to use.\n");
				(void) fclose(fp);
				(void) fclose(log.fp);
				usage(argv[0]);
				}
			(void) printf("Binary log file found.\n");
			log.type = BINARY;
			}
		else if ((log.fp = fopen(COLL_LOG_TXT,"r"))) {
			(void) printf("Text log file found.\n");
			log.type = TEXT;
			}
		else {
			(void) fprintf(stderr,"Unable to find/open log file.\n");
			return 1;
			}
		break;
	default:
		assert(0); /* this can't happen */
		}

	assert(log.type != UNKNOWN);

	/* Initialize from first file */

	(void) init_map(argv[optind],&time0,&n0,&map);

	/* Open XDR stream if binary log */

	if (log.type == BINARY)
		xdrstdio_create(&log.xdrs,log.fp,XDR_DECODE);

	/* Map to second file */

	(void) make_map(argv[optind + 1],time0,n0,map,&log);

	/* Close log */

	if (log.type == BINARY)
		xdr_destroy(&log.xdrs);

	(void) fclose(log.fp);

	/* Write map */

	write_map(argv[optind],argv[optind + 1],n0,map);

	/* All done */

	return 0;
	}

/* ssg.c */
