#ifndef SS_HINCLUDED
#define SS_HINCLUDED

/*
 ** ss.h -- DCR 97-08-06
 ** ====
 ** Header file specific to Solar System utilities.
 */

#include <math.h>
#include <ssdefs.h>

#include "getopt.h"
/* for getopt() with intel compiler */
#ifdef __INTEL_COMPILER
extern int optind;
extern char* optarg;
#endif

/* Math macros */

#undef PI
#ifdef M_PI
# define PI M_PI
#endif

#ifndef PI
# error Unable to find a value for PI!
#endif

#undef TWO_PI
#define TWO_PI (2.0*(PI))

#define DEG_TO_RAD ((PI)/180.0)
#define RAD_TO_DEG (180.0/(PI))

/* From http://asa.usno.navy.mil/ and wikipedia */

#define AU		1.495978707e11	/* One au in metres */
#define JD		86400.0			/* One Julian day in seconds (24*60*60) */
#define SID_YR	3.15581497635e7	/* One sidereal year in seconds (2000.0) */
#define TROP_YR 3.15569252507e7	/* One tropical year in seconds (2000.0) */
#define JUL_YR	(365.25*JD)		/* One Julian year in seconds */
#define M_SUN	1.9884e30		/* Solar mass in kilograms */
#define R_SUN	6.96e8			/* Solar radius in metres */
#define M_EARTH	5.9722e24		/* Earth mass in kilograms */
#define G		6.6743e-11		/* Gravitation constant in mks */
#define GAUSS_K	0.01720209895	/* Gaussian gravitational constant (not used) */

/* Unit scaling */

#define M_SCALE M_SUN
#define L_SCALE AU
#define T_SCALE (SID_YR/TWO_PI) /* was (JD/GAUSS_K) */
#define V_SCALE (L_SCALE/T_SCALE)
#define W_SCALE (1.0/(TWO_PI*T_SCALE))
#define D_SCALE (M_SCALE/(L_SCALE*L_SCALE*L_SCALE))

/* Currently support the four giant planets */

#define MAX_NUM_PLANETS 4 /* Jupiter, Saturn, Uranus, Neptune */

/* Definitions for ssg.c, ssh.c, and sso.c */

#define MAP_EXT ".map"
#define SSH_EXT ".his"
#define SSO_EXT ".osc"

#endif
