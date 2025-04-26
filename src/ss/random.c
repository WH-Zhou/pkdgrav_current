/*
** Random number generation routines and helper functions.
** Includes a standalone test driver.
*/

/*#define STANDALONE*/ /* uncomment to compile standalone test */

#ifdef NUMREC
/* compile with -DNUMREC */
#endif

#ifdef GSL
/* compile with -DGSL -I/path/to/headers -L/path/to/lib -lgsl */
#endif

#include <stdio.h>
#include <stdlib.h>
#include "random.h"

Ullong randReadUrandom(void)
{
    /*
    ** Returns sizeof(Ullong) bytes from /dev/urandom as a Ullong.
    ** Since this uses device access, it may be best for seeding.
    */

    FILE *fp;
    Ullong ullNum;

    fp = fopen("/dev/urandom","rb");
    if (fp == NULL) {
	fprintf(stderr, "Unable to open /dev/urandom for reading.\n");
	exit(1);
	}
    fread(&ullNum, sizeof(ullNum), 1, fp);
    fclose(fp);
    return ullNum;
    }

#ifdef NUMREC

#include <math.h>

/*
** Following routines adapted from Numerical Recipes (3rd ed).
*/

/* used by randSeedGenerator() -- don't use this value as the seed! */
static Ullong ullState = 4101842887655102017LL;

static Ullong rand_uniform_int64(void)
/*
** Based on Ranq1.int64(), NR(3e) 7.1.3.
**
** Recommended generator for everyday use.  The period is ~1.8e19.
*/
{
    ullState ^= ullState >> 21;
    ullState ^= ullState << 35;
    ullState ^= ullState >> 4;

    return ullState*2685821657736338717LL;
    }

void randSeedGenerator(Ullong ullSeed)
/*
** Based on Ranq1 constructor, NR(3e) 7.1.3.
*/
{
    /* seeds this processor's random number generator */

    ullState ^= ullSeed;
    ullState = rand_uniform_int64();
    }

double randUniform(void)
/*
** Based on Ranq1.doub(), in NR(3e) 7.1.3.
*/
{
    return 5.42101086242752217e-20*(double)rand_uniform_int64();
    }

/*
** Following routines adapted from Numerical Recipes in C (2nd ed).
*/

double randRayleigh(void) /* based on expdev(), NRiC(2e) 7.2 */
{
    double dum;
    
    do
        dum = randUniform();
    while (dum == 0.0);
    return sqrt(-2.0*log(dum));
}

double randGaussian(void) /* based on gasdev(), NRiC(2e) 7.2 */
/*
** Returns a normally distributed deviate with zero mean and unit
** variance, using randUniform() as the source of uniform deviates
** (assumed already seeded).
*/
{
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (iset == 0) {
	/*
	** We don't have an extra deviate handy, so pick two uniform
	** numbers in the square extending from -1 to +1 in each
	** direction, see if they are in the unit circle, and if they
	** are not, try again.
	*/
	do {
	    v1 = 2.0*randUniform() - 1.0;
	    v2 = 2.0*randUniform() - 1.0;
	    rsq = v1*v1 + v2*v2;
	    } while (rsq >= 1.0 || rsq == 0.0);
	fac = sqrt(-2.0*log(rsq)/rsq);
	/*
	** Now make the Box-Muller transformation to get two normal
	** deviates.  Return one and save the other for next time.
	*/
	gset = v1*fac;
	iset = 1; /* set flag */
	return v2*fac;
	}
    else {
	/*
	** We have an extra deviate handy, so unset the flag, and
	** return it.
	*/
	iset = 0;
	return gset;
	}
    }

static double gammln(double xx) /* based on gammln(), NRiC(2e) 6.1 */
/*
** Returns the value ln[Gamma(xx)] for xx > 0.
*/
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.12086505973866179e-2,-0.5395239384953e-5};
    int j;
    
    y = x = xx;
    tmp = x + 5.5;
    ser = 1.000000000190015;
    for (j=0;j<5;j++)
	ser += cof[j]/++y;
    return -tmp + log(2.5066282746310005*ser/x);
    }

double randPoisson(double dMean) /* based on poidev(), NRiC(2e) 7.3 */
/*
** Returns as a floating-point number an integer value that is a
** random deviate drawn from a Poisson distribution of mean "dMean",
** using randUniform() as a source of uniform random deviates (assumed
** already seeded).
*/
{
    static double sq,alxm,g;
    static double oldm=(-1.0); /* flag if dMean has changed since last call */
    double em,t,y;

    if (dMean < 12.0) { /* use direct method */
	if (dMean != oldm) {
	    oldm = dMean;
	    g = exp(-dMean);
	    }
	em = -1.0;
	t = 1.0;
	do {
	    /*
	    ** Instead of adding exponential deviates it is equivalent
	    ** to multiply uniform deviates.  We never actually have
	    ** to take the log, merely compare to the pre-computed
	    ** exponential.
	    */
	    ++em;
	    t *= randUniform();
	    } while (t > g);
	}
    else { /* use rejection method */
	if (dMean != oldm) {
	    /*
	    ** If dMean has changed since the last call, then
	    ** precompute some functions that occur below.
	    */
	    oldm = dMean;
	    sq = sqrt(2.0*dMean);
	    alxm = log(dMean);
	    g = dMean*alxm - gammln(dMean + 1.0);
	    }
	do {
	    do {
		/* y is a deviate from a Lorentzian comparison function */
		y = tan(M_PI*randUniform());
		em = sq*y + dMean; /* em is y, shifted and scaled */
		} while (em < 0.0); /* reject if in regime of 0 probability */
	    em = floor(em); /* the trick for integer-valued distributions */
	    t = 0.9*(1.0 + y*y)*exp(em*alxm - gammln(em + 1.0) - g);
	    /*
	    ** Above, "t" is the ratio of the desired distribution to
	    ** the comparison function; we accept or reject by
	    ** comparing it to another uniform deviate.  The factor of
	    ** 0.9 is chosen so that t never exceeds 1.
	    */
	    } while (randUniform() > t);
	}
    return em;
    }

#endif /* NUMREC */

#ifdef GSL

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng; /* pointer to memory for random number generator */

void randSeedGenerator(Ullong ullSeed)
{
    /* allocates and seeds this processor's random number generator */
    
    /* note: should really call gsl_rng_free() when done... */
    rng = gsl_rng_alloc(gsl_rng_default); /* MT19937 generator */
    if (rng == NULL) {
        fprintf(stderr, "randSeedGenerator(): Unable to allocate generator.\n");
        exit(1);
    }
    gsl_rng_set(rng, ullSeed % UINT32_MAX); /* mod to 32 bits (unsigned long) */
}

double randUniform(void)
{
    return gsl_rng_uniform(rng); /* range [0,1) */
}

double randRayleigh(void)
{
    return gsl_ran_rayleigh(rng, 1.0); /* unit scale parameter */
}

double randGaussian(void)
{
    return gsl_ran_gaussian(rng, 1.0); /* zero mean, unit variance */
}

double randPoisson(double dMean)
{
    return gsl_ran_poisson(rng, dMean); /* returns as double */
}

#endif /* GSL */

#ifdef STANDALONE
int main(void) {
    Ullong ullSeed = randReadUrandom();
    randSeedGenerator(ullSeed);
    printf("Uniform random number = %g\n", randUniform());
    printf("Rayleigh random number (scale parameter 1) = %g\n", randRayleigh());
    printf("Gaussian random number (mean 0, width 1) = %g\n", randGaussian());
    printf("Poisson random number (mean 100) = %g\n", randPoisson(100.0));
    return 0;
}
#endif /* STANDALONE */
