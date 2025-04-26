/*
** Routines for random number generation (normally used only for
** generating initial conditions, but in some circumstances needed
** during a run) are collected here for convenience.
*/

#include <stdio.h>
#include <unistd.h> /* for getpid() and getppid() */
#include <time.h> /* for time() */
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "random.h"

#ifdef SSIO_USE_MPI
#include "mpi.h"
#endif

/* MWC256 from Usenet posting by G. Marsaglia - Period 2^8222 */
static unsigned int Q[256], c=362436;

unsigned int MWC256(void)
{
    unsigned long long t;
    static unsigned char i=255;
    t = 809430660ULL * Q[++i] + c;
    c = (t>>32);
    return (Q[i]=t);
}

/*
    David Jones, UCL Bioinformatics Group
    Good Practice in (Pseudo) Random Number Generation for
    Bioinformatics Applications
*/

double randUniform
(void)
{
    double x;
    unsigned int a, b;
    a = MWC256() >> 6; /* Upper 26 bits */
    b = MWC256() >> 5; /* Upper 27 bits */
    x = (a * 134217728.0 + b) / 9007199254740992.0;
    return x;
}

/* following routines adapted from NRiC(2e) */

double randGaussian(void) /* gasdev() in NRiC(2e) */
/*
** Returns a normally distributed deviate with zero mean and unit
** variance, using randUniform() as the source of uniform deviates
** (assumed already seeded).
*/
{
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (iset == 0)
    {
        /*
        ** We don't have an extra deviate handy, so pick two uniform
        ** numbers in the square extending from -1 to +1 in each
        ** direction, see if they are in the unit circle, and if they
        ** are not, try again.
        */
        do
        {
            v1 = 2.0*randUniform() - 1.0;
            v2 = 2.0*randUniform() - 1.0;
            rsq = v1*v1 + v2*v2;
        }
        while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0*log(rsq)/rsq);
        /*
        ** Now make the Box-Muller transformation to get two normal
        ** deviates.  Return one and save the other for next time.
        */
        gset = v1*fac;
        iset = 1; /* set flag */
        return v2*fac;
    }
    else
    {
        /*
        ** We have an extra deviate handy, so unset the flag, and
        ** return it.
        */
        iset = 0;
        return gset;
    }
}

static double gammln(double xx)
/*
** Returns the value ln[Gamma(xx)] for xx > 0.
*/
{
    double x,y,tmp,ser;
    static double cof[6]= {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.12086505973866179e-2,-0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    ser = 1.000000000190015;

    for (j=0; j<5; j++)
    {
        ser += cof[j]/++y;
    }

    return -tmp + log(2.5066282746310005*ser/x);
}

double randPoisson(double dMean) /* poidev() in NRiC(2e) */
/*
** Returns as a floating-point number an integer value that is a
** random deviate drawn from a Poisson distribution of mean "dMean",
** using randUniform() as a source of uniform random deviates (assumed
** already seeded).
*/
{
    static double sq,alxm,g,oldm=(-1.0); /* oldm is a flag for whether dMean has changed since last call */
    double em,t,y;

    if (dMean < 12.0)   /* use direct method */
    {
        if (dMean != oldm)
        {
            oldm = dMean;
            g = exp(-dMean);
        }

        em = -1.0;
        t = 1.0;

        do
        {
            /*
            ** Instead of adding exponential deviates it is equivalent
            ** to multiply uniform deviates.  We never actually have
            ** to take the log, merely compare to the pre-computed
            ** exponential.
            */
            ++em;
            t *= randUniform();
        }
        while (t > g);
    }
    else   /* use rejection method */
    {
        if (dMean != oldm)
        {
            /*
            ** If dMean has changed since the last call, then
            ** precompute some functions that occur below.
            */
            oldm = dMean;
            sq = sqrt(2.0*dMean);
            alxm = log(dMean);
            g = dMean*alxm - gammln(dMean + 1.0);
        }

        do
        {
            do
            {
                y = tan(M_PI*randUniform()); /* y is a deviate from a Lorentzian comparison function */
                em = sq*y + dMean; /* em is y, shifted and scaled */
            }
            while (em < 0.0);   /* reject if in regime of zero probability */

            em = floor(em); /* the trick for integer-valued distributions */
            t = 0.9*(1.0 + y*y)*exp(em*alxm - gammln(em + 1.0) - g);
            /*
            ** Above, "t" is the ratio of the desired distribution to
            ** the comparison function; we accept or reject by
            ** comparing it to another uniform deviate.  The factor of
            ** 0.9 is chosen so that t never exceeds 1.
            */
        }
        while (randUniform() > t);
    }

    return em;
}

double randRayleigh
(void)
{
    /* based on expdev(), but for Rayleigh distribution */

    double dum;

    do
    {
        dum = randUniform();
    }
    while (dum == 0);

    return sqrt(-2.0*log(dum));
}

void seedRandom(int iProcID) /* based on Ranq1 constructor in NR(3e) */
{
    /* seeds this processor's random number generator */

    unsigned int seed;

    /*
    // seed should be 32 bit uint otherwise implicit type conversion below could remove useful bits
    seed = time(NULL) % getpid() + getppid(); // random-ish seed
    */

    // 'better', and process of filling Q is defined below so can always reproduce
    FILE * fp;
    assert(fp = fopen("/dev/urandom", "r"));
    assert(fread(&seed, sizeof(seed), 1, fp));
    assert(fclose(fp) != EOF);

    int ii;

    // has to be > 0 - any arbitrary value would do
    if (seed == 0)
    {
        seed = 1;
    }

// make absolutely sure none have the same seed
#ifdef SSIO_USE_MPI
    int initd;
    MPI_Initialized(&initd);

    if (!initd)
    {
        MPI_Init(NULL, NULL);
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Bcast(&seed, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    seed += rank;
    // stops small number of processes having similar seed
    // node counts of >32 will start to circle round again, but will still be different (slightly)
    seed = (seed >> rank) | (seed << sizeof(seed)*8 - rank);
#endif

    // set c and fill Q with seed
    c = seed;

    for (ii = 0; ii < 256; ii++)
    {
        Q[ii] = c + c*ii;
    }

    printf("seedRandom(pid): processor %i random seed = %llu\n", iProcID, seed);

    // discard some numbers to advance it a bit
    for (ii = 0; ii < 2024; ii++)
    {
        MWC256();
    }

#ifdef SSIO_USE_MPI
    // make sure it's different again
    c = Q[rank % 256];
    c += rank;

    for (ii = 0; ii < 256; ii++)
    {
        Q[ii] = c + c*ii;
    }

    // discard some numbers to advance it a bit
    for (ii = 0; ii < 2024; ii++)
    {
        MWC256();
    }
#endif
}

