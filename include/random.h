#ifndef RANDOM_HINCLUDED
#define RANDOM_HINCLUDED

#include <stdint.h> /* for uint64_t and UINT32_MAX */

typedef uint64_t Ullong; /* unsigned long long requires C99 standard or later */

Ullong randReadUrandom(void);
void randSeedGenerator(Ullong ullSeed);
double randUniform(void);
double randRayleigh(void);
double randGaussian(void);
double randPoisson(double dMean);

#endif /* !RANDOM_HINCLUDED */
