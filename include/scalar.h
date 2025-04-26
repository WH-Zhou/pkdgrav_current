/*
 * scalar.h -- DCR 94-04-19
 * ========================
 * Useful macros for common operations on scalars.
 *
 */

/* Basic math macros */

#define POW2(n)  (1 << (n))                             /* 2 to nth power  */
#define EXP10(n) (pow(10.0, (double) (n)))              /* 10 to nth power */
#define SQ(x)    ((x) * (x))                            /* Square of x */
#define CUBE(x)  ((x) * (x) * (x))                      /* Cube of x */
#define ABS(x)   ((x) < 0 ? (- (x)) : (x))              /* Abs val of x */
#define SGN(x)   ((x) == 0 ? 0 : ((x) < 0 ? (-1) : 1))  /* Signum of x */

#ifndef MIN
# define MIN(x,y) ((x) < (y) ? (x) : (y))               /* Min of x and y */
#endif

#ifndef MAX
# define MAX(x,y) ((x) > (y) ? (x) : (y))               /* Max of x and y */
#endif

/*
 * Macros for rough comparisons of double-precision numbers.
 *
 * NOTE: APPROX_EQ() cannot be used to check for near equality with zero.
 *       It may be appropriate to use (ABS(value) < PRECISION) instead,
 *       depending on the situation.
 *
 */

#define PRECISION 1.0e-9 /* This is the best value found so far... */

#define APPROX_EQ(x,y) (ABS((x) - (y)) <= PRECISION * MAX(ABS(x), ABS(y)))
#define APPROX_LT(x,y) ((y) - (x) > PRECISION * ABS(y))
#define APPROX_GT(x,y) ((x) - (y) > PRECISION * ABS(x))
#define APPROX_LE(x,y) (APPROX_EQ((x),(y)) || APPROX_LT((x),(y)))
#define APPROX_GE(x,y) (APPROX_EQ((x),(y)) || APPROX_GT((x),(y)))

/* scalar.h */
