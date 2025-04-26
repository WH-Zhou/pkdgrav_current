/*
 ** rp2.c -- DCR 98-12-08
 ** =====
 ** Utility for rps script. Given 2 ss files (each should contain
 ** a single rubble pile), returns:
 **
 **   M   B   D   R   V   W   RR
 **
 ** where M = total mass
 **       B = sum of radii
 **       D = effective (volume-weighted) bulk density
 **       R = effective radius
 **       V = critical speed (was effective escape speed)
 **       W = effective maximum spin magnitude
 **      RR = effective Roche radius
 **      VE = mutual escape speed in units of V
 */

#include <rpu.h>
#include <math.h>
#include <assert.h>

static int
get_rp(const char *filename,RUBBLE_PILE *rp)
{
	SSIO ssio;
	SSHEAD h;
	int i;

	assert(filename != NULL && rp != NULL);

	if (ssioOpen(filename,&ssio,SSIO_READ)) {
		(void) fprintf(stderr,"Unable to open \"%s\"\n",filename);
		return 1;
		}

	if (ssioHead(&ssio,&h) || h.n_data < 0) {
		(void) fprintf(stderr,"Corrupt header\n");
		(void) ssioClose(&ssio);
		return 1;
		}

	if (h.n_data == 0) {
		(void) fprintf(stderr,"No data found!");
		(void) ssioClose(&ssio);
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

	rp->n_particles = h.n_data;
	rpuMalloc(rp);

	for (i=0;i<rp->n_particles;i++)
		if (ssioData(&ssio,&rp->data[i])) {
			(void) fprintf(stderr,"Corrupt data\n");
			(void) ssioClose(&ssio);
			return 1;
			}

	(void) ssioClose(&ssio);

	rpuAnalyze(rp);

	return 0;
	}

int
main(int argc,char *argv[])
{
	RUBBLE_PILE rp1,rp2;
	double m,mu,b,d,x,r,v,w,rr,ve;

	setbuf(stdout,(char *)NULL);

	if (argc != 3) {
		(void) fprintf(stderr,"Usage: %s ss-file1 ss-file2\n",argv[0]);
		return 1;
		}

	if (get_rp(argv[1],&rp1)) return 1;
	if (get_rp(argv[2],&rp2)) return 1;

	m = rp1.mass + rp2.mass;
	mu = rp1.mass*rp2.mass/m;
	b = rp1.radius + rp2.radius;
	d = m/(rpuVolEll(rp1.axis_len) + rpuVolEll(rp2.axis_len));
	x = 4*PI/3*d;
	r = pow(m/x,1.0/3);
/*	v = sqrt(2*m/r);*/
	v = sqrt(1.2*SQ(m)/(mu*r)); /* kinetic energy = binding energy */
	w = sqrt(x);
	rr = 1.52*pow(m/d,1.0/3);
	ve = sqrt(2*m/b);

	(void) printf("%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n",
		m,b,d,r,v,w,rr,ve/v);

	rpuFree(&rp1);
	rpuFree(&rp2);

	return 0;
	}

/* rp2.c */
