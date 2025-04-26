/*
 ** rotate.c -- DCR 01-08-09
 ** ========
 ** Simple program to rotate a vector.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <vector.h>

int
main(int argc,char *argv[])
{
	/* rotates input vector by specified angle around given rotation axis */

	VECTOR r,n,nxr,rp;
	double phi,cphi,sphi,ndotr;

	if (argc != 8) {
		(void) fprintf(stderr,"Usage: %s x y z phi rx ry rz\n",argv[0]);
		return 1;
		}

	SET_VEC(r,atof(argv[1]),atof(argv[2]),atof(argv[3]));
	assert(MAG(r) > 0.0);
	phi = atof(argv[4]);
	SET_VEC(n,atof(argv[5]),atof(argv[6]),atof(argv[7]));
	assert(MAG(n) > 0.0);

	ndotr = DOT(n,r);
	CROSS(n,r,nxr);
	cphi = cos(phi);
	sphi = sin(phi);

	SCALE_VEC(r,cphi);
	SCALE_VEC(n,ndotr*(1.0-cphi));
	SCALE_VEC(nxr,sphi);

	ADD_VEC(r,n,rp);
	ADD_VEC(rp,nxr,rp);

	(void) printf("%g %g %g\n",rp[X],rp[Y],rp[Z]);

	return 0;
	}
