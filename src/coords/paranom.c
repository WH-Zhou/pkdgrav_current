/*
 ** This function solves Barker's equation, M = (2^1.5)*(D + D*D*D/3), for D.
 ** It is mainly called by conversion from the Delaunay elements to 
 ** cartesian coordinates.
 **
 ** Joachim Stadel, Jan. 11, 1995
 */
#include <math.h>
#include <assert.h>
#include "coords.h"

double dParAnom(double M)
{
	double sgm,A,B,desc,D,m,m1;
	int i;

	if (M < 0.0) {
		M = -M;
		sgm = -1.0;
		}
	else if (M > 0.0) {
		sgm = 1.0;
		}
	else return(0.0);
	m = 0.75*M/sqrt(2.0);
	desc = sqrt(1.0 + 1.0/(m*m));
	A = pow(m*(desc + 1.0),1.0/3.0);
	B = pow(m*(desc - 1.0),1.0/3.0);
	D = A - B;
	/*
	 ** Add two root polishing steps, which looks like it is enough.
	 */
	for (i=0;i<2;++i) {
		m = D*D*D/3.0 + D - 0.5*M/sqrt(2.0);
		m1 = D*D + 1.0;
		D -= m/m1;
		}
	return(sgm*D);
	}
