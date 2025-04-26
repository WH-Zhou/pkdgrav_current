/*
 ** This function converts from heliocentric cartesian coordinates
 ** to Delaunay orbital elements. It can handle elliptic, hyperbolic and
 ** parabolic orbits. This function currently cannot deal with the
 ** rectilinear orbits, this ability may be added later. For a reference
 ** to this code consult ORBITAL MOTION by A. E. Roy (sec 4.11).
 **
 **	Joachim Stadel, Jan. 11, 1995
 */
#include <math.h>
#include <assert.h>
#include "constants.h"
#include "coords.h"

#define TINY 1.0e-14 /*DEBUG it's crashed with 1.0e-15 on Linux-alpha...*/

void heltodel(double mu,struct helio *h,struct delaunay *d)
{
	double r,v2,rv,hx,hy,hz,hxy,h2,c,rcoswf,rsinwf,wf,f,w,E,F,D;

	r = sqrt(h->x*h->x + h->y*h->y + h->z*h->z);
	v2 = h->vx*h->vx + h->vy*h->vy + h->vz*h->vz;
	rv = h->x*h->vx + h->y*h->vy + h->z*h->vz;
	hx = h->y*h->vz - h->z*h->vy;
	hy = h->z*h->vx - h->x*h->vz;
	hz = h->x*h->vy - h->y*h->vx;
	hxy = hx*hx + hy*hy;
	h2 = hxy + hz*hz;
	/*
	 ** Make sure we are not dealing with a rectilinear orbit.
	 */
	assert(h2 > 0.0);
	/*
	 ** Roy's Book is not correct on the calculation of d->lan, no
	 ** changes of sign are required for -ve or +ve values of hz!
	 */
	d->lan = atan2(hx,-hy);
	/*
	 ** Convert from -Pi <= d->lan <= Pi to 0 <= d->lan < 2*Pi.
	 */
	if (d->lan < 0.0) d->lan += 2.0*CONST_PI;
	/*
	 ** This could have been written d->inc = acos(hz/sqrt(h2)).
	 ** However, we should avoid using acos and asin where atan2
	 ** can be used. In the case of small inclination acos gets a 
	 ** value of nearly 1.0 which limits the precision, using asin
	 ** in this case does not suffer this precision problem. However,
	 ** atan2 takes care of all these considerations.
	 */
	d->inc = atan2(sqrt(hxy),hz);
	/*
	 ** Calculate, (angle of perihelion) + (true anomaly). wf for short.
	 */
	rcoswf = h->x*cos(d->lan) + h->y*sin(d->lan);
	if (d->inc > 0.0) {
		rsinwf = h->z/sin(d->inc);
		}
	else {
		rsinwf = h->y*cos(d->lan) - h->x*sin(d->lan);
		}
	wf = atan2(rsinwf,rcoswf);
	if (wf < 0.0) wf += 2.0*CONST_PI;
	/*
	 ** Find the energy.
	 */
	c = 0.5*v2 - mu/r;
	if (c < 0.0) {
		double x;
		/*
		 ** Elliptical orbit.
		 */
		d->sma = -0.5*mu/c;
		/* DCR: Added precision-limit checks -- DCR 97-10-09 */
		x = 1.0 - h2/(mu*d->sma);
		if (fabs(x) < TINY)
			d->ecc = TINY;
		else
			d->ecc = sqrt(x);
		x = h2/(mu*r) - 1.0;
		if (fabs(x) > d->ecc)
			f = (x < 0 ? CONST_PI : 0);
		else
			f = acos(x/d->ecc);
		x = 1.0 - r/d->sma;
		if (fabs(x) > d->ecc)
			E = (x < 0 ? CONST_PI : 0);
		else
			E = acos(x/d->ecc);
		if (rv < 0.0) {
			w = wf + f;
			if (w > 2.0*CONST_PI) w -= 2.0*CONST_PI;
			E = -E;
			}
		else {
			w = wf - f;
			if (w < 0.0) w += 2.0*CONST_PI;
			}
		/*
		 ** Now w should be between 0 and 2*Pi!
		 */
		d->lop = d->lan + w;
		d->mea = E - d->ecc*sin(E);
		}
	else if (c > 0.0) {
		/*
		 ** Hyperbolic orbit.
		 */
		d->sma = 0.5*mu/c;
		d->ecc = sqrt(1.0 + h2/(mu*d->sma));
		f = acos((h2/(mu*r) - 1.0)/d->ecc);	 /* note: same as elliptic case */
		F = acosh((1.0 + r/d->sma)/d->ecc);
		if (rv < 0.0) {
			w = wf + f;
			if (w > 2.0*CONST_PI) w -= 2.0*CONST_PI;
			F = -F;
			}
		else {
			w = wf - f;
			if (w < 0.0) w += 2.0*CONST_PI;
			}
		/*
		 ** Now w should be between 0 and 2*Pi!
		 */
		d->lop = d->lan + w;
		d->mea = d->ecc*sinh(F) - F;
		}
	else {
		/*
		 ** Parabolic orbit.
		 */
		d->sma = 2.0*mu/v2;
		d->ecc = 1.0;
		f = acos(h2/(mu*r) - 1.0);
		D = sqrt(r/d->sma - 1.0);
		if (rv < 0.0) {
			w = wf + f;
			if (w > 2.0*CONST_PI) w -= 2.0*CONST_PI;
			D = -D;
			}
		else {
			w = wf - f;
			if (w < 0.0) w += 2.0*CONST_PI;
			}
		/*
		 ** Now w should be between 0 and 2*Pi!
		 */
		d->lop = d->lan + w;
		d->mea = pow(2.0,1.5)*D*(1.0 + D*D/3.0);
		}
	}
