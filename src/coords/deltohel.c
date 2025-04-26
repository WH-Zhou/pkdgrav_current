/*
 ** This function converts from Delaunay orbital elements to heliocentric
 ** cartesian coordinates. It can handle elliptic, hyperbolic and
 ** parabolic orbits. This function currently cannot deal with the
 ** rectilinear orbits, this ability may be added later.
 ** The code requires the three functions for solving the "Eccentric
 ** Anomallies" E, F and D for elliptic, hyperbolic and parabolic orbits.
 ** For a reference to this code consult ORBITAL MOTION by
 ** A. E. Roy (sec 4.11).
 **
 **	Joachim Stadel, Jan. 13, 1995
 */
#include <math.h>
#include <assert.h>
#include "constants.h"
#include "coords.h"

void deltohel(double mu,struct delaunay *d,struct helio *h)
{
	double et,nu,det,dnu;
	double caop,saop,clan,slan,cinc,sinc,l1,m1,n1,l2,m2,n2;

	assert(d->ecc >= 0.0);
	if (d->ecc < 1.0) {
		/*
		 ** Elliptic orbit. Careful, may not be true for rectilinear ellipse!
		 ** Solve Kepler's Equation.
		 */
		double E,cE,sE,rt,vs;

		E = dEccAnom(d->mea,d->ecc);
		/*
		 ** Must solve for et,nu and det and dnu.
		 */
		cE = cos(E);
		sE = sin(E);
		rt = sqrt(1.0 - d->ecc*d->ecc);
		et = d->sma*(cE - d->ecc);
		nu = d->sma*rt*sE;
		vs = sqrt(mu*d->sma)/(d->sma*(1.0 - d->ecc*cE));
		det = -vs*sE;
		dnu = vs*rt*cE;
		}
	else if (d->ecc > 1.0) {
		/*
		 ** Hyperbolic orbit.
		 */
		double F,cF,sF,rt,vs;

		F = dHypAnom(d->mea,d->ecc);
		/*
		 ** Must solve for et,nu and det and dnu.
		 */
		cF = cosh(F);
		sF = sinh(F);
		rt = sqrt(d->ecc*d->ecc - 1.0);
		et = d->sma*(d->ecc - cF);
		nu = d->sma*rt*sF;
		vs = sqrt(mu*d->sma)/(d->sma*(d->ecc*cF - 1.0));
		det = -vs*sF;
		dnu = vs*rt*cF;
		}
	else {
		/*
		 ** Parabolic orbit.
		 */
		double D,p,vs;

		D = dParAnom(d->mea);
		/*
		 ** Must solve for et,nu and det and dnu.
		 */
		p = 2.0*d->sma;
		nu = p*D;
		et = nu*nu + 0.5*p;
		vs = 1.0/(1.0 + D*D);
		det = 2.0*sqrt(mu*p)*vs*D;
		dnu = sqrt(mu/p)*vs;
		}
	/*
	 ** Compute direction cosines.
	 */
	caop = cos(d->lop - d->lan);	/* lop - lan = argument of perihelion */
	saop = sin(d->lop - d->lan);
	clan = cos(d->lan);
	slan = sin(d->lan);
	cinc = cos(d->inc);
	sinc = sin(d->inc);
	l1 = clan*caop - slan*saop*cinc;
	m1 = slan*caop + clan*saop*cinc;
	n1 = saop*sinc;
	l2 = -clan*saop - slan*caop*cinc;
	m2 = -slan*saop + clan*caop*cinc;
	n2 = caop*sinc;
	/*
	 ** Calculate heliocentric coordinates.
	 */
	h->x = l1*et + l2*nu;
	h->y = m1*et + m2*nu;
	h->z = n1*et + n2*nu;
	h->vx = l1*det + l2*dnu;
	h->vy = m1*det + m2*dnu;
	h->vz = n1*det + n2*dnu;
	}
