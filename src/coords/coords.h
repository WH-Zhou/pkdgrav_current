#ifndef COORDS_HINCLUDED
#define COORDS_HINCLUDED

struct delaunay {
	double mea;		/* mean anomaly 'M' */
	double lop;		/* longitude of perihelion 'pomega' */
	double lan;		/* longitude of ascending node 'Omega' */
	double sma;		/* semi-major axis 'a' */
	double ecc;		/* eccentricity 'e' */
	double inc;		/* inclination 'i' */
	};

struct helio {
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	};

double dEccAnom(double,double);
double dHypAnom(double,double);
double dParAnom(double);

void deltohel(double mu,struct delaunay *d,struct helio *h);
void heltodel(double mu,struct helio *h,struct delaunay *d);

#endif
