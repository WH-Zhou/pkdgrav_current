#include <stdio.h>
#include <strings.h> /* for strcasecmp() */
#include <math.h> /* for M_SQRT2 and M_PI */
#include <assert.h>
#include <ss.h>

#define OUTFILENAME "shape.ss"

enum {TetraNF=4,CubeNF=6,OctaNF=8,DodecaNF=12,IcosaNF=20}; /* number of faces */

enum {TetraNV=4,CubeNV=8,OctaNV=6,DodecaNV=20,IcosaNV=12}; /* number of vertices */

typedef double Vector[3];

#define SQRT2 M_SQRT2
#define SQRT3 1.73205080756887729353
#define SQRT6 (SQRT2*SQRT3)

/* vertex positions from Wolfram-Alpha: center at origin, unit sides */

const Vector vTetra[] = { /* 4 vertices (enclosing radius ~ 0.612372) */
  {0.0,0.0,SQRT2/SQRT3 - 0.5/SQRT6},
  {-0.5/SQRT3,-0.5,-0.5/SQRT6},
  {-0.5/SQRT3,0.5,-0.5/SQRT6},
  {1.0/SQRT3,0.0,-0.5/SQRT6}
};

const Vector vCube[] = { /* 8 vertices (enclosing radius ~ 0.866025) */
  {-0.5,-0.5,-0.5},
  {-0.5,-0.5,0.5},
  {-0.5,0.5,-0.5},
  {-0.5,0.5,0.5},
  {0.5,-0.5,-0.5},
  {0.5,-0.5,0.5},
  {0.5,0.5,-0.5},
  {0.5,0.5,0.5}
};

const Vector vOcta[] = { /* 6 vertices (enclosing radius ~ 0.707107) */
  {-1.0/SQRT2,0,0},
  {0.0,1.0/SQRT2,0.0},
  {0.0,0.0,-1.0/SQRT2},
  {0.0,0.0,1.0/SQRT2},
  {0.0,-1.0/SQRT2,0.0},
  {1.0/SQRT2,0.0,0.0}
};

const Vector vDodeca[] = { /* 20 vertices (enclosing radius ~ 1.40126) */
  {-1.37638192047117353821,0.0,2.62865556059566803015e-1},
  {1.37638192047117353821,0.0,-2.62865556059566803015e-1},
  {-4.2532540417601996609e-1,-1.30901699437494742410,2.62865556059566803015e-1},
  {-4.2532540417601996609e-1,1.30901699437494742410,2.62865556059566803015e-1},
  {1.11351636441160673519,-8.09016994374947424102e-1,2.62865556059566803015e-1},
  {1.11351636441160673519,8.09016994374947424102e-1,2.62865556059566803015e-1},
  {-2.62865556059566803015e-1,-8.09016994374947424102e-1,1.11351636441160673519},
  {-2.62865556059566803015e-1,8.09016994374947424102e-1,1.11351636441160673519},
  {-6.88190960235586769105e-1,-0.5,-1.11351636441160673520},
  {-6.88190960235586769105e-1,0.5,-1.11351636441160673520},
  {6.881909602355867691e-1,-0.5,1.11351636441160673519},
  {6.881909602355867691e-1,0.5,1.11351636441160673519},
  {8.5065080835203993218e-1,0.0,-1.11351636441160673520},
  {-1.11351636441160673520,-8.09016994374947424102e-1,-2.62865556059566803015e-1},
  {-1.11351636441160673520,8.09016994374947424102e-1,-2.62865556059566803015e-1},
  {-8.5065080835203993218e-1,0.0,1.11351636441160673519},
  {2.62865556059566803015e-1,-8.09016994374947424102e-1,-1.11351636441160673520},
  {2.62865556059566803015e-1,8.09016994374947424102e-1,-1.11351636441160673520},
  {4.2532540417601996609e-1,-1.30901699437494742410,-2.62865556059566803015e-1},
  {4.2532540417601996609e-1,1.30901699437494742410,-2.62865556059566803015e-1}
};

const Vector vIcosa[] = { /* 12 vertices (enclosing radius ~ 0.951057) */
  {0.0,0.0,-9.51056516295153572116e-1},
  {0.0,0.0,9.51056516295153572116e-1},
  {-8.5065080835203993218e-1,0.0,-4.25325404176019966092e-1},
  {8.5065080835203993218e-1,0.0,4.25325404176019966092e-1},
  {6.88190960235586769105e-1,-0.5,-4.25325404176019966092e-1},
  {6.88190960235586769105e-1,0.5,-4.25325404176019966092e-1},
  {-6.88190960235586769105e-1,-0.5,4.25325404176019966092e-1},
  {-6.88190960235586769105e-1,0.5,4.25325404176019966092e-1},
  {-2.62865556059566803014e-1,-8.090169943749474241e-1,-4.25325404176019966092e-1},
  {-2.62865556059566803014e-1,8.090169943749474241e-1,-4.25325404176019966092e-1},
  {2.62865556059566803014e-1,-8.090169943749474241e-1,4.25325404176019966092e-1},
  {2.62865556059566803014e-1,8.090169943749474241e-1,4.25325404176019966092e-1}
};
  
static int write_data(int nVert,SSDATA *data)
{
	SSIO ssio;
	SSHEAD head;
	int i;

	if (ssioOpen(OUTFILENAME,&ssio,SSIO_WRITE)) {
		fprintf(stderr,"Unable to open \"%s\" for writing.\n",OUTFILENAME);
		return 1;
		}

	head.time = 0.0;
	head.n_data = nVert;
	head.iMagicNumber = SSIO_MAGIC_STANDARD;
	if (ssioHead(&ssio,&head)) {
		fprintf(stderr,"Unable to write header.\n");
		ssioClose(&ssio);
		return 1;
		}

	for (i=0;i<nVert;i++)
		if (ssioData(&ssio,&data[i])) {
			fprintf(stderr,"Error writing data (particle %i).\n",i);
			ssioClose(&ssio);
			return 1;
			}

	ssioClose(&ssio);

	return 0;
	}

static void construct_solid(int nVert,int nFaces,double dSideLength,double dRadius,double dMass,int iColor,SSDATA *d)
{
	const Vector *vVert;
	int i;

	for (i=0;i<nVert;i++) {
		d[i].mass = dMass;
		d[i].radius = dRadius;
		switch (nFaces) {
		case TetraNF:
			vVert = &vTetra[i];
			break;
		case CubeNF:
			vVert = &vCube[i];
			break;
		case OctaNF:
			vVert = &vOcta[i];
			break;
		case DodecaNF:
			vVert = &vDodeca[i];
			break;
		case IcosaNF:
			vVert = &vIcosa[i];
			break;
		default:
			assert(0); /* shouldn't be here! */
			}
		d[i].pos[0] = (*vVert)[0]*dSideLength;
		d[i].pos[1] = (*vVert)[1]*dSideLength;
		d[i].pos[2] = (*vVert)[2]*dSideLength;
		d[i].vel[0] = d[i].vel[1] = d[i].vel[2] = 0.0;
		d[i].spin[0] = d[i].spin[1] = d[i].spin[2] = 0.0;
		d[i].color = iColor;
		d[i].org_idx = i;
		}
	}

static void usage(const char *achProgName)
{
	fprintf(stderr,
			"Usage: %s n-faces side-length radius units density color\n"
			"where n-faces is: 4=tetra; 6=cube; 8=octa; 12=dodec; 20=icosa;\n"
			"side-length is length of any side (they're all the same);\n"
			"radius is particle radius;\n"
			"units is \"cm\", \"m\", \"km\", or \"au\";\n"
			"density is in g/cc only and refers to the particle density;\n"
			"and color is an ssdraw color (integer).\n",achProgName);
	exit(1);
	}

int main(int argc,char *argv[])
{
	SSDATA *data;
	double dSideLength,dRadius,dScale,dDensity,dMass;
	int nFaces,nVert,iColor;

	if (argc != 7)
		usage(argv[0]);

	nFaces = atoi(argv[1]);
	nVert = 0; /* to suppress compiler warning */
	switch (nFaces) {
	case TetraNF:
		puts("Tetrahedron selected.");
		nVert = TetraNV;
		break;
	case CubeNF:
		puts("Cube selected.");
		nVert = CubeNV;
		break;
	case OctaNF:
		puts("Octahedron selected.");
		nVert = OctaNV;
		break;
	case DodecaNF:
		puts("Dodecahedron selected.");
		nVert = DodecaNV;
		break;
	case IcosaNF:
		puts("Icosahedron selected.");
		nVert = IcosaNV;
		break;
	default:
		fprintf(stderr,"Invalid number of faces (%i).\n",nFaces);
		usage(argv[0]);
		}

	dSideLength = atof(argv[2]);
	if (dSideLength <= 0.0) {
		fprintf(stderr,"Invalid side length (%g).\n",dSideLength);
		usage(argv[0]);
		}

	dRadius = atof(argv[3]);
	if (dRadius <= 0.0) {
		fprintf(stderr,"Invalid particle radius (%g).\n",dRadius);
		usage(argv[0]);
		}
	if (dRadius < 0.5*dSideLength)
		puts("Note: particles will not be touching.");
	else if (dRadius > 0.5*dSideLength)
		puts("Note: particles will be overlapping.");

	dScale = 1.0;
	if (strcasecmp(argv[4],"cm") == 0) {
		puts("Using centimeters.");
		dScale = 0.01/L_SCALE; /* cm -> au */
		}
	else if (strcasecmp(argv[4],"m") == 0) {
		puts("Using meters.");
		dScale = 1.0/L_SCALE; /* m -> au */
		}
	else if (strcasecmp(argv[4],"km") == 0) {
		puts("Using meters.");
		dScale = 1000.0/L_SCALE; /* km -> au */
		}
	else if (strcasecmp(argv[4],"au") == 0)
		puts("Using astronomical units.");
	else {
		fprintf(stderr,"Invalid length unit (\"%s\").\n",argv[4]);
		usage(argv[0]);
		}
	dSideLength *= dScale;
	dRadius *= dScale;

	dDensity = atof(argv[5]);
	if (dDensity <= 0.0) {
		fprintf(stderr,"Invalid density (%g).\n",dDensity);
		usage(argv[0]);
		}
	dMass = 4.0/3.0*M_PI*dRadius*dRadius*dRadius*dDensity*1000.0/D_SCALE;

	iColor = (int) atof(argv[6]);
	if (iColor < 0 || iColor > 255) {
		fprintf(stderr,"Invalid color (%i).\n",iColor);
		usage(argv[0]);
		}
	if (iColor == 0)
		puts("WARNING: chosen particle color is black.");

	data = (SSDATA *) malloc((size_t) nVert*sizeof(SSDATA));
	assert(data != NULL);
	construct_solid(nVert,nFaces,dSideLength,dRadius,dMass,iColor,data);
	if (write_data(nVert,data)) {
		free((void *) data);
		return 1;
		}

	free((void *) data);

	printf("Data written to %s.\n",OUTFILENAME);

	return 0;
	}
