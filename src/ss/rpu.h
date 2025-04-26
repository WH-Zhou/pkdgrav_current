#ifndef RPU_HINCLUDED
#define RPU_HINCLUDED

/*
 ** rpu.h -- DCR 98-09-23
 ** =====
 ** Header file specific to rubble pile utilities.
 */

#include <ss.h>
#include <boolean.h>
#include <vector.h>
#include <colors.h>

#define MAJOR(prp) ((prp)->axis_ord[0])
#define INTER(prp) ((prp)->axis_ord[1])
#define MINOR(prp) ((prp)->axis_ord[2])

typedef struct {
	/* traditional quantities */
	double mass;
	double radius;
	VECTOR pos;
	VECTOR vel;
	VECTOR spin;
	int color;
	/* special quantities */
	int agg_id;
	MATRIX inertia;
	VECTOR moments;
	int mom_ord[N_DIM];
	MATRIX axes; /* principal axes are ROWS of this matrix, i.e., applying this matrix will convert from space coords to body coords */
        VECTOR DEEVE; /* DEEVE dimensions */
	VECTOR ang_mom;
	double ang_mom_mag;
	double kin_energy;
	double dyn_inertia;
	double eff_spin;
	double rot_idx;
	VECTOR axis_len; /* semi-axes *//*DEBUG! rename so it's OBVIOUS these are semi-axes!!!*/
	int axis_ord[N_DIM];
	double density;
	double packing;
	/* particle data */
	int n_particles;
	SSDATA *data;
	} RUBBLE_PILE;

/* aggregate definitions (patterned after macros in pkdgrav aggs.h) */

int AGG_IDX(const SSDATA *d);
#define AGG_IDX(d) (-1 - (d)->org_idx)

int AGG_SET_IDX(SSDATA *d,int IDX);
#define AGG_SET_IDX(d,IDX) {assert((IDX) >= 0); (d)->org_idx = -1 - (IDX);}

/* rpu function prototypes */

double rpuVolSph(double r);
double rpuVolEll(const VECTOR a);
BOOLEAN rpuInEllipsoid(const VECTOR r0, double R, const VECTOR a);
void rpuMalloc(RUBBLE_PILE *rp);
void rpuRealloc(RUBBLE_PILE *rp);
void rpuFree(RUBBLE_PILE *rp);
void rpuCopyData(const RUBBLE_PILE *rps, RUBBLE_PILE *rpd, const size_t offset);
void rpuCalcMass(RUBBLE_PILE *rp);
void rpuCalcRadius(RUBBLE_PILE *rp);
void rpuCalcPos(RUBBLE_PILE *rp);
void rpuCalcVel(RUBBLE_PILE *rp);
void rpuCalcInertia(RUBBLE_PILE *rp);
void rpuCalcAxisLengths(RUBBLE_PILE *rp);
void rpuCalcAxes(RUBBLE_PILE *rp);
void rpuCalcSpin(RUBBLE_PILE *rp);
void rpuCalcColor(RUBBLE_PILE *rp);
void rpuCalcAggID(RUBBLE_PILE *rp);
void rpuCalcDensity(RUBBLE_PILE *rp);
void rpuCalcPacking(RUBBLE_PILE *rp);
void rpuScaleMass(RUBBLE_PILE *rp, double f);
void rpuScaleRadius(RUBBLE_PILE *rp, double f, BOOLEAN just_particles);
void rpuApplyPos(RUBBLE_PILE *rp);
void rpuApplyVel(RUBBLE_PILE *rp);
void rpuApplyColor(RUBBLE_PILE *rp);
void rpuApplyAggID(RUBBLE_PILE *rp);
void rpuRotate(RUBBLE_PILE *rp, MATRIX rot, BOOLEAN body);
void rpuAddSpin(RUBBLE_PILE *rp, BOOLEAN body);
void rpuAddAngMom(RUBBLE_PILE *rp);
void rpuAnalyze(RUBBLE_PILE *rp);

#endif
