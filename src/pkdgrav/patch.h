#ifndef PATCH_HINCLUDED
#define PATCH_HINCLUDED

#ifdef SLIDING_PATCH

enum {STRIP_LEFT_ONLY,STRIP_RIGHT_ONLY,STRIP_BOTH};

typedef struct {
  int bPatch; /* non-zero if sliding patch enabled */
  double dCentMass;
  double dOrbDist;
  double dWidth; /* same as dxPeriod */
  double dLength; /* same as dyPeriod */
  /* following values computed from other parameters */
  double dOrbFreq;
  double dOrbFreqZ2;
  /* external perturber parameters */
  int bExtPert;
  double dPertOrbDist;
  double dPertMass;
  double dPertMaxZ; /* maximum vertical displacement of perturber from orbital plane */
  double dPertOrbFreqZ; /* vertical frequency of perturber oscillation */
  double dPertPhase; /* initial phase of perturber's orbit, relative to patch */
  double dPertPhaseZ; /* initial phase of perturber's vertical oscillation, relative to orbital plane */
  double dPertOrbFreq; /* computed */
  /* azimuthal wrap randomization parameters */
  int bRandAzWrap;
  int bNoRandomX; /* diable randomization of x component? */
  int nWrapAttempts; /* max number of wrap attempts (without overlap) before giving up */
  int iStripOption; /* 0=left only,1=right only,2=both */
  double dStripInner; /* between 0 and dWidth/2 */
  double dStripOuter; /* ditto (dStripOuter > dStripInner) */
  double dVelDispX,dVelDispY; /* desired velocity dispersion components after wrap */
  double dAvgVertAmp; /* desired mean vertical oscillation amplitude after wrap */
  double dAvgMass; /* should not include any perturbers */
  double dDelta; /* timestep */
} PATCH_PARAMS;

#endif /* SLIDING_PATCH */

#endif /* !PATCH_HINCLUDED */
