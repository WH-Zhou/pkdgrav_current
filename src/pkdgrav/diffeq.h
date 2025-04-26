#ifndef DIFFEQ_HINCLUDED
#define DIFFEQ_HINCLUDED

#include "floattype.h"

typedef int (*diffeqDerivsT)(FLOAT t, const FLOAT fVars[], FLOAT fDerivs[],
                             void *pUserData);

void diffeqIntegrate(FLOAT fVars[], int nVars, FLOAT fStart, FLOAT fStop,
                     FLOAT fAccuracy, FLOAT fStepTry, FLOAT fStepMin,
                     diffeqDerivsT funcDerivs, void *pUserData);

#endif /* !DIFFEQ_HINCLUDED */
