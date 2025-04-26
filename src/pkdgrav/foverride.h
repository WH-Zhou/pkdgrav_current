#ifndef FOVERRIDE_HINCLUDED
#define FOVERRIDE_HINCLUDED

#ifdef SPRINGS
#include "springs.h"
#endif

#ifdef DEM
#include "dem.h"
#endif

enum {FO_NONE,FO_VANDERWAALS,FO_STRENGTH};

typedef struct {
  int iForceOverrideOption; /* see enum above */
  /* following for FO_STRENGTH option */
#ifdef SPRINGS
  SPRING_PARAMS SP;
#endif
#ifdef DEM
  DEM_PARAMS DP;
#endif
  } FOVERRIDE_PARAMS;

#endif /* !FOVERRIDE_HINCLUDED */
