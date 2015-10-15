#include "lfa.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
  {"lfa_threshold", (DL_FUNC) &lfa_threshold, 2},
  {"lfa_scaling", (DL_FUNC) &lfa_scaling, 2},
  {"centerscale", (DL_FUNC) &centerscale, 1},
  {"center", (DL_FUNC) &center, 1},
  {"lreg", (DL_FUNC) &lreg, 4},
  {"mv", (DL_FUNC) &mv, 2},
  {"tmv", (DL_FUNC) &tmv, 2},
  {NULL, NULL, 0}
};

void
R_init_lfa(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
