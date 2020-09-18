#include "lfa.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
  {"lfa_threshold", (DL_FUNC) &lfa_threshold, 2},
  {"lfa_scaling", (DL_FUNC) &lfa_scaling, 2},
  {"centerscale_c", (DL_FUNC) &centerscale_c, 1},
  {"lreg_c", (DL_FUNC) &lreg_c, 4},
  {"mv_c", (DL_FUNC) &mv_c, 2},
  {"tmv_c", (DL_FUNC) &tmv_c, 2},
  {NULL, NULL, 0}
};

void R_init_lfa(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

