#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP))

SEXP lfa_threshold(SEXP, SEXP);
SEXP lfa_scaling(SEXP, SEXP);
SEXP centerscale_c(SEXP);
SEXP lreg_c(SEXP, SEXP, SEXP, SEXP);
SEXP mv_c(SEXP, SEXP);
SEXP tmv_c(SEXP, SEXP);
