#include "lfa.h"

SEXP mv(SEXP RA, SEXP Rv){
    int *dimA;
    double *v, *A;

    dimA = getDims(RA);
    PROTECT(RA=coerceVector(RA, REALSXP));
    PROTECT(Rv=coerceVector(Rv, REALSXP));
    A = REAL(RA);
    v = REAL(Rv);

    SEXP Rret;
    double *ret;
    PROTECT(Rret = allocVector(REALSXP, dimA[0]));
    ret = REAL(Rret);

    double alpha = 1.0;
    double zero = 0.0;
    char tr = 'N';
    int one = 1;
    F77_CALL(dgemv)(&tr,dimA,dimA+1,&alpha,A,dimA,v,&one,&zero,ret,&one);

    UNPROTECT(3);

    return Rret;
}

SEXP tmv(SEXP RA, SEXP Rv){
    int *dimA;
    double *v, *A;

    dimA = getDims(RA);
    PROTECT(RA=coerceVector(RA, REALSXP));
    PROTECT(Rv=coerceVector(Rv, REALSXP));
    A = REAL(RA);
    v = REAL(Rv);

    SEXP Rret;
    double *ret;
    PROTECT(Rret = allocVector(REALSXP, dimA[1]));
    ret = REAL(Rret);

    double alpha = 1.0;
    double zero = 0.0;
    char tr = 'T';
    int one = 1;
    F77_CALL(dgemv)(&tr,dimA,dimA+1,&alpha,A,dimA,v,&one,&zero,ret,&one);

    UNPROTECT(3);

    return Rret;
}
