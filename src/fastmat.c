#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP))


SEXP mv(SEXP RA, SEXP Rv){
    int *dimA, *dimv;
    double *v, *A;

    dimA = getDims(RA); //266118 443
    dimv = getDims(Rv); //443    1
    RA=coerceVector(RA, REALSXP);
    Rv=coerceVector(Rv, REALSXP);
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

    UNPROTECT(1);

    return Rret;
}

SEXP tmv(SEXP RA, SEXP Rv){
    int *dimA, *dimv;
    double *v, *A;

    dimA = getDims(RA); //266118 443
    dimv = getDims(Rv); //443    1
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


SEXP minv(SEXP RA){
    int *dimA;
    double *A;
    
    dimA = getDims(RA);
    if(dimA[0] != dimA[1]) error("expecting sq matrix\n");
    PROTECT(RA=coerceVector(RA, REALSXP));
    A = REAL(RA);

    SEXP Rret = PROTECT(duplicate(RA));
    double *ret = REAL(Rret);
    
    int dimsq;
    if(dimA[0] < 1000)
        dimsq = dimA[0]*dimA[0];
    else if(dimA[0] < 10000)
        dimsq = dimA[0]*2;
    else 
        dimsq = dimA[0];
    int* ipiv = (int*) malloc(sizeof(int)*dimA[0]);
    double *wo = (double*) malloc(sizeof(double)*dimsq);
    int  info = 0;
    
    F77_CALL(dgetrf)(dimA,dimA,ret,dimA,ipiv,&info);
    if(info != 0) error("dgetrf error\n");
    F77_CALL(dgetri)(dimA,ret,dimA,ipiv,wo,&dimsq,&info);
    if(info != 0) error("dgetri error\n");
    
    free(wo);
    free(ipiv);
    UNPROTECT(2);
    return Rret;
}


/* So, given that we have a tall matrix X, we construct a small
 * symmetrized square version of it with X'X. Then we compute the SVD 
 * X'X, calling the components U*, D*, and V*. To recover the SVD of X,
 * components called U, D, and V, we have to do D=sqrt(D*), V=V*, and
 * U = XVD^{-1}. This function performs the last matrix multiplication,
 * taking advantage of memory saving tricks in dgemm.
 * 
 * Expected inputs are X, V, and 1/D. 
 * 
 * Output is U (m by n)
 */
SEXP tall_svd_mm(SEXP RX, SEXP RV, SEXP RD){
    int *dimX, *dimV, *dimD;
    double *X,    *V,    *D;
    
    //get dimensions
    dimX = getDims(RX); // m by n
    dimV = getDims(RV); // n by n
    dimD = getDims(RD); // n by n
    
    //check dimensions
    if(dimV[1] != dimD[0] || dimX[1] != dimV[0])
        error("dimension issue\n");
    
    //do the rest of the R stuff
    PROTECT(RX = coerceVector(RX, REALSXP));
    PROTECT(RV = coerceVector(RV, REALSXP));
    PROTECT(RD = coerceVector(RD, REALSXP));
    X = REAL(RX);
    V = REAL(RV);
    D = REAL(RD);
    
    SEXP RU; //return value
    double *U;
    PROTECT(RU = allocMatrix(REALSXP, dimX[0], dimD[1]));
    U = REAL(RU);
    
    double zero = 0.0, one = 1.0;
    char tr = 'n';
    
    //tmp square matrix that stores V 1/D
    double *T = (double*) malloc(sizeof(double)*dimV[0]*dimD[1]); 
    
    //OK, let's do T = V 1/D
    F77_CALL(dgemm)(&tr,&tr,dimV,dimD+1,dimV+1,&one,V,dimV,D,dimD,&zero,T,dimV);
    
    //now let's do X times the earlier result T
    F77_CALL(dgemm)(&tr,&tr,dimX,dimD+1,dimX+1,&one,X,dimX,T,dimV,&zero,U,dimX);
    
    //OK we're done
    UNPROTECT(4);
    free(T);
    return RU;
}


/* This function does t(X) %*% X in really fast manner using a symmetric
 * update from BLAS
 */
SEXP crossprod(SEXP RX){
    int *dimX, m, n, i, j, ind1, ind2;
    double *X;
    
    dimX = getDims(RX); // m by n
    m = dimX[0];
    n = dimX[1];
    
    PROTECT(RX = coerceVector(RX, REALSXP));
    X = REAL(RX);
    
    SEXP RU; //return value
    double *U;
    PROTECT(RU = allocMatrix(REALSXP, n, n));
    U = REAL(RU);
    
    char tr = 't', uplo = 'u';
    double zero = 0.0, one = 1.0;
    
    F77_CALL(dsyrk)(&uplo,&tr,&n,&m,&one,X,&m,&zero,U,&n);

    //symmetize U, column major format
    for(i = 0; i < n; i++){
        for(j = 0; j < i; j++){
            ind1 = j*n + i; // i-th row
            ind2 = i*n + j;
            U[ind1] = U[ind2];
        }
    }
    
    UNPROTECT(2);
    return RU;
}
