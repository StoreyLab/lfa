#include "lfa.h"

SEXP lfa_threshold(SEXP RX, SEXP Rthresh){
    int *dimX, n, i, ind;
    double *X, max, min;
    double thresh = (double)(*REAL(Rthresh));
    
    dimX = getDims(RX);
    PROTECT(RX = coerceVector(RX, REALSXP));
    X = REAL(RX);

    if(dimX[1] <1) 
        Rprintf("dimension problem in lfa_threshold...");

    SEXP Rret; //returns boolean list of valid rows
    double *ret;
    PROTECT(Rret = allocVector(REALSXP, dimX[0]));
    ret = REAL(Rret);

    for(n = 0; n < dimX[0]; n++){
        min = X[n]; //set min/max to first element
        max = X[n];
        ind = n + dimX[0]; //start from second element
        for(i = 1; i < dimX[1]; i++){
            if(X[ind] > max)
                max = X[ind];
            else if(X[ind] < min)
                min = X[ind];
            ind += dimX[0]; //iterate across loops of course
        }
        //Rprintf("%f %f\n", max, min);
        if((max < (1-thresh)) && (min > thresh))
            ret[n] = 1;
        else
            ret[n] = 0;
    }
    
    UNPROTECT(2);
    return Rret;
}


//This function seeks to do the following lines of R code:
//    mean_x = apply(x,1,mean)
//    sd_x = apply(x,1,sd)
//    z = (z*sd_x) + mean_x
//    z = z/2
//except be really efficient by taking full advantage of passing by
//reference.
SEXP lfa_scaling(SEXP RX, SEXP RZ){
    int *dimX, n, i, ind;
    double *X, *Z, mean, sd;
    
    dimX = getDims(RX);
    PROTECT(RX = coerceVector(RX, REALSXP));
    X = REAL(RX);
    
    PROTECT(RZ = coerceVector(RZ, REALSXP));
    Z = REAL(RZ);
    
    for(n = 0; n < dimX[0]; n++){
        mean = 0;
        sd = 0;
        
        ind = n;
        for(i = 0; i < dimX[1]; i++){
            mean += X[ind];
            ind += dimX[0]; //looping over rows...
        }
        mean = mean/dimX[1];
        
        ind=n;
        for(i = 0; i < dimX[1]; i++){
            Z[ind] *= sd;
            Z[ind] += mean;
            Z[ind] /= 2;
            ind += dimX[0]; //looping over rows...
        }
    }
    
    UNPROTECT(2);
    return R_NilValue;
}


//two utility functions for centerscale
double sd(double* A, int n, int inc){
    int i, ind=0;
    double sum = 0;
    for(i = 0; i < n; i++){
        sum += A[ind];
        ind += inc;
    }
    
    double mean = sum/n;
    sum = 0;
    ind = 0;
    
    for(i = 0; i < n; i++) {
        sum += (A[ind]-mean) * (A[ind]-mean);
        ind += inc;
    }
    
    return sqrt(sum/(n-1));
}

double mean(double* A, int n, int inc){
    int i, ind = 0;
    double sum = 0;
    
    for(i = 0; i < n; i++){
        sum += A[ind];
        ind += inc;
    }
    
    return sum/n;    
}

SEXP centerscale(SEXP RA){
    int *dimA;
    double *A;
    
    dimA = getDims(RA);
    if(dimA[0] <= 1) error("er, first dimension is 1? that's weird.");
    if(dimA[1] <= 1) error("er, second dimension is 1? that's weird.");
    PROTECT(RA=coerceVector(RA, REALSXP));
    A = REAL(RA);

    SEXP Rret = PROTECT(duplicate(RA));
    double *ret = REAL(Rret);
    
    int i, j, ind;
    double m, s;
    for(i = 0; i < dimA[0]; i++){
        ind = i;
        m = mean(A+i, dimA[1], dimA[0]);
        s = sd(A+i, dimA[1], dimA[0]);
        
        for(j = 0; j < dimA[1]; j++){
            if (s != 0) {
                ret[ind] = (A[ind] - m)/s;
                ind += dimA[0];
            }
            else {
                ret[ind] = 0;
                ind += dimA[0];
            }
        }
    }
    
    UNPROTECT(2);
    return Rret;
}

SEXP center(SEXP RA){
    int *dimA;
    double *A;
    
    dimA = getDims(RA);
    if(dimA[0] <= 1) error("er, first dimension is 1? that's weird.");
    if(dimA[1] <= 1) error("er, second dimension is 1? that's weird.");
    PROTECT(RA=coerceVector(RA, REALSXP));
    A = REAL(RA);

    SEXP Rret = PROTECT(duplicate(RA));
    double *ret = REAL(Rret);
    
    int i, j, ind;
    double m;
    for(i = 0; i < dimA[0]; i++){
        ind = i;
        m = mean(A+i, dimA[1], dimA[0]);
        
/*        if( i ==0) {
            printf("%f %f\n", m, s);
        }*/
        
        for(j = 0; j < dimA[1]; j++){
            ret[ind] = A[ind] - m;
            ind += dimA[0];
        }
    }
    
    UNPROTECT(2);
    return Rret;
}

