#include "lfa.h"

//logistic regression
//you MUST add the constant before this
SEXP lreg(SEXP RX, SEXP Ry, SEXP Rmi, SEXP Rtol){
    int *dimX, maxiter = (int)(*REAL(Rmi));
    double *X, *y, tol = (double)(*REAL(Rtol));
    int i, j, k, ind1, ind2, ind;
    int flag;

    dimX = getDims(RX);
    PROTECT(Ry = coerceVector(Ry, REALSXP));
    PROTECT(RX = coerceVector(RX, REALSXP));
    y = REAL(Ry);
    X = REAL(RX);
    
    SEXP Rret;
    double *ret;
    PROTECT(Rret = allocVector(REALSXP, dimX[1]));
    ret = REAL(Rret);
    
    int numblock = 64*dimX[1];
    double *b  = (double*) malloc(sizeof(double)*dimX[1]); //beta
    double *bl = (double*) malloc(sizeof(double)*dimX[1]); //beta last
    double *f  = (double*) malloc(sizeof(double)*dimX[1]); //tmp
    double *p  = (double*) malloc(sizeof(double)*dimX[0]); //mle

    double *w  = (double*) malloc(sizeof(double)*dimX[1]*dimX[1]);
    int* ipiv  = (int*)    malloc(sizeof(int)   *dimX[1]); //for inverting
    double *wo = (double*) malloc(sizeof(double)*numblock);
    double max; // check convergence

    int iter = 1;
    double alpha = -1.0, zero = 0.0, one = 1.0;
    int ione = 1;
    int info=0;
    char tr = 'n';
    double tmp;
    for(i = 0; i < dimX[1]; i++) {
        b[i]  = 0;
        bl[i] = 0;
    }

    //IRLS
    flag = 0;
    while(iter <= maxiter){
        ///////////////////////////////////////////////////////////////////////
        //p <- as.vector(1/(1 + exp(-X %*% b)))
        F77_CALL(dgemv)(&tr,dimX,dimX+1,&alpha,X,dimX,b,&ione,&zero,p,&ione);
        for(i = 0; i < dimX[0]; i++) 
            p[i] = 1/(1+exp(p[i]));
        
        ///////////////////////////////////////////////////////////////////////
        //var.b <- solve(crossprod(X, p * (1 - p) * X))
        //
        //here, solve is inverting the matrix. 
        //p*(1-p) is applied to cols of X.
        //at the moment I am manually computing the crossprod
        //which is guaranteed to be symmetric
        for(i = 0; i < dimX[1]; i++){ //rows
            for(j = i; j < dimX[1]; j++){ //columns
                ind1 = i*dimX[0]; //i-th col of X
                ind2 = j*dimX[0]; //j-th col of X
                ind  = dimX[1]*i + j; //position on w
                w[ind] = 0;
                for(k = 0; k < dimX[0]; k++){ //loop over X'p(1-p)X
                    w[ind]+=X[ind1]*X[ind2]*p[k]*(1-p[k]);
                    ind1++;
                    ind2++;
                }
                if(i != j) //reflect it
                    w[dimX[1]*j+i] = w[ind];
            }
        }
        
        //actually inverting here. remember to pay attention to includes
        F77_CALL(dgetrf)(dimX+1,dimX+1,w,dimX+1,ipiv,&info);
        if(info != 0) {
            //Rprintf("warning: dgetrf error, NA used\n");
            //Rprintf("info:%i iter:%i\n", info, iter);
            //error("dgetrf error\n");
            flag = 1;
        }
        F77_CALL(dgetri)(dimX+1,w,dimX+1,ipiv,wo,&numblock,&info);
        if(info != 0) {
            //Rprintf("warning: dgetri error, NA used\n");
            //Rprintf("info:%i iter:%i\n", info, iter);
            //error("dgetri error\n");
            flag = 1;
        }
        
        //if a failure, skip outta here.
        if(flag == 1){
            for(i = 0; i < dimX[1]; i++) ret[i] = R_NaReal;
            free(b);
            free(bl);
            free(f);
            free(p);
            free(w);
            free(ipiv);
            free(wo);
            UNPROTECT(3);
            return Rret;
        }


        ///////////////////////////////////////////////////////////////////////
        //b <- b + var.b %*% crossprod(X, y - p)
        //use f to calculate crossprod(X,y-p) first.
        //then use dgemv
        ind  = 0; //since we are iterating over X in order
        for(i = 0; i < dimX[1]; i++){ //cols of X, values of f
            f[i] = 0;
            for(j = 0; j < dimX[0]; j++){ //rows of X, values of y-p
                f[i] += X[ind] * (y[j] - p[j]);
                ind++;
            }
        }

        F77_CALL(dgemv)(&tr,dimX+1,dimX+1,&one,w,dimX+1,f,&ione,&one,b,&ione);
        
        
        ///////////////////////////////////////////////////////////////////////
        //if (max(abs(b - b.last)/(abs(b.last) + 0.01*tol)) < tol) break
        //check to see if we need to break
        max = 0.0;
        for(i = 0; i < dimX[1]; i++) {
            tmp = fabs(b[i] - bl[i])/(fabs(bl[i]) + 0.01*tol);
            if(tmp > max) max = tmp;
        }

        if(max < tol)
            break;

        
        ///////////////////////////////////////////////////////////////////////
        //b.last <- b
        //it <- it + 1
        for(i = 0; i < dimX[1]; i++) bl[i] = b[i];
        
        iter++;
    }

    //if(iter > maxiter) printf("warning: max iterations exceeded\n");
    
    //set the return...
    for(i = 0; i < dimX[1]; i++) ret[i] = b[i];
    
    free(b);
    free(bl);
    free(f);
    free(p);
    free(w);
    free(ipiv);
    free(wo);
    UNPROTECT(3);
    return Rret;
}


