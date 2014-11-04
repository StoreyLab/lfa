#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.


#' @title Logistic factor analysis.
#' 
#' @details
#' This function performs logistic factor analysis on SNP data. As it 
#' stands, we follow the convention where \eqn{d=1} is intercept only,
#' and for \eqn{d>1} we compute \eqn{d-1} singular vectors and postpend 
#' the intercept. 
#'
#' @note Genotype matrix is expected to be a matrix of integers with 
#' values 0, 1, and 2. Currently no support for missing values. Note 
#' that the coding of the SNPs does not affect the algorithm.
#'
#' @param X a matrix of SNP genotypes, i.e. an integer matrix of 0's, 
#' 1's, and 2's. Sparse matrices of class Matrix are not supported 
#' (yet).
#' @param d number of logistic factors, including the intercept
#' @param override optional boolean to bypass Lanczos bidiagonalization
#' SVD. Usually not advised unless encountering a bug in the SVD code.
#' @param safety optional boolean to bypass checks on the genotype 
#' matrices, which require a non-trivial amount of computation.
#' @return matrix of logistic factors, with the intercept at the end.
#' @export
#' @useDynLib lfa
lfa = function(X, d, override=FALSE, safety=FALSE){
    if(safety)
        check.geno(X)
    
    m = nrow(X)
    n = ncol(X)

    #check for d validity
    if(d != as.integer(d)){
        stop("d should be integer")
    } else if(d < 1){
        stop("d should be at least 1")
    } else if(d == 1){
        return(matrix(1, n, 1))
    } else if(d >1){
        d = d-1 #for the svd stuff
    }
    
    adjust = 8
    if((n-d ) < 10) adjust=n-d-1
    
    norm_X = center(X)
    mysvd = trunc.svd(norm_X, d=d, adjust=adjust, tol=1e-13, override=override)
    
    mean_X = apply(X,1,mean)
    sd_X = apply(X,1,sd)
    
    rm(norm_X)
    D = mysvd$d
    U = mysvd$u
    V = mysvd$v
    rm(mysvd)
    
    z = U %*% diag(D, d, d)  %*% t(V)
        
    #z = (z*sd_X) + mean_X
    z = z + mean_X
    z = z/2
    rm(U); rm(D); rm(V)

    #The .Call() is equivalent to the following lines of R code:
    #
    #zmin = apply(z, 1, min)
    #zmax = apply(z, 1, max)
    #ind  = (zmax<(1-2/n)) & (zmin>(2/n))
    ind = as.logical(.Call("lfa_threshold", z, 1/(2*n)))
    z = z[ind,]
    z = log(z/(1-z))

    norm_z = centerscale(z)
    v = trunc.svd(norm_z, d=d, adjust=adjust, tol=1e-13, override=override)$v
    v = cbind(v,1)
    return(v)
}

#' @title Matrix centering and scaling
#' 
#' @description
#' C routine to row-center and scale a matrix
#'
#' @param A matrix
#' @return matrix same dimensions \code{A} but row centered and scaled
centerscale <- function(A){
    as.matrix(.Call("centerscale", A))
}

#' @title Matrix centering
#' 
#' @description
#' C routine to row-center a matrix
#'
#' @param A matrix
#' @return matrix same dimensions \code{A} but row centered
center <- function(A){
    as.matrix(.Call("center", A))
}


check.geno <- function(X){
    if(class(X) != "matrix")
        stop("expecting genotypes in class matrix")
    
    if(class(X[1]) != "integer")
        stop("elements of genotype matrix should be integer")
    
    classes = names(table(as.vector(X)))
    if(length(classes) != 3)
        stop("expecting genotypes to be 0, 1, and 2 (no missing values)")
    if(class != c("0", "1", "2"))
        stop("expecting genotypes to be 0, 1, and 2")
    
    m = nrow(X)
    n = ncol(X)
    
    if(m <=n )
        stop("genotype matrix shoudl be tall")
    
}

check.LF <- function(LF){
    
}
