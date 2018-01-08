#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.


#' @title Logistic factor analysis
#'
#' @description
#' Fit a factor model of dimension \eqn{d} for binomial data. Returns 
#' logistic factors.
#'
#' @details
#' This function performs logistic factor analysis on SNP data. As it
#' stands, we follow the convention where \eqn{d=1} is intercept only,
#' and for \eqn{d>1} we compute \eqn{d-1} singular vectors and postpend
#' the intercept.
#'
#' @note Genotype matrix is expected to be a matrix of integers with
#' values 0, 1, and 2. Note
#' that the coding of the SNPs does not affect the algorithm.
#'
#' @param X a matrix of SNP genotypes, i.e. an integer matrix of 0's,
#' 1's, and 2's. Sparse matrices of class Matrix are not supported
#' (yet).
#' @param d number of logistic factors, including the intercept
#' @param adjustments a matrix of adjustment variables to hold fixed 
#' during estimation. 
#' @param override optional boolean to bypass Lanczos bidiagonalization
#' SVD. Usually not advised unless encountering a bug in the SVD code.
#' @param safety optional boolean to bypass checks on the genotype
#' matrices, which require a non-trivial amount of computation.
#' @return matrix of logistic factors, with the intercept at the end.
#' @export
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' dim(LF)
#' head(LF)
#' @importFrom corpcor fast.svd
#' @useDynLib lfa
lfa <- function(X, d, adjustments=NULL, override=FALSE, safety=FALSE){
    if(safety)
        check.geno(X)
    
    m <- nrow(X)
    n <- ncol(X)

    # check for d validity
    if(d != as.integer(d)){
        stop("d should be integer")
    } else if(d < 1){
        stop("d should be at least 1")
    } else if(d == 1){
        return(matrix(1, n, 1))
    } else if(d >1){
        d <- d-1 #for the svd stuff
    }
    
    # check adjustments vars
    if(!is.null(adjustments)){
        if (nrow(adjustments) != ncol(X)){
          stop("adjustments needs to have same number of rows as individuals")
        }
        if (ncol(adjustments) >= d){
          stop("need to estimate at least one non-adjustment logistic factor")
        }
        if (sum(!complete.cases(adjustments)) > 0){
          stop("no missing values in adjustments")
        }
    }
    
    adjust <- 8
    if((n-d ) < 10) adjust <- n-d-1
    
    # index the missing values
    NA_IND <- is.na(X)
    
    #center the matrix...
    mean_X <- rowMeans(X, na.rm=TRUE)
    norm_X <- X - mean_X

    # ...then 'impute'
    norm_X[NA_IND] <- 0
    
    # first SVD
    mysvd <- trunc.svd(norm_X, d=d, adjust=adjust, tol=1e-13, override=override)

    rm(norm_X)
    D <- diag(mysvd$d, d, d)
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    
    # form projection
    z <- U %*% D %*% t(V)

    z <- z + mean_X
    z <- z/2
    rm(U); rm(D); rm(V)

    #The .Call() is equivalent to the following lines of R code:
    #
    #zmin <- apply(z, 1, min)
    #zmax <- apply(z, 1, max)
    #ind  <- (zmax<(1-2/n)) & (zmin>(2/n))
    ind <- as.logical(.Call("lfa_threshold", z, 1/(2*n)))
    z <- z[ind,]
    z <- log(z/(1-z))

    norm_z <- centerscale(z)
    
    # regress out adjustment vars, if relevant
    if(!is.null(adjustments)){
      norm_z <- t(residuals(lm(t(norm_z)~adjustments-1)))
      d <- d - ncol(adjustments)
    }

    v <- trunc.svd(norm_z, d=d, adjust=adjust, tol=1e-13, override=override)$v
    v <- cbind(v,1)
    if(!is.null(adjustments)){
      v <- cbind(adjustments, v)
    }
    return(v)
}

#' @title Matrix centering and scaling
#'
#' @description
#' C routine to row-center and scale a matrix. Doesn't work with missing data.
#'
#' @param A matrix
#' @examples
#' centerscale(hgdp_subset)
#' @return matrix same dimensions \code{A} but row centered and scaled
#' @export
centerscale <- function(A){
    as.matrix(.Call("centerscale", A))
}

#' @title Matrix centering
#'
#' @description
#' C routine to row-center a matrix
#'
#' @param A matrix
#' @examples
#' center(hgdp_subset)
#' @return matrix same dimensions \code{A} but row centered
#' @export
center <- function(A){
    as.matrix(.Call("center", A))
}

# returns T/F for missing values or not
check.geno <- function(X){
    ret <- FALSE
    if(class(X) != "matrix")
        stop("The input must be genotypes in a matrix class.")

    if(class(X[1]) != "integer")
        stop("Elements of the genotype matrix should be integer.")

    classes <- names(table(as.vector(X)))
    if(length(classes) != 3)
        stop("Expecting genotypes to be 0, 1, and 2.")
    if(classes != c("0", "1", "2"))
        stop("Expecting genotypes to be 0, 1, and 2.")

    unique.length <- apply(X,1,function(x) length(unique(x)))
    if(sum(unique.length==1) > 1) {
        stop(paste0("Remove ", unique.length," loci without any variation across samples."))
    }

    m <- nrow(X)
    n <- ncol(X)

    if(m <= n)
        stop("The genotype matrix should be tall.")

}

#' @title Hardy-Weinberg Equilibrium in structure populations
#' 
#' @description
#' Compute structural Hardy-Weinberg Equilibrium (sHWE) p-values
#' on a SNP-by-SNP basis. These p-values can be aggregated to 
#' determine genome-wide goodness-of-fit for a particular value
#' of \eqn{d}. See \url{https://doi.org/10.1101/240804} for more
#' details.
#'
#' @param LF matrix of logistic factors
#' @param B number of null datasets to generate - \eqn{B=1} is usualy
#' sufficient. If computational time/power allows, a few extra
#' \eqn{B} could be helpful
#' @inheritParams lfa
#' @inheritParams sHWE
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' gof_4 <- sHWE(hgdp_subset, LF, 3)
#' LF <- lfa(hgdp_subset, 10)
#' gof_10 <- sHWE(hgdp_subset, LF, 3)
#' hist(gof_4)
#' hist(gof_10)
#' @return a vector of p-values for each SNP.
#' @export
sHWE <- function(X, LF, B){
    obs_stat <- apply(X, 1, gof.stat.snp, LF)
    d <- ncol(LF)
    AF <- af(X,LF)
    rm(X)

    stat0 <- compute.nulls(AF, d, B)

    obs_stat_order <- order(obs_stat)
    obs_stat <- sort(obs_stat)
    stat0 <- sort(as.vector(stat0))

    p <- rep(0, length(obs_stat))
    m <- length(obs_stat)
    B0 <- length(stat0)
    i1 <- 1
    i2 <- 1

    while(i1 <= m) {
        while((i2 <= B0) & (obs_stat[i1] >= stat0[i2])) {
            i2 <- i2+1
        }
        p[obs_stat_order[i1]] <- 1 - ((i2-1)/B0)
        i1 <- i1+1
    }

    p
}

#' @title LFA model goodness of fit
#' 
#' @description
#' Compute SNP-by-SNP goodness-of-fit when compared to population 
#' structure. This can be aggregated to determine genome-wide 
#' goodness-of-fit for a particular value of \eqn{d}.
#' 
#' @details
#' This function returns p-values for LFA model goodness of fit based
#' on a simulated null.
#'
#' @note Genotype matrix is expected to be a matrix of integers with
#' values 0, 1, and 2. Currently no support for missing values. Note
#' that the coding of the SNPs does not affect the algorithm.
#'
#' @inheritParams lfa
#' @inheritParams sHWE
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' gof_4 <- model.gof(hgdp_subset, LF, 3)
#' LF <- lfa(hgdp_subset, 10)
#' gof_10 <- model.gof(hgdp_subset, LF, 3)
#' hist(gof_4)
#' hist(gof_10)
#' @return vector of p-values for each SNP.
#' @export
model.gof <- sHWE

inverse_2x2 <- function(X) {
    denom <- X[1]*X[4]-X[2]*X[3]
    stopifnot(denom != 0)
    matrix( c(X[4], -X[2],
              -X[3], X[1]), 2, 2) / denom
}

gof.stat.snp <- function(snp, LF){
    NA_IND <- is.na(snp)
    snp <- snp[!is.na(snp)]
    LF  <- LF[!NA_IND, ,drop=FALSE]

    p <- af_snp(snp, LF)

    p0 <- (1-p)^2
    p1 <- 2*p*(1-p)

    est <- c(sum(p0), sum(p1))
    N   <- c(sum( snp==0 ), sum( snp==1 ))

    SIGMA <- matrix( c(sum(p0*(1-p0)), -1*sum(p0*p1),
                      -1*sum(p0*p1), sum(p1*(1-p1))), 2, 2)

    stat <- t(N-est) %*% inverse_2x2(SIGMA) %*%(N-est)

    return(stat)
}

compute.nulls <- function(AF, d, B) {
    m <- nrow(AF)
    n <- ncol(AF)

    stat0 <- matrix(0, m, B)
    for(i in 1:B) {
        suppressWarnings(X0 <- matrix(rbinom(m*n, 2, as.numeric(AF)), m, n))
        LF <- lfa(X0, d)
        stat0[,i] <- apply(X0, 1, gof.stat.snp, LF)
    }

    return(stat0)
}


