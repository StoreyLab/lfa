#' @title Allele frequencies
#' @description Compute matrix of individual-specific allele frequencies
#' @inheritParams lfa
#' @param LF Matrix of logistic factors, with intercept. Pass in the 
#' return value from \code{lfa}!
#' @details Computes the matrix of individual-specific allele 
#' frequencies, which has the same dimensions of the genotype matrix.
#' Be warned that this function could use a ton of memory, as the 
#' return value is all doubles. It could be wise to pass only a 
#' selection of the SNPs in your genotype matrix to get an idea for
#' memory usage. Use \code{gc} to check memory usage!
#' @return Matrix of individual-specific allele frequencies.
#' @export
af = function(X, LF, safety=FALSE){
    if(safety)
        check.geno(X)
    t(apply(X, 1, af_snp, LF))
}

#' @title Allele frequencies for SNP
#' @description Computes individual-specific allele frequencies for a 
#' single SNP.
#' @return vector of allele frequencies
af_snp = function(snp, LF){
    #b_0 = clreg(snp, LF)
    #est = .Call("mv", cbind(1, LF), b_0)
    b_0    = lreg(snp,LF) #coefficients from logreg assuming HWE
    est    = .Call("mv", LF, b_0)
    
    exp(est)/(1+exp(est))
}

#C based logistic regression 
lreg <- function(y,X){
    if(is.null(X) || !hasIntercept(X)){
        stop("you must supply the intercept!")
    }
    X = rbind(X, X)
    y1 = as.numeric((y==1) | (y==2))
    y2 = as.numeric(y==2)
    y=c(y1,y2)
    .Call("lreg", X, y, 20, 1e-10)
}

hasIntercept <- function(x) {
    return(TRUE) #3/6/2014
    counts = apply(x, 2, function(y) {length(unique(y))})
    counts = (counts == 1)
    if(sum(counts) != 1) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

