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
#' @examples
#' LF = lfa(hgdp_subset, 4)
#' allele_freqs = af(hgdp_subset, LF)
#' @return Matrix of individual-specific allele frequencies.
#' @export
af <- function(X, LF, safety=FALSE){
    if(safety)
        check.geno(X)
    t(apply(X, 1, af_snp, LF))
}

#' @title Allele frequencies for SNP
#' @description Computes individual-specific allele frequencies for a 
#' single SNP.
#' @inheritParams af
#' @param snp vector of 0's, 1's, and 2's
#' @return vector of allele frequencies
#' @export
af_snp <- function(snp, LF){
    NA_IND <- is.na(snp)
    af <- rep(NA, length(snp))
    
    snp <- snp[!NA_IND]
    LF  <- LF[!NA_IND, , drop=FALSE]
    b_0 <- lreg(snp, LF) #coefficients from logreg 
    est <- .Call("mv", LF, b_0)
    
    af[!NA_IND] <- exp(est)/(1+exp(est))
    af
}

#C based logistic regression 
lreg <- function(y, X){
    if(is.null(X) || !hasIntercept(X)){
        stop("you must supply the intercept!")
    }
    X <- rbind(X, X)
    y1 <- as.numeric((y==1) | (y==2))
    y2 <- as.numeric(y==2)
    y <- c(y1,y2)
    b <- .Call("lreg", X, y, 20, 1e-10)
    
    #if coefficients are NA, use glm
    if(sum(is.na(b)) > 0) {
        # print(paste("harmless warning:"))
        b <- glm(cbind(y, 2-y) ~ -1 + X, family="binomial")$coef
        names(b) <- NULL
    }
    b
}

hasIntercept <- function(x) {
    return(TRUE) #3/6/2014
    counts <- apply(x, 2, function(y) {length(unique(y))})
    counts <- (counts == 1)
    if(sum(counts) != 1) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' @title PCA Allele frequencies
#' @description Compute matrix of individual-specific allele frequencies
#' via PCA
#' @inheritParams lfa
#' @details This corresponds to algorithm 1 in the paper. Only used for
#' comparison purposes.
#' @return Matrix of individual-specific allele frequencies.
#' @examples
#' LF = lfa(hgdp_subset, 4)
#' allele_freqs_lfa = af(hgdp_subset, LF)
#' allele_freqs_pca = pca_af(hgdp_subset, 4, LF)
#' summary(abs(allele_freqs_lfa-allele_freqs_pca))
#' @export
pca_af <- function(X, d, override=FALSE){
    m <- nrow(X)
    n <- ncol(X)
    
    #check for d validity
    if(d != as.integer(d)){
        stop("d should be integer")
    } else if(d < 1){
        stop("d should be at least 1")
    } else if(d == 1){
        return(matrix(1, n, 1))
    } else if(d >1){
        d <- d-1 #for the svd stuff
    }
    
    adjust <- 8
    if((n-d ) < 10) adjust=n-d-1

    #index the missing values
    NA_IND <- is.na(X)
    
    #center the matrix...
    mean_X <- rowMeans(X, na.rm=TRUE)
    norm_X <- X - mean_X
    #sd_X <- apply(X, 1, function(snp){sd(snp, na.rm=TRUE)})
    
    #...then 'impute'
    norm_X[NA_IND] <- 0
    
    mysvd <- trunc.svd(norm_X, d=d, adjust=adjust, tol=1e-13, override=override)
    
    rm(norm_X)
    D <- mysvd$d
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    
    z <- U %*% diag(D, d, d)  %*% t(V)
        
    z <- z + mean_X
    z <- z/2
    
    AF <- t(apply(z, 1, function(x) {
            IND <- x > 1 - 1/(2*n)
            x[IND] <- 1 - 1/(2*n)
            IND <- x < 1/(2*n)
            x[IND] <- 1/(2*n)
            x
           }))
    
    rm(z)
    
    AF[NA_IND] <- NA
    return(AF)
}
