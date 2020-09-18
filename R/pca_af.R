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
#' allele_freqs_pca = pca_af(hgdp_subset, 4)
#' summary(abs(allele_freqs_lfa-allele_freqs_pca))
#' @export
pca_af <- function(X, d, override=FALSE){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( is.null( d ) )
        stop( "Principal components number `d` is required!" )
    
    m <- nrow(X)
    n <- ncol(X)
    
    #check for d validity
    if(d != as.integer(d)){
        stop("d should be integer")
    } else if(d < 1){
        stop("d should be at least 1")
    } else if(d >= 1){
        d <- d-1 #for the svd stuff
    }
    
    adjust <- 8
    if( n-d < 10 )
        adjust <- n-d-1
    
    # index the missing values
    NA_IND <- is.na(X)
    
    # center the matrix...
    mean_X <- rowMeans(X, na.rm=TRUE)

    if ( d == 0 ) {
        # this is "intercept only"
        # all allele frequencies are just the mean
        AF <- matrix(
            rep.int( mean_X, n ),
            nrow = m,
            ncol = n
        )
        return( AF )
    }
    # else continue
    
    norm_X <- X - mean_X
    # sd_X <- apply(X, 1, function(snp){sd(snp, na.rm=TRUE)})
    
    #...then 'impute'
    norm_X[NA_IND] <- 0
    
    mysvd <- trunc_svd(norm_X, d=d, adjust=adjust, tol=1e-13, override=override)
    
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

    # WHY???
    #AF[NA_IND] <- NA
    
    return(AF)
}
