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
pca_af <- function(
                   X,
                   d,
                   override = FALSE,
                   ploidy = 2,
                   tol = 1e-13,
                   m_chunk = 1000 # gave good performance in tests
                   ){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( is.null( d ) )
        stop( "Principal components number `d` is required!" )
    
    # check class
    is_BEDMatrix <- FALSE
    if ( "BEDMatrix" %in% class(X) ) {
        is_BEDMatrix <- TRUE
    } else if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    
    #check for d validity
    if(d != as.integer(d)){
        stop("d should be integer")
    } else if(d < 1){
        stop("d should be at least 1")
    } else if(d >= 1){
        d <- d-1 #for the svd stuff
    }

    if ( is_BEDMatrix ) {
        # data dimensions (transposed for BEDMatrix)
        m <- ncol(X)
        n <- nrow(X)

        # NOTE: this is inefficient for the edge case d=0 (no PCs, just mean), but that's a dumb case nobody should actually want in practice.
        # calculate covariance matrix and loci means
        obj <- covar_BEDMatrix( X, m_chunk = m_chunk )
        covar <- obj$covar
        X_mean <- obj$X_mean
        
        if ( d == 0 ) {
            # this is "intercept only"
            # all allele frequencies are just the mean
            P <- matrix(
                rep.int( X_mean, n ),
                nrow = m,
                ncol = n
            )
            return( P )
        }
        # else continue
        
        # get truncated eigendecomposition
        obj <- RSpectra::eigs_sym( covar, d )
        V <- obj$vectors
        # turn eigenvector matrix into a projection matrix
        P_V <- tcrossprod( V )
        
        # most painful part is forming P in parts, so X is not in memory all at once.
        # P will be in memory at once though, can be too much memory usage in some cases and even negate the BEDMatrix advantage
        # initialize desired matrix and vector
        P <- matrix( 0, nrow = m, ncol = n )
        
        # navigate chunks
        i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
        while (TRUE) { # start an infinite loop, break inside as needed
            # here m is known...

            # this means all SNPs have been covered!
            if (i_chunk > m)
                break

            # range of SNPs to extract in this chunk
            indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m)

            # transpose for our usual setup (makes centering easiest)
            Xi <- t( X[, indexes_loci_chunk, drop = FALSE] )

            # update starting point for next chunk! (overshoots at the end, that's ok)
            i_chunk <- i_chunk + m_chunk

            # get row means from precomputed data
            Xi_mean <- X_mean[ indexes_loci_chunk ]
            
            # center data
            Xi <- Xi - Xi_mean
            
            # set NAs to zero ("impute")
            if ( anyNA(Xi) )
                Xi[ is.na(Xi) ] <- 0

            # use projection trick explained in covar_logit_BEDMatrix
            # use eigenvectors V to project data (via P_V, transformed V, see above)
            Pi <- ( Xi %*% P_V + Xi_mean ) / ploidy
            # cap allele frequencies (up to now, they could have been out of range)
            Pi <- af_cap( Pi )
            # save to output matrix
            P[ indexes_loci_chunk, ] <- Pi
        }

        # done with BEDMatrix version of processing!
        return( P )
    }
    # else below is regular X (R matrix, not BEDMatrix)

    # data dimensions
    m <- nrow(X)
    n <- ncol(X)
    
    adjust <- 8
    if( n-d < 10 )
        adjust <- n-d-1
    
    # index the missing values
    NA_IND <- is.na(X)
    
    # center the matrix...
    X_mean <- rowMeans(X, na.rm=TRUE)

    if ( d == 0 ) {
        # this is "intercept only"
        # all allele frequencies are just the mean
        P <- matrix(
            rep.int( X_mean, n ),
            nrow = m,
            ncol = n
        )
        return( P )
    }
    # else continue
    
    norm_X <- X - X_mean
    # sd_X <- apply(X, 1, function(snp){sd(snp, na.rm=TRUE)})
    
    #...then 'impute'
    norm_X[NA_IND] <- 0
    
    mysvd <- trunc_svd(norm_X, d = d, adjust = adjust, tol = tol, override = override)
    
    rm(norm_X)
    D <- mysvd$d
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    
    P <- U %*% diag(D, d, d)  %*% t(V)
        
    P <- P + X_mean
    P <- P / ploidy

    # cap allele frequencies (up to now, they could have been out of range)
    P <- af_cap(P)
    
    # WHY???
    #P[NA_IND] <- NA
    
    return( P )
}
