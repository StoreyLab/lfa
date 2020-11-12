# based on covar_BEDMatrix, itself based on popkinsuppl::kinship_std
# - computes effective covariance matrix that is used in the second SVD step of LFA, which involves projection into the first SVD step, truncation, and logit transformation
covar_logit_BEDMatrix <- function(
                                  X,
                                  X_mean,
                                  V,
                                  ploidy = 2,
                                  m_chunk = 1000 # gave good performance in tests
                                  ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( X_mean ) )
        stop( 'Mean locus frequency `X_mean` is required!' )
    if ( missing( V ) )
        stop( 'Truncated eigenvector matrix `V` is required!' )
    if ( !('BEDMatrix' %in% class(X) ))
        stop( '`X` must be a BEDMatrix object!' )

    # get dimensions
    n <- nrow(X)
    m <- ncol(X)

    # turn eigenvector matrix into a projection matrix
    P_V <- tcrossprod( V )
    
    # initialize desired matrix and vector
    covar <- matrix( 0, nrow = n, ncol = n )
    
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

        # projection trick
        #
        # full rank version: (X is centered matrix though!)
        # X = U D t(V)
        # X V = U D t(V) V = U D
        # X V D^(-1) = U
        #
        # limited rank now:
        # Z = U_r D_r t(V_r)
        # Z = (X V D^(-1))_r D_r t(V_r)
        # Z = (X V)_r D_r^(-1) D_r t(V_r)
        # Z = X V_r t(V_r)
        
        # use first pass eigenvectors V to project data (via P_V, transformed V, see above)
        Zi <- Xi %*% P_V + Xi_mean
        Zi <- Zi / ploidy

        # apply LFA threshold to this subset, will remove some loci
        indexes_loci_keep <- as.logical( .Call( "lfa_threshold", Zi, 1 / ( ploidy * n ) ) )
        # for small enough chunks, nobody may have passed, just move on to next chunk!
        if ( !any( indexes_loci_keep ) )
            next
        # subset loci
        Zi <- Zi[ indexes_loci_keep, , drop = FALSE ]
        # logit transformation of whole matrix
        Zi <- log( Zi / ( 1 - Zi ) )

        # center and scale this reduced matrix
        Zi <- centerscale( Zi )
        
        # cross product matrix at this SNP, add to running sum.
        covar <- covar + crossprod( Zi )
    }
    
    return( covar )
}

