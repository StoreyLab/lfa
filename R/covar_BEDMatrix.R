# based on popkinsuppl::kinship_std
# - limited to BEDMatrix case
# - returns locus mean values too (needed by LFA)
# - does not use popkin for memory control
# - does not normalize by p(1-p), to match how LFA does it
# - similarly does not "average" non-NA cases, just sums (doesn't average for fixed m either)
covar_BEDMatrix <- function(
                            X,
                            m_chunk = 1000 # gave good performance in tests
                            ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( !('BEDMatrix' %in% class(X) ))
        stop( '`X` must be a BEDMatrix object!' )

    # get dimensions
    n <- nrow(X)
    m <- ncol(X)
    
    # initialize desired matrix and vector
    covar <- matrix( 0, nrow = n, ncol = n )
    X_mean <- vector( 'numeric', m )
    
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

        # standard mean
        Xi_mean <- rowMeans(Xi, na.rm = TRUE)
        # store for output too
        X_mean[ indexes_loci_chunk ] <- Xi_mean
        
        # center before cross product...
        Xi <- Xi - Xi_mean
        
        # set NAs to zero ("impute")
        if ( anyNA(Xi) )
            Xi[ is.na(Xi) ] <- 0
        
        # cross product matrix at this SNP, add to running sum.
        covar <- covar + crossprod(Xi)
    }
    
    return(
        list(
            covar = covar,
            X_mean = X_mean
        )
    )
}

