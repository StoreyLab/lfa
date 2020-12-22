# compute empirical p-values by brute force (clear implementation)
# used only to validate `pvals_empir`, which is way faster on large data, but code is much harder to understand
# both versions handle NAs in inputs
pvals_empir_brute <- function( stats1, stats0 ) {

    # first remove NAs in stats0
    if ( anyNA( stats0 ) )
        stats0 <- stats0[ !is.na( stats0 ) ]
    # NAs in stats1 get preserved though
    
    # for loop and output length
    m <- length( stats1 )
    # for p-value normalization
    m0 <- length( stats0 )

    # create output vector
    pvals <- rep.int( NA, m )
    for ( i in 1 : m ) {
        # NAs are preserved
        if ( !is.na( stats1[ i ] ) )
            pvals[i] <- sum( stats0 > stats1[ i ] ) / m0
    }
    return( pvals )
}

