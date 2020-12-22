pvals_empir <- function( stats1, stats0 ) {
    if ( missing( stats1 ) )
        stop( '`stats1` observed statistics vector is required!' )
    if ( missing( stats0 ) )
        stop( '`stats0` null statistics (vector or matrix) is required!' )
    # NOTE: values in either stats1 and stats0 can be NA
    # the observed cases must be kept (their p-values are NA too)
    
    # helps us map back to original position after sort below
    # NOTE: preserves NA, orders them last!
    stats1_order <- order( stats1 )
    # NAs go last, same length as input
    stats1_sorted <- stats1[ stats1_order ]
    
    # stats0 can be a matrix, but we need it as a vector now
    # original order doesn't matter, so let's sort and forget about original
    # NAs are also removed
    stats0 <- sort( as.vector( stats0 ) )
    # number of non-NA values in stats0
    m0 <- length( stats0 )

    # begin calculating p-values!
    m <- length( stats1 )
    # initialize to NAs, to preserve stats1 NAs
    pvals <- rep(NA, m)
    i0 <- 1 # index on stats0
    for ( i1 in 1 : m ) {
        # i1 is index in stats1_sorted
        # look at i1'th observed statistic
        stats1_sorted_i1 <- stats1_sorted[ i1 ]
        
        # if there are any NAs, the appear in the end
        # in that case, time to stop!
        if ( is.na( stats1_sorted_i1 ) )
            break
        
        # i0 is the number of null stats that are *smaller or equal* than the current observed statistic
        # so below the p-value is the proportion of null stats that is *strictly larger* than the current observed statistic
        
        # i0 gets incremented until the null statistic i0 exceeds the observed statistic i1 (both in ascending order)
        while ( ( i0 <= m0 ) && ( stats1_sorted_i1 >= stats0[ i0 ] ) ) {
            i0 <- i0 + 1
        }
        # the p-value is 1 minus the proportion of null statistics smaller than the observed statistic
        # stats1_order[ i1 ] puts the value in the original position (verified in a toy case)
        pvals[ stats1_order[ i1 ] ] <- 1 - ( ( i0 - 1 ) / m0 )
        # equals
        # 1 - ( ( i0 - 1 ) / m0 )
        # = ( m0 - i0 + 1 ) / m0
        # extrema:
        # i0 == 1, p = 1
        # i0 == m0, p = 1 / m0
    }
    
    return( pvals )
}
