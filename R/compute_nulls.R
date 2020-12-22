compute_nulls <- function(P, d, B) {
    m <- nrow( P )
    n <- ncol( P )
    
    # since m and n are integers, multiplying them causes a buffer overflow
    # let's multiply them as doubles, overcomes the problem
    n_data <- ( n + 0.0 ) * ( m + 0.0 )
    
    stats0 <- matrix( 0, m, B )
    for(i in 1:B) {
        X0 <- matrix(
            stats::rbinom( n_data, 2, P ),
            nrow = m,
            ncol = n
        )
        LF0 <- lfa( X0, d )
        # this calculates stats correctly, even when X0 is BEDMatrix!
        stats0[ , i ] <- gof_stat( X0, LF0 )
    }

    return( stats0 )
}
