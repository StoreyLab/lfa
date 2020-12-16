compute_nulls <- function(P, d, B) {
    m <- nrow(P)
    n <- ncol(P)

    # since m and n are integers, multiplying them causes a buffer overflow
    # let's multiply them as doubles, overcomes the problem
    n_data <- ( n + 0.0 ) * ( m + 0.0 )
    
    stat0 <- matrix(0, m, B)
    for(i in 1:B) {
        X0 <- matrix(
            stats::rbinom( n_data, 2, P),
            m,
            n
        )
        LF <- lfa(X0, d)
        stat0[ , i ] <- apply(X0, 1, gof_stat_snp, LF)
    }

    return(stat0)
}
