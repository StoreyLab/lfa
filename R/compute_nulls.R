compute_nulls <- function(P, d, B) {
    m <- nrow(P)
    n <- ncol(P)

    stat0 <- matrix(0, m, B)
    for(i in 1:B) {
        X0 <- matrix(
            stats::rbinom(m*n, 2, P),
            m,
            n
        )
        LF <- lfa(X0, d)
        stat0[ , i ] <- apply(X0, 1, gof_stat_snp, LF)
    }

    return(stat0)
}
