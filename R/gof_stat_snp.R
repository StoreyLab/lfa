.gof_stat_snp <- function(snp, LF, max_iter = 100, tol = 1e-10) {
    # remove NAs before calculating GOF statistics
    keep <- !is.na(snp)
    snp <- snp[keep]
    LF <- LF[keep, , drop = FALSE]
    # get vector of allele frequencies at this SNP
    p <- af_snp(snp, LF, max_iter = max_iter, tol = tol)
    # some intermediate calcs
    p0 <- (1 - p)^2
    p1 <- 2 * p * (1 - p)
    est <- c(sum(p0), sum(p1))
    N <- c(sum(snp == 0), sum(snp == 1))
    # construct Sigma and inverse
    sigma11 <- sum(p0 * (1 - p0))
    sigma12 <- -sum(p0 * p1)
    sigma22 <- sum(p1 * (1 - p1))
    # determinant
    determ <- sigma11 * sigma22 - sigma12^2
    # not invertible if this is zero
    if (determ == 0)
        return(NA)
    # else continue
    Sigma_inv <- c(sigma22, -sigma12, -sigma12, sigma11)
    Sigma_inv <- matrix(Sigma_inv, nrow=2, ncol=2)/determ
    stat <- t(N - est) %*% Sigma_inv %*% (N - est)
    return(stat)
}
