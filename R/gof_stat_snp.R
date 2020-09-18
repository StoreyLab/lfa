gof_stat_snp <- function(snp, LF){
    # remove NAs before calculating GOF statistics
    NA_IND <- is.na(snp)
    snp <- snp[!is.na(snp)]
    LF  <- LF[!NA_IND, ,drop=FALSE]

    # get vector of allele frequencies at this SNP
    p <- af_snp(snp, LF)

    p0 <- (1-p)^2
    p1 <- 2*p*(1-p)

    est <- c(sum(p0), sum(p1))
    N   <- c(sum( snp==0 ), sum( snp==1 ))

    SIGMA <- matrix(
        c(sum(p0*(1-p0)), -1*sum(p0*p1),
          -1*sum(p0*p1), sum(p1*(1-p1))),
        nrow = 2,
        ncol = 2
    )
    # invert matrix first (it may be NA)
    SIGMA_inv <- inverse_2x2(SIGMA)
    # stat is NA if that matrix wasn't invertible
    stat <- if ( is.null( SIGMA_inv ) ) {
                NA
            } else {
                t(N-est) %*% SIGMA_inv %*% (N-est)
            }
    return(stat)
}

