#' @title Allele frequencies for SNP
#' @description Computes individual-specific allele frequencies for a 
#' single SNP.
#' @inheritParams af
#' @param snp vector of 0's, 1's, and 2's
#' @return vector of allele frequencies
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' # pick one SNP only
#' snp <- hgdp_subset[ 1, ]
#' # allele frequency vector for that SNO only
#' allele_freqs_snp <- af_snp(snp, LF)
#' @seealso [af()]
#' @export
af_snp <- function(snp, LF, max_iter = 100, tol = 1e-10) {
    if (missing(snp))
        stop("`snp` is required!")
    if (missing(LF))
        stop("`LF` matrix is required!")

    # dimensions should agree
    if (length(snp) != nrow(LF))
        stop("Number of individuals in `snp` and `LF` disagree!")

    # can only regress with non-NA individuals
    indexes_keep <- !is.na(snp)
    snp <- snp[indexes_keep]  # overwrite
    LF2 <- LF[indexes_keep, , drop = FALSE]  # don't overwite LF
    # get coefficients from logistic regression
    betas <- .lreg(snp, LF2, max_iter, tol)

    # get allele frequencies.  Though `snp` may be NA, no `LF` or `beta` are
    # NA, so this imputes the missing genotypes!
    est <- .Call("mv_c", LF, betas)
    # very large `est` can result in NaN's (i.e. est==1000).  Oddly, very large
    # negative are no problem
    af <- ifelse(est > 100, 1, exp(est)/(1 + exp(est)))
    return(af)
}
