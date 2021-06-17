#' @title Hardy-Weinberg Equilibrium in structure populations
#' 
#' @description
#' Compute structural Hardy-Weinberg Equilibrium (sHWE) p-values
#' on a SNP-by-SNP basis. These p-values can be aggregated to 
#' determine genome-wide goodness-of-fit for a particular value
#' of `d`. See <https://doi.org/10.1101/240804> for more
#' details.
#'
#' @param LF matrix of logistic factors
#' @param B number of null datasets to generate, `B = 1` is usually
#' sufficient. If computational time/power allows, a few extra
#' `B` could be helpful
#' @inheritParams lfa
#' @inheritParams af
#' @examples
#' # get LFs
#' LF <- lfa(hgdp_subset, 4)
#' # look at a small (300) number of SNPs for rest of this example:
#' hgdp_subset_small <- hgdp_subset[ 1:300, ]
#' gof_4 <- sHWE(hgdp_subset_small, LF, 3)
#' LF <- lfa(hgdp_subset, 10)
#' gof_10 <- sHWE(hgdp_subset_small, LF, 3)
#' hist(gof_4)
#' hist(gof_10)
#' @return a vector of p-values for each SNP.
#' @export
sHWE <- function(X, LF, B, max_iter = 100, tol = 1e-10) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(LF))
        stop("`LF` matrix is required!")
    if (missing(B))
        stop("`B` scalar is required!")

    # check class
    if (!is.matrix(X)) # BEDMatrix returns TRUE
        stop("`X` must be a matrix!")

    # get dimensions
    if (methods::is(X, "BEDMatrix")) {
        m <- ncol(X)
        n <- nrow(X)
    } else {
        n <- ncol(X)
        m <- nrow(X)
    }

    # dimensions should agree
    if (n != nrow(LF))
        stop("Number of individuals in `X` and `LF` disagrees!")

    # calculate observed stats across matrix
    stats1 <- .gof_stat(X, LF, max_iter=max_iter, tol=tol)

    # to create null statistics, get P matrix, then simulate data from it
    d <- ncol(LF)
    # this already works on BEDMatrix, but produces this large matrix!
    P <- af(X, LF)
    rm(X)
    stats0 <- .compute_nulls(P, d, B, max_iter=max_iter, tol=tol)

    # calculate empirical p-values based on these distributions
    pvals <- .pvals_empir(stats1, stats0)

    return(pvals)
}
