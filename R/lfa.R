#' Logistic factor analysis
#'
#' Fit logistic factor model of dimension `d` to binomial data.
#' Computes `d - 1` singular vectors followed by intercept.
#'
#' Genotype matrix should have values in 0, 1, 2, or `NA`.
#' The coding of the SNPs (which case is 0 vs 2) does not change the output.
#'
#' @param X A matrix of SNP genotypes, i.e. an integer matrix of 0's,
#' 1's, 2's and `NA`s.
#' BEDMatrix is supported.
#' Sparse matrices of class Matrix are not supported (yet).
#' @param d Number of logistic factors, including the intercept
#' @param adjustments A matrix of adjustment variables to hold fixed during 
#' estimation.  Number of rows must equal number of individuals in `X`.
#' These adjustments take the place of LFs in the output, so the number of
#' columns must not exceed `d-2` to allow for the intercept and at least one
#' proper LF to be included.
#' When present, these adjustment variables appear in the first columns of the
#' output.
#' Not supported when `X` is a BEDMatrix object.
#' @param rspectra If `TRUE`, use
#' [RSpectra::svds()] instead of default
#' [trunc_svd()] or
#' [corpcor::fast.svd()] options.
#' Ignored if `X` is a BEDMatrix object.
#' @param override Optional boolean passed to [trunc_svd()] 
#' to bypass its Lanczos bidiagonalization SVD, instead using
#' [corpcor::fast.svd()].
#' Usually not advised unless encountering a bug in the SVD code.
#' Ignored if `X` is a BEDMatrix object.
#' @param safety Optional boolean to bypass checks on the genotype
#' matrices, which require a non-trivial amount of computation.
#' Ignored if `X` is a BEDMatrix object.
#' @param ploidy Ploidy of data, defaults to 2 for bi-allelic unphased SNPs
#' @param tol Tolerance value passed to [trunc_svd()]
#' Ignored if `X` is a BEDMatrix object.
#' @param m_chunk If `X` is a BEDMatrix object, number of loci to read per
#' chunk (to control memory usage).
#' 
#' @return The matrix of logistic factors, with individuals along rows and
#' factors along columns.
#' The intercept appears at the end of the columns, and adjustments in the
#' beginning if present.
#' 
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' dim(LF)
#' head(LF)
#' @useDynLib lfa, .registration = TRUE
#' @export
lfa <- function(X, d, adjustments = NULL, override = FALSE, safety = FALSE,
    rspectra = FALSE, ploidy = 2, tol = .Machine$double.eps, m_chunk = 1000) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(d))
        stop("Dimension number `d` is required!")
    # check class
    if (!is.matrix(X)) # BEDMatrix returns TRUE
        stop("`X` must be a matrix!")
    # data dimensions (BEDMatrix is transposed)
    n <- if (methods::is(X, "BEDMatrix")) nrow(X) else ncol(X)
    # check for d validity
    if (!is.numeric(d)) {
        stop("d must be numeric")
    } else if (d != as.integer(d)) {
        stop("d should be integer")
    } else if (d < 1) {
        stop("d should be at least 1")
    } else if (d == 1) {
        return(matrix(1, n, 1)) # return intercept column vector only
    } else if (d > 1) {
        d <- d - 1  #for the svd stuff
    }
    # check adjustments vars
    if (!is.null(adjustments)) {
        if (methods::is(X, "BEDMatrix"))
            stop("`adjustments` are not supported when `X` is class BEDMatrix!")
        if (!is.matrix(adjustments))
            stop("`adjustments` must be a matrix!")
        if (nrow(adjustments) != n)
            stop("`adjustments` and `X` number of individuals disagree!")
        if (ncol(adjustments) >= d)
            stop("need to estimate at least one non-adjustment logistic factor")
        if (anyNA(adjustments))
            stop("`adjustments` must not have missing values!")
    }
    if (methods::is(X, "BEDMatrix"))
        return(.lfa_BEDMatrix(X, d, ploidy=ploidy, m_chunk=m_chunk))
    # else continue
    if (safety)
        .check_geno(X) # check data if asked to
    # now use 'R matrix' version of code, return those LFs
    return(.lfa_matrix(X, d, adjustments, override, rspectra, ploidy, tol))
}
