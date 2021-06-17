# internal version of lfa for BEDMatrix data.  `d` should be `d-1` from `lfa`
# input!  `adjustments` not supported yet.  `lfa` checks not repeated
.lfa_BEDMatrix <- function(X, d, ploidy = 2, m_chunk = 1000) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(d))
        stop("Dimension number `d` is missing!")
    if (!("BEDMatrix" %in% class(X)))
        stop("`X` must be a BEDMatrix object!")

    # calculate covariance matrix and loci means
    obj <- .covar_BEDMatrix(X, m_chunk = m_chunk)
    covar <- obj$covar
    X_mean <- obj$X_mean

    # get truncated eigendecomposition
    obj <- RSpectra::eigs_sym(covar, d)
    V <- obj$vectors

    # calculate covariance matrix for second step (after logit
    # filter/transform)
    covar <- .covar_logit_BEDMatrix(X, X_mean, V)

    # get truncated eigendecomposition of second level, which yields the
    # logistic factors
    obj <- RSpectra::eigs_sym(covar, d)
    V <- obj$vectors

    # add intercept column last
    V <- cbind(V, 1)

    return(V)
}


