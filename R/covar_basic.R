# basic covariance formula for an R genotype matrix X, to checks the more
# elaborate .covar_BEDMatrix against this
.covar_basic <- function(X) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (!is.matrix(X))
        stop("`X` must be a matrix!")

    # standard mean
    X_mean <- rowMeans(X, na.rm=TRUE)

    # center before cross product...
    X <- X - X_mean

    # before applying cross product, to prevent NA errors, just set those
    # values to zero and it works out!
    if (anyNA(X))
        X[is.na(X)] <- 0

    # cross product matrix is what we desire
    covar <- crossprod(X)

    return(covar)
}

