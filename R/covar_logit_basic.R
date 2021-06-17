# basic covariance for second (logit) SVD of X, to check
# .covar_logit_BEDMatrix against
.covar_logit_basic <- function(X, V, ploidy = 2) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(V))
        stop("Truncated eigenvector matrix `V` is required!")
    if (!is.matrix(X))
        stop("`X` must be a matrix!")

    # get dimensions
    n <- ncol(X)

    # standard mean
    X_mean <- rowMeans(X, na.rm=TRUE)

    # center data
    X <- X - X_mean

    # set NAs to zero ('impute')
    if (anyNA(X))
        X[is.na(X)] <- 0

    # project data to rowspace of V (first-pass eigenvectors)
    Z <- X %*% tcrossprod(V) + X_mean
    Z <- Z/ploidy
    
    # apply LFA thhreshold to this subset
    ind <- as.logical(.Call("lfa_threshold", Z, 1/(ploidy * n)))
    # subset loci
    Z <- Z[ind, , drop=FALSE]
    # logit transformation of whole matrix
    Z <- log(Z/(1 - Z))

    # center and scale this reduced matrix
    Z <- centerscale(Z)

    # cross product matrix
    covar <- crossprod(Z)

    return(covar)
}

