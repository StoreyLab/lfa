.pca_af_BEDMatrix <- function(X, d, ploidy, m_chunk) {
    # data dimensions (transposed for BEDMatrix)
    m <- ncol(X)
    n <- nrow(X)
    # Calculate covariance matrix and loci means.
    # NOTE: inefficient for d=0 (no PCs, just mean).
    obj <- .covar_BEDMatrix(X, m_chunk=m_chunk)
    covar <- obj$covar
    X_mean <- obj$X_mean
    # this is 'intercept only', all allele frequencies are just the mean
    if (d == 0)
        return(matrix(rep.int(X_mean, n), nrow=m, ncol=n))
    # get truncated eigendecomposition
    obj <- RSpectra::eigs_sym(covar, d)
    V <- obj$vectors
    P_V <- tcrossprod(V) # turn eigenvectors into projection matrix
    # Form P in parts, so X is not in memory all at once.
    # P is fully in memory, potentially negating the BEDMatrix advantage
    P <- matrix(0, nrow = m, ncol = n) # initialize
    # navigate chunks
    i_chunk <- 1 # start of first chunk
    while (TRUE) {
        # start an infinite loop, break inside as needed.
        if (i_chunk > m)
            break # reached end
        # range of SNPs to extract in this chunk
        indexes_loci_chunk <- i_chunk:min(i_chunk + m_chunk - 1, m)
        # transpose for our usual setup (makes centering easiest)
        Xi <- t(X[, indexes_loci_chunk, drop = FALSE])
        # update for next chunk! (overshoots at end, that's ok)
        i_chunk <- i_chunk + m_chunk
        # get row means from precomputed data
        Xi_mean <- X_mean[indexes_loci_chunk]
        Xi <- Xi - Xi_mean # center
        if (anyNA(Xi))
            Xi[is.na(Xi)] <- 0 # set NAs to zero ('impute')
        # project data to V rowspace
        Pi <- (Xi %*% P_V + Xi_mean)/ploidy
        # cap allele frequencies, store in output matrix
        P[indexes_loci_chunk, ] <- .af_cap(Pi)
    }
    return(P)
}
