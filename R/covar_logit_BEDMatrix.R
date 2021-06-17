# based on .covar_BEDMatrix / popkinsuppl::kinship_std.  Computes
# covariance for second SVD step of LFA.
.covar_logit_BEDMatrix <- function(X, X_mean, V, ploidy = 2, m_chunk = 1000) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(X_mean))
        stop("Mean locus frequency `X_mean` is required!")
    if (missing(V))
        stop("Truncated eigenvector matrix `V` is required!")
    if (!("BEDMatrix" %in% class(X)))
        stop("`X` must be a BEDMatrix object!")
    # get dimensions
    n <- nrow(X)
    m <- ncol(X)
    # turn eigenvector matrix into a projection matrix
    P_V <- tcrossprod(V)
    # initialize desired matrix and vector
    covar <- matrix(0, nrow = n, ncol = n)
    # navigate chunks
    i_chunk <- 1  # start of first chunk
    while (TRUE) {
        # start an infinite loop, break inside as needed
        if (i_chunk > m)
            break  # reached end
        # range of SNPs to extract in this chunk
        loci_chunk <- i_chunk:min(i_chunk + m_chunk - 1, m)
        # transpose for our usual setup (makes centering easiest)
        Xi <- t(X[, loci_chunk, drop = FALSE])
        # update for next chunk! (overshoots at end, that's ok)
        i_chunk <- i_chunk + m_chunk
        Xi_mean <- X_mean[loci_chunk]  # precomputed data
        Xi <- Xi - Xi_mean  # center
        if (anyNA(Xi))
            Xi[is.na(Xi)] <- 0 # set NAs to zero ('impute')
        # project data using first-pass eigenvecs (P_V)
        Zi <- (Xi %*% P_V + Xi_mean)/ploidy
        # apply LFA threshold to this subset, will remove some loci
        loci_keep <- as.logical(.Call("lfa_threshold", Zi, 1/(ploidy * n)))
        if (!any(loci_keep))
            next  # move on if nothing passed
        Zi <- Zi[loci_keep, , drop = FALSE]  # subset loci
        Zi <- log(Zi/(1 - Zi))  # logit transform whole matrix
        Zi <- centerscale(Zi)
        covar <- covar + crossprod(Zi) # add to running sum.
    }
    return(covar)
}

# projection trick.  Full rank version: (X is centered matrix though!): X = U D
# t(V); X V = U D t(V) V = U D; X V D^(-1) = U; limited rank now: Z = U_r D_r
# t(V_r); Z = (X V D^(-1))_r D_r t(V_r); Z = (X V)_r D_r^(-1) D_r t(V_r); Z = X
# V_r t(V_r)
