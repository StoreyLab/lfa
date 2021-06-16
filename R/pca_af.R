#' @title PCA Allele frequencies
#' @description Compute matrix of individual-specific allele frequencies
#' via PCA
#' @inheritParams lfa
#' @details This corresponds to algorithm 1 in the paper. Only used for
#' comparison purposes.
#' @return Matrix of individual-specific allele frequencies.
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' allele_freqs_lfa <- af(hgdp_subset, LF)
#' allele_freqs_pca <- pca_af(hgdp_subset, 4)
#' summary(abs(allele_freqs_lfa-allele_freqs_pca))
#' @export
pca_af <- function(X, d, override = FALSE, ploidy = 2, tol = 1e-13,
        m_chunk = 1000) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(d))
        stop("Principal components number `d` is required!")
    # check class
    if (!is.matrix(X)) # returns true for BEDMatrix
        stop("`X` must be a matrix!")
    # check for d validity
    if (d != as.integer(d)) {
        stop("d should be integer")
    } else if (d < 1) {
        stop("d should be at least 1")
    } else if (d >= 1) {
        d <- d - 1  #for the svd stuff
    }
    if (methods::is(X, "BEDMatrix"))
        return(.pca_af_BEDMatrix(X, d, ploidy, m_chunk))
    # else below is regular X (R matrix, not BEDMatrix)
    m <- nrow(X) # data dimensions
    n <- ncol(X)
    adjust <- 8
    if (n - d < 10)
        adjust <- n - d - 1
    X_mean <- rowMeans(X, na.rm=TRUE)
    if (d == 0) {
        # this is 'intercept only' all allele frequencies are just the mean
        P <- matrix(rep.int(X_mean, n), nrow=m, ncol=n)
        return(P)
    }
    norm_X <- X - X_mean # center
    norm_X[is.na(X)] <- 0 # then 'impute'
    mysvd <- trunc_svd(norm_X, d=d, adjust=adjust, tol=tol, override=override)
    rm(norm_X)
    D <- mysvd$d
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    P <- U %*% diag(D, d, d) %*% t(V)
    P <- P + X_mean
    P <- P/ploidy
    P <- .af_cap(P) # cap allele frequencies (they could be out of range)
    return(P)
}
