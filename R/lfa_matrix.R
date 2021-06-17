# LFA for in-memory R matrices only (as opposed to BEDMatrix).
# Skips validations already performed in [lfa()].
.lfa_matrix <- function(X, d, adjustments, override, rspectra, ploidy, tol) {
    n <- ncol(X) # number of individuals
    if (!rspectra) {
        adjust <- 8 # a mysterious param for trunc_svd
        if (n - d < 10)
            adjust <- n - d - 1
    }
    NA_IND <- is.na(X) # index the missing values
    mean_X <- rowMeans(X, na.rm=TRUE)
    norm_X <- X - mean_X # center
    norm_X[NA_IND] <- 0 # then 'impute'
    # first SVD
    if (rspectra) {
        mysvd <- RSpectra::svds(norm_X, d)
    } else {
        mysvd <- trunc_svd(norm_X, d, adjust, tol, override=override)
    }
    rm(norm_X)
    D <- diag(mysvd$d, nrow=d) # pass `d` so `diag` gets `d=1` right
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    # form projection
    z <- U %*% D %*% t(V)
    z <- z + mean_X
    z <- z/ploidy
    rm(U); rm(D); rm(V)
    # remove rows that exceed logit (0,1) domain
    z <- z[as.logical(.Call("lfa_threshold", z, 1/(ploidy * n))), ]
    z <- log(z/(1 - z)) # logit
    norm_z <- centerscale(z) # center/scale in logit scale now
    # regress out adjustment vars, if relevant
    if (!is.null(adjustments)) {
        norm_z <- t(stats::residuals(stats::lm(t(norm_z) ~ adjustments - 1)))
        d <- d - ncol(adjustments)
    }
    # second SVD yields the logistic factors
    if (rspectra) {
        v <- RSpectra::svds(norm_z, d)$v
    } else {
        v <- trunc_svd(norm_z, d, adjust, tol, override=override)$v
    }
    v <- cbind(v, 1) # add intercept column last
    if (!is.null(adjustments))
        v <- cbind(adjustments, v) # add adjustment variables first
    return(v)
}
