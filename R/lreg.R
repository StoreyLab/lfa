# C based logistic regression
.lreg <- function(x, LF, max_iter = 100, tol = 1e-10) {
    if (missing(x))
        stop("Genotype vector `x` is required!")
    if (is.null(LF))
        stop("`LF` matrix is required!")

    # make sure the data is NA-free.  Focus on x only, that's a more common
    # issue (it'd be wasteful to test LFs repeatedly (for each locus))
    if (anyNA(x))
        stop("Genotype vector `x` must not have NA values!")

    # why weird doubling of everything?
    LF2 <- rbind(LF, LF)
    x1 <- as.numeric((x == 1) | (x == 2))
    x2 <- as.numeric(x == 2)
    x2 <- c(x1, x2)
    # get the desired coefficients
    betas <- .Call("lreg_c", LF2, x2, max_iter, tol)

    # if coefficients are NA, use glm
    if (anyNA(betas)) {
        # `-1` is because LF already has intercept.  NOTE: this reduces betas
        # by 1 as well, we don't match `lreg_c` otherwise!
        # suppressWarnings: because sometimes we get warning 'glm.fit: fitted
        # probabilities numerically 0 or 1 occurred'. Occurs on randomly
        # simulated data, nothing particularly problematic, so meh
        suppressWarnings(betas <- stats::glm(cbind(x, 2 - x) ~ -1 + LF,
            family = "binomial")$coef)
        names(betas) <- NULL
    }
    return(betas)
}
