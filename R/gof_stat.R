.gof_stat <- function(X, LF, max_iter = 100, tol = 1e-10) {
    # wrapper around .gof_stat_snp, applying it correctly across matrix whether
    # input is a regular R matrix or a BEDMatrix object

    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(LF))
        stop("`LF` matrix is required!")

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
        stop("Number of individuals in `X` and `LF` disagree!")
    
    if (!methods::is(X, "BEDMatrix")) {
        # usual R object behavior
        gof_stats <- apply(X, 1, .gof_stat_snp, LF, max_iter=max_iter, tol=tol)
    } else {
        # BEDMatrix case: write an explicit loop around the genotype matrix.
        # Questions: is it better to write a simple loop (one locus at the
        # time) or to read big chunks (1000 loci at the time)?  Since `af_snp`
        # is the bottleneck, maybe the difference is small

        # output vector
        gof_stats <- vector("numeric", m)
        for (i in seq_len(m)) {
            # get locus i genotype vector
            xi <- X[, i]
            # calculate and store result
            gof_stats[i] <- .gof_stat_snp(xi, LF, max_iter=max_iter, tol=tol)
        }
    }
    return(gof_stats)
}
