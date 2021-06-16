# for PCA method.  Caps IAFs to 1/(2n) of [0,1] boundary.  Preserves input type
# (matrix vs vector).
.af_cap <- function(P) {
    if (missing(P))
        stop("Individual-specific allele frequency matrix `P` is required!")

    # getting sample size is only step that varies between vectors and matrices
    n_ind <- if (is.matrix(P))
        ncol(P) else length(P)

    # calculate allele frequency caps according to sample size
    p_cap_lo <- 1/(2 * n_ind)
    p_cap_hi <- 1 - p_cap_lo  # symmetric capping

    # apply caps throughout the matrix or vector
    P[P > p_cap_hi] <- p_cap_hi
    P[P < p_cap_lo] <- p_cap_lo

    # return modified individual-specific allele frequency matrix
    return(P)
}

