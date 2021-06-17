.pvals_empir <- function(stats1, stats0) {
    if (missing(stats1))
        stop("`stats1` observed statistics vector is required!")
    if (missing(stats0))
        stop("`stats0` null statistics (vector or matrix) is required!")
    # NOTE: values in `stats1`, `stats0` can be NA. Observed cases must be kept
    # (their p-values are NA too)

    # helps us map back to original position after sort below. NOTE: preserves
    # NA, orders them last!
    stats1_order <- order(stats1)
    # NAs go last, same length as input
    stats1_sorted <- stats1[stats1_order]

    # Flatten stats0 to vector.  Original order doesn't matter (sort and forget
    # about original).  NAs are removed
    stats0 <- sort(as.vector(stats0))
    # number of non-NA values in stats0
    m0 <- length(stats0)

    # begin calculating p-values!
    m <- length(stats1)
    # initialize to NAs, to preserve stats1 NAs
    pvals <- rep(NA, m)
    i0 <- 1  # index on stats0
    for (i1 in seq_len(m)) {
        # i1 is index in stats1_sorted; look at i1'th observed statistic
        stats1_sorted_i1 <- stats1_sorted[i1]

        # since NAs appear in the end, stop if we see one
        if (is.na(stats1_sorted_i1))
            break

        # i0 = |null stats <= current observed stat|, so p-value is proportion
        # of null stats *strictly larger* than the current observed stat.
        # Increment i0 until null stat i0 exceeds obs stat i1 (both ascending)
        while ((i0 <= m0) && (stats1_sorted_i1 >= stats0[i0])) {
            i0 <- i0 + 1
        }
        # pval = 1 - prop null stats smaller than obs stat.  stats1_order[ i1 ]
        # puts value in orig pos
        pvals[stats1_order[i1]] <- 1 - ((i0 - 1)/m0)
    }

    return(pvals)
}
