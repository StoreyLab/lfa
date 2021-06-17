.convtests <- function(Bsz, tol, d_org, residuals, d, Smax) {
    Len_res <- sum(residuals[seq_len(d_org)] < tol * Smax)
    # if this happens, we've converged!
    if (Len_res == d_org)
        return(list(converged=TRUE, d=d))

    # ... otherwise not converged.
    d <- max(d, d_org + Len_res)
    if (d > Bsz - 3)
        d <- Bsz - 3
    return(list(converged=FALSE, d=d))
}
