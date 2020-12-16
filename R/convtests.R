convtests <- function(
                      Bsz,
                      tol,
                      d_org,
                      residuals, 
                      d,
                      SVTol,
                      Smax
                      ) {
    # criterion based on `tol`
    Len_res <- sum( residuals[ 1 : d_org ] < tol * Smax )
    # if this happens, we've converged!
    if (Len_res == d_org)
        return(
            list(
                converged = TRUE,
                d = d
            )
        )

    # ... otherwise not converged
    # alter `d` now, for some reason
    # first redo criterion, but based on `SVTol` (why?)
    Len_res <- sum( residuals[ 1 : d_org ] < SVTol * Smax )
    d <- max( d, d_org + Len_res )
    if ( d > Bsz - 3 )
        d <- Bsz - 3
    return(
        list(
            converged = FALSE,
            d = d
        )
    )
}
