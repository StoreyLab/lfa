convtests <- function(Bsz, tol, d_org, U_B, S_B, V_B, residuals, 
                      d, SVTol, Smax) {
    Len_res <- sum(residuals[1:d_org] < tol * Smax)
    
    if (Len_res == d_org)
        return(
            list(
                converged = TRUE,
                U_B = U_B[ , 1:d_org, drop = FALSE ],
                S_B = S_B[ 1:d_org, drop = FALSE ], 
                V_B = V_B[ , 1:d_org, drop = FALSE ],
                d = d
            )
        )
    
    Len_res <- sum(residuals[1:d_org] < SVTol * Smax)
    d <- max(d, d_org+Len_res)
    if (d > Bsz - 3) 
        d <- Bsz - 3
    return(
        list(
            converged = FALSE,
            d = d
        )
    )
}
