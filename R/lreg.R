# C based logistic regression 
lreg <- function(x, LF, max_iter = 20, tol = 1e-10){
    if ( missing(x) )
        stop( 'Genotype vector `x` is required!' )
    if ( is.null(LF) )
        stop( "`LF` matrix is required!" )

    # question: weird doubling of everything, appears to be for modeling dominance effects
    LF <- rbind(LF, LF)
    x1 <- as.numeric( ( x == 1 ) | ( x == 2 ) )
    x2 <- as.numeric( x == 2 )
    x <- c(x1, x2)
    # get the desired coefficients
    betas <- .Call("lreg_c", LF, x, max_iter, tol)
    
    # if coefficients are NA, use glm
    if ( anyNA(betas) ) {
        # print(paste("harmless warning:"))
        # `-1` is because LF already has intercept
        # questions:
        # - is there a second doubling of data here?
        betas <- stats::glm(
            cbind(x, 2-x) ~ -1 + LF,
            family = "binomial"
        )$coef
        names(betas) <- NULL
    }
    return( betas )
}
