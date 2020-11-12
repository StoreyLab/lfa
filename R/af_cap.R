# for PCA method
# input is individual-specific allele frequency matrix estimate, with values potentially out of range
# this function caps data to be within a reasonable range (not setting any frequencies to 0 or 1, but rather to the minimum according to the sample size
# if input is matrix, returns matrix
# if input is vector, returns vector
af_cap <- function ( P ) {
    if ( missing( P ) )
        stop( 'Individual-specific allele frequency matrix `P` is required!' )
    
    # getting sample size is only step that varies between vectors and matrices
    n_ind <- if ( is.matrix( P ) ) ncol( P ) else length( P )
    
    # calculate allele frequency caps according to sample size
    p_cap_lo <- 1 / ( 2 * n_ind )
    p_cap_hi <- 1 - p_cap_lo # symmetric capping

    # apply caps throughout the matrix or vector
    P[ P > p_cap_hi ] <- p_cap_hi
    P[ P < p_cap_lo ] <- p_cap_lo
    
    # return modified individual-specific allele frequency matrix
    return( P )
}

