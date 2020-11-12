#' @title Allele frequencies
#' @description Compute matrix of individual-specific allele frequencies
#' @inheritParams lfa
#' @param LF Matrix of logistic factors, with intercept. Pass in the 
#' return value from \code{lfa}!
#' @details Computes the matrix of individual-specific allele 
#' frequencies, which has the same dimensions of the genotype matrix.
#' Be warned that this function could use a ton of memory, as the 
#' return value is all doubles. It could be wise to pass only a 
#' selection of the SNPs in your genotype matrix to get an idea for
#' memory usage. Use \code{gc} to check memory usage!
#' @examples
#' LF = lfa(hgdp_subset, 4)
#' allele_freqs = af(hgdp_subset, LF)
#' @return Matrix of individual-specific allele frequencies.
#' @export
af <- function( X, LF, safety = FALSE ){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( is.null(LF) )
        stop( "`LF` matrix is required!" )
    
    # check class
    is_BEDMatrix <- FALSE
    if ( "BEDMatrix" %in% class(X) ) {
        is_BEDMatrix <- TRUE
    } else if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )

    if (! is_BEDMatrix ) {
        # usual R object behavior
        if (safety)
            check_geno(X)
        
        return(
            t( apply( X, 1, af_snp, LF ) )
        )
    } else {
        # BEDMatrix case
        # write an explicit loop around the genotype matrix
        # questions: is it better to write a simple loop (one locus at the time) or to read big chunks (1000 loci at the time)?  Since `af_snp` is the bottleneck, maybe the difference is small
        m <- ncol(X)
        n <- nrow(X)
        # output matrix
        P <- matrix( 0, m, n )
        for ( i in 1 : m ) {
            # get locus i genotype vector
            xi <- X[ , i ]
            # calculate and store result
            P[ i, ] <- af_snp( xi, LF )
        }
        # done!
        return( P )
    }
}

