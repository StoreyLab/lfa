gof_stat <- function(X, LF) {
    # wrapper around gof_stat_snp, applying it correcto across matrix whether input is a regular R matrix or a BEDMatrix object

    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( LF ) )
        stop( '`LF` matrix is required!' )

    # check class
    is_BEDMatrix <- FALSE
    if ( "BEDMatrix" %in% class(X) ) {
        is_BEDMatrix <- TRUE
    } else if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )

    # get dimensions
    if ( is_BEDMatrix ) {
        m <- ncol(X)
        n <- nrow(X)
    } else {
        n <- ncol(X)
        m <- nrow(X)
    }

    # dimensions should agree
    if ( n != nrow(LF) )
        stop( 'Number of individuals in `X` must equal number of individuals (rows) in `LF`' )
    
    if (! is_BEDMatrix ) {
        # usual R object behavior
        gof_stats <- apply(X, 1, gof_stat_snp, LF)
    } else {
        # BEDMatrix case
        # write an explicit loop around the genotype matrix
        # questions: is it better to write a simple loop (one locus at the time) or to read big chunks (1000 loci at the time)?  Since `af_snp` is the bottleneck, maybe the difference is small
        
        # output vector
        gof_stats <- vector( 'numeric', m )
        for ( i in 1 : m ) {
            # get locus i genotype vector
            xi <- X[ , i ]
            # calculate and store result
            gof_stats[ i ] <- gof_stat_snp( xi, LF )
        }
    }
    return( gof_stats )
}
