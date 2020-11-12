#' @title Hardy-Weinberg Equilibrium in structure populations
#' 
#' @description
#' Compute structural Hardy-Weinberg Equilibrium (sHWE) p-values
#' on a SNP-by-SNP basis. These p-values can be aggregated to 
#' determine genome-wide goodness-of-fit for a particular value
#' of \eqn{d}. See \url{https://doi.org/10.1101/240804} for more
#' details.
#'
#' @param LF matrix of logistic factors
#' @param B number of null datasets to generate - \eqn{B=1} is usualy
#' sufficient. If computational time/power allows, a few extra
#' \eqn{B} could be helpful
#' @inheritParams lfa
#' @inheritParams sHWE
#' @examples
#' # get LFs
#' LF <- lfa(hgdp_subset, 4)
#' # look at a small (300) number of SNPs for rest of this example:
#' hgdp_subset_small <- hgdp_subset[ 1:300, ]
#' gof_4 <- sHWE(hgdp_subset_small, LF, 3)
#' LF <- lfa(hgdp_subset, 10)
#' gof_10 <- sHWE(hgdp_subset_small, LF, 3)
#' hist(gof_4)
#' hist(gof_10)
#' @return a vector of p-values for each SNP.
#' @export
sHWE <- function(X, LF, B){
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( LF ) )
        stop( '`LF` matrix is required!' )
    if ( missing( B ) )
        stop( '`B` scalar is required!' )

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
        obs_stat <- apply(X, 1, gof_stat_snp, LF)
    } else {
        # BEDMatrix case
        # write an explicit loop around the genotype matrix
        # questions: is it better to write a simple loop (one locus at the time) or to read big chunks (1000 loci at the time)?  Since `af_snp` is the bottleneck, maybe the difference is small
        
        # output vector
        obs_stat <- vector( 'numeric', m )
        for ( i in 1 : m ) {
            # get locus i genotype vector
            xi <- X[ , i ]
            # calculate and store result
            obs_stat[ i ] <- gof_stat_snp( xi, LF )
        }
    }
    
    d <- ncol(LF)
    # this already works on BEDMatrix, but produces this large matrix!
    P <- af(X, LF)
    rm(X)
    
    stat0 <- compute_nulls(P, d, B)

    # NOTE: values in either obs_stat and stat0 can be NA
    # the observed cases must be kept

    # helps us map back to original position after sort below
    obs_stat_order <- order(obs_stat)
    # sorted obs_stat has NAs removed
    obs_stat <- sort(obs_stat)
    # stat0 was a matrix, but we need it as a vector now (and sorted)
    # also, sorting this way removes NAs (it's perfect for null stats)
    stat0 <- sort(as.vector(stat0))

    p <- rep(NA, m) # initialize to NAs, to preserve obs_stat NAs
    B0 <- length(stat0)
    i1 <- 1 # index on p and obs_stat
    i2 <- 1 # index on stat0

    # i1 gets incremented in every iteration (a normal for loop?)
    # m (stopping value) is for non-NA obserations (obs_stat NA's are never in loop, so in p they stay NAs as desired)
    while(i1 <= m) {
        # i2 gets incremented until the null statistic i2 exceeds the observed statistic i1 (both in ascending order)
        while( ( i2 <= B0 ) & ( obs_stat[i1] >= stat0[i2] ) ) {
            i2 <- i2 + 1
        }
        # the p-value is 1 minus the proportion of null statistics smaller than the observed statistic
        # obs_stat_order[ i1 ] puts the value in the original position (verified in a toy case)
        p[ obs_stat_order[ i1 ] ] <- 1 - ( ( i2 - 1 ) / B0 )
        i1 <- i1 + 1
    }

    return( p )
}
