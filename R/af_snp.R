#' @title Allele frequencies for SNP
#' @description Computes individual-specific allele frequencies for a 
#' single SNP.
#' @inheritParams af
#' @param snp vector of 0's, 1's, and 2's
#' @return vector of allele frequencies
#' @export
af_snp <- function(snp, LF){
    if ( missing( snp ) )
        stop( '`snp` is required!' )
    if ( missing( LF ) )
        stop( "`LF` matrix is required!" )
    
    # dimensions should agree
    if ( length(snp) != nrow(LF) )
        stop( 'Number of individuals in `snp` must equal number of individuals (rows) in `LF`' )
    
    # can only regress with non-NA individuals
    indexes_keep <- !is.na(snp)
    # get coefficients from logreg 
    betas <- lreg(
        snp[ indexes_keep ],
        LF[ indexes_keep, , drop = FALSE ]
    )

    # get allele frequencies for all individuals
    # NOTE: though some `snp` values may be NA, no `LF` are NA, and no `betas` either, so this effectively imputes the missing genotypes!
    est <- .Call("mv_c", LF, betas)
    # this is a tricky step, in that (without this trick) very large values can result in NaN's (i.e. est==1000)
    # oddly, very large negative values present no problem
    af <- ifelse( est > 100, 1, exp(est)/(1+exp(est)) )
    return( af )
}
