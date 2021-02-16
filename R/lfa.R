#' Logistic factor analysis
#'
#' Fit a factor model of dimension `d` for binomial data, returning logistic factors (LFs).
#' Note that `d = 1` is intercept only, and for `d > 1` we compute `d - 1` singular vectors and postpend the intercept.
#'
#' Genotype matrix is expected to be a matrix of integers with values 0, 1, 2, or `NA`.
#' Note that the coding of the SNPs (which case is 0 vs 2) does not change the output.
#'
#' @param X A matrix of SNP genotypes, i.e. an integer matrix of 0's,
#' 1's, 2's and `NA`s.
#' BEDMatrix is supported.
#' Sparse matrices of class Matrix are not supported (yet).
#' @param d Number of logistic factors, including the intercept
#' @param adjustments A matrix of adjustment variables to hold fixed during estimation.
#' Number of rows must equal number of individuals in `X`.
#' These adjustments take the place of LFs in the output, so the number of columns must not exceed `d-2` to allow for the intercept and at least one proper LF to be included.
#' When present, these adjustment variables appear in the first columns of the output.
#' Not supported when `X` is a BEDMatrix object.
#' @param rspectra If `TRUE`, use [RSpectra::svds()] instead of default [trunc_svd()] or [corpcor::fast.svd()] options.
#' Ignored if `X` is a BEDMatrix object.
#' @param override Optional boolean passed to [trunc_svd()] to bypass its Lanczos bidiagonalization SVD, instead using [corpcor::fast.svd()].
#' Usually not advised unless encountering a bug in the SVD code.
#' Ignored if `X` is a BEDMatrix object.
#' @param safety Optional boolean to bypass checks on the genotype
#' matrices, which require a non-trivial amount of computation.
#' Ignored if `X` is a BEDMatrix object.
#' @param ploidy Ploidy of data, defaults to 2 for bi-allelic unphased SNPs
#' @param tol Tolerance value passed to [trunc_svd()]
#' Ignored if `X` is a BEDMatrix object.
#' @param m_chunk If `X` is a BEDMatrix object, number of loci to read per chunk (to control memory usage).
#' 
#' @return The matrix of logistic factors, with individuals along rows and factors along columns.
#' The intercept appears at the end of the columns, and adjustments in the beginning if present.
#' 
#' @examples
#' LF <- lfa(hgdp_subset, 4)
#' dim(LF)
#' head(LF)
#' @useDynLib lfa, .registration = TRUE
#' @export
lfa <- function(
                X,
                d,
                adjustments = NULL,
                override = FALSE,
                safety = FALSE,
                rspectra = FALSE,
                ploidy = 2,
                tol = 1e-13,
                m_chunk = 1000 # gave good performance in tests
                ) {
    if ( missing( X ) )
        stop( 'Genotype matrix `X` is required!' )
    if ( missing( d ) )
        stop( 'Dimension number `d` is required!' )

    # check class
    is_BEDMatrix <- FALSE
    if ( "BEDMatrix" %in% class(X) ) {
        is_BEDMatrix <- TRUE
    } else if ( !is.matrix( X ) )
        stop( '`X` must be a matrix!' )
    
    # data dimensions
    if ( is_BEDMatrix ) {
        # matrix is transposed in this case
        # only this one is needed until the pure BEDMatrix code kicks in
        n <- nrow(X)
    } else {
        m <- nrow(X)
        n <- ncol(X)
    }

    # check for d validity
    if( !is.numeric(d) ) {
        stop("d must be numeric")
    } else if( d != as.integer(d) ) {
        stop("d should be integer")
    } else if( d < 1 ) {
        stop("d should be at least 1")
    } else if( d == 1 ) {
        # return intercept column vector only
        return( matrix(1, n, 1) )
    } else if( d > 1 ) {
        d <- d - 1 #for the svd stuff
    }

    # check adjustments vars
    if ( !is.null( adjustments ) ){
        if ( is_BEDMatrix )
            stop( '`adjustments` are not supported when `X` is a BEDMatrix object!' )
        if ( !is.matrix( adjustments ) )
            stop( '`adjustments` must be a matrix!' )
        if ( nrow( adjustments ) != n )
            stop( "`adjustments` needs to have same number of rows as individuals" )
        if ( ncol( adjustments ) >= d )
            stop("need to estimate at least one non-adjustment logistic factor")
        if ( anyNA( adjustments ) )
            stop( '`adjustments` must not have missing values!' )
    }

    if ( is_BEDMatrix )
        return( lfa_BEDMatrix( X, d, ploidy = ploidy, m_chunk = m_chunk ) )
    # else continue
    
    # check data if asked to do it
    if ( safety )
        check_geno(X)

    if (!rspectra) {
        # a mysterious parameter for trunc_svd
        adjust <- 8
        if ( n - d < 10 )
            adjust <- n-d-1
    }
    
    # index the missing values
    NA_IND <- is.na(X)
    
    # center the matrix...
    mean_X <- rowMeans(X, na.rm=TRUE)
    norm_X <- X - mean_X

    # ...then 'impute'
    norm_X[NA_IND] <- 0

    # first SVD
    if (rspectra) {
        mysvd <- RSpectra::svds( norm_X, d )
    } else {
        mysvd <- trunc_svd(
            norm_X,
            d = d,
            adjust = adjust,
            tol = tol,
            override = override
        )
    }
    rm(norm_X)
    
    # for diag, pass dimension so it doesn't make mistakes
    # problem is diag( 5.5 ) returns a matrix with 5 rows/cols, only diag( 5.5, d = 1 ) is as desired
    D <- diag( mysvd$d, nrow = d )
    U <- mysvd$u
    V <- mysvd$v
    rm(mysvd)
    
    # form projection
    z <- U %*% D %*% t(V)
    z <- z + mean_X
    z <- z / ploidy
    rm(U); rm(D); rm(V) 

    #The .Call() is equivalent to the following lines of R code:
    #
    #zmin <- apply(z, 1, min)
    #zmax <- apply(z, 1, max)
    #indexes_loci_keep  <- (zmax<(1-ploidy/n)) & (zmin>(ploidy/n))
    indexes_loci_keep <- as.logical(.Call("lfa_threshold", z, 1/(ploidy*n)))
    z <- z[ indexes_loci_keep, ]
    z <- log( z / (1-z) )

    norm_z <- centerscale( z )
    
    # regress out adjustment vars, if relevant
    if ( !is.null( adjustments ) ) {
        norm_z <- t( stats::residuals( stats::lm( t(norm_z) ~ adjustments - 1 ) ) )
        d <- d - ncol( adjustments )
    }

    # second SVD yields the logistic factors
    if (rspectra) {
        v <- RSpectra::svds( norm_z, d )$v
    } else {
        v <- trunc_svd(
            norm_z,
            d = d,
            adjust = adjust,
            tol = tol,
            override = override
        )$v
    }
    # add intercept column last
    v <- cbind( v, 1 )
    # add adjustment variables first
    if ( !is.null( adjustments ) )
        v <- cbind( adjustments, v )
    
    return(v)
}


