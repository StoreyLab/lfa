# generate random data for tests

# data dimensions
n_ind <- 10
m_loci <- 300
# total data size
n_data <- n_ind * m_loci
# add missingness
miss <- 0.1

# completely unstructured genotypes
# create ancestral allele frequencies
p_anc <- runif( m_loci )
# create genotypes
X <- rbinom( n_data, 2, p_anc )
# add missing values
X[ sample( n_data, n_data * miss ) ] <- NA
# turn into matrix
X <- matrix( X, nrow = m_loci, ncol = n_ind )

# to have a reasonable dataset always, remove fixed loci and all NA loci
# first remove loci that are entirely NA (with just 10 indiviuals, very possible)
loci_keep <- rowSums( !is.na(X) ) > 0
X <- X[ loci_keep, ]
# now identify fixed loci
p_anc_hat <- rowMeans( X, na.rm = TRUE )
loci_keep <- (0 < p_anc_hat) & (p_anc_hat < 1)
X <- X[ loci_keep, ]
# update number of loci and data size
m_loci <- nrow( X )
n_data <- n_ind * m_loci

# also create a matrix A != X without missingness, can be continuous values
# same dimensions as X
A <- matrix(
    rnorm( n_data ),
    nrow = m_loci,
    ncol = n_ind
)

test_that( "trunc_svd works, matches base::svd", {
    # expect errors when things are missing (both X and d are required)
    expect_error( trunc_svd() )
    expect_error( trunc_svd( A = A ) )
    expect_error( trunc_svd( d = 1 ) )

    # NOTE: since all dimensions are small, internally this defaults to fast.svd
    # test all d values, for completeness
    for ( force in c(FALSE, TRUE) ) {
        # Lanczos works best for small d (accuracy declines dramatically as d gets closer to n_ind)
        d_max <- if (force) n_ind / 2 else n_ind
        
        for ( d in 1 : d_max ) {
            # try to run successfully
            expect_silent( 
                obj <- trunc_svd(
                    A = A,
                    d = d,
                    force = force
                )
            )
            # test return values
            expect_true( is.list(obj) )
            expect_equal( length(obj), 4 )
            expect_equal( names(obj), c('d', 'u', 'v', 'iter') )
            # these must be matrices
            expect_true( is.matrix( obj$u ) )
            expect_true( is.matrix( obj$v ) )
            # dimensions, these are all different but obviously related
            expect_equal( length( obj$d ), d )
            expect_equal( nrow( obj$u ), m_loci )
            expect_equal( ncol( obj$u ), d )
            expect_equal( nrow( obj$v ), n_ind )
            expect_equal( ncol( obj$v ), d )
            
            # ultimate test is to compare to R's vanilla SVD (must agree!)
            obj2 <- svd( A, nu = d, nv = d )
            # svd's d is always length n_ind, must subset
            expect_equal( obj$d, obj2$d[ 1:d ] )
            # signs differ randomly, just compare absolute values
            expect_equal( abs(obj$u), abs(obj2$u) )
            expect_equal( abs(obj$v), abs(obj2$v) )
            
            # NOTE: though this would have been more precise, for some reason sign alignments didn't work well
            # signs differ randomly, align using first column of `u`
            ## # sgn has length m_loci
            ## sgn <- sign( obj$u[ , 1 ] * obj2$u[ , 1 ] )
            ## sgn[ sgn == 0 ] <- 1 # never use zeroes, just preserve (probably extremely rare)
            ## # this fixes signs, multiplies down columns, which is what we want
            ## expect_equal( obj$u, sgn * obj2$u )
            ## # sign flips are the same here, but only for a smaller number of rows
            ## expect_equal( obj$v, sgn[ 1 : n_ind ] * obj2$v )
        }
    }
})

test_that("lfa works", {
    # expect errors when things are missing (both X and d are required)
    expect_error( lfa() )
    expect_error( lfa( X = X ) )
    expect_error( lfa( d = 3 ) )
    # and when d is invalid
    expect_error( lfa( X = X, d = 'a' ) )
    expect_error( lfa( X = X, d = 5.9 ) ) # d must be integer
    expect_error( lfa( X = X, d = 0 ) ) # require d >= 1
    
    # test several d values, for completeness
    # NOTES:
    # - due to `lfa_threshold` removing too many SNPs in our toy examples, d can't be too large
    # - there's no "force" version here for `trunc_svd` (essentially only fast.svd outputs are tested, though they've all been verified to agree
    for ( d in 1 : (n_ind/2) ) {
        # test run overall
        expect_silent(
            LFs <- lfa( X = X, d = d )
        )
        expect_true( is.matrix( LFs ) )
        # test dimensions
        expect_equal( nrow( LFs ), n_ind )
        expect_equal( ncol( LFs ), d )
        # last column should always be intercept
        expect_equal( LFs[, d], rep.int(1, n_ind) )
        # nothing should be NA
        expect_true( !anyNA( LFs ) )
    }
})

test_that("lfa works with adjustments", {
    # weird thing is that adjustments take the place of LFs, so d >= ncol(adjustments) + 2!
    # this ensures there is at least the intercept and one proper LF)
    # (below we try 1 and 2 adjustments, so smallest d to test is 4)
    d <- 4
    
    # trigger errors when adjustments are the wrong type/dimensions
    # adjustments must be a matrix
    expect_error( lfa( X = X, d = d, adjustments = 1:n_ind ) )
    # adjustments rows must equal n_ind
    expect_error( lfa( X = X, d = d, adjustments = cbind( 2:n_ind ) ) )
    # adjustments columns must not equal or exceed `d-1`
    expect_error( lfa( X = X, d = d, adjustments = cbind( 1:n_ind, 1:n_ind, 1:n_ind ) ) )
    # adjustments aren't allowed to have NAs
    expect_error( lfa( X = X, d = d, adjustments = cbind( c(2:n_ind, NA) ) ) )
    
    # create random data for test
    # adjustments are matrices in general
    # try 1-column adjustments 
    adjustments1 <- cbind( rnorm( n_ind ) )
    # and 2 columns
    adjustments2 <- cbind( adjustments1, rnorm( n_ind ) )

    # repeat all tests for both
    for (adjustments in list( adjustments1, adjustments2 ) ) {
        # test run overall
        expect_silent(
            LFs <- lfa( X = X, d = d, adjustments = adjustments )
        )
        expect_true( is.matrix( LFs ) )
        # test dimensions
        expect_equal( nrow( LFs ), n_ind )
        expect_equal( ncol( LFs ), d ) # always d columns, regardless of adjustments size
        # last column should always be intercept
        expect_equal( LFs[ , d ], rep.int(1, n_ind) )
        # adjustment variables are repeated in first columns
        # (attributes differ, so use *_equivalent instead of *_equal)
        expect_equivalent( LFs[ , 1:ncol(adjustments) ], adjustments )
        # nothing should be NA
        expect_true( !anyNA( LFs ) )
    }
})

test_that( "lreg works", {
    # this core function is for data without missingness only!

    # get LFs from the full data with missingness (that's ok)
    d <- 3
    LFs <- lfa( X = X, d = d )
    # now generate a new unstructured genotype vector without missingness
    p_anc <- 0.5
    # create genotypes
    x <- rbinom( n_ind, 2, p_anc )

    # expect errors when key data is missing
    expect_error( lreg( ) )
    expect_error( lreg( x = x ) )
    expect_error( lreg( LF = LFs ) )
    
    # begin test!
    expect_silent(
        betas <- lreg( x = x, LF = LFs )
    )
    # test that coefficients are as expected
    expect_true( is.numeric( betas ) )
    expect_equal( length( betas ), d )
    expect_true( !anyNA( betas ) )
})

test_that( "af_snp works", {
    # like lreg, except NAs are handled and returns allele frequencies instead of coefficients

    # get LFs from the full data
    d <- 3
    LFs <- lfa( X = X, d = d )

    # expect errors when key data is missing
    expect_error( af_snp( ) )
    expect_error( af_snp( snp = X[ 1, ] ) )
    expect_error( af_snp( LF = LFs ) )
    # expect errors for mismatched dimensions
    # here number of individuals disagrees
    expect_error( af_snp( snp = X[ 1, ], LF = LFs[ 2:n_ind, ] ) )
    
    # begin test!
    # test a few SNPs in the same data (not all, that'd be overkill)
    m_loci_max <- 10
    for ( i in 1 : m_loci_max ) {
        xi <- X[ i, ]
        expect_silent(
            af <- af_snp( snp = xi, LF = LFs )
        )
        # test that AFs are as expected
        expect_true( is.numeric( af ) )
        expect_equal( length( af ), n_ind )
        expect_true( !anyNA( af ) )
    }
})

test_that( "af works", {
    # this is a boring wrapper around af_snp, applying it to the whole genome
    
    # get LFs from the full data
    d <- 3
    LFs <- lfa( X = X, d = d )
    
    # expect errors when key data is missing
    expect_error( af( ) )
    expect_error( af( X = X ) )
    expect_error( af( LF = LFs ) )
    
    # begin test!
    expect_silent(
        P <- af( X = X, LF = LFs )
    )
    # test that AFs are as expected
    expect_true( is.numeric( P ) )
    expect_true( is.matrix( P ) )
    expect_equal( nrow( P ), m_loci )
    expect_equal( ncol( P ), n_ind )
    expect_true( !anyNA( P ) )
})

test_that( "pca_af works", {
    # expect errors when key data is missing
    expect_error( pca_af( ) )
    expect_error( pca_af( X = X ) )
    expect_error( pca_af( d = d ) )
    
    # in all these cases dimensions are so small only fast.svd version is run, so all d possible values should work
    for ( d in 1 : n_ind ) {
        # try a successful run
        expect_silent(
            P <- pca_af( X = X, d = d )
        )
        # test that AFs are as expected
        expect_true( is.numeric( P ) )
        expect_true( is.matrix( P ) )
        expect_equal( nrow( P ), m_loci )
        expect_equal( ncol( P ), n_ind )
        expect_true( !anyNA( P ) )
    }
})

test_that( "centerscale works", {
    # use this function
    # NOTE: only works for data without missingness!
    expect_silent(
        A_cs <- centerscale(A)
    )
    # compare to expected value
    # first compute means
    x_m <- rowMeans( A )
    # now compute standard deviation, scale by it
    x_sd <- sqrt( rowSums( ( A - x_m )^2 ) / (n_ind-1) )
    A_cs2 <- ( A - x_m ) / x_sd
    expect_equal( A_cs, A_cs2 )
})

test_that( "check_geno works", {
    # our simulated data should pass this check
    expect_silent( check_geno( X ) )
    
    # now creater expected failures
    # this tests all cases implemented
    # ... if encoding is different this way
    expect_error( check_geno( X - 1 ) )
    # ... if matrix is not tall
    expect_error( check_geno( t(X) ) )
    # ... if it's a vector instead of a matrix
    expect_error( check_geno( 0:2 ) )
    # ... if there is a fixed locus
    # (create a 4x3 matrix, so it is tall, and with data in correct range otherwise)
    expect_error( check_geno( rbind(0:2, c(0,0,0), 2:0, 0:2 ) ) )
    # ... with the other continuous matrix
    expect_error( check_geno( A ) )
})

test_that( "inverse_2x2 works", {
    # this should match `solve`, except behavior is different when matrices are singular
    # create a random 2x2 matrix
    # keep drawing random data until it is not singular!
    rcondM <- 0 # condition number
    while ( rcondM < 1e-10 ) {
        M <- matrix(
            rnorm(4),
            nrow = 2,
            ncol = 2
        )
        rcondM <- rcond( M )
    }
    # invert this way
    expect_silent(
        inv1 <- inverse_2x2( M )
    )
    # invert the other way
    inv2 <- solve( M )
    # compare
    expect_equal( inv1, inv2 )

    # now make sure than an actual singular matrix returns NA as expected
    M <- matrix(
        c(rnorm(2), 0, 0),
        nrow = 2,
        ncol = 2
    )
    expect_silent(
        inv1 <- inverse_2x2( M )
    )
    expect_true( is.null( inv1 ) )
})

test_that( "gof_stat_snp works", {
    # get LFs for test
    d <- 3
    LFs <- lfa( X = X, d = d )
    
    # begin test!
    # test a few SNPs in the same data (not all, that'd be overkill)
    m_loci_max <- 10
    for ( i in 1 : m_loci_max ) {
        xi <- X[ i, ]
        expect_silent(
            stat <- gof_stat_snp( snp = xi, LF = LFs )
        )
        # validate features of the stat, which should be a scalar
        expect_equal( length(stat), 1 )
    }
})

test_that( "compute_nulls works", {
    d <- 3
    B <- 2
    # first compute LFs
    LFs <- lfa( X = X, d = d )
    # then compute allele frequencies
    P <- af( X = X, LF = LFs )
    # now test begins
    expect_silent(
        stat0 <- compute_nulls(AF = P, d = d, B = B)
    )
    # test return value
    expect_true( is.matrix( stat0 ) )
    expect_equal( nrow( stat0 ), m_loci )
    expect_equal( ncol( stat0 ), B )
})

test_that( "sHWE works", {
    # get LFs from the full data
    d <- 3
    LFs <- lfa( X = X, d = d )
    # just use default suggestion
    B <- 1

    # expect errors when key data is missing
    expect_error( sHWE( ) )
    expect_error( sHWE( X = X ) )
    expect_error( sHWE( LF = LFs ) )
    expect_error( sHWE( B = B ) )
    expect_error( sHWE( LF = LFs, B = B ) )
    expect_error( sHWE( X = X, B = B ) )
    expect_error( sHWE( X = X, LF = LFs ) )

    # now a successful run
    expect_silent(
        pvals <- sHWE( X = X, LF = LFs, B = B )
    )
    # test output dimensions, etc
    expect_equal( length( pvals ), m_loci )
    expect_true( max( pvals, na.rm = TRUE ) <= 1 )
    expect_true( min( pvals, na.rm = TRUE ) >= 0 )
})

### BEDMatrix tests

# require external packages for this...

if (
    suppressMessages(suppressWarnings(require(BEDMatrix))) &&
    suppressMessages(suppressWarnings(require(genio)))
) {
    context('lfa_BEDMatrix')
    
    # write the same data we simulated onto a temporary file
    file_bed <- tempfile('delete-me-random-test') # output name without extensions!
    genio::write_plink( file_bed, X )

    # load as a BEDMatrix object
    X_BEDMatrix <- suppressMessages(suppressWarnings( BEDMatrix( file_bed ) ))

    test_that( "covar_BEDMatrix and covar_logit_BEDMatrix work", {
        # computes not only covariance structure, but also mean vector

        # first compute data from ordinary R matrix, standard methods
        covar_direct <- covar_basic( X )
        X_mean <- rowMeans(X, na.rm = TRUE)

        # now compute from BEDMatrix object!
        expect_silent(
            obj <- covar_BEDMatrix(X_BEDMatrix)
        )
        # used "equivalent" because attributes differ, doesn't matter
        expect_equivalent( covar_direct, obj$covar )
        expect_equal( X_mean, obj$X_mean )

        # get eigendecomposition, make sure it agrees as expected with vanilla SVD
        # this is a test for whether the last `obj$covar` is scaled correctly or not
        d <- 3
        obj2 <- RSpectra::eigs_sym( obj$covar, d )
        V <- obj2$vectors
        # ultimate test is to compare to R's vanilla SVD (must agree!)
        # but have to transform X the same way as is normal
        Xc <- X - X_mean
        Xc[ is.na(Xc) ] <- 0
        obj3 <- svd( Xc, nu = d, nv = d )
        # sqrt(eigenvalues) should be singular values
        expect_equal( sqrt(obj2$values), obj3$d[ 1:d ] )
        # signs differ randomly, just compare absolute values
        expect_equal( abs(V), abs(obj3$v) )
        
        ## # this is a test of recovering U when it's not available
        ## expect_equal( abs( Xc %*% V %*% diag( 1/sqrt(obj2$values), d )), abs(obj3$u) )
        
        ## # construct projected data with proper SVD (truncated)
        ## Z <- obj3$u %*% diag( obj3$d[ 1:d ], d ) %*% t( obj3$v )
        ## # match it up with my prediction
        ## Z2 <- Xc %*% tcrossprod( V )
        ## expect_equal( Z, Z2 )

        # now test that subsequent step is also as desired
        expect_silent(
            covar_Z <- covar_logit_BEDMatrix( X_BEDMatrix, X_mean, V )
        )
        expect_silent(
            covar_Z_basic <- covar_logit_basic( X, V )
        )
        expect_equal( covar_Z, covar_Z_basic )
        
        # repeat with edge case m_chunk
        expect_silent(
            obj <- covar_BEDMatrix(X_BEDMatrix, m_chunk = 1)
        )
        expect_equivalent( covar_direct, obj$covar )
        expect_equal( X_mean, obj$X_mean )
        expect_silent(
            covar_Z <- covar_logit_BEDMatrix( X_BEDMatrix, X_mean, V, m_chunk = 1 )
        )
        expect_equal( covar_Z, covar_Z_basic )
    })

    test_that( "lfa works with BEDMatrix", {
        # large d doesn't work in toy data (see first `lfa` tests above for notes)
        for ( d in 1 : (n_ind/2) ) {
            # essentially the previously-tested version, no need to retest
            LFs <- lfa( X = X, d = d )
            # new version for BEDMatrix
            expect_silent(
                LFs2 <- lfa( X = X_BEDMatrix, d = d )
            )
            # signs vary randomly, but otherwise should match!
            expect_equal( abs(LFs), abs(LFs2) )
        }
    })

    test_that( "af works with BEDMatrix", {
        for ( d in 1 : (n_ind/2) ) {
            # setup data
            #d <- 3
            LFs <- lfa( X = X, d = d )
            # get ordinary `af` output
            P_basic <- af( X = X, LF = LFs )
            # get BEDMatrix version
            expect_silent(
                P_BM <- af( X = X_BEDMatrix, LF = LFs )
            )
            expect_equal( P_basic, P_BM )
        }
    })

    # delete temporary data when done
    genio::delete_files_plink( file_bed )
}
