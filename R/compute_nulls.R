.compute_nulls <- function(
                           X, LF, B, max_iter = 100, tol = 1e-10,
                           m_chunk = 1000
                           ) {
    # get dimensions
    if ( methods::is( X, "BEDMatrix" ) ) {
        m_loci <- ncol( X )
        n_ind <- nrow( X )
    } else {
        n_ind <- ncol( X )
        m_loci <- nrow( X )
    }

    d <- ncol( LF )

    # initialize output now, both ways
    stats0 <- matrix(0, m_loci, B)
    
    if ( methods::is( X, "BEDMatrix" ) ) {
        # BEDMatrix trigger low-memory version
        # this low-memory version based on code from: jackstraw_BEDMatrix

        # random files to create
        # fit IAFs once, requires each of these to be created at the same time
        files_rand <- vector( 'character', B )
        for ( b in seq_len( B ) )
            files_rand[b] <- tempfile( paste0( 'sHWE_BEDMatrix_rand_b', b, '_' ), fileext = '.bed' )

        # first pass fits IAFs, simulates genotypes to write to file as we go
        # start writing in chunks for balance of memory usage and speed
        i_chunk <- 1
        while ( i_chunk <= m_loci ) {
            # calculate end of current chunk, careful not to exceed total number of loci
            i_chunk_max <- i_chunk + m_chunk - 1
            if ( i_chunk_max > m_loci )
                i_chunk_max <- m_loci
            # range of data to process
            indexes <- i_chunk : i_chunk_max
            # actual length of chunk
            m_chunk2 <- length( indexes )

            # read chunk out of BEDMatrix object, orient immediately as needed internally
            X_chunk <- t( X[ , indexes, drop = FALSE ] )

            # calculate AFs for this chunk
            P_chunk <- af( X_chunk, LF )

            for ( b in seq_len( B ) ) {
                # simulated data for this b
                X0_chunk <- matrix(
                    stats::rbinom( m_chunk2 * n_ind, 2, P_chunk ),
                    nrow = m_chunk2,
                    ncol = n_ind
                )

                # write this random data to plink BED
                genio::write_bed(
                           files_rand[b],
                           X0_chunk,
                           verbose = FALSE,
                           append = TRUE
                       )
            }

            # increment for next round
            i_chunk <- i_chunk_max + 1
        }

        # now that the files are complete, can finish processing each b
        for ( b in seq_len( B ) ) {
            # load the random genotypes
            X0 <- BEDMatrix::BEDMatrix( files_rand[b], n = n_ind, p = m_loci )
            # calculate its LFs
            LF0 <- lfa( X0, d )
            # this calculates stats correctly, even when X0 is BEDMatrix!
            stats0[, b] <- .gof_stat( X0, LF0, max_iter = max_iter, tol = tol )

            # we can delete this temp file once done!
            # NOTE: on Windows there's a peculiar issue, that these temporary files cannot be removed because BEDMatrix left them "open", silence those warnings!
            X0 <- NULL
            # try to delete, ignore warnings if it failed
            invisible( suppressWarnings( file.remove( files_rand[b] ) ) )
        }
        
    } else {
        # original, high memory solution

        # since m and n are integers, multiplying them causes a buffer overflow
        # let's multiply them as doubles, overcomes the problem
        n_data <- (n_ind + 0) * (m_loci + 0)
        # this already works on BEDMatrix, but produces this large matrix!
        P <- af(X, LF)

        for (i in seq_len(B)) {
            X0 <- matrix(stats::rbinom(n_data, 2, P), nrow=m_loci, ncol=n_ind)
            LF0 <- lfa(X0, d)
            # this calculates stats correctly, even when X0 is BEDMatrix!
            stats0[, i] <- .gof_stat(X0, LF0, max_iter=max_iter, tol=tol)
        }
    }

    return(stats0)
}
