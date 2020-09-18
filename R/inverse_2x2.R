inverse_2x2 <- function(M) {
    # determinant
    denom <- M[1] * M[4] - M[2] * M[3]
    # not invertible if this is zero
    if ( denom == 0 ) {
        return( NULL )
    } else {
        # final inverse matrix
        inv <- matrix(
            c(M[4], -M[2],
              -M[3], M[1]),
            nrow = 2,
            ncol = 2
        ) / denom
        return( inv )
    }
}
