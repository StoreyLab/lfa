#' @title Allele frequencies
#' @description Compute matrix of individual-specific allele frequencies
#' @inheritParams lfa
#' @param LF Matrix of logistic factors, with intercept.
#' Pass in the return value from [lfa()]!
#' @param max_iter Maximum number of iterations for logistic regression
#' @param tol Numerical tolerance for convergence of logistic regression
#' @details Computes the matrix of individual-specific allele 
#' frequencies, which has the same dimensions of the genotype matrix.
#' Be warned that this function could use a ton of memory, as the 
#' return value is all doubles. It could be wise to pass only a 
#' selection of the SNPs in your genotype matrix to get an idea for
#' memory usage. Use [gc()] to check memory usage!
#' @examples
#' LF <- lfa( hgdp_subset, 4 )
#' allele_freqs <- af( hgdp_subset, LF )
#' @return Matrix of individual-specific allele frequencies.
#' @export
af <- function(X, LF, safety = FALSE, max_iter = 100, tol = 1e-10) {
    if (missing(X))
        stop("Genotype matrix `X` is required!")
    if (missing(LF))
        stop("`LF` matrix is required!")

    # check class
    if (!is.matrix(X)) # BEDMatrix returns TRUE here
        stop("`X` must be a matrix!")

    # get dimensions
    if (methods::is(X, "BEDMatrix")) {
        m <- ncol(X)
        n <- nrow(X)
    } else {
        n <- ncol(X)
        # m not used in this case
    }

    # dimensions should agree
    if (n != nrow(LF))
        stop("Number of individuals in `X` and `LF` disagree!")
    
    if (!methods::is(X, "BEDMatrix")) {
        # usual R object behavior
        if (safety)
            .check_geno(X)
        return(t(apply(X, 1, af_snp, LF, max_iter=max_iter, tol=tol)))
    } else {
        # BEDMatrix case.
        P <- matrix(0, m, n) # init output matrix
        for (i in seq_len(m)) {
            # get locus i genotype vector
            xi <- X[, i]
            # calculate and store result
            P[i, ] <- af_snp(xi, LF, max_iter=max_iter, tol=tol)
        }
        # done!
        return(P)
    }
}

