#' @title LFA model goodness of fit
#' 
#' @description
#' Compute SNP-by-SNP goodness-of-fit when compared to population 
#' structure. This can be aggregated to determine genome-wide 
#' goodness-of-fit for a particular value of `d`.
#' 
#' @details
#' This function returns p-values for LFA model goodness of fit based
#' on a simulated null.
#'
#' @note Genotype matrix is expected to be a matrix of integers with
#' values 0, 1, and 2. Currently no support for missing values. Note
#' that the coding of the SNPs does not affect the algorithm.
#'
#' @param X A matrix of SNP genotypes, i.e. an integer matrix of 0's,
#' 1's, 2's and `NA`s.
#' BEDMatrix is supported.
#' @param LF matrix of logistic factors
#' @param B number of null datasets to generate, `B = 1` is usually
#' sufficient. If computational time/power allows, a few extra
#' `B` could be helpful
#' @return vector of p-values for each SNP.
#' @name model.gof-deprecated
#' @usage model.gof(X, LF, B)
#' @seealso [lfa-deprecated()]
#' @keywords internal
NULL

#' @rdname lfa-deprecated
#' @section `model.gof`:
#' For `model.gof`, use [sHWE()].
#' @export
model.gof <- function(X, LF, B) {
    .Deprecated('sHWE')
    sHWE(X, LF, B)
}
