#' @title Matrix centering and scaling
#'
#' @description
#' C routine to row-center and scale a matrix. Doesn't work with missing data.
#'
#' @param A matrix
#' @examples
#' centerscale(hgdp_subset)
#' @return matrix same dimensions `A` but row centered and scaled
#' @export
centerscale <- function(A){
    as.matrix(.Call("centerscale_c", A))
}
