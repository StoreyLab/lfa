#' @title Matrix centering
#'
#' @description
#' C routine to row-center a matrix
#'
#' @param A matrix
#' @return `A` but row centered
#' @name center-deprecated
#' @usage center(A)
#' @seealso [lfa-deprecated()]
#' @keywords internal
NULL

#' @rdname lfa-deprecated
#' @section `center`:
#' For `center`, use `function(x) x - rowMeans(x)`.
#' @export
center <- function(A) {
    .Deprecated('function(x) x - rowMeans(x)')
    return(A - rowMeans(A))
}
