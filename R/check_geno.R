# stops if X is not good on some way or another
.check_geno <- function(X) {
    ret <- FALSE
    if (!is.matrix(X))
        stop("The input must be genotypes in a matrix class.")

    if (!is.integer(X[1]))
        stop("Elements of the genotype matrix should be integer.")

    classes <- names(table(as.vector(X)))
    if (!all(classes %in% c("0", "1", "2")))
        stop("Expecting genotypes to be 0, 1, and 2.")

    uniqLen <- apply(X, 1, function(x) length(unique(x)))
    if (sum(uniqLen == 1) > 1)
        stop("Remove ", uniqLen, " loci without variation across samples.")
    
    m <- nrow(X)
    n <- ncol(X)

    if (m <= n)
        stop("The genotype matrix should be tall.")
}
