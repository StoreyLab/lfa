#' @title Truncated singular value decomposition
#'
#' @description
#' Truncated SVD
#'
#' @details
#' Performs singular value decomposition but only returns the first `d` 
#' singular vectors/values.
#' The truncated SVD utilizes Lanczos bidiagonalization.
#' See references.
#' 
#' This function was modified from the package irlba 1.0.1 under GPL.
#' Replacing the [crossprod()] calls with the C wrapper to
#' `dgemv` is a dramatic difference in larger datasets.
#' Since the wrapper is technically not a matrix multiplication function, it
#' seemed wise to make a copy of the function.
#' 
#' @param A matrix to decompose
#' @param d number of singular vectors
#' @param adjust extra singular vectors to calculate for accuracy
#' @param tol convergence criterion
#' @param override `TRUE` means we use
#' [corpcor::fast.svd()] instead of the
#' iterative algorithm (useful for small data or very high `d`).
#' @param force If `TRUE`, forces the Lanczos algorithm to be used on all
#' datasets (usually
#' [corpcor::fast.svd()]
#' is used on small datasets or large `d`)
#' @param maxit Maximum number of iterations
#' @return list with singular value decomposition.  Has elements 'd', 'u', 'v',
#' and 'iter'
#' @examples
#' obj <- trunc_svd( hgdp_subset, 4 )
#' obj$d
#' obj$u
#' obj$v
#' obj$iter
#' @export
trunc_svd <- function(A, d, adjust = 3, tol = .Machine$double.eps,
    override = FALSE, force = FALSE, maxit = 1000) {
    if (missing(A))
        stop("Input matrix `A` is required!")
    if (missing(d))
        stop("Dimension number `d` is required!")
    if (d <= 0)
        stop("d must be positive")
    m <- nrow(A)
    n <- ncol(A)
    if (d > min(m, n))
        stop("d must be less than min(m,n)")
    if (!force) { # uses fast.svd() for small matrices or large `d`
        if ((log10(m) + log10(n)) <= 6 || m < 1000 || n < 100 || d > n/20 ||
            override) {
            mysvd <- corpcor::fast.svd(A)
            indexes <- seq_len(d)
            return(list(d = mysvd$d[indexes], u = mysvd$u[, indexes,
                drop = FALSE], v = mysvd$v[, indexes, drop = FALSE], iter = NA))
        }
    }
    d_org <- d  # remember original value
    d <- d + adjust # *adjust* d
    if (m < n)
        stop("expecting tall or sq matrix")
    if (d > min(m, n))
        stop("d must be less than min(m,n)-adjust")
    W <- matrix(0, m, d + 1)
    V <- matrix(0, n, d + 1)
    V <- .new_col_ortho_unit(V, 1)
    dat <- list()
    iter <- 1
    while (iter <= maxit) {
        dat <- .trunc_svd_iter(A, V, W, dat$B, dat$Smax, d, d_org, iter, tol)
        V <- dat$V
        W <- dat$W
        Bsvd <- dat$Bsvd
        if (dat$converged || iter >= maxit) break  # break criterion
        d <- dat$d
        V[, seq_len(d + 1)] <- cbind(V[, seq_len(nrow(Bsvd$v))] %*%
            Bsvd$v[, seq_len(d)], dat$G)
        dat$B <- cbind(diag(Bsvd$d[seq_len(d)]), dat$R[seq_len(d)])
        W[, seq_len(d)] <- W[, seq_len(nrow(Bsvd$u))] %*% Bsvd$u[, seq_len(d)]
        iter <- iter + 1
    }
    d <- Bsvd$d[seq_len(d_org)]
    u <- W[, seq_len(nrow(Bsvd$u))] %*% Bsvd$u[, seq_len(d_org)]
    v <- V[, seq_len(nrow(Bsvd$v))] %*% Bsvd$v[, seq_len(d_org)]
    return(list(d = d, u = u, v = v, iter = iter))
}

.trunc_svd_iter <- function(A, V, W, B, Smax, d, d_org, iter, tol) {
    j <- 1
    if (iter != 1) 
        j <- d + 1
    W[, j] <- .Call("mv_c", A, V[, j]) # W=AV
    if (iter != 1)
        W[, j] <- .orthog(W[, j], W[, seq_len(j) - 1])
    S <- .norm(W[, j])
    if (S < tol) { # normalize W and check for dependent vectors
        W <- .new_col_ortho_unit(W, j)
        S <- 0
    } else W[, j] <- W[, j]/S
    # lanczos steps
    while (j <= ncol(W)) {
        G <- .Call("tmv_c", A, W[, j]) - S * V[, j]
        G <- .orthog(G, V[, seq_len(j)])
        if (j + 1 <= ncol(W)) { # while not the 'edge' of the bidiag matrix
            R <- .norm(G)
            if (R <= tol) { # check for dependence
                V <- .new_col_ortho_unit(V, j + 1)
                G <- V[, j + 1]
                R <- 0
            } else V[, j + 1] <- G/R
            if (is.null(B)) {
                B <- cbind(S, R) # make block diag matrix
            } else B <- rbind(cbind(B, 0), c(rep(0, j - 1), S, R))
            W[, j + 1] <- .Call("mv_c", A, V[, j + 1]) - W[, j] * R
            if (iter == 1)
                W[, j + 1] <- .orthog(W[, j + 1], W[, seq_len(j)])
            S <- .norm(W[, j + 1])
            if (S <= tol) {
                W <- .new_col_ortho_unit(W, j + 1)
                S <- 0
            } else W[, j + 1] <- W[, j + 1]/S
        } else B <- rbind(B, c(rep(0, j - 1), S)) # add block
        j <- j + 1
    }
    Bsz <- nrow(B)
    R_F <- .norm(G)
    G <- G/R_F
    Bsvd <- svd(B) # SVD of bidiag matrix
    Smax <- max(Smax, Bsvd$d[1], tol^(2/3))
    R <- R_F * Bsvd$u[Bsz, ] # compute residuals
    ct <- .convtests(Bsz, tol, d_org, abs(R), d, Smax) # check convergence
    return(c(ct, list(V=V, W=W, B=B, G=G, R=R, Bsvd=Bsvd, Smax=Smax)))
}

# replace column with random data!
.new_col_ortho_unit <- function(W, j) {
    # new column with random data
    Wj <- stats::rnorm(nrow(W))
    # remove projection to existing data in W (cols < j).
    # Nothing to do if j==1
    if (j > 1)
        Wj <- .orthog(Wj, W[, seq_len(j-1)])
    # unit normalize and store in desired column
    W[, j] <- Wj/.norm(Wj)
    return(W) # return whole matrix
}

# these work just fine if x/X/Y are dropped to vectors
.norm <- function(x) return(sqrt(drop(crossprod(x))))
.orthog <- function(Y, X) return(Y - X %*% crossprod(X, Y))
