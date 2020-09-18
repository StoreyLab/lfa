#' @title Truncated singular value decomposition
#'
#' @description
#' Truncated SVD
#'
#' @details
#' Performs singular value decomposition but only returns the first 
#' \code{d} singular vectors/values. The truncated SVD utilizes Lanczos
#' bidiagonalization. See references.
#' 
#' This function was modified from the package irlba 1.0.1 under 
#' GPL. Replacing the \code{crossprod()} calls with the C wrapper to 
#' \code{dgemv} is a dramatic difference in larger datasets. Since the 
#' wrapper is technically not a matrix multiplication function, it seemed 
#' wise to make a copy of the function.
#' @param A matrix to decompose
#' @param d number of singular vectors
#' @param adjust extra singular vectors to calculate for accuracy
#' @param tol convergence criterion
#' @param V optional initial guess
#' @param seed seed
#' @param ltrace debugging output
#' @param override TRUE means we use fast.svd instead of the iterative
#' @param force If TRUE, forces the Lanczos algorithm to be used on all datasets (usually fast.svd is used on small datasets and other conditions)
#' algorithm (useful for small data or very high d).
#' @return list with singular value decomposition.
trunc_svd <- function (A, d, adjust = 3, tol = 1e-10, V = NULL, 
                       seed=NULL, ltrace = FALSE, override=FALSE,
                       force = FALSE
                       ) {
    if ( missing( A ) )
        stop( 'Input matrix `A` is required!' )
    if ( missing( d ) )
        stop( 'Dimension number `d` is required!' )
    
    if (!is.null(seed))
        set.seed(seed)
    maxit <- 1000
    eps <- .Machine$double.eps
    
    m <- nrow(A)
    n <- ncol(A)

    if ( !force ) {
        # uses fast.svd() instead if approximate conditions are satisified
        # this mostly assumes small matrices (memory not an issue anyway)
        # d > n/20: this could be for large matrices, but large d too (when n is very large, this is unlikely)
        if((log10(m)+log10(n)) <= 6 || m < 1000 || n < 100 || d > n/20 || override){
            # this check clarifies problems in toy cases, where m becomes too small after `lfa_threshold`
            if ( d > min(m, n) )
                stop("d must be less than min(m,n)")
            mysvd <- corpcor::fast.svd(A)
            return(
                list(
                    d = mysvd$d[ 1:d ],
                    u = mysvd$u[ , 1:d, drop = FALSE ],
                    v = mysvd$v[ , 1:d, drop = FALSE ],
                    iter = 0
                )
            )
        }
    }
    
    # *adjust* d
    d_org <- d # remember original value
    d <- d+adjust
    if (m < n)
        stop("expecting tall or sq matrix")
    if (d <= 0) 
        stop("d must be positive")
    if (d > min(m, n)) 
        stop("d must be less than min(m,n)-adjust")
    if (tol < 0) 
        stop("tol must be non-negative")
    if (maxit <= 0) 
        stop("maxit must be positive")
    
    m_b <- 3
    if (m_b >= min(n, m)) {
        m_b <- floor(min(n, m) - 0.1)
        if (m_b - d - 1 < 0) {
            adjust <- 0
            d <- m_b - 1
        }
    }
    if (m_b - d - 1 < 0) 
        m_b <- ceiling(d+1+0.1)
    if (m_b >= min(m, n)) {
        m_b <- floor(min(m, n) - 0.1)
        adjust <- 0
        d <- m_b - 1
    }

    if (tol < eps) 
        tol <- eps
    
    W <- matrix(0, m, m_b)
    G <- matrix(0, n, 1)
    if (is.null(V)) {
        V <- matrix(0, n, m_b)
        V[, 1] <- stats::rnorm(n)
    }
    else {
        V <- cbind(V, matrix(0, n, m_b - ncol(V)))
    }
    B <- NULL
    Bsz <- NULL
    eps23 <- eps^(2/3)
    I <- NULL
    J <- NULL
    iter <- 1
    R_F <- NULL
    Smax <- 1
    Smin <- NULL
    SVTol <- min(sqrt(eps), tol)

    while (iter <= maxit) {
        j <- 1
        #normalize starting vector
        if (iter == 1) 
            V[, 1] <- V[, 1, drop = FALSE]/norm(V[, 1, drop = FALSE])
        else j <- d+1
        
        # compute W=AV using mv
        W[,j] <- as.matrix(.Call("mv_c", A, V[,j,drop=FALSE]))
        
        #orthogonalize W
        if (iter != 1) {
            W[, j] <- orthog(W[, j, drop = FALSE], W[, 1:j - 
                                                       1, drop = FALSE])
        }
        
        #normalize W and check for dependent vectors
        S <- norm(W[, j, drop = FALSE]) #L_2 norm
        if ((S < SVTol) && (j == 1)) 
            stop("starting vector near the null space")
        if (S < SVTol) { #check if enters??
            W[, j] <- stats::rnorm(nrow(W))
            W[, j] <- orthog(W[, j, drop = FALSE], W[, 1:j - 
                                                       1, drop = FALSE])
            W[, j] <- W[, j, drop = FALSE]/norm(W[, j, drop = FALSE])
            S <- 0
        }
        else W[, j] <- W[, j, drop = FALSE]/S
        
        #lanczos steps
        while (j <= m_b) {
            G <- as.matrix(.Call("tmv_c", A, W[,j,drop=FALSE]))

            G <- G - S * V[, j, drop = FALSE]
            
            #orthog
            G <- orthog(G, V[, 1:j, drop = FALSE])
            
            #while not the 'edge' of the bidiagonal matrix
            if (j+1 <= m_b) {
                R <- norm(G)
                #check for dependence
                if (R <= SVTol) {
                    G <- matrix(stats::rnorm(nrow(V)),nrow(V),1)
                    G <- orthog(G, V[, 1:j, drop = FALSE])
                    V[, j+1] <- G/norm(G)
                    R <- 0
                }
                else V[, j+1] <- G/R
                
                #make block diag matrix
                if (is.null(B)) B <- cbind(S, R)
                else B <- rbind(cbind(B, 0), c(rep(0, j-1), S, R))

                W[,j+1] <- as.matrix(.Call("mv_c", A, V[,j+1,drop = FALSE]))
                
                #reorthog
                W[, j+1] <- W[, j+1, drop = FALSE] - W[, j, drop = FALSE] * R
                
                #orthog
                if (iter == 1) 
                    W[, j+1] <- orthog(W[, j+1, drop = FALSE], W[, 1:j, drop = FALSE])
                S <- norm(W[, j+1, drop = FALSE])

                if (S <= SVTol) {
                    W[, j+1] <- stats::rnorm(nrow(W))
                    W[, j+1] <- orthog(W[, j+1, drop = FALSE], W[, 1:j, drop = FALSE])
                    W[, j+1] <- W[, j+1, drop = FALSE]/norm(W[, j+1, drop = FALSE])
                    S <- 0
                }
                else W[, j+1] <- W[, j+1, drop = FALSE]/S
            }
            else {
                #add block
                B <- rbind(B, c(rep(0, j - 1), S))
            }
            j <- j+1
        }
        
        #compute SVD of bidiag matrix 
        Bsz <- nrow(B)
        R_F <- norm(G)
        G <- G/R_F
        Bsvd <- svd(B)
        
        #print(rev(sort(sapply(ls(), function (object.name) object.size(get(object.name))))))

        if (iter == 1) {
            Smax <- Bsvd$d[1]
            Smin <- Bsvd$d[Bsz]
        }
        else {
            Smax <- max(Smax, Bsvd$d[1])
            Smin <- min(Smin, Bsvd$d[Bsz])
        }
        Smax <- max(eps23, Smax)
        
        #compute residuals
        R <- R_F * Bsvd$u[Bsz, , drop = FALSE]
        
        #check for convergence
        ct <- convtests(Bsz, tol, d_org, Bsvd$u, Bsvd$d, Bsvd$v, 
                        abs(R), d, SVTol, Smax)
        d <- ct$d
        
        #break criterion
        if (ct$converged) 
            break
        if (iter >= maxit) 
            break
        
        #next step in iteration --- re-initialize starting V and B
        V[, 1:(d+dim(G)[2])] <- cbind(V[, 1:(dim(Bsvd$v)[1]), drop = FALSE] %*% Bsvd$v[, 1:d, drop = FALSE], G)
        B <- cbind(diag(Bsvd$d[1:d, drop = FALSE]), R[1:d, drop = FALSE])
        
        #update left SVd
        W[, 1:d] <- W[, 1:(dim(Bsvd$u)[1]), drop = FALSE] %*% 
            Bsvd$u[, 1:d, drop = FALSE]
        iter <- iter+1
    }
    d <- Bsvd$d[1:d_org]
    u <- W[, 1:(dim(Bsvd$u)[1]), drop = FALSE] %*% Bsvd$u[, 1:d_org, 
                                                          drop = FALSE]
    v <- V[, 1:(dim(Bsvd$v)[1]), drop = FALSE] %*% Bsvd$v[, 1:d_org, 
                                                          drop = FALSE]
    return(list(d = d, u = u, v = v, iter = iter))
}

norm <- function(x)
    return( sqrt( drop( crossprod( x ) ) ) )

orthog <- function(Y, X)
    return( Y - X %*% crossprod( X, Y ) )

