cen <- function (X, w = rep(1, nrow(X))/nrow(X)) 
{
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    wm <- t(w) %*% X
    e <- matrix(rep(1, n), ncol = 1)
    Xc <- X - e %*% wm
    return(Xc)
}
