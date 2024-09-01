dens.mvnormSisI <-
function(X,m) {
  n <- nrow(X)
  p <- ncol(X)
  X <- X - rep(1,n)%o%m
  X2 <- exp(-0.5*rowSums(X*X))
  X2 <- X2/((sqrt(2*pi))^p)
  return(X2)
}
