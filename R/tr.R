`tr` <-
function(X) {
  if(is.matrix(X)) y <- sum(diag(X))
  if(is.vector(X)) {
    cat("warning: X is not a matrix\n")
    y <- sum(X)
  }
  return(y)
}

