lrcco <- function(X,Y) {
  nnonpositiveX <- sum(X<=0)
  nnonpositiveY <- sum(Y<=0)
  if(nnonpositiveX > 0) {
    stop("X contains non-positive elements.")
  }
  if(nnonpositiveY > 0) {
    stop("X contains non-positive elements.")
  }
  X.clr <- clrmat(X)
  Y.clr <- clrmat(Y)
  colnames(X.clr) <- paste("X",1:ncol(X),sep="")
  colnames(Y.clr) <- paste("Y",1:ncol(Y),sep="")
  out.lrcco <- canocov(X.clr,Y.clr)
  return(out.lrcco)
}
