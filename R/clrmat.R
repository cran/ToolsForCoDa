clrmat <- function(X) {
  lX <- log(X)
  Xclr <- cen(t(cen(t(lX))))
  return(Xclr)
}
