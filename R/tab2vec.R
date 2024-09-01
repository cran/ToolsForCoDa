tab2vec <- function(x) {
  yl <- names(x)
  y <- as.vector(x)
  names(y) <- yl
  return(y)
}
