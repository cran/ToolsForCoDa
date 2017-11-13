canocov <- function (X, Y) 
{
  Xs <- cen(X)
  Ys <- cen(Y)
  Rxx <- cov(X)
  Ryy <- cov(Y)
  Rxy <- cov(X, Y)
  d <- diag(eigen(Rxx)$values)
  v <- eigen(Rxx)$vectors
  Rxxmh <- v %*% sqrt(ginv(d)) %*% t(v)
  d <- diag(eigen(Ryy)$values)
  v <- eigen(Ryy)$vectors
  Ryymh <- v %*% sqrt(ginv(d)) %*% t(v)
  K <- Rxxmh %*% Rxy %*% Ryymh
  D <- diag(svd(K)$d)
  Ah <- svd(K)$u
  Bh <- svd(K)$v
  A <- Rxxmh %*% Ah
  B <- Ryymh %*% Bh
  U <- Xs %*% A
  V <- Ys %*% B
  Fs <- Rxx %*% A
  Fp <- Rxx %*% A %*% D
  Gs <- Ryy %*% B
  Gp <- Ryy %*% B %*% D
  Rxu <- cor(Xs,U)
  Rxv <- cor(Xs,V)
  Ryv <- cor(Ys,V)
  Ryu <- cor(Ys,U)
  Sxu <- cov(Xs,U)
  Sxv <- cov(Xs,V)
  Syv <- cov(Ys,V)
  Syu <- cov(Ys,U)
  lamb <- diag(D^2)
  frac <- lamb/sum(lamb)
  cumu <- cumsum(frac)
  fitRxy <- rbind(lamb, frac, cumu)
  
  AdeXcov <- diag(Fs%*%t(Fs))
  AdeXcovr <- AdeXcov/tr(Rxx)
  fitXc <- rbind(AdeXcov,AdeXcovr)
  
  AdeX <- apply(Rxu * Rxu, 2, mean)
  AdeY <- apply(Ryv * Ryv, 2, mean)
  RedX <- apply(Rxv * Rxv, 2, mean)
  RedY <- apply(Ryu * Ryu, 2, mean)
  cAdeX <- cumsum(AdeX)
  cAdeY <- cumsum(AdeY)
  cRedX <- cumsum(RedX)
  cRedY <- cumsum(RedY)
  fitXs <- rbind(AdeX, cAdeX)
  fitXp <- rbind(RedX, cRedX)
  fitYs <- rbind(AdeY, cAdeY)
  fitYp <- rbind(RedY, cRedY)
  return(list(ccor = D, A = A, B = B, U = U, V = V, Fs = Fs, 
              Gs = Gs, Fp = Fp, Gp = Gp, Rxu = Rxu, Rxv = Rxv, Ryu = Ryu, 
              Sxu = Sxu, Sxv = Sxv, Syu = Syu, Syv = Syv,
              Ryv = Ryv, fitRxy = fitRxy, fitXs = fitXs, 
              fitXp = fitXp, fitYs = fitYs, fitYp = fitYp, fitXc = fitXc))
}
