lrlda <- function (Xtrain, group, Xtest = NULL, divisorn = FALSE, verbose = FALSE) 
{
  ng <- length(table(group))
  gsizes <- tab2vec(table(group))
  if (verbose) {
    cat("Group sizes:\n\n")
    print(gsizes)
  }
  p <- ncol(Xtrain)
  n <- nrow(Xtrain)
  dim.sol <- min(p - 1, ng - 1)
  ntrain <- nrow(Xtrain)
  Xl <- log(Xtrain)
  Xlc <- t(scale(t(Xl), scale = FALSE))
  m <- colMeans(Xlc)
  d <- data.frame(Xlc, group)
  Mclr <- aggregate(d[, 1:p], list(d$group), mean)
  Xlc <- Xlc - rep(1, n) %o% m
  if (verbose) {
    cat("Mean vectors of CLR coordinates:\n\n")
    print(Mclr)
  }
  M <- rep(1, ng) %o% m
  Mc <- as.matrix(Mclr[, -1] - M)
  rownames(Mc) <- Mclr[, 1]
  w <- tab2vec(table(group))
  w <- w/sum(w)
  Dw <- diag(w)
  Dwh <- sqrt(Dw)
  S.list <- lapply(levels(d$group), function(x) cov(d[d$group == 
                                                        x, -ncol(d)], use = "na.or.complete"))
  Sp <- matrix(0, ncol = p, nrow = p)
  if (divisorn) {
    for (i in 1:ng) {
      Sp <- Sp + w[i] * ((gsizes[i] - 1)/gsizes[i]) * S.list[[i]]
    }
  }
  else {
    for (i in 1:ng) {
      Sp <- Sp + ((gsizes[i] - 1)/(sum(gsizes) - ng)) * 
        S.list[[i]]
    }
  }
  if (verbose) {
    cat("Pooled covariance matrix:\n\n")
    print(Sp)
  }
  out <- eigen(Sp)
  Dl <- diag(out$values)
  V <- out$vectors
  Spmh <- V %*% sqrt(ginv(Dl)) %*% t(V)
  K <- Dwh %*% Mc %*% Spmh
  out.svd <- svd(K)
  U <- out.svd$u[, 1:dim.sol]
  D <- diag(out.svd$d)
  D <- D[1:dim.sol, 1:dim.sol]
  V <- out.svd$v[, 1:dim.sol]
  if(dim.sol > 1) {
    la <- diag(D * D)  
  } else {
    la <- D
  }
  laf <- la/sum(la)
  lac <- cumsum(laf)
  if(dim.sol > 1) {
    decom <- rbind(la, laf, lac)  
    colnames(decom) <- paste("LD", 1:dim.sol, sep = "")  
    rownames(decom) <- c("lambda","fraction","cumulative")
  } else if(dim.sol==1) {
    decom <- c(la, laf, lac)
    decom <- matrix(decom,ncol=1)
    colnames(decom) <- paste("LD", 1:dim.sol, sep = "")  
    rownames(decom) <- c("lambda","fraction","cumulative")
  }
  Fp <- solve(Dwh) %*% U %*% D
  rownames(Fp) <- levels(d$group)
  Gs <- V
  if (is.matrix(Gs)) {
    rownames(Gs) <- colnames(Xtrain)
  }
  else if (is.vector(Gs)) {
    names(Gs) <- colnames(Xtrain)
  }
  else {
    stop("unknown type for Gs")
  }
  LD <- Xlc %*% Spmh %*% Gs
  colnames(LD) <- paste("LD", 1:dim.sol, sep = "")
  rownames(LD) <- rownames(Xlc)
 
 
  Mld <- aggregate(LD, list(group), mean)
  rownames(Mld) <- Mld[, 1]
  Mld <- Mld[, -1]
  Mld <- as.matrix(Mld)
  if (verbose) {
    cat("Group means on the linear discriminants:\n\n")
    print(Mld)
  }
  prob.posterior <- matrix(0, nrow = nrow(LD), ncol = ng)
  for (i in 1:ng) {
    if (ncol(LD) > 1) {
      prob.posterior[, i] <- dens.mvnormSisI(LD, m = Mld[i, 
      ])
    }
    else {
      prob.posterior[, i] <- dnorm(LD, mean = Mld[i, ], 
                                   sd = 1)
    }
  }
  PostTot <- rep(0, nrow(LD))
  for (i in 1:ng) {
    PostTot <- PostTot + w[i] * prob.posterior[, i]
  }
  for (i in 1:ng) {
    prob.posterior[, i] <- w[i] * prob.posterior[, i]/PostTot
  }
  rownames(prob.posterior) <- rownames(Xlc)
  colnames(prob.posterior) <- names(gsizes)

  pred <- rep(NA, nrow(LD))

  for (i in 1:nrow(LD)) {
    pred[i] <- which.max(prob.posterior[i,])
  }

  CM <- table(group, pred)
  if (verbose) {
    cat("Confusion matrix\n\n")
    print(CM)
  }

  if (!is.null(Xtest)) {
    ntest <- nrow(Xtest)
    Xtestl <- log(Xtest)
    Xtestlc <- t(scale(t(Xtestl), scale = FALSE))
    Xtestlc <- Xtestlc - rep(1, ntest) %o% m
    LD.test <- Xtestlc %*% Spmh %*% Gs
    colnames(LD.test) <- paste("LD", 1:dim.sol, sep = "")
    rownames(LD.test) <- rownames(Xtestlc)
   
    prob.posterior.test <- matrix(0, nrow = nrow(LD.test), 
                                  ncol = ng)
    for (i in 1:ng) {
      if (ncol(LD.test) > 1) {
        prob.posterior.test[, i] <- dens.mvnormSisI(LD.test, 
                                                    m = Mld[i, ])
      }
      else {
        prob.posterior.test[, i] <- dnorm(LD.test, mean = Mld[i, 
        ], sd = 1)
      }
    }
    PostTot <- rep(0, nrow(LD.test))
    for (i in 1:ng) {
      PostTot <- PostTot + w[i] * prob.posterior.test[, 
                                                      i]
    }
    for (i in 1:ng) {
      prob.posterior.test[, i] <- w[i] * prob.posterior.test[, 
                                                             i]/PostTot
    }
    rownames(prob.posterior.test) <- rownames(Xtestlc)
    colnames(prob.posterior.test) <- names(gsizes)

    pred.test <- rep(NA, nrow(LD.test))

    for (i in 1:nrow(LD.test)) {
      pred.test[i] <- which.max(prob.posterior.test[i,])
    }

    LD   <- LD.test
    pred <- pred.test
    prob.posterior <- prob.posterior.test
  }
  return(list(LD = LD, Fp = Fp, Gs = Gs, Sp = Sp, Mc = Mc, 
              S.list = S.list, la = la, pred = pred, CM = CM, gsizes = gsizes, 
              Mclr = Mclr, prob.posterior = prob.posterior, decom = decom))
}
