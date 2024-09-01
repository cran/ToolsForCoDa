lrpca <- function (Xcom) 
{
    n <- nrow(Xcom)
    p <- ncol(Xcom)
    Xl <- log(Xcom)
    Xclr <- cen(t(cen(t(Xl))))
    result <- svd(Xclr/sqrt(n))
    U <- result$u
    V <- result$v
    D <- diag(result$d)
    La <- D %*% D
    la <- diag(La)
    laf <- la/sum(la)
    lac <- cumsum(laf)
    decom <- rbind(la, laf, lac)
    colnames(decom) <- paste("PC", 1:p, sep = "")
    Fp <- sqrt(n) * U %*% D
    Gs <- V
    Fs <- sqrt(n) * U
    Gp <- V %*% D
    colnames(Fp) <- paste("PC", 1:p, sep = "")
    colnames(Fs) <- paste("PC", 1:p, sep = "")
    colnames(Gp) <- paste("PC", 1:p, sep = "")
    colnames(Gs) <- paste("PC", 1:p, sep = "")
    casenames <- rownames(Xcom)
    if (!is.null(casenames)) {
        rownames(Fp) <- casenames
        rownames(Fs) <- casenames
    }
    varnames <- colnames(Xcom)
    if (!is.null(varnames)) {
        rownames(Gp) <- varnames
        rownames(Gs) <- varnames
    }
    lac <- la[-length(la)]
    cn <- lac[1]/lac
    cn <- round(cn, digits = 2)
    lac <- round(lac, digits = 4)
    pc <- 1:length(lac)
    Gs <- Gs[, 1:(p - 1)]
    Gss <- round(Gs, digits = 2)
    kappalist <- data.frame(pc, cn, lac, t(Gss))
    i <- order(cn, decreasing = TRUE)
    kappalist <- kappalist[i, ]
    list(Fp = Fp, Fs = Fs, Gp = Gp, Gs = Gs, La = La, D = D, 
        decom = decom, kappalist = kappalist)
}
