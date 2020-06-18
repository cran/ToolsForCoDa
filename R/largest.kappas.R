largest.kappas <- function(Xcom,nparts=3,sizetoplist=10) {
  n <- nrow(Xcom)
  if(is.null(colnames(Xcom))) colnames(Xcom) <- paste("X",1:ncol(Xcom),sep="")
  Z <- t(combn(ncol(Xcom),nparts))
  ncomb <- nrow(Z)
  CNvec <- numeric(ncomb)
  EV <- matrix(NA,nrow=ncomb,ncol=nparts)
  PartNames <- matrix(NA,nrow=ncomb,ncol=nparts)
  for(i in 1:nrow(Z)) {
    cols <- Z[i,]
    lXsub <- Xcom[,cols]  
    Xclrsub <- cen(t(cen(t(lXsub))))
    Sclr <- (1/(n-1))*t(Xclrsub)%*%Xclrsub
    out.sd <- eigen(Sclr)
    V <- out.sd$vectors
    la <- out.sd$values
    Dl <- diag(la)
    CNvec[i] <- sqrt(la[1]/la[nparts-1])
    EV[i,] <- V[,nparts-1]
    PartNames[i,] <- colnames(Xclrsub)
  }
  EV <- round(EV,digits=2)
  ii <- order(CNvec,decreasing = TRUE)
  CNvec <- CNvec[ii]
  Z <- Z[ii,]
  EV <- EV[ii,]
  PartNames <- PartNames[ii,]
  TopTen <- data.frame(CNvec,PartNames,EV)
  if(!is.null(sizetoplist)) TopTen <- TopTen[1:sizetoplist,]
  return(TopTen)
}
