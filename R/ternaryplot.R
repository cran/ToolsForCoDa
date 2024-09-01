ternaryplot <- function(X, vertexlab=colnames(X), vertex.cex = 1, pch = 19, addpoints = TRUE, grid = FALSE, gridlabels = TRUE,  ...) {
#
# Plot a ternary diagram that represents all rows of X as points.
#
  if(is.vector(X)) {
    if(length(X)!=3) {
      stop("X must have three elements")
    } else {
       X <- matrix(X,ncol=3,dimnames=list(c("1"),names(X)))
    }
  }
  X <- as.matrix(X)
  nr <- nrow(X)
  nc <- ncol(X)
  if(any(X<0)) stop("X must be non-negative")
  if(nc != 3) stop("X must have three columns")
  Xc <- X/rowSums(X)
  
  M <- matrix(c(-1/sqrt(3),0,0,1,1/sqrt(3),0),ncol=2,byrow=T)

  opar <- par(pty = "m", xpd = TRUE)
  on.exit(par(opar))
  plot(M[, 1], M[, 2], type = "n", axes = FALSE, xlab = "", 
            ylab = "", pch = 19, asp = 1, cex.main = 2, ...)
  polygon(M)
  eps <- 0.04 * vertex.cex
  Mlab <- M + matrix(c(-eps, 0, 0, eps, eps, 0), ncol = 2, 
            byrow = T)
  text(Mlab[, 1], Mlab[, 2], vertexlab, cex = vertex.cex)
  #
  # coordinates for grid
  #
  gridstyle <- "dashed"
  B1 <- rbind(c(0.8,0.2,0.0),
              c(0.6,0.4,0.0),
              c(0.4,0.6,0.0),
              c(0.2,0.8,0.0))
  E1 <- rbind(c(0.0,0.2,0.8),
              c(0.0,0.4,0.6),
              c(0.0,0.6,0.4),
              c(0.0,0.8,0.2))
  B2 <- B1[,c(3,1,2)]
  E2 <- E1[,c(3,1,2)]
  B3 <- B1[,c(2,3,1)]
  E3 <- E1[,c(2,3,1)]
  B <- rbind(B1,B2,B3)
  E <- rbind(E1,E2,E3)
  Bc <- B%*%M
  Ec <- E%*%M

  if (grid) {
    for(i in 1:nrow(Bc)) {
      segments(Bc[i,1],Bc[i,2],Ec[i,1],Ec[i,2],lty=gridstyle)
    }

    if(gridlabels) {
    
      text(Ec[1,1],Ec[1,2],"0.2",pos=4,cex=0.75)
      text(Ec[2,1],Ec[2,2],"0.4",pos=4,cex=0.75)
      text(Ec[3,1],Ec[3,2],"0.6",pos=4,cex=0.75)
      text(Ec[4,1],Ec[4,2],"0.8",pos=4,cex=0.75)
    
      text(Ec[5,1],Ec[5,2],"0.2",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[6,1],Ec[6,2],"0.4",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[7,1],Ec[7,2],"0.6",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[8,1],Ec[8,2],"0.8",pos=1,cex=0.75,srt=-120,offset=1)
    
      text(Ec[9,1],Ec[9,2],"0.2",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[10,1],Ec[10,2],"0.4",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[11,1],Ec[11,2],"0.6",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[12,1],Ec[12,2],"0.8",pos=3,cex=0.75,srt=+120,offset=1)
    }
  }

  if(addpoints) {
    Xc <- X%*%M # cartesian coordinates
    points(Xc[,1],Xc[,2],pch=pch,...)
    text(Xc[,1],Xc[,2],rownames(X))
  }
}

