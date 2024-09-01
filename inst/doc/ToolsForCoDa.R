## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(warnings = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(fig.width = 6, fig.height = 6) 
rm(list=ls())

## ----preinstall---------------------------------------------------------------
#install.packages("ToolsForCoDa")
library(calibrate)
library(Correlplot)
library(ToolsForCoDa)

## ----pinotnoir----------------------------------------------------------------
data(PinotNoir)
head(PinotNoir)
Aroma <- PinotNoir[,c("Aroma")]

## ----closure------------------------------------------------------------------
X  <- as.matrix(PinotNoir[,1:17])
Xc <- X/rowSums(X)
out.lrpca <- lrpca(Xc)

## ----decomposition------------------------------------------------------------
round(out.lrpca$decom,2)
plot(1:ncol(out.lrpca$decom),
     out.lrpca$decom[1,],type="b",main="Scree-plot",
     xlab="PC",ylab="Eigenvalue")

## ----biplot-------------------------------------------------------------------
lims <- jointlim(out.lrpca$Fs,2.5*out.lrpca$Gp)
bplot(out.lrpca$Fs,2.5*out.lrpca$Gp,rowlab="",collab=colnames(Xc),rowch=1,colch=NA,
      xl=lims$xlim,yl=lims$ylim,main="Covariance biplot")

pc1lab <- paste("PC1 (",toString(round(100*out.lrpca$decom[2,1],1)),"%)",sep="")
pc2lab <- paste("PC2 (",toString(round(100*out.lrpca$decom[2,2],1)),"%)",sep="")

text(1,-0.1,pc1lab,cex=0.75)
text(0.1,1.5,pc2lab,cex=0.75,srt=90)

## ----aroma--------------------------------------------------------------------
cor(Aroma,out.lrpca$Fs[,1:2])

## ----correlation--------------------------------------------------------------
logScSr <- log(Xc[,c("Cr")]/Xc[,c("Sr")])
cor(Aroma,logScSr)

## ----artificial---------------------------------------------------------------
data("Artificial")
Xsim.com <- Artificial$Xsim.com
Ysim.com <- Artificial$Ysim.com
colnames(Xsim.com) <- paste("X",1:3,sep="")
colnames(Ysim.com) <- paste("Y",1:3,sep="")

## ----ternaries, fig.width = 6, fig.height = 3---------------------------------
opar <- par(mfrow=c(1,2),mar=c(2,2,1,0)+0.5,pty="s")
par(mfg=c(1,1))
  ternaryplot(Xsim.com,pch=1)
par(mfg=c(1,2))
  ternaryplot(Ysim.com,pch=1)
par(opar)


## ----lrcco--------------------------------------------------------------------
out.lrcco <- lrcco(Xsim.com,Ysim.com)

## ----cancortab----------------------------------------------------------------
round(diag(out.lrcco$ccor),digits=3)

## ----canweights---------------------------------------------------------------
out.lrcco$A
out.lrcco$B

## ----canloadings--------------------------------------------------------------
out.lrcco$Rxu
out.lrcco$Ryv

## ----canadequacy--------------------------------------------------------------
out.lrcco$fitXs
out.lrcco$fitYs

## ----canredundancey-----------------------------------------------------------
out.lrcco$fitXp
out.lrcco$fitYp

## ----panelbiplots-------------------------------------------------------------
opar <- par(mfrow=c(2,2),mar=c(2,2,1,0)+0.5,mgp=c(2,1,0))
par(mfg=c(1,1))
#
# Figure A
#
Z <- rbind(out.lrcco$Fs,out.lrcco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1,xlab="",ylab="",main="A")
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red",angle=10,length=0.1)
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue",angle=10,length=0.1)
text(out.lrcco$Fs[,1],out.lrcco$Fs[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(out.lrcco$Gp[,1],out.lrcco$Gp[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.15
points(fa*out.lrcco$U[,1],fa*out.lrcco$U[,2])

par(mfg=c(1,2))
#
# Figure B
#
Z <- rbind(out.lrcco$Fp,out.lrcco$Gs)
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1,xlab="",ylab="",main="B")
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red",angle=10,length=0.1)
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue",angle=10,length=0.1)
text(out.lrcco$Fp[,1],out.lrcco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(out.lrcco$Gs[,1],out.lrcco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*out.lrcco$V[,1],fa*out.lrcco$V[,2])

#
# Standardizing the clr transformed data
#
Xstan.clr <- scale(clrmat(Xsim.com))
Ystan.clr <- scale(clrmat(Ysim.com))
res.stan.cco <- canocov(Xstan.clr,Ystan.clr)
par(mfg=c(2,1))
#
# Figure C
#
Z <- rbind(res.stan.cco$Fs,res.stan.cco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1,xlab="",ylab="",main="C")
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red",angle=10,length=0.1)
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue",angle=10,length=0.1)
text(res.stan.cco$Fs[,1],res.stan.cco$Fs[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.stan.cco$Gp[,1],res.stan.cco$Gp[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.2
points(fa*res.stan.cco$U[,1],fa*res.stan.cco$U[,2])
circle()
par(mfg=c(2,2))
#
# Figure D
#
Z <- rbind(res.stan.cco$Fp,res.stan.cco$Gs)
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1,xlab="",ylab="",main="D")
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red",angle=10,length=0.1)
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue",angle=10,length=0.1)
text(res.stan.cco$Fp[,1],res.stan.cco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.stan.cco$Gs[,1],res.stan.cco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*res.stan.cco$V[,1],fa*res.stan.cco$V[,2])
circle()
par(opar)

## ----bentonites---------------------------------------------------------------
data("bentonites")
head(bentonites)

## ----clrbentonites------------------------------------------------------------
X <- bentonites[,1:9]
X <- X[,-4]
X <- X/rowSums(X)
Y <- scale(bentonites[,10:11])
Xclr <- clrmat(X)
cco <- canocov(Xclr,Y)

## ----twocancor----------------------------------------------------------------
diag(cco$ccor)

## ----biplotbentonites---------------------------------------------------------
plot(cco$Fs[,1],cco$Fs[,2],col="red",pch=NA,xlab="First principal axis",
     ylab="Second principal axis",xlim=c(-1,1),ylim=c(-1,1),asp=1)
textxy(cco$Fs[,1],cco$Fs[,2],colnames(X),cex=0.75)
arrows(0,0,cco$Gp[,1],cco$Gp[,2],angle=10,col="blue",length=0.1)
arrows(0,0,cco$Fs[,1],cco$Fs[,2],angle=10,col="red",length=0.1)
text(cco$Gp[,1],cco$Gp[,2],colnames(Y),pos=c(3,1))
fa <- 0.45
points(fa*cco$U[,1],fa*cco$U[,2])
textxy(fa*cco$U[,1],fa*cco$U[,2],1:14)

## ----lognaca------------------------------------------------------------------
lnNaCa <- log(X[,c("Na")]/X[,c("Mg")])
cor(Y[,c("d18O")],lnNaCa)

## ----tubbdata-----------------------------------------------------------------
data(Tubb)
head(Tubb)
site             <- factor(Tubb$site)
Oxides           <- as.matrix(Tubb[,2:10])
rownames(Oxides) <- Tubb$Sample
Oxides           <- Oxides/rowSums(Oxides)

## ----lrldatubbdata------------------------------------------------------------
out.lrlda <- lrlda(Oxides,site)

## ----groupsizes---------------------------------------------------------------
out.lrlda$gsizes

## ----groupmeans---------------------------------------------------------------
out.lrlda$Mclr

## ----ldscores-----------------------------------------------------------------
head(out.lrlda$LD)

## ----confusion----------------------------------------------------------------
out.lrlda$CM

## ----posteriorprobs-----------------------------------------------------------
head(round(out.lrlda$prob.posterior,4))

## ----biplotcoordinates--------------------------------------------------------
Fp <- out.lrlda$Fp
Gs <- out.lrlda$Gs
LD <- out.lrlda$LD

colvec <- rep(NA,nrow(Oxides))
colvec[site=="G"]  <- "red"
colvec[site=="NF"] <- "green"
colvec[site=="W"]  <- "blue"


lims <- jointlim(LD,Fp)
opar <- par(bty="n",xaxt="n",yaxt="n")   
plot(Fp[,1],Fp[,2],asp=1,pch=17,xlab="",ylab="",col=c("red","green","blue"),
     xlim=lims$xlim,ylim=lims$ylim,cex=1.25)
points(LD[,1],LD[,2],col=colvec)
origin()
arrows(0,0,10*Gs[,1],10*Gs[,2],angle = 10, length = 0.1)
textxy(10*Gs[,1],10*Gs[,2],colnames(Oxides))
par(opar)
legend("topleft",c("G","NF","W"),col=c("red","green","blue"),pch=1,cex=0.5)

## ----lrmgoal2o3---------------------------------------------------------------
lrMgOAl2O2 <- Oxides[,c("MgO")]/Oxides[,c("Al2O3")]
boxplot(lrMgOAl2O2~site,col=c("red","green","blue"))

