### R code from vignette source 'ToolsForCoDa.Rnw'

###################################################
### code chunk number 1: ToolsForCoDa.Rnw:71-72
###################################################
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: ToolsForCoDa.Rnw:75-77 (eval = FALSE)
###################################################
## install.packages("ToolsForCoDa")
## library("ToolsForCoDa")


###################################################
### code chunk number 3: ToolsForCoDa.Rnw:84-85 (eval = FALSE)
###################################################
## vignette("ToolsForCoDa")


###################################################
### code chunk number 4: ToolsForCoDa.Rnw:95-102
###################################################
library(HardyWeinberg) # needed for making some ternary diagrams
library(ToolsForCoDa)
data("Artificial")
Xsim.com <- Artificial$Xsim.com
Ysim.com <- Artificial$Ysim.com
colnames(Xsim.com) <- paste("X",1:3,sep="")
colnames(Ysim.com) <- paste("Y",1:3,sep="")


###################################################
### code chunk number 5: ToolsForCoDa.Rnw:109-115
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,0)+0.5,mgp=c(2,1,0),pty="s")
par(mfg=c(1,1))
  HWTernaryPlot(Xsim.com,n=100,region=0,hwcurve=FALSE,vbounds=FALSE)
par(mfg=c(1,2))
  HWTernaryPlot(Ysim.com,n=100,region=0,hwcurve=FALSE,vbounds=FALSE)
par(opar)


###################################################
### code chunk number 6: ToolsForCoDa.Rnw:124-128
###################################################
Xsub.clr <- clrmat(Xsim.com)
Ysub.clr <- clrmat(Ysim.com)
colnames(Xsub.clr) <- paste("X",1:3,sep="")
colnames(Ysub.clr) <- paste("Y",1:3,sep="")


###################################################
### code chunk number 7: ToolsForCoDa.Rnw:133-135
###################################################
res.cco <- canocov(Xsub.clr,Ysub.clr)
res.cco$ccor


###################################################
### code chunk number 8: ToolsForCoDa.Rnw:141-142
###################################################
round(diag(res.cco$ccor),digits=3)


###################################################
### code chunk number 9: ToolsForCoDa.Rnw:147-149
###################################################
res.cco$A
res.cco$B


###################################################
### code chunk number 10: ToolsForCoDa.Rnw:156-158
###################################################
res.cco$Rxu
res.cco$Ryv


###################################################
### code chunk number 11: ToolsForCoDa.Rnw:163-165
###################################################
res.cco$fitXs
res.cco$fitYs


###################################################
### code chunk number 12: ToolsForCoDa.Rnw:170-172
###################################################
res.cco$fitXp
res.cco$fitYp


###################################################
### code chunk number 13: ToolsForCoDa.Rnw:181-260
###################################################

opar <- par(mfrow=c(2,2),mar=c(3,3,2,0)+0.5,mgp=c(2,1,0))
par(mfg=c(1,1))
#
# Figure A
#
Z <- rbind(res.cco$Fs,res.cco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.cco$Fs[,1],res.cco$Fs[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.cco$Gp[,1],res.cco$Gp[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.15
points(fa*res.cco$U[,1],fa*res.cco$U[,2])

par(mfg=c(1,2))
#
# Figure B
#

Z <- rbind(res.cco$Fp,res.cco$Gs)
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.cco$Fp[,1],res.cco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.cco$Gs[,1],res.cco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*res.cco$V[,1],fa*res.cco$V[,2])

par(mfg=c(2,1))
#
# Standardizing the transformed data
#

Xstan.clr <- scale(Xsub.clr)
Ystan.clr <- scale(Ysub.clr)
res.stan.cco <- canocov(Xstan.clr,Ystan.clr)

#
# Figure C
#

Z <- rbind(res.stan.cco$Fs,res.stan.cco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
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
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.stan.cco$Fp[,1],res.stan.cco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.stan.cco$Gs[,1],res.stan.cco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*res.stan.cco$V[,1],fa*res.stan.cco$V[,2])
circle()
par(opar)


###################################################
### code chunk number 14: ToolsForCoDa.Rnw:267-346
###################################################

opar <- par(mfrow=c(2,2),mar=c(3,3,2,0)+0.5,mgp=c(2,1,0))
par(mfg=c(1,1))
#
# Figure A
#
Z <- rbind(res.cco$Fs,res.cco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.cco$Fs[,1],res.cco$Fs[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.cco$Gp[,1],res.cco$Gp[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.15
points(fa*res.cco$U[,1],fa*res.cco$U[,2])

par(mfg=c(1,2))
#
# Figure B
#

Z <- rbind(res.cco$Fp,res.cco$Gs)
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.cco$Fp[,1],res.cco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.cco$Gs[,1],res.cco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*res.cco$V[,1],fa*res.cco$V[,2])

par(mfg=c(2,1))
#
# Standardizing the transformed data
#

Xstan.clr <- scale(Xsub.clr)
Ystan.clr <- scale(Ysub.clr)
res.stan.cco <- canocov(Xstan.clr,Ystan.clr)

#
# Figure C
#

Z <- rbind(res.stan.cco$Fs,res.stan.cco$Gp)
plot(Z[,1],Z[,2],type="n",xlim=c(-1,1),ylim=c(-1,1),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
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
plot(Z[,1],Z[,2],type="n",xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),asp=1)
arrows(0,0,Z[1:3,1],Z[1:3,2],col="red")
arrows(0,0,Z[4:6,1],Z[4:6,2],col="blue")
text(res.stan.cco$Fp[,1],res.stan.cco$Fp[,2],
     c(expression(X[1]),expression(X[2]),expression(X[3])))
text(res.stan.cco$Gs[,1],res.stan.cco$Gs[,2],
     c(expression(Y[1]),expression(Y[2]),expression(Y[3])),pos=c(4,3,1))
grid()
fa <- 0.25
points(fa*res.stan.cco$V[,1],fa*res.stan.cco$V[,2])
circle()
par(opar)


###################################################
### code chunk number 15: ToolsForCoDa.Rnw:357-359
###################################################
data("bentonites")
head(bentonites)


###################################################
### code chunk number 16: ToolsForCoDa.Rnw:364-369
###################################################
X <- bentonites[,1:9]
X <- X[,-4]
Y <- scale(bentonites[,10:11])
Xclr <- clrmat(X)
cco <- canocov(Xclr,Y)


###################################################
### code chunk number 17: ToolsForCoDa.Rnw:374-375
###################################################
diag(cco$ccor)


###################################################
### code chunk number 18: bentobiplot
###################################################
plot(cco$Fs[,1],cco$Fs[,2],col="red",pch=19,xlab="First principal axis",
     ylab="Second principal axis",xlim=c(-1,1),ylim=c(-1,1),asp=1)
textxy(cco$Fs[,1],cco$Fs[,2],colnames(X),cex=0.75)
points(cco$Gp[,1],cco$Gp[,2],col="blue",pch=19)
arrows(0,0,cco$Gp[,1],cco$Gp[,2])
text(cco$Gp[,1],cco$Gp[,2],colnames(Y),pos=c(3,1))
fa <- 0.45
points(fa*cco$U[,1],fa*cco$U[,2])
textxy(fa*cco$U[,1],fa*cco$U[,2],1:14)


