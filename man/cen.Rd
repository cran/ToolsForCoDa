\name{cen}
\alias{cen}
\title{centring of a data matrix}
\description{
  centres the columns of a matrix to mean zero.
}
\usage{cen(X,w=rep(1,nrow(X))/nrow(X))}
\arguments{
  \item{X}{a raw data matrix.}
  \item{w}{a vector of case weights.}
}
\value{returns a matrix}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
X<-matrix(runif(10),ncol=2)
Y<-cen(X)
print(Y)
}
\keyword{multivariate}

