\name{largest.kappas}
\alias{largest.kappas}
\title{
  Calculate condition indices for subcompositions
}
\description{
  Function \code{largest.kappas} calculates the condition numbers for all subcompositions of a given size, for a particular compositional data set.
}
\usage{
largest.kappas(Xcom, nparts = 3, sizetoplist = 10)
}
\arguments{
  \item{Xcom}{A data matrix with compositions in rows}
  \item{nparts}{The number of parts for the subcompositions to be analysed.}
  \item{sizetoplist}{The length of the list of the "best" subcompositions}
}
\details{
Log-ratio PCA is executed for each subcompostion, and the resulting eigenvalues and eigenvectors are stored.
}
\value{A data frame with an ordered list of subcompositions}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
X <- matrix(runif(600),ncol=6)
Xcom <- X/rowSums(X)
Results <- largest.kappas(Xcom)
}
\keyword{multivariate}

