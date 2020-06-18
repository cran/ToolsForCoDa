\name{lrpca}
\alias{lrpca}
\title{
  Logratio principal component analysis with condition indices
}
\description{
  Function \code{lrpca} performs logratio principal component analysis. It returns the variance decomposition, principal components, biplot coordinates and a table with condition indices.
}
\usage{
lrpca(Xcom)
}
\arguments{
  \item{Xcom}{A matrix with compositions in its rows}
}
\details{
  Calculations are based on the singular value decompositon of the clr transformed compositions.
}
\value{
  \item{Fp}{matrix with principal components}
  \item{Fs}{matrix with standardized principal components}
  \item{Gp}{biplot markers for parts (principal coordinates)}
  \item{Gs}{biplot markers for parts (standard coordinates)}
  \item{La}{eigenvalues}
  \item{D}{singular values}
  \item{decom}{table with variance decomposition}
  \item{kappalist}{table with condition indices and eigenvectors}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
  \code{\link{princomp}},
}
\examples{
data(bentonites)
Ben <- bentonites[,1:8]
Ben.com <- Ben/rowSums(Ben)
out.lrpca <- lrpca(Ben.com)
}
\keyword{multivariate}
