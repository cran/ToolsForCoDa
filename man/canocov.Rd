\name{canocov}
\alias{canocov}
\title{
Canonical correlation analysis.
}
\description{
Function \code{canocov} performs a canonical correlation analysis. It
operates on raw data matrices, which are only centered in the
program. It uses generalized inverses and can deal with structurally
singular covariance matrices.
}
\usage{
canocov(X, Y)
}
\arguments{
  \item{X}{The n times p X matrix of observations}
  \item{Y}{The n times q Y matrix of observations}
}
\details{
  \code{canocov} computes the solution by a singular value 
  decomposition of the transformed between set covariance matrix.
}
\value{
Returns a list with the following results
\item{ccor }{ the canonical correlations }
\item{A }{ canonical weights of the X variables } 
\item{B }{ canonical weights of the Y variables } 
\item{U }{ canonical X variates } 
\item{V }{ canonical Y variates } 
\item{Fs }{ biplot markers for X variables (standard coordinates) }
\item{Gs }{ biplot markers for Y variables (standard coordinates) }
\item{Fp }{ biplot markers for X variables (principal coordinates) }
\item{Gp }{ biplot markers for Y variables (principal coordinates) }
\item{Rxu }{ canonical loadings, (correlations X variables, canonical X variates)}
\item{Rxv }{ canonical loadings, (correlations X variables, canonical Y variates)}
\item{Ryu }{ canonical loadings, (correlations Y variables, canonical X variates)}
\item{Ryv }{ canonical loadings, (correlations Y variables, canonical Y variates)}
\item{Sxu }{ covariance X variables, canonical X variates}
\item{Sxv }{ covariance X variables, canonical Y variates}
\item{Syu }{ covariance Y variables, canonical X variates}
\item{Syv }{ covariance Y variables, canonical Y variates}
\item{fitRxy }{ goodness of fit of the between-set correlation matrix }
\item{fitXs }{ adequacy coefficients of X variables } 
\item{fitXp }{ redundancy coefficients of X variables } 
\item{fitYs }{ adequacy coefficients of Y variables } 
\item{fitYp }{ redundancy coefficients of Y variables} 
}
\references{ 
   Hotelling, H. (1935) The most predictable criterion. Journal of Educational 
   Psychology (26) pp. 139-142.

   Hotelling, H. (1936) Relations between two sets of variates. Biometrika
   (28) pp. 321-377.

   Johnson, R. A. and Wichern, D. W. (2002) Applied Multivariate Statistical Analysis.
   New Jersey: Prentice Hall.   
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{cancor}}}
\examples{
set.seed(123)
X <- matrix(runif(75),ncol=3)
Y <- matrix(runif(75),ncol=3)
cca.results <- canocov(X,Y)
}
\keyword{multivariate}

