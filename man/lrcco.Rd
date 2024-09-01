\name{lrcco}
\alias{lrcco}
\title{
Logratio Canonical Correlation Analysis
}
\description{
  Function \code{lrcco} is a wrapper function around \code{canocov}. It performs logratio canonical correlation analysis (LR-CCO) accepting
  two compositional data matrices as input.
}
\usage{
lrcco(X, Y)
}
\arguments{
  \item{X}{The matrix of X compositions}
  \item{Y}{The matrix of Y compositions}
}
\details{
  Matrices \code{X} and \code{Y} are assumed to contain positive elements only, and there rows sum to one. 
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

Graffelman, J. and Pawlowsky-Glahn, V. and Egozcue, J.J. and Buccianti, A. (2018)
Exploration of geochemical data with compositional canonical biplots, Journal of
Geochemical Exploration 194, pp. 120--133. \doi{10.1016/j.gexplo.2018.07.014}

}

\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }

\seealso{\code{\link{cancor}},\code{\link{canocov}}}

\examples{
set.seed(123)
X  <- matrix(runif(75),ncol=3)
Y  <- matrix(runif(75),ncol=3)
Xc <- X/rowSums(X) # create compositions by closure
Yc <- Y/rowSums(Y)
out.lrcco <- lrcco(X,Y)
}

\keyword{multivariate}

