\name{clrmat}
\alias{clrmat}
\title{
Centred log-ratio transformation
}
\description{
Program \code{clrmat} calculates the centred log-ratio transformation
for a matrix of compositions. 
}
\usage{
clrmat(X)
}
\arguments{
  \item{X}{A matrix of compositions}
}
\value{A matrix containing the transformed data}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\examples{
data(Artificial)
Xsim.com <- Artificial$Xsim.com
Xclr <- clrmat(Xsim.com)
}
\keyword{multivariate}

