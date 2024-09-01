\name{ternaryplot}
\alias{ternaryplot}
\title{
  Create a Ternary Plot for three-part Compositions 
}
\description{
  Function \code{ternaryplot} accepts a matrix of three part compositions or non-negative
  counts and presents these in a ternary diagram.
}
\usage{
ternaryplot(X, vertexlab = colnames(X), vertex.cex = 1, pch = 19, addpoints = TRUE,
            grid = FALSE, gridlabels = TRUE, ...)
}
\arguments{
  \item{X}{A matrix of counts or compositions with three columns}
  \item{vertexlab}{Labels for the vertices of the tenary diagram}
  \item{vertex.cex}{Character expansion factor for vertex labels}
  \item{pch}{Plotting character for the compositions}
  \item{addpoints}{Show the compositions \code{addpoints=TRUE} or not}
  \item{grid}{Place a grid over the ternary diagram}
  \item{gridlabels}{Place grid labels or not}
  \item{\dots}{Additional arguments for the \code{points} function}
}
\value{NULL}
\author{Jan Graffelman (jan.graffelman@upc.edu)}
\examples{
data("Artificial")
Xsim.com <- Artificial$Xsim.com
colnames(Xsim.com) <- paste("X",1:3,sep="")
ternaryplot(Xsim.com)
}
\keyword{aplot}
