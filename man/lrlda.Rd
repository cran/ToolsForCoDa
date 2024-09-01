\name{lrlda}
\alias{lrlda}
\title{Logratio Linear Discriminant Analysis}
\description{
   Function \code{lrlda} implements logratio linear discriminant analysis for compositional data, using the centred logratio
   transformation (clr) 
}
\usage{
lrlda(Xtrain, group, Xtest = NULL, divisorn = FALSE, verbose = FALSE)
}
\arguments{
  \item{Xtrain}{A compositional data set, the training data for logratio-LDA.}
  \item{group}{A categorical variable defining the groups.}
  \item{Xtest}{A compositional data set for which group prediction is sought (the test data). If no test data is supplied, the
  training data itself is classified.}
  \item{divisorn}{Use divisor "n" (\code{divisorn=TRUE}) in the calculation of covariance or use "n-1" (\code{divisorn=TRUE})}
  \item{verbose}{Print output (\code{verbose = TRUE}) or not.}
}
\details{
  Function \code{lrlda} uses the centred logratio transformation, which produces a singular covariance matrix. This singularity is
  dealt with by using a generalized inverse. When test data is supplied via argument \code{Xtest}, the scores of the linear classifier,
  the poster probabilities and the predicted classes are calculated for the test data. If no test data is supplied, these quantities
  are calculated for the training data.
}
\value{
  \item{LD}{Scores on the linear classifier for the test observations. These are also the biplot coordinates of the individuals.}
  \item{Fp}{Biplot coordinates of the group means.}
  \item{Gs}{Biplot coordinates of the variables.}  
  \item{Sp}{Pooled covariance matrix.}
  \item{Mc}{Matrix of centred clr mean vectors, one row for each group.}
  \item{S.list}{Covariance matrices of each group.}
  \item{la}{Vector of eigenvalues.}
  \item{pred}{Predicted class for the test observations.}
  \item{CM}{The confusion matrix.}
  \item{gsize}{Sample size of each group.}
  \item{Mclr}{Matrix of mean vectors for clr coordinates, one row for each group.}
  \item{prob.posterior}{Vector of posterior probabilities.}
  \item{decom}{Table with decomposition of variability as expressed by the eigenvalues.}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
   \code{\link{lrpca}},\code{\link{lrlda}}
}
\examples{
  data(Tubb)
  sampleid <- Tubb$Sample
  site     <- factor(Tubb$site)
  Oxides   <- as.matrix(Tubb[,2:10])
  rownames(Oxides) <- sampleid
  Oxides   <- Oxides/rowSums(Oxides)
  out.lda  <- lrlda(Oxides,site,verbose=FALSE)
}
\keyword{multivariate}
