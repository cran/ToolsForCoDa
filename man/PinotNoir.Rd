\name{PinotNoir}
\alias{PinotNoir}
\docType{data}
\title{
Chemical composition of Pinot Noir wines
}
\description{
Dataframe \code{PinotNoir} contains the composition of 17 chemical components for 37 Pinot Noir wines, as well as
an Aroma evaluation.
}
\usage{data("PinotNoir")}
\format{
  A data frame with 37 observations on the following 18 variables.
  \describe{
    \item{\code{Cd}}{Cadmium}
    \item{\code{Mo}}{Molybdenum}
    \item{\code{Mn}}{Manganese}
    \item{\code{Ni}}{Nickel}
    \item{\code{Cu}}{Copper}
    \item{\code{Al}}{Aluminium}
    \item{\code{Ba}}{Barium}
    \item{\code{Cr}}{Chromium}
    \item{\code{Sr}}{Strontium}
    \item{\code{Pb}}{Lead}
    \item{\code{B}}{Boron}
    \item{\code{Mg}}{Magnesium}
    \item{\code{Si}}{Silicon}
    \item{\code{Na}}{Sodium}
    \item{\code{Ca}}{Calcium}
    \item{\code{P}}{Phosphorus}
    \item{\code{K}}{Potassium}
    \item{\code{Aroma}}{Aroma evaluation}
  }
}
\source{
\doi{10.1016/S0003-2670(00)84245-2}
}
\references{
Frank, I.E. and Kowalski, B.R. (1984) Prediction of Wine Quality and Geographic Origin from Chemical Measurements by Partial Least-Squares Regression Modeling. Analytica Chimica Acta 162, pp. 241--251 \doi{10.1016/S0003-2670(00)84245-2}
}
\examples{
data(PinotNoir)
}
\keyword{datasets}
