\name{bentonites}
\alias{bentonites}
\docType{data}
\title{Isotopic and chemical compositions of bentonites}
\description{
The data consists of 14 geological samples from the US with their major oxide
composition (SiO2, Al2O3, Fe2O3, MnO, MgO, CaO, K2O, Na2O and H2O+) and
delta Deuterium and delta-18-Oxysgen (dD,d18O). 
}
\usage{data("bentonites")}
\format{
  A data frame with 14 observations on the following 11 variables.
  \describe{
    \item{\code{Si}}{a numeric vector}
    \item{\code{Al}}{a numeric vector}
    \item{\code{Fe}}{a numeric vector}
    \item{\code{Mn}}{a numeric vector}
    \item{\code{Mg}}{a numeric vector}
    \item{\code{Ca}}{a numeric vector}
    \item{\code{K}}{a numeric vector}
    \item{\code{Na}}{a numeric vector}
    \item{\code{H20}}{a numeric vector}
    \item{\code{dD}}{a numeric vector}
    \item{\code{d18O}}{a numeric vector}
  }
}
\source{
  Cadrin, A.A.J (1995), Tables 1 and 2.
  Reyment, R. A. and Savazzi, E. (1999), pp. 220-222.
}
\references{
Cadrin, A.A.J., Kyser, T.K., Caldwell, W.G.E. and Longstaffe,
F.J. (1995) Isotopic and chemical compositions of bentonites as
paleoenvironmental indicators of the Cretaceous Western Interior Seaway
Palaeogeography, Palaeoclimatology, Palaeoecology 119 pp. 301--320.
  
Reyment, R. A. and Savazzi, E. (1999) Aspects of Multivariate Statistical Analysis in Geology, Elsevier
Science B.V., Amsterdam. https://doi.org/10.1016/B978-044482568-1/50012-4
}
\examples{
data(bentonites)
}
\keyword{datasets}
