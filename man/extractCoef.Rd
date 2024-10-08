% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractCoef.R
\name{extractCoef}
\alias{extractCoef}
\alias{extractCoef,Fourier_series-method}
\title{Method of class Fourier_series to extract Fourier coefficients}
\usage{
extractCoef(object)

\S4method{extractCoef}{Fourier_series}(object)
}
\arguments{
\item{object}{an object of \code{\link{Fourier_series}} class.}
}
\value{
A list that contains the coefficients.
}
\description{
Method of class Fourier_series to extract Fourier coefficients
}
\examples{
fsc = Fourier_series(
           double_constant = 0.5,
           cos = c(0,0.3),
           sin = c(1,0.7),
           k_cos = 1:2,
           )
 extractCoef(fsc)
}
\author{
Heyang Ji
}
