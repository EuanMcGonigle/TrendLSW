\name{Atau.mat.calc}
\alias{Atau.mat.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Lagged Autocorrelation Wavelet Inner Product Calculation
}
\description{
Function that computes the matrix of lagged autocorrelation wavelet inner products.
}
\usage{
Atau.mat.calc(J, filter.number = 1, family = "DaubExPhase", lag = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
The dimension of the matrix required. Should be a positive integer.
}
  \item{filter.number}{
The index of the wavelet used to compute the inner product matrix.
}
  \item{family}{
The family of wavelet used to compute the inner product matrix.
}
  \item{lag}{
The lag of matrix to calculate. A lag of 0 corresponds to the standard A matrix.
}
}
\details{
This function computes the lagged inner product matrix of the discrete non-decimated autocorrelation wavelets. This matrix is used in the calculation to correct the wavelet periodogram of the differenced time series. The matrix returned is the one called A^tau in McGonigle et al. (2020).
}
\value{
A J-dimensional square matrix giving the lagged inner product autocorrelation wavelet matrix.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling time-varying first and second-order structure of time series via wavelets and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448. 
}
\author{
Euan T. McGonigle
}
\examples{
Atau.mat.calc(J=5, filter.number = 1, family = "DaubExPhase", lag = 1)
}
