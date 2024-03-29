% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Atau.mat.calc.R
\name{Atau.mat.calc}
\alias{Atau.mat.calc}
\title{Lagged Autocorrelation Wavelet Inner Product Calculation}
\usage{
Atau.mat.calc(J, filter.number = 1, family = "DaubExPhase", lag = 1)
}
\arguments{
\item{J}{The dimension of the matrix required. Should be a positive integer.}

\item{filter.number}{The index of the wavelet used to compute the inner
product matrix.}

\item{family}{The family of wavelet used to compute the inner product
matrix.}

\item{lag}{The lag of matrix to calculate. A lag of 0 corresponds to the
matrix \eqn{A} defined in Nason et al. (2000).}
}
\value{
A J-dimensional square matrix giving the lagged inner product
autocorrelation wavelet matrix.
}
\description{
Internal function for computing the matrix of lagged autocorrelation wavelet inner
products. This is not intended for general use by regular users of the package.
}
\details{
Computes the lagged inner product matrix of the discrete
non-decimated autocorrelation wavelets. This matrix is used in the
calculation to correct the wavelet periodogram of the differenced time
series. With \code{lag} \eqn{= \tau}, the matrix returned is the matrix \eqn{A^\tau} in McGonigle et al. (2022).
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
time-varying first and second-order structure of time series via wavelets
and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.

Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and
adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of
the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.
}
\seealso{
\link{TLSW}
}
\keyword{internal}
