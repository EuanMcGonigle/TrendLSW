\name{Cmat.calc}
\alias{Cmat.calc}

\title{
  Cross Autocorrelation Wavelet Inner Product Matrix Calculation
}
\description{
Function that computes the cross autocorrelation matrix of inner products.
}
\usage{
Cmat.calc(J, gen.filter.number = 1, an.filter.number = 1,
                      gen.family = "DaubExPhase", an.family = "DaubExPhase")
}

\arguments{
  \item{J}{
 The dimension of the matrix required. Should be a positive integer.
}
\item{gen.filter.number}{
   The index of the generating wavelet used to compute the inner product matrix.
}
\item{an.filter.number}{
   The index of the analysing wavelet used to compute the inner product matrix.
}
\item{gen.family}{
   The family of generating wavelet used to compute the inner product matrix.
}
\item{an.family}{
   The family of analysing wavelet used to compute the inner product matrix.
}
}
\details{
This function computes the cross inner product matrix of the discrete non-decimated autocorrelation wavelets. This matrix is used to correct the wavelet periodogram analysed using a diffferent wavelet to the wavelet that is assumed to generate the time series. The matrix returned is the one called C in McGonigle et al. (2020).

}
\value{
A J-dimensional square matrix giving the cross inner product autocorrelation wavelet matrix.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. \emph{Journal of Time Series Analysis}.
}
\author{
Euan T. McGonigle
}

\seealso{
\code{\link{ewspec.trend}}, \code{\link{wav.diff.trend.est}}
}
\examples{
Cmat.calc(J=5, gen.filter.number = 1, an.filter.number = 2,
gen.family = "DaubExPhase", an.family = "DaubExPhase")
}

