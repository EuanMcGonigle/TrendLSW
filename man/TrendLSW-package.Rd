\name{TrendLSW-package}
\alias{TrendLSW-package}
\alias{TrendLSW}
\docType{package}
\title{
Wavelet Methods for Analysing Locally Stationary Time Series
}
\description{
Provides wavelet-based methods for trend, spectrum and autocovariance estimation of locally
    stationary time series.}
\details{
\tabular{ll}{
Package: \tab TrendLSW\cr
Type: \tab Package\cr
Version: \tab 0.2.0 \cr
Date: \tab 2023-11-21\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}

}
\author{
Euan T. McGonigle <e.t.mcgonigle@soton.ac.uk>, Rebecca Killick <r.killick@lancs.ac.uk>, and Matthew Nunes <m.a.nunes@bath.ac.uk>

Maintainer: Euan T. McGonigle <e.t.mcgonigle@soton.ac.uk>
}

\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
time-varying first and second-order structure of time series via wavelets
and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.

McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
locally stationary wavelet processes. \emph{Journal of Time Series
Analysis}, 43(6), 895-917.
}
\keyword{time series}
\keyword{nonstationary}
\seealso{
\code{\link{TLSW}},\code{\link{TLSWsim}},\code{\link{plot.TLSW}}
}
\examples{
# simulates an example time series and estimates its trend and evolutionary wavelet spectrum

spec <- wavethresh::cns(512)
spec <- wavethresh::putD(spec, level = 8, 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2)
# note that 8 in wavethresh is 1 in our re-numbering

trend <- seq(from = 0, to = 5, length = 512)

set.seed(1)

x <- TLSWsim(trend = trend, spec = spec)

x.TLSW <- TLSW(x)

summary(x.TLSW)
}
