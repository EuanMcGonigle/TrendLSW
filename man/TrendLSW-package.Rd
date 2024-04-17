\name{TrendLSW-package}
\alias{TrendLSW-package}
\alias{TrendLSW}
\docType{package}
\title{
Wavelet Methods for Analysing Locally Stationary Time Series
}
\description{
Provides wavelet-based methods for trend, spectrum and autocovariance estimation of locally
    stationary time series.  See \code{\link{TLSW}} for the main estimation function.}
\details{
\tabular{ll}{
Package: \tab TrendLSW\cr
Type: \tab Package\cr
Version: \tab 1.0.0 \cr
Date: \tab 2024-04-17\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Euan T. McGonigle <e.t.mcgonigle@soton.ac.uk>, Rebecca Killick <r.killick@lancs.ac.uk>, and Matthew Nunes <m.a.nunes@bath.ac.uk>

Maintainer: Euan T. McGonigle <e.t.mcgonigle@soton.ac.uk>
}

\references{
Spectral estimation with differencing/nonlinear trend estimator: McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling
time-varying first and second-order structure of time series via wavelets
and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.

Spectral estimation in presence of trend/linear trend estimator: McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend
locally stationary wavelet processes. \emph{Journal of Time Series
Analysis}, 43(6), 895-917.

LSW processes without trend: Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and
adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of
the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.

lacf estimation without trend: Nason, G. P. (2013). A test for second-order stationarity and approximate
confidence intervals for localized autocovariances for locally stationary
time series. \emph{Journal of the Royal Statistical Society: Series B
(Statistical Methodology)}, \bold{75(5)}, 879--904.
}
\keyword{time series}
\keyword{nonstationary}
\keyword{wavelet}
\seealso{
\code{\link{TLSW}},\code{\link{TLSWsim}},\code{\link{plot.TLSW}}
}
\examples{
# simulates an example time series and estimates its trend and evolutionary wavelet spectrum

spec <- matrix(0, nrow = 9, ncol = 512)
spec[1,] <- 1 + sin(seq(from = 0, to = 2 * pi, length = 512))^2

trend <- seq(from = 0, to = 5, length = 512)

set.seed(1)

x <- TLSWsim(trend = trend, spec = spec)

x.TLSW <- TLSW(x)

summary(x.TLSW)

plot(x.TLSW)
}
