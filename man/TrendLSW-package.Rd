\name{TrendLSW-package}
\alias{TrendLSW-package}
\alias{TrendLSW}
\docType{package}
\title{
\packageTitle{TrendLSW}
}
\description{
\packageDescription{TrendLSW}
}
\details{
\tabular{ll}{
Package: \tab TrendLSW \cr
Type: \tab Package\cr
Version: \tab 0.1.0 \cr
Date: \tab 2022-06-01 \cr
License: \tab GPL (>=3) \cr
LazyLoad: \tab yes\cr
}
}
\author{
\packageAuthor{TrendLSW}

Maintainer: \packageMaintainer{TrendLSW}
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. \emph{Journal of Time Series Analysis}.

McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling time-varying first and second-order structure of time series via wavelets and differencing. \emph{arXiv preprint} \href{https://arxiv.org/abs/2108.07550}{arXiv:2108.07550}.

Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{62(2)}, 271--292.
}
\seealso{

\code{\link{ewspec.trend}}, \code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}, \code{\link{wav.diff.trend.est}}

}
\keyword{wavelet}
\examples{
# simulates an example time series and estimates its evolutionary wavelet spectrum

spec = wavethresh::cns(512)
spec = wavethresh::putD(spec, level = 8, 1+sin(seq(from = 0, to = 2*pi,length = 512))^2)

noise = wavethresh::LSWsim(spec)
trend = seq(from = 0, to = 5,length = 512)

x = trend+noise

spec.est = ewspec.trend(x, an.filter.number = 4, an.family = "DaubExPhase",
gen.filter.number = 1, gen.family = "DaubExPhase")

simple.spec.plot(spec.est$S)

# estimate trend function

trend.est = wav.trend.est(x, filter.number = 4, family = "DaubLeAsymm", boundary.handle =TRUE)

plot.ts(x, lty = 1, col = 8)
lines(trend, col = 2, lwd = 2)
lines(trend.est,col = 4, lwd = 2, lty = 2)


}
