\name{wav.trend.est}
\alias{wav.trend.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Wavelet Thresholding Trend Estimation of Second-Order Nonstationary Time Series
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

This function computes the linear wavelet thresholding trend estimate for a time series that may be second-order nonstationary. The function calculates the wavelet transform of the time series, sets to zero the non-boundary wavelet coffifients, then inverts the transform to give the final estimate.

}
\usage{
wav.trend.est(data, filter.number = 4, family = "DaubLeAsymm",
              max.scale = floor(log2(length(data))*0.7),
              type = "dec", boundary.handle =FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 The time series you want to esitmate the trend function of.
}
\item{filter.number}{
  Selects the index of the wavelet used in the estimation procedure. For Daubechies compactly supported wavelets the filter number is the number of vanishing moments.
}
\item{family}{
  Selects the wavelet family to use. Recommended to only use the Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
}
\item{max.scale}{
Selects the coarsest scale of the wavelet transform to analyse to. Should be a value from 1 (finest) to J-1 (coarsest), where T=2^J is the length of the time series.
}
\item{type}{
The type of wavelet transform used. By default, it is "dec" which is the standard discrete wavelet transform. Can also be "nondec", which uses a non-decimated wavelet transform.
}
\item{boundary.handle}{
 Can be TRUE or FALSE. If TRUE, the time series is boundary corrected, to get a less variable trend estimate at the boundaries of the times series. If FALSE, no boundary correction is applied.
}
}
\value{
A vector of length length(data) containing the trend estimate.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. \emph{Journal of Time Series Analysis}.
}
\author{
Euan T. McGonigle
}
\seealso{
\code{\link{wav.diff.trend.est}}
}
\examples{
# computes trend estimator of simulated linear trend time series

set.seed(1)

noise = rnorm(512)
trend = seq(from = 0, to = 5,length = 512)
x = trend+noise

trend.est = wav.trend.est(x, filter.number = 4, family = "DaubLeAsymm", boundary.handle =TRUE)

plot.ts(x, lty = 1, col = 8)
lines(trend, col = 2, lwd = 2)
lines(trend.est,col = 4, lwd = 2, lty = 2)

}

