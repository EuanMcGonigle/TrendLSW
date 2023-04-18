\name{wav.diff.trend.est}
\alias{wav.diff.trend.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wavelet Thresholding Trend Estimation of Second-Order Nonstationary Time Series
}
\description{
This function computes the wavelet thresholding trend estimate for a time series that may be second-order nonstationary. The function calculates the wavelet transform of the time series, thresholds the coefficients based on an estimate of their variance, and inverts to give the trend estimate.
}
\usage{
wav.diff.trend.est(data, spec, filter.number = 4, thresh.type = "soft",
normal = TRUE, family = "DaubLeAsymm", max.scale = floor(0.7 * log2(length(data))),
boundary.handle = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The time series you want to estimate the trend function of.
}
  \item{spec}{
You must supply the estimate of the evolutionary wavelet spectrum of the time series. This is calculated using the function ewspec.diff, and selecting the S component from the returned list object.
}
  \item{filter.number}{
 Selects the index of the wavelet used in the estimation procedure. For Daubechies compactly supported wavelets the filter number is the number of vanishing moments.
}
  \item{thresh.type}{
The type of thresholding function used. Currently only "soft" and "hard" are available. Recommended to use "soft".
}
  \item{normal}{
If TRUE, uses a threshold assuming the data are normally distributed. If FALSE, uses a larger threshold to reflect non-normality.
}
  \item{family}{
  Selects the wavelet family to use. Recommended to only use the Daubechies compactly supported wavelets DaubExPhase and DaubLeAsymm.
}
  \item{max.scale}{
Selects the number of scales of the wavelet transform to apply thresholding to. Should be a value from 1 (finest) to J-1 (coarsest), where T=2^J is the length of the time series. Recommended to use 2J/3 scales.
}
\item{boundary.handle}{
Logical variable, decides if boundary handling should be applied to the time series before estimation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function estimates the trend function of a locally stationary time series, by incorporating the evolutionary wavelet spectrum estimate in a wavelet thresholding procedure. To use this function, first compute the spectral estimate of the time series, using the function ewspec.diff.

The function works as follows:

1. The wavelet transform of the time series is calculated.

2. The wavelet coefficients are individually thresholded using the universal threshold sqrt(var*2*log T) and an estimate of their variance. The variance estimate is calculated using the spectral estimated, supplied by the user in the spec argument.

3. The inverse wavelet transform is applied to obtain the final estimate.



}
\value{
A vector of length length(data) containing the trend estimate.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Modelling time-varying first and second-order structure of time series via wavelets and differencing. \emph{Electronic Journal of Statistics}, 6(2), 4398-4448.
}
\author{
Euan T. McGonigle
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}

}
\examples{

spec = wavethresh::cns(1024, filter.number = 4)
spec = wavethresh::putD(spec, level = 8, 1+sin(seq(from = 0, to = 2*pi,length = 1024))^2)

set.seed(120)

noise = wavethresh::LSWsim(spec)
sine_trend = -2*sin(seq(from=0,to=2*pi,length=1024))-
1.5*cos(seq(from=0,to=pi,length=1024))

x = sine_trend+noise

spec.est = ewspec.diff(data = x,family = "DaubExPhase", filter.number = 4,max.scale = 7)

trend.est = wav.diff.trend.est(data = x, spec = spec.est$S)

plot.ts(x, lty = 1, col = 8)
lines(sine_trend, col = 2, lwd = 2)
lines(trend.est,col = 4, lwd = 2, lty = 2)


}

