\name{ewspec.trend}
\alias{ewspec.trend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of Evolutionary Wavelet Spectrum for Non-Zero Mean Time Series
%%  ~~function to do ... ~~
}
\description{
This function computes the evolutionary wavelet spectrum (EWS) estimate from a time series that may include a trend component. The estimate is computed by taking the non-decimated wavelet transform of the time series data, squaring it; smoothing using a running mean and then correction for bias using the appropriate correction matrix. Inherits the smoohing functionality from the ewspec3 function in thr R package wavethresh.
}
\usage{
ewspec.trend(data, an.filter.number  = 10, an.family = "DaubLeAsymm",
                        gen.filter.number = an.filter.number, gen.family = an.family,
                        binwidth = floor(2*sqrt(length(data))),
                        max.scale = floor(log2(length(data))*0.7), WP.smooth = TRUE,
                        AutoReflect = TRUE, supply.mat = FALSE, mat = NULL,
                        boundary.handle = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The time series you wish to analyse.
}
\item{an.filter.number}{
The index number for the wavelet used to analyse the time series. For the "DaubExPhase" family, the filter number can be bewtween 1 to 10. For the "DaubLeAsymm" family, the filter number can be between 4 to 10. Similarly for gen.filter.number.
}
\item{an.family}{
The family of the analysing wavelet. It is recommended to use either the Daubechies Extremal Phase family, or the Daubechies Least Asymmetric family, corresponding to the "DaubExPhase" or the "DaubLeAsymm" options. Similarly for gen.family.
}
\item{gen.filter.number}{
The index number for the wavelet that generates the stochastic component of the time series.
}
\item{gen.family}{
The family of the generating wavelet.
}
\item{binwidth}{
The bin width of the running mean smoother used to smooth the raw wavelet periodogram.
}
\item{max.scale}{
The coarsest level to which the time series is analysed to. Should be a positive integer less than J, where T=2^J is the length of the time series. The default setting is 0.7J, to control for bias from the trend and boundary effects.
}
\item{WP.smooth}{
Argument that dictates if smoothing is performed on the raw wavelet periodogram.
}
\item{AutoReflect}{
As in wavethresh. Decides whether or not the time series is reflected when computing the wavelet transform. Helps estimation at the boundaries.
}
\item{supply.mat}{
Not intended for use. If TRUE, user must supply the appropriate correction matrix
}
\item{mat}{
If supply.mat is TRUE, user must supply the appropriate correction matrix used to correct the raw wavelet periodogram. Equal to solve(C).
}
\item{boundary.handle}{
 Can be TRUE or FALSE. If TRUE, the time series is boundary corrected, to get a more accurate spectrum estimate at the boundaries of the times series. If FALSE, no boundary correction is applied. Recommned to use TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

This function computes an estimate of the evolutionary wavelet spectrum of a time series that displays a smooth mean and nonstationary autocovariance. The estimation procedure is as follows:

1. The squared modulus of the non-decimated wavelet transform is computed, known as the raw wavelet periodogram. This is returned by the function.

2. The raw wavelet periodogram is smoothed using a running mean smoother.

3. The smoothed periodogram is bias corrected using the inverse of the bias matrix. The correction is applied across the finest max.scale scales. If the analyising wavelet and generating wavelet are different, this is given by the iverse of the C matrix defined in McGonigle et al. (2020). If they are the same, this is the inverse of the A matrix, defined in Nason et al. (2000). If you are unsure on the filter and wavelet choices, it is recommended to use the same wavelet for generating and analysing purposes.

The final estimate, stored in the S component, can be plotted using the plot function, please see the example below.
}
\value{
A list object, containing the following objects:
\item{S}{The evolutionary wavelet spectral estimate of the input data. This object is of class wd and so can be plotted and printed in the usual way using wavethresh functionality.
}
\item{WavPer}{
The raw wavelet periodogram of the input data. The EWS estimate (above) is the smoothed corrected version of the raw wavelet periodgram.
}
\item{SmoothWavPer}{
The smoothed, un-corrected raw wavelet periodogram of the input data.
}
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. Journal of Time Series Analysis.

Nason, G. P., von Sachs, R., and Kroisandt, G. (2000). Wavelet processes and adaptive estima- tion of the evolutionary wavelet spectrum. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 62(2):271–292.
}
\author{
E T McGonigle
}
\examples{

##---- simulates an example time series and estimates its EWS


spec = wavethresh::cns(1024)
spec = wavethresh::putD(spec, level = 8, 1+sin(seq(from = 0, to = 2*pi,length = 1024))^2)

noise = wavethresh::LSWsim(spec)
trend = seq(from = 0, to = 10,length = 1024)

x = trend+noise

spec.est = ewspec.trend(x, an.filter.number = 4, an.family = "DaubExPhase",
gen.filter.number = 1, gen.family = "DaubExPhase", max.scale = 7)

simple.spec.plot(spec.est$S)


}

