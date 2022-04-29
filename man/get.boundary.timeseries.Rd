\name{get.boundary.timeseries}
\alias{get.boundary.timeseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Boundary Extended Time Series
%%  ~~function to do ... ~~
}
\description{
A function to calclate the boundary extended time series, to be used within the ewspec.trend and wav.trend.est functions.
}
\usage{
get.boundary.timeseries(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The time series used to calculate the boundary extended version.
}
}
\value{
A vector of 4 times the length of the input vector.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. Journal of Time Series Analysis.
}
\author{
E T McGonigle
}
\examples{

  trend = seq(from = 0, to = 2, length = 400)
  set.seed(1)

  x = rnorm(400) + trend

  x.b = get.boundary.timeseries(x)

  plot.ts(x.b)
}

