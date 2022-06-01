\name{get.boundary.timeseries}
\alias{get.boundary.timeseries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Boundary Extended Time Series
%%  ~~function to do ... ~~
}
\description{
A function to calclate the boundary extended time series, to be used within the \code{\link{ewspec.trend}}, \code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}, and \code{\link{wav.diff.trend.est}} functions. Not recommended for general usage.
}
\usage{
get.boundary.timeseries(data, type = "TLSW")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The time series used to calculate the boundary extended version.
}
\item{type}{
The type of boundary handling used. Either "TLSW" for boundary handling as described in McGonigle, E. T., Killick, R., and Nunes, M. (2022) Trend locally stationary wavelet processes, or "LSW.diff" for a periodic version of this, used for the differencing-based functions.
}
}
\value{
A vector of 4 times the length of the input vector.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2022). Trend locally stationary wavelet processes. \emph{Journal of Time Series Analysis}.
}
\author{
Euan T. McGonigle
}
\seealso{
\code{\link{ewspec.trend}}, \code{\link{ewspec.diff}}, \code{\link{wav.trend.est}}, \code{\link{wav.diff.trend.est}}
}
\examples{

  trend = seq(from = 0, to = 2, length = 400)
  set.seed(1)

  x = rnorm(400) + trend

  x.b = get.boundary.timeseries(x)

  plot.ts(x.b)
}

