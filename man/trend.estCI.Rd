\name{trend.estCI}
\alias{trend.estCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Confidene Intervals Of Wavelet-Based Trend Estimate
%%  ~~function to do ... ~~
}
\description{
A function to calculate appropriate confidence intervals for the trend estimate, based on a given trend estimate and covariance matrix of the time series.
}
\usage{
trend.estCI(trend_est, cov_mat, filter.number, family, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trend_est}{
The trend estimate of the time series.
}
  \item{cov_mat}{
The covariance matrix estimate of the time series.
}
\item{filter.number}{
  The index of the wavelet used in the estimation procedure.
}
\item{family}{
  The wavelet family used in the estimation procedure.
}
\item{alpha}{
  The significance level of the confidence interval to calculate. Default is 0.95.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list object, containing the following objects:
\item{trend.est}{The trend estimate of the input data.
}
\item{trend.var}{
The variance estimate of the trend estimate.
}
\item{upper.conf}{
The upper pointwise confidence interval for the trend estimate.
}
\item{lower.conf}{
The lower pointwise confidence interval for the trend estimate.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
E T McGonigle
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

\code{\link{wav.trend.est}}, \code{\link{lacf.calc}}, \code{\link{create.covmat}}
}
\examples{

}

