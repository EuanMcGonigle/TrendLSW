\name{Atau.mat.calc}
\alias{Atau.mat.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Lagged Autocorrelation Wavelet Inner Product Calculation
}
\description{
Function that computes the matrix of lagged autocorrelation wavelet inner products.
}
\usage{
Atau.mat.calc(J, filter.number = 1, family = "DaubExPhase", lag = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{
The dimension of the matrix required. Should be a positive integer.
}
  \item{filter.number}{
The index of the wavelet used to compute the inner product matrix.
}
  \item{family}{
The family of wavelet used to compute the inner product matrix.
}
  \item{lag}{
The lag of matrix to calculate. A lage of 0 corresponds to the standrad A matrix.
}
}
\details{
This function computes the lagged inner product matrix of the discrete non-decimated autocorrelation wavelets. This matrix is used in the calculation to correct the wavelet periodogram of the differenced time series. The matrix returned is the one called A^tau in McGonigle et al. (2020).
}
\value{
A J-dimensional square matrix giving the lagged inner product autocorrelation wavelet matrix.
}
\references{
McGonigle, E. T., Killick, R., and Nunes, M. (2021). Modelling Time-Varying First and Second-Order Structure of Time Series via Wavelets and Differencing. arXiv preprint arXiv:2108.07550.
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

}
\examples{
Atau.mat.calc(J=5, filter.number = 1, family = "DaubExPhase", lag = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
