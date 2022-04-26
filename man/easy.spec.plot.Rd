\name{easy.spec.plot}
\alias{easy.spec.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quick Plotting of Evolutionary Wavelet Spectrum
%%  ~~function to do ... ~~
}
\description{
Plots a standard version of the EWS plot.
}
\usage{
easy.spec.plot(spec, bylev = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
A spectrum object of class wd which you wish to plot
}
\item{bylev}{
If TRUE, plots each level of the spectrum on its own individual scaling.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

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
}
\examples{

spec = wavethresh::cns(1024)
spec = wavethresh::putD(spec, level = 8, 1+sin(seq(from = 0, to = 2*pi,length = 1024))^2)

easy.spec.plot(spec)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
