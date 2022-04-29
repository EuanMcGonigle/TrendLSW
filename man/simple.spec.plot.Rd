\name{simple.spec.plot}
\alias{simple.spec.plot}
\title{
Quick Plotting of Evolutionary Wavelet Spectrum
}
\description{
Plots a standard version of the EWS plot.
}
\usage{
simple.spec.plot(spectrum, bylev = FALSE)
}
\arguments{
  \item{spectrum}{
A spectrum object of class wd which you wish to plot
}
\item{bylev}{
If TRUE, plots each level of the spectrum on its own individual scaling.
}
}
\author{
E T McGonigle
}

\examples{

spec = wavethresh::cns(1024)
spec = wavethresh::putD(spec, level = 8, 1+sin(seq(from = 0, to = 2*pi,length = 1024))^2)

simple.spec.plot(spec)

}

