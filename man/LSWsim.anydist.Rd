\name{LSWsim.anydist}
\alias{LSWsim.anydist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate locally stationary wavelet processes with random innovations not necessarily Gaussian.
}
\description{
Simulates a locally stationary wavelet process given a spectrum and distribution for the random innovations. Extension of the LSWsim function from wavethresh.
}
\usage{
LSWsim.anydist(spec, distribution = "Normal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
An object of class wd which contains the spectrum for simulating an LSW process.
}
  \item{distribution}{
The distribution of the random variables used to simulate the process. Can be "Normal", "Exponential", "Chisquare" or "Poisson".
}
}
\value{
A vector simulated from the spectral description given in the spec description. The returned vector will exhibit the spectral characteristics defined by spec.
}
\author{
E T McGonigle
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~

\code{\link{LSWsim}}

}
\examples{

spec = wavethresh::cns(1024)

spec = wavethresh::putD(spec,level = 8, seq(from = 2, to = 8,length = 1024))

x = LSWsim.anydist(spec, distribution = "Exponential")

plot.ts(x)

}

