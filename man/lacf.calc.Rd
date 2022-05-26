\name{lacf.calc}
\alias{lacf.calc}
\title{
Compute localised autocovariance estimate from spectrum estimate.~
}
\description{
Function that produces an lacf object that contains the local autocovariance and autocorrelation estimates, given an input of a spectrum estimate. Provides the same functionality as the function \code{lacf} from the locits package, but user specifies the spectrum estimate in the argument.
}
\usage{
lacf.calc(x, filter.number = 10, family = "DaubLeAsymm",
                      spec.est, lag.max = NULL)
}
\arguments{
\item{x}{
	The time series you wish to analyse.
}
\item{filter.number}{
  Wavelet filter number that generated the time series.
}
\item{family}{
Wavelet family that generated the time series.
}
\item{spec.est}{
Estiamted spectrum from which the lacf estimate will be calculated.
}
\item{lag.max}{
The maximum lag of acf required. If NULL then the same default as in the regular acf function is used.
}

}
\value{
An object of class \code{lacf} which contains the autocovariance.
}

\author{
E T McGonigle
}

\seealso{

\code{\link{lacf}}

}


