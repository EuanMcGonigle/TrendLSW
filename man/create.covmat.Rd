\name{create.covmat}
\alias{create.covmat}
\title{ Create Covariance Matrix from lacf Object
Function used to create the autocovariance matrix from an lacf object input.
}
\description{
This function used to create the autocovariance matrix from an lacf object input, which is used to calculate the confidence interval for the trend estimate.
}
\usage{
create.covmat(lacf)
}
\arguments{
  \item{lacf}{
An lacf object, which can be created from a spectral estimate using the function lacf.calc.
}
}
\value{
A T x T autocovariance matrix.
}
\author{
E T McGonigle
}
\seealso{

\code{\link{lacf.calc}}, \code{\link{lacf}}

}


