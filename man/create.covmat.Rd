\name{create.covmat}
\alias{create.covmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create Covariance Matrix from lacf Object
Function used to create the autocovariance matrix from an lacf object input.
}
\description{
This function used to create the autocovariance matrix from an lacf object input, which is used to calculate the confidence interval for the trend estimate.
}
\usage{
create.covmat(lacf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lacf}{
An lacf object, which can be created from a spectral estimate using the function lacf.calc.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A T x T autocovariance matrix.
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
\code{\link{lacf.calc}}, \code{\link{lacf}}

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
