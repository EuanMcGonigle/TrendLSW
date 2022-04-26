\name{calc.final.spec}
\alias{calc.final.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Auxilliary Function to Convert wd Objects
}
\description{
Auxiliary function used as a result of boundary handling. Not for general use.
}
\usage{
calc.final.spec(spec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spec}{
The spectrum object based on the boundary corrected time series of length 4*length(data).
}
}
\details{
Converts the boundary converted spectrum object back to its correct size.
}
\value{
Spectrum object of class wd.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
