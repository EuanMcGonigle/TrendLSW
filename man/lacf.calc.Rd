\name{lacf.calc}
\alias{lacf.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create lacf Object from Spectrum Estimate
%%  ~~function to do ... ~~
}
\description{
Function that produces an lacf object that contains the local autocovariance and autocorrelation estimates, given an input of a spectrum estimate. Same functionality as lacf function from locits, but user specifies the spectrum estimate in the argument.
}
\usage{
lacf.calc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Only additional argument is
  \item{altered_spec}{
Spectrum from which the lacf estimate will be calculated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{lacf}}

}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
