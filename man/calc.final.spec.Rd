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
calc.final.spec(spec, dyadic=TRUE,data.len)
}

\arguments{
  \item{spec}{
The spectrum object based on the boundary corrected time series.
}
  \item{dyadic}{
Whether or not the time serie is a power of 2.
}
  \item{data.len}{
Length of time series.
}
}
\details{
Converts the boundary converted spectrum object back to its correct size.
}
\value{
Spectrum object of class wd.
}

\author{
E T McGonigle
}

