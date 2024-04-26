\name{z.labels}
\docType{data}
\alias{z.labels}
\title{
  Activity Labels for Human Activity Monitoring
}
\description{

This dataset gives the labelled activities recorded during the time period of observations given in the data
object \code{z.acc}.

}
\usage{
  z.labels
}
\format{ A data frame with 6 rows and 3 variables:
 \describe{
   \item{activity}{The activity recorded, either \code{"downstairs"} or \code{"upstairs"}, corresponding to walking downstairs and upstairs respectively.}
   \item{start}{the starting time of the activity.}
   \item{end}{the ending time of the activity.}
 }
}
\source{Kelly M., Longjohn R., and Nottingham, K. (2024). The UCI Machine Learning Repository. \doi{10.24432/C54G7M}.
}
\references{
  Reyes-Ortiz, J. L., Oneto, L., Sama, A., Parra, X., and Anguita, D. (2016). Transition-Aware Human Activity Recognition Using Smartphones. \emph{Neurocomputing}, 171, 754–767.
}

\keyword{datasets}
\seealso{
\code{\link{z.acc}}
}
