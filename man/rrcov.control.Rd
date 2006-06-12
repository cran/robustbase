\name{rrcov.control}
\alias{rrcov.control}
%%% FIXME --- naming and more --------------
%%% -----
\title{Control object for the estimation parameters }
\description{
     Auxilary function for passing the estimation options as parameters to the
     estimation functions.

     \bold{NOTE: The name  WILL change !!!!}
}
\usage{
rrcov.control(alpha = 1/2, nsamp = 500, seed = NULL, trace = FALSE,
              use.correction = TRUE, adjust = FALSE)
}
\arguments{
  \item{alpha}{This parameter controls the size of the subsets over
    which the determinant is minimized, i.e., \code{alpha*n} observations
    are used for computing the determinant.  Allowed values are between 0.5
    and 1 and the default is 0.5. }
  \item{nsamp}{number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}. Default is \code{nsamp = 500}.
    If \code{nsamp="best"} exhaustive enumeration is done, as far as
    the number of trials do not exceed 5000. If \code{nsamp="exact"}
    exhaustive enumeration will be attempted however many samples
    are needed. In this case a warning message will be displayed
    saying that the computation can take a very long time. }
  \item{seed}{initial seed for R's random number generator; see
    \code{\link{.Random.seed}} and the description of the \code{seed}
    argument in \code{\link{lmrob.control}}.}
%% FIXME: This  REALLY should be renamed to  'trace' !!
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{use.correction}{ whether to use finite sample correction factors. Default is \code{use.correction=TRUE}}
  \item{adjust}{whether to perform intercept adjustment at each step. This could be quite
time consuming, therefore the default is \code{adjust = FALSE} }
}
\details{
For details about the estimation options see the corresponding estimation functions.
}
\value{
  A list with components, as the parameters passed  by the invocation

}
\author{Valentin Todorov}
% \references{
% }

\examples{
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])

ctrl <- rrcov.control(alpha=0.75, trace=TRUE)
covMcd(hbk.x,      control = ctrl)
covMcd(log(brain), control = ctrl)
}
\keyword{robust}
\keyword{multivariate}
