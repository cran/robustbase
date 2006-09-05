\name{summarizeRobWeights}
\alias{summarizeRobWeights}
\title{Print a Nice "summary" of Robustness Weights}
\description{
  Print a nice \dQuote{summary} about a numeric vector of robustness
  weights.  Observations with weights around zero are marked as outliers.
}
\usage{
summarizeRobWeights(w, digits = getOption("digits"),
                    header = "Robustness weights:",
                    eps = 1e-04, eps1 = eps, \dots)
}
\arguments{
  \item{w}{numeric vector of robustness weigths.}
  \item{digits}{digits to be used for \code{\link{print}}ing.}
  \item{header}{string to be printed as header line.}
  \item{eps}{numeric tolerance \eqn{\epsilon}{eps}: values of \code{w}
    with \eqn{\left|w_i\right| < \epsilon/n}{|w[i]| < eps/n} are said to
    be outliers.}
  \item{eps1}{numeric tolerance: values of \code{w} with
    \eqn{\left|1 - w_i\right| < eps1}{|1 - w[i]| < eps1} are said to
    have weight \sQuote{\code{~= 1}}.}
  \item{\dots}{potential further arguments, passed to
    \code{\link{print}()}.}
}
\value{
  none; the function is used for its side effect of printing.
}
\author{Martin Maechler}
\seealso{The \code{\link{summary}} methods for \code{\link{lmrob}}
  and \code{\link{glmrob}} make use of \code{summarizeRobWeights()}.
}
\examples{
%% TODO: find an explicit one, particularly for "playing with" eps
}
\keyword{utilities}
