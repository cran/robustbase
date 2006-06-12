\name{lmrob.MM}
\alias{lmrob.MM}
\title{ Function to compute MM-regression estimators }
\description{
  This function performs RWLS iterations to find an
  MM-regression estimator
}
\usage{
lmrob.MM(x, y, beta.initial, scale, control)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{beta.initial}{numeric vector (of length \eqn{p}) of initial
    estimate.  Usually the result of an S-regression estimator.}
  \item{scale}{ A robust residual scale estimate. Usually
    an S-scale estimator.}
  \item{control}{A list of control parameters as returned
    by \code{\link{lmrob.control}}.}
}
\details{
  This function is used by \code{\link{lmrob.fit.MM}} and not intended
  to be used on its own.
}
\value{A list with the following elements:
  \item{coef }{ The MM-regression estimator }
  \item{cov  }{ The covariance matrix of the MM-regression estimator}
  \item{control}{ The control list used}
  \item{scale }{ The residual scale estimate}
  \item{seed  }{ The random number generator seed}
  \item{converged}{ \code{TRUE} if the RWLS iterations converged,
  \code{FALSE} otherwise}
}
\references{
  Yohai, 1987
}
\seealso{
  \code{\link{lmrob.fit.MM}}, \code{\link{lmrob}}
}
\author{ Matias Salibian-Barrera }
\keyword{robust}
\keyword{regression}
