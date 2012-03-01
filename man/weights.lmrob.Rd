\name{weights.lmrob}
\alias{weights.lmrob}
\title{
  Extract Robustness and Model Weights
}
\description{
  \code{weights} which extracts robustness weights or fitting
  (or prior) weights from a \code{lmrob} object.
}
\usage{
\method{weights}{lmrob}(object, type = c("prior", "robustness"), ...)
}
\arguments{
  \item{object}{
    an object of class \code{lmrob}, typically the result of
    a call to \code{\link{lmrob}}.
}
  \item{type}{
    the type of weights which should be returned.
    The alternatives are: \code{"prior"} (default), and
    \code{"robustness"}.
}
  \item{\dots}{
    not used currently.
}
}
\details{
  The \dQuote{prior weights} correspond to the weights specified using
  the \dQuote{weights} argument when calling \code{lmrob}. The
  \dQuote{robustness weights} are the weights assigned by the
  M-estimator of regression, \eqn{\psi(r_i/S) / (r_i/S)}. The robust
  coefficient estimate then numerically corresponds to a weighted least
  squares fit using the product of both types of weights as weights.
}
\value{
  Weights extracted from the object \code{object}.
}
\author{
  Manuel Koller
}

\seealso{
  \code{\link{lmrob}} and \code{\link{weights}}
}
