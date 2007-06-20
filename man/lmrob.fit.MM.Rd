\name{lmrob.fit.MM}
\alias{lmrob.fit.MM}
\title{ MM-regression estimators  }
\description{
  Fit function for MM-regression estimators
}
\usage{
lmrob.fit.MM(x, y, control,
             init.S = lmrob.S(x = x, y = y, control = control))
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{control}{A list of control parameters as returned
    by \code{\link{lmrob.control}}.}
  \item{init.S}{a \code{\link{list}} with components \code{coef},
    \code{scale}, and \code{cov} specifying an initial S-estimate, such
    as resulting from \code{\link{lmrob.S}(..)}.}
  %% FIXME: used in practice, should also contain convergence info
  %% ---- *AND* that should be propagated lmrob.S() --> lmrob.fit.MM() --> lmrob()
}
\details{This function is used by \code{\link{lmrob}} and typically not
  to be used on its own.
}
\value{
  A list with components
  \item{fitted.values}{}
  \item{residuals}{}
  \item{weights}{}
  \item{rank}{}
  \item{degree.freedom}{\code{n - rank}}% more!
  \item{coefficients}{regression coefficient estimators}
  \item{initial.coefficients}{}% = iCoef,
  \item{scale}{the robustly estimated error standard deviation}% = final.MM$scale,
  \item{cov}{}% = cov.matrix,
  \item{control}{}% = control,
  \item{iter}{}% = final.MM$iter, <<<<<< also 'init.S' !
  \item{converged}{}% = final.MM$converged) <<<< also 'init.S'
}
\seealso{
  \code{\link{lmrob}},
  \code{\link{lmrob.MM}},
  \code{\link{lmrob.S}}
}
\keyword{robust}
\keyword{regression}

