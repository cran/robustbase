\name{lmrob.fit.MM}
\alias{lmrob.fit.MM}
\title{ MM-estimator for regression }
\description{
  Compute MM-estimators of regression:  An S-estimator is
  used as starting value, and an M-estimator with fixed scale and
  redescending psi-function is used from there.
}
\usage{
lmrob.fit.MM(x, y, control)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{control}{A list of control parameters as returned
    by \code{\link{lmrob.control}}, used for both the initial S-estimate
    and the subsequent M-estimate.}
}
\details{This function is the basic fitting function for MM-estimation,
  called by \code{\link{lmrob}} and typically not to be used on its own.

  It calls \code{\link{lmrob.S}(..)} and uses it as initial estimator.
  Note that the inference used (covariance matrix) depends crucially on
  the S-estimator used, and hence it is currently no longer possible to
  specify the S-estimator at this level.
}
\value{
  A list with components
  \item{fitted.values}{\eqn{X \beta}{X beta}, i.e. \code{X \%*\% coefficients}.}
  \item{residuals}{the raw residuals, \code{y - fitted.values}}
  \item{weights}{robustness weights derived from the final M-estimator
    residuals (even when not converged).}
  \item{rank}{}
  \item{degree.freedom}{\code{n - rank}}% more!
  \item{coefficients}{estimated regression coefficient vector}
  \item{initial.coefficients}{}% = iCoef,
  \item{scale}{the robustly estimated error standard deviation}% = final.MM$scale,
  \item{cov}{variance-covariance matrix of \code{coefficients}, if the
    RWLS iterations have converged, otherwise the vcov-matrix of the
    initial estimator.}
  \item{control}{}% = control,
  \item{iter}{}% = final.MM$iter, <<<<<< also 'init.S' !
  \item{converged}{logical indicating if the RWLS iterations have converged.}
  \item{init.S}{the whole initial S-estimator result, including its own
    \code{converged} flag, see \code{\link{lmrob.S}}.}
}
\seealso{
  \code{\link{lmrob}},
  \code{\link{lmrob..M..fit}},
  \code{\link{lmrob.S}}
}
\keyword{robust}
\keyword{regression}

