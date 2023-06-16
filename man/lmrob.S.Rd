\name{lmrob.S}
\alias{lmrob.S}
\title{ S-regression estimators }
\description{
  Computes an S-estimator for linear regression,
  using the \dQuote{fast S} algorithm.% of Matias Salibian & Victor Yohai ...
}
\usage{
lmrob.S(x, y, control,
        trace.lev = control$trace.lev,
        only.scale = FALSE, mf)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n * p})}
  \item{y}{numeric vector of responses (or residuals for \code{only.scale=TRUE}).}
  \item{control}{ list as returned by \code{\link{lmrob.control}}; the
    following components are used for \code{lmrob.S()}:
    \code{"trace.lev"}, % as default for direct arg {trace.lev}, see below
    \code{"nResample"},
    \code{"groups"},
    \code{"n.group"},
    \code{"fast.s.large.n"},
    \code{"seed"},
    \code{"bb"},
    \code{"psi"}, \code{"tuning.chi"},% c.chi
    \code{"best.r.s"},
    \code{"k.fast.s"},
    \code{"k.max"},
    \code{"maxit.scale"},
    \code{"refine.tol"}, \code{"solve.tol"}, \code{"scale.tol"},
    \code{"mts"},
    \code{"subsampling"}.
  }
  \item{trace.lev}{integer indicating if the progress of the algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
  \item{only.scale}{\code{\link{logical}} indicating if only the scale
    of \code{y} should be computed.  In this case, \code{y} will
    typically contain \emph{residuals}.}%% FIXME: explain
  %% namely, s, fulfilling the \sum_i chi(i) = 1/2  equation.
  \item{mf}{defunct.}
}
\details{
  This function is used by \code{\link{lmrob.fit}} and typically not
  to be used on its own (because an S-estimator has too low
  efficiency \sQuote{on its own}).

  By default, the subsampling algorithm uses a customized LU
  decomposition which ensures a non singular subsample (if this is at
  all possible). This makes the Fast-S algorithm also feasible for
  categorical and mixed continuous-categorical data.

  One can revert to the old subsampling scheme by setting the parameter
  \code{subsampling} in \code{control} to \code{"simple"}.
}
\value{
  By default (when \code{only.scale} is false), a list with components
  \item{coefficients}{numeric vector (length \eqn{p}) of S-regression coefficient estimates.}
  \item{scale}{the S-scale residual estimate}% 'residual estimate' ?? % resid. VAR !?
  % \item{cov}{covariance matrix (\eqn{p \times p}{p x p}) of the
  %   coefficient estimates.}
  \item{fitted.values}{numeric vector (length \eqn{n}) of the fitted
    values.}
  \item{residuals}{numeric vector (length \eqn{n}) of the residuals.}
  \item{rweights}{numeric vector (length \eqn{n}) of the robustness weights.}
  \item{k.iter}{(maximal) number of refinement iterations used.}
  \item{converged}{logical indicating if \bold{all} refinement
    iterations had converged.}
  \item{control}{the same list as the \code{control} argument.}

  If \code{only.scale} is true, the computed scale (a number) is returned.
}
\seealso{\code{\link{lmrob}}, also for references.
}
\author{ Matias Salibian-Barrera and Manuel Koller; Martin Maechler for
  minor new options and more documentation.
}
\examples{
set.seed(33)
x1 <- sort(rnorm(30)); x2 <- sort(rnorm(30)); x3 <- sort(rnorm(30))
X. <- cbind(x1, x2, x3)
y <-  10 + X. \%*\% (10*(2:4)) + rnorm(30)/10
y[1] <- 500   # a moderate outlier
X.[2,1] <- 20 # an X outlier
X1  <- cbind(1, X.)

(m.lm <- lm(y ~ X.))
set.seed(12)
m.lmS <- lmrob.S(x=X1, y=y,
                 control = lmrob.control(nRes = 20), trace.lev=1)
m.lmS[c("coefficients","scale")]
all.equal(unname(m.lmS$coef), 10 * (1:4), tolerance = 0.005)
stopifnot(all.equal(unname(m.lmS$coef), 10 * (1:4), tolerance = 0.005),
          all.equal(m.lmS$scale, 1/10, tolerance = 0.09))

## only.scale = TRUE:  Compute the S scale, given residuals;
s.lmS <- lmrob.S(X1, y=residuals(m.lmS), only.scale = TRUE,
                 control = lmrob.control(trace.lev = 3))
all.equal(s.lmS, m.lmS$scale) # close: 1.89e-6 [64b Lnx]
}
\keyword{robust}
\keyword{regression}
