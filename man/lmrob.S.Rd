\name{lmrob.S}
\alias{lmrob.S}
\title{ S-regression estimators }
\description{
  Computes an S-estimator for linear regression,
  using the \dQuote{fast S} algorithm.% of Matias Salibian & Victor Yohai ...
}
\usage{
lmrob.S(x, y, control, trace.lev = 0)
}
\arguments{
  \item{x}{ design matrix }
  \item{y}{ response vector }
  \item{control}{ list as returned by \code{\link{lmrob.control}} }
  \item{trace.lev}{integer indicating if the progress of the algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
}
\details{
  This function is used by \code{\link{lmrob.fit.MM}} and not
  intended to be used on its own (because an S-estimator has too low
  efficiency \sQuote{on its own}).
}
\value{
  A list with components
  \item{coef}{numeric vector (length \eqn{p}) of S-regression coefficient estimates.}
  \item{scale}{the S-scale residual estimate}% 'residual estimate' ?? % resid. VAR !?
  \item{cov}{covariance matrix (\eqn{p \times p}{p x p}) of the
    coefficient estimates.}
  \item{k.iter}{(maximal) number of refinement iterations used.}
  \item{converged}{logical indicating if \bold{all} refinement
    iterations had converged.}
  \item{control}{the same list as the \code{control} argument.}
}
\seealso{\code{\link{lmrob}}, also for references.
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
m.lmS[c("coef","scale")]
stopifnot(all.equal(m.lmS$coef,
                    10 * (1:4),        tol = 0.005),
          all.equal(m.lmS$scale, 1/10, tol = 0.09))
}
\author{ Matias Salibian-Barrera (and Martin Maechler for minor details) }
\keyword{robust}
\keyword{regression}
