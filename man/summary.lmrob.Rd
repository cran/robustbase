\name{summary.lmrob}
\title{Summary Method for "lmrob" Objects}
%
\alias{summary.lmrob}
\alias{vcov.lmrob}
\alias{print.summary.lmrob}
\alias{model.matrix.lmrob}
%
\description{
  Summary method for \R object of class \code{"lmrob"} and
  \code{\link{print}} method for the summary object.

  Further, methods \code{\link{fitted}()}, \code{\link{residuals}()} or
  \code{\link{weights}()} work (via the default methods), and
  \code{\link{predict}()} (see \code{\link{predict.lmrob}},
  \code{\link{vcov}()}, \code{\link{model.matrix}()} have explicitly
  defined \code{lmrob} methods.
}
\usage{
\method{summary}{lmrob}(object, correlation = FALSE,
        symbolic.cor = FALSE, \dots)
\method{print}{summary.lmrob}(x, digits = max(3, getOption("digits") - 3),
      symbolic.cor= x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), \dots)

\method{vcov}{lmrob}(object, cov = object$control$cov, \dots)
\method{model.matrix}{lmrob}(object, \dots)
}
\arguments{
  \item{object}{an \R object of class \code{lmrob}, typically created by
    \code{\link{lmrob}}.}
  \item{correlation}{logical variable indicating whether
    to compute the correlation matrix of the estimated coefficients.}
  \item{symbolic.cor}{logical indicating whether
    to use symbols to display the above correlation matrix.}

  \item{x}{an \R object of class \code{summary.lmrob}, typically
    resulting from \code{summary(\link{lmrob}(..),..)}.}
  \item{digits}{number of digits for printing, see \code{digits} in
    \code{\link{options}}.}
  \item{signif.stars}{logical variable indicating
    whether to use stars to display different levels of
    significance in the individual t-tests.}
  \item{cov}{covariance estimation function to use.}
  \item{\dots}{potentially more arguments passed to methods.}
}
\seealso{\code{\link{lmrob}}, \code{\link{predict.lmrob}},
  \code{\link{summary.lm}},
  \code{\link{print}}, \code{\link{summary}}.
}
\examples{
mod1 <- lmrob(stack.loss ~ ., data = stackloss)
sa <- summary(mod1)  # calls summary.lmrob(....)
sa                   # dispatches to call print.summary.lmrob(....)

## correlation between estimated coefficients:
cov2cor(vcov(mod1))

cbind(fit = fitted(mod1), resid = residuals(mod1),
      wgts= weights(mod1),
      predict(mod1, interval="prediction"))

data(heart)
sm2 <- summary( m2 <- lmrob(clength ~ ., data = heart) )
sm2
}
\keyword{robust}
\keyword{regression}
