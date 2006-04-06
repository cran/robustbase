\name{summary.lmrob}
\alias{summary.lmrob}
\alias{print.summary.lmrob}
\title{Summary Method for "lmrob" Objects}
\description{
  Summary method for \R object of class \code{"lmrob"} and
  \code{\link{print}} method for the summary object.
}
\usage{
\method{summary}{lmrob}(object, correlation = FALSE,
        symbolic.cor = FALSE, \dots)
\method{print}{summary.lmrob}(x, digits = max(3, getOption("digits") - 3),
      symbolic.cor= x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), \dots)
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
  \item{\dots}{potentially more arguments passed to methods.}
}
\seealso{\code{\link{lmrob}}, \code{\link{summary.lm}},
    \code{\link{print}}, \code{\link{summary}}.
}
\examples{
mod1 <- lmrob(stack.loss ~ ., data = stackloss)
sa <- summary(mod1)  # calls summary.lmrob(....)
sa                   # dispatches to call print.summary.lmrob(....)

data(heart)
sm2 <- summary( m2 <- lmrob(clength ~ ., data = heart) )
sm2
}
\keyword{robust}
\keyword{regression}
