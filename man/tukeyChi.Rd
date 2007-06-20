\name{tukeyChi}
\alias{tukeyChi}
\title{Tukey's "Chi", the Bi-square Loss (Rho) Function}
\description{
  Computes Tukey's bi-square loss function, \code{chi(x)} and its
  first two derivatives.  Note that in the general context of
  \eqn{M}-estimators, these loss functions are called
  \eqn{\rho (rho)}{rho}-functions.
}
\usage{
tukeyChi(x, cc, deriv = 0)
}
\arguments{
  \item{x}{numeric vector.}
  \item{cc}{ tuning constant }
  \item{deriv}{integer in \eqn{\{0,1,2\}} specifying the order of the
    derivative; the default, \code{deriv = 0} computes the chi- (or
    rho-)function.}
}
\value{
  a numeric vector of the same length as \code{x}.
}
\note{\code{tukeyChi(x, d)} and \code{\link{tukeyPsi1}(x, d-1)} are just
  re-scaled versions of each other (for \code{d in 0:2}).
}
\author{Matias Salibian-Barrera and Martin Maechler}
\seealso{
  \code{\link{lmrob}} and \code{\link{tukeyPsi1}}.
}
\examples{
op <- par(mfrow = c(3,1), oma = c(0,0, 2, 0),
          mgp = c(1.5, 0.6, 0), mar= .1+c(3,4,3,2))
x <- seq(-2.5, 2.5, length = 201)
cc <- 1.55 # as set by default in lmrob.control()
plot. <- function(...) { plot(...); abline(h=0,v=0, col="gray", lty=3)}
plot.(x, tukeyChi(x, cc), type = "l", col = 2)
plot.(x, tukeyChi(x, cc, deriv = 1), type = "l", col = 2)
plot.(x, tukeyChi(x, cc, deriv = 2), type = "l", col = 2)
%                               \ is escape for Rd
mtext(sprintf("tukeyChi(x, c = \%g, deriv),  deriv = 0,1,2", cc),
      outer = TRUE, font = par("font.main"), cex = par("cex.main"))
par(op)
}
\keyword{robust}
