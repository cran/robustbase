\name{lmrob.Psi}
\alias{lmrob.Psi}
\title{Tukey's Bi-square Score (Psi) Function and Derivative}
\description{
  Compute Tukey's bi-square score (psi) function and its first
  derivative.
}
\usage{
lmrob.Psi(x, cc, deriv = 0)
}
\arguments{
  \item{x}{numeric vector.}
  \item{cc}{ tuning constant }
  \item{deriv}{integer in \eqn{\{0,1,2\}} specifying the order of the
    derivative; the default, \code{deriv = 0} computes the chi-function.}
}
\value{
  a numeric vector of the same length as \code{x}.
}
\seealso{
  \code{\link{lmrob}} and \code{\link{lmrob.Chi}}.
}
\author{Matias Salibian-Barrera and Martin Maechler}
\examples{
op <- par(mfrow = 2:1, oma = c(0,0, 2, 0),
          mgp = c(1.5, 0.6, 0), mar= .1+c(3,4,3,2))
x <- seq(-5, 5, length = 201)
cc <- 4.69 # as set by default in lmrob.control()
plot. <- function(...) { plot(...); abline(h=0,v=0, col="gray", lty=3)}
plot.(x, lmrob.Psi(x, cc), type = "l", col = 2)
plot.(x, lmrob.Psi(x, cc, deriv = 1), type = "l", col = 2)
%                               \ is escape for Rd
mtext(sprintf("lmrob.Psi(x, c = \%g, deriv),  deriv = 0,1", cc),
      outer = TRUE, font = par("font.main"), cex = par("cex.main"))
par(op)
}
\keyword{robust}
