\name{tukeyPsi1}
\alias{tukeyPsi1}
\title{Tukey's Bi-square Score (Psi) Function and Derivative}
\description{
  Compute Tukey's bi-square score (psi) function, its first
  derivative or it's integral/\dQuote{principal function}.  This is
  scaled such that \eqn{\psi'(0) = 1}, i.e.,
  \eqn{\psi(x) \approx x}{psi(x) ~= x} around 0.
}
\usage{
tukeyPsi1(x, cc, deriv = 0)
}
\arguments{
  \item{x}{numeric vector.}
  \item{cc}{ tuning constant }
  \item{deriv}{integer in \eqn{\{-1,0,1\}} specifying the order of the
    derivative; the default, \code{deriv = 0} computes the psi-function.}
}
\value{
  a numeric vector of the same length as \code{x}.
}
\note{\code{tukeyPsi1(x, d)} and \code{\link{tukeyChi}(x, d+1)} are just
  re-scaled versions of each other (for \code{d in -1:1}), i.e.,
  \deqn{\chi^{(\nu)}(x, c) = (6/c^2) \psi^{(\nu-1)}(x,c),}
  for \eqn{\nu = 0,1,2}.

  We use the name \sQuote{tukeyPsi\bold{1}}, because \code{tukeyPsi} is
  reserved for a future \dQuote{Psi Function} class object, see
  \code{\link{psiFunc}}.
}
\seealso{
  \code{\link{lmrob}} and \code{\link{tukeyChi}}; further
  \code{\link{anova.lmrob}} which needs the \code{deriv = -1}.
}
\author{Matias Salibian-Barrera, Martin Maechler and Andreas Ruckstuhl}
\examples{
op <- par(mfrow = c(3,1), oma = c(0,0, 2, 0),
          mgp = c(1.5, 0.6, 0), mar= .1+c(3,4,1,1))
x <- seq(-5, 5, length = 201)
cc <- 4.69 # as set by default in lmrob.control()
plot. <- function(...) { plot(..., asp = 1); abline(h=0,v=0, col="gray", lty=3)}
plot.(x, tukeyPsi1(x, cc), type = "l", col = 2)
abline(0:1, lty = 3, col = "light blue")
plot.(x, tukeyPsi1(x, cc, deriv = -1), type = "l", col = 2)
plot.(x, tukeyPsi1(x, cc, deriv =  1), type = "l", col = 2); abline(h=1,lty=3)
%                               \ is escape for Rd
mtext(sprintf("tukeyPsi1(x, c = \%g, deriv),  deriv = 0, -1, 1", cc),
      outer = TRUE, font = par("font.main"), cex = par("cex.main"))
par(op)
}
\keyword{robust}
