\name{Qn}
\alias{Qn}
\title{Robust Location-Free Scale Estimate More Efficient than MAD}
\description{
  Compute the robust scale estimator \eqn{Q_n}{Qn}, an efficient
  alternative to the MAD.
}
\usage{
Qn(x, constant = 2.2219, finite.corr = missing(constant))
}
\arguments{
  \item{x}{numeric vector of observations.}
  \item{constant}{number by which the result is multiplied; the default
    achieves consisteny for normally distributed data.}
  \item{finite.corr}{logical indicating if the finite sample bias
    correction factor should be applied.  Default to \code{TRUE} unless
    \code{constant} is specified.}
}
\value{
  a number, the \eqn{Q_n}{Qn} robust scale estimator, scaled to be
  consistent for \eqn{\sigma^2} and i.i.d. Gaussian observatsions,
  optionally bias corrected for finite samples.
}
\details{
  ............  FIXME ........
}
\references{
  Rousseeuw, P.J. and Croux, C. (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{Journal of the American Statistical Association} \bold{88}, 1273--1283.

  Christophe Croux and Peter J. Rousseeuw (1992)
  Time-Efficient Algorithms for Two Highly Robust Estimators of Scale,
  \emph{Computational Statistics, Vol.\ 1}, ed. Dodge and Whittaker,
  Physica-Verlag Heidelberg, 411--428;\cr
  also available from
  \url{http://win-www.uia.ac.be/u/statis/abstract/Timeff92.htm}.
}
\seealso{\code{\link{mad}} for the `most robust' but much less efficient
  scale estimator; \code{\link{Sn}} for another alternative.
}
\author{Original Fortran code:
  Christophe Croux and Peter Rousseeuw \email{rousse@wins.uia.ac.be}.
  \cr
  Port to C and R: Martin Maechler, \email{maechler@R-project.org}
}
\examples{
set.seed(153)
x <- sort(c(rnorm(80), rt(20, df = 1)))
Qn(x)
Qn(x, finite.corr = FALSE)
}
\keyword{robust}
\keyword{univar}
