\name{scaleTau2}
\alias{scaleTau2}
\title{Robust Tau-Estimate of Scale}
\description{
  Computes the robust \eqn{\tau}-estimate of univariate scale, as
  proposed by Maronna and Zamar (2002); improved by a consistency factor.
%% FIXME: and a finite sample correction by Martin Maechler !!
}
\usage{
scaleTau2(x, c1 = 4.5, c2 = 3.0, consistency = TRUE,
          mu.too = FALSE, ...)
}
\arguments{
  \item{x}{numeric vector}
  \item{c1,c2}{non-negative numbers, specifying cutoff values for the
    biweighting of the mean and the rho function respectively.}
  \item{mu.too}{logical indicating if both location and scale should be
    returned or just the scale (when \code{mu.too=FALSE} as by default).}
  \item{consistency}{logical indicating if the consistency correction
    factor (for the scale) should be applied.}
  \item{\dots}{potentially additional arguments which are not used.}
}
\details{
  First, \eqn{s_0}{s0} := MAD, i.e. the equivalent of \code{\link{mad}(x,
    constant=1)} is computed.  Robustness weights
  \eqn{w_i := w_{c1}((x_i - med(X))/ s_0)} are computed, where
  \eqn{w_c(u) = max(0, (1 - (u/c)^2)^2)}.  The robust location
  measure is defined as \eqn{\mu(X) := (\sum_i w_i x_i)/(\sum_i w_i)},
  and the robust \eqn{\tau (tau)}{tau}-estimate is \eqn{s(X)^2 :=
    s_0^2 * (1/n) \sum_i \rho_{c2}((x_i - \mu(X))/s_0)},
  where \eqn{\rho_c(u) = min(c^2, u^2)}.
  \cr
  \code{scaleTau2(*, consistency=FALSE)} returns \eqn{s(X)}, whereas
  this value is divided by its asymptotic limit when \code{consistency =
    TRUE} as by default.

  Note that for \code{n = length(x) == 2}, all equivariant scale estimates are
  proportional, and specifically, \code{scaleTau2(x, consistency=FALSE)
    == mad(x, constant=1)}.  See also the reference.
}
\value{
  numeric vector of length one (if \code{mu.too} is \code{FALSE} as by
  default) or two (when \code{mu.too = TRUE}) with robust scale or
  (location,scale) estimators
  \eqn{\hat\sigma(x)}{s^(x)} or
  \eqn{(\hat\mu(x),\hat\sigma(x))}{(m^(x), s^(x))}.
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.
}
\author{Original by Kjell Konis with substantial modifications by
  Martin Maechler.
}
\seealso{\code{\link{Sn}}, \code{\link{Qn}}, \code{\link{mad}};
  further \code{\link{covOGK}} for which \code{scaleTau2} was designed.
}
\examples{
x <- c(1:7, 1000)
sd(x) # non-robust std.deviation
scaleTau2(x)
scaleTau2(x, mu.too = TRUE)
}
\keyword{robust}
\keyword{univar}
