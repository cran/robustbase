\name{huberize}
\alias{huberize}
\title{Huberization -- Bringing Outliers In}
\description{
  Huberization (named after Peter Huber's M-estimation algorithm for
  location originally) replaces outlying values in a sample \code{x} by
  their respective boundary: when \eqn{x_j < c_1} it is replaced by \eqn{c_1}
  and when \eqn{x_j > c_2}  it is replaced by \eqn{c_2}.  Consequently,
  values inside the interval \eqn{[c_1, c_2]} remain unchanged.

  Here, \eqn{c_j = M \pm c\cdot s}{c1,c2 = M +/- c*s} where \eqn{s := s(x)} is
  the \emph{robust} scale estimate \code{\link{Qn}(x)} if that is positive,
  and by default, \eqn{M} is the robust huber estimate of location
  \eqn{\mu} (with tuning constant \eqn{k}).

  In the degenerate case where \code{\link{Qn}(x) == 0}, trimmed means of
  \code{abs(x - M)} are tried as scale estimate \eqn{s}, with decreasing
  trimming proportions specified by the decreasing \code{trim} vector.
}
\usage{
huberize(x, M = huberM(x, k = k)$mu, c = k,
         trim = (5:1)/16,
         k = 1.5,
         warn0 = getOption("verbose"), saveTrim = TRUE)
}
\arguments{
  \item{x}{numeric vector which is to be huberized.}
  \item{M}{a number; defaulting to \code{\link{huberM}(x, k)}, the robust
    Huber M-estimator of location.}
  \item{c}{a positive number, the tuning constant for huberization of the
    sample \code{x}.}
  \item{trim}{a \emph{decreasing} vector of trimming proportions in
    \eqn{[0, 0.5]}, only used to trim the absolute deviations from \code{M}
    in case \code{\link{Qn}(x)} is zero.
  }
  \item{k}{used if \code{M} is not specified as huberization center
    \code{M}, and so, by default is taken as Huber's M-estimate
    \code{\link{huberM}(x, k)}.}
  \item{warn0}{\code{\link{logical}} indicating if a warning should be
    signalled in case \code{\link{Qn}(x)} is zero and the trimmed means for
    all trimming proportions \code{trim} are zero as well.}
  \item{saveTrim}{a \code{\link{logical}} indicating if the last tried
    \code{trim[j]} value should be stored if \code{\link{Qn}(x)} was zero.}
}
\details{
  \itemize{
    \item In regular cases, \code{s = \link{Qn}(x)} is positive and used to
    huberize values of \code{x} outside \code{[M - c*s, M + c*s]}.

    \item In degenerate cases where \code{\link{Qn}(x) == 0}, we search for
    an \eqn{s > 0} by trying the trimmed mean \code{s := mean(abs(x-M), trim =
      trim[j])} with less and less trimming (as the trimming
    proportions \code{trim[]} must decrease).
    If even the last, \code{trim[length(trim)]}, leads to \eqn{s = 0}, a
    warning is printed when \code{warn0} is true.
  }
}
\value{
  a numeric vector as \code{x}; in case \code{\link{Qn}(x)} was zero and
  \code{saveTrim} is true, also containing the (last) \code{trim}
  proportion used (to compute the scale \eqn{s}) as attribute \code{"trim"}
  (see \code{\link{attr}()}, \code{\link{attributes}}).
}
\author{Martin Maechler}
\note{
  For the use in \code{\link{mc}()} and similar cases where mainly numerical
  stabilization is necessary, a large \code{c = 1e12} will lead to \emph{no}
  huberization, i.e., all \code{y == x} for \code{y <- huberize(x, c)}
  for typical non-degenerate samples.
}
\seealso{
  \code{\link{huberM}} and \code{\link{mc}} which is now stabilized by
  default via something like \code{huberize(*, c=1e11)}.
}
\examples{
## For non-degenerate data and large c, nothing is huberized,
## as there are *no* really extreme outliers :
set.seed(101)
x <- rnorm(1000)
stopifnot(all.equal(x, huberize(x, c=100)))
## OTOH, the "extremes" are shrunken towards the boundaries for smaller c:
xh <- huberize(x, c = 2)
table(x != xh)
## 45 out of a 1000:
table(xh[x != xh])# 26 on the left boundary -2.098 and 19 on the right = 2.081
## vizualization:
stripchart(x); text(0,1, "x {original}", pos=3); yh <- 0.9
stripchart(xh, at = yh, add=TRUE, col=2)
text(0, yh, "huberize(x, c=2)",   col=2, pos=1)
arrows( x[x!=xh], 1,
       xh[x!=xh], yh, length=1/8, col=adjustcolor("pink", 1/2))
}
\keyword{robust}
\keyword{univar}
