\name{scaleTau2}
\alias{scaleTau2}
\title{Robust Tau-Estimate of Scale}
\description{
  Computes the robust \eqn{\tau}-estimate of univariate scale, as
  proposed by Maronna and Zamar (2002); improved by a consistency factor,
  %% FIXME: TODO: add a finite sample correction by Martin Maechler
  %% (currently have 'n-2' but can even be better !!!!
}
\usage{%--> ../R/OGK.R <<<
scaleTau2(x, c1 = 4.5, c2 = 3.0, na.rm = FALSE, consistency = TRUE,
          mu0 = median(x),
          sigma0 = median(x.), mu.too = FALSE, iter = 1, tol.iter = 1e-7)
}
\arguments{
  \item{x}{numeric vector}
  \item{c1,c2}{non-negative numbers, specifying cutoff values for the
    biweighting of the mean and the rho function respectively.}
  \item{na.rm}{a logical value indicating whether \code{NA}
    values should be stripped before the computation proceeds.}
  \item{consistency}{logical indicating if the consistency correction
    factor (for the scale) should be applied.}
  \item{mu0}{the initial location estimate \eqn{\mu_0}{mu0}, defaulting to
    the \code{\link{median}}.}
  \item{sigma0}{the initial scale estimate \eqn{s_0}{s0}, defaulting to
    the MAD; may be set to a positive value when the MAD is zero.}
  \item{mu.too}{logical indicating if both location and scale should be
    returned or just the scale (when \code{mu.too=FALSE} as by default).}
  \item{iter}{positive integer or logical indicating if and how many
    iterations should be done.  The default, \code{iter = 1} computes the
    \dQuote{traditional} tau-estimate of scale.}
  \item{tol.iter}{if \code{iter} is true, or \code{iter > 1}, stop the
    iterations when \eqn{|s_n - s_o| \le \epsilon s_n}, where
    \eqn{\epsilon :=}\code{tol.iter}, and  \eqn{s_o} and \eqn{s_n} are the
    previous and current estimates of \eqn{\sigma}.}
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
  When \code{iter=TRUE} or \code{iter > 1}, the above estimate is
  \emph{iterated} in a fixpoint iteration, setting \eqn{s_0} to the current
  estimate \eqn{s(X)} and iterating until the number of iterations is
  larger than \code{iter} or the fixpoint is found in the sense that
  \
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
  % MM: ~/save/papers/robust-diverse/Maronna-Zamar-OGK_2002.pdf

  Yohai, V.J., and Zamar, R.H. (1988).
  High breakdown-point estimates of regression by means of the
  minimization of an efficient scale.
  \emph{Journal of the American Statistical Association} \bold{83}, 406--413.
  % MM: ~/save/papers/robust-diverse/Yohai-Zamar-tau_JASA1988.pdf
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
(sI  <- scaleTau2(c(x,Inf), mu.too = TRUE))
(sIN <- scaleTau2(c(x,Inf,NA), mu.too = TRUE, na.rm=TRUE))
stopifnot({
  identical(sI, sIN)
  all.equal(scaleTau2(c(x, 999), mu.too = TRUE), sIN,
            tol = 1e-15)
})

if(doExtras <- robustbase:::doExtras()) {
 set.seed(11)
 ## show how much faster this is, compared to Qn
 x <- sample(c(rnorm(1e6), rt(5e5, df=3)))
 (system.time(Qx <- Qn(x)))         ## 2.04 [2017-09, lynne]
 (system.time(S2x <- scaleTau2(x))) ## 0.25    (ditto)
 cbind(Qn = Qx, sTau2 = S2x)
}##       Qn    sTau2
##  1.072556 1.071258
}
\keyword{robust}
\keyword{univar}
