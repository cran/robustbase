\name{Qn}
\alias{Qn}
\alias{Qn.old}
\alias{s_Qn}
%
\title{Robust Location-Free Scale Estimate More Efficient than MAD}
\description{
  Compute the robust scale estimator \eqn{Q_n}{Qn}, an efficient
  alternative to the MAD.

  By default, \eqn{Q_n(x_1, \ldots, x_n)}{Qn(x1, .., xn)} is the \eqn{k}-th
  order statistic (a quantile) of the \code{choose(n, 2)} absolute
  differences \eqn{|x_i - x_j|}{abs(x[i] - x[j])},
  (for \eqn{1 \le i < j \le n}{1 <= i < j <= n}),
  where by default (originally only possible value) \eqn{k = choose(n\%/\% 2 + 1, 2)}
  which is about the first quartile (25\% quantile) of these
  pairwise differences.  See the references for more.
}
\usage{
Qn(x, constant = NULL, finite.corr = is.null(constant) && missing(k),
   na.rm = FALSE, k = choose(n \%/\% 2 + 1, 2), warn.finite.corr = TRUE)

s_Qn(x, mu.too = FALSE, \dots)
}
\arguments{% >> ../R/qnsn.R <<
  \item{x}{numeric vector of observations.}
  \item{constant}{number by which the result is multiplied; the default
    achieves consistency for normally distributed data.  Note that until
    Nov. 2010, \dQuote{thanks} to a typo in the very first papers, a slightly
    wrong default constant, 2.2219, was used instead of the correct one
    which is equal to \code{1 / (sqrt(2) * qnorm(5/8))} (as mentioned
    already on p.1277, after (3.7) in Rousseeuw and Croux (1993)).

    If you need the old slightly off version for historical
    reproducibility, you can use \code{Qn.old()}.

    Note that the relative difference is only about 1 in 1000, and that
    the correction should not affect the finite sample corrections for
    \eqn{n \le 9}{n <= 9}.
  }
  \item{finite.corr}{logical indicating if the finite sample bias
    correction factor should be applied.  Defaults to \code{TRUE} unless
    \code{constant} is specified.  Note the for non-default \code{k}, the
    consistency \code{constant} already depends on \code{n} leading to
    \emph{some} finite sample correction, but no simulation-based
    small-\code{n} correction factors are available.}
  \item{na.rm}{logical specifying if missing values (\code{\link{NA}})
    should be removed from \code{x} before further computation.  If false
    as by default, and if there are \code{NA}s, i.e., \code{if(anyNA(x))},
    the result will be \code{NA}.}
  \item{k}{integer, typically half of n, specifying the \dQuote{quantile}, i.e., rather the
    order statistic that \code{Qn()} should return; for the Qn() proper,
    this has been hard wired to \code{choose(n\%/\%2 +1, 2)}, i.e.,
    \eqn{\lfloor\frac{n}{2}\rfloor +1}{floor(n/2) +1}.  Choosing a large \code{k} is less robust but
    allows to get non-zero results in case the default \code{Qn()} is zero.}
  \item{warn.finite.corr}{logical indicating if a \code{\link{warning}}
    should be signalled when \code{k} is non-default, in which case specific
    small-\eqn{n} correction is not yet provided.}
  \item{mu.too}{logical indicating if the \code{\link[stats]{median}(x)} should
    also be returned for \code{s_Qn()}.}
  \item{\dots}{potentially further arguments for \code{s_Qn()} passed to
    \code{Qn()}.}
}
\value{
  \code{Qn()} returns a number, the \eqn{Q_n}{Qn} robust scale
  estimator, scaled to be consistent for \eqn{\sigma^2} and
  i.i.d. Gaussian observations, optionally bias corrected for finite
  samples.

  \code{s_Qn(x, mu.too=TRUE)} returns a length-2 vector with location
  (\eqn{\mu}) and scale; this is typically only useful for
  \code{\link{covOGK}(*, sigmamu = s_Qn)}.
}
\details{
  As the (default, consistency) constant needed to be corrected,
  the finite sample correction has been based on a much more extensive
  simulation, and on a 3rd or 4th degree polynomial model in \eqn{1/n}
  for odd or even n, respectively.
}
\references{
  Rousseeuw, P.J. and Croux, C. (1993)
  Alternatives to the Median Absolute Deviation,
  \emph{Journal of the American Statistical Association} \bold{88}, 1273--1283.
  \doi{10.2307/2291267}% JSTOR
  % MM: ~/save/papers/Rousseeuw/93/R+Croux_MAD_Sn_Qn.pdf

  Christophe Croux and Peter J. Rousseeuw (1992)
  A class of high-breakdown scale estimators based on subranges ,
  \emph{Communications in Statistics - Theory and Methods} \bold{21}, 1935--1951;
  \doi{10.1080/03610929208830889}

  Christophe Croux and Peter J. Rousseeuw (1992)
  Time-Efficient Algorithms for Two Highly Robust Estimators of Scale,
  \emph{Computational Statistics, Vol. 1}, ed. Dodge and Whittaker,
  Physica-Verlag Heidelberg, 411--428; available via Springer Link.
  % MM: ~/save/papers/robust-diverse/Croux-Rousseeuw-Timeff_Scale_1992.pdf
  %% no longer \url{http://win-www.uia.ac.be/u/statis/abstract/Timeff92.htm}.

  About the typo in the \code{constant}:\cr
  Christophe Croux (2010)
  Private e-mail, Fri Jul 16, w/ Subject
  \emph{Re: Slight inaccuracy of Qn implementation \dots\dots}.
}
\seealso{\code{\link[stats]{mad}} for the \sQuote{most robust} but much less efficient
  scale estimator; \code{\link{Sn}} for a similar faster but less
  efficient alternative.  Finally, \code{\link{scaleTau2}} which some
  consider \dQuote{uniformly} better than Qn or competitors.
}
\author{Original Fortran code:
  Christophe Croux and Peter Rousseeuw \email{rousse@wins.uia.ac.be}.
  \cr
  Port to C and R: Martin Maechler, \email{maechler@R-project.org}
}
\examples{
set.seed(153)
x <- sort(c(rnorm(80), rt(20, df = 1)))
s_Qn(x, mu.too = TRUE)
Qn(x, finite.corr = FALSE)

## A simple pure-R version of Qn() -- slow and memory-rich for large n: O(n^2)
Qn0R <- function(x, k = choose(n \%/\% 2 + 1, 2)) { %
    n <- length(x <- sort(x))
    if(n == 0) return(NA) else if(n == 1) return(0.)
    stopifnot(is.numeric(k), k == as.integer(k), 1 <= k, k <= n*(n-1)/2)
    m <- outer(x,x,"-")# abs not needed as x[] is sorted
    sort(m[lower.tri(m)], partial = k)[k]
}
(Qx1 <- Qn(x, constant=1)) # 0.5498463
## the C-algorithm "rounds" to 'float' single precision ..
stopifnot(all.equal(Qx1, Qn0R(x), tol = 1e-6))


(qn <- Qn(c(1:4, 10, Inf, NA), na.rm=TRUE))
stopifnot(is.finite(qn), all.equal(4.075672524, qn, tol=1e-10))

## -- compute for different 'k' :

n <- length(x) # = 100 here
(k0 <- choose(floor(n/2) + 1, 2)) # 51*50/2 == 1275
stopifnot(identical(Qx1, Qn(x, constant=1, k=k0)))
nn2 <- n*(n-1)/2
all.k <- 1:nn2
system.time(Qss <- sapply(all.k, function(k) Qn(x, 1, k=k)))
system.time(Qs  <- Qn  (x, 1, k = all.k))
system.time(Qs0 <- Qn0R(x,    k = all.k) )
stopifnot(exprs = {
   Qs[1]   == min(diff(x))
   Qs[nn2] == diff(range(x))
   all.equal(Qs,  Qss, tol = 1e-15) # even exactly
   all.equal(Qs0, Qs, tol = 1e-7) # see 2.68e-8, as Qn() C-code rounds to (float)
})

plot(2:nn2, Qs[-1], type="b", log="y", main = "Qn(*, k),  k = 2..n(n-1)/2")
}
\keyword{robust}
\keyword{univar}
