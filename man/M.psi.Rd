\name{Mpsi}
\alias{Mchi}
\alias{Mpsi}
\alias{Mwgt}
\alias{MrhoInf}
\alias{.Mchi}
\alias{.Mpsi}
\alias{.Mwgt}
\alias{.MrhoInf}
\title{Psi / Chi / Wgt / Rho Functions for *M-Estimation}
\description{
  Compute Psi / Chi / Wgt / Rho functions for M-estimation,
  i.e., including MM, etc.
  %% TODO: More, notably definitions ... but they are all nicely in the
  %% vignette.. How can we link from here to there ???

  \code{MrhoInf(x)} computes \eqn{\rho(\infty)}{rho(Inf)}, i.e., the
  normalizing or scaling constant for the transformation
  from \eqn{\rho(\cdot)}{rho(.)} to
  \eqn{\tilde\rho(\cdot)}{rho~(.)}, where the latter, aka as
  \eqn{\chi()}{chi()} fulfills \eqn{\tilde\rho(\infty) = 1}{rho~(Inf) = 1}
  which makes only sense for \dQuote{redescending} psi functions, i.e.,
  not for \code{"huber"}.

  \code{Mwgt(x, *)} computes \eqn{\psi(x)/x}  (fast and numerically accurately).
}
\usage{
Mpsi(x, cc, psi, deriv = 0)
Mchi(x, cc, psi, deriv = 0)
Mwgt(x, cc, psi)
MrhoInf(cc, psi)
}
\arguments{
  \item{x}{numeric (\dQuote{abscissa} values) vector.}
  \item{cc}{numeric tuning constant, for some \code{psi} of length
    \eqn{> 1}.}
  \item{psi}{a string specifying the psi / chi / rho / wgt function;
    either \code{"huber"}, or one of the same possible specifiers as for
    \code{psi} in \code{\link{lmrob.control}}, i.e. currently,
    \code{"bisquare"}, \code{"lqq"}, \code{"welsh"}, \code{"optimal"},
    \code{"hampel"}, or \code{"ggw"}.}
  \item{deriv}{an integer, specifying the \emph{order} of derivative to
    consider; particularly, \code{Mpsi(x, *, deriv = -1)} is the
    principal function of \eqn{\psi()}{psi()}, typically denoted
    \eqn{\rho()}{rho()} in the literature.}
  %% FIXME: mention that  deriv = 2 is *partially* implemented
}
\details{
  Theoretically, \code{Mchi()} would not be needed explicitly as it can be computed
  from \code{Mpsi()} and \code{MrhoInf()}, namely, by
  \preformatted{Mchi(x, *, deriv = d)  ==  Mpsi(x, *, deriv = d-1) / MrhoInf(*)}
  for \eqn{d = 0, 1, 2}  (and \sQuote{*} containing \code{par, psi}, and
  equality is in the sense of \code{\link{all.equal}(x,y, tol)} with a
  small \code{tol}.

  Similarly, \code{Mwgt} would not be needed strictly, as it could be
  defined via \code{Mpsi}), but the explicit definition takes care of
  0/0 and typically is of a more simple form.

  For experts, there are slightly even faster versions,
  \code{.Mpsi()}, \code{.Mwgt()}, etc.
}
\value{
  a numeric vector of the same length as \code{x}, with corresponding
  function (or derivative) values.
}
\references{
  See the vignette about %% \link{ .. } ????
  \dQuote{\eqn{\psi}{psi}-Functions Available in Robustbase}.
  %% ../inst/doc/psi_functions.Rnw
}
\author{
  Manuel Koller, notably for the original C implementation;
  tweaks and speedup via \code{\link{.Call}} by Martin Maechler.
}
\seealso{
 \code{\link{psiFunc}} and the \code{\linkS4class{psi_func}} class, both
 of which provide considerably more on the \R side, but are less
 optimized for speed.
}
\examples{
x <- seq(-5,7, by=1/8)
matplot(x, cbind(Mpsi(x, 4, "biweight"),
                 Mchi(x, 4, "biweight"),
                 Mwgt(x, 4, "biweight")), type = "l")
abline(h=0, v=0, lty=2, col=adjustcolor("gray", 0.6))

hampelPsi
(ccHa <- hampelPsi @ xtras $ tuningP $ k)
psHa <- hampelPsi@psi(x)
% FIXME: interesting as long as hampelPsi does not use Mpsi(... "hampel") !
## using Mpsi():
Mp.Ha <- Mpsi(x, cc = ccHa, psi = "hampel")
stopifnot(all.equal(Mp.Ha, psHa, tol = 1e-15))
}
\keyword{robust}
