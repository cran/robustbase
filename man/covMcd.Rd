\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\name{covMcd}
\alias{covMcd}
\alias{print.mcd}
%
\title{Robust Location and Scatter Estimation via MCD}
\concept{High breakdown point}
\description{
     Compute a robust multivariate location and scale estimate with a high
     breakdown point, using the \sQuote{Fast MCD} (Minimum Covariance
     Determinant) estimator.
}
\usage{
covMcd(x, cor = FALSE, raw.only = FALSE,
       alpha =, nsamp =, nmini =, seed =, tolSolve=, trace =,
       use.correction =, wgtFUN =, control = rrcov.control())
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{cor}{should the returned result include a correlation matrix?
    Default is \code{cor = FALSE}}.
  \item{raw.only}{should only the \dQuote{raw} estimate be returned,
    i.e., no weighted reweit the returned result include a correlation matrix?
    Default is \code{cor = FALSE}}.
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{nsamp}{number of subsets used for initial estimates or \code{"best"}
    or \code{"exact"}.  Default is \code{nsamp = 500}.  For
    \code{nsamp = "best"} exhaustive enumeration is done, as long as the
    number of trials does not exceed 100'000 (\code{= nLarge}).  For \code{"exact"},
    exhaustive enumeration will be attempted however many samples are
    needed.  In this case a warning message may be displayed saying
    that the computation can take a very long time.}
  \item{nmini}{for large \eqn{n}, the algorithm splits the data into
    maximally \eqn{k_{\mathrm{rep}} = 5}{krep = 5} subsets of size
    \code{nmini}. %% FIXME? -- more accurately and in detail ?!
    The original algorithm had \code{nmini = 300} hard coded.}
  \item{seed}{initial seed for random generator, see \code{\link{rrcov.control}}.}
  \item{tolSolve}{numeric tolerance to be used for inversion
    (\code{\link{solve}}) of the covariance matrix in \code{\link{mahalanobis}}.}
  \item{trace}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{FALSE}; values \eqn{\ge 2}{>= 2}
    also produce print from the internal (Fortran) code.}
  \item{use.correction}{ whether to use finite sample correction
    factors; defaults to \code{TRUE}.}
  \item{wgtFUN}{a character string or \code{\link{function}}, specifying
    how the weights for the reweighting step should be computed.
    Up to April 2013, the only option has been the original proposal in
    (1999), now specified by \code{wgtFUN = "01.original"} (or via \code{control}).}
  \item{control}{a list with estimation options - this includes those
    above provided in the function specification, see
    \code{\link{rrcov.control}} for the defaults.  If \code{control} is
    supplied, the parameters from it will be used.  If parameters are
    passed also in the invocation statement, they will override the
    corresponding elements of the control object.}
}
\details{
  The minimum covariance determinant estimator of location and scatter
  implemented in \code{covMcd()} is similar to \R function
  \code{\link[MASS]{cov.mcd}()} in \pkg{MASS}.  The MCD method looks for
  the \eqn{h (> n/2)} (\eqn{h = h(\alpha,n,p) =}
  \code{\link{h.alpha.n}(alpha,n,p)}) observations (out of \eqn{n})
  whose classical covariance matrix has the lowest possible determinant.

  The raw MCD estimate of location is then the average of these \eqn{h} points,
  whereas the raw MCD estimate of scatter is their covariance matrix,
  multiplied by a consistency factor and (if \code{use.correction} is
  true) a finite sample correction factor, to make it consistent at the
  normal model and unbiased at small samples.  Both rescaling factors
  (consistency and finite sample) are returned in the length-2 vector
  \code{raw.cnp2}.

  The implementation of \code{covMcd} uses the Fast MCD algorithm of
  Rousseeuw and Van Driessen (1999) to approximate the minimum
  covariance determinant estimator.

  Based on these raw MCD estimates, (unless argument \code{raw.only} is
  true), a reweighting step is performed, i.e., \code{V <- \link{cov.wt}(x,w)},
  where \code{w} are weights determined by \dQuote{outlyingness} with
  the respect to the scaled raw MCD.  Again, a consistency factor and
  (if \code{use.correction} is true) a finite sample correction factor
  are applied.
  The reweighted covariance is typically considerably more efficient
  than the raw one, see Pison et al. (2002).

  The two rescaling factors for the reweighted estimates are returned in
  \code{cnp2}.  Details for the computation of the finite sample
  correction factors can be found in Pison et al. (2002).
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}, based on
  work written for S-plus by Peter Rousseeuw and Katrien van Driessen
  from University of Antwerp.% no E-mails for spam-protection

  Visibility of (formerly internal) tuning parameters, notably
  \code{wgtFUN()}: Martin Maechler
}
\value{
  An object of class \code{"mcd"} which is basically a
  \code{\link{list}} with components
  \item{center}{the final estimate of location.}
  \item{cov}{the final estimate of scatter.}
  \item{cor}{the (final) estimate of the correlation matrix (only if
    \code{cor = TRUE}).}
  \item{crit}{the value of the criterion, i.e. the determinant.}
  \item{best}{the best subset found and used for computing the raw
    estimates, with \code{\link{length}(best) == quan =
      \link{h.alpha.n}(alpha,n,p)}.}
  \item{mah}{mahalanobis distances of the observations using the final
    estimate of the location and scatter.}
  \item{mcd.wt}{weights of the observations using the final estimate of
    the location and scatter.}
  \item{cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the final estimate of the covariance matrix.}
  \item{raw.center}{the raw (not reweighted) estimate of location.}
  \item{raw.cov}{the raw (not reweighted) estimate of scatter.}
  \item{raw.mah}{mahalanobis distances of the observations based on the
    raw estimate of the location and scatter.}
  \item{raw.weights}{weights of the observations based on the raw
    estimate of the location and scatter.}
  \item{raw.cnp2}{a vector of length two containing the consistency
    correction factor and the finite sample correction factor of
    the raw estimate of the covariance matrix.}
  \item{X}{the input data as numeric matrix, without \code{\link{NA}}s.}
  \item{n.obs}{total number of observations.}
  \item{alpha}{the size of the subsets over which the determinant is
    minimized (the default is \eqn{(n+p+1)/2}).}
  \item{quan}{the number of observations, \eqn{h}, on which the MCD is
    based.  If \code{quan} equals \code{n.obs}, the MCD is the classical
    covariance matrix.}
  \item{method}{character string naming the method (Minimum Covariance
    Determinant).}
  \item{call}{the call used (see \code{\link{match.call}}).}
}
\references{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection.} Wiley.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.

  Pison, G., Van Aelst, S., and Willems, G. (2002),
  Small Sample Corrections for LTS and MCD,
  \emph{Metrika} \bold{55}, 111--123.% ~/save/papers/robust-diverse/Pison_VanAelst_Willems.pdf
}
\seealso{
  \code{\link[MASS]{cov.mcd}} from package \CRANpkg{MASS};
  \code{\link{covOGK}} as cheaper alternative for larger dimensions.

  \code{\link[robustX]{BACON}} and \code{\link[robustX]{covNNC}},
  from package \CRANpkg{robustX};
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
covMcd(hbk.x)

## the following three statements are equivalent
c1 <- covMcd(hbk.x, alpha = 0.75)
c2 <- covMcd(hbk.x, control = rrcov.control(alpha = 0.75))
## direct specification overrides control one:
c3 <- covMcd(hbk.x, alpha = 0.75,
             control = rrcov.control(alpha=0.95))
c1
}
\keyword{robust}
\keyword{multivariate}
