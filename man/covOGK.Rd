\name{covOGK}
\alias{covOGK}
\alias{covGK}
\alias{hard.rejection}
%
\title{Orthogonalized Gnanadesikan-Kettenring (OGK) Covariance Matrix Estimation}
\description{
  Computes the orthogonalized pairwise covariance matrix estimate described in
  in Maronna and Zamar (2002).  The pairwise proposal goes back to
  Gnanadesikan and Kettenring (1972).
}
\usage{
covOGK(X, n.iter, sigmamu, rcov = covGK, weight.fn,
       keep.data = FALSE, ...)
}
\arguments{
  \item{X}{data in something that can be coerced into a numeric matrix.}
  \item{n.iter}{number of orthogonalization iterations.  Usually 1 or 2;
    values greater than 2 are unlikely to have any significant effect on
    the estimate (other than increasing the computing time).}
  \item{sigmamu}{a function that computes univariate robust location and
    scale estimates.  By default \code{sigmamu} should return a single
    numeric value containing the robust scale (standard deviation)
    estimate.  When \code{mu.too} is true, \code{sigmamu()} should
    return a numeric vector of length 2 containing robust location and
    scale estimates.  See \code{\link{scaleTau2}} for an example.}
  \item{rcov}{function that computes a robust covariance estimate
    between two vectors.  The default, Gnanadesikan-Kettenring's
    \code{covGK}, is simply \eqn{(s^2(X+Y) - s^2(X-Y))/4} where
    \eqn{s()} is the scale estimate \code{sigmamu()}.}
  \item{weight.fn}{a function of the robust distances and the number of
    variables \eqn{p} to compute the weights used in the reweighting step.}
  \item{keep.data}{logical indicating if the (untransformed) data matrix
    \code{X} should be kept as part of the result.}
  \item{\dots}{additional arguments to be passed to \code{sigmamu()} and
    \code{weight.fn()}.}
}
\details{
  Typical default values for the \emph{function} arguments
  \code{sigmamu}, \code{rcov}, and \code{weight.fn}, are
  available as well, see the \emph{Examples} below,
  \bold{but} their names and calling sequences are
  still subject to discussion and may be changed in the future.
}
\value{
  currently a list with components
  \item{center}{robust location: numeric vector of length \eqn{p}.}
  \item{cov}{robust covariance matrix estimate: \eqn{p\times p}{p x p}
    matrix.}
  \item{wcenter, wcov}{re-\bold{w}eighted versions of \code{center} and
    \code{cov}.}
  \item{weights}{the robustness weights used.}
  \item{distances}{the mahalanobis distances computed using
    \code{center} and \code{cov}.}

  \dots\dots

  \bold{but note that this might be radically changed to returning an
    S4 classed object!}
}
\references{
  Maronna, R.A. and Zamar, R.H. (2002)
  Robust estimates of location and dispersion of high-dimensional datasets;
  \emph{Technometrics} \bold{44}(4), 307--317.

  Gnanadesikan, R. and John R. Kettenring (1972)
  Robust estimates, residuals, and outlier detection with multiresponse data.
  \emph{Biometrics} \bold{28}, 81--124.
}
\author{Kjell Konis \email{konis@stats.ox.ac.uk}, with modifications by
    Martin Maechler.}
\seealso{\code{\link{scaleTau2}},
  \code{\link{covMcd}}, \code{\link[MASS]{cov.rob}}.
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])

cO1 <- covOGK(hbk.x, n.iter = 2,
              sigmamu = scaleTau2, weight.fn = hard.rejection)
%% FIXME: Add this!
% scaleQn <- function (...)  { .... ; c(Qn, median) }
% cO2 <- covOGK(hbk.x, n.it = 2, sigmamu = scaleQn, weight.fn = hard.rejection)

%% FIXME: Add time comparison, here or in "vignette", "demo", "...
}
\keyword{robust}
\keyword{multivariate}
