\name{lmrob.control}
\title{Tuning Parameters for lmrob() and Auxiliaries}
\encoding{utf8}
\alias{lmrob.control}
\alias{lmrob.control.sfs}
\description{
  Tuning parameters for \code{\link{lmrob}}, the MM-type regression
  estimator and the associated S-, M- and D-estimators.  Using
  \code{setting="KS2011"} sets the defaults as suggested by
  Koller and Stahel (2011).
}
\usage{
lmrob.control(setting, seed = NULL, nResample = 500,
              tuning.chi = NULL, bb = 0.5, tuning.psi = NULL,
              max.it = 50, groups = 5, n.group = 400,
              k.fast.s = 1, best.r.s = 2,
              k.max = 200, maxit.scale = 200, k.m_s = 20,
              refine.tol = 1e-7, rel.tol = 1e-7, solve.tol = 1e-7,
              trace.lev = 0,
              mts = 1000, subsampling = c("nonsingular", "simple"),
              compute.rd = FALSE, method = 'MM',
              psi = c('bisquare', 'lqq', 'welsh', 'optimal', 'hampel', 'ggw'),
              numpoints = 10, cov = NULL,
              split.type = c("f", "fi", "fii"), fast.s.large.n = 2000, ...)
}
\arguments{
  \item{setting}{a string specifying alternative default values.  Leave
    empty for the defaults or use \code{"KS2011"} for the defaults
    suggested by Koller and Stahel (2011).  See \emph{Details}.}
  \item{seed}{\code{NULL} or an integer vector compatible with
    \code{\link{.Random.seed}}: the seed to be used for random
    re-sampling used in obtaining candidates for the initial
    S-estimator.  The current value of \code{.Random.seed} will be
    preserved if \code{seed} is set, i.e. non-\code{NULL};
    otherwise, as by default, \code{.Random.seed} will be used and
    modified as usual from calls to \code{\link{runif}()} etc.
  }
  \item{nResample}{number of re-sampling candidates to be
    used to find the initial S-estimator.  Currently defaults to 500
    which works well in most situations (see references).}
  \item{tuning.chi}{tuning constant vector for the S-estimator.  If
    \code{NULL}, as by default, sensible defaults are set (depending on
    \code{psi}) to yield a 50\% breakdown estimator.  See \emph{Details}.}
  \item{bb}{expected value under the normal model of the
    \dQuote{chi} (rather \eqn{\rho (rho)}{rho}) function with tuning
    constant equal to \code{tuning.chi}.  This is used to compute the
    S-estimator.}
  \item{tuning.psi}{tuning constant vector for the redescending
    M-estimator.  If \code{NULL}, as by default, this is set (depending
    on \code{psi}) to yield an estimator with asymptotic efficiency of
    95\% for normal errors.  See \emph{Details}.}
  \item{max.it}{integer specifying the maximum number of IRWLS iterations.}
  \item{groups}{(for the fast-S algorithm): Number of
    random subsets to use when the data set is large.}
  \item{n.group}{(for the fast-S algorithm): Size of each of the
    \code{groups} above.  Note that this must be at least \eqn{p}.}
  \item{k.fast.s}{(for the fast-S algorithm): Number of
    local improvement steps (\dQuote{\emph{I-steps}}) for each
    re-sampling candidate.}
  \item{k.m_s}{(for the M-S algorithm): specifies after how many
    unsucessful refinement steps the algorithm stops.}
  \item{best.r.s}{(for the fast-S algorithm): Number of
    of best candidates to be iterated further (i.e.,
    \dQuote{\emph{\bold{r}efined}}); is denoted \eqn{t} in
    Salibian-Barrera & Yohai(2006).}
  \item{k.max}{(for the fast-S algorithm): maximal number of
    refinement steps for the \dQuote{fully} iterated best candidates.}
  \item{maxit.scale}{integer specifying the maximum number of C level
    \code{find_scale()} iterations.}
  \item{refine.tol}{(for the fast-S algorithm): relative convergence
    tolerance for the fully iterated best candidates.}
  \item{rel.tol}{(for the RWLS iterations of the MM algorithm): relative
    convergence tolerance for the parameter vector.}
  \item{solve.tol}{(for the S algorithm): relative
    tolerance for inversion.  Hence, this corresponds to
    \code{\link{solve.default}()}'s \code{tol}.}
  \item{trace.lev}{integer indicating if the progress of the MM-algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
  \item{mts}{maximum number of samples to try in subsampling
    algorithm.}
  \item{subsampling}{type of subsampling to be used, a string:
    \code{"simple"} for simple subsampling (default prior to version 0.9),
    \code{"nonsingular"} for nonsingular subsampling.  See also
    \code{\link{lmrob.S}}.}
  \item{compute.rd}{logical indicating if robust distances (based on
    the MCD robust covariance estimator \code{\link{covMcd}}) are to be
    computed for the robust diagnostic plots.  This may take some
    time to finish, particularly for large data sets, and can lead to
    singularity problems when there are \code{\link{factor}} explanatory
    variables (with many levels, or levels with \dQuote{few}
    observations). Hence, is \code{FALSE} by default.}
  \item{method}{string specifying the estimator-chain. \code{MM}
    is interpreted as \code{SM}. See \emph{Details} of
    \code{\link{lmrob}} for a description of the possible values.}
  \item{psi}{string specifying the type \eqn{\psi}-function
    used.  See \emph{Details} of \code{\link{lmrob}}.  Defaults to
    \code{"bisquare"} for S and MM-estimates, otherwise \code{"lqq"}.}
  \item{numpoints}{number of points used in Gauss quadrature.}
  \item{cov}{function or string with function name to be used to
    calculate covariance matrix estimate.  The default is
    \code{if(method \%in\% c('SM', 'MM')) ".vcov.avar1" else ".vcov.w"}.
    See \emph{Details} of \code{\link{lmrob}}.}
  \item{split.type}{determines how categorical and continuous variables
    are split. See \code{\link{splitFrame}}.}
  \item{fast.s.large.n}{minimum number of observations required to
    switch from ordinary \dQuote{fast S} algorithm to an efficient
    \dQuote{large n} strategy.}
  \item{...}{further arguments to be added as \code{\link{list}}
    components to the result.}
}
\value{
  a named \code{\link{list}} with over twenty components, corresponding
  to the arguments, where \code{tuning.psi} and \code{tuning.chi} are
  typically computed, see above.
}
\details{The option \code{setting="KS2011"} alters the default
  arguments.  They are changed to \code{method = 'SMDM', psi = 'lqq',
    max.it = 500, k.max = 2000, cov = '.vcov.w'}.  The defaults of all
  the remaining arguments are not changed.

  By default, \code{tuning.chi} and \code{tuning.psi} are set to
  yield an MM-estimate with break-down point \eqn{0.5} and efficiency of
  \eqn{95\%}  at the normal. They are:
  \tabular{rll}{
    \code{psi} \tab \code{tuning.chi} \tab \code{tuning.psi} \cr
    \code{bisquare} \tab \code{1.54764} \tab \code{4.685061} \cr
    \code{welsh} \tab \code{0.5773502} \tab \code{2.11} \cr
    \code{ggw} \tab \code{c(-0.5, 1.5, NA, 0.5)} \tab
    \code{c(-0.5, 1.5, 0.95, NA)} \cr
    \code{lqq} \tab \code{c(-0.5, 1.5, NA, 0.5)} \tab
    \code{c(-0.5, 1.5, 0.95, NA)} \cr
    \code{optimal} \tab \code{0.4047} \tab \code{1.060158} \cr
    \code{hampel} \tab \code{c(1.5, 3.5, 8)*0.2119163} \tab
    \code{c(1.5, 3.5, 8)*0.9014}
  }
  The values for the tuning constant for the \code{ggw} psi function are
  hard coded. The constants vector has four elements: minimal slope, b
  (controlling the bend at the maximum of the curve), efficiency,
  break-down point. Use \code{NA} for an unspecified value, see examples
  in the tables.

  The constants for the \code{hampel} psi function are chosen to have a
  redescending slope of \eqn{-1/3}. Constants for a slope of \eqn{-1/2}
  would be
  \tabular{rll}{
    \code{psi} \tab \code{tuning.chi} \tab \code{tuning.psi} \cr
    \code{hampel} \tab \code{c(2, 4, 8) * 0.1981319} \tab
    \code{c(2, 4, 8) * 0.690794}
  }

  Alternative coefficients for an efficiency of \eqn{85\%} at the normal
  are given in the table below.
  \tabular{rl}{
    \code{psi} \tab \code{tuning.psi} \cr
    \code{bisquare} \tab \code{3.443689} \cr
    \code{welsh} \tab \code{1.456} \cr
    \code{ggw}, \code{lqq} \tab \code{c(-0.5, 1.5, 0.85, NA)} \cr
    \code{optimal} \tab \code{0.8684} \cr
    \code{hampel} (-1/3) \tab \code{c(1.5, 3.5, 8)* 0.5704545} \cr
    \code{hampel} (-1/2) \tab \code{c( 2,  4,  8) * 0.4769578}
  }
}
\references{
  Koller, M. and Stahel, W.A. (2011), Sharpening Wald-type inference in
  robust regression for small samples, \emph{Computational Statistics &
  Data Analysis} \bold{55}(8), 2504--2515.
}
\author{ Matias Salibian-Barrera, Martin Maechler and Manuel Koller}
\seealso{ \code{\link{lmrob}}, also for references and examples.
}
\examples{
## Show the default settings:
str(lmrob.control())

## Artificial data for a  simple  "robust t test":
set.seed(17)
y <- y0 <- rnorm(200)
y[sample(200,20)] <- 100*rnorm(20)
gr <- as.factor(rbinom(200, 1, prob = 1/8))
lmrob(y0 ~ 0+gr)

## Use  Koller & Stahel(2011)'s recommendation but a larger  'max.it':
str(ctrl <- lmrob.control("KS2011", max.it = 1000))
}
\keyword{robust}
\keyword{regression}
