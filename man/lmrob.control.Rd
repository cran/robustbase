\name{lmrob.control}
\title{Tuning Parameters for lmrob() and Auxiliaries}
\encoding{utf8}
\alias{lmrob.control}
\alias{update.lmrobCtrl}
\alias{.Mchi.tuning.default}
\alias{.Mpsi.tuning.default}
\alias{.Mchi.tuning.defaults}
\alias{.Mpsi.tuning.defaults}
\description{
  Tuning parameters for \code{\link{lmrob}}, the MM-type regression
  estimator and the associated S-, M- and D-estimators.  Using
  \code{setting="KS2011"} sets the defaults as suggested by
  Koller and Stahel (2011) and analogously for \code{"KS2014"}.

  The \code{.M*.default} \code{\link{function}}s and
  \code{.M*.defaults} \code{\link{list}}s contain default tuning
  parameters for all the predefined \eqn{\psi}{psi} functions, see also
  \code{\link{Mpsi}}, etc.
}
\usage{
lmrob.control(setting, seed = NULL, nResample = 500,
              tuning.chi = NULL, bb = 0.5, tuning.psi = NULL,
              max.it = 50, groups = 5, n.group = 400,
              k.fast.s = 1, best.r.s = 2,
              k.max = 200, maxit.scale = 200, k.m_s = 20,
              refine.tol = 1e-7, rel.tol = 1e-7, scale.tol = 1e-10, solve.tol = 1e-7,
              zero.tol = 1e-10,
              trace.lev = 0,
              mts = 1000, subsampling = c("nonsingular", "simple"),
              compute.rd = FALSE, method = "MM", psi = "bisquare",
              numpoints = 10, cov = NULL,
              split.type = c("f", "fi", "fii"), fast.s.large.n = 2000,
              # only for outlierStats() :
              eps.outlier = function(nobs) 0.1 / nobs,
              eps.x = function(maxx) .Machine$double.eps^(.75)*maxx,
              compute.outlier.stats = method,
              warn.limit.reject = 0.5,
              warn.limit.meanrw = 0.5, \dots)

\method{update}{lmrobCtrl}(object, \dots)

.Mchi.tuning.defaults
.Mchi.tuning.default(psi)
.Mpsi.tuning.defaults
.Mpsi.tuning.default(psi)
}
\arguments{
  \item{setting}{a string specifying alternative default values.  Leave
    empty for the defaults or use \code{"KS2011"} or \code{"KS2014"}
    for the defaults suggested by Koller and Stahel (2011, 2017).
    See \emph{Details}.}
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
    unsuccessful refinement steps the algorithm stops.}
  \item{best.r.s}{(for the fast-S algorithm): Number of
    of best candidates to be iterated further (i.e.,
    \dQuote{\emph{\bold{r}efined}}); is denoted \eqn{t} in
    Salibian-Barrera & Yohai(2006).}
  \item{k.max}{(for the fast-S algorithm): maximal number of
    refinement steps for the \dQuote{fully} iterated best candidates.}
  \item{maxit.scale}{integer specifying the maximum number of C level
    \code{find_scale()} iterations (in fast-S and M-S algorithms).}
  \item{refine.tol}{(for the fast-S algorithm): relative convergence
    tolerance for the fully iterated best candidates.}
  \item{rel.tol}{(for the RWLS iterations of the MM algorithm): relative
    convergence tolerance for the parameter vector.}
  \item{scale.tol}{(for the scale estimation iterations of the S algorithm): relative
    convergence tolerance for the \code{scale} \eqn{\sigma(.)}.}
  \item{solve.tol}{(for the S algorithm): relative
    tolerance for inversion.  Hence, this corresponds to
    \code{\link{solve.default}()}'s \code{tol}.}
  \item{zero.tol}{for checking 0-residuals in the S algorithm, non-negative number
    \eqn{\epsilon_z}{ez} such that
    \eqn{\{i; \left|\tilde{R}_i\right| \le \epsilon_z\}}{{i; |R~[i]| <= ez}}
    correspond to \eqn{0}-residuals, where \eqn{\tilde{R}_i}{R~[i]} are standardized residuals,
    \eqn{\tilde{R}_i = R_i/s_y}{R~[i] = R[i]/sy} and
    \eqn{s_y = \frac{1}{n} \sum_{i=1}^n \left|y_i\right|}{%
          sy = ave_i |y[i]| = 1/n sum(i=1..n, |y[i]|)}.}
  \item{trace.lev}{integer indicating if the progress of the MM-algorithm
    and the fast-S algorithms, see \code{\link{lmrob.S}},
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
    observations).  Hence, is \code{FALSE} by default.}
  \item{method}{string specifying the estimator-chain. \code{MM}
    is interpreted as \code{SM}.  See \emph{Details} of
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
    are split.  See \code{\link{splitFrame}}.}
  \item{fast.s.large.n}{minimum number of observations required to
    switch from ordinary \dQuote{fast S} algorithm to an efficient
    \dQuote{large n} strategy.}
  \item{eps.outlier}{limit on the robustness weight below which an observation
    is considered to be an outlier.
    Either a \code{numeric(1)} or a function that takes the number of observations as
    an argument.  Used only in \code{\link{summary.lmrob}} and
    \code{\link{outlierStats}}.}
  \item{eps.x}{limit on the absolute value of the elements of the design matrix
    below which an element is considered zero.
    Either a \code{numeric(1)} or a function that takes the maximum absolute value in
    the design matrix as an argument.}
  \item{compute.outlier.stats}{vector of \code{\link{character}}
    strings, each valid to be used as \code{method} argument.  Used to
    specify for which estimators outlier statistics (and warnings)
    should be produced.  Set to empty (\code{NULL} or \code{character(0)})
    if none are required.
    \cr Note that the default is \code{method} which by default is either
    \code{"MM"}, \code{"SM"}, or \code{"SMDM"}; hence using
    \code{compute.outlier.stats = "S"} provides \code{\link{outlierStats}()}
    to a \code{\link{lmrob.S}()} result.}
  \item{warn.limit.reject}{limit of ratio
    \eqn{\#\mbox{rejected} / \#\mbox{obs in level}}{# rejected / # obs in level}
    above (\eqn{\geq}{>=}) which a warning is produced.
    Set to \code{NULL} to disable warning.}
  \item{warn.limit.meanrw}{limit of the mean robustness per factor level
    below which (\eqn{\leq}{<=}) a warning is produced.
    Set to \code{NULL} to disable warning.}

  \item{object}{an \code{"lmrobCtrl"} object, as resulting from a
    \code{lmrob.control(*)} or an \code{update(<lmrobCtrl>, *)} call.}
  \item{\dots}{for \describe{
      \item{\code{lmrob.control()}:}{further arguments to be added as
	\code{\link{list}} components to the result, e.g., those to be used in
	\code{.vcov.w()}.}
      \item{\code{update(object, *)}:}{(named) components from
	\code{object}, to be \emph{modified}, \bold{not} \code{setting = *}.}
  }}
}
\value{
  \code{.Mchi.tuning.default(psi)} and \code{.Mpsi.tuning.default(psi)}
  return a short \code{\link{numeric}} vector of tuning constants which
  are defaults for the corresponding psi-function, see the \emph{Details}.
  They are based on the named \code{\link{list}}s
  \code{.Mchi.tuning.defaults} and \code{.Mpsi.tuning.defaults},
  respectively.

  \code{lmrob.control()} returns a named \code{\link{list}} with over
  twenty components, corresponding to the arguments, where
  \code{tuning.psi} and \code{tuning.chi} are typically computed, as
  \code{.Mpsi.tuning.default(psi)} or \code{.Mchi.tuning.default(psi)},
  respectively.
  It is of \code{\link{class}} \code{"lmrobCtrl"} and we provide
  \code{print()}, \code{\link{update}()} and \code{\link{within}} methods.

  \code{update(<lmrobCtrl>, ....)} does \emph{not} allow a
  \code{setting="<...>"} in \code{....}.
}
\details{The option \code{setting="KS2011"} alters the default
  arguments.  They are changed to \code{method = "SMDM"}, \code{psi = "lqq"},
  \code{max.it = 500}, \code{k.max = 2000}, \code{cov = ".vcov.w"}.
  The defaults of all the remaining arguments are not changed.

  The option \code{setting="KS2014"} builds upon \code{setting="KS2011"}.
  More arguments are changed to \code{best.r.s = 20, k.fast.s = 2,
  nResample = 1000}.  This setting should produce more stable estimates
  for designs with \code{\link{factor}}s.

  By default, and in \code{.Mpsi.tuning.default()} and \code{.Mchi.tuning.default()},
  \code{tuning.chi} and \code{tuning.psi} are set to yield an
  MM-estimate with breakdown point \eqn{0.5} and efficiency of 95\% at
  the normal.

  If numeric \code{tuning.chi} or \code{tuning.psi} are specified, say
  \code{cc}, for \code{psi = "ggw"} or \code{"lqq"},
  \code{\link{.psi.const}(cc, psi)} is used, see its help page.

  To get the defaults, e.g., \code{.Mpsi.tuning.default(psi)} is
  equivalent to but more efficient than the formerly widely used
  \code{lmrob.control(psi = psi)$tuning.psi}.

  These defaults are:
  \tabular{rll}{
    \code{psi}     \tab\code{tuning.chi}               \tab\code{tuning.psi} \cr
    \code{bisquare}\tab\code{1.54764}                  \tab\code{4.685061} \cr
    \code{welsh}   \tab\code{0.5773502}                \tab\code{2.11} \cr
    \code{ggw} 	   \tab\code{c(-0.5, 1.5, NA, 0.5)}    \tab\code{c(-0.5, 1.5, 0.95, NA)} \cr
    \code{lqq}     \tab\code{c(-0.5, 1.5, NA, 0.5)}    \tab\code{c(-0.5, 1.5, 0.95, NA)} \cr
    \code{optimal} \tab\code{0.4047}                   \tab\code{1.060158} \cr
    \code{hampel}  \tab\code{c(1.5, 3.5, 8)*0.2119163} \tab\code{c(1.5, 3.5, 8)*0.9014}
  }
  The values for the tuning constant for the \code{ggw} and \code{lqq}
  psi functions are specified differently here by a vector with four
  elements: minimal slope, b (controlling the bend at the maximum of the curve),
  efficiency, breakdown point.
  Use \code{NA} for an unspecified value of either efficiency or
  breakdown point, see examples in the tables (above and below).
  For these table examples, the respective \dQuote{inner constants} are
  stored precomputed, see \code{\link{.psi.lqq.findc}} for more.

  The constants for the \code{"hampel"} psi function are chosen to have a
  redescending slope of \eqn{-1/3}.  Constants for a slope of \eqn{-1/2}
  would be
  \tabular{rll}{
    \code{psi}     \tab\code{tuning.chi}             \tab\code{tuning.psi} \cr
    \code{"hampel"}\tab\code{c(2, 4, 8) * 0.1981319} \tab\code{c(2, 4, 8) * 0.690794}
  }

  Alternative coefficients for an efficiency of 85\%
  at the normal are given in the table below.
  \tabular{rl}{
    \code{psi}		  \tab\code{tuning.psi} \cr
    \code{bisquare}	  \tab\code{3.443689} \cr
    \code{welsh}	  \tab\code{1.456} \cr
    \code{ggw}, \code{lqq}\tab\code{c(-0.5, 1.5, 0.85, NA)} \cr
    \code{optimal}        \tab\code{0.8684} \cr
    \code{hampel} (-1/3)  \tab\code{c(1.5, 3.5, 8)* 0.5704545} \cr
    \code{hampel} (-1/2)  \tab\code{c( 2,  4,  8) * 0.4769578}
  }
}
\references{
  Koller, M. and Stahel, W.A. (2011)
  Sharpening Wald-type inference in robust regression for small samples.
  \emph{Computational Statistics & Data Analysis} \bold{55}(8), 2504--2515.

  Koller, M. and Stahel, W.A. (2017)
  Nonsingular subsampling for regression S estimators with categorical predictors,
  \emph{Computational Statistics} \bold{32}(2): 631--646.
  \doi{10.1007/s00180-016-0679-x}.
  Referred as \code{"KS2014"} everywhere in \pkg{robustbase};  A shorter first
  version, Koller (2012) has been available from \url{https://arxiv.org/abs/1208.5595}.
}
\author{Matias Salibian-Barrera, Martin Maechler and Manuel Koller}
\seealso{  \code{\link{Mpsi}}, etc, for the (fast!) psi function computations;
  \code{\link{lmrob}}, also for references and examples.
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

str(.Mpsi.tuning.defaults)
stopifnot(identical(.Mpsi.tuning.defaults,
                   sapply(names(.Mpsi.tuning.defaults),
                          .Mpsi.tuning.default)))
## Containing (names!) all our (pre-defined) redescenders:
str(.Mchi.tuning.defaults)

## Difference between settings:
Cdef <- lmrob.control()
C11 <- lmrob.control("KS2011")
C14 <- lmrob.control("KS2014")
str(C14)
## Differences:
diffD <- names(which(!mapply(identical, Cdef,C11, ignore.environment=TRUE)))
diffC <- names(which(!mapply(identical, C11, C14, ignore.environment=TRUE)))
## KS2011 vs KS2014:  Apart from `setting` itself, they only differ in three places:
cbind(KS11 = unlist(C11[diffC[-1]]),
      KS14 = unlist(C14[diffC[-1]]))
##           KS11 KS14
## nResample  500 1000
## best.r.s     2   20
## k.fast.s     1    2
## default vs KS2011: a bit more: setting + 8
str2simpLang <-  function(x) {
    r <- if(is.null(x)) quote((NULL)) else str2lang(deparse(x))
    if(is.call(r)) format(r) else r
}
cbind(deflt= lapply(Cdef[diffD], str2simpLang),
      KS11 = lapply(C11 [diffD], str2simpLang))

## update()ing a lmrob.control() , e.g.,
C14mod <- update(C14, trace.lev = 2) # the same as
C14m.d <- C14; C14m.d$trace.lev <- 2
stopifnot(identical(C14mod, C14m.d))
## changing psi --> updates tuning.{psi,chi}:
C14mp <- update(C14, psi = "hampel", seed=101)
## updating 'method' is "smart" :
C.SMDM <- update(Cdef, method="SMDM")
all.equal(Cdef, C.SMDM) # changed also psi, tuning.{psi,chi} and cov !
chgd <- c("method", "psi", "tuning.chi",  "tuning.psi", "cov")
str(Cdef  [chgd])
str(C.SMDM[chgd])
C14m <- update(C14, method="SMM")
(ae <- all.equal(C14, C14mp))# changed tuning.psi & tuning.chi, too
stopifnot(exprs = {
    identical(C14, update(C14, method="SMDM")) # no change!
    identical(c("psi", "seed", "tuning.chi", "tuning.psi"),
              sort(gsub("[^.[:alpha:]]", "", sub(":.*", "", sub("^Component ", "", ae)))))
    identical(C14m, local({C <- C14; C$method <- "SMM"; C}))
})
##
try( update(C14, setting="KS2011") ) #--> Error: .. not allowed
\dontshow{tools::assertError(update(C14, setting="KS2011"))}
}
\keyword{robust}
\keyword{regression}
