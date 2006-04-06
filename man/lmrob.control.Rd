\name{lmrob.control}
\alias{lmrob.control}
\title{ Tuning parameters for lmrob }
\description{
  Tuning parameters for, \code{\link{lmrob}}, the MM-regression
  estimator and the associated S-estimator.
}
\usage{
lmrob.control(seed = 37, Nres = 500,
	tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061,
	max.it = 50, groups = 5, n.group = 400, k.fast.s = 1,
        compute.rd = TRUE)
}
\arguments{
  \item{seed}{random seed for the re-samples used in obtaining candiates
    for the initial S-estimator.  The default, \code{37} used to be
    frozen in the underlying C code.}
  \item{Nres}{number of re-sampling candidates to be
    used to find the initial S-estimator.  Currently defaults to 500
    which works well in most situations (see References below).
    User-choice capability will be added in future releases.}
  \item{tuning.chi}{tuning constant for the S-estimator.
    The default, \code{1.54764}, yields a 50\% breakdown estimator.}
  \item{bb}{expected value under the normal model of the
    "chi" function with tuning constant equal to
    \code{tuning.chi}.  This is used to compute the S-estimator.}
  \item{tuning.psi}{tuning constant for the re-descending M-estimator.
    The choice \code{4.685061} yields an estimator with asymptotic
    efficiency of 95\% for normal errors.}
  \item{max.it}{integer specifying the maximum number of IRWLS iterations.}
  \item{groups}{This parameter is for the fast-S algorithm.  Number of
    random subsets to use when the data set is large.}
  \item{n.group}{This parameter is for the fast-S algorithm.  Size of
    each of the \code{groups} above.}
  \item{k.fast.s}{This parameter is for the fast-S algorithm.  Number of
    local improvement steps for each re-sampling candidate.}
  \item{compute.rd}{logical indicating if robust distances (based on
    the MCD robust covariance estimator \code{\link{covMcd}}) are to be
    computed for the robust diagnostic plots.  This may take some
    time to finish, particularly for large data sets.}
}
\author{ Matias Salibian-Barrera }
\seealso{ \code{\link{lmrob}}, also for references and examples.
}
\keyword{robust}
\keyword{regression}

