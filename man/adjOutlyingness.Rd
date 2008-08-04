\name{adjOutlyingness}
\alias{adjOutlyingness}
\title{Compute Skewness-adjusted Multivariate Outlyingness}
\description{
  For an \eqn{n \times p}{n * p} data matrix (or data frame) \code{x},
  compute the \dQuote{\emph{outlyingness}} of all \eqn{n} observations.
  Outlyingness here is a generalization of the Donoho-Stahel
  outlyingness measure, where skewness is taken into account via the
  medcouple, \code{\link{mc}()}.
}
\usage{
adjOutlyingness(x, ndir = 250, clower = 3, cupper = 4,
                alpha.cutoff = 0.75, coef = 1.5, qr.tol = 1e-12)
}
\arguments{
  \item{x}{a numeric \code{\link{matrix}} or \code{\link{data.frame}}.}
  \item{ndir}{positive integer specifying the number of directions that
    should be searched.}
  \item{clower, cupper}{the constant to be used for the lower and upper
    tails, in order to transform the data towards symmetry.}
  %%% MM: How can I modify it to the  *NON*-adjusted case ?
  \item{alpha.cutoff}{number in (0,1) specifying the quantiles
    \eqn{(\alpha, 1-\alpha)} which determine the \dQuote{outlier}
    cutoff.}
  \item{coef}{positive number specifying the factor with which the
    interquartile range (\code{\link{IQR}}) is multiplied to determine
    \sQuote{boxplot hinges}-like upper and lower bounds.}
  \item{qr.tol}{positive tolerance to be used for \code{\link{qr}} and
    \code{\link{solve.qr}} for determining the \code{ndir} directions
    each determined by a random sample of \eqn{p} (out of \eqn{n})
    observations.}
}
\details{

  \bold{FIXME}:  Details in the comment of the Matlab code;
  also in the reference(s).
%% SEE ~/R/MM/STATISTICS/robust/MC/mcmatl/adjoutlyingness.m
%% ---- which has notes about input/output etc of the corresponding
%%      Matlab code

 The method as described can be useful as preprocessing in
 FASTICA (\url{http://www.cis.hut.fi/projects/ica/fastica/}
}
\value{
  a list with components
  \item{adjout}{numeric of \code{length(n)} giving the adjusted
    outlyingness of each observation.}
  \item{cutoff}{cutoff for \dQuote{outlier} with respect to the adjusted
    outlyingnesses, and depending on \code{alpha.cutoff}.}
  \item{nonOut}{logical of \code{length(n)}, \code{TRUE} when the
    corresponding observation is \bold{non}-outlying with respect to the
    cutoff and the adjusted outlyingnesses.}
}
\references{
  Brys, G., Hubert, M., and Rousseeuw, P.J. (2005)
  A Robustification of Independent Component Analysis;
  \emph{Journal of Chemometrics}, \bold{19}, 1--12.

  For the up-to-date reference, please consult
  \url{http://wis.kuleuven.be/stat/robust.html}
}
\author{Guy Brys; help page and improvements by Martin Maechler}
\seealso{the adjusted boxplot, \code{\link{adjbox}} and the medcouple,
  \code{\link{mc}}.
}
\examples{
## An Example with bad condition number and "border case" outliers

if(FALSE) {## Not yet ok, because of bug in adjOutl
  dim(longley)
  set.seed(1) ## result is random %% and there's a bug - FIXME! -- try set.seed(3)
  ao1 <- adjOutlyingness(longley)
  ## which are not outlying ?
  table(ao1$nonOut)  ## all of them
  stopifnot(all(ao1$nonOut))
}

## An Example with outliers :

dim(hbk)
set.seed(1)
ao.hbk <- adjOutlyingness(hbk)
str(ao.hbk)
hist(ao.hbk $adjout)## really two groups
table(ao.hbk$nonOut)## 14 outliers, 61 non-outliers:
## outliers are :
which(! ao.hbk$nonOut) # 1 .. 14   --- but not for all random seeds!

## here, they are the same as found by (much faster) MCD:
cc <- covMcd(hbk)
stopifnot(all(cc$mcd.wt == ao.hbk$nonOut))

## This is revealing (about 1--2 cases, where outliers are *not* == 1:14
##  but needs almost 1 [sec] per call:
if(interactive()) {
  for(i in 1:30) {
    print(system.time(ao.hbk <- adjOutlyingness(hbk)))
    if(!identical(iout <- which(!ao.hbk$nonOut), 1:14)) {
	 cat("Outliers:\n"); print(iout)
    }
  }
}

}
\keyword{robust}
\keyword{multivariate}
