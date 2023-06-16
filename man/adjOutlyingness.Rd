\name{adjOutlyingness}
\alias{adjOutlyingness}
\title{Compute (Skewness-adjusted) Multivariate Outlyingness}
\description{
  For an \eqn{n \times p}{n * p} data matrix (or data frame) \code{x},
  compute the \dQuote{\emph{outlyingness}} of all \eqn{n} observations.
  Outlyingness here is a generalization of the Donoho-Stahel
  outlyingness measure, where skewness is taken into account via the
  medcouple, \code{\link{mc}()}.
}
\usage{
adjOutlyingness(x, ndir = 250, p.samp = p, clower = 4, cupper = 3,
                IQRtype = 7,
                alpha.cutoff = 0.75, coef = 1.5,
                qr.tol = 1e-12, keep.tol = 1e-12,
                only.outlyingness = FALSE, maxit.mult = max(100, p),
                trace.lev = 0,
                mcReflect = n <= 100, mcScale = TRUE, mcMaxit = 2*maxit.mult,
                mcEps1 = 1e-12, mcEps2 = 1e-15,
                mcTrace = max(0, trace.lev-1))
}
\arguments{
  \item{x}{a numeric \eqn{n \times p}{n * p} \code{\link{matrix}} or
    \code{\link{data.frame}},
    which must be of full rank \eqn{p}.}
  \item{ndir}{positive integer specifying the number of directions that
    should be searched.}
  \item{p.samp}{the sample size to use for finding good random
    directions, must be at least \code{p}.  The default, \code{p} had
    been hard coded previously.}
  \item{clower, cupper}{the constant to be used for the lower and upper
    tails, in order to transform the data towards symmetry.  You can set
    \code{clower = 0, cupper = 0} to get the \emph{non}-adjusted,
    i.e., classical (\dQuote{central} or \dQuote{symmetric})
    outlyingness.  In that case, \code{\link{mc}()} is not used.}
  \item{IQRtype}{a number from \code{1:9}, denoting type of empirical
    quantile computation for the \code{\link{IQR}()}.  The default 7
    corresponds to \code{\link{quantile}}'s and \code{\link{IQR}}'s
    default.  MM has evidence that \code{type=6} would be a bit more stable
    for small sample size.}
  \item{alpha.cutoff}{number in (0,1) specifying the quantiles
    \eqn{(\alpha, 1-\alpha)} which determine the \dQuote{outlier}
    cutoff.  The default, using quartiles, corresponds to the definition
    of the medcouple (\code{\link{mc}}), but there is no stringent
    reason for using the same alpha for the outlier cutoff.}
  \item{coef}{positive number specifying the factor with which the
    interquartile range (\code{\link{IQR}}) is multiplied to determine
    \sQuote{boxplot hinges}-like upper and lower bounds.}
  \item{qr.tol}{positive tolerance to be used for \code{\link{qr}} and
    \code{\link{solve.qr}} for determining the \code{ndir} directions,
    each determined by a random sample of \eqn{p} (out of \eqn{n})
    observations.  Note that the default \eqn{10^{-12}} is rather small,
    and \code{\link{qr}}'s default \code{= 1e-7} may be more appropriate.}
  \item{keep.tol}{positive tolerance to determine which of the sample
    direction should be kept, namely only those for which
    \eqn{\|x\| \cdot \|B\|}{||x|| * ||B||} is larger than \code{keep.tol}.}
  \item{only.outlyingness}{logical indicating if the final outlier
    determination should be skipped.  In that case, a vector is
    returned, see \sQuote{Value:} below.}
  \item{maxit.mult}{integer factor; \code{maxit <- maxit.mult * ndir}
    will determine the maximal number of direction searching
    iterations.  May need to be increased for higher dimensional data,
    though increasing \code{ndir} may be more important.}
  \item{trace.lev}{an integer, if positive allows to monitor the
    direction search.}
  %% new (Aug-Dec 2020), related to mc():  see >>> ./mc.Rd <<<
  \item{mcReflect}{passed as \code{doReflect} to \code{\link{mc}()}.}
  \item{mcScale}{passed as \code{doScale} to \code{\link{mc}()}.}
  \item{mcMaxit}{passed as \code{maxit} to \code{\link{mc}()}.}
  \item{mcEps1}{passed as \code{eps1} to \code{\link{mc}()}; the default is slightly
    looser (100 larger) than the default for \code{mc}().}
  \item{mcEps2}{passed as \code{eps2} to \code{\link{mc}()}.}
  \item{mcTrace}{passed as \code{trace.lev} to \code{\link{mc}()}.}
}
\note{
  % During Lukas Graz' Master's thesis (Spring 2021), it finally became clear to MM:
  If there are too many degrees of freedom for the projections, i.e., when
  \eqn{n \le 4p}{n <= 4*p}, the current definition of adjusted outlyingness
  is ill-posed, as one of the projections may lead to a denominator
  (quartile difference) of zero, and hence formally an adjusted
  outlyingness of infinity.
  The current implementation avoids \code{Inf} results, but will return
  seemingly random \code{adjout} values of around \eqn{10^{14} -- 10^{15}} which may
  be completely misleading, see, e.g., the \code{longley} data example.

  The result is \emph{random} as it depends on the sample of
  \code{ndir} directions chosen; specifically, to get sub samples the algorithm uses
  \code{\link{sample.int}(n, p.samp)}
  which from \R version 3.6.0 depends on
  \code{\link{RNGkind}(*, sample.kind)}.  Exact reproducibility of results
  from \R versions 3.5.3 and earlier, requires setting
  \code{\link{RNGversion}("3.5.0")}.% same text in ./glmrob.Rd ("MT")
  In any case, do use \code{\link{set.seed}()} yourself
  for reproducibility!

 Till Aug/Oct. 2014, the default values for \code{clower} and \code{cupper} were
 accidentally reversed, and the signs inside \code{exp(.)} where swapped
 in the (now corrected) two expressions \preformatted{
 tup <- Q3 + coef * IQR * exp(.... + clower * tmc * (tmc < 0))
 tlo <- Q1 - coef * IQR * exp(.... - cupper * tmc * (tmc < 0))
}
 already in the code from Antwerpen (\file{mcrsoft/adjoutlingness.R}),
 contrary to the published reference.

 Further, the original algorithm had not been scale-equivariant in the
 direction construction, which has been amended in 2014-10 as well.

 The results, including diagnosed outliers, therefore have changed,
 typically slightly, since \pkg{robustbase} version 0.92-0.
}
\details{
  \bold{FIXME}:  Details in the comment of the Matlab code;
  also in the reference(s).
  %% SEE /u/maechler/R/MM/STATISTICS/robust/MC/mcmatl/adjoutlyingness.m
  %% ---- which has notes about input/output etc of the corresponding
  %%      Matlab code

  The method as described can be useful as preprocessing in
  FASTICA (\url{http://research.ics.aalto.fi/ica/fastica/}
  see also the \R package \CRANpkg{fastICA}.
}
\value{
  If \code{only.outlyingness} is true, a vector \code{adjout},
  otherwise, as by default, a list with components
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

  Hubert, M., Van der Veeken, S. (2008)
  Outlier detection for skewed data;
  \emph{Journal of Chemometrics} \bold{22}, 235--246;
  \doi{10.1002/cem.1123}.
  %% preprint \url{http://wis.kuleuven.be/stat/robust/papers/2008/outlierdetectionskeweddata-revision.pdf}
  %%MM: Journal-pdf  ~/save/papers/robust-diverse/Hubert_VdV_skewed-Chemom_2008.pdf
  %%MM: Compstat 2010: Slides (of talk) and paper of Mia H:
  %% ~/save/papers/robust-diverse/Hubert_skewed-CS2010-slides.pdf  and
  %% ~/save/papers/robust-diverse/Hubert_skewed-CS2010-paper.pdf (slides are better !!)

  For the up-to-date reference, please consult
 %\url{https://wis.kuleuven.be/stat/robust}% to sound modern:
  \url{https://wis.kuleuven.be/statdatascience/robust}
}
\author{Guy Brys; help page and improvements by Martin Maechler}
\seealso{the adjusted boxplot, \code{\link{adjbox}} and the medcouple,
  \code{\link{mc}}.
}
\examples{
## An Example with bad condition number and "border case" outliers

dim(longley) # 16 x 7  // set seed, as result is random :
set.seed(31)
ao1 <- adjOutlyingness(longley, mcScale=FALSE)
## which are outlying ?
which(!ao1$nonOut) ## for this seed, two: "1956", "1957"; (often: none)
## For seeds 1:100, we observe (Linux 64b)
if(FALSE) {
  adjO <- sapply(1:100, function(iSeed) {
            set.seed(iSeed); adjOutlyingness(longley)$nonOut })
  table(nrow(longley) - colSums(adjO))
}
## #{outl.}:  0  1  2  3
## #{cases}: 74 17  6  3


## An Example with outliers :

dim(hbk)
set.seed(1)
ao.hbk <- adjOutlyingness(hbk)
str(ao.hbk)
hist(ao.hbk $adjout)## really two groups
table(ao.hbk$nonOut)## 14 outliers, 61 non-outliers:
## outliers are :
which(! ao.hbk$nonOut) # 1 .. 14   --- but not for all random seeds!

## here, they are(*) the same as found by (much faster) MCD:
## *) only "almost", since the 2023-05 change to covMcd() 
cc <- covMcd(hbk)
table(cc = cc$mcd.wt, ao = ao.hbk$nonOut)# one differ..:
stopifnot(sum(cc$mcd.wt != ao.hbk$nonOut) <= 1)

## This is revealing: About 1--2 cases, where outliers are *not* == 1:14
## (2023: ~ 1/8 [sec] per call)
if(interactive()) {
  for(i in 1:30) {
    print(system.time(ao.hbk <- adjOutlyingness(hbk)))
    if(!identical(iout <- which(!ao.hbk$nonOut), 1:14)) {
	 cat("Outliers:\n"); print(iout)
    }
  }
}

## "Central" outlyingness: *not* calling mc()  anymore, since 2014-12-11:
trace(mc)
out <- capture.output(
  oo <- adjOutlyingness(hbk, clower=0, cupper=0)
)
untrace(mc)
stopifnot(length(out) == 0)

## A rank-deficient case
T <- tcrossprod(data.matrix(toxicity))
try(adjOutlyingness(T, maxit. = 20, trace.lev = 2)) # fails and recommends:
T. <- fullRank(T)
aT <- adjOutlyingness(T.)
plot(sort(aT$adjout, decreasing=TRUE), log="y")
plot(T.[,9:10], col = (1:2)[1 + (aT$adjout > 10000)])
## .. (not conclusive; directions are random, more 'ndir' makes a difference!)
}
\keyword{robust}
\keyword{multivariate}
