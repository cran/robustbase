\name{mc}
\alias{mc}
\title{Medcouple, a Robust Measure of Skewness}
\description{
  Compute the \sQuote{medcouple}, a \emph{robust} concept and estimator
  of skewness.  The medcouple is defined as a scaled median difference
  of the left and right half of distribution, and hence \emph{not} based
  on the third moment as the classical skewness.
}
\usage{
mc(x, na.rm = FALSE, doReflect = (length(x) <= 100),
   doScale = FALSE,     # was hardwired=TRUE, then default=TRUE
   c.huberize = 1e11,   # was implicitly = Inf originally
   eps1 = 1e-14, eps2 = 1e-15,   # << new in 0.93-2 (2018-07..)
   maxit = 100, trace.lev = 0, full.result = FALSE)
}
\arguments{
  \item{x}{a numeric vector}
  \item{na.rm}{logical indicating how missing values (\code{\link{NA}}s)
    should be dealt with.}
  \item{doReflect}{logical indicating if the internal MC should also be
    computed on the \emph{reflected} sample \code{-x}, with final result
    \code{(mc.(x) - mc.(-x))/2}.  This makes sense since the internal
    MC, \code{mc.()} computes the himedian() which can differ slightly from
    the median.}%% only whenever sum(x <= med) * sum(x >= med) is even
  \item{doScale}{logical indicating if the internal algorithm should
    also \emph{scale} the data (using the most distant value from the
    median which is unrobust and numerically dangerous); scaling has been
    hardwired in the original algorithm and \R's \code{mc()} till summer
    2018, where it became the default.  Since \pkg{robustbase} version 0.95-0,
    March 2022, the default is \code{FALSE}.  As this may change the
    result, a message is printed about the new default, once per \R
    session.  You can suppress the message by specifying \code{doScale = *}
    explicitly, or, by setting \code{\link{options}(mc_doScale_quiet=TRUE)}.}
  \item{c.huberize}{a positive number (default: \code{1e11}) used to
    stabilize the sample via \code{x <- \link{huberize}(x, c = c.huberize)}
    for the \code{mc()} computations in the case of a nearly degenerate
    sample (many observations practically equal to the median) or very
    extreme outliers.  In previous versions of \pkg{robustbase} no such
    huberization was applied which is equivalent to \code{c.huberize = Inf}.}
  \item{eps1, eps2}{tolerance in the algorithm; \code{eps1} is used as a  for
    convergence tolerance, where \code{eps2} is only used in the internal
    \code{h_kern()} function to prevent underflow to zero, so could be
    considerably smaller.  The original code implicitly \emph{hard
      coded} in C \code{eps1 := eps2 := 1e-13};  only change with care!}
  \item{maxit}{maximal number of iterations; typically a few should be
    sufficient.}
  \item{trace.lev}{integer specifying how much diagnostic output the
    algorithm (in C) should produce.  No output by default, most output
    for \code{trace.lev = 5}.}
  \item{full.result}{logical indicating if the full return values (from
    C) should be returned as a list via \code{attr(*, "mcComp")}.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  a number between -1 and 1, which is the medcouple, \eqn{MC(x)}.
  For \code{r <- mc(x, full.result = TRUE, ....)}, then
  \code{attr(r, "mcComp")} is a list with components
  \item{medc}{the medcouple  \eqn{mc.(x)}.}
  \item{medc2}{the medcouple \eqn{mc.(-x)} if \code{doReflect=TRUE}.}
  \item{eps}{tolerances used.}
  \item{iter,iter2}{number of iterations used.}
  \item{converged,converged2}{logical specifying \dQuote{convergence}.}
}
\section{Convergence Problems}{
  For extreme cases there were convergence problems which should not
  happen anymore as we now use \code{doScale=FALSE} and huberization (when
  \code{c.huberize < Inf}).

  %% Some of them can be alleviated by \dQuote{loosening} the tolerances
  %% \code{eps1} and \code{eps2}.
  %% \cr

  The original algorithm and \code{mc(*, doScale=TRUE)} not only centers
  the data around the median but
  also scales them by the extremes which may have a negative effect
  e.g., when changing an extreme outlier to even more extreme, the
  result changes wrongly; see the 'mc10x' example.
}
\references{
  Guy Brys, Mia Hubert and Anja Struyf (2004)
  A Robust Measure of Skewness;
  \emph{JCGS} \bold{13} (4), 996--1017.

  Hubert, M. and Vandervieren, E. (2008).
  An adjusted boxplot for skewed distributions,
  \emph{Computational Statistics and Data Analysis} \bold{52}, 5186--5201.

  Lukas Graz (2021). Improvement of the Algorithms for the Medcoule and the
  Adjusted Outlyingness; unpublished BSc thesis, supervised by M.Maechler, ETH Zurich.
}
\author{Guy Brys; modifications by Tobias Verbeke and bug fixes and
  extensions by Manuel Koller and Martin Maechler.

  The new default \code{doScale=FALSE}, and the new \code{c.huberize} were
  introduced as consequence of Lukas Graz' BSc thesis.
}
\seealso{\code{\link{Qn}} for a robust measure of scale (aka
  \dQuote{dispersion}), ....
}
\examples{
mc(1:5) # 0 for a symmetric sample

x1 <- c(1, 2, 7, 9, 10)
mc(x1) # = -1/3

data(cushny)
mc(cushny) # 0.125

stopifnot(mc(c(-20, -5, -2:2, 5, 20)) == 0,
          mc(x1, doReflect=FALSE) ==  -mc(-x1, doReflect=FALSE),
          all.equal(mc(x1, doReflect=FALSE), -1/3, tolerance = 1e-12))

## Susceptibility of the current algorithm to large outliers :
dX10 <- function(X) c(1:5,7,10,15,25, X) # generate skewed size-10 with 'X'
x <- c(10,20,30, 100^(1:20))
## (doScale=TRUE, c.huberize=Inf)  were (implicit) defaults in earlier {robustbase}:
(mc10x <- vapply(x, function(X) mc(dX10(X), doScale=TRUE, c.huberize=Inf), 1))
## limit X -> Inf  should be 7/12 = 0.58333...  but that "breaks down a bit" :
plot(x, mc10x, type="b", main = "mc( c(1:5,7,10,15,25, X) )", xlab="X", log="x")
## The new behavior is much preferable {shows message about new 'doScale=FALSE'}:
(mc10N <- vapply(x, function(X) mc(dX10(X)), 1))
lines(x, mc10N, col=adjustcolor(2, 3/4), lwd=3)
mtext("mc(*, c.huberize=1e11)",  col=2)
stopifnot(all.equal(c(4, 6, rep(7, length(x)-2))/12, mc10N))
## Here, huberization already solves the issue:
mc10NS <- vapply(x, function(X) mc(dX10(X), doScale=TRUE), 1)
stopifnot(all.equal(mc10N, mc10NS))
}
\keyword{robust}
\keyword{univar}
