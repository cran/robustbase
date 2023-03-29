\name{lmc}
\alias{lmc}
\alias{rmc}
\title{Left and Right Medcouple, Robust Measures of Tail Weight}
\description{
  Compute the left and right \sQuote{medcouple}, \emph{robust} estimators
  of tail weight, in some sense robust versions of the kurtosis, the very
  unrobust centralized 4th moment.
}
\usage{
lmc(x, mx = median(x, na.rm=na.rm), na.rm = FALSE, doReflect = FALSE, ...)
rmc(x, mx = median(x, na.rm=na.rm), na.rm = FALSE, doReflect = FALSE, ...)
}
\arguments{
  \item{x}{a numeric vector}
  \item{mx}{number, the \dQuote{center} of \code{x} wrt which the left and
    right parts of \code{x} are defined: \preformatted{
    lmc(x, mx, *) :=  mc(x[x <= mx], *)
    rmc(x, mx, *) :=  mc(x[x >= mx], *)}
  }
  \item{na.rm}{logical indicating how missing values (\code{\link{NA}}s)
    should be dealt with.}
  \item{doReflect}{logical indicating if \code{\link{mc}} should also be
    computed on the \emph{reflected} sample \code{-x}.  Setting
    \code{doReflect=TRUE} makes sense for mathematical strictness reasons,
    as the internal MC computes the himedian() which can differ slightly from
    the median.
    Note that \code{\link{mc}()}'s own default is true iff \code{length(x) <= 100}.
  }
  \item{\dots}{further arguments to \code{\link{mc}()}, see its help page.}
}
\value{
  each a number (unless \code{\dots} contains \code{full.result = TRUE}).
}
\references{
  Brys, G., Hubert, M. and Struyf, A. (2006).
  Robust measures of tail weight,
  \emph{Computational Statistics and Data Analysis} \bold{50(3)}, 733--759.

  and those in \sQuote{References} of \code{\link{mc}}.
}
\examples{
mc(1:5)  # 0 for a symmetric sample
lmc(1:5) # 0
rmc(1:5) # 0

x1 <- c(1, 2, 7, 9, 10)
mc(x1) # = -1/3
c( lmc( x1),  lmc( x1, doReflect=TRUE))#   0  -1/3
c( rmc( x1),  rmc( x1, doReflect=TRUE))# -1/3 -1/6
c(-rmc(-x1), -rmc(-x1, doReflect=TRUE)) # 2/3  1/3

data(cushny)
lmc(cushny) # 0.2
rmc(cushny) # 0.45

isSym_LRmc <- function(x, tol = 1e-14)
    all.equal(lmc(-x, doReflect=TRUE),
              rmc( x, doReflect=TRUE), tolerance = tol)

sym <- c(-20, -5, -2:2, 5, 20)
stopifnot(exprs = {
    lmc(sym) == 0.5
    rmc(sym) == 0.5
    isSym_LRmc(cushny)
    isSym_LRmc(x1)
})

## Susceptibility to large outliers:
## "Sensitivity Curve" := empirical influence function
dX10 <- function(X) c(1:5,7,10,15,25, X) # generate skewed size-10 with 'X'
x <- c(26:40, 45, 50, 60, 75, 100)
(lmc10N <- vapply(x, function(X) lmc(dX10(X)), 1))
(rmc10N <- vapply(x, function(X) rmc(dX10(X)), 1))
cols <- adjustcolor(2:3, 3/4)

plot(x, lmc10N, type="o", cex=1/2, main = "lmc & rmc( c(1:5,7,10,15,25, X) )",
     xlab=quote(X), log="x", col=cols[1])
lines(x, rmc10N, col=cols[2], lwd=3)
legend("top", paste0(c("lmc", "rmc"), "(X)"), col=cols, lty=1, lwd=c(1,3), pch = c(1, NA), bty="n")

n <- length(x)
stopifnot(exprs = {
    all.equal(current = lmc10N, target = rep(0, n))
    all.equal(current = rmc10N, target = c(3/19, 1/5, 5/21, 3/11, 7/23, rep(1/3, n-5)))
    ## and it stays stable with outlier  X --> oo :
    lmc(dX10(1e300)) == 0
    rmc(dX10(1e300)) == rmc10N[6]
})
}
\keyword{robust}
\keyword{univar}
