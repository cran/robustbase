\name{tolEllipsePlot}
\alias{tolEllipsePlot}
\title{Tolerance Ellipse Plot}
\description{
  Plots the 0.975 tolerance ellipse of the bivariate data set \code{x}.
  The ellipse is defined by those data points whose distance
  is equal to the squareroot of the 0.975 chisquare quantile
  with 2 degrees of freedom.
}
\usage{
tolEllipsePlot(x, mcd = NULL, cutoff = NULL, id.n = NULL,
	       classic = FALSE, tol = 1e-07,
	       xlab = "", ylab = "",
               main = "Tolerance ellipse (97.5\%)",
	       txt.leg = c("robust", "classical"),
	       col.leg = c("red", "blue"),
	       lty.leg = c("solid","dashed"))
}

\arguments{
  \item{x}{a two dimensional matrix or data frame. }
  \item{mcd}{an object similar to those of class \code{"mcd"}, only its
    components \code{center} and \code{cov} will be used.  If missing, the
    MCD will be computed (via \code{\link{covMcd}()}).}
  \item{cutoff}{ Distance needed to flag data points outside the ellipse  }
  \item{id.n}{number of observations to be identified by a label.  If
    not supplied, the number of observations with distance larger than
    \code{cutoff} is used.}
  \item{classic}{whether to plot the classical distances too. Default is \code{FALSE}.}.
  \item{tol}{tolerance to be used for computing the inverse, see
    \code{\link{solve}}.  Defaults to \code{tol = 1e-7}.}
  \item{xlab, ylab, main}{passed to \code{\link{plot.default}}.}
  \item{txt.leg, col.leg, lty.leg}{character vectors of length 2 for the
    legend, only if used if \code{class = TRUE}.}
}

\author{Peter Filzmoser, Valentin Todorov and Martin Maechler}
%\details{}
%\value{}
%\references{ }
\seealso{
  \code{\link[cluster]{ellipsoidhull}} and
  \code{\link[cluster]{predict.ellipsoid}} from package \pkg{cluster}.
}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
mcd <- covMcd(hbk.x)       # compute mcd in advance
## must be a 2-dimensional data set: take the first two columns :
tolEllipsePlot(hbk.x[,1:2])
}
\keyword{hplot}
\keyword{robust}

