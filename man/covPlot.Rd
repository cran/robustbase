\name{plot.mcd}
\alias{covPlot}
\alias{plot.mcd}
%NO \alias{ddplot}
%NO \alias{distplot}
%NO \alias{chi2qqplot}
%NO \alias{ellipse}
\title{Robust Distance Plots}
\description{
  Shows the Mahalanobis distances based on robust and classical estimates
  of the location and the covariance matrix in different plots.
  The following plots are available:
  \itemize{
    \item index plot of the robust and mahalanobis distances
    \item distance-distance plot
    \item Chisquare QQ-plot of the robust and mahalanobis distances
    \item plot of the tolerance ellipses (robust and classic)
    \item Scree plot - Eigenvalues comparison plot
  }
}
\usage{
\method{plot}{mcd}(x,
     which = c("all", "dd", "distance", "qqchi2", "tolEllipsePlot", "screeplot"),
     classic = FALSE, ask = (which=="all" && dev.interactive()),
     cutoff, id.n, tol = 1e-7, \dots)

%    covPlot(x, which = c("all","dd","distance","qqchi2","tolEllipsePlot","screeplot"), classic=FALSE, ask=(which=="all" && dev.interactive()), mcd, cutoff, id.n, tol = 1e-7)

%% ddplot(x, \dots)
%% distplot(x, \dots)
%% chi2qqplot(x, \dots)
%% ellipse(x, \dots)
}

\arguments{
  \item{x}{a \code{mcd} object, typically result of \code{covMcd}.}
  \item{which}{string indicating which plot to show.  See the
    \emph{Details} section for a description of the options.  Defaults to \code{"all"}.}.
  \item{classic}{whether to plot the classical distances too. Default is \code{FALSE}.}.
  \item{ask}{logical indicating if the user should be \emph{ask}ed
    before each plot, see \code{\link{par}(ask=.)}.  Defaults to
    \code{which == "all" && \link{dev.interactive}()}.
  }
  \item{cutoff}{the cutoff value for the distances.}
  \item{id.n}{number of observations to be identified by a label.  If
    not supplied, the number of observations with distance larger than
    \code{cutoff} is used.}
  \item{tol}{tolerance to be used for computing the inverse, see
    \code{\link{solve}}.  Defaults to \code{tol = 1e-7}.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
}
\details{
  This function produces several plots based on the robust and classical
  location and covariance matrix. Which of them to select is specified
  by the attribute  \code{which}.
  The possible options are:

 \code{distance} - index plot of the robust distances;

 \code{dd} - distance-distance plot;

 \code{qqchi2} - a qq-plot of the robust distances versus the quantiles of the chi-squared distribution

 \code{tolEllipsePlot} - a tolerance ellipse

 \code{screeplot} - an eigenvalues comparison plot - screeplot

 The Distance-Distance Plot, introduced by
 Rousseeuw and van Zomeren (1990), displays the robust distances
 versus the classical Mahalanobis distances. The dashed line is the set of
 points where the robust distance is equal to the classical distance.
 The horizontal and vertical lines are drawn at values equal to the cutoff
 which defaults to square root of the 97.5\% quantile of a chi-squared
 distribution with p degrees of freedom. Points beyond these lines can
 be considered outliers.
}
%\value{}
\references{
  P. J. Rousseeuw and van Zomeren, B. C. (1990).
  Unmasking Multivariate Outliers and Leverage Points.
  \emph{Journal of the American Statistical Association} \bold{85}, 633--639.

  P. J. Rousseeuw and K. van Driessen (1999)
  A fast algorithm for the minimum covariance determinant estimator.
  \emph{Technometrics} \bold{41}, 212--223.
}
\seealso{
  \code{\link{tolEllipsePlot}}
}
\examples{
data(Animals, package ="MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
mcd <- covMcd(log(brain))

plot(mcd, which = "distance", classic = TRUE)# 2 plots
plot(mcd, which = "dd")
plot(mcd, which = "tolEllipsePlot", classic = TRUE)
op <- par(mfrow = c(2,3))
plot(mcd) ## -> which = "all" (5 plots)
par(op)
\dontshow{
op <- par(mfrow = c(3,3))
plot(mcd, classic = TRUE) ## -> which = "all" (7 plots)
par(op)
}
}
\keyword{hplot}
\keyword{robust}
\keyword{multivariate}
