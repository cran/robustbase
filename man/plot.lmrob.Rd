\name{plot.lmrob}
\alias{plot.lmrob}
\title{Plot Method for "lmrob" Objects}
\description{
  Diagnostic plots for elements of class lmrob
}
\usage{
\method{plot}{lmrob}(x, which = 1:5,
     caption = c("Standardized residuals vs. Robust Distances",
   "Normal Q-Q vs. Residuals", "Response vs. Fitted Values",
   "Residuals vs. Fitted Values" , "Sqrt of abs(Residuals) vs. Fitted Values"),
   panel = points, sub.caption = deparse(x$call), main = "",
   ask = prod(par("mfcol")) < length(which) && dev.interactive(),
   \dots, p=0.025)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object as created by \code{lmrob} }
  \item{which}{ integer number between 1 and 5 to specify which
    plot is desired }
  \item{caption}{Caption for the different plots}
  \item{panel}{Panel}
  \item{main}{Main title}
  \item{sub.caption}{sub titles}
  \item{ask}{waits for user input before displaying each plot }
  \item{...}{optional arguments for \code{par}}
  \item{p}{threshold for distance-distance plot}
}
\references{ Robust diagnostic plots as in Rousseuw and van Zomeren
(1990) }

\seealso{ \code{\link{lmrob}}, also for examples, \code{\link{plot.lm}}.
}
\examples{
data(starsCYG)
## Plot simple data and fitted lines
plot(starsCYG)
 lmST <-    lm(log.light ~ log.Te, data = starsCYG)
RlmST <- lmrob(log.light ~ log.Te, data = starsCYG)
RlmST
abline(lmST, col = "red")
abline(RlmST, col = "blue")

op <- par(mfrow = c(2,2), mgp = c(1.5, 0.6, 0), mar= .1+c(3,3,3,1))
plot(RlmST, which = c(1:2, 4:5))
par(op)
}
\keyword{robust}
\keyword{regression}

