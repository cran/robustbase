\name{adjboxStats}
\alias{adjboxStats}
\title{Statistics for Skewness-adjusted Boxplots}
\description{
  Computes the \dQuote{statistics} for producing boxplots adjusted for
  skewed distributions as proposed in Hubert and Vandervieren (2004).
}
\usage{
adjboxStats(x, coef = 1.5, a = -4, b = 3, do.conf = TRUE, do.out = TRUE)
}
\arguments{
  \item{x}{a numeric vector for which adjusted boxplot statistics are computed.}
  \item{coef}{number determining how far \sQuote{whiskers} extend out
    from the box, see \code{\link{boxplot.stats}}.}
  \item{a, b}{scaling factors multiplied by the medcouple
    \code{\link{mc}()} to determine outlyer boundaries; see the references.}
  \item{do.conf,do.out}{logicals; if \code{FALSE}, the \code{conf} or
    \code{out} component respectively will be empty in the result.}
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{list}} with the components

  \item{stats}{a matrix, each column contains the extreme of the lower
    whisker, the lower hinge, the median, the upper hinge and the extreme of
    the upper whisker for one group/plot. If all the inputs have the same
    class attribute, so will this component.}
  \item{n}{a vector with the number of observations in each group.}
  \item{coef}{a matrix where each column contains the lower and upper
    extremes of the notch.}
  \item{out}{the values of any data points which lie beyond the extremes
    of the whiskers.}
}
\author{ R Core Development Team, slightly adapted by Tobias Verbeke }
\note{The code only slightly modifies the code of \R's
  \code{\link{boxplot.stats}}.
}
\seealso{\code{\link{adjbox}()}, also for references,
  the function which mainly uses this one;
  further \code{\link{boxplot.stats}}.
}
\examples{
data(condroz)
astat <- adjboxStats(condroz[,"Ca"])
astat
}
\keyword{robust}
\keyword{univar}
