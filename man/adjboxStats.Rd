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
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{coef}{ ~~Describe \code{coef} here~~ }
  \item{a}{ ~~Describe \code{a} here~~ }
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{do.conf}{ ~~Describe \code{do.conf} here~~ }
  \item{do.out}{ ~~Describe \code{do.out} here~~ }
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
\note{The code slightly modifies the code of \code{\link{boxplot.stats}}.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
\keyword{univar}
