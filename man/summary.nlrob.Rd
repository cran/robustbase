\name{summary.nlrob}
\alias{summary.nlrob}

\title{Summarizing Robust Fits of Nonlinear Regression Models }
\description{
  \code{summary} method for class \code{"nlrob"}.
}
\usage{
\method{summary}{nlrob}(object, ...)
}

\arguments{
  \item{object}{an object of class \code{"nlrob"}, usually, a result of
    a call to \code{\link{nlrob}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

% \details{
%   \code{\link{summary.nlrob}} tries to be smart about formatting the
%   coefficients, standard errors, etc. and additionally gives
%   \dQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
% }

\value{
  The current implementation is trival and returns the same object.
  This will change in the future.
%   The function \code{\link{summary.nlrob}} computes and returns a list
%   of summary statistics of the robustly fitted linear model given in
%   \code{object}. The following elements are in the list:
%   \item{FixIt}{???}
}

%%\references{ ~put references to the literature/web site here ~ }
\author{Andreas Ruckstuhl}

\examples{
## To do
}
\keyword{regression}
\keyword{nonlinear}
\keyword{robust}
