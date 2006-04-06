\name{lmrob.MM}
\alias{lmrob.MM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to compute MM-regression estimators }
\description{
  This function performs RWLS iterations to find an
  MM-regression estimator
}
\usage{
lmrob.MM(x, y, beta.initial, scale, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Design matrix }
  \item{y}{ Response vector }
  \item{beta.initial}{ Vector of initial estimators. Usually
  an S-regression estimator}
  \item{scale}{ A robust residual scale estimate. Usually
  an S-scale estimator.}
  \item{control}{ A list of control parameters as returned
  by \code{lmrob.control} }
}
\details{
  This function is used by \code{lmrob.fit.MM} and not intended
  to be used on its own.
}
\value{
	A list with the following elements:
  \item{coef }{ The MM-regression estimator }
  \item{cov  }{ The covariance matrix of the MM-regression estimator} 
  \item{control}{ The control list used}
  \item{scale }{ The residual scale estimate}
  \item{seed  }{ The random number generator seed}
  \item{converged}{ \code{TRUE} if the RWLS iterations converged,
  \code{FALSE} otherwise}
}
\references{ Yohai, 1987 }
\author{ Matias Salibian-Barrera }
\keyword{ robust }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
