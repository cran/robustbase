\name{summary.mcd}
\alias{summary.mcd}
\title{Summary Method for MCD objects}
\usage{
\method{summary}{mcd}(object, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{object}{an object of class \code{"mcd"}, usually, a result of a call to \code{\link{covMcd}}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  \code{summary} method for class \code{"mcd"}.
}
\details{

%%MM: FIXME:  summary() should return summary object

  This function prints summary statistics for the weighted covariance
  matrix and location estimates with weights based on MCD estimates.
  While the function \code{print.mcd} prints only the robust estimates
  of the location and the covariance matrix, \code{summary.mcd} will
  print also the correlation matrix (if requested in the call to
  \code{covMcd} with \code{cor=TRUE}), the eigenvalues of the covariance
  or the correlation matrix and the robust distances.
}
\value{
  The \code{mcd} object is returned.
}
\seealso{
  \code{\link{summary}}
}
\examples{
data(Animals, package = "MASS")
brain <- Animals[c(1:24, 26:25, 27:28),]
lbrain <- log(brain)
summary(cLB <- covMcd(lbrain))
}
\keyword{multivariate}
\keyword{robust}
