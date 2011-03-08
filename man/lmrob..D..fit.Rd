\name{lmrob..D..fit}
\alias{lmrob..D..fit}
\title{Compute Design Adaptive Scale estimate}
\description{This function calculates a Design Adaptive Scale estimate
  for a given MM-estimate. This is supposed to be a part of a chain of
  estimates like \code{SMD} or \code{SMDM}.
}
\usage{
lmrob..D..fit(obj, x=obj$x, control = obj$control)
}
\arguments{
  \item{obj}{\code{lmrob}-object based on which the estimate is to be
    calculated.}
  \item{x}{The design matrix, if missing the method tries to get it from
    \code{obj$x}and if this fails from \code{obj$model}.}
  \item{control}{list of control parameters, as returned
    by \code{\link{lmrob.control}}.}  
}
\details{
  This function is used by \code{\link{lmrob.fit}} and typically not to
  be used on its own.
}
\value{The given \code{lmrob}-object with the following elements updated:
  \item{scale}{The Design Adaptive Scale estimate}
  \item{converged}{ \code{TRUE} if the scale calculation converged,
    \code{FALSE} other.}
}
\references{
  Koller and Stahel, 2010 (tbp).
}
\seealso{
  \code{\link{lmrob.fit}}, \code{\link{lmrob}}
}
\examples{
data(stackloss)
## Compute manual SMD-estimate:
## 1) MM-estimate
m1 <- lmrob(stack.loss ~ ., data = stackloss)
## 2) Add Design Adaptive Scale estimate
m2 <- lmrob..D..fit(m1)
print(c(m1$scale, m2$scale))

summary(m1)
summary(m2) ## the covariance matrix estimate is also updated
}
\author{Manuel Koller}
\keyword{robust}
\keyword{regression}
