\name{lmrob..M..fit}
\alias{lmrob..M..fit}
\title{Compute M-estimators of regression}
\description{
  This function performs RWLS iterations to find an
  M-estimator of regression with Tukey's bisquare psi-function.  When
  started from an S-estimated \code{beta.initial}, this results in an
  MM-estimator.
}
\usage{
lmrob..M..fit(x, y, beta.initial, scale, c.psi, control)
}
\arguments{
  \item{x}{design matrix (\eqn{n \times p}{n x p}) typically including a
    column of \code{1}s for the intercept.}
  \item{y}{numeric response vector (of length \eqn{n}).}
  \item{beta.initial}{numeric vector (of length \eqn{p}) of initial
    estimate.  Usually the result of an S-regression estimator.}
  \item{scale}{robust residual scale estimate. Usually
    an S-scale estimator.}
  \item{c.psi}{the tuning constant of the \eqn{\psi /}{} psi-function of
    the M-estimator used.}
  \item{control}{list of control parameters, as returned
    by \code{\link{lmrob.control}}.  Currently, only the components
    \code{c("max.it", "rel.tol","trace.lev")} are accessed.}
}
\details{
  This function is used by \code{\link{lmrob.fit.MM}} %% and *fit.M() %% << TODO
  and typically not to be used on its own.
}
\value{A list with the following elements:
  \item{coef}{the M-estimator (or MM-estim.) of regression}
  \item{control}{the \code{control} list input used}
  \item{scale}{ The residual scale estimate}
  \item{seed }{ The random number generator seed}
  \item{converged}{ \code{TRUE} if the RWLS iterations converged,
    \code{FALSE} otherwise}
}
\references{
  Yohai, 1987
}
\seealso{
  \code{\link{lmrob.fit.MM}}, \code{\link{lmrob}};
  for M-estimators with more flexibility, e.g., choice of rho/psi:
  \code{\link[MASS]{rlm}} from package \pkg{MASS}.
}
\author{Matias Salibian-Barrera and Martin Maechler}
\examples{
data(stackloss)
X <- model.matrix(stack.loss ~ . , data = stackloss)
y <- stack.loss
## Compute manual MM-estimate:
## 1) initial LTS:
m0 <- ltsReg(X[,-1], y)
## 2) M-estimate started from LTS:
m1 <- lmrob..M..fit(X, y, beta.initial = coef(m0),
                    scale = m0$scale, c.psi = 1.6,
                    control = lmrob.control())
cbind(m0$coef, m1$coef)
## the scale is kept fixed:
stopifnot(identical(unname(m0$scale), m1$scale))

##  robustness weights: are
r.s <- with(m1, resid/scale) # scaled residuals
m1.wts <- robustbase:::tukeyPsi1(r.s, cc = 1.6) / r.s
summarizeRobWeights(m1.wts)
##--> outliers 1,3,4,13,21
which(m0$lts.wt == 0) # 1,3,4,21 but not 13

}
\keyword{robust}
\keyword{regression}
