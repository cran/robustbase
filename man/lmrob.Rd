\name{lmrob}
\alias{lmrob}
\title{MM-Estimators for Linear Regression}
\description{
  Computes fast MM-estimators for linear (regression) models.
}
\usage{
lmrob(formula, data, subset, weights, na.action,
      model = TRUE, x = !control$compute.rd, y = FALSE, singular.ok = TRUE,
      contrasts = NULL, offset = NULL, control = lmrob.control(...), ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit.  See
    \code{\link{lm}} and \code{\link{formula}} for more details.}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{lmrob} is called.}

  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used
    in the fitting process. %%% If specified, weighted least squares is used
    %%% with weights \code{weights} (that is, minimizing \code{sum(w*e^2)});
    %%% otherwise ordinary least squares is used.
    }
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}

  \item{model, x, y}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response) are returned.}

  \item{singular.ok}{logical. If \code{FALSE} (the default in S but
    not in \R) a singular fit is an error.}

  \item{contrasts}{an optional list.  See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.}

  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  An \code{\link{offset}} term can be included in the
    formula instead or as well, and if both are specified their sum is used.}

  \item{control}{a \code{\link{list}} specifying control parameters; use
    the function \code{\link{lmrob.control}(.)} and see its help page.}
  \item{\dots}{can be used to specify control parameters directly
    instead of via \code{control}.}
}
\details{
  This function computes an MM-regression estimator
  as described in Yohai (1987).  It uses a bi-square re-desceding
  score function, and by default it returns a highly
  robust and highly efficient estimator (with
  50\% breakdown point and 95\% asymptotic efficiency
  for normal errors).  The computation is carried out by a call to
  \code{\link{lmrob.fit.MM}()}.

  It uses an S-estimator (Rousseeuw
  and Yohai, 1984) for the errors which
  is also computed with a bi-square score function.
  The S-estimator is computed using the
  Fast-S algorithm of Salibian-Barrera and Yohai (2006), calling the
  function \code{\link{lmrob.S}}.

  Standard errors are computed using the
  formulas of Croux, Dhaene and Hoorelbeke (2003).
}
\value{
  An object of class \code{lmrob}. A list that includes the
  following components:
  \item{coefficients}{The MM-estimator of the coefficient vector}
  \item{initial.coefficients}{the S-estimator}
  \item{scale}{The scale as returned by S-estimator and used in the M one.}
  \item{cov}{The estimated covariance matrix of the regression coefficients}
  \item{residuals}{Residuals associated with the MM-estimator}
  \item{fitted.values}{Fitted values associated with the MM-estimator}
  \item{weights}{the \dQuote{robustness weights} \eqn{\psi(r_i/S) / (r_i/S)}.}
  \item{converged}{\code{TRUE} if the IRWLS iterations have converged}
}
\references{
  Croux, C., Dhaene, G. and Hoorelbeke, D. (2003)
  \emph{Robust standard errors for robust estimators},
  Discussion Papers Series 03.16, K.U. Leuven, CES.

  Rousseeuw, P.J. and Yohai, V.J. (1984)
  Robust regression by means of S-estimators,
  In \emph{Robust and Nonlinear Time Series},
  J. Franke, W. H\"ardle and R. D. Martin (eds.).
  Lectures Notes in Statistics 26, 256--272,
  Springer Verlag, New York.

  Salibian-Barrera, M. and Yohai, V.J. (2006)
  A fast algorithm for S-regression estimates,
  \emph{Journal of Computational and Graphical Statistics},
  in press. %%% <<< FIXME

  Yohai, V.J. (1987)
  High breakdown-point and high efficiency estimates for regression.
  \emph{The Annals of Statistics} \bold{15}, 642--65.
}
\author{ Matias Salibian-Barrera }
\seealso{
  \code{\link{lmrob.control}};
  for the algorithms \code{\link{lmrob.S}} and \code{\link{lmrob.fit.MM}};
  and for methods,
  \code{\link{summary.lmrob}},
  \code{\link{print.lmrob}}, and \code{\link{plot.lmrob}}.
}
\examples{
data(coleman)
summary( m1 <- lmrob(Y ~ ., data=coleman) )

data(starsCYG, package = "robustbase")
## Plot simple data and fitted lines
plot(starsCYG)
 lmST <-    lm(log.light ~ log.Te, data = starsCYG)
(RlmST <- lmrob(log.light ~ log.Te, data = starsCYG))
abline(lmST, col = "red")
abline(RlmST, col = "blue")
summary(RlmST)
vcov(RlmST)
stopifnot(all.equal(fitted(RlmST),
                    predict(RlmST, newdata = starsCYG),
                    tol = 1e-14))
}
\keyword{robust}
\keyword{regression}
