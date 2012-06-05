\name{lmrob}
\encoding{utf8}
\alias{lmrob}
% "Link to here", even those are not exported:
\alias{.vcov.avar1}
\alias{.vcov.w}
\title{MM-type Estimators for Linear Regression}
\description{
  Computes fast MM-type estimators for linear (regression) models.
}
\usage{
lmrob(formula, data, subset, weights, na.action, method = "MM",
      model = TRUE, x = !control$compute.rd, y = FALSE,
      singular.ok = TRUE, contrasts = NULL, offset = NULL,
      control = NULL, init = NULL, ...)
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

  \item{method}{string specifying the estimator-chain. \code{MM}
    is interpreted as \code{SM}. See \emph{Details}.}

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

  \item{init}{an optional argument to specify or supply the initial
  estimate. See \emph{Details}.}

  \item{\dots}{can be used to specify control parameters directly
    instead of via \code{control}.}
}
\details{
  \describe{
    \item{Overview}{
      This function computes an MM-type regression estimator as
      described in Yohai (1987) and Koller and Stahel (2011).  By
      default it uses a bi-square redescending score function, and it
      returns a highly robust and highly efficient estimator (with 50\%
      breakdown point and 95\% asymptotic efficiency for normal errors).
      The computation is carried out by a call to \code{\link{lmrob.fit}()}.

      The argument \code{setting} of \code{\link{lmrob.control}} is provided
      to set alternative defaults as suggested in Koller and Stahel (2011)
      (use \code{setting='KS2011'}).  For details, see
      \code{\link{lmrob.control}}.
    }
    \item{Initial Estimator \code{init}:}{
      The initial estimator may be specified using the argument
      \code{init}. This argument can either be a string, a function or a
      list.  A \emph{string} can be used to specify built in internal
      estimators (currently \code{S} and \code{M-S}, see \emph{See also}
      below).  A \emph{function} taking arguments \code{x, y, control,
      mf} (where \code{mf} stands for \code{model.frame}) and returning
      a list containing at least the initial coefficients as
      \code{coefficients} and the initial scale estimate
      \code{scale}. Or a \emph{list} giving the initial coefficients and
      scale as \code{coefficients} and \code{scale}. See also \emph{Examples}.

      Note that if the init argument supplied is a function or list, the
      method argument must \emph{not} contain the initial estimator, e.g.,
      use \code{MDM} instead of \code{SMDM}.

      The default, equivalent to \code{init = "S"}, uses as initial
      estimator an S-estimator (Rousseeuw and Yohai, 1984) which is
      computed using the Fast-S algorithm of Salibian-Barrera and Yohai
      (2006), calling \code{\link{lmrob.S}()}.  That function, since
      March 2012, uses \emph{nonsingular} subsampling which makes the
      Fast-S algorithm feasible for categorical data as well.
    }
    \item{Method \code{method}:}{
      The following chain of estimates is customizable via the
      \code{method} argument. % of \code{\link{lmrob.control}}.
      There are currently two types of estimates available,
      \describe{
	\item{\code{"M"}:}{corresponds to the standard M-regression
	  estimate.}
	\item{\code{"D"}:}{stands for the Design Adaptive Scale estimate
	  as proposed in Koller and Stahel (2011).}
      }
      The \code{method} argument takes a string that specifies the
      estimates to be calculated as a chain.  Setting
      \code{method='SMDM'} will result in an intial S-estimate, followed
      by an M-estimate, a Design Adaptive Scale estimate and a final
      M-step.  For methods involving a \code{D}-step, the default value
      of \code{psi} (see \code{\link{lmrob.control}}) is changed to
      \code{"lqq"}.

      By default, standard errors are computed using the formulas of
      Croux, Dhaene and Hoorelbeke (2003) (\code{\link{lmrob.control}}
      option \code{cov=".vcov.avar1"}).  This method, however, works only
      for MM-estimates. For other \code{method} arguments, the
      covariance matrix estimate used is based on the asymptotic
      normality of the estimated coefficients (\code{cov=".vcov.w"}) as
      described in Koller and Stahel (2011).
    }
  }%end {describe}
}
\value{
  An object of class \code{lmrob}. A list that includes the
  following components:
  \item{coefficients}{The estimate of the coefficient vector}
  \item{init.S}{The list returned by \code{\link{lmrob.S}} or
  \code{\link{lmrob.M.S}} (for MM-estimates only)}
  \item{init}{A similar list that contains the results of intermediate
  estimates (not for MM-estimates).}
  \item{scale}{The scale as used in the M estimator.}
  \item{cov}{The estimated covariance matrix of the regression coefficients}
  \item{residuals}{Residuals associated with the estimator}
  \item{fitted.values}{Fitted values associated with the estimator}
  \item{weights}{the \dQuote{robustness weights} \eqn{\psi(r_i/S) / (r_i/S)}.}
  \item{converged}{\code{TRUE} if the IRWLS iterations have converged}
}
\references{
  Croux, C., Dhaene, G. and Hoorelbeke, D. (2003)
  \emph{Robust standard errors for robust estimators},
  Discussion Papers Series 03.16, K.U. Leuven, CES.

  Koller, M. and Stahel, W.A. (2011), Sharpening Wald-type inference in
  robust regression for small samples, \emph{Computational Statistics &
    Data Analysis} \bold{55}(8), 2504--2515.

  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.

  Rousseeuw, P.J. and Yohai, V.J. (1984)
  Robust regression by means of S-estimators,
  In \emph{Robust and Nonlinear Time Series},
  J. Franke, W. Härdle and R. D. Martin (eds.).
  Lectures Notes in Statistics 26, 256--272,
  Springer Verlag, New York.

  Salibian-Barrera, M. and Yohai, V.J. (2006)
  A fast algorithm for S-regression estimates,
  \emph{Journal of Computational and Graphical Statistics},
  \bold{15}(2), 414--427.

  Yohai, V.J. (1987)
  High breakdown-point and high efficiency estimates for regression.
  \emph{The Annals of Statistics} \bold{15}, 642--65.
}
\author{ Matias Salibian-Barrera and Manuel Koller}
\seealso{
  \code{\link{lmrob.control}};
  for the algorithms \code{\link{lmrob.S}}, \code{\link{lmrob.M.S}} and
  \code{\link{lmrob.fit}};
  and for methods,
  \code{\link{predict.lmrob}}, \code{\link{summary.lmrob}},
  \code{\link{print.lmrob}}, and \code{\link{plot.lmrob}}.
}
\examples{
data(coleman)
set.seed(0)
summary( m1 <- lmrob(Y ~ ., data=coleman) )
summary( m2 <- lmrob(Y ~ ., data=coleman, setting = 'KS2011') )

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

## --- init argument
## string
set.seed(0)
m3 <- lmrob(Y ~ ., data=coleman, init = "S")
stopifnot(all.equal(m1[-17], m3[-17]))
## function
initFun <- function(x, y, control, mf) {
    init.S <- lmrob.S(x, y, control)
    list(coefficients=init.S$coef, scale = init.S$scale)
}
set.seed(0)
m4 <- lmrob(Y ~ ., data=coleman, method = "M", init = initFun)
## list
m5 <- lmrob(Y ~ ., data=coleman, method = "M",
            init = list(coefficients = m3$init$coef, scale = m3$scale))
stopifnot(all.equal(m4[-17], m5[-17]))

}
\keyword{robust}
\keyword{regression}
