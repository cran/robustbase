\name{glmrob}
\alias{glmrob}

\title{Robust Fitting of Generalized Linear Models}

\description{
  \code{glmrob} is used to fit generalized linear models by robust
  methods.  The models are specified by giving a symbolic description of
  the linear predictor and a description of the error distribution.
  Currently, robust methods are implemented only for discrete response
  distributions, i.e. \code{\link{family} = binomial} or \code{ = poisson}.
}

\usage{
glmrob(formula, family, data, weights, subset, na.action,
       start = NULL, offset, method = "Mqle",
       weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL,
       model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, trace = FALSE, ...)
}

\arguments{
  \item{formula}{a \code{\link{formula}}, i.e., a symbolic description
    of the model to be fit (cf. \code{\link{glm}} or \code{\link{lm}}).}

  \item{family}{a description of the error distribution and link function to
    be used in the model. This can be a character string naming a
    family function, a family function or the result of a call to
    a family function.  (See \code{\link{family}} for details of
    family functions.)  As mentioned, currently this must be
    \code{binomial} or \code{poisson}.}

  \item{data}{an optional data frame containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from
    which \code{glmrob} is called.}

  \item{weights}{an optional vector of weights to be used in the fitting
    process.}

  \item{subset}{an optional vector specifying a subset of observations to be
    used in the fitting process.}

  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting in \code{\link{options}}. The
    \dQuote{factory-fresh} default is \code{\link{na.omit}}.}

  \item{start}{starting values for the parameters in the linear predictor.}

  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.}

  \item{method}{a character string specifying the robust fitting method.
    The details of method specification are given below. }

  \item{weights.on.x}{
    character string (can be abbreviated) specifying how points
    (potential outliers) in x-space are downweighted.  If \code{"hat"},
    weights on the design of the form \eqn{\sqrt{1-h_{ii}}} are used,
    where \eqn{h_{ii}} are the diagonal elements of the hat matrix.  If
    \code{"robCov"}, weights based on the robust Mahalanobis distance of the
    design matrix (intercept excluded) are used where the covariance
    matrix and the centre is estimated by \code{\link[MASS]{cov.rob}}
    from the package \pkg{MASS}.\cr
    Similarly,  if \code{"covMcd"}, robust weights are computed using
    \code{\link{covMcd}}.  The default is \code{"none"}.}

  \item{control}{a list of parameters for controlling the fitting process.
    See the documentation for \code{\link{glmrobMqle.control}} for
    details.}

  \item{model}{a logical value indicating whether \emph{model frame}
    should be included as a component of the returned value.}

  \item{x, y}{logical values indicating whether the response
    vector and model matrix used in the fitting process should be
    returned as components of the returned value.}

  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{model.matrix.default}.}
  \item{trace}{logical (or integer) indicating if intermediate results
    should be printed; defaults to \code{FALSE}.}
  \item{\dots}{arguments passed to \code{\link{glmrobMqle.control}} when
    \code{control} is \code{NULL} (as per default).}
}

\details{
  \code{method="Mqle"} fits a generalized linear
  model using Mallows or Huber type robust estimators, as described in
  Cantoni and Ronchetti (2001).  In contrast to the implementation
  described in Cantoni (2004), the pure influence algorithm is
  implemented.
  \cr
  Currently no other method is implemented.

  \code{weights.on.x= "robCov"} makes sense if all explanatory variables
  are continuous.
}

\value{
  \code{glmrob} returns an object of class \code{"glmrob"} and is also
  inheriting from \code{\link{glm}}.
  \cr
  The function \code{\link{summary}} (i.e., \code{\link{summary.glmrob}}) can
  be used to obtain or print a summary of the results.
  \cr
  The generic accessor functions \code{\link{coefficients}},
  \code{effects}, \code{fitted.values} and \code{residuals} can be used to
  extract various useful features of the value returned by \code{glmrob()}.

  An object of class \code{"glmrob"} is a list with at least the
  following components:
  \item{coefficients}{a named vector of coefficients}
  \item{residuals}{the \emph{Pearson} residuals}
  \item{fitted.values}{the fitted mean values, obtained by transforming
    the linear predictors by the inverse of the link function.}
  \item{w.r}{robustness weights for each observations; i.e.,
    \code{residuals} \eqn{\times}{*} \code{w.r} equals the psi-function of the
    Preason's residuals.}
  \item{w.x}{weights used to down-weight observations based on the
    position of the observation in the design space.}
  \item{cov}{the estimated asymptotic covariance matrix of the estimated
    coefficients.}
  \item{tcc}{the tuning constant c in Huber's psi-function.}
  \item{family}{the \code{\link{family}} object used.}
  \item{linear.predictors}{the linear fit on link scale.}
  \item{deviance}{up to a constant, minus twice the maximized
    log-likelihood.  Where sensible, the constant is chosen so that a
    saturated model has deviance zero.}
  \item{null.deviance}{The deviance for the null model, comparable with
    \code{deviance}. The null model will include the offset, and an
    intercept if there is one in the model}
  \item{iter}{the number of iterations used by the influence algorithm.}
  \item{converged}{logical. Was the IWLS algorithm judged to have converged?}
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{terms}{the \code{\link{terms}} object used.}
  \item{data}{the \code{data argument}.}
  \item{offset}{the offset vector used.}
  \item{control}{the value of the \code{control} argument used.}
  \item{method}{the name of the robust fitter function used.}
  \item{contrasts}{(where relevant) the contrasts used.}
  \item{xlevels}{(where relevant) a record of the levels of the factors
    used in fitting.}

  If a \code{\link{binomial}} is specified by giving a two-column
  response ???
%%   ??? FIXME
%%   If a \code{\link{binomial}} \code{glm} model is specified by giving a
%%   two-column response, the weights returned by \code{prior.weights} are
%%   the total numbers of cases (multipied by the supplied case weights) and
%%   the component \code{y} of the result is the proportion of successes.
}

\references{
  E. Cantoni and E. Ronchetti (2001)
  Robust Inference for Generalized Linear Models.
  \emph{JASA} \bold{96} (455), 1022--1030.

  E.Cantoni (2004)
  Analysis of Robust Quasi-deviances for Generalized Linear Models.
  \emph{Journal of Statistical Software}, \bold{10},
  \url{http://www.jstatsoft.org}
}

\author{Andreas Ruckstuhl}

%%\note{ }

\seealso{\code{\link{glmrobMqle.control}}}
\examples{
## Binomial response --------------
data(carrots)

Cfit1 <- glm(cbind(success, total-success) ~ logdose + block,
             data = carrots, family = binomial)
summary(Cfit1)

Rfit1 <- glmrob(cbind(success, total-success) ~ logdose + block,
                family = binomial, data = carrots, method= "Mqle",
                control= glmrobMqle.control(tcc=1.2))
summary(Rfit1)

Rfit2 <- glmrob(success/total ~ logdose + block, weights = total,
                family = binomial, data = carrots, method= "Mqle",
                control= glmrobMqle.control(tcc=1.2))
coef(Rfit2)  ## The same as Rfit1


## Binary response --------------
data(vaso)

Vfit1 <- glm(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso)
coef(Vfit1)

Vfit2 <- glmrob(Y ~ log(Volume) + log(Rate), family=binomial, data=vaso,
                method="Mqle", control = glmrobMqle.control(tcc=3.5))
## Note the problems with  tcc <= 3 %% FIXME algorithm ???
coef(Vfit2) # c = 3.5 ==> not much different from classical



## Poisson response --------------
data(epilepsy)

Efit1 <- glm(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy)
summary(Efit1)

Efit2 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
                data = epilepsy, method= "Mqle",
                control = glmrobMqle.control(tcc= 1.2))
summary(Efit2)

## 'x' weighting:
(Efit3 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
 	         data = epilepsy, method= "Mqle", weights.on.x = "hat",
		 control = glmrobMqle.control(tcc= 1.2)))

try( # gives singular cov matrix: 'Trt' is binary factor -->
     # affine equivariance and subsampling are problematic
Efit4 <- glmrob(Ysum ~ Age10 + Base4*Trt, family = poisson,
                data = epilepsy, method= "Mqle", weights.on.x = "covMcd",
                control = glmrobMqle.control(tcc=1.2, maxit=100))
)

}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
