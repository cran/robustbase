\name{predict.lmrob}
\alias{predict.lmrob}
\title{Predict method for Robust Linear Model ("lmrob") Fits}
\description{
  Predicted values based on robust linear model object.
}
\usage{
\method{predict}{lmrob}(object, newdata, se.fit = FALSE,
       scale = NULL, df = NULL,
       interval = c("none", "confidence", "prediction"), level = 0.95,
       type = c("response", "terms"), terms = NULL,
       na.action = na.pass, pred.var = res.var/weights, weights = 1, ...)
}
\arguments{
%% the following is +- copy-pasted from  predict.lm.Rd:
  \item{object}{object of class inheriting from \code{"lmrob"}}
  \item{newdata}{an optional data frame in which to look for variables with
    which to predict.  If omitted, the fitted values are used.}
  \item{se.fit}{a switch indicating if standard errors are required.}
  \item{scale}{scale parameter for std.err. calculation}
  \item{df}{degrees of freedom for scale}
  \item{interval}{type of interval calculation.}
  \item{level}{tolerance/confidence level}
  \item{type}{Type of prediction (response or model term).}
  \item{terms}{if \code{type="terms"}, which terms (default is all terms)}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to predict \code{NA}.}
  \item{pred.var}{the variance(s) for future observations to be assumed
    for prediction intervals.  See \sQuote{Details}.}
  \item{weights}{variance weights for prediction. This can be a numeric
    vector or a one-sided model formula. In the latter case, it is
    interpreted as an expression evaluated in \code{newdata}}
  \item{\dots}{further arguments passed to or from other methods.}
}
% \details{
% }
\value{
%% the following is +- copy-pasted from  predict.lm.Rd:
  \code{predict.lmrob} produces a vector of predictions or a matrix of
  predictions and bounds with column names \code{fit}, \code{lwr}, and
  \code{upr} if \code{interval} is set.  If \code{se.fit} is
  \code{TRUE}, a list with the following components is returned:
  \item{fit}{vector or matrix as above}
  \item{se.fit}{standard error of predicted means}
  \item{residual.scale}{residual standard deviations}
  \item{df}{degrees of freedom for residual}
}
% \references{

% }
\author{Andreas Ruckstuhl}
\seealso{
  \code{\link{lmrob}} and the (non-robust) traditional
  \code{\link{predict.lm}} method.
}
% \examples{
% }
\keyword{robust}
\keyword{regression}
