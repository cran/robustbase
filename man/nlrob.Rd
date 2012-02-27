\name{nlrob}
\title{Robust Fitting of Nonlinear Regression Models}
\alias{nlrob}
\alias{fitted.nlrob}
\alias{residuals.nlrob}
\alias{predict.nlrob}
\description{
  \code{nlrob} fits a nonlinear regression model by robust
  M-estimators, using iterated reweighted least squares (IWLS).
}
\usage{
nlrob(formula, data, start, weights = NULL, na.action = na.fail,
      psi = psi.huber, test.vec = c("resid", "coef", "w"), maxit = 20,
      acc = 1e-06, algorithm = "default", control = nls.control(),
      trace = FALSE, ...)

\method{fitted}{nlrob}(object, ...)
\method{residuals}{nlrob}(object, type = , ...)% FIXME: more 'type's + DOCU
\method{predict}{nlrob}(object, newdata, ...)
}
\arguments{
  \item{formula}{a nonlinear \code{\link{formula}} including variables
    and parameters of the model, such as \code{y ~ f(x, theta)} (cf. \code{\link{nls}}).
    (For some checks: if \eqn{f(.)} is linear, then we need
    parentheses, e.g., \code{y ~ (a + b * x)}.)
%% FIXME in code -- long overdue, as nls() is more flexible *SINCE* R 2.2.1
    Do not use \code{w} as variable or parameter name!
  }
  \item{data}{an optional data frame containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from
    which \code{nlrob} is called.}
  \item{start}{a named numeric vector of starting parameters estimates.}
  \item{weights}{an optional vector of weights to be used in the fitting
    process (for intrinsic weights, not the weights \code{w} used in the
    iterative (robust) fit). I.e.,
    \code{sum(w * e^2)} is minimized with \code{e} = residuals,
    \eqn{e_i = y_i - f(xreg_i, \theta)}{e[i] = y[i] - f(xreg[i], theta)},
    where \eqn{f(x,\theta)}{f(x, theta)} is the nonlinear function,
    and \code{w} are the robust weights from \code{resid * weights}.}
  \item{na.action}{a function which indicates what should happen when the data
    contain \code{NA}s.  The default action is for the procedure to
    fail.  If NAs are present, use \code{na.exclude} to have residuals with
    \code{length == nrow(data) == length(w)}, where \code{w} are the
    weights used in the iterative robust loop.  This is better if the explanatory
    variables in
    \code{formula} are time series (and so the NA location is important).
    For this reason, \code{na.omit}, which leads to omission of cases
    with missing values on any required variable, is not suitable
    here since the residuals length is different from
    \code{nrow(data) == length(w)}.
  }
  \item{psi}{a function (possibly by name) of the form \code{g(x, 'tuning
      constant(s)', deriv)} that for \code{deriv=0} returns
    \eqn{\psi(x)/x}{psi(x)/x} and for \code{deriv=1} returns
    \eqn{\psi'(x)}{psi'(x)}.  Tuning constants will be passed in via one
    or several arguments, depending on the psi function. (See also
    \code{\link[MASS]{rlm}}).}
  \item{test.vec}{character string specifying the convergence
    criterion. The relative change is tested for residuals with a value
    of \code{"resid"} (the default), for coefficients with
    \code{"coef"}, and for weights with \code{"w"}.}
  \item{maxit}{maximum number of iterations in the robust loop.}
  \item{acc}{convergence tolerance for the robust fit.}
  \item{algorithm}{character string specifying the algorithm to use for
    \code{\link{nls}}, see there. The default algorithm is a
    Gauss-Newton algorithm.}
  \item{control}{an optional list of control settings for \code{\link{nls}}.
    See \code{\link{nls.control}} for the names of the settable control
    values and their effect.}
  \item{trace}{logical value indicating if a \dQuote{trace} of
    the \code{nls} iteration progress should be printed.  Default is
    \code{FALSE}. \cr
    If \code{TRUE}, in each robust iteration, the residual
    sum-of-squares and the parameter values are printed at the
    conclusion of each \code{nls} iteration.
    When the \code{"plinear"} algorithm is used, the conditional
    estimates of the linear parameters are printed after the nonlinear
    parameters.}
  \item{\dots}{optional arguments (tuning constants) to be passed to the
    psi function (see \code{psi} above).}
  \item{object}{an \R object of class \code{"nlrob"}, typically
    resulting from \code{nlrob(..)}.}
  \item{type}{a string specifying the \emph{type} of residuals desired.
    Currently, \code{"response"} and \code{"working"} are supported.
    %% FIXME: 1. document these (here)  2. write and support more types
  }
  \item{newdata}{a data frame (or list) with the same names as the
    original \code{data}, see e.g., \code{\link{predict.nls}}.}
}
% \details{

% }
\value{
  \code{nlrob()} returns an object of S3 class \code{"nlrob"}, also
  inheriting from class \code{"nls"}, (see \code{\link{nls}}).

  It is a list with components
  \item{FIXME }{???}

  There are methods (at least) for the generic accessor functions
  \code{\link{summary}}, \code{\link{coefficients}},
  \code{fitted.values} and \code{residuals}.

  \code{residuals(.)}, by default \code{type = "response"}, returns
  the residuals \eqn{e_i}, defined above as
  \eqn{e_i = Y_i - f_(x_i, \theta)}{e[i] = Y[i] - f(x[i],theta)}.
  These differ from the standardized or weighted residuals which, e.g.,
  are assumed to be normally distributed, and a version of which is
  returned in \code{working.residuals} component.
  %% and another is working.residuals/Scale
}

\author{
  Andreas Ruckstuhl (inspired by \code{\link[MASS]{rlm}}() and
  \code{\link{nls}}()), in July 1994 for S-plus.\cr
  Christian Sangiorgio did the update to \R and corrected some errors,
  from June 2002 to January 2005, and Andreas contributed slight changes
  and the first methods in August 2005.\cr
  Help page, testing, more cleanup, methods: Martin Maechler.
}
\note{
  This function used to be named \code{\link[sfsmisc]{rnls}}
  and has been in package \pkg{sfsmisc}, but will be deprecated and dropped
  there, eventually.
}

\seealso{ \code{\link{nls}}, \code{\link[MASS]{rlm}}.
}
\examples{
DNase1 <- DNase[ DNase$Run == 1, ]

## note that selfstarting models don't work yet % <<< FIXME !!!

##--- without conditional linearity ---

## classical
fmNase1 <- nls( density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DNase1,
                start = list( Asym = 3, xmid = 0, scal = 1 ),
                trace = TRUE )
summary( fmNase1 )

## robust
RmN1  <- nlrob( density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DNase1, trace = TRUE,
                start = list( Asym = 3, xmid = 0, scal = 1 ))
summary( RmN1 )

##--- using conditional linearity ---

## classical
fm2DNase1 <- nls( density ~ 1/(1 + exp(( xmid - log(conc) )/scal ) ),
                  data = DNase1,
                  start = c( xmid = 0, scal = 1 ),
                  alg = "plinear", trace = TRUE )
summary( fm2DNase1 )

## robust
if(FALSE) { # currently fails %% FIXME error in nls's nlsModel.plinear()
frm2DNase1 <- nlrob(density ~ 1/(1 + exp(( xmid - log(conc) )/scal ) ),
                  data = DNase1, start = c( xmid = 0, scal = 1 ),
                  alg = "plinear", trace = TRUE )
summary( frm2DNase1 )
} # not yet

### -- new examples -- "moderate outlier":
DN2 <- DNase1
DN2[10,"density"] <- 2*DN2[10,"density"]

fm3DN2 <- nls(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
              data = DN2, trace = TRUE,
              start = list( Asym = 3, xmid = 0, scal = 1 ))

## robust
Rm3DN2 <- nlrob(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                data = DN2, trace = TRUE,
                start = list( Asym = 3, xmid = 0, scal = 1 ))
Rm3DN2
summary(Rm3DN2)

## utility function sfsmisc::lseq() :
lseq <- function (from, to, length)
  2^seq(log2(from), log2(to), length.out = length)
## predict() {and plot}:
h.x <- lseq(min(DN2$conc), max(DN2$conc), length = 100)
nDat <- data.frame(conc = h.x)

h.p  <- predict(fm3DN2, newdata = nDat)# classical
h.rp <- predict(Rm3DN2, newdata = nDat)# robust

plot(density ~ conc, data=DN2, log="x",
     main = format(formula(Rm3DN2)))
lines(h.x, h.p,  col="blue")
lines(h.x, h.rp, col="magenta")
legend("topleft", c("classical nls()", "robust nlrob()"),
       lwd = 1, col= c("blue", "magenta"), inset = 0.05)
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
