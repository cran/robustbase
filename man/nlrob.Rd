\name{nlrob}
\alias{nlrob}
\title{Robust Fitting of Nonlinear Regression Models}
\description{
  \code{nlrob} fits a nonlinear regression model by robust
  M-estimators. The fitting is done by iterated reweighted least squares
  (IWLS).
}
\usage{
nlrob(formula, data, start, weights = NULL, na.action = na.fail,
      psi = psi.huber, test.vec = c("resid", "coef", "w"), maxit = 20,
      acc = 1e-06, algorithm = "default", control = nls.control(),
      trace = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a nonlinear \code{\link{formula}} including variables
    and parameters of the model, such as \code{y ~ f(x, alpha)} (cf. \code{\link{nls}}).
    (For some checks: if \eqn{f(.)} is linear, then we need
    parentheses, e.g., \code{y ~ (a + b * x)}.)
%% FIXME in code:
    Do not use as \code{w} as variable or parameter name!
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
    \eqn{e_i = y_{i} - f(xreg_{i}, alpha's)}{e[i] = y[i] - f(xreg[i], alpha's)},
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
    constant(s)', deriv)} that for deriv=0 returns psi(x)/x and for
    deriv=1 returns psi'(x). Tuning constants will be passed in via one
    or several arguments, depending on the psi function. (see also
    \code{\link[MASS]{rlm}}).}
  \item{test.vec}{character string specifying the convergence
    criterion. The relative change is tested for residuals with a value
    of \code{"resid"} (the default), for coefficients with
    \code{"coef"}, and for weights with \code{"w"}.}
  \item{maxit}{maximum number of iterations in the robust loop.}
  \item{acc}{convergence tolerance for the robust fit.}
  \item{algorithm}{character string specifying the algorithm to use for
    \code{nls}. The default algorithm is a Gauss-Newton algorithm. The other
    alternative is \code{"plinear"}, the Golub-Pereyra algorithm for
    partially linear least-squares models.}
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
  \item{\dots}{potentially arguments to be passed to the psi function
    (see above).}
}
% \details{

% }
\value{
  An object of S3 class \code{"nlrob"}, also inheriting from class
  "nls", (see \code{\link{nls}}).

  There methods (at least) for the generic accessor functions
  \code{\link{summary}}, \code{\link{coefficients}},
  \code{fitted.values} and \code{residuals}.

  It is a list with components
  \item{FIXME }{???}
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
fm3DNase1 <- nls( density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                  data = DNase1,
                  start = list( Asym = 3, xmid = 0, scal = 1 ),
                  trace = TRUE )
summary( fm3DNase1 )

## robust
Rm3DNase1 <- nlrob(density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                   data = DNase1, trace = TRUE,
                   start = list( Asym = 3, xmid = 0, scal = 1 ))
summary( Rm3DNase1 )

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

### -- new examples
DNase1[10,"density"] <- 2*DNase1[10,"density"]

fm3DNase1 <- nls(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                       data = DNase1, trace = TRUE,
                       start = list( Asym = 3, xmid = 0, scal = 1 ))

## robust
Rm3DNase1 <- nlrob(density ~  Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
                   data = DNase1, trace = TRUE,
                   start = list( Asym = 3, xmid = 0, scal = 1 ))
Rm3DNase1 ## summary() is not yet there {FIXME}

## utility function sfsmisc::lseq() :
lseq <- function (from, to, length)
  2^seq(log2(from), log2(to), length.out = length)
## predict() {and plot}:
h.x <- lseq(min(DNase1$conc), max(DNase1$conc), length = 100)
nDat <- data.frame(conc = h.x)

h.p  <- predict(fm3DNase1, newdata = nDat)# classical
h.rp <- predict(Rm3DNase1, newdata= nDat)# robust

plot(density ~ conc, data=DNase1, log="x",
     main = deparse(Rm3DNase1$call$formula))
lines(h.x, h.p,  col="blue")
lines(h.x, h.rp, col="magenta")
legend("topleft", c("classical nls()", "robust nlrob()"),
       lwd = 1, col= c("blue", "magenta"), inset = 0.05)
}
\keyword{robust}
\keyword{regression}
\keyword{nonlinear}
