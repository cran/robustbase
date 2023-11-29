\name{nlrob.control}
\alias{nlrob.control}
\title{Control Nonlinear Robust Regression Algorithms}
\description{
  Allow the user to specify details for the different nonlinear robust
  regression algorithms in \code{\link{nlrob}}.
}
\usage{
nlrob.control(method,
              psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
              init = c("S", "lts"),
              optimizer = "JDEoptim", optArgs  = list(),
              ...)
}
\arguments{
  \item{method}{\code{\link{character}} string specifying the method}
  \item{psi}{string specifying the psi-function which defines the estimator.}
  \item{init}{for some methods, currently, \code{"MM"} only, a string
    specifying the initial estimator.
  }
  \item{optimizer}{currently only \code{"JDEoptim"} from package \CRANpkg{DEoptimR}.}
  \item{optArgs}{
    a \code{\link{list}} of optional arguments to the optimizer.
    Currently, that is \code{\link[DEoptimR]{JDEoptim}} from package
    \CRANpkg{DEoptimR}.
  }
  \item{\dots}{optional arguments depending on \code{method}, such as
    \code{fnscale}, \code{tuning.chi} or both \code{tuning.chi.tau}  and
    \code{tuning.chi.scale}; for \code{method = "MM"} also \code{optim.control}
    to be passed to the \code{\link{optim}(.., hessian=TRUE)} call.
    %%
    Internally, \code{nlrob.control()} will choose (or check) defaults for
    the psi/rho/chi related tuning parameters, also depending on the \code{method}
    chosen; see e.g., the \sQuote{Examples}.
  }
}
%% Code in >>> ../R/nlregrob.R <<<
%%
\value{
  a \code{\link{list}} with several named components.
  The contents depend quite a bit on the \code{method}.
}
%\author{Martin Maechler}
\seealso{
  \code{\link{nlrob}}; for some details, \code{\link{nlrob.algorithms}}.
}
\examples{
## Show how the different 'method's  have different smart defaults :
str(nlrob.control("MM"))
str(nlrob.control("MM", psi = "hampel"))# -> other tuning.psi.M and tuning.chi.scale defaults
str(nlrob.control("MM", psi = "lqq", tol = 1e-10))# other tuning.psi.M & tuning.chi.scale defaults
str(nlrob.control("tau"))
str(nlrob.control("tau",psi= "lqq"))
str(nlrob.control("CM")) # tuning.chi undefined, unneeded
str(nlrob.control("CM", psi= "optimal"))
str(nlrob.control("mtl"))
}
\keyword{utilities}
