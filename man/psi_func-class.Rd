\name{psi_func-class}
\docType{class}
\alias{psi_func-class}
\title{Class of "Psi Functions" for M-Estimation}
%% 'psi_func' cannot be used in title!
\description{
  The class \code{"psi_func"} is used to store \eqn{\psi \ (psi)}{psi}
  functions for M-estimation.  In particular, an object of the class
  contains \eqn{\rho(x) \ (\code{rho})}{rho(x)}, its derivative
  \eqn{\psi(x) \ (psi)}{psi(x)}, the weight function \eqn{\psi(x)/x}, and
  first derivative of \eqn{\psi}, \code{Dpsi = } \eqn{\psi'(x)}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("psi_func", ...)},
  but preferably by \code{\link{psiFunc}(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{rho}:}{the \eqn{\rho()}{rho()} function, an object of
      class \code{"functionX"}. This is used to formulate the
      objective function; \eqn{\rho()} can be regarded as generalized
      negative log-likelihood.}
    \item{\code{psi}:}{\eqn{\psi()}{psi()} is the derivative of \eqn{\rho},
      \eqn{\psi(x) = \frac{d}{dx} \rho(x)}{psi(x) = d/dx rho(x)};
      also of class \code{"functionX"}.}
    \item{\code{wgt}:}{The weight function \eqn{\psi(x)/x},
      of class \code{"functionX"}.}
    \item{\code{Dpsi}:}{the derivative of \eqn{\psi},
      \eqn{Dpsi(x) = psi'(x)}; of class \code{"functionX"}.}
    \item{\code{tDefs}:}{\emph{named} numeric vector of \bold{t}uning
      parameter \bold{Def}ault values.}
    %%% FIXME : Replace these by *methods*
    \item{\code{Erho}:}{Object of class \code{"functionXal"} ~~ }
    \item{\code{Epsi2}:}{Object of class \code{"functionXal"} ~~ }
    \item{\code{EDpsi}:}{Object of class \code{"functionXal"} ~~ }

    \item{\code{xtras}:}{Potentially further information.}
  }
}
\section{Methods}{
  No methods defined with class \code{"psi_func"} in the signature.
}
\author{Martin Maechler}
% \note{ ~~further notes~~ }

%  ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{psiFunc}}.
}
\examples{
str(huberPsi, give.attr = FALSE)
}
\keyword{classes}
\keyword{robust}
