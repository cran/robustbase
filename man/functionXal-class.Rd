\name{functionXal-class}
\docType{class}
\alias{functionXal-class}
\title{Class "functionXal" of Functionals (of Psi-like functions)}
\description{
  The class \code{"functionXal"} is a class of functionals (typically
  integrals) typically of \code{\link{functionX-class}} functions.

  Since the \code{functionX} functions typically also depend on tuning
  parameters, objects of this class (\code{"functionXal"}) are functions
  of these tuning parameters.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Directly extends class \code{"function"}.}
  }
}
\section{Extends}{
  Class \code{"function"}, from data part.
  Class \code{"OptionalFunction"}, by class \code{"function"}.
  Class \code{"PossibleMethod"}, by class \code{"function"}.
}
% \section{Methods}{
%   No methods defined with class "functionXal" in the signature.
% }
\seealso{
  \code{\link{functionX-class}};
  \code{\link{psiFunc}}, and
  \code{\link{psi_func-class}} which has several \code{functionXal} slots.
}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% }
\keyword{classes}
