\name{glmrobMqle.control}
\alias{glmrobMqle.control}
\title{Auxiliary for Controlling Robust GLM Fitting by the Method "Mqle"}
\description{
  Auxiliary function as user interface for \code{\link{glmrob}} fitting
  when method "Mqle" is applied.  Typically only used when calling
  \code{\link{glmrob}}.
}

\usage{
glmrobMqle.control(acc = 1e-04, test.acc = "coef", maxit = 50, tcc = 1.345)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{acc}{positive convergence tolerance;
    the iterations converge when ???}
  \item{test.acc}{Only "coef" is currently implemented}
  \item{maxit}{integer giving the maximum number of iterations. }
  \item{tcc}{tuning constant c for Huber's psi-function}
}

%% \details{
%% }

\value{
  A list with the arguments as components.
}

%% \references{ ~put references to the literature/web site here ~ }

\author{Andreas Ruckstuhl}

%% \\note{ ~~further notes~~ }
%%
%%  ~Make other sections like Warning with \section{Warning }{....} ~
%%

\seealso{\code{\link{glmrob}}}

%% \examples{
%% }

\keyword{robust}
\keyword{regression}
\keyword{nonlinear}

