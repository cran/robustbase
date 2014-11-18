\name{internals}
\alias{internals}
\title{Internal Functions of Package 'robustbase'}
%
\alias{robMD}
\alias{mahalanobisD}
%
\alias{r6pack}
\alias{doScale}
\alias{.signflip}
%
\description{
  These functions are for internal use \emph{or} not yet documented.
}
\usage{
%% for glmrob
\method{print}{glmrob}(x, digits = max(3, getOption("digits") - 3), ...)

%% MM: These are not *needed*, since they are namespace-hidden
%% --  OTOH, it shows what we might want to export {and then document}
glmrobMqle(X, y, weights = NULL, start = NULL, offset = NULL,
           family, weights.on.x = "none",
           control = glmrobMqle.control(), intercept = TRUE, trace = FALSE)
% glmrobMqleEpsiB     (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiPois  (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiSB    (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiSPois (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsi2B    (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsi2Pois (mu, Vmu, ni, H, K, tcc)

glmrobMqleDiffQuasiDevB   (mu, mu0, y, ni, w.x, phi, tcc)
glmrobMqleDiffQuasiDevPois(mu, mu0, y, ni, w.x, phi, tcc)

% rflmFixedScale(x, y, weights, w = rep(1, nrow(x)), init, scale,
%                psi = psi.huber, maxit = 20, acc = 1e-04,
%                test.vec = "resid", wt.method = "inv.var")


%% For modsel.* aka anova.glmrob(); later for nlrob() only;
%% eliminated, 2013-07-26 {for robustbase 0.9-11}, MM :
%% psi.bisquare(u, c = 4.685, deriv = 0, rho = FALSE)
%% psi.hampel(u, a = 2, b = 4, c = 8, deriv = 0, rho = FALSE)
%% psi.huber(u, k = 1.345, deriv = 0, rho = FALSE)

%%--------From here on, are exported (and MM would want minimal doc, evt.) ------
robMD(x, intercept, wqr, ...)

mahalanobisD(x, center, sd)

## Utilities currently for the deterministic MCD only:
## subject to change / be renamed ?
r6pack(x, h, hsets.init, full.h, scaled=TRUE,
       scalefn = rrcov.control()$scalefn)
doScale(x, center, scale)
.signflip(loadings)
}

\keyword{internal}

