\name{internals}
\alias{internals}

\title{Internal Functions}
%% MM: This page is not *needed*, since all the functions are namespace-hidden
%% --  OTOH, it shows what we might want to export {and then document}
\description{
  These functions are for internal use only or not yet documented.
}
\usage{
%% for glmrob
\method{print}{glmrob}(x, digits = max(3, getOption("digits") - 3), ...)

glmrobMqle(X, y, weights = NULL, start = NULL, offset = NULL,
           family, weights.on.x = "none",
           control = glmrobMqle.control(), intercept = TRUE, trace = FALSE)
% glmrobMqleEpsiB     (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiPois  (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiSB    (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsiSPois (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsi2B    (mu, Vmu, ni, H, K, tcc)
% glmrobMqleEpsi2Pois (mu, Vmu, ni, H, K, tcc)

glmrobMqleDiffQuasiDevB   (mu, mu0, y, ni, w.x, tcc)
glmrobMqleDiffQuasiDevPois(mu, mu0, y, ni, w.x, tcc)

wts_HiiDist(X)
wts_RobDist(X, intercept, covFun)

% rflmFixedScale(x, y, weights, w = rep(1, nrow(x)), init, scale,
%                psi = psi.huber, maxit = 20, acc = 1e-04,
%                test.vec = "resid", wt.method = "inv.var")

%% for nlrob
\method{print}{nlrob}(x, \dots)

%% for modsel.* aka anova.glmrob():
psi.bisquare(u, c = 4.685, deriv = 0, rho = FALSE)
psi.hampel(u, a = 2, b = 4, c = 8, deriv = 0, rho = FALSE)
psi.huber(u, k = 1.345, deriv = 0, rho = FALSE)
}

\author{Andreas Ruckstuhl}

\keyword{internal}

