\name{modsel.glmrob}
\alias{modsel.glmrob}

%% Title must fit on one line:
\title{ Robust Model Selection for "glmrob" Objects}
\description{
  Compute a robust Wald-type or deviance-type Test, i.e., a robust test
  for multiple coefficients of an \code{\link{glmrob}} object.
}
\usage{
modsel.glmrob(full.mfit, reduced.mfit, quad.form = FALSE,
             test = c("Wald", "Quasi-Deviance"))
}

\arguments{
  \item{full.mfit}{ full model: object of class \code{\link{glmrob}} }
  \item{reduced.mfit}{ reduced model: object of class \code{\link{glmrob}} }
  \item{quad.form}{ indicates if quadratic approximation is used for
    calculating the test statistic - only for  \code{"Quasi-Deviance"}-type. }
  \item{test}{ a character string specifying the test statistic to be used.
    Can be one of \code{"Wald"}-type or \code{"Quasi-Deviance"}-type, with
	  partial matching allowed. }
}

%% \details{
%% }

\value{
  Basically, an object of class \code{"htest"}, hence with the standard print
  methods for hypothesis tests.  This is basically a list with components
  \item{statistic}{the test statistic, according to ???}% FIXME
  \item{df}{numerator and denominator degrees of freedom.}
  \item{data.name}{extracted from input \code{full.mfit}.}
  \item{alternative}{\code{"two.sided"}, always.}
  \item{p.value}{the P-value, using an F-test on \code{statistic} and
    \code{df[1:2]}.}
}
\references{
  E. Cantoni and E. Ronchetti (2001)
  Robust Inference for Generalized Linear Models.
  \emph{JASA} \bold{96} (455), 1022--1030.

  E.Cantoni (2004)
  Analysis of Robust Quasi-deviances for Generalized Linear Models.
  \emph{Journal of Statistical Software} \bold{10},
  \url{http://www.jstatsoft.org}
}

\author{ Andreas Ruckstuhl }

%% \note{ ~~further notes~~ }

\seealso{ \code{\link{glmrob}}}%, \code{\link{modsel.rflm}}}
\examples{
## Binomial response
data(carrots)
Cfit2 <- glmrob(cbind(success, total-success) ~ logdose + block,
                family=binomial, data=carrots, method="Mqle",
                control=glmrobMqle.control(tcc=1.2))
summary(Cfit2)

Cfit4 <- glmrob(cbind(success, total-success) ~ logdose * block,
                family=binomial, data=carrots, method="Mqle",
                control=glmrobMqle.control(tcc=1.2))

modsel.glmrob(Cfit4, Cfit2, test="Wald")

modsel.glmrob(Cfit4, Cfit2, test="Quasi-Deviance")

modsel.glmrob(Cfit4, Cfit2, test="Quasi-Deviance", quad.form=TRUE)

## Poisson response
data(epilepsy)

Efit2 <- glmrob(Ysum ~ Age10 + Base4*Trt, family=poisson, data=epilepsy,
               method="Mqle", control=glmrobMqle.control(tcc=1.2,maxit=100))
summary(Efit2)

Efit3 <- glmrob(Ysum ~ Age10 + Base4 + Trt, family=poisson, data=epilepsy,
               method="Mqle", control=glmrobMqle.control(tcc=1.2,maxit=100))

modsel.glmrob(Efit2, Efit3, test = "Wald")

modsel.glmrob(Efit2, Efit3, test = "Quasi-Deviance")

modsel.glmrob(Efit2, Efit3, test = "Quasi-Deviance", quad.form = TRUE)
}
\keyword{robust}
\keyword{models}
\keyword{htest}

