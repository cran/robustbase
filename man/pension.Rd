\name{pension}
\alias{pension}
\docType{pension}
\title{Pension Funds Data}
\description{
 The total 1981 premium income of pension funds of Dutch firms,
 for 18 Professional Branches, from de Wit (1982).
}
\usage{data(pension)}
\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{Income}}{Premium Income (in millions of guilders)}
    \item{\code{Reserves}}{Premium Reserves (in millions of guilders)}
  }
}
\source{
  P. J. Rousseeuw and A. M. Leroy (1987)
  \emph{Robust Regression and Outlier Detection};
  Wiley, p.76, table 13.
}
\examples{
data(pension)
summary(lm.pension <- lm(Reserves ~., data=pension))
%%summary(lts.pension <- ltsReg(Reserves ~., data=pension))
}
\keyword{datasets}
