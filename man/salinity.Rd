\name{salinity}
\alias{salinity}
\docType{data}
\title{Salinity Data}
\description{
  This is a data set consisting of measurements of water salinity (i.e.,
  its salt concentration) and river discharge taken in North Carolina's
  Pamlico Sound; This dataset was listed by Ruppert and Carroll
  (1980).  In Carrol and Ruppert (1985) the physical background of the
  data is described.  They indicated that observations 5 and 16
  correspond to periods of very heavy discharge and showed that the
  discrepant observation 5 was masked by observations 3 and 16, i.e.,
  only after deletion of these observations it was possible to identify
  the influential observation 5.

  This data set is a prime example of the masking effect.
}
\usage{data(salinity)}
\format{
  A data frame with 28 observations on the following 4 variables.
  \describe{
    \item{\code{X1}}{Lagged Salinity}
    \item{\code{X2}}{Trend}
    \item{\code{X3}}{Discharge}
    \item{\code{Y}}{Salinity}
  }
}
\source{
 P. J. Rousseeuw and A. M. Leroy (1987)
 \emph{Robust Regression and Outlier Detection};
 Wiley, p.82, table 5.
}
\examples{
data(salinity)
summary(lm.sali  <-        lm(Y ~ . , data = salinity))
summary(rlm.sali <- MASS::rlm(Y ~ . , data = salinity))
summary(lts.sali <-    ltsReg(Y ~ . , data = salinity))

salinity.x <- data.matrix(salinity[, 1:3])
c_sal <- covMcd(salinity.x)
plot(c_sal, "tolEllipsePlot")
}
\keyword{datasets}
