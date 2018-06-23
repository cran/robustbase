\name{steamUse}
\title{Steam Usage Data (Excerpt)}
\alias{steamUse}
\docType{data}
\encoding{utf8}
\description{
  The monthly use of steam (\code{Steam}) in a factory may be
  modeled and described as function of the
  operating days per month (\code{Operating.Days}) and
  mean outside temperature per month (\code{Temperature}).
}
\usage{data("steamUse", package="robustbase")}
\format{
  A data frame with 25 observations on the following 9 variables.
  \describe{
    \item{\code{Steam}:}{regression response \eqn{Y}, the poinds of
      steam used monthly.}
    \item{\code{fattyAcid}:}{pounds of Real Fatty Acid in storage per month.}
    \item{\code{glycerine}:}{pounds of crude glycerine made.}
    \item{\code{wind}:}{average wind velocity in miles per hour (a numeric vector).}
    \item{\code{days}:}{an integer vector with number of days of that
      month, i.e., in \eqn{28..31}.}
    \item{\code{op.days}:}{the number of operating days for the given
      month (integer).}
    \item{\code{freeze.d}:}{the number of days below
      32 degrees Fahrenheit (\eqn{= 0}\enc{°C}{'C} (C=Celsius) \eqn{=}
      freezing temperature of water).}
    \item{\code{temperature}:}{a numeric vector of average outside temperature in
      Fahrenheit (F).}
    \item{\code{startups}:}{the number of startups (of production in that month).}
  }
}
\details{
  Nor further information is given in Draper and Smith, about the place
  and exacts years of the measurements, though some educated guesses
  should be possible, see the examples.
}
\source{
  Data from Draper and Smith, 1st ed, 1966; appendix A.

  A version of this has been used in teaching at SfS ETH Zurich, since at least 1996,
  \url{https://stat.ethz.ch/Teaching/Datasets/NDK/dsteam.dat}

  The package \CRANpkg{aprean3} contains all data sets from the 3rd
  edition of Draper and Smith (1998), and this data set with variable
  names \code{x1 .. x10} (\code{x9} being \code{wind^2}, hence extraneous).
}
\references{
  Draper and Smith (1981) Applied Regression Analysis (2nd ed., p. 615 ff)
}
\examples{
\dontrun{
if(require("aprean3")) { # show how  'steamUse'  is related to 'dsa01a'
  stm <- dsa01a
  names(stm) <- c("Steam", "fattyAcid", "glycerine", "wind",
		  "days", "op.days", "freeze.d",
		  "temperature", "wind.2", "startups")
  ## prove that wind.2 is  wind^2, "traditionally" rounded to 1 digit:
  stopifnot(all.equal(floor(0.5 + 10*stm[,"wind"]^2)/10,
                      stm[,"wind.2"], tol = 1e-14))
  ## hence drop it
  steamUse <- stm[, names(stm) != "wind.2"]
}
}% dont
data(steamUse)
str(steamUse)
## Looking at this,
cbind(M=rep_len(month.abb, 25), steamUse[,5:8, drop=FALSE])
## one will conjecture that these were 25 months, Jan--Jan in a row,
## starting in a leap year (perhaps 1960 ?).

plot(steamUse)

summary(fm1 <- lmrob(Steam ~ temperature + op.days, data=steamUse))
## diagnoses 2 outliers: month of July, maybe company-wide summer vacations
%% no longer visible: summary(fmF <- lmrob(Steam ~ ., data=steamUse))
## KS2014 alone seems not robust enough:
summary(fm.14 <- lmrob(Steam ~ temperature + op.days, data=steamUse,
         setting="KS2014"))
pairs(Steam ~ temperature+op.days, steamUse)
}
\keyword{datasets}
