\name{radarImage}
\alias{radarImage}
\docType{data}
\title{Satellite Radar Image Data from near Munich}
\description{
  The data were supplied by A. Frery.  They are a part of a synthetic aperture
  satellite radar image corresponding to a suburb of Munich.  Provided
  are coordinates and values corresponding to three frequency bands for
  each of 1573 pixels.
}
\usage{data(radarImage)}
\format{
  A data frame with 1573 observations on the following 5 variables.
  \describe{
    \item{\code{X.coord}}{a numeric vector}
    \item{\code{Y.coord}}{a numeric vector}
    \item{\code{Band.1}}{a numeric vector}
    \item{\code{Band.2}}{a numeric vector}
    \item{\code{Band.3}}{a numeric vector}
  }
}
% \details{
%   ~~ If necessary, more details than the __description__ above ~~
% }
\source{
  The website accompanying the MMY-book:
  \url{http://www.wiley.com/legacy/wileychi/robust_statistics}
}
% \references{
%   ~~ possibly secondary sources and usages ~~
% }
\examples{
data(radarImage)
plot(Y.coord ~ X.coord, data = radarImage)
## see outliers
pairs(radarImage[, 3:5], main = "radarImage (n = 1573)")
}
\keyword{datasets}
