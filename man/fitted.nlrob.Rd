\name{fitted.nlrob}
\alias{fitted.nlrob}

\title{Extract Fitted Values from "nlrob" Model}
\description{
  \code{\link{fitted}} is a generic function which extracts fitted values from
  objects returned by modeling functions. \code{fitted.values} is an
  alias for it.
}
\usage{
\method{fitted}{nlrob}(object, ...)
}

\arguments{
  \item{object}{object of class \code{\link{nlrob}} }
  \item{\dots}{ additional arguments }
}

\value{
  A numeric vector of fitted values extracted from \code{object}.
}

\author{ Andreas Ruckstuhl }

\seealso{\code{\link[stats]{fitted}},
  \code{\link{nlrob}}
}

\keyword{ regression }
