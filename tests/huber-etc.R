library(robustbase)

source(system.file("xtraR/styleData.R", package = "robustbase"))  # -> smallD  list of small datasets
str(smallD,, 20)

lx <- lapply(smallD,
             function(x) {
                 m <- mad(x)
                 hx <-
                     if(!is.na(m) && m > 0 && m != Inf) # in all these cases, MASS::huber() fails
                         MASS::huber(x)
                     else list(m=NA, s=NA)
                 hMx <- huberM(x)
                 list(loc =
                      c(median = median(x),
                        huber  =  hx$m,
                        huberM = hMx$m),
                      scale=
                      c(mad    = m,
                        huber  =  hx$s,
                        huberM = hMx$s))
             })


r <- list(mu = sapply(lx, function(x) x$loc),
          s  = sapply(lx, function(x) x$scale))
r

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
