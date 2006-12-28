library(robustbase)

source(system.file("test_MCD.R", package = "robustbase"))
##          ../inst/test_MCD.R

## -- now do it:
options(digits = 5)
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed
doMCDdata()
##                vvvv no timing for 'R CMD Rdiff' outputs
doMCDdata(nrep = 12, time=FALSE)
doMCDdata(nrep = 12, time=FALSE, method = "MASS")

###--- now the "close to singular mahalanobis case:
(c3 <- covMcd(mort3))
## rescale variables:
scaleV <- c(0.1, 0.1, 1, 1, .001, 0.1, 0.1, 100)
mm <- data.matrix(mort3) * rep(scaleV, each = nrow(mort3))
C3 <- covMcd(mm)
stopifnot(C3$mcd.wt == c3$mcd.wt)
try(## error: with "old default tolerance:
  covMcd(mm, control= rrcov.control(tol = 1e-10))
)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
