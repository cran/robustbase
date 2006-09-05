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

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
