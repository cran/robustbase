library(robustbase)

## minimal testing only
data(ruspini, package = "cluster")

### NOTA BENE:  scale.tau2() is *not* consistent {by constant factor}
rub1 <- covOGK(ruspini, 1, scaleTau2, covGK, hard.rejection)
rub2 <- covOGK(ruspini, 2, scaleTau2, covGK, hard.rejection)

AE <- function(x,y) all.equal(x,y, tol = 2e-15)
## The following test is already fulfilled by Kjell Konis'  original code:
stopifnot(AE(c(rub1$wcov)[c(1,3:4)],
             c(917.99893333333, 94.9232, 2340.319288888888)),
          all.equal(rub1$wcov, rub2$wcov, tol=0)
          ,
          AE(c(rub1$cov)[c(1,3:4)],
             c(923.5774514441657, 91.5385216376565, 2342.4556232436971))
          ,
          AE(c(rub2$cov)[c(1,3:4)],
             c(927.2465953711782, 91.8009184487779, 2346.5790105548940))
          )

## More tests:
## set.seed(101)
