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

###--- now the "close to singular" mahalanobis case:
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

## Now, some small sample cases:

## maximal values:
n. <- 10
p. <-  8
set.seed(44)
(X. <- cbind(1:n., round(10*rt(n.,3)), round(10*rt(n.,2)),
             matrix(round(10*rnorm(n. * (p.-3)), 1),  nrow = n., ncol = p.-3)))

## 2 x 1 ---> Error
r <- try(covMcd(X.[1:2, 2, drop=FALSE]), silent=TRUE)
stopifnot(inherits(r, "try-error"),
          grep("too small sample size", r) == 1)

## 3 x 2 --- ditto
r <- try(covMcd(X.[1:3, 2:3]), silent=TRUE)
stopifnot(inherits(r, "try-error"),
          grep("too small sample size", r) == 1)

## 5 x 3  [ n < 2 p  ! ]  --- also works for MASS
X <- X.[1:5, 1:3]
set.seed(101)
## the finite-sample correction is definitely doubtful:
(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 2 3 4 5
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tol = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.673549282206, 3*3)))

## p = 4 -- 6 x 4 & 7 x 4  [ n < 2 p  ! ]
p <- 4
n <- 7
X <- X.[1:n, 1+(1:p)]
stopifnot(dim(X) == c(n,p))
(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 1 2 4 5 6 7
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tol = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7782486992881, p*p)))
n <- 6
X <- X[1:n,]
(cc <- covMcd(X, use.correction = FALSE))
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tol = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7528695976179, p*p)))

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

## nsamp = "exact" -- here for p=7
coleman.x <- data.matrix(coleman[, 1:6])
system.time(CcX <- covMcd(coleman.x, nsamp="exact")) # ~ 3 sec.
stopifnot(all.equal(CcX$best,
                    c(2, 5:9, 11,13, 14:16, 19:20), tol=0))
