### test subsample
### LU decomposition and singular subsamples handling
require(robustbase)
source(system.file("xtraR/subsample-fns.R", package = "robustbase", mustWork=TRUE))

A <- matrix(c(0.001, 1, 1, 2), 2)
set.seed(11)
str(sa <- subsample(A))

A <- matrix(c(3, 2, 6, 17, 4, 18, 10, -2, 12), 3)
subsample(A)

## test some random matrix
set.seed(1002)
A <- matrix(rnorm(100), 10)
subsample(A)

## test singular matrix handling
A <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1), 4, byrow=TRUE)
subsample(A)


## test subsample with mts > 0
data <- data.frame(y = rnorm(9), expand.grid(A = letters[1:3], B = letters[1:3]))
x <- model.matrix(y ~ ., data)
y <- data$y
## this should produce a warning and return status == 2
z <- Rsubsample(x, y, mts=2)
stopifnot(z$status == 2)


## test equilibration
## columns only
X <- matrix(c(1e-7, 2, 1e-10, 0.2), 2)
y <- 1:2
subsample(t(X), y)

## rows only
X <- matrix(c(1e-7, 2, 1e10, 0.2), 2)
y <- 1:2
subsample(X, y)

## both
X <- matrix(c(1e-7, 1e10, 2, 2e12), 2)
y <- 1:2
subsample(X, y)


## test real data example
data(possumDiv)## 151 * 9; the last two variables are factors
with(possumDiv, table(eucalyptus, aspect))

mf <- model.frame(Diversity ~ .^2, possumDiv)
X <- model.matrix(mf, possumDiv)
y <- model.response(mf)
stopifnot(qr(X)$rank == ncol(X))

## this used to fail: different pivots in step 37
str(s1 <- subsample(X, y))
s2 <- subsample(X / max(abs(X)), y / max(abs(X)))
s3 <- subsample(X * 2^-50, y * 2^-50)
## all components *BUT*  x, y, lu, Dr, Dc, rowequ, colequ :
nm <- names(s1); nm <- nm[is.na(match(nm, c("x","y","lu", "Dr", "Dc", "rowequ", "colequ")))]
stopifnot(all.equal(s1[nm], s2[nm], tol=1e-10),
	  all.equal(s1[nm], s3[nm], tol=1e-10))

## test subsampling
testSubSampling <- function(X, y) {
    lX <- X[sample(nrow(X)), ]
    ## C version
    zc <- Rsubsample(lX, y)
    ## R version
    zR <- LU.gaxpy(t(lX))
    if (as.logical(zc$status)) {
        ## singularity in C detected
        if (!zR$singular)
            stop("singularity in C but not in R")
    } else {
        ## no singularity detected
        if (zR$singular)
            stop("singularity in R but not in C")
    }
    zR$singular
}

set.seed(10)
nsing <- sum(replicate(200, testSubSampling(X, y)))
stopifnot(nsing == 0)

## test example with many categorical predictors
set.seed(10)
r1 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none")
## lmrob.S used to fail for this seed:
set.seed(108)
lmrob(Diversity ~ .^2 , data = possumDiv, cov="none") #, trace=4)

## investigate problematic subsample:
idc <- 1 + c(140, 60, 12, 13, 89, 90, 118, 80, 17, 134, 59, 94, 36,
         43, 46, 93, 107, 62, 57, 116, 11, 45, 35, 38, 120, 34, 29,
         33, 147, 105, 115, 92, 61, 91, 104, 141, 138, 129, 130, 84,
         119, 132, 6, 135, 112, 16, 67, 41, 102, 76, 111, 82, 148, 24,
         131, 10, 96, 0, 87, 21, 127, 56, 124)

rc <- lm(Diversity ~ .^2 , data = possumDiv, subset = idc)

X <- model.matrix(rc)
y <- possumDiv$Diversity[idc]
subsample(X, y)

lu <- LU.gaxpy(t(X))
stopifnot(lu$sing)
zc <- Rsubsample(X, y)
stopifnot(zc$status > 0)
## column 52 is linearly dependent and should have been discarded
## qr(t(X))$pivot

image(as(round(zc$lu - (lu$L + lu$U - diag(nrow(lu$U))), 10), "Matrix"))
image(as(sign(zc$lu) - sign(lu$L + lu$U - diag(nrow(lu$U))), "Matrix"))


## test equilibration
## colequ only
X <- matrix(c(1e-7, 2, 1e-10, 0.2), 2)
y <- 1:2
subsample(t(X), y)

## rowequ only
X <- matrix(c(1e-7, 2, 1e10, 0.2), 2)
y <- 1:2
subsample(X, y)

## both
X <- matrix(c(1e-7, 1e10, 2, 2e12), 2)
y <- 1:2
subsample(X, y)


### real data, see MM's ~/R/MM/Pkg-ex/robustbase/hedlmrob.R
##  close to singular cov():
attach(system.file("external", "d1k27.rda", package="robustbase", mustWork=TRUE))

fm1 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27)
##     ^^^^^ gave error, earlier, now with a warning -- use ".vcov.w"
## --> cov = ".vcov.w"
fm2 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27,
             cov = ".vcov.w", trace = TRUE)

## Q: does it change to use numeric instead of binary factors ?
## A: not really ..
d1k.n <- d1k27
d1k.n[-(1:5)] <- lapply(d1k27[,-(1:5)], as.numeric)

fm1.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n)
fm2.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n,
             cov = ".vcov.w", trace = 2)

summary(weights(fm1))
hist(weights(fm1), main="robustness weights of fm1"); rug(weights(fm1))

##
fmc <- lm   (y ~ poly(a,2)-a + poly(tf, 2)-tf + poly(A, 2)-A + . , data = d1k27)
summary(fmc)
## -> has NA's for  'a, tf, A'  --- bad that it did *not* work to remove them

nform <- update(formula(fm1), ~ .
                +poly(A,2)  -A  -I(A^2)
                +poly(a,2)  -a  -I(a^2)
                +poly(tf,2) -tf -I(tf^2))

fm1. <- lmrob(nform, data = d1k27)# now w/o warning !? !!
fm2. <- lmrob(nform, data = d1k27, cov = ".vcov.w", trace = TRUE)


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
