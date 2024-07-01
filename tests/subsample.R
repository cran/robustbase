### test subsample
### LU decomposition and singular subsamples handling
require(robustbase)
source(system.file("xtraR/subsample-fns.R", package = "robustbase", mustWork=TRUE))
## instead of relying on  system.file("test-tools-1.R", package="Matrix"):
source(system.file("xtraR/test-tools.R", package = "robustbase")) # assert.EQ(), showProc.time() ..
options(nwarnings = 4e4, warnPartialMatchArgs = FALSE)

cat("doExtras:", doExtras <- robustbase:::doExtras(),"\n")
showProc.time()

A <- rbind(c(0.001, 1),
           c(1,     2))
set.seed(11)
## IGNORE_RDIFF_BEGIN
sa <- tstSubsample(A) # (now typically also shows Matrix version ..)
## IGNORE_RDIFF_END
str(sa)

A <- rbind(c(3, 17, 10),
           c(2,  4, -2),
           c(6, 18, 12))
tstSubsample(A)

## test some random matrix
set.seed(1002)
A <- matrix(rnorm(100), 10)
tstSubsample(A)

## test singular matrix handling
A <- rbind(c(1, 0, 0),
           c(0, 1, 0),
           c(0, 1, 0),
           c(0, 0, 1))
tstSubsample(A)


## test subsample with mts > 0
data <- data.frame(y = rnorm(9), expand.grid(A = letters[1:3], B = letters[1:3]))
x <- model.matrix(y ~ ., data)
y <- data$y
## this should produce a warning and return status == 2
showSys.time(z <- Rsubsample(x, y, mts=2))
stopifnot(identical(2L, z$status)) # (z$status may be NULL; stopifnot(NULL) does not trigger)


## test equilibration
## columns only
X <- rbind(c(1e-7, 1e-10),
           c(2   , 0.2))
y <- 1:2
tstSubsample(t(X), y)

## rows only
X <- rbind(c(1e-7, 1e+10),
           c(2   , 0.2))
y <- 1:2
tstSubsample(X, y)

## both
X <- rbind(c(1e-7,  2  ),
           c(1e10, 2e12))
y <- 1:2
tstSubsample(X, y)
showProc.time()


## test real data example
data(possumDiv)## 151 * 9; the last two variables are factors
with(possumDiv, table(eucalyptus, aspect))

mf <- model.frame(Diversity ~ .^2, possumDiv)
X <- model.matrix(mf, possumDiv)
ncol(X) # 63
y <- model.response(mf)
stopifnot(identical(qr(X)$rank, ncol(X)))

## this used to fail: different pivots in step 37
str(s1 <- tstSubsample(X, y))
s2 <- tstSubsample(X / max(abs(X)), y / max(abs(X)))
s3 <- tstSubsample(X * 2^-50, y * 2^-50)
## all components *BUT*  x, y, lu, Dr, Dc, rowequ, colequ :
nm <- names(s1); nm <- nm[is.na(match(nm, c("x","y","lu", "Dr", "Dc", "rowequ", "colequ")))]
stopifnot(all.equal(s1[nm], s2[nm], tolerance=1e-10),
	  all.equal(s1[nm], s3[nm], tolerance=1e-10))
showProc.time()

set.seed(10)
nsing <- sum(replicate(if(doExtras) 200 else 20, tstSubsampleSing(X, y)))
stopifnot(nsing == 0)
showProc.time()

## test example with many categorical predictors - 2 different random seeds:
set.seed(10) ; r1 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none")
set.seed(108); r2 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none")# lmrob.S() failed
(i1 <- r1$init) # print(<lmrob.S>)
(i2 <- r1$init) # ... and they are "somewhat" close:
stopifnot(all.equal(r1[names(r1) != "init.S"],
                    r2[names(r2) != "init.S"], tol = 0.40))
c1 <- coef(r1)
c2 <- coef(r2)
relD <- (c1-c2)*2/(c1+c2)
xCf <- which(abs(relD) >= 10)
stopifnot(exprs = {
    identical(xCf, c(`Bark:aspectSW-NW` = 46L))
    all.equal(c1[-xCf], c2[-xCf], tol = 0.35) # 0.3418
    sign(c1[-xCf]) == sign(c2[-xCf])
})
showProc.time()

## investigate problematic subsample:
idc <- 1 + c(140, 60, 12, 13, 89, 90, 118, 80, 17, 134, 59, 94, 36,
         43, 46, 93, 107, 62, 57, 116, 11, 45, 35, 38, 120, 34, 29,
         33, 147, 105, 115, 92, 61, 91, 104, 141, 138, 129, 130, 84,
         119, 132, 6, 135, 112, 16, 67, 41, 102, 76, 111, 82, 148, 24,
         131, 10, 96, 0, 87, 21, 127, 56, 124)
rc <- lm(Diversity ~ .^2 , data = possumDiv, subset = idc)

X <- model.matrix(rc)
y <- possumDiv$Diversity[idc]
tstSubsample(X, y)## have different pivots ... could not find non-singular

lu <- LU.gaxpy(t(X))
stopifnot(length(lusi <- lu$sing) >= 1, lusi)
zc <- Rsubsample(X, y)
stopifnot(length(st <- zc$status) > 0, st > 0)
## column 52 is linearly dependent and should have been discarded
## qr(t(X))$pivot

image(as(round(zc$lu -      (lu$L + lu$U - diag(nrow(lu$U))), 10), "Matrix"))
image(as( sign(zc$lu) - sign(lu$L + lu$U - diag(nrow(lu$U))),      "Matrix"))
showProc.time()

## test equilibration
## colequ only
X <- matrix(c(1e-7, 2, 1e-10, 0.2), 2)
y <- 1:2
tstSubsample(t(X), y)

## rowequ only
X <- matrix(c(1e-7, 2, 1e10, 0.2), 2)
y <- 1:2
tstSubsample(X, y)

## both
X <- matrix(c(1e-7, 1e10, 2, 2e12), 2)
y <- 1:2
tstSubsample(X, y)
showProc.time()

### real data, see MM's ~/R/MM/Pkg-ex/robustbase/hedlmrob.R
##  close to singular cov():
attach(system.file("external", "d1k27.rda", package="robustbase", mustWork=TRUE))

fm1 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27)
##     ^^^^^ gave error, earlier, now with a warning -- use ".vcov.w"
## --> cov = ".vcov.w"
fm2 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27,
             cov = ".vcov.w", trace = TRUE)
showProc.time()# 2.77

if(doExtras) withAutoprint({##---------------------------------------------------------

## Q: does it change to use numeric instead of binary factors ?
## A: not really ..
d1k.n <- d1k27
d1k.n[-(1:5)] <- lapply(d1k27[,-(1:5)], as.numeric)

fm1.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n)
fm2.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n,
             cov = ".vcov.w", trace = 2)

summary(weights(fm1, type="robustness"))
   hist(weights(fm1, type="robustness"), main="robustness weights of fm1")
rug(weights(fm1, type="robustness"))
showProc.time()## 2.88

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

## now lmrob takes care of NA coefficients automatically
lmrob(y ~ poly(a,2)-a + poly(tf, 2)-tf + poly(A, 2)-A + . , data = d1k27)
showProc.time() ## 4.24
}) ## only if(doExtras) -----------------------------------------------------

## test exact fit property
set.seed(20)
data <- data.frame(y=c(rep.int(0, 20), round(100*rnorm(5))),
                   group=rep(letters[1:5], each=5))
x <- model.matrix(y ~ group, data)
(ini <- lmrob.S(x, data$y, lmrob.control()))
(ret <- lmrob(y ~ group, data))
summary(ret)
showProc.time() ## 4.24

##--- continuous x -- exact fit -- inspired by Thomas Mang's real data example
mkD9 <- function(iN, dN = 1:m) {
    stopifnot((length(iN) -> m) == length(dN), 1 <= m, m <= 5,
              iN == as.integer(iN), is.numeric(dN), !is.na(dN))
    x <- c(-3:0,0:1,1:3) # {n=9; sorted; x= 0, 1  are "doubled"}
    y <- x+5
    y[iN] <- y[iN] + dN
    data.frame(x,y)
}

mkRS <- function(...) { set.seed(...); .Random.seed }

d <- mkD9(c(1L,3:4, 7L))
rs2 <- mkRS(2)
Se <- tryCatch(error = identity,
               with(d, lmrob.S(cbind(1,x), y, lmrob.control("KS2014", seed=rs2))))
## gave DGELS rank error {for lmrob.c+wg..}

if(inherits(Se, "error")) {
    cat("Caught ")
    print(Se)
} else withAutoprint({ ## no error
    coef(Se)
    stopifnot(coef(Se) == c(5, 1)) # was (0 0)
    residuals(Se) # was == y  ---- FIXME
})

## try 100 different seeds
repS <- lapply(1:100, function(ii) tryCatch(error = identity,
                with(d, lmrob.S(cbind(1,x), y, lmrob.control("KS2014", seed = mkRS(ii))))))
if(FALSE)
 ## was
 str(unique(repS))## ==> 100 times the same error
## now completely different: *all* returned properly
str(cfS <- t(sapply(repS, coef))) # all numeric -- not *one* error --
## even all the *same* (5 1) solution:
(ucfS <- unique(cfS))
stopifnot(identical(ucfS, array(c(5, 1), dim = 1:2, dimnames = list(NULL, c("", "x")))))

## *Not* "KS2014" but the defaults works *all the time* (!)
repS0 <- lapply(1:100, function(ii) tryCatch(error = identity,
               with(d, lmrob.S(cbind(1,x), y, lmrob.control(seed = mkRS(ii))))))
summary(warnings())
## 100 identical warnings:
## In lmrob.S(cbind(1, x), y, lmrob.control(seed = mkRS(ii))) :
##   S-estimated scale == 0:  Probably exact fit; check your data

str(cfS0 <- t(sapply(repS0, coef))) # all numeric -- not *one* error
## even all the same *and* the same as "KS2014"
(ucfS0 <- unique(cfS0))
stopifnot(nrow(ucfS0) == 1L,
          ucfS0 == c(5,1))



d9L <- list(
    mkD9(c(1L,3L, 5L, 7L))
  , mkD9(c(1L,3L, 8:9))
  , mkD9(2L*(1:4))
)

if(doExtras) {
sfsmisc::mult.fig(length(d9L)); invisible(lapply(d9L, function(d) plot(y ~ x, data=d)))
}

dorob <- function(dat, control=lmrob.control(...), meth = c("S", "MM"),
                  doPl=interactive(), cex=1, ...) {
    meth <- match.arg(meth)
    stopifnot(is.data.frame(dat), c("x","y") %in% names(dat), is.list(control))
    if(doPl) plot(y ~ x, data=dat) ## with(dat, n.plot(x, y, cex=cex))
    ans <- tryCatch(error = identity,
                  switch(meth
                       , "S" = with(dat, lmrob.S(cbind(1,x), y, control))
                       , "MM"= lmrob(y ~ x, data = dat, control=control)
                       , stop("invalid 'meth'")))
    if(!doPl)
        return(ans)
    ## else
    if(!inherits(ans, "error")) {
        abline(coef(ans))
    } else { # error
        mtext(paste(paste0("lmrob.", meth), "Error:", conditionMessage(ans)))
    }
    invisible(ans)
}

## a bad case -- much better new robustbase >= 0.99-0
Se <- dorob(d9L[[1]], lmrob.control("KS2014", mkRS(2), trace.lev=4))
## was really bad -- ended returning  coef = (0 0); fitted == 0, residuals == 0 !!

if(doExtras) sfsmisc::mult.fig(length(d9L))
r0 <- lapply(d9L, dorob, seed=rs2, doPl=doExtras) # 3 x ".. exact fit" warning
if(doExtras) print(r0)
## back to 3 identical fits: (5 1)
(cf0 <- sapply(r0, coef))
stopifnot(cf0 == c(5,1))

if(doExtras) sfsmisc::mult.fig(length(d9L))
### Here, all 3 were "0-models"
r14 <- lapply(d9L, dorob, control=lmrob.control("KS2014", seed=rs2), doPl=doExtras)
## --> 3 (identical) warnings:   In lmrob.S(cbind(1, x), y, control) :#
##                        S-estimated scale == 0:  Probably exact fit; check your data
## now *does* plot
if(doExtras) print(r14)
## all 3 are "identical"
(cf14 <- sapply(r14, coef))
identical(cf0, cf14) # see TRUE; test a bit less:
stopifnot(all.equal(cf0, cf14, tol=1e-15))

## use "large n"
ctrl.LRG.n <- lmrob.control("KS2014", seed=rs2, trace.lev = if(doExtras) 2 else 1, # 3: too much (for now),
                            nResample = 60,
                            fast.s.large.n = 7, n.group = 3, groups = 2)
rLrg.n <- lapply(d9L, \(d) lmrob.S(cbind(1,d$x), d$y, ctrl.LRG.n))
summary(warnings())
sapply(rLrg.n, coef)
## currently ... ....  really would want always (5 1)
##      [,1] [,2]     [,3]
## [1,]    5    5 7.333333
## [2,]    1    1 1.666667


## ==> use  lmrob() instead of  lmrob.S():

mm0 <- lapply(d9L, dorob, meth = "MM", seed=rs2, doPl=doExtras) # looks all fine -- no longer: error in [[3]]
if(doExtras) print(mm0)
## now, the 3rd one errors (on Linux, not on M1 mac!)
(cm0 <- sapply(mm0, function(.) if(inherits(.,"error")) noquote(paste("Caught", as.character(.))) else coef(.)))

## no longer needed
c0.12 <- rbind(`(Intercept)` = c(5.7640215, 6.0267156),
               x             = c(0.85175883, 1.3823841))
if(is.list(cm0)) { ## after error {was on Linux+Win, not on M1 mac}:
    ## NB: This does *not* happen on Macbuilder -- there the result it cf = (5 1)  !!
    stopifnot(all.equal(tol = 1e-8, # seen 4.4376e-9
                        c0.12, simplify2array(cm0[1:2])))
    print(cm0[[3]])
    ## FIXME?:   Caught Error in eigen(ret, symmetric = TRUE): infinite or missing values in 'x'\n
} else if(is.matrix(cm0)) { # when no error happened
    k <- ncol(cm0)
    stopifnot(all.equal(tol = 1e-8, rbind(`(Intercept)` = rep(5,k), "x" = rep(1,k)), cm0))
} else warning("not yet encountered this case {and it should not happen}")


se3 <- lmrob(y ~ x, data=d9L[[3]], init = r0[[3]], seed=rs2, trace.lev=6)


if(doExtras) sfsmisc::mult.fig(length(d9L))
### Here, all 3 were "0-models"
##  now, have 3 *different* cases {with this seed}
## [1] : init fails (-> r14[[1]] above)
## [2] : init s=0, b=(5,1) .. but  residuals(),fitted() wrong
## [3] : init s=0, b=(5,1) ..*and* residuals(),fitted() are good

cm14 <- lapply(d9L, dorob, meth = "MM", control=lmrob.control("KS2014", seed=rs2), doPl=doExtras)
## now, first is error; for others, coef = (5, 1) are correct:
stopifnot(exprs = {
    sapply(cm14[-1], coef)  == c(5,1)
    sapply(cm14[-1], sigma) == 0
})

m2 <- cm14[[2]]
summary(m2) # prints quite nicely; and this is perfect (for scale=0), too:
## {residual != 0 <==> weights = 0}:
cbind(rwgt = weights(m2, "rob"), res = residuals(m2), fit = fitted(m2), y = d9L[[2]][,"y"])

sapply(cm14, residuals) ## now, [2] is good; [3] still wrong - FIXME
sapply(cm14, fitted)
sapply(cm14, weights, "robust")## [2]: 0 1 0 1 1 1 1 0 0;  [3]: all 0

## (unfinished ... do *test* once we've checked platform consistency)

summary(warnings())
showProc.time()

