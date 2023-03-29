
#### Testing  medcouple mc() and related functions

### here, we do "strict tests" -- hence no *.Rout.save
### hence, can also produce non-reproducible output such as timing

library(robustbase)
for(f in system.file("xtraR", c("mcnaive.R", # -> mcNaive()
                                "styleData.R",  # -> smallD  list of small datasets
			      "platform-sessionInfo.R"), # -> moreSessionInfo()
                     package = "robustbase", mustWork=TRUE)) {
    cat("source(",f,"):\n", sep="")
    source(f)
}

## instead of relying on  system.file("test-tools-1.R", package="Matrix"):
source(system.file("xtraR/test-tools.R", package = "robustbase")) # assert.EQ() etc

assertEQm12 <- function(x,y, giveRE=TRUE, ...)
    assert.EQ(x,y, tol = 1e-12, giveRE=giveRE, ...)
## ^^ shows *any* difference ("tol = 0") unless there is no difference at all
##
c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
DO <- function(...) S.time(stopifnot(...))

## from {sfsmisc}:
lseq <- function(from, to, length) exp(seq(log(from), log(to), length.out = length))


mS <- moreSessionInfo(print.=TRUE)

(doExtras <- robustbase:::doExtras())# TRUE if interactive() or activated by envvar

if(!dev.interactive(orNone=TRUE)) pdf("mc-strict.pdf")

tools::assertCondition(mc(1:11), "message") # change of default to  doScale=FALSE

smlMC  <- vapply(smallD, mc, pi)
smlMCo <- vapply(smallD, mc, pi, doScale=TRUE, c.huberize=Inf)
yI <- c("yI", "yI."); notI  <- setdiff(names(smallD), yI)
yI2 <- c(yI, "x3I");  notI2 <- setdiff(names(smallD), yI2)
assert.EQ(smlMC [notI],
          smlMCo[notI], tol = 4e-11, giveRE=TRUE)
## above small diff. is from 'x3I';  dropping that, too, leaves no differences
table(smlMC [notI2] == smlMCo[notI2])

n.set <- c(1:99, 1e5L+ 0:1) # large n gave integer overflow in earlier versions
DO(0 == sapply(n.set, function(n) mc(seq_len(n))))
DO(0 == sapply(n.set, function(n) mc(seq_len(n), doRefl=FALSE)))

DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "simple")))
DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "h.use" )))


x1 <- c(1, 2, 7, 9, 10)
mcNaive(x1) # = -1/3
assertEQm12(-1/3, mcNaive(x1, "simple"))
assertEQm12(-1/3, mcNaive(x1, "h.use"))
assertEQm12(-1/3, mc(x1))

x2 <- c(-1, 0, 0, 0, 1, 2)
mcNaive(x2, meth="simple") # = 0 - which is wrong
mcNaive(x2, meth="h.use")  # = 1/6 = 0.16666
assertEQm12(1/6, mc(x2))
assertEQm12(1/6, mcNaive(x2, "h.use"))

x4 <- c(1:5,7,10,15,25, 1e15) ## - bombed in original algo
mcNaive(x4,"h.use") # 0.5833333
assertEQm12( 7/12, mcNaive(x4, "h.use"))
assertEQm12( 7/12, mcNaive(x4, "simple"))
assertEQm12( 7/12, mc( x4, doRefl= FALSE))
assertEQm12(-7/12, mc(-x4, doRefl= FALSE))

xx <- c(-3, -3, -2, -2, -1, rep(0, 6), 1, 1, 1, 2, 2, 3, 3, 5)
stopifnot(exprs = {
    mc(xx, doScale=TRUE , c.huberize = Inf) == 0 ## old mc()
    mc(xx) == 0
    mc(xx,  doReflect=FALSE) == 0
   -mc(-xx, doReflect=FALSE) == 0
    mcNaive(xx, "h.use" ) == 0
    mcNaive(xx, "simple") == 0
})

set.seed(17)
for(n in 3:50) {
    cat(" ")
    for(k in 1:5) {
	x <- rlnorm(n)
	mc1 <- mc(x)
	mc2 <- mcNaive(x, method = "simple")
	mc3 <- mcNaive(x, method = "h.use" )
	stopifnot(all.equal(mc1, mc3, tolerance = 1e-10),# 1e-12 not quite ok
		  mc2 == mc3)
	cat(".")
    }
};  cat("\n")


###----  Strict tests of adjOutlyingness():
###                      ================= changed after long-standing bug fix in Oct.2014

## For longley, note n < 4p and hence "random nonsense" numbers
set.seed(1);  S.time(a1.1 <- adjOutlyingness(longley))
set.seed(11); S.time(a1.2 <- adjOutlyingness(longley))
##
set.seed(2); S.time(a2 <- adjOutlyingness(hbk)) # 75 x 4
set.seed(3); S.time(a3 <- adjOutlyingness(hbk[, 1:3]))# the 'X' space
set.seed(4); S.time(a4 <- adjOutlyingness(milk)) # obs.63 = obs.64
set.seed(5); S.time(a5 <- adjOutlyingness(wood)) # 20 x 6  ==> n < 4p
set.seed(6); S.time(a6 <- adjOutlyingness(wood[, 1:5]))# ('X' space) 20 x 5: n = 4p (ok!)

## 32-bit <-> 64-bit different results {tested on Linux only}
is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()[["sysname"]] == "Darwin"
isSun <- Sys.info()[["sysname"]] == "SunOS"
Rnk <- function(u) rank(unname(u), ties.method = "first")
## to use for testing below:
cat("\nRnk(a3 $ adjout): "); dput(Rnk(a3$adjout), control= {})
cat("\nRnk(a4 $ adjout): "); dput(Rnk(a4$adjout), control= {})

(i.a4Out <- which( ! a4$nonOut)) # the outliers -- varies "wildly"
stopifnot(70 %in% i.a4Out)
{
    if(is32 && !isMac)
        all.equal(i.a4Out, c(1, 2, 41, 70))
    ## and this is "typically" true, but not for a 64-bit Linux version bypassing BLAS in matprod
    else if(isSun || isMac)
        TRUE
    else if(grepl("^Fedora", osVersion) && !is32)
        identical(i.a4Out, 70L) # since Dec 2020 (F 32)
    else
        all.equal(i.a4Out, c(9:19, 23:27,57, 59, 70, 77))
}

## only for ATLAS (BLAS/Lapack), not all are TRUE; which ones [but n < 4p]
if(!all(a5$nonOut))
  print(which(!a5$nonOut)) # if we know, enable check below

stopifnot(exprs = {
    which(!a2$nonOut) == 1:14
    which(!a3$nonOut) == 1:14
    ## 'longley', 'wood' have no outliers in the "adjOut" sense:
    if(doExtras && !isMac) { ## longley also has n < 4p (!)
        if(mS$ strictR)
            sum(a1.2$nonOut) >= 15 # sum(.) = 16 [nb-mm3, Oct.2014]
        else ## however, openBLAS Fedora Linux /usr/bin/R gives sum(a1.2$nonOut) = 13
            sum(a1.2$nonOut) >= 13
    } else TRUE
    if(doExtras) { ## have n < 4p (!)
        if(mS$ strictR) a5$nonOut
        else ## not for ATLAS
            sum(a5$nonOut) >= 18 # 18: OpenBLAS
    } else TRUE
    a6$nonOut[-20]
    ## hbk (n = 75, p = 3) should be "stable" (but isn't quite)
    abs(Rnk(a3$adjout) -
        c(62, 64, 69, 71, 70,    66, 65, 63, 68, 67,    73, 75, 72, 74, 35,
          60, 55,  4, 22, 36,     6, 33, 34, 28, 53,    16, 13,  9, 27, 31,
          49, 39, 20, 50, 14,     2, 24, 40, 54, 21,    17, 37, 52, 23, 58,
          19, 61, 11, 25,  8,    46, 59, 48, 47, 29,    44, 43, 42,  7, 30,
          18, 51, 41, 15, 10,    38,  3, 56, 57,  5,     1, 12, 26, 32, 45)
        ) <= 3 ## all 0 on 64-bit (F 32) Linux
})

## milk (n = 86) : -- Quite platform dependent!
r <- Rnk(a4$adjout)
r64 <- ## the 64-bit (ubuntu 14.04, nb-mm3) values:
    c(65, 66, 61, 56, 47,   51, 19, 37, 74, 67,   79, 86, 83, 84, 85,
      82, 81, 73, 80, 55,   27,  3, 70, 68, 78,   76, 77, 53, 48,  8,
      29, 33,	 6, 32, 28,   31, 36, 40, 22, 58,   64, 52, 39, 63, 44,
      30, 57, 46, 43, 45,   25, 54, 12,  1,  9,    2, 71, 14, 75, 23,
      4, 10, 34, 35, 17,   24, 15, 20, 38, 72,   42, 13, 50, 60, 62,
      26, 69, 18,  5, 21,    7, 49, 11, 41, 59,   16)
r32 <- ## Linux 32bit (florence: 3.14.8-100.fc19.i686.PAE)
    c(78, 79, 72, 66, 52,   61, 22, 41, 53, 14,   74, 85, 82, 83, 84,
      80, 81, 56, 73, 65,   30,  3, 16, 17, 68,   57, 58, 63, 54,  8,
      32, 37,  6, 36, 31,   35, 40, 44, 25, 69,   77, 62, 43, 76, 48,
      34, 67, 51, 47, 49,   28, 64, 12,  1,  9,    2, 33, 15, 59, 26,
      4, 10, 38, 39, 20,   27, 18, 23, 42, 86,   46, 13, 60, 71, 75,
      29, 50, 21,  5, 24,    7, 55, 11, 45, 70,   19)
d <- (r - if (is32) r32 else r64)
cbind(r, d)
       table(abs(d))
cumsum(table(abs(d))) # <=> unscaled ecdf(d)

## For the biggest part (79 out of 86), the ranks are "close":
## 2014: still true, but in a different sense..
##       ^ typically, but e.g., *not* when using non-BLAS matprod():
sum(abs(d) <= 17) >= 78
sum(abs(d) <= 13) >= 75


## check of adjOutlyingness *free* bug
## reported by Kaveh Vakili <Kaveh.Vakili@wis.kuleuven.be>
set.seed(-37665251)
X <- matrix(rnorm(100*5),   100, 5)
Z <- matrix(rnorm(10*5)/100, 10, 5)
Z[,1] <- Z[,1] + 5
X[91:100,] <- Z # if anything these should be outliers, but ...
for (i in 1:10) {
    ## this would produce an error in the 6th iteration
    aa <- adjOutlyingness(x=X, ndir=250)
    if(any(is.out <- !aa$nonOut))
        cat("'outliers' at obs.", paste(which(is.out), collapse=", "),"\n")
    stopifnot(1/4 < aa$adjout & aa$adjout < 16)
}

## Check "high"-dimensional Noise ... typically mc() did *not* converge for some re-centered columns
## Example by Valentin Todorov:
n <- 50
p <- 30
set.seed(1) # MM
a <- matrix(rnorm(n * p), nrow=n, ncol=p)
str(a)
kappa(a) # 20.42 (~ 10--20 or so; definitely not close to singular)
a.a <- adjOutlyingness(a, mcScale=FALSE, # <- my own recommendation
                       trace.lev=1)
a.s <- adjOutlyingness(a, mcScale=TRUE, trace.lev=1)
## a.a :
str(a.a) # high 'adjout' values "all similar" -> no outliers .. hmm .. ???
(hdOut <- which( ! a.a$nonOut)) ## indices of "outlier" -- very platform dependent !
a.a$MCadjout; all.equal(a.a$MCadjout, 0.136839766177,
          tol = 1e-12) # seen 7.65e-14   and "big" differences on non-default platforms
## a.s :
which(! a.s$nonOut ) # none  [all TRUE]
a.s$MCadjout # platform dependent; saw
all.equal(a.s$MCadjout, 0.32284906741568, tol = 1e-13) # seen 2.2e-15 ..
                                        # and big diffs on non-default platforms
##
## The adjout values are all > 10^15  !!!   why ??
## Now (2021) I know: n < 4*p ==> can find 1D-projection where 1 of the 2 {Q3-Q2, Q2-Q1} is 0 !
##---------------------------------------------------------------------------------------------


###-- Back to  mc()  checks for "hard" cases
###           =====  -----------------------

## "large n" (this did overflow sum_p, sum_q  earlier ==> had inf.loop):
set.seed(3); x <- rnorm(2e5)
(mx <- mc(x, trace.lev=3))
stopifnot(print(abs(mx - -0.000772315846101988)) < 1e-15)
					# 3.252e-19, 64b Linux
					# 1.198e-16, 32b Windows

### Some platform info :
local({ nms <- names(Si <- Sys.info())
        dropNms <- c("nodename", "machine", "login")
        structure(Si[c("nodename", nms[is.na(match(nms, dropNms))])],
                  class="simple.list") })

if(identical(1L, grep("linux", R.version[["os"]]))) { ##----- Linux - only ----
    ##
    Sys.procinfo <- function(procfile)
    {
        l2 <- strsplit(readLines(procfile),"[ \t]*:[ \t]*")
        r <- sapply(l2[sapply(l2, length) == 2],
                    function(c2)structure(c2[2], names= c2[1]))
        attr(r,"Name") <- procfile
        class(r) <- "simple.list"
        r
    }
    ##
    Scpu <- Sys.procinfo("/proc/cpuinfo")
    Smem <- Sys.procinfo("/proc/meminfo")
    print(Scpu[c("model name", "cpu MHz", "cache size", "bogomips")])
    print(Smem[c("MemTotal", "SwapTotal")])
}

##' Checking the breakdown point of mc() --- Hubert et al. theory said : 25%
##' using non-default  doReflect=FALSE  as that corresponds to original Hubert et al.
##'
##' @title Medcouple mc() checking
##' @param x
##' @param Xfun
##' @param eps
##' @param NAiferror
##' @param doReflect
##' @param ...
##' @return mc(*,..) or NaN in case mc() signals an error [non-convergence]
##' @author Martin Maechler
mcX <- function(x, Xfun, eps=0, NAiferror=FALSE, doReflect=FALSE, ...) {
    stopifnot(is.numeric(x), is.function(Xfun), "eps" %in% names(formals(Xfun)))
    myFun <-
	if(NAiferror)
	    function(u) tryCatch(mc(Xfun(u, eps=eps), doReflect=doReflect, ...),
				 error = function(e) NaN)
	else
	    function(u) mc(Xfun(u, eps=eps), doReflect=doReflect, ...)
    vapply(x, myFun, 1.)
}

X1. <- function(u, eps=0) c(1,2,3, 7+(-10:10)*eps, u + (-1:1)*eps)
## ==> This *did* breakdown [but points not "in general position"]:
## but now is stable:
r.mc1 <- curve(mcX(x, X1.), 10, 1e35, log="x", n=1001)
stopifnot(r.mc1$y == 0) # now stable
if(FALSE) {
rt1 <- uniroot(function(x) mcX(exp(x), X1.) - 1/2, lower=0, upper=500)
exp(rt1$root) #  4.056265e+31
}

## eps > 0  ==> No duplicated points ==> theory says breakdown point = 0.25
## -------  but get big numerical problems:
if(FALSE) { # ==> convergence problem [also in maxit = 1e5] .. really an *inf* loop!
r.mc1.1  <- curve(mcX(x, X1., eps= .1  ), 10, 1e35, log="x", n=1001)
r.mc1.2  <- curve(mcX(x, X1., eps= .01 ), 10, 1e35, log="x", n=1001)
r.mc1.3  <- curve(mcX(x, X1., eps= .001), 10, 1e35, log="x", n=1001)
r.mc1.5  <- curve(mcX(x, X1., eps= 1e-5), 10, 1e35, log="x", n=1001)
r.mc1.8  <- curve(mcX(x, X1., eps= 1e-8), 10, 1e35, log="x", n=1001)
r.mc1.15 <- curve(mcX(x, X1., eps=1e-15), 10, 1e35, log="x", n=1001)# still!
}
## practically identical to  eps = 0 where we have breakdown (see above)
r.mc1.16 <- curve(mcX(x, X1., eps=1e-16), 10, 1e35, log="x", n=1001)
all.equal(r.mc1, r.mc1.16, tol=1e-15)#-> TRUE

## Quite bad case: Non convergence
X2. <- function(u) c(1:3, seq(6, 8, by = 1/8), u, u, u)
## try(mc(X2.(4.3e31)))## -> error: no convergence
## but now
stopifnot(exprs = {
    all.equal(1/30, mc(X2.(4.3e31)), tol=1e-12)
    all.equal(1/30, mc(X2.(4.3e31), eps1=1e-7, eps2=1e-100), tol=1e-12)
})
## related, more direct:
X3. <- function(u) c(10*(1:3), 60:80, (4:6)*u)
stopifnot(0 == mc(X3.(1e31), trace=5)) # fine convergence in one iter.
stopifnot(0 == mc(X3.(1e32), trace=3)) # did *not* converge

### TODO : find example with *smaller* sample size -- with no convergence
X4. <- function(u, eps, ...) c(10, 70:75, (2:3)*u)
mc(X4.(1e34))# "fine"
## now stable too:
r.mc4 <- curve(mcX(x, X4.), 100, 1e35, log="x", n=2^12)
stopifnot(abs(1/3 - r.mc4$y) < 1e-15)

X5. <- function(u) c(10*(1:3), 70:78, (4:6)*u)
stopifnot(all.equal(4/15, mc(X5.(1e32), maxit=1000)))

X5. <- function(u, eps,...) c(5*(1:12), (4:6)*u)
str(r.mc5 <- mc(X5.(1e32), doReflect=FALSE, full.result = TRUE))
## Now, stable:
stopifnot(all.equal(1/5, c(r.mc5))) ## was 1; platform dependent ..
stopifnot(all.equal(4/15, mc(X5.(5e31)))) # had  no convergence w/ maxit=10000
r.mc5Sml <- curve(mcX(x, X5.), 1,  100, log="x", n=1024) ## quite astonishing
x <- lseq(1, 1e200, 2^11)
mc5L <- mcX(x, X5.)
table(err <- abs(0.2 - mc5L[x >= 24])) # I see all 0!
stopifnot(abs(err) < 1e-15)

c.time(proc.time())

summary(warnings()) # seen 15 x  In mcComp(....) :
## maximal number of iterations (100 =? 100) reached prematurely

