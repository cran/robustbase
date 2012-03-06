
#### Testing  medcouple	 and related functions

### here, we do "strict tests" -- hence no *.Rout.save
### hence, can also produce non-reproducible output such as timing

library(robustbase)
source(system.file("xtraR/mcnaive.R", package = "robustbase"))# mcNaive()

allEQ <- function(x,y) all.equal(x,y, tol = 1e-12)
##
c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
DO <- function(...) S.time(stopifnot(...))

n.set <- c(1:99, 1e5L+ 0:1) # large n gave integer overflow in earlier versions
DO(0 == sapply(n.set, function(n) mc(seq_len(n))))
DO(0 == sapply(n.set, function(n) mc(seq_len(n), doRefl=FALSE)))

DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "simple")))
DO(0 == sapply(1:100, function(n) mcNaive(seq_len(n), "h.use" )))


x1 <- c(1, 2, 7, 9, 10)
mcNaive(x1) # = -1/3
stopifnot(allEQ(-1/3, mcNaive(x1)),
	  allEQ(-1/3, mcNaive(x1, "h.use")),
	  allEQ(-1/3, mc(x1)))

x2 <- c(-1, 0, 0, 0, 1, 2)
mcNaive(x2, meth="simple") # = 0 - which is wrong
mcNaive(x2, meth="h.use")  # = 1/6 = 0.16666
stopifnot(allEQ(1/6, mc(x2)),
	  allEQ(1/6, mcNaive(x2, "h.use")))

x4 <- c(1:5,7,10,15,25, 1e15) ## - bombed in orignal algo
mcNaive(x4,"h.use") # 0.5833333
stopifnot(allEQ( 7/12, mcNaive(x4, "h.use")),
	  allEQ( 7/12, mc( x4, doRefl= FALSE)),
	  allEQ(-7/12, mc(-x4, doRefl= FALSE)))


set.seed(17)
for(n in 3:50) {
    cat(" ")
    for(k in 1:5) {
	x <- rlnorm(n)
	mc1 <- mc(x)
	mc2 <- mcNaive(x, method = "simple")
	mc3 <- mcNaive(x, method = "h.use" )
	stopifnot(all.equal(mc1, mc3, tol = 1e-10),# 1e-12 not quite ok
		  mc2 == mc3)
	cat(".")
    }
};  cat("\n")

###----  Strict tests of adjOutlyingness():

set.seed(1);  S.time(a1.1 <- adjOutlyingness(longley))
set.seed(11); S.time(a1.2 <- adjOutlyingness(longley))
##
set.seed(2); S.time(a2 <- adjOutlyingness(hbk))
set.seed(3); S.time(a3 <- adjOutlyingness(hbk[, 1:3]))# the 'X' space
set.seed(4); S.time(a4 <- adjOutlyingness(milk))
set.seed(5); S.time(a5 <- adjOutlyingness(wood))
set.seed(6); S.time(a6 <- adjOutlyingness(wood[, 1:5]))# the 'X' space

## 32-bit <-> 64-bit different results {tested on Linux only}
is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()["sysname"] == "Darwin"
isSun <- Sys.info()["sysname"] == "SunOS"
Rnk <- function(u) rank(unname(u), ties.method = "first")
## for later testing:
dput(Rnk(a3$adjout),, {})
dput(Rnk(a4$adjout),, {})

stopifnot(which(!a2$nonOut) == 1:14,
	  which(!a3$nonOut) == 1:14,
	  if(isSun) TRUE else
	  which(!a4$nonOut) == if(is32 && !isMac) c(1, 2, 41, 70) else c(12, 70),
	  ## 'longley', 'wood' have no outliers in the "adjOut" sense:
	  ## FIXME: longley is platform dependent too
	  if(isMac) TRUE else sum(a1.2$nonOut) >= 15,
	  a5$nonOut, a6$nonOut,
	  ## hbk (n = 75) :
	  abs(Rnk(a3$adjout) -
	     c(62, 64, 68, 71, 70,   65, 66, 63, 69, 67,   73, 75, 72, 74, 18,
	       52, 44,	4, 12, 19,    6, 20, 15, 21, 59,   14, 22, 16, 45, 39,
	       49, 33,	9, 54, 23,    2, 24, 50, 56, 10,   32, 40, 43, 39, 60,
	       36, 61, 25, 13, 11,   48, 55, 47, 42, 17,   30, 51, 26,	7, 37,
	       27, 58, 39, 28, 29,   34,  3, 53, 57,  5,    1,	8, 31, 35, 46)
	      ) <= 3,
	  ## milk (n = 86) : -- Quite platform dependent!
      {
	  r <- Rnk(a4$adjout)
	  r64 <- ## the 64-bit (FC13 Linux) values:
	      c(79, 83, 54, 57,	 9,   20, 27, 15, 59, 23,   74, 85, 80, 82, 84,
		78, 77, 58, 75, 64,   30, 10, 36, 35, 72,   69, 70, 61, 46, 14,
		41, 51,	 6, 16, 26,   42, 40, 29, 11, 49,   81, 63, 45, 76, 12,
		13, 73, 66, 65, 60,   48, 62, 33,  5,  3,    1, 37, 22, 68, 17,
		 4, 47, 38, 39, 31,   19,  7, 52, 34, 86,   24, 25, 43, 67, 71,
		21, 53, 2, 28, 32,     8, 18, 55, 56, 44,   50)
	  ## for the biggest part (79 out of 86), the ranks are "close":
	  table(d <- (r - r64)[-c(9, 24:27, 59, 71)]) # 32b Linux: 0 .. 6
	  abs(d) <= 7
      })


## check of adjOutlyingness *free* bug
## reported by Kaveh Vakili <Kaveh.Vakili@wis.kuleuven.be>
set.seed(-37665251)
X<-matrix(rnorm(100*5),100,5)
Z<-matrix(rnorm(100*5,0,1/100),10,5)
Z<-sweep(Z,2,c(5,rep(0,4)),FUN="+")
X[91:100,]<-Z
for (i in 1:10) {
    ## this would produce an error in the 6th iteration
    aa <- adjOutlyingness(x=X,ndir=250)
}

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

c.time(proc.time())
