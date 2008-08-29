
#### Testing  medcouple	 and related functions

### here, we do "strict tests" -- hence no *.Rout.save
### hence, can also produce non-reproducible output such as timing

library(robustbase)
source(system.file("mcnaive.R", package = "robustbase"))# mcNaive()

allEQ <- function(x,y) all.equal(x,y, tol = 1e-12)
DO <- function(...) system.time(stopifnot(...))

DO(0 == sapply(1:100, function(n) mc(seq_len(n))))
DO(0 == sapply(1:100, function(n) mc(seq_len(n), doRefl=FALSE)))


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
	  allEQ( 7/12, mc( x4, doRefl= FALSE, eps1=.Machine$double.xmin)),
	  allEQ(-7/12, mc(-x4, doRefl= FALSE, eps1=.Machine$double.xmin)))


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


DO(0 == sapply(1:100, function(n)
   mc(seq_len(n), doRefl=FALSE, eps1=.Machine$double.xmin)))

###----  Strict tests of adjOutlyingness():

## works for this (!) seed:
set.seed(1); system.time(a1 <- adjOutlyingness(longley))
set.seed(2); system.time(a2 <- adjOutlyingness(hbk))
set.seed(3); system.time(a3 <- adjOutlyingness(hbk[, 1:3]))# the 'X' space
set.seed(4); system.time(a4 <- adjOutlyingness(milk))
set.seed(5); system.time(a5 <- adjOutlyingness(wood))
set.seed(6); system.time(a6 <- adjOutlyingness(wood[, 1:5]))# the 'X' space

## FIXME:  32-bit <-> 64-bit different results {tested on Linux only}
is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()["sysname"] == "Darwin"
stopifnot(which(!a2$nonOut) == 1:14,
	  which(!a3$nonOut) == 1:14,
	  which(!a4$nonOut) == if(is32 && !isMac) c(1, 2, 41, 70) else c(12, 70),
	  ## 'longley', 'wood' have no outliers in the "adjOut" sense:
	  if(isMac) TRUE else a1$nonOut,
          a5$nonOut, a6$nonOut,
          ## milk (n = 86) :
	  if(is32 && !isMac) ## FIXME: This is platform (32 <-> 64) dependent!
          rank(a4$adjout) ==
          c(83, 85, 59, 62, 11,   26, 27, 15, 43, 24,   73, 82, 78, 79, 81,
            76, 77, 63, 72, 68,   30, 11, 36, 18, 56,  47.5, 51, 65, 49, 14,
            42, 55,  6, 16, 22,   41, 40, 29, 11, 53,   84, 67, 46, 80, 11,
            11, 75, 70, 69, 64,   52, 66, 35,  5,  3,    1, 33, 23, 47.5, 17,
            4, 50, 38.5,38.5,31,  20,  7, 57, 37, 86,   34, 25, 44, 71, 74,
            21, 58,  2, 28, 32,    8, 19, 60, 61, 45,   54) else TRUE,
          ## hbk (n = 75) :
          rank(a3$adjout) ==
          c(62, 64, 68, 71, 70,   65, 66, 63, 69, 67,   73, 75, 72, 74, 18,
            52, 44,  4, 12, 24,    6, 24, 15, 24, 59,   14, 24, 16, 45, 39,
            49, 33,  9, 54, 24,    2, 24, 50, 56, 10,   32, 41, 43, 37, 60,
            36, 61, 24, 13, 11,   48, 55, 47, 42, 17,   30, 51, 24,  7, 38,
            24, 58, 40, 24, 24,   34,  3, 53, 57,  5,    1,  8, 31, 35, 46)
          )



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

proc.time()
