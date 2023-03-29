
## Left Medcouple
lmc <- function(x, mx = median(x, na.rm=na.rm), na.rm = FALSE, doReflect = FALSE, ...) {
    -mc(x[x <= mx], na.rm=na.rm, doReflect=doReflect, ...)
}

## Right Medcouple
rmc <- function(x, mx = median(x, na.rm=na.rm), na.rm = FALSE, doReflect = FALSE, ...) {
    mc(x[x >= mx], na.rm=na.rm, doReflect=doReflect, ...)
}

.optEnv$mc_doScale_msg <- TRUE # initially
mc <- function(x, na.rm = FALSE, doReflect = (length(x) <= 100)
             , doScale = FALSE # was hardwired=TRUE (till 2018-07); then default=TRUE (till 2022-03)
             , c.huberize = 1e11 # was implicitly = Inf originally
             , eps1 = 1e-14, eps2 = 1e-15 # << new in 0.93-2 (2018-07..)
             , maxit = 100, trace.lev = 0
             , full.result = FALSE
               )
{
    x <- as.numeric(x)
    ina <- is.na(x)
    if (na.rm)
        x <- x[!ina]
    else if (any(ina))
        return(NA_real_)
    ## ==> x is NA-free from here on
    stopifnot(length(c.huberize) == 1L, c.huberize >= 0)

    ## For robustbase 0.95-0 (March 2022); drop the message eventually:
    if(missing(doScale) && .optEnv$mc_doScale_msg && !getOption("mc_doScale_quiet", FALSE)) {
        message("The default of 'doScale' is FALSE now for stability;\n  ",
                "set options(mc_doScale_quiet=TRUE) to suppress this (once per session) message")
        .optEnv$mc_doScale_msg <- FALSE
    }

    if(is.finite(c.huberize))
        x <- huberize(x, c=c.huberize, warn0 = trace.lev > 0, saveTrim=FALSE)
    rr <- mcComp(x, doReflect, doScale=doScale, eps1=eps1, eps2=eps2,
                 maxit=maxit, trace.lev=trace.lev)

    if(!(conv1 <- rr[["converged"]]) |
       (doReflect && !(conv2 <- rr[["converged2"]]))) {
        stop("mc(): not 'converged' ",
             if(!conv1) paste("in", rr[["iter"]], "iterations"),
             if(doReflect && !conv2)
             paste(if(!conv1)" *and*",
                   "'reflect part' in", rr[["iter2"]], "iterations"),
             "; try enlarging eps1, eps2 !?\n")
    }

    m <- if (doReflect) (rr[["medc"]] - rr[["medc2"]]) / 2  else rr[["medc"]]
    structure(m, mcComp = if(full.result) rr)
}

## eps1 = 1e-13, eps2 = eps1  <==>  original code which only had  'eps = 1e-13'
##                                  hardcoded in C code.
## These defaults do *not* make sense here, but in mc().
## However, currently they are used in ../tests/mc-etc.R
mcComp <- function(x, doReflect, doScale, eps1, eps2, maxit = 1000, trace.lev = 1)
{
    stopifnot(is.logical(doReflect), length(doReflect) == 1L, !is.na(doReflect),
              is.logical(doScale),   length(doScale)   == 1L, !is.na(doScale),
              is.1num(eps1), eps1 >= 0,
              is.1num(eps2), eps2 >= 0,
              length(maxit     <- as.integer(maxit)) == 1,
              length(trace.lev <- as.integer(trace.lev)) == 1
              )
    ## Assumption [from caller, = mc()]: 'x' has no NAs (but can have +-Inf)
    x <- as.numeric(x)
    n <- as.integer(length(x))
    eps <- as.double(c(eps1, eps2))
    c.iter <- c(maxit, trace.lev)
    ## NAOK=TRUE: to allow  +/- Inf to be passed
    ans <- .C(mc_C, x, n,
              eps = eps, iter = c.iter,
              medc = double(1)
            , doScale = doScale
            , NAOK=TRUE)[c("medc", "eps", "iter")]
    it <- ans[["iter"]]
    ans[["converged"]] <- it[2] == 1
    ans[["iter"]] <- it[1]

    if (doReflect) { ## also compute on reflected data
	a2 <- .C(mc_C, -x, n,
                 eps2 = eps, iter2 = c.iter,
                 medc2 = double(1)
               , doScale = doScale
               , NAOK=TRUE)[c("medc2", "iter2", "doScale")]
        it <- a2[["iter2"]]
        a2[["converged2"]] <- it[2] == 1
        a2[["iter2"]] <- it[1]

	c(ans, a2)
    }
    else ans
}

