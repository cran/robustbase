## Left Medcouple
lmc <- function(x, na.rm = FALSE, ...) {
    -mc(x[x <= median(x, na.rm = na.rm)], na.rm = na.rm, ...)
}

## Right Medcouple
rmc <- function(x, na.rm = FALSE, ...) {
    mc(x[x >= median(x, na.rm = na.rm)], na.rm = na.rm, ...)
}


## Generic function
mc <- function (x, ...)
      UseMethod("mc")

## Default method (for numeric vectors):
mc.default <- function(x, na.rm = FALSE,
		       doReflect = (length(x) <= 100),
## experiment:                       eps1 = 1e-13, eps2 = eps1, maxit = 100,
                       eps1 = .Machine$double.xmin, eps2 = eps1, maxit = 100,
                       trace.lev = 0, full.result = FALSE,
                       ...)
{
    x <- as.numeric(x)
    ina <- is.na(x)
    if (na.rm)
        x <- x[!ina]
    else if (any(ina))
        return(as.numeric(NA))

    if(length(l.. <- list(...)))
        stop("In mc(): invalid argument(s) : ",
             paste(sQuote(names(l..)), collapse=","), call. = FALSE)
    rr <- mcComp(x, doReflect, eps1=eps1, eps2=eps2,
                 maxit=maxit, trace.lev = trace.lev)

    if(!(conv1 <- rr[["converged"]]) |
       (doReflect && !(conv2 <- rr[["converged2"]]))) {
        stop("mc(): not 'converged' ",
             if(!conv1) paste("in", rr[["iter"]], "iterations"),
             if(doReflect && !conv2)
             paste(if(!conv1)" *and*",
                   "'reflect part' in", rr[["iter2"]], "iterations"))
    }

    m <- if (doReflect) (rr[["medc"]] - rr[["medc2"]]) / 2  else rr[["medc"]]
    structure(m, mcComp = if(full.result) rr)
}

mcComp <- function(x, doReflect, eps1 = 1e-13, eps2 = eps1, maxit = 1000,
                   trace.lev = 1)

{
    stopifnot(is.logical(doReflect), length(doReflect) == 1,
              is.numeric(eps1), length(eps1) == 1, eps1 >= 0,
              is.numeric(eps2), length(eps2) == 1, eps2 >= 0,
              length(maxit     <- as.integer(maxit)) == 1,
              length(trace.lev <- as.integer(trace.lev)) == 1
              )
    x <- as.numeric(x)
    n <- as.integer(length(x))
    eps <- as.double(c(eps1, eps2))
    c.iter <- c(maxit, trace.lev)
    ans <- .C(mc_C, x, n,
              eps = eps, iter = c.iter,
              medc = double(1))[c("medc", "eps", "iter")]
    it <- ans[["iter"]]
    ans[["converged"]] <- it[2] == 1
    ans[["iter"]] <- it[1]

    if (doReflect) { ## also compute on reflected data
        a2 <- .C(mc_C, (max(x) - x), n,
                 eps2 = eps, iter2 = c.iter,
                 medc2 = double(1))[c("medc2", "eps2", "iter2")]
        it <- a2[["iter2"]]
        a2[["converged2"]] <- it[2] == 1
        a2[["iter2"]] <- it[1]

	c(ans, a2)
    }
    else ans
}

