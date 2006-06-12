### FIXME:
### ----- MM wants to change

### 3) allow the 'control' entries to enter via "..." as well -- Done

### 4) lmrob() should really behave like lm() {in R; not S !}
###		--> 'subset' etc			      -- Done

### 5) There are still quite a few things hard-coded in ../src/lmrob.c
###    E.g., 'nResample' is used, but MAX_NO_RESAMPLES = 500 cannot be changed.

### 6) Use ' method = "MM" ' and a general scheme for "plugin" of other estimators!!


### The first part of lmrob()  much cut'n'paste from lm() - on purpose!
lmrob <-
    function(formula, data, subset, weights, na.action,
	     model = TRUE, x = !control$compute.rd, y = FALSE,
	     singular.ok = TRUE, contrasts = NULL, offset = NULL,
	     control = lmrob.control(...), ...)
{
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms") # allow model.frame to update it
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != NROW(y))
	stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
		      length(offset), NROW(y)), domain = NA)

    if (is.empty.model(mt)) {
	x <- NULL
	z <- list(coefficients = if (is.matrix(y))
		    matrix(,0,3) else numeric(0), residuals = y,
		  fitted.values = 0 * y, weights = w, rank = 0,
		  df.residual = NROW(y), converged = TRUE)
	if(!is.null(offset)) z$fitted.values <- offset
    }
    else {
	x <- model.matrix(mt, mf, contrasts)
	if (!singular.ok)
	    warning("only 'singular.ok = TRUE' is currently implemented.")
	if (!is.null(w))
	    stop("Weights are not yet implemented for this estimator")
	if(!is.null(offset))
	    stop("'offset' not yet implemented for this estimator")
	z <- lmrob.fit.MM(x, y, control = control)
    }

    class(z) <- "lmrob"
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if( control$compute.rd && !is.null(x)) {
	x0 <- if(attr(mt, "intercept") == 1) x[, -1, drop=FALSE] else x
	if(ncol(x0) >= 1) {
	    rob <- covMcd(x0)
	    z$MD <- sqrt( mahalanobis(x0, rob$center, rob$cov) )
	}
    }
    if (model)
	z$model <- mf
    if (ret.x)
	z$x <- x
    if (ret.y)
	z$y <- y
    z$control <- control
    z
}


print.lmrob <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if(length((cf <- coef(x)))) {
	if( x$converged )
	    cat("Coefficients:\n")
	else
	    cat("Algorithm did not converge\n\n",
		"Coefficients of the *initial* S-estimator:\n")
	print(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}


summary.lmrob <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    z <- object
    if (is.null(z$terms))
	stop("invalid 'lmrob' object:  no terms component")
    p <- z$rank
    df <- z$degree.freedom
    if (p > 0) {
	n <- p + df
	se <- sqrt(diag(z$cov))
	est <- z$coefficients
	tval <- est/se
	ans <- z[c("call", "terms", "residuals", "scale", "weights",
		   "converged", "iter", "control")]
	ans$df <- df
	ans$coefficients <-
	    if( ans$converged )
		cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
	    else cbind(est, NA, NA, NA)
	dimnames(ans$coefficients) <-
	    list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$cov.unscaled <- z$cov
	dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
	if (correlation) {
	    ans$correlation <- ans$cov.unscaled / outer(se, se)
	    ans$symbolic.cor <- symbolic.cor
	}
    } else { ## p = 0: "null model"
	ans <- z
	ans$coefficients <- matrix(, 0, 4)
	ans$df <- df
	ans$cov.unscaled <- z$cov
    }
    class(ans) <- "summary.lmrob"
    ans
}


print.summary.lmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
	"\n\n", sep = "")
    resid <- x$residuals
    df <- x$df
    ##df
    cat(if (!is.null(x$w) && diff(range(x$w))) "Weighted ",
	"Residuals:\n", sep = "")
    if (df > 5) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	if (NCOL(resid) > 1)
	    rq <- structure(apply(t(resid), 1, quantile),
			    dimnames = list(nam, dimnames(resid)[[2]]))
	else rq <- structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else print(resid, digits = digits, ...)
    if( length(x$coef) ) {
	if( !(x$converged) ) {
	    cat("\nAlgorithm did not converge\n")
	    cat("\nCoefficients of *initial* S-estimator:\n")
	    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
			 ...)
	} else {
	    cat("\nCoefficients:\n")
	    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
			 ...)
	    cat("\nRobust residual standard error:",
		format(signif(x$scale, digits)),"\n")
	    correl <- x$correlation
	    if (!is.null(correl)) {
		p <- NCOL(correl)
		if (p > 1) {
		    cat("\nCorrelation of Coefficients:\n")
		    if (is.logical(symbolic.cor) && symbolic.cor) {
			print(symnum(correl), abbr.col = NULL)
		    }
		    else { correl <- format(round(correl, 2), nsmall = 2,
					    digits = digits)
			   correl[!lower.tri(correl)] <- ""
			   print(correl[-1, -p, drop = FALSE], quote = FALSE)
		       }
		}
	    }
	    cat("Convergence in", x$iter, "IRWLS iterations\n")
	}
	cat("\n")

        summarizeRobWeights(x$weights, digits = digits, ...)

    } else cat("\nNo Coefficients\n")

    printControl(x$control, digits = digits)

    invisible(x)
}

## hidden in namespace:
printControl <-
    function(ctrl, digits = getOption("digits"),
	     str.names = "seed",
	     header = "Algorithmic parameters:",
	     ...)
{
    ## Purpose: nicely and sensibly print a 'control' structure
    ## Author: Martin Maechler, Date: 31 May 2006
    cat(header,"\n")
    is.str <- (nc <- names(ctrl)) %in% str.names
    real.ctrl <- sapply(ctrl, function(x) length(x) > 0 && x != round(x))
    print(unlist(ctrl[!is.str & real.ctrl]), digits = digits, ...)
    ## non-real ones, but dropping 0-length ones
    print(unlist(ctrl[!is.str & !real.ctrl]), ...)
    if(any(is.str))
	for(n in nc[is.str]) {
	    cat(n,":")
	    str(ctrl[[n]], vec.len = 2)
	    ## 'vec.len = 2' is smaller than normal, but nice for Mersenne seed
	}
}

summarizeRobWeights <-
    function(w, digits = getOption("digits"),
             header = "Robustness weights:", ...)
{
    ## Purpose: nicely print a "summary" of robustness weights
    stopifnot(is.numeric(w))
    cat(header,"\n")
    n <- length(w)
    if(n <= 10) print(w, digits = digits, ...)
    else {
	n1 <- sum(w1 <- abs(w - 1) < 1e-4)
	n0 <- sum(w0 <- abs(w) < 1e-4 / n)
	if(n0 > 0 || n1 > 0) {
	    if(n0 > 0)
		cat(n0, " observations c(",
		    strwrap(paste(which(w0),collapse=",")),
		    ")\n  are outliers with |weights| < ", formatC(1e-4 / n),".\n",
		    sep='')
	    if(n1 > 0)
		cat(n1, "weights are ~= 1.\n")
	    cat("The remaining", n - n0 - n1,
		" ones are summarized as\n")
	    w <- w[!w1 & !w0]
	}
	print(summary(w, digits = digits), digits = digits, ...)
    }
}

