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
    function(formula, data, subset, weights, na.action, method = 'MM',
             model = TRUE, x = !control$compute.rd, y = FALSE,
             singular.ok = TRUE, contrasts = NULL, offset = NULL,
             control = NULL, ...)
{
    ## to avoid problems with setting argument
    ## call lmrob.control here either with or without method arg.
    if (missing(control))
        control <- if (missing(method))
            lmrob.control(...) else lmrob.control(method = method, ...)
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
    if (!missing(control) && !missing(method) && method != control$method) {
        warning("Methods argument set by method is different from method in control\n",
                "Using method = ", method)
        control$method <- method
    }

    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(,0,3) else numeric(0),
                  residuals = y, fitted.values = 0 * y,
                  cov = matrix(,0,0), weights = w, rank = 0,
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

        z <- lmrob.fit(x, y, control)
    }

    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    if(control$compute.rd && !is.null(x))
        z$MD <- robMD(x, attr(mt, "intercept"))
    if (model)
        z$model <- mf
    if (ret.x)
        z$x <- x
    if (ret.y)
        z$y <- y
    z
}

## internal function, used in lmrob() and maybe plot.lmrob()
robMD <- function(x, intercept, ...) {
    if(intercept == 1) x <- x[, -1, drop=FALSE]
    if(ncol(x) >= 1) {
	rob <- covMcd(x, ...)
	sqrt( mahalanobis(x, rob$center, rob$cov) )
    }
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


vcov.lmrob <- function (object, cov=object$control$cov, ...) {
  if (!is.null(object$cov) && identical(cov, object$control$cov))
    return(object$cov)
  else {
    lf.cov <- if (!is.function(cov))
      get(cov, mode = "function")
    else cov
    lf.cov(object, ...)
  }
}


## residuals.default works for "lmrob"  {unless we'd allow re-weighted residuals
## fitted.default works for "lmrob"

## This seems to work - via lm
model.matrix.lmrob <- function (object, ...) {
    stats::model.matrix.lm(object, ...)
}

if(FALSE) ## now replaced with more sophsticated in ./lmrobPredict.R
## learned from MASS::rlm() : via "lm" as well
predict.lmrob <- function (object, newdata = NULL, scale = NULL, ...)
{
    class(object) <- c(class(object), "lm")
    object$qr <- qr(sqrt(object$weights) * object$x)
    predict.lm(object, newdata = newdata, scale = object$s, ...)
}


summary.lmrob <- function(object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    if (is.null(object$terms))
	stop("invalid 'lmrob' object:  no terms component")
    p <- object$rank
    df <- object$degree.freedom
    if (p > 0) {
	n <- p + df
	se <- sqrt(diag(object$cov))
	est <- object$coefficients
	tval <- est/se
	ans <- object[c("call", "terms", "residuals", "scale", "weights",
			"converged", "iter", "control")]
	ans$df <- df
	ans$coefficients <-
	    if( ans$converged )
		cbind(est, se, tval, 2 * pt(abs(tval), df, lower.tail = FALSE))
	    else cbind(est, NA, NA, NA)
	dimnames(ans$coefficients) <-
	    list(names(est), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$cov.unscaled <- object$cov
	dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
	if (correlation) {
	    ans$correlation <- ans$cov.unscaled / outer(se, se)
	    ans$symbolic.cor <- symbolic.cor
	}
    } else { ## p = 0: "null model"
	ans <- object
	ans$coefficients <- matrix(, 0, 4)
	ans$df <- df
	ans$cov.unscaled <- object$cov
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
			print(symnum(correl), abbr.colnames = NULL)
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

    if (x$control$method == 'SM') x$control$method <- 'MM'
    printControl(x$control, digits = digits)

    invisible(x)
}

## hidden in namespace
printControl <-
    function(ctrl, digits = getOption("digits"),
	     str.names = "seed",
	     header = "Algorithmic parameters:",
	     ...)
{
    ## Purpose: nicely and sensibly print a 'control' structure
    ##		currently  for lmrob(), glmrob()
    ## Author: Martin Maechler, Date: 31 May 2006
    PR <- function(LST, ...) if(length(LST)) print(unlist(LST), ...)

    cat(header,"\n")
    is.str <- (nc <- names(ctrl)) %in% str.names
    is.ch <- sapply(ctrl, is.character)
    real.ctrl <- sapply(ctrl, function(x)
			length(x) > 0 && is.numeric(x) && x != round(x))
    PR(ctrl[!is.str & real.ctrl], digits = digits, ...)
    ## non-real, non-char ones (typically integers), but dropping 0-length ones
    PR(ctrl[!is.str & !is.ch & !real.ctrl], ...)
    ## char ones
    PR(ctrl[!is.str & is.ch], ...)
    if(any(is.str))
	for(n in nc[is.str]) {
	    cat(n,":")
	    str(ctrl[[n]], vec.len = 2)
	    ## 'vec.len = 2' is smaller than normal, but nice for Mersenne seed
	}
}

summarizeRobWeights <-
    function(w, digits = getOption("digits"), header = "Robustness weights:",
	     eps = 0.1 / length(w), eps1 = 1e-3, ...)
{
    ## Purpose: nicely print a "summary" of robustness weights
    stopifnot(is.numeric(w))
    cat(header,"\n")
    cat0 <- function(...) cat('', ...)
    n <- length(w)
    if(n <= 10) print(w, digits = digits, ...)
    else {
	n1 <- sum(w1 <- abs(w - 1) < eps1)
	n0 <- sum(w0 <- abs(w) < eps)
        if(any(w0 & w1))
            warning("weights should not be both close to 0 and close to 1!\n",
                    "You should use different 'eps' and/or 'eps1'")
	if(n0 > 0 || n1 > 0) {
	    if(n0 > 0) {
                formE <- function(e) formatC(e, digits = max(2, digits-3), width=1)
		i0 <- which(w0)
		maxw <- max(w[w0])
		c3 <- paste("with |weight| ",
			    if(maxw == 0) "= 0" else paste("<=", formE(maxw)),
			    " ( < ", formE(eps), ");", sep='')
		cat0(if(n0 > 1) {
		       cc <- sprintf("%d observations c(%s)",
				     n0, strwrap(paste(i0, collapse=",")))
		       c2 <- " are outliers"
		       paste(cc,
			     if(nchar(cc)+ nchar(c2)+ nchar(c3) > getOption("width"))
			     "\n	", c2, sep='')
		     } else
		       sprintf("observation %d is an outlier", i0),
		     c3, "\n")
	    }
	    if(n1 > 0)
		cat0(ngettext(n1, "one weight is",
			     sprintf("%s%d weights are",
				     if(n1 == n)"All " else '', n1)), "~= 1.")
	    n.rem <- n - n0 - n1
	    if(n.rem <= 0) { # < 0 possible if w0 & w1 overlap
		if(n1 > 0) cat("\n")
		return(invisible())
	    }
	    cat0("The remaining",
		 ngettext(n.rem, "one", sprintf("%d ones", n.rem)), "are")
	    if(is.null(names(w)))
		names(w) <- as.character(seq(along = w))
	    w <- w[!w1 & !w0]
	    if(n.rem <= 10) {
		cat("\n")
		print(w, digits = digits, ...)
		return(invisible())
	    }
	    else cat(" summarized as\n")
	}
	print(summary(w, digits = digits), digits = digits, ...)
    }
}

