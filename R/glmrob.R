glmrob <-
function (formula, family, data, weights, subset,
	  na.action, start = NULL, offset, method = "Mqle",
	  weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL,
	  model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    if (is.character(family))
	family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
	family <- family()
    if (is.null(family$family)) {
	print(family)
	stop("'family' not recognized")
    }
    if (!(family$family %in% c("binomial", "poisson"))) {
	if(family$family == "gaussian") {
	    if(weights.on.x != "none")
		stop("Use lmrob(formula, ...) for the Gaussian case;\n",
		     " 'weights.on.x' needs to be changed")
	    return(lmrob(formula, data= data, weights= weights, subset= subset,
			 na.action= na.action, offset = offset,
			 control = control, model = model, x = x, y = y,
			 contrasts = contrasts, ...) )
	}
	else
	    stop("Robust fitting method not yet implemented for this family")
    }
    if(is.null(control)) # -> use e.g., glmrobMqle.control()
	control <- get(paste("glmrob", method, ".control", sep = ""))(...)
    if (missing(data))
	data <- environment(formula)
    ##
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
	       names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if(method == "model.frame") return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "numeric")
    X <- if (!is.empty.model(mt))
	model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0))
	stop("'weights' must be non-negative")
    if (!is.null(offset) && length(offset) != NROW(Y))
	stop("Number of offsets is ", length(offset), ", should equal ",
	     NROW(Y), " (number of observations)")
    weights.on.x <- match.arg(weights.on.x)
    fit <- switch(method,
		  "cubif" =
		  glmrobCubif(X = X, y = Y, weights = weights, start = start,
			      offset = offset, family = family,
			      weights.on.x = weights.on.x, control = control,
			      intercept = attr(mt, "intercept") > 0),
		  "Mqle" =
		  glmrobMqle(X = X, y = Y, weights = weights, start = start,
			     offset = offset, family = family,
			     weights.on.x = weights.on.x, control = control,
			     intercept = attr(mt, "intercept") > 0),
		  stop("invalid 'method': ", method))
    ##-	    if (any(offset) && attr(mt, "intercept") > 0) {
    ##-		fit$null.deviance <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
    ##-		    y = Y, weights = weights, offset = offset, family = family,
    ##-		    control = control, intercept = TRUE)$deviance
    ##-	    }
    fit$na.action <- attr(mf, "na.action")
    if (model)
	fit$model <- mf
    if (x)
	fit$x <- X
    if (!y)
	fit$y <- NULL
    fit <- c(fit,
	     list(call = call, formula = formula, terms = mt, data = data,
		  offset = offset, control = control, method = method,
		  contrasts = attr(X, "contrasts"),
		  xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("glmrob", "glm")
    fit
}



summary.glmrob <- function(object, correlation=FALSE, symbolic.cor=FALSE, ...)
{
    dispersion <- 1
    coefs <- object$coefficients
    aliased <- is.na(coefs)# needs care; also used in print method
    if(any(aliased))
	coefs <- coefs[!aliased]
    covmat <- object$cov
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    zvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind("Estimate" = coefs, "Std. Error" = s.err,
			"z-value" = zvalue, "Pr(>|z|)" = pvalue)

    ans <- c(object[c("call", "terms", "family", "iter", "control", "method",
		      "residuals", "fitted.values", "w.r", "w.x")],
	     ## MM: should rather keep more from 'object' ?
	     ##	    currently, cannot even print the asympt.efficiency!
	     list(deviance=NULL, df.residual=NULL, null.deviance=NULL,
		  df.null= NULL, df= NULL, ## (because of 0 weights; hmm,...)
		  aliased = aliased,
		  coefficients = coef.table, dispersion = dispersion,
		  cov.unscaled = covmat))
    if (correlation) {
	ans$correlation <- cov2cor(covmat)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glmrob"
    return(ans)
}

## almost a copy of vcov.glm() [if that didn't have summmary.glm() explicitly]
vcov.glmrob <- function (object, ...)
{
    so <- summary(object, corr = FALSE, ...)
    so$dispersion * so$cov.unscaled
}


print.glmrob <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(coef(x))) {
	cat("Coefficients")
	if (is.character(co <- x$contrasts))
	    cat("  [contrasts: ", apply(cbind(names(co), co),
					1, paste, collapse = "="), "]")
	cat(":\n")
	print.default(format(x$coefficients, digits = digits),
		      print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nNumber of observations:", length(x$residuals),
	"\nFitted by method ", sQuote(x$method), "\n")
    invisible(x)
}

print.summary.glmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")
    if (length(cf <- coef(x))) {
	if(nsingular <- sum(x$aliased)) # glm has   df[3] - df[1]
	    cat("\nCoefficients: (", nsingular,
		" not defined because of singularities)\n", sep = "")
	else cat("\nCoefficients:\n")
	printCoefmat(cf, digits = digits, signif.stars = signif.stars,
		     na.print = "NA", ...)

	summarizeRobWeights(x$w.r * x$w.x, digits = digits,
			    header = "Robustness weights w.r * w.x:", ...)
    }
    else cat("No coefficients\n\n")

    n <- length(x$residuals)
    cat("\nNumber of observations:", n,
	"\nFitted by method", sQuote(x$method)," (in", x$iter, "iterations)\n")

    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",sep = "")
    if(any(!is.null(unlist(x[c("null.deviance", "deviance")]))))
	cat(apply(cbind(paste(format(c("Null", "Residual"), justify="right"),
			      "deviance:"),
			format(unlist(x[c("null.deviance", "deviance")]),
			       digits=max(5, digits + 1)), " on",
			format(unlist(x[c("df.null", "df.residual")])),
			" degrees of freedom\n"), 1, paste, collapse=" "),
	    "\n", sep = "")
    else
	cat("No deviance values available \n")
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if (isTRUE(symbolic.cor)) {
		print(symnum(correl, abbr.col=NULL))
	    }
	    else {
		correl <- format(round(correl, 2), nsmall=2, digits=digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote=FALSE)
	    }
	}
    }

    printControl(x$control, digits = digits)

    cat("\n")
    invisible(x)
}
