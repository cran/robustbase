glmrob <-
function (formula, family = binomial, data, weights, subset,
		   na.action, start = NULL, offset, method = "Mqle",
		   weights.on.x = c("none", "hat", "robCov"), control = NULL,
		   model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,	 ...)
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
	if(family$family == "gaussian")
	    stop("Use rlm(formula, ...) from library MASS for the Gaussian case")
	else
	    stop("Robust fitting method not yet implemented for this family")
    }
    if(is.null(control)) # -> use e.g., glmrobMqle.control()
	control <- get(paste("glmrob", method, ".control", sep = ""))()
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
	stop("Negative wts not allowed")
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
    if (model)
	fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
	fit$x <- X
    if (!y)
	fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
		       data = data, offset = offset, control = control, method = method,
		       contrasts = attr(X, "contrasts"),
		       xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("glmrob", "glm")
    fit
}



summary.glmrob <- function (object, correlation=FALSE, symbolic.cor=FALSE, ...)
{
    dispersion <- 1
    coefs <- object$coefficients
    covmat <- object$cov
    dimnames(covmat) <- list(names(coefs), names(coefs))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    zvalue <- coefs/s.err
    dn <- c("Estimate", "Std.Error")

    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coefs, s.err, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coefs),
				 c(dn, "z-value", "Pr(>|z|)"))
    ans <- c(object[c("call", "terms", "family", "iter", "method")],
	     ## MM: maybe we should rather keep the full object?
	     list(deviance=NULL, df.residual=NULL, null.deviance=NULL,
		  df.null=NULL, deviance.resid=NULL, coefficients=coef.table,
		  dispersion=dispersion, df=NULL , cov.unscaled=covmat))
    if (correlation) {
	ans$correlation <- cov2cor(covmat)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glmrob"
    return(ans)
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
    cat("\nNumber of oberservations:", length(x$y),
	"\nFitted by method = ", x$method, "\n")
    invisible(x)
}


print.summary.glmrob <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
		 na.print = "NA", ...)
    cat("\n(Dispersion parameter for ", x$family$family,
	" family taken to be ", format(x$dispersion), ")\n\n",sep = "")
    if(any(!is.null(unlist(x[c("null.deviance", "deviance")]))))
	cat(apply(cbind(paste(format(c("Null", "Residual"), justify="right"),
			      "deviance:"),
			format(unlist(x[c("null.deviance", "deviance")]),
			       digits=max(5, digits + 1)), " on",
			format(unlist(x[c("df.null", "df.residual")])),
			" degrees of freedom\n"), 1, paste, collapse=" "),
	    "\n\n", sep = "")
    else
	cat("No deviance values available \n\n")
    cat("Fitted by method = ", x$method, "\n",
	"Number of iterations: ", x$iter, "\n", sep = "")
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if (is.logical(symbolic.cor) && symbolic.cor) {
		print(symnum(correl, abbr.col=NULL))
	    }
	    else {
		correl <- format(round(correl, 2), nsmall=2, digits=digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote=FALSE)
	    }
	}
    }
    cat("\n")
    invisible(x)
}
