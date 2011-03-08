nlrob <-
    function (formula, data, start, weights = NULL, na.action = na.fail,
	      psi = psi.huber, test.vec = c("resid", "coef", "w"),
	      maxit = 20, acc = 1e-06, algorithm = "default",
	      control = nls.control(), trace = FALSE, ...)
{
    ## Purpose:
    ##	Robust fitting of nonlinear regression models. The fitting is
    ##	done by iterated reweighted least squares (IWLS) as in rlm() of
    ##	the package MASS. In addition, see also 'nls'.
    ##
    ## --> see the help file,  ?nlrob  (or ../man/nlrob.Rd in the source)
    ## -------------------------------------------------------------------------

    ##- some checks
    mf <- match.call() # << and more as in nls()  [FIXME or drop]
    formula <- as.formula(formula)
    if (length(formula) != 3)
	stop("'formula' should be a formula of the type 'y  ~ f(x, alpha)'")
    test.vec <- match.arg(test.vec)
    varNames <- all.vars(formula)
    dataName <- substitute(data)
    data <- as.data.frame(data)

    ## FIXME:  nls() allows  a missing 'start';	 we don't :
    if (length(pnames <- names(start)) != length(start))
	stop("'start' must be fully named (list or numeric vector)")
    if (!((is.list(start) && all(sapply(start, is.numeric))) ||
	  (is.vector(start) && is.numeric(start)))
	|| any(is.na(match(pnames, varNames))))
	stop("'start' must be a list or numeric vector named with parameters in 'formula'")
    if ("w" %in% varNames || "w" %in% pnames || "w" %in% names(data))
	stop("Do not use 'w' as a variable name or as a parameter name")
    if (!is.null(weights)) {
	if (length(weights) != nrow(data))
	    stop("'length(weights)' must equal the number of observations")
	if (any(weights < 0) || any(is.na(weights)))
	    stop("'weights' must be nonnegative and not contain NAs")
    }

    ## if (any(is.na(data)) & options("na.action")$na.action == "na.omit")
    ##	 stop("if NAs are present, use 'na.exclude' to preserve the residuals length")

    irls.delta <- function(old, new) sqrt(sum((old - new)^2, na.rm = TRUE)/
					  max(1e-20, sum(old^2, na.rm = TRUE)))

    ##- initialize testvec  and update formula with robust weights
    coef <- start
    fit <- eval(formula[[3]], c(as.list(data), start))
    y <- eval(formula[[2]], as.list(data))
    resid <- y - fit
    w <- rep(1, nrow(data))
    if (!is.null(weights))
	w <- w * weights
    ## The following "put everything on the right in order to use weights"
    ## is only necessary for R versions <= 2.2.1  (FIXME eventually)
    oform <- formula
    formula <- as.formula(substitute(~(LHS-RHS) * w, list(LHS = formula[[2]],
							  RHS = formula[[3]])))
    ##- robust loop (IWLS)
    converged <- FALSE
    status <- "converged"
    method.exit <- FALSE
    for (iiter in 1:maxit) {
	if (trace)
	    cat("robust iteration", iiter, "\n")
	previous <- get(test.vec)
	Scale <- median(abs(resid), na.rm = TRUE)/0.6745
	if (Scale == 0) {
	    convi <- 0
	    method.exit <- TRUE
	    warning(status <- "could not compute scale of residuals")
	    ## FIXME : rather use a "better" Scale in this case, e.g.,
	    ## -----   Scale <- min(abs(resid)[resid != 0])
	}
	else {
	    w <- psi(resid/Scale, ...)
	    if (!is.null(weights))
		w <- w * weights
	    data$w <- sqrt(w)
	    out <- nls(formula, data = data, start = start, algorithm = algorithm,
		       trace = trace, na.action = na.action, control = control)

	    ## same sequence as in start! Ok for test.vec:
	    coef <- coefficients(out)
	    start <- coef
	    resid <- -residuals(out)/sqrt(w) ## == - (y - f(x))*sqrt(w)
	    convi <- irls.delta(previous, get(test.vec))
	}
	converged <- convi <= acc
	if (converged)
	    break
    }

    if(!converged && !method.exit)
	warning(status <- paste("failed to converge in", maxit, "steps"))

    if(!is.null(weights)) {
	tmp <- weights != 0
	w[tmp] <- w[tmp]/weights[tmp]
    }

    ## --- Estimated asymptotic covariance of the robust estimator
   if (!converged && !method.exit) {
        asCov <- NA
    } else {
        AtWAinv <- chol2inv(out$m$Rmat())
        dimnames(AtWAinv) <- list(names(coef), names(coef))
        tau <- (mean(psi(resid/Scale, ...)^2) /
                (mean(psi(resid/Scale, deriv=TRUE, ...))^2))
        asCov <- AtWAinv * Scale^2 * tau
    }

    ## returned object:
    out <- list(m = out$m, call = match.call(), formula = oform,
                new.formula = formula,
		coefficients = coef, residuals = resid,
		fitted.values = y - out$residuals,
		Scale = Scale, w = w, w.r = psi(resid/Scale, ...),
                cov=asCov, status = status, iter=iiter,
                psi = psi, data = dataName,
		dataClasses = attr(attr(mf, "terms"), "dataClasses"))
    ##MM: Where would this "label" really make sense?
    ##MM: attr(out$fitted.values, "label") <- "Fitted values"
    ##-	  names(out$residuals) <- rownames(data)
    ##-	  names(out$fitted.values) <- rownames(data)
    class(out) <- c("nlrob", "nls")
    out
}


fitted.nlrob <- function (object, ...)
{
    val <- as.vector(object$fitted.values)
    if (!is.null(object$na.action))
	val <- napredict(object$na.action, val)
    ##MM: attr(val, "label") <- "Fitted values"
    val
}


## formula() works "by default"

predict.nlrob <- function (object, newdata, ...)
{
    if (missing(newdata))
	return(as.vector(fitted(object)))
    if (!is.null(cl <- object$dataClasses))
	.checkMFClasses(cl, newdata)
    eval(object$formula[[3]], c(as.list(newdata), coef(object)))
}


print.nlrob <- function (x, ...)
{
    cat("Robustly fitted nonlinear regression model\n")
    cat("  model: ", deparse(formula(x)), "\n")
    cat("   data: ", deparse(x$data), "\n")
    print(coef(x), ...)
    cat(" status: ", x$status, "\n")
    invisible(x)
}


residuals.nlrob <- function (object, ...)
{
    val <- as.vector(object$residuals)
    if (!is.null(object$na.action))
	val <- naresid(object$na.action, val)
    ##MM: attr(val, "label") <- "Residuals"
    val
}


summary.nlrob <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    w <- object$w ## weights * w.r, scaled such that sum(w)=1
    n <- sum(w > 0)
    param <- coef(object)
    p <- length(param)
    rdf <- n - p
    ans <- object[c("formula", "residuals", "Scale", "w", "w.r", "cov",
		    "call", "status", "iter", "control")]
    ans$df <- c(p, rdf)
    cf <-
	if(ans$status == "converged") {
	    se <- sqrt(diag(object$cov))
	    tval <- param/se
	    cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
	} else cbind(param, NA, NA, NA)
    dimnames(cf) <- list(names(param),
			 c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$coefficients <- cf
    if(correlation && rdf > 0 && ans$status == "converged") {
	ans$correlation <- object$cov / outer(se, se)
	ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.nlrob"
    ans
}

print.summary.nlrob <-
    function (x, digits = max(3, getOption("digits") - 3),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nFormula: ")
    cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
    df <- x$df
    rdf <- df[2L]
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
		 ...)
    if(x$status == "converged") {
	cat("\nRobust residual standard error:",
	    format(signif(x$Scale, digits)), "\n")
	correl <- x$correlation
	if (!is.null(correl)) {
	    p <- NCOL(correl)
	    if (p > 1) {
		cat("\nCorrelation of Parameter Estimates:\n")
		if(is.logical(symbolic.cor) && symbolic.cor) {
		    print(symnum(correl, abbr.colnames = NULL))
		} else {
		    correl <- format(round(correl, 2), nsmall = 2, digits = digits)
		    correl[!lower.tri(correl)] <- ""
		    print(correl[-1, -p, drop=FALSE], quote = FALSE)
		}
	    }
	}
	cat("Convergence in", x$iter, "IRWLS iterations\n\n")
	summarizeRobWeights(x$w.r, digits = digits, ...)
    }
    else
	cat("** IRWLS iterations did *not* converge!\n\n")
    invisible(x)
}
