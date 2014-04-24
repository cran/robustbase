nlrob <-
    function (formula, data, start, lower, upper,
              weights = NULL, na.action = na.fail,
	      method = c("M", "MM", "tau", "CM", "mtl"),
	      psi = .Mwgt.psi1("huber", cc=1.345),
	      test.vec = c("resid", "coef", "w"),
	      maxit = 20, tol = 1e-06, acc,
	      algorithm = "default", doCov = FALSE,
	      control = if(method == "M") nls.control() else
			nlrob.control(method, optArgs = list(trace=trace), ...),
              trace = FALSE, ...)
{
    ## Purpose:
    ##	Robust fitting of nonlinear regression models. The fitting is
    ##	done by iterated reweighted least squares (IWLS) as in rlm() of
    ##	the package MASS. In addition, see also 'nls'.
    ##
    ## --> see the help file,  ?nlrob  (or ../man/nlrob.Rd in the source)
    ## -------------------------------------------------------------------------

    ##- some checks
    call <- match.call() # << and more as in nls()
    formula <- as.formula(formula)
    if (length(formula) != 3)
	stop("'formula' should be a formula of the type 'y  ~ f(x, alpha)'")
    ## Had 'acc'; now use 'tol' which is more universal; 'acc' should work for a while
    if(!missing(acc) && is.numeric(acc)) {
        if(!missing(tol)) stop("specifying both 'acc' and 'tol' is invalid")
        tol <- acc
        message("The argument 'acc' has been renamed to 'tol'; do adapt your code.")
    }
    method <- match.arg(method)
    dataName <- substitute(data)
    dataCl <- attr(attr(call, "terms"), "dataClasses")

    if(method != "M") {
      if(!is.null(weights))
          stop("specifying 'weights' is not yet supported for method ", method)
      if(!missing(start))
	  warning("Starting values will not be used for method ", method)
      force(control)
      fixAns <- function(mod) {
          mod$call <- call # the nlrob() one, not nlrob.<foo>()
          mod$data <- dataName
          mod$dataClasses <- dataCl
          mod
      }
      switch(method,
	     "MM" = {
		 return(fixAns(nlrob.MM (formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "tau" = {
		 return(fixAns(nlrob.tau(formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "CM"	 = {
		 return(fixAns(nlrob.CM (formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     },
	     "mtl" = {
		 return(fixAns(nlrob.mtl(formula, data, lower=lower, upper=upper,
					 tol=tol, ctrl= control)))
	     })
    }

    ## else: method == "M", original method, the only one based on 'nls' :
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    test.vec <- match.arg(test.vec)
    if(missing(lower)) lower <- -Inf
    if(missing(upper)) upper <- +Inf
    ## FIXME:  nls() allows  a missing 'start'; we don't :
    if(length(pnames <- names(start)) != length(start))
        stop("'start' must be fully named (list or numeric vector)")
    if (!((is.list(start) && all(sapply(start, is.numeric))) ||
	  (is.vector(start) && is.numeric(start))))
	stop("'start' must be a named list or numeric vector")
    if(any(is.na(match(pnames, varNames))))
	stop("parameter names must appear in 'formula'")
    nm <- "._nlrob.w"
    if (nm %in% c(varNames, pnames, names(data)))
	stop(gettextf("Do not use '%s' as a variable name or as a parameter name",
		      nm), domain=NA)
    if (!is.null(weights)) {
	if (length(weights) != nobs)
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
    fit <- eval(formula[[3]], c(data, start))
    y <- eval(formula[[2]], data)
    nobs <- length(y)
    resid <- y - fit
    w <- rep.int(1, nobs)
    if (!is.null(weights))
	w <- w * weights
    ##- robust loop (IWLS)
    converged <- FALSE
    status <- "converged"
    method.exit <- FALSE
    for (iiter in seq_len(maxit)) {
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
	    w <- psi(resid/Scale)
	    if (!is.null(weights))
		w <- w * weights
	    data$._nlrob.w <- w ## use a variable name the user "will not" use
	    ._nlrob.w <- NULL # workaround for codetools "bug"
            ## Case distinction against "wrong warning" as long as
            ## we don't require R > 3.0.2:
	    if(identical(lower, -Inf) && identical(upper, Inf))
	    out <- nls(formula, data = data, start = start,
		       algorithm = algorithm, trace = trace,
		       weights = ._nlrob.w,
		       na.action = na.action, control = control)
	    else
	    out <- nls(formula, data = data, start = start,
		       algorithm = algorithm, trace = trace,
		       lower=lower, upper=upper,
		       weights = ._nlrob.w,
		       na.action = na.action, control = control)

	    ## same sequence as in start! Ok for test.vec:
	    coef <- coefficients(out)
	    start <- coef
            resid <- residuals(out)
	    convi <- irls.delta(previous, get(test.vec))
	}
	converged <- convi <= tol
	if (converged)
	    break
    }

    if(!converged || method.exit) {
	warning(status <- paste("failed to converge in", maxit, "steps"))
        if(method.exit) converged <- FALSE
    }

    if(!is.null(weights)) { ## or just   out$weights  ??
	tmp <- weights != 0
	w[tmp] <- w[tmp]/weights[tmp]
    }

    ## --- Estimated asymptotic covariance of the robust estimator
    rw <- psi(res.sc <- resid/Scale)
    asCov <- if(!converged || !doCov) NA else {
        ## a version of  .vcov.m(.) below
	AtWAinv <- chol2inv(out$m$Rmat())
	dimnames(AtWAinv) <- list(names(coef), names(coef))
	tau <- mean(rw^2) / mean(psi(res.sc, d=TRUE))^2
	AtWAinv * Scale^2 * tau
    }

    ## returned object:	 ==  out$m$fitted()  [FIXME?]
    fit <- eval(formula[[3]], c(data, coef))
    names(fit) <- obsNames
    structure(class = c("nlrob", "nls"),
	      list(m = out$m, call = call, formula = formula,
		   new.formula = formula, nobs = nobs,
		   coefficients = coef,
		   working.residuals = as.vector(resid),
		   fitted.values = fit, residuals = y - fit,
		   Scale = Scale, w = w, rweights = rw,
		   cov=asCov, status = status, iter=iiter,
		   psi = psi, data = dataName, dataClasses = dataCl,
		   control = control))
}

.vcov.m <- function(m, nms.coef, psi, Scale, resid, res.sc = resid/Scale) {
    AtWAinv <- chol2inv(m$Rmat())
    stopifnot(length(Scale) == 1, Scale >= 0,
              is.character(nms.coef), length(nms.coef) == nrow(AtWAinv))
    dimnames(AtWAinv) <- list(nms.coef, nms.coef)
    rw <- psi(res.sc)
    tau <- mean(rw^2) / mean(psi(res.sc, d=TRUE))^2
    AtWAinv * Scale^2 * tau
}


## The 'nls' method is *not* correct
formula.nlrob <- function(x, ...) x$formula

sigma.nlrob <- function(object, ...)
    if(!is.null(s <- object$Scale)) s else object$coefficients[["sigma"]]

estimethod <- function(object, ...) UseMethod("estimethod")
estimethod.nlrob <- function(object, ...)
    if(is.list(object$m) && inherits(object, "nls")) "M" else object$ctrl$method

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
    eval(formula(object)[[3]], c(as.list(newdata), coef(object)))
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


residuals.nlrob <- function (object, type = c("response", "working", "pearson"), ...)
{
    type <- match.arg(type)
    R <- switch(type,
                "pearson"=
            {
                stop("type 'pearson' is not yet implemented")
                ## as.vector(object$working.residuals)
            },
                "working"=
            {
                object$working.residuals
            },
                "response"=
            {
                object$residuals
            },
                stop("invalid 'type'"))# ==> programming error, as we use match.arg()
    if (!is.null(object$na.action))
        R <- naresid(object$na.action, R)
    ## FIXME: add 'names'!
    ##MM no labels; residuals.glm() does neither: attr(val, "label") <- "Residuals"
    R
}


vcov.nlrob <- function (object, ...) {
    if(!is.na(cv <- object$cov)) cv
    else {
        sc <- object$Scale
        .vcov.m(object$m, names(coef(object)), psi=object$psi, Scale= sc,
                res.sc = object$working.residuals / sc)
    }
}

summary.nlrob <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    w <- object$w ## weights * rweights, scaled such that sum(w)=1
    n <- sum(w > 0)
    param <- coef(object)
    p <- length(param)
    rdf <- n - p
    ans <- object[c("formula", "residuals", "Scale", "w", "rweights", "cov",
		    "call", "status", "iter", "control")]
    conv <- ans$status == "converged"
    sc   <- ans$Scale
    if(is.na(ans[["cov"]]) && conv)
	ans$cov <- .vcov.m(object$m, names(param), psi=object$psi, Scale= sc,
			   res.sc = object$working.residuals / sc)
    ans$df <- c(p, rdf)
    cf <-
	if(conv) {
	    se <- sqrt(diag(ans$cov))
	    tval <- param/se
	    cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
	} else cbind(param, NA, NA, NA)
    dimnames(cf) <- list(names(param),
			 c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$coefficients <- cf
    if(correlation && rdf > 0 && conv) {
	ans$correlation <- ans$cov / outer(se, se)
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
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
	"\n\n", sep = "")
    ## cat("\nFormula: ")
    ## cat(paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
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
	summarizeRobWeights(x$rweights, digits = digits, ...)
    }
    else
	cat("** IRWLS iterations did *not* converge!\n\n")
    invisible(x)
}
