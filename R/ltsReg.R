#### This is from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.


ltsReg <- function(x, ...) UseMethod("ltsReg")

ltsReg.formula <- function(formula, data, ...,
	model = TRUE, x.ret = FALSE, y.ret = FALSE)
{
    ##	  method <- match.arg(method)

    mf <- match.call(expand.dots = FALSE)
    mf$method <- mf$contrasts <- mf$model <- mf$x.ret <- mf$y.ret <- mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)

    ##	  if (method == "model.frame") return(mf)

    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf, contrasts)

    ##	 Check if there is an Intercept in the model - default.
    ##	 A formula without intercept looks like this: Y~.-1
    ##	 If so, remove the column named Intercept and call ltsReg with intercept=TRUE.
    ##	 Otherwise call ltsReg with intercept=FALSE
    ##
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if(xint) x <- x[, -xint, drop = FALSE]

    fit <- ltsReg(x, y, intercept = (xint > 0), ...)

    fit$terms <- mt
    fit$call <- match.call()
    fit$contrasts <- attr(x, "contrasts")
    fit$xlevels <- .getXlevels(mt, mf)

    ##	  fit$na.action <- attr(mf, "na.action")

    if(model) fit$model <- mf
    if(x.ret) fit$x <- x
    if(y.ret) fit$y <- y

    fit
}


ltsReg.default <- function (x, y,
		    intercept = TRUE,
		    alpha = NULL,
		    nsamp = 500,
		    adjust = FALSE,
		    mcd = TRUE,
		    qr.out = FALSE,
		    yname = NULL,
		    seed = 0,
		    control,
		    ...)
{

    ##	 Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    if(!missing(control)) {
	defcontrol <- rrcov.control()	# default control
	if(length(alpha) == 0 && control$alpha != defcontrol$alpha)
	    alpha <- control$alpha
	if(nsamp == defcontrol$nsamp)
	    nsamp <- control$nsamp
	if(seed == defcontrol$seed)
	    seed <- control$seed
	##	  if(print.it == defcontrol$print.it)
	##	      print.it <- control$print.it
	if(adjust == defcontrol$adjust)
	    adjust <- control$adjust
    }


    ##cat("++++++ Entering ltsReg() ...\n")
    if (is.vector(y) || (is.matrix(y) && !is.data.frame(y))) {
	if (!is.numeric(y))
	    stop("y is not a numeric dataframe or vector.")
    }
    if ((!is.matrix(y) && !is.vector(y)) || is.data.frame(y)) {
	if ((!is.data.frame(y) && !is.numeric(y)) ||
            (!all(sapply(y, data.class) == "numeric")))
	    stop("y is not a numeric dataframe or vector.")
    }

    y <- as.matrix(y)
    if (dim(y)[2] != 1)
	stop("y is not onedimensional.")

    if (missing(x)) {
	##cat("++++++ Prepare: x is missing...\n")
	x <- rep(1, nrow(y))
	if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
	    if (!is.numeric(x))
		stop("x is not a numeric dataframe or matrix.")
	}
	if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
	    if ((!is.data.frame(x) && !is.numeric(x)) ||
                (!all(sapply(x, data.class) == "numeric")))
		stop("x is not a numeric dataframe or matrix.")
	}
	if (!is.matrix(x))
	    x <- array(x, c(length(x), 1), list(names(x), deparse(substitute(x))))
	x <- as.matrix(x)
	dimny <- dimnames(y)
	dimnx <- dimnames(x)
	na.x <- !is.finite(x %*% rep(1, ncol(x)))
	na.y <- !is.finite(y)
	if (nrow(na.x) != nrow(na.y))
	    stop("Number of observations in x and y not equal")
	ok <- !(na.x | na.y)
	y <- y[ok, , drop = FALSE]
	dy <- nrow(y)
	rownames <- dimny[[1]]
	yn <- if (!is.null(yname)) yname else dimny[[2]]
	if (!length(yn))
	    yn <- "Y"
	storage.mode(y) <- "double"
	x <- x[ok, , drop = FALSE]
	storage.mode(x) <- "double"
	dx <- dim(x)
	if (!length(dx))
	    stop("All observations have missing values!")
	n <- dx[1]
    } else {

	##cat("++++++ Prepare: x is present...\n")
	if (is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
	    if (!is.numeric(x))
		stop("x is not a numeric dataframe or matrix.")
	}
	if ((!is.matrix(x) && !is.vector(x)) || is.data.frame(x)) {
	    if ((!is.data.frame(x) && !is.numeric(x)) ||
                (!all(sapply(x, data.class) == "numeric")))
		stop("x is not a numeric dataframe or matrix.")
	}

	##VT:: if the data is supplied as a data.frame, the following expressions results in an error
	## as workaround convert the data.frame to a matrix
	if(is.data.frame(x))
	    x <- as.matrix(x)
	else if(!is.matrix(x)) {
	    x <- array(x, c(length(x), 1),
		       list(names(x), deparse(substitute(data))))
	    x <- as.matrix(x)
	}
	dimny <- dimnames(y)
	dimnx <- dimnames(x)
	na.x <- !is.finite(x %*% rep(1, ncol(x)))
	na.y <- !is.finite(y)
	if (nrow(na.x) != nrow(na.y))
	    stop("Number of observations in x and y not equal")
	ok <- !(na.x | na.y)
	y <- y[ok, , drop = FALSE]
	dy <- nrow(y)
	rownames <- dimny[[1]]
	yn <- if (!is.null(yname))
	    yname
	else dimny[[2]]
	if (!length(yn))
	    yn <- "Y"
	storage.mode(y) <- "double"
	x <- x[ok, , drop = FALSE]
	storage.mode(x) <- "double"
	dx <- dim(x)
	if (!length(dx))
	    stop("All observations have missing values!")
	n <- dx[1]
	constantcolom <- function(x) {
	    c1 <- range(x)
	    c1[1] == c1[2]
	}
	if (sum(apply(x, 2, constantcolom)) > 0)
	    stop("There is at least one constant column. Remove this column and set intercept=T")
    }

    ##cat("++++++ Prepare: Ready.\n")

    dn <- dimnames(x)
    xn <- dn[[2]]
    if (!length(xn))
	if (dx[2] > 1)
	    xn <- paste("X", 1:dx[2], sep = "")
	else xn <- "X"
    X <- x
    dimnames(X) <- list(NULL, xn)
    y <- as.vector(y)


    if (all(x == 1)) {

	##cat("++++++ A - all x == 1...\n")
	if (length(alpha)) {
	    if (alpha < 1/2)
		stop("alpha is out of range!")
	    else if (alpha > 1)
		stop("alpha is greater than 1")
	    quan <- quan.f(alpha, n, dx[2])
	}
	else {
	    alpha <- 1/2
	    quan <- quan.f(alpha, n, dx[2])
	}

	p <- 1
	xbest <- NULL
	if (alpha == 1) {
	    scale <- sqrt(cov.wt(x)$cov)
	    center <- as.vector(mean(x))
	} else {
	    sh <- .fastmcd(as.matrix(y), as.integer(quan), nsamp = 0, seed)

	    y <- as.vector(y)
	    center <- as.double(sh$initmean)
	    qalpha <- qchisq(quan/n, 1)
	    calphainvers <- pgamma(qalpha/2, 1/2 + 1)/(quan/n)
	    calpha <- 1/calphainvers
	    correct <- LTScorrfactor.s(1, intercept = intercept, n, alpha)
	    scale <- sqrt(as.double(sh$initcovariance)) * sqrt(calpha) * correct
	    xbest <- sort(as.vector(sh$inbest))
	}
	resid <- y - center
	ans <- list()
	ans$best <- xbest
	ans$coefficients <- center
	ans$alpha <- alpha
	ans$quan <- quan
	ans$raw.resid <- resid/scale
	weights <- rep(NA, n)
	if (abs(scale) < 1e-07) {
	    weights <- ifelse(abs(resid) < 1e-07, 1, 0)
	    ans$scale <- ans$raw.scale <- 0
	    ans$crit <- 0
	    ans$coefficients <- ans$raw.coefficients <- center
	}
	if (abs(scale) >= 1e-07) {
	    ans$raw.scale <- scale
	    ans$raw.coefficients <- center
	    quantiel <- qnorm(0.9875)
	    weights <- ifelse(abs(resid/scale) <= quantiel, 1, 0)
	    reweighting <- cov.wt(y, wt = weights)
	    ans$coefficients <- reweighting$center
	    ans$scale <- sqrt(sum(weights)/(sum(weights) - 1) * reweighting$cov)
	    resid <- y - ans$coefficients
	    ans$crit <- sum(sort((y - center)^2, quan)[1:quan])
	    if (sum(weights) == n) {
		cdelta.rew <- 1
		correct.rew <- 1
	    }
	    else {
		qdelta.rew <- qchisq(sum(weights)/n, 1)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, 1/2 + 1)/(sum(weights)/n)
		cdelta.rew <- sqrt(1/cdeltainvers.rew)
		correct.rew <- LTScorrfactor.rew.s(1, intercept = intercept, n, alpha)
	    }
	    ans$scale <- ans$scale * cdelta.rew * correct.rew
	    weights <- ifelse(abs(resid/ans$scale) <= quantiel, 1, 0)
	}
	ans$resid <- resid/ans$scale
	ans$rsquared <- 0
	ans$residuals <- rep(NA, length(na.y))
	ans$residuals[ok] <- resid
	ans$lts.wt <- rep(NA, length(na.y))
	ans$lts.wt[ok] <- weights
	ans$intercept <- intercept
	ans$method <- paste("Univariate location and scale estimation.")
	if (abs(scale) < 1e-07)
	    ans$method <- paste(ans$method, "\nMore than half of the data are equal!")
	names(ans$coefficients) <- names(ans$raw.coefficients) <- yn
	names(ans$scale) <- names(ans$raw.scale) <- yn
	names(ans$rsquared) <- yn
	names(ans$crit) <- yn
	names(ans$residuals) <- rownames
	names(ans$lts.wt) <- rownames
	ans$X <- x
	ans$Y <- y	# VT:: 01.09.2004 - add y to the result object
	if (length(rownames))
	    dimnames(ans$X)[[1]] <- rownames[ok]
	else {
	    xx <- seq(1, length(na.x))
	    dimnames(ans$X) <- list(NULL, NULL)
	    dimnames(ans$X)[[1]] <- xx[ok]
	}
	class(ans) <- "lts"
	attr(ans, "call") <- sys.call()

	##cat("++++++ A - all x == 1...Ready and Return.\n")
	return(ans)
    }

    ans <- list()
    if (intercept) {
	dx <- dx + c(0, 1)
	xn <- c(xn, "Intercept")
	x <- array(c(x, rep(1, n)), dx, dimnames = list(dn[[1]], xn))
    }
    p <- dx[2]
    if (n <= 2 * p)
	stop("Need more than twice as many observations as variables.")
    if (length(alpha)) {
	if (alpha > 1)
	    stop("alpha is greater than 1")
	if (alpha == 1) { ## alpha == 1 -----------------------

	    ##cat("++++++ B - alpha == 1...\n")
	    z <- lsfit(x, y, intercept = FALSE)

	    ## VT:: 26.12.2004
	    ## Reorder the coeficients,so that the intercept moves to the beginning of the array
	    ## Skip this if p == 1 (i.e. p=1 and intercept=FALSE).
	    ## Do the same for the names and for ans$coef - see below

	    if(p > 1 && intercept) {
		ans$raw.coefficients[2:p] <- z$coef[1:(p - 1)]
		ans$raw.coefficients[1] <- z$coef[p]
		names(ans$raw.coefficients)[2:p] <- xn[1:(p - 1)]
		names(ans$raw.coefficients)[1] <- xn[p]
	    } else {
		ans$raw.coefficients <- z$coef
		names(ans$raw.coefficients) <- xn
	    }

	    ans$alpha <- alpha
	    ans$quan <- quan <- n   # VT:: 01.09.2004 - bug in alpha=1
	    ## (ans$quan was not set)

	    s0 <- sqrt((1/(n - p)) * sum(z$residuals^2))

	    ##cat("++++++ B - alpha == 1... - s0=",s0,"\n")
	    if(abs(s0) < 1e-07) {
		fitted <- x %*% z$coef
		weights <- ifelse(abs(z$residuals) <= 1e-07, 1, 0)
		ans$scale <- ans$raw.scale <- 0
		ans$coefficients <- ans$raw.coefficients
	    }
	    else {
		ans$raw.scale <- s0
		ans$raw.resid <- ans$residuals/ans$raw.scale
		weights <- ifelse(abs(z$residuals/s0) <= qnorm(0.9875), 1, 0)

		## vt:: weights has to be a vector instead of a matrix -
		##	to avoid "Error in x * wtmult : non-conformable arrays"
		##
		weights <- as.vector(weights)
		sum.w <- sum(weights)
		z <- lsfit(x, y, wt = weights, intercept = FALSE)

		## VT:: 26.12.2004
		ans$coefficients <-
		    if(p > 1 && intercept) z$coef[c(p, 1:(p - 1))] else z$coef

		fitted <- x %*% z$coef
		ans$scale <- sqrt(sum(weights * z$residuals^2)/(sum.w - 1))
		if (sum.w == n) {
		    cdelta.rew <- 1
		}
		else {
		    qn.w <- qnorm((sum.w + n)/(2 * n))
		    cdelta.rew <- 1/sqrt(1 - (2 * n)/(sum.w/qn.w) * dnorm(qn.w))
		}
		ans$scale <- ans$scale * cdelta.rew
		weights <- ifelse(abs(z$residuals/ans$scale) <= qnorm(0.9875), 1, 0)
		ans$resid <- z$residuals/ans$scale
	    }

	    ## VT:: 26.12.2004
	    names(ans$coefficients) <-
		if(p > 1 && intercept) xn[c(p, 1:(p - 1))] else xn

	    ans$crit <- sum(z$residuals^2)
	    if (intercept) {
		s1 <- sum(z$residuals^2)
		center <- mean(y)
		sh <- sum((y - center)^2)
	    }
	    else {
		s1 <- sum(z$residuals^2)
		sh <- sum(y^2)
	    }
	    ans$rsquared <- max(0, min(1, 1 - (s1/sh)))
	    ans$residuals <- rep(NA, length(na.y))
	    ans$residuals[ok] <- z$residuals
	    ans$lts.wt <- matrix(NA, length(na.y))
	    ans$lts.wt[ok] <- weights
	    ans$intercept <- intercept
	    ans$method <- paste("Least Squares Regression.")
	    if (abs(s0) < 1e-07)
		ans$method <- paste(ans$method, "\nAn exact fit was found!")
	    if (mcd) {
		## VT:: changed name of the function from 'cov.mcd.default' to 'covMcd'
		mcd <- covMcd(X, alpha = 1)
		if(-(determinant(mcd$cov, log = TRUE)$modulus[1])/p > 50) {
		    ans$RD[1] <- "singularity"
		} else {
		    ans$RD <- rep(NA, length(na.y))
		    ans$RD[ok] <- sqrt(mahalanobis(X, mcd$center, mcd$cov))
		    names(ans$RD) <- rownames
		}
	    }
	    names(ans$residuals) <- rownames
	    names(ans$lts.wt) <- rownames
	    names(ans$scale) <- names(ans$raw.scale) <- yn
	    names(ans$rsquared) <- yn
	    names(ans$crit) <- yn
	    ans$X <- x
	    ans$Y <- y	# VT:: 01.09.2004 - add y to the result object
	    if (length(rownames))
		dimnames(ans$X)[[1]] <- rownames[ok]
	    else {
		xx <- seq(1, length(na.x))
		dimnames(ans$X) <- list(NULL, NULL)
		dimnames(ans$X)[[1]] <- xx[ok]
	    }
	    ans$fitted.values <- rep(NA, length(na.y))
	    ans$fitted.values[ok] <- fitted
	    names(ans$fitted.values) <- rownames
	    if (qr.out)
		ans$qr <- z$qr
	    class(ans) <- "lts"
	    attr(ans, "call") <- sys.call()

	    ##cat("+++++ B - alpha == 1...Ready and return\n")
	    return(ans)
	}
    }


    coefs <- rep(NA, p)
    names(coefs) <- xn
    if(qr.out)
	qrx <- qr(x)
    else
	qrx <- qr(x)[c("rank", "pivot")]

    rk <- qrx$rank
    if (rk < p) {
	stop("x is singular")
    }
    else
	piv <- 1:p

    if (!length(alpha)) {
	alpha <- 1/2
	quan <- quan.f(alpha, n, rk)
    } else {
	if (alpha < 1/2)
	    stop("alpha is out of range!")
	quan <- quan.f(alpha, n, rk)
    }

    z <- .fastlts(x, y, quan, nsamp, intercept, adjust, seed)

    ## vt:: lm.fit.qr == lm.fit(...,method=qr,...)
    ##	cf <- lm.fit.qr(x[z$inbest, , drop = FALSE], y[z$inbest])$coef
    cf <- lm.fit(x[z$inbest, , drop = FALSE], y[z$inbest])$coef
    ans$best <- sort(as.vector(z$inbest))
    fitted <- x %*% cf
    resid <- y - fitted
    coefs[piv] <- cf

    ## VT:: 26.12.2004
    if(p > 1 && intercept) {
	ii <- c(p, 1:(p - 1))
	ans$raw.coefficients <- coefs[ii]
	names(ans$raw.coefficients) <- names(coefs)[ii] ## < unneeded ?
    } else {
	ans$raw.coefficients <- coefs
	names(ans$raw.coefficients) <- names(coefs) ## < unneeded ?
    }

    ans$alpha <- alpha
    ans$quan <- quan
    correct <- LTScorrfactor.s(p, intercept = intercept, n, alpha)
    s0 <- sqrt((1/quan) * sum(sort(resid^2, quan)[1:quan]))
    sh0 <- s0
    qn.q <- qnorm((quan + n)/ (2 * n))
    s0 <- s0 / sqrt(1 - (2 * n)/(quan / qn.q) * dnorm(qn.q)) * correct

    if (abs(s0) < 1e-07) {
	weights <- ifelse(abs(resid) <= 1e-07, 1, 0)
	ans$scale <- ans$raw.scale <- 0
	ans$coefficients <- ans$raw.coefficients
    }
    else {
	ans$raw.scale <- s0
	ans$raw.resid <- resid/ans$raw.scale
	quantiel <- qnorm(0.9875)
	weights <- ifelse(abs(resid/s0) <= quantiel, 1, 0)

	## vt:: weights has to be a vector instead of a matrix -
	##	to avoid "Error in x * wtmult : non-conformable arrays"
	##
	weights <- as.vector(weights)
	sum.w <- sum(weights)

	z1 <- lsfit(x, y, wt = weights, intercept = FALSE)

	## VT:: 26.12.2004
	if(p > 1) {
	    ans$coefficients[2:p] <- z1$coef[1:(p - 1)]
	    ans$coefficients[1] <- z1$coef[p]
	} else {
	    ans$coefficients <- z1$coef
	}

	fitted <- x %*% z1$coef
	resid <- z1$residuals
	ans$scale <- sqrt(sum(weights * resid^2)/(sum.w - 1))
	if (sum.w == n) {
	    cdelta.rew <- 1
	    correct.rew <- 1
	}
	else {
	    qn.w <- qnorm((sum.w + n)/(2 * n))
	    cdelta.rew <- 1 / sqrt(1 - (2 * n)/(sum.w / qn.w) * dnorm(qn.w))
	    correct.rew <- LTScorrfactor.rew.s(p, intercept = intercept, n, alpha)
	}
	ans$scale <- ans$scale * cdelta.rew * correct.rew
	ans$resid <- resid/ans$scale
	quantiel <- qnorm(0.9875)
	weights <- ifelse(abs(resid/ans$scale) <= quantiel, 1, 0)
    }
    names(ans$coefficients) <- names(ans$raw.coefficients)
    ans$lts.wt <- matrix(NA, length(na.y))
    ans$lts.wt[ok] <- weights
    ans$crit <- z$objfct
    if (intercept) {
	sh <- .fastmcd(as.matrix(y), as.integer(quan), nsamp = 0, seed)
	y <- as.vector(y)
	sh <- as.double(sh$adjustcov)
	ans$rsquared <- 1 - (sh0/sh)^2
    }
    else {
	s1 <- sum(sort(resid^2, quan)[1:quan])
	sh <- sum(sort(y^2, quan)[1:quan])
	ans$rsquared <- 1 - (s1/sh)
    }

    ## VT:: 03.11.04 - consider the case when sh=0 (i.e. rsquared=NaN or -Inf)
    if(is.nan(ans$rsquared) || is.infinite(ans$rsquared)) {
	ans$rsquared <- 0
    } else if (ans$rsquared > 1) {
	ans$rsquared <- 1
    } else if (ans$rsquared < 0) {
	ans$rsquared <- 0
    }

    attributes(resid) <- attributes(fitted) <- attributes(y)
    ans$residuals <- rep(NA, length(na.y))
    ans$residuals[ok] <- resid
    ans$intercept <- intercept
    ans$method <- paste("Least Trimmed Squares Robust Regression.")
    if(abs(s0) < 1e-07)
	ans$method <- paste(ans$method, "\nAn exact fit was found!")
    if (mcd) {
	mcd <- covMcd(X, alpha = alpha)
	if(-(determinant(mcd$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
	    ans$RD[1] <- "singularity"
	}
	else {
	    ans$RD <- rep(NA, length(na.y))
	    ans$RD[ok] <- sqrt(mahalanobis(X, mcd$center, mcd$cov))
	    names(ans$RD) <- rownames
	}
    }
    names(ans$residuals) <- rownames
    names(ans$lts.wt) <- rownames
    names(ans$scale) <- names(ans$raw.scale) <- yn
    names(ans$rsquared) <- yn
    names(ans$crit) <- yn
    ans$X <- x
    ans$Y <- y		# VT:: 01.09.2004 - add y to the result object
    if (length(rownames))
	dimnames(ans$X)[[1]] <- rownames[ok]
    else {
	xx <- seq(1, length(na.x))
	dimnames(ans$X) <- list(NULL, NULL)
	dimnames(ans$X)[[1]] <- xx[ok]
    }
    ans$fitted.values <- rep(NA, length(na.y))
    ans$fitted.values[ok] <- fitted
    names(ans$fitted.values) <- rownames
    if (qr.out)
	ans$qr <- qrx
    class(ans) <- "lts"
    attr(ans, "call") <- sys.call()
    return(ans)
}

##predict.lts <- function (object, newdata, na.action = na.pass, ...)
##{
##    if (missing(newdata)) return(fitted(object))
##    ## work hard to predict NA for rows with missing data
##    Terms <- delete.response(terms(object))
##    m <- model.frame(Terms, newdata, na.action = na.action,
##		       xlev = object$xlevels)
##    if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
##    X <- model.matrix(Terms, m, contrasts = object$contrasts)
##    drop(X %*% object$coefficients)
##}

summary.lts <- function (object, correlation = FALSE, ...)
{
    z <- object
    r <- z$residuals
    f <- z$fitted
    w <- as.vector(z$lts.wt)
    n <- sum(w)

##    Qr <- qr(object$X)
    Qr <- qr(t(t(object$X) %*% diag(as.vector(w))))
    p <- Qr$rank
    p1 <- 1:p

    rdf <- n - p

    mss <-  if(z$intercept) {
		m <- sum(w * f /sum(w))
		sum(w * (f - m)^2)
	    } else
		sum(w * f^2)
    rss <- sum(w * r^2)

    r <- sqrt(w) * r
    resvar <- rss/rdf

    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])

    ## Reorder R, so that the intercept (if any) moves
    ## to the beginning. Skip this if p == 1 or intercept=FALSE.
	RR <- R
	RR[2:p, 2:p] <- R[1:(p - 1), 1:(p-1)]
	rr <- R[p,]
	rr[2:p] <- R[p, 1:(p - 1)]
	rr[1] <- R[p,p]
	RR[1,] <- rr
	RR[,1] <- rr
	R <- RR
    se <- sqrt(diag(R) * resvar)

    est <- z$coefficients
    tval <- est/se

    ans <- z[c("call", "terms")]
    attr(ans, "call") <- attr(z,"call")
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2*pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients),
				 c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))

    df.int <- if(z$intercept) 1 else 0
    if(p - df.int > 0) {
	ans$r.squared <- mss/(mss + rss)
	ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
	ans$fstatistic <- c(value = (mss/(p - df.int))/resvar, numdf = p - df.int, dendf = rdf)
    } else
	ans$r.squared <- ans$adj.r.squared <- 0

    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]

    if (correlation) {
	ans$correlation <- (R * resvar)/outer(se, se)
	dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    }

    class(ans) <- "summary.lts"
    ans
}

print.summary.lts <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    xcall <- x$call
    if(length(xcall) == 0)
	xcall <- attr(x,"call")
    cat(paste(deparse(xcall), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    resid <- x$residuals
    df <- x$df
    rdf <- df[2]

    if(rdf > 5) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <-	if(length(dim(resid)) == 2)
		    structure(apply(t(resid), 1, quantile),
			      dimnames = list(nam, dimnames(resid)[[2]]))
		else
		    structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else if(rdf > 0) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }

    if (nsingular <- df[3] - df[1])
	cat("\nCoefficients: (", nsingular,
	    " not defined because of singularities)\n", sep = "")
    else
	cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE, na.print = "NA", ...)

    cat("\nResidual standard error:",
    format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")

    if(!is.null(x$fstatistic)) {
	cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
	cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits = digits),
	    "\nF-statistic:", formatC(x$fstatistic[1], digits = digits),
	    "on", x$fstatistic[2], "and",
	    x$fstatistic[3], "DF,  p-value:",
	    format.pval(pf(x$fstatistic[1], x$fstatistic[2],
			   x$fstatistic[3], lower.tail = FALSE), digits = digits),
	"\n")
    }

    correl <- x$correlation
    if(!is.null(correl)) {
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    correl <- format(round(correl, 2), nsmall = 2, digits = digits)
	    correl[!lower.tri(correl)] <- ""
	    print(correl[-1, -p, drop = FALSE], quote = FALSE)
	}
    }
    cat("\n")
    invisible(x)
}


### --- Namespace hidden (but parsed once and for all) : -------------

LTScorrfactor.s <- function(p, intercept = intercept, n, alpha)
{
    stopifnot(0.5 <= alpha, alpha <= 1)
    if (intercept)
	p <- p - 1
    stopifnot(p == as.integer(p), p >= 0)
    if (p == 0) {
	fp.500.n <- 1 - exp( 0.262024211897096) / n^ 0.604756680630497
	fp.875.n <- 1 - exp(-0.351584646688712) / n^ 1.01646567502486
	if ((0.5 <= alpha) && (alpha <= 0.875)) {
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	    fp.alpha.n <- sqrt(fp.alpha.n)
	}
	if ((0.875 < alpha) && (alpha < 1)) {
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
	    fp.alpha.n <- sqrt(fp.alpha.n)
	}
    }
    else { ## p >= 1
        if (p == 1) {
            if (intercept) {
                fp.500.n <- 1 - exp( 0.630869217886906 ) / n^ 0.650789250442946
                fp.875.n <- 1 - exp( 0.565065391014791 ) / n^ 1.03044199012509
            }
            else {
                fp.500.n <- 1 - exp(-0.0181777452315321) / n^ 0.697629772271099
                fp.875.n <- 1 - exp(-0.310122738776431 ) / n^ 1.06241615923172
            }
        } else { ## --- p > 1 ---
            if (intercept) {
		##                           "alfaq"            "betaq"    "qwaarden"
                coefgqpkwad875 <- matrix(c(-0.458580153984614, 1.12236071104403, 3,
                                           -0.267178168108996, 1.1022478781154,  5), ncol = 2)
                coefeqpkwad500 <- matrix(c(-0.746945886714663, 0.56264937192689,  3,
                                           -0.535478048924724, 0.543323462033445, 5), ncol = 2)
            }
            else {
		##                           "alfaq"            "betaq"    "qwaarden"
                coefgqpkwad875 <- matrix(c(-0.251778730491252, 0.883966931611758, 3,
                                           -0.146660023184295, 0.86292940340761,  5), ncol = 2)
                coefeqpkwad500 <- matrix(c(-0.487338281979106, 0.405511279418594, 3,
                                           -0.340762058011,    0.37972360544988,  5), ncol = 2)
            }

            y.500 <- log(- coefeqpkwad500[1, ] / p^ coefeqpkwad500[2, ])
            y.875 <- log(- coefgqpkwad875[1, ] / p^ coefgqpkwad875[2, ])

            A.500 <- cbind(1, - log(coefeqpkwad500[3, ] * p^2))
            coeffic.500 <- solve(A.500, y.500)
            A.875 <- cbind(1, - log(coefgqpkwad875[3, ] * p^2))

            coeffic.875 <- solve(A.875, y.875)
            fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
            fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
        }

        if(alpha <= 0.875)
            fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
        else ##	 0.875 < alpha <= 1
            fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
    }## else (p >= 1)

    return(1/fp.alpha.n)
} ##  LTScorrfactor.s

LTScorrfactor.rew.s <- function(p, intercept = intercept, n, alpha)
{
    stopifnot(0.5 <= alpha, alpha <= 1)
    if (intercept)
	p <- p - 1
    stopifnot(p == as.integer(p), p >= 0)

    if (p == 0) {
	fp.500.n <- 1 - exp( 1.11098143415027) / n^ 1.5182890270453
	fp.875.n <- 1 - exp(-0.66046776772861) / n^ 0.88939595831888

	if(alpha <= 0.875)
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	else ##	 0.875 < alpha <= 1
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
        ## MM: sqrt() {below} is ``different logic'' than below.. (??)
        fp.alpha.n <- sqrt(fp.alpha.n)
    }
    else {
	if (p == 1) {
	    if (intercept) {
		fp.500.n <- 1 - exp(1.58609654199605 ) / n^ 1.46340162526468
		fp.875.n <- 1 - exp(0.391653958727332) / n^ 1.03167487483316
	    }
	    else {
		fp.500.n <- 1 - exp( 0.6329852387657)	/ n^ 1.40361879788014
		fp.875.n <- 1 - exp(-0.642240988645469) / n^ 0.926325452943084
	    }
	}
	else { ##  --- p > 1 ---
	    if (intercept) {
		##                           "alfaq"            "betaq"    "qwaarden"
		coefqpkwad875 <- matrix(c(-0.474174840843602, 1.39681715704956, 3,
                                          -0.276640353112907, 1.42543242287677, 5), ncol = 2)
		coefqpkwad500 <- matrix(c(-0.773365715932083, 2.02013996406346, 3,
                                          -0.337571678986723, 2.02037467454833, 5), ncol = 2)
	    }
	    else {
		##                           "alfaq"            "betaq"    "qwaarden"
		coefqpkwad875 <- matrix(c(-0.267522855927958, 1.17559984533974, 3,
                                          -0.161200683014406, 1.21675019853961, 5), ncol = 2)
		coefqpkwad500 <- matrix(c(-0.417574780492848, 1.83958876341367, 3,
                                          -0.175753709374146, 1.8313809497999, 5), ncol = 2)
	    }
	    y.500 <- log( - coefqpkwad500[1, ] / p^ coefqpkwad500[2, ])
	    y.875 <- log( - coefqpkwad875[1, ] / p^ coefqpkwad875[2, ])
	    A.500 <- cbind(1, - log(coefqpkwad500[3, ] * p^2))
	    coeffic.500 <- solve(A.500, y.500)
	    A.875 <- cbind(1, - log(coefqpkwad875[3, ] * p^2))
	    coeffic.875 <- solve(A.875, y.875)
	    fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
	    fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
        }

	if(alpha <= 0.875)
	    fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
	else ##	 0.875 < alpha <= 1
	    fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    }## else (p >= 1)

    return(1/fp.alpha.n)
} ## LTScorrfactor.rew.s

.fastlts <- function(x, y, quan, nsamp, intercept, adjust, seed) {
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    y <- as.matrix(y)
    x1 <- matrix(0, ncol = p + 1, nrow = n)
    x1 <- cbind(x, y)
    x1 <- as.matrix(x1)

    objfct <- 0
    interc <- ifelse(intercept, 1, 0)
    intadjust <- ifelse(adjust, 1, 0)

    storage.mode(x1) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(quan) <- "integer"
    storage.mode(nsamp) <- "integer"
    storage.mode(seed) <- "integer"
    storage.mode(objfct) <- "double"
    storage.mode(interc) <- "integer"
    storage.mode(intadjust) <- "integer"

    inbest <- matrix(10000, nrow = quan, ncol = 1)

    storage.mode(inbest) <- "integer"

    datt <- matrix(0, ncol = p + 1, nrow = n)
    storage.mode(datt) <- "double"
    nvad <- p + 1
    storage.mode(nvad) <- "integer"

##   Allocate temporary storage for the fortran implementation

    weights <- matrix(0, nrow = n, ncol = 1)
    temp <- matrix(0, nrow = n, ncol = 1)
    index1 <- matrix(0, nrow = n, ncol = 1)
    index2 <- matrix(0, nrow = n, ncol = 1)
    aw2 <- matrix(0, nrow = n, ncol = 1)
    aw <- matrix(0, nrow = n, ncol = 1)
    residu <- matrix(0, nrow = n, ncol = 1)
    yy <- matrix(0, nrow = n, ncol = 1)
    nmahad <- matrix(0, nrow = n, ncol = 1)
    ndist <- matrix(0, nrow = n, ncol = 1)
    am <- matrix(0, nrow = n, ncol = 1)
    am2 <- matrix(0, nrow = n, ncol = 1)
    slutn <- matrix(0, nrow = n, ncol = 1)

    storage.mode(weights) <- "double"
    storage.mode(temp) <- "integer"
    storage.mode(index1) <- "integer"
    storage.mode(index2) <- "integer"
    storage.mode(aw2) <- "double"
    storage.mode(aw) <- "double"
    storage.mode(residu) <- "double"
    storage.mode(yy) <- "double"
    storage.mode(nmahad) <- "double"
    storage.mode(ndist) <- "double"
    storage.mode(am) <- "double"
    storage.mode(am2) <- "double"
    storage.mode(slutn) <- "double"

    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

##   integer jmiss(nvad)		 --> p+1
    jmiss <- matrix(0, nrow = p+1, ncol = 1)
    storage.mode(jmiss) <- "integer"
##   double precision xmed(nvad)	 --> p+1
    xmed <- matrix(0, nrow = p+1, ncol = 1)
    storage.mode(xmed) <- "double"
##   double precision xmad(nvad)	 p+1
    xmad <- matrix(0, nrow = p+1, ncol = 1)
    storage.mode(xmad) <- "double"
##   double precision a(nvad)		 p+1
    a <- matrix(0, nrow = p+1, ncol = 1)
    storage.mode(a) <- "double"

##   double precision da(nvad)		 p+1
    da <- matrix(0, nrow = p+1, ncol = 1)
    storage.mode(da) <- "double"
##   double precision h(nvar,nvad)	     p*(p+1)
    h <- matrix(0, nrow = p*(p+1), ncol = 1)
    storage.mode(h) <- "double"
##   double precision hvec(nvar*nvad)	     p*(p+1)
    hvec <- matrix(0, nrow = p*(p+1), ncol = 1)
    storage.mode(hvec) <- "double"
##   double precision c(nvar,nvad)	     p*(p+1)
    c <- matrix(0, nrow = p*(p+1), ncol = 1)
    storage.mode(c) <- "double"
##   double precision cstock(10,nvar*nvar)   10*p*p
    cstock <- matrix(0, nrow = 10*p*p, ncol = 1)
    storage.mode(cstock) <- "double"
##   double precision mstock(10,nvar)	     10*p
    mstock <- matrix(0, nrow = 10*p, ncol = 1)
    storage.mode(mstock) <- "double"
##   double precision c1stock(km10,nvar*nvar)	 km10*p*p
    c1stock <- matrix(0, nrow = km10*p*p, ncol = 1)
    storage.mode(c1stock) <- "double"
##    double precision m1stock(km10,nvar)	  km10*p
    m1stock <- matrix(0, nrow = km10*p, ncol = 1)
    storage.mode(m1stock) <- "double"
##    double precision dath(nmaxi,nvad)		  nmaxi*(p+1)
    dath <- matrix(0, nrow = nmaxi*(p+1), ncol = 1)
    storage.mode(dath) <- "double"
##   double precision sd(nvar)		     p
    sd <- matrix(0, nrow = p, ncol = 1)
    storage.mode(sd) <- "double"
##   double precision means(nvar)	     p
    means <- matrix(0, nrow = p, ncol = 1)
    storage.mode(means) <- "double"
##   double precision bmeans(nvar)	     p
    bmeans <- matrix(0, nrow = p, ncol = 1)
    storage.mode(bmeans) <- "double"

    zlts <- .Fortran("rfltsreg",
	    x1 = x1,
	    n,
	    p,
	    quan,
	    nsamp,
	    inbest = inbest,
	    objfct = objfct,
	    interc,
	    intadjust,
	    nvad,
	    datt,
	    seed,
	    weights,
	    temp,
	    index1,
	    index2,
	    aw2,
	    aw,
	    residu,
	    yy,
	    nmahad,
	    ndist,
	    am,
	    am2,
	    slutn,
	    jmiss,xmed,xmad,a,da,h,hvec,c,cstock,mstock,c1stock,
	    m1stock,dath,sd,means,bmeans,
	    PACKAGE = "robustbase")
    zlts
}

