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
plot.lts <- function(x,
		   which = c("all", "rqq","rindex", "rfit", "rdiag"),
		   classic = FALSE,
		   ask = (which == "all" && dev.interactive()),
		   id.n, ...) {
    if (!inherits(x, "lts"))
	stop("Use only with 'lts' objects")

    ltsPlot(x, which, classic, ask, id.n, ...)
}

ltsPlot <- function(x,
		   which = c("all", "rqq","rindex", "rfit", "rdiag"),
		   classic = FALSE,
		   ask = FALSE,
		   id.n, ...)
{
    ##@bdescr
    ##	Make plots for model checking and outlier detection based on
    ##	    the LTS regression estimates:
    ##	rqq	 -  normal quantile plot of the LTS and LS residuals
    ##	rindex	 -  standardized LTS/LS Residuals versus index
    ##	rfit	 -  standardized LTS/LS Residuals versus fitted values
    ##	rdiag	 -  regression diagnostic plot
    ##
    ##@edescr
    ##
    ##@in  x		   : [object] An lts object
    ##@in  which	  : [character] A plot option, one of:
    ##				  rqq:
    ##				  rdiag:
    ##				  rfit:
    ##				  rindex:
    ##				default is "rqq"
    ##@in  classic	     : [logical] If true the classical plot will be displayed too
    ##					 default is classic=FALSE
    ##@in  id.n		      : [number] number of observations to be identified with a label.

    label <- function(x, y, ord, lab, id.n, ...)
    {
	if(id.n) {
	    n <- length(y)
	    which <- order(ord)[(n - id.n + 1):n]
	    if(missing(lab))
		lab <- which
	    else
		lab <- lab[which]
	    ## how to adjust the labels?
	    ## a) adj=0.1
	    ## b) x=x+xrange
	    ## c) pos=4 (to the left of the observation)
	    ## d) additionaly to pos specify offset=0.2 (fraction of a character)
	    xrange <- par("usr")
	    xrange <- (xrange[2] - xrange[1])/50
	    text(x[which], y[which], pos = 4, offset = 0.2, lab, ...)
	}
    }

    ## The R function 'qqline' (package::stats) adds a line to a
    ## normal quantile-quantile plot which passes through the
    ## first and third quartiles. In S this function returns the
    ## slope and intercept of the line, but not in R.
    ## Here we need the slope and intercept in order to sort the
    ## residuals according to their distance from the line.

    myqqline <- function(y, datax = FALSE, ...) {
	y <- quantile(y[!is.na(y)],c(0.25, 0.75))
	x <- qnorm(c(0.25, 0.75))
	if(datax) {
	    slope <- diff(x)/diff(y)
	    int <- x[1] - slope*y[1]
	} else {
	    slope <- diff(y)/diff(x)
	    int <- y[1]-slope*x[1]
	}
	abline(int, slope, ...)
	invisible(list(int = int, slope = slope))
    }

    myqqplot <- function(r, classic = FALSE, lab, id.n, ...) {
	##  Normal QQ-plot of residuals:
	##  Produces a Quantile-Quantile plot in which the vector r is plotted
	##  against the quantiles of a standard normal distribution.
	##

	mgp = c(2.5, 1, 0) # set the margin line (in 'mex' units) for the:
	## - axis title,
	## - axis labels and
	## - axis line.
	## The default is 'c(3, 1, 0)'.
	xlab <- "Quantiles of the standard normal distribution"
	ylab <- "Standardized LTS residual"
	if(classic)
	    ylab <- "Standardized LS residual"

	qq <- qqnorm(r, mgp = mgp, xlab = xlab, ylab = ylab, ...)
	ll <- myqqline(r, lty = 2, ...)
	ord <- abs(qq$y - ll$int - ll$slope * qq$x)
	label(qq$x, qq$y, ord, lab, id.n, ...)
	title(main = "Normal Q-Q plot")
    }

    indexplot <- function(r, scale, classic = FALSE, lab, id.n, ...) {
	##  Index plot of standardized residuals:
	##  Plot the vector r (LTS or LS residuals) against
	##  the observation indexes. Identify by a label the id.n
	##  observations with largest value of r.
	##  Use classic=FALSE/TRUE to choose the label of the vertical axes

	## VT:: 26.12.2004
	if(scale == 0)
	    stop("Index plot of standardized residuals is not avalable if scale = 0")

	mgp = c(2.5, 1, 0) # set the margin line (in 'mex' units) for the:
	## - axis title,
	## - axis labels and
	## - axis line.
	## The default is 'c(3, 1, 0)'.
	xlab <- "Index"
	ylab <- "Standardized LTS residual"
	if(classic)
	    ylab <- "Standardized LS residual"

	x <- 1:length(r)
	y <- r/scale
	ylim <- c(min(-3, min(y)), max(3, max(y)))

	plot(x, y, ylim = ylim, mgp = mgp, xlab = xlab, ylab = ylab, ...)
	label(x, y, ord = abs(y), lab, id.n, ...)
	abline(h = -2.5, ...)
	abline(h = 0, lty = 4, ...)
	abline(h = 2.5, ...)

	mtext("-2.5", side = 4, line = 1.2, at = -2.5, ...)
	mtext("2.5", side = 4, line = 1.2, at = 2.5, ...)

	title(main = "Residuals vs Index")
    }

    fitplot <- function(obj, classic = FALSE, lab, id.n, ...) {
	##  Standardized residuals vs Fitted values plot:
	##  Plot the vector r (LTS or LS residuals) against
	##  the corresponding fitted values. Identify by a
	##  label the id.n observations with largest value of r.
	##  Use classic=FALSE/TRUE to choose the label of the vertical axes

	## VT:: 26.12.2004
	if(obj$scale == 0)
	    stop("Standardized residuals vs Fitted values plot is not avalable if scale = 0")

	mgp = c(2.5, 1, 0) # set the margin line (in 'mex' units) for the:
	## - axis title,
	## - axis labels and
	## - axis line.
	## The default is 'c(3, 1, 0)'.

	##    x <- obj$X %*% as.matrix(obj$coef)
	x <- obj$fitted.values
	y <- obj$residuals/obj$scale
	ylim <- c(min(-3, min(y)), max(3, max(y)))
	yname <- names(obj$scale)
	xlab <- paste("Fitted :", yname)
	ylab <- "Standardized LTS residual"
	if(classic)
	    ylab <- "Standardized LS residual"

	plot(x, y, ylim = ylim, mgp = mgp, xlab = xlab, ylab = ylab, ...)
	label(x, y, ord = abs(y), lab, id.n, ...)
	abline(h = -2.5, ...)
	abline(h = 0, lty = 4, ...)
	abline(h = 2.5, ...)

	mtext("-2.5", side = 4, line = 1.2, at = -2.5, ...)
	mtext("2.5", side = 4, line = 1.2, at = 2.5, ...)

	title(main = "Residuals vs Fitted")

    }


    rdiag <- function(obj, classic = FALSE, lab, id.n, ...) {
	##  Regression diagnostic plot:
	##  Plot the vector of the standardized residuals against
	##  the robust distances of the predictor variables
	##  Identify by a label the id.n observations with largest value of r.
	##  Use classic=FALSE/TRUE to choose the label of the vertical axes

	p <- if(obj$intercept) length(obj$coef) - 1 else length(obj$coef)
	if(p <= 0)
	    warning("Diagnostic plot is not available for univar\niate location and scale estimation")

	## VT:: 26.12.2004
	if(obj$scale <= 0)
	    stop("Regression Diagnostic plot is not avalable if scale = 0")

	if(is.null(obj$RD))
	    stop("Regression Diagnostic plot is not avalable: option mcd=F was set in ltsReg().")
	if(obj$RD[1] == "singularity")
	    stop("The MCD covariance matrix was singular.")

	mgp = c(2.5, 1, 0) # set the margin line (in 'mex' units) for the:
	## - axis title,
	## - axis labels and
	## - axis line.
	## The default is 'c(3, 1, 0)'.

	xlab <- "Robust distance computed by MCD"
	ylab <- "Standardized LTS residual"
	if(classic) {
	    xlab <- "Mahalanobis distance"
	    ylab <- "Standardized LS residual"
	}

	## VT:: 18.01.20045
	## set id.n to the number of all outliers:
	##  regression outliers (weight==0)+ leverage points (RD > cutoff)
	if(missing(id.n)) {
	    id.n <- length(unique(c(which(obj$RD > sqrt(qchisq(0.975, p))), which(obj$lts.wt == 0))))
	}

	quant <- max(c(sqrt(qchisq(0.975, p)), 2.5))
	x <- obj$RD
	y <- obj$residuals/obj$scale
	xlim <- c(0, max(quant + 0.1, max(x)))
	ylim <- c(min(-3, min(y)), max(3, max(y)))

	plot(x, y, ylim = ylim, mgp = mgp, xlab = xlab, ylab = ylab, ...)
	ord = apply(abs(cbind(x/2.5, y/quant)), 1, max)
	label(x, y, ord = ord, lab, id.n, ...)
	abline(h = -2.5, ...)
	abline(h = 2.5, ...)
	abline(v = quant, ...)

	mtext("-2.5", side = 4, line = 1.2, at = -2.5, ...)
	mtext("2.5", side = 4, line = 1.2, at = 2.5,...)

	title(main = "Regression Diagnostic Plot")
    }

    ##	parameters and preconditions

    which <- match.arg(which)
    r <- residuals(x)
    n <- length(r)

    id.n.default <- TRUE # if id.n is missing, it will be set to a default for
    ## for each plot.
    if(!missing(id.n) && !is.null(id.n)) {
	id.n.default <- FALSE
	id.n <- as.integer(id.n)
	if(id.n < 0 || id.n > n)
	    stop("`id.n' must be in {1,..,",n,"}")
    }

    if(!classic)
	par(mfrow = c(1,1), pty = "m")
    else {
	par(mfrow = c(1,2), pty = "m")

	## calculate the LS regression (using LTS with alpha = 1)
	## if intercept, obj$X is augmented with a column of 1s - remove it

	if(x$intercept &&		# model with intercept
	   length(dim(x$X)) == 2 &&	# X is 2-dimensional
	   dim(x$X)[2] > 1 &&		# X has more than 1 column
	   all(x$X[,dim(x$X)[2]] == 1)) # the last column of X is all 1s
	    X <- x$X[,1:(dim(x$X)[2]-1)]
	else
	    X <- x$X
	obj.cl <- ltsReg(X, x$Y, intercept = x$intercept, alpha = 1)
    }

    if (ask) {
	op <- par(ask = TRUE)
	on.exit(par(op))
    }

    if(which == "all" || which == "rqq") {
	nx <- if(id.n.default) length(which(x$lts.wt == 0)) else id.n # set id.n to the number of regression outliers (weight==0)
	myqqplot(x$residuals, id.n = nx, ...) # normal QQ-plot of the LTS residuals
	if(classic)
	    myqqplot(obj.cl$residuals, classic = TRUE, id.n = nx, ...) # normal QQ-plot of the LS residuals
    }

    if(which == "all" || which == "rindex") {
	nx <- if(id.n.default) length(which(x$lts.wt == 0)) else id.n # set id.n to the number of regression outliers (weight==0)
	indexplot(x$residuals, x$scale, id.n = nx, ...) # index plot of the LTS residuals
	if(classic)
	    indexplot(obj.cl$residuals, obj.cl$scale, classic = TRUE, id.n = nx, ...) # index plot of the LS residuals
    }

    if(which == "all" || which == "rfit") {
	nx <- if(id.n.default) length(which(x$lts.wt == 0)) else id.n # set id.n to the number of regression outliers (weight==0)
	fitplot(x, id.n = nx, ...)
	if(classic)
	    fitplot(obj.cl, classic = TRUE, id.n = nx, ...)
    }

    if(which == "all" || which == "rdiag") {
	rdiag(x, id.n = id.n, ...)
	if(classic)
	    rdiag(obj.cl, classic = TRUE, id.n = id.n, ...)
    }
}
