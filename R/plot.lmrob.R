
plot.lmrob <-
function (x, which = 1:5,
          caption = c("Standardized residuals vs. Robust Distances",
          "Normal Q-Q vs. Residuals", "Response vs. Fitted Values",
          "Residuals vs. Fitted Values" ,
          "Sqrt of abs(Residuals) vs. Fitted Values"),
          panel = points, sub.caption = deparse(x$call), main = "",
          compute.MD = TRUE, # maybe  (n < 1000 && p < 20)
          ask = prod(par("mfcol")) < length(which) && dev.interactive(),
          ..., p = 0.025)
{
    if (!inherits(x, "lmrob"))
        stop("Use only with 'lmrob' objects")
    show <- rep(FALSE, 5)
    if (!is.numeric(which) || any(which < 1) || any(which > 5))
        stop("'which' must be in 1:5")
    show[which] <- TRUE
    r <- residuals(x)
    n <- length(r)
    sr <- r/x$scale
    yh <- fitted(x)
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }
    if (show[1]) {
	if(is.null(x[['MD']]) && compute.MD) {
	    message("recomputing robust Mahalanobis distances")
	    x$MD <- ## need to recompute
		robMD(x = if(!is.null(x[['x']])) x$x else
		      if(!is.null(x[['model']])) model.matrix(x, x$model)
		      else stop("need 'model' or 'x' component for robust Mahalanobis distances"),
		      intercept = attr(x$terms,"intercept"))
	    ## try to "cache" them with the object
	    if(identical(parent.frame(), .GlobalEnv) &&
	       exists((cnx <- as.character(nx <- match.call()[["x"]])),
		      .GlobalEnv)) {
		assign(cnx, x, envir = .GlobalEnv)
		message("saving the robust distances 'MD' as part of ",
			sQuote(cnx))
	    }
	}
        if(!is.null(x[['MD']])) {
            if (p < 0 || p > 1)
                stop ("Tolerance range must be between 0% to 100%")
            else chi <- sqrt( qchisq(p = 1-p, df = x$rank) )
            plot(x$MD,sr,
                 xlab = "Robust Distances",
                 ylab = "Robust Standardized residuals",
                 main = main, type = "p", ...)
            mtext(caption[1], 3, 0.25)
            if (one.fig)
                title(sub = sub.caption, ...)
            abline(h = c(2.5,-2.5), lty = 3)
            abline(v = chi, lty = 3)
        }
    }
    if (show[2]) {
        qqnorm(r, ylab = "Residuals", main = main,...)
        qqline(r)
        mtext(caption[2], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
    }
    if (show[3]) {
        y <- if(!is.null(x[['model']])) model.response(x$model) else yh + r
        m1 <- min(yh,y)
        m2 <- max(yh,y)
        plot(yh, y, xlab = "Fitted Values", ylab = "Response",
             xlim = c(m1,m2), ylim = c(m1,m2), main = main, ...)
        mtext(caption[3], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
        abline(a = 0,b = 1)
    }
    if (show[4]) {
        plot(yh, r, xlab = "Fitted Values", ylab = "Residuals", main = main, ...)
        mtext(caption[4], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
        abline(h = c(2.5*x$scale,0,-2.5*x$scale), lty = 3)
    }
    if (show[5]) {
        sqrtabsr <- sqrt(abs(r))
        plot(yh, sqrtabsr, xlab = "Fitted Values", ylab = "Sqrt of abs(Residuals)",
             main = main, ...)
        mtext(caption[5], 3, 0.25)
        if (one.fig)
            title(sub = sub.caption, ...)
    }
    invisible()
}

