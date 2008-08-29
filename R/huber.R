
##  A modified "safe" (and more general) Huber estimator:
huberM <-
    function(x, k = 1.5, weights = NULL,
	     tol = 1e-06,
	     mu = if(is.null(weights)) median(x) else wgt.himedian(x, weights),
	     s = if(is.null(weights)) mad(x, center=mu)
		 else wgt.himedian(abs(x - mu), weights),
	     warn0scale = getOption("verbose"))
{
    ## Author: Martin Maechler, Date: 6 Jan 2003

    ## implicit 'na.rm = TRUE':
    if(any(i <- is.na(x))) {
        x <- x[!i]
        if(!is.null(weights)) weights <- weights[!i]
    }
    n <- length(x)
    sum.w <-
        if(!is.null(weights)) {
            stopifnot(is.numeric(weights), weights >= 0, length(weights) == n)
            sum(weights)
        } else n
    it <- 0:0
    if(sum.w == 0) # e.g 'x' was all NA
	return(list(mu = NA, s = NA, it = it)) # instead of error

    if (s <= 0) {
        if(s < 0) stop("negative scale 's'")
        if(warn0scale && n > 1)
            warning("scale 's' is zero -- returning initial 'mu'")
    }
    else {
        wsum <- if(is.null(weights)) sum else function(u) sum(u * weights)
	repeat {
	    it <- it + 1:1
            y <- pmin(pmax(mu - k * s, x), mu + k * s)
	    mu1 <- wsum(y) / sum.w
	    if (abs(mu - mu1) < tol * s)
		break
	    mu <- mu1
	}
    }
    list(mu = mu, s = s, it = it)
}

## this is a compatible improvement of MASS' huber() :
## 1) returning median() if mad()=0
## 2)	"	NA when y has only NAs (or length 0)

if(FALSE)
huber <- function (y, k = 1.5, tol = 1e-06)
{
    y <- y[!is.na(y)]
    n <- length(y)
    if(n == 0) # e.g 'y' was all na
	return(list(mu = NA, s = NA))# instead of error
    mu <- median(y)
    s <- mad(y)
    if (s == 0) { # FIXME?  make this warning optional
	if(n > 1) warning("scale MAD is zero for this sample")
    }
    else repeat {
	yy <- pmin(pmax(mu - k * s, y), mu + k * s)
	mu1 <- sum(yy)/n
	if (abs(mu - mu1) < tol * s)
	    break
	mu <- mu1
    }
    list(mu = mu, s = s)
}
