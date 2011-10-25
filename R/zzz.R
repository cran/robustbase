## Here and in NAMESPACE:
if(getRversion() < "2.13") {
    nobs <- function (object, ...) UseMethod("nobs")
    ## also used for mlm fits  *and* lmrob :
    nobs.lm <- function(object, ...)
	if(!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)
    ## for glmrob :
    nobs.glm <- function(object, ...) sum(!is.na(object$residuals))

    if(getRversion() < "2.12") {
	droplevels <- function(x, ...) UseMethod("droplevels")
	droplevels.factor <- function(x, ...) factor(x)
	droplevels.data.frame <- function(x, except = NULL, ...)
	{
	    ix <- vapply(x, is.factor, NA)
	    if (!is.null(except)) ix[except] <- FALSE
	    x[ix] <- lapply(x[ix], factor)
	    x
	}

    }
}
