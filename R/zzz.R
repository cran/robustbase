## Here and in NAMESPACE:
## if(getRversion() < "2.13") {
##     nobs <- function (object, ...) UseMethod("nobs")
##     ## also used for mlm fits  *and* lmrob :
##     nobs.lm <- function(object, ...)
## 	if(!is.null(w <- object$weights)) sum(w != 0) else NROW(object$residuals)
##     ## for glmrob :
##     nobs.glm <- function(object, ...) sum(!is.na(object$residuals))
## }

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_robustbase_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
