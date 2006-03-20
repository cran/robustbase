.onLoad <- function(libname, pkgname) {
    require("methods")
}

if(paste(R.version$major, R.version$minor, sep=".") < 2.3) {
    ## These are substitutes such that newer code still runs in older R
    tcrossprod <- function(x, y = NULL) x %*% t(if(is.null(y)) x else y)
}

