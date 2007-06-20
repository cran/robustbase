.onLoad <- function(libname, pkgname) {
    require("methods")
}

if(getRversion() < "2.3") {
    ## These are substitutes such that newer code still runs in older R
    tcrossprod <- function(x, y = NULL) x %*% t(if(is.null(y)) x else y)
}

if(getRversion() < "2.6") {
    nzchar <- function(x) nchar(x) > 0
}
