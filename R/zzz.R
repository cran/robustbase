.onLoad <- function(libname, pkgname) {
    require("methods")
}

if(getRversion() < "2.9") {
    ## The "real" grepl() function has more arguments and is more powerful
    ## This is sufficient for a workaround for old R versions:
    grepl <- function (pattern, x, ...) {
	seq_along(x) %in% grep(pattern, x, ...)
    }
}
