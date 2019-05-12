## <---> sync with ~/R/Pkgs/CLA/inst/xtraR/platform-sessionInfo.R
##                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##' return 'x' unless it is NULL where you'd use 'orElse'
`%||%` <- function(x, orElse) if(!is.null(x)) x else orElse

##'  not %in%  :
`%nin%` <- function (x, table) is.na(match(x, table))

##' Derive more sessionInfo() like information, notably about BLAS, LAPACK, arithmetic, etc
moreSessionInfo <- function(print. = FALSE) {
    .M <- .Machine
    if(print.) str(.M[grep("^sizeof", names(.M))]) ## differentiate long-double..
    b64 <- .M$sizeof.pointer == 8
    onWindows <- .Platform$OS.type == "windows"
    ## Do we have 64bit but no-long-double ?
    arch <- Sys.info()[["machine"]]
    b64nLD <- (arch == "x86_64" && .M$sizeof.longdouble != 16)
    if(b64nLD) arch <- paste0(arch, "--no-long-double")
    if(print.)
        cat(sprintf("%d bit platform type '%s'  ==> onWindows: %s\narch: %s\n",
                    if(b64) 64 else 32, .Platform$OS.type, onWindows, arch))
    sInfo <- sessionInfo()
    if(!exists("osVersion")) osVersion <- sInfo$running
    if(print.) cat("osVersion (0):", osVersion, "\n")
    if(is.null(osVersion)) osVersion <- "Fedora" # very last resort
    if(!length(BLAS.is.LAPACK <- sInfo$BLAS == sInfo$LAPACK))
        BLAS.is.LAPACK <- NA # R versions <= 3.3.x
    ## A cheap check (that works on KH's debian-gcc setup, 2019-05):
    if(!length(BLAS.is.openBLAS <- grepl("openblas", sInfo$BLAS, ignore.case=TRUE)))
        BLAS.is.openBLAS <- NA
    if(!length(Lapack.is.openBLAS <- grepl("openblas", sInfo$LAPACK, ignore.case=TRUE)))
        Lapack.is.openBLAS <- NA
    if(print.)
        cat("osVersion:", osVersion, "\n"
          ,'+  BLAS "is" Lapack:', BLAS.is.LAPACK
          , '| BLAS=OpenBLAS:', BLAS.is.openBLAS
          , '| Lapack=OpenBLAS:', Lapack.is.openBLAS
          , "\n")
    ## NB: sessionInfo() really gets these:
    if(getRversion() >= "3.4") local({
        is.BLAS.LAPACK <- exists("La_library", mode="function") && ## R 3.4.0 and newer
            identical(La_library(), extSoftVersion()[["BLAS"]])
        stopifnot(isTRUE(is.BLAS.LAPACK == BLAS.is.LAPACK))
    })
    ## also TRUE for Windows [since both are "" !!]

    ## Find out if we are running Micrsoft R Open
    is.MS.Ropen <- {
        file.exists(Rpr <- file.path(R.home("etc"), "Rprofile.site")) &&
            length(lnsRpr <- readLines(Rpr)) &&
            ## length(grep("[Mm]icrosoft", lnsRpr)) > 3 # MRO 3.5.1 has '20' times "[Mm]icrosoft"
            length(grep("Microsoft R Open", lnsRpr, fixed=TRUE, value=TRUE)) > 0 ## MRO 3.5.1 has it twice
    }
    if(print. && is.MS.Ropen) cat("We are running 'Microsoft R Open'\n")

    ## I'd really would want to know which of (OpenBLAS | ATLAS | MKL | R's own BLAS+LAPACK)
    ##
    ## Next best, I'd really like
    ##
    ##    strictR <- we_are_using_Rs_own_BLAS_and_Lapack()  [ ==> BLAS != Lapack ]
    ##
    ## Actually the following aims to be equivalent to {and *is* for MM on Fedora, 2019-03}
    ## strictR <- !(using ATLAS || OpenBLAS || MKL )
    if(TRUE) {
        strictR <-
            !BLAS.is.LAPACK   && !is.MS.Ropen		&&
            !BLAS.is.openBLAS && !Lapack.is.openBLAS	&&
            TRUE
    } else { ## workaround:
        strictR <- print(Sys.info()[["user"]]) == "maechler"# actually
        ## but not when testing with /usr/bin/R [OpenBLAS on Fedora!] (as "maechler"):
        if(strictR && substr(osVersion, 1,6) == "Fedora" && R.home() == "/usr/lib64/R")
            strictR <- FALSE
    }
    if(print.) cat("strictR:", strictR, "\n")
    structure(class = "moreSessionInfo",
              list(
                  arch = arch
                , b64 = b64 # 64-bit (:<==> sizeof.pointer == 8 )
                , b64nLD = b64nLD # 64-bit, but --no-long-double (sizeof.longdouble != 16)
                , BLAS.is.LAPACK = BLAS.is.LAPACK
                , BLAS.is.openBLAS = BLAS.is.openBLAS
                , Lapack.is.openBLAS = Lapack.is.openBLAS
                , is.MS.Ropen = is.MS.Ropen # is R a version of Microsoft R Open (==> MKL-linked BLAS)
                , onWindows = onWindows
                , osVersion = osVersion
                , strictR = strictR # are BLAS & Lapack from R's source, and "otherwise known safe platform"
              ))
}

if(getRversion() < "3.4.0") withAutoprint <- function(x, ...) x
if(isTRUE(getOption("chk.moreSessionInfo"))) withAutoprint({
    ms1 <- moreSessionInfo()
    ms. <- moreSessionInfo(print. = TRUE)
    stopifnot(is.list(ms1), length(ms1) > 1, identical(ms1, ms.) )
})
