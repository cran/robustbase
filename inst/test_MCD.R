#### Utility functions for testing covMCD()
#### -------------------------------------- ../tests/tmcd.R

repMCD <- function(x, nrep = 1, method = c("FASTMCD","MASS"))
{
    stopifnot(length(nrep) == 1, nrep >= 1)
    method <- match.arg(method)
    if(method == "MASS")
         for(i in 1:nrep) MASS::cov.mcd(x)
    else for(i in 1:nrep) covMcd(x)
}

doMCDdata <- function(nrep = 1, time = nrep >= 3, short = time, full = !short,
		   method = c("FASTMCD", "MASS"))
{
    ##@bdescr
    ## Test the function covMcd() on the literature datasets:
    ##
    ## Call covMcd() for "all" regression datasets available in robustbase
    ## and print:
    ##  - execution time (if time)
    ##  - objective function
    ##  - best subsample found (if not short)
    ##  - outliers identified (with cutoff 0.975) (if not short)
    ##  - estimated center and covariance matrix (if full)
    ##
    ##@edescr
    ##
    ##@in  nrep              : [integer] number of repetitions to use for estimating the
    ##                                   (average) execution time
    ##@in  time              : [boolean] whether to evaluate the execution time
    ##@in  short             : [boolean] whether to do short output (i.e. only the
    ##                                   objective function value). If short == FALSE,
    ##                                   the best subsample and the identified outliers are
    ##                                   printed. See also the parameter full below
    ##@in  full              : [boolean] whether to print the estimated center and covariance matrix
    ##@in  method            : [character] select a method: one of (FASTMCD, MASS)


    domcd <- function(x, xname, nrep = 1)
    {
        n <- dim(x)[1]
        p <- dim(x)[2]
        if(method == "MASS") {
            mcd <- MASS::cov.mcd(x)
            quan <- as.integer(floor((n + p + 1)/2)) #default: floor((n+p+1)/2)
        }
        else {
            mcd <- covMcd(x) # trace = FALSE
            quan <- as.integer(mcd$quan)
        }

        crit <- if(method == "MASS") mcd$crit else log(mcd$crit)

        xres <- sprintf("%*s %3d %3d %3d %12.6f", lname, xname, n, p, quan, crit)
        if(time) {
            xtime <- system.time(repMCD(x, nrep, method))[1]/nrep
            xres <- sprintf("%s %10.1f", xres, 1000 * xtime)
        }
        cat(xres, "\n")

        if(!short) {
            cat("Best subsample: \n")
            print(mcd$best)

            ibad <- which(mcd$mcd.wt == 0)
            names(ibad) <- NULL
            nbad <- length(ibad)
            cat("Outliers: ",nbad,"\n")
            if(nbad > 0)
                print(ibad)
            if(full) {
                cat("-------------\n")
                print(mcd)
            }
            cat("--------------------------------------------------------\n")
        }
    }

    lname <- 20
    method <- match.arg(method)

    data(heart)
    data(phosphor)
    data(starsCYG)
    data(stackloss)
    data(coleman)
    data(salinity)
    data(wood)
    data(hbk)

    data(Animals, package = "MASS")
    brain <- Animals[c(1:24, 26:25, 27:28),]
    data(milk)
    data(bushfire)

    ##    data(x1000)
    ##    data(x5000)

    tmp <- sys.call()
    cat("\nCall: ", deparse(substitute(tmp)),"\n")

    cat("Data Set               n   p  Half   LOG(obj)  Time [ms]\n")
    cat("========================================================\n")
    domcd(heart[, 1:2], data(heart), nrep)
    domcd(data.matrix(subset(phosphor, select = -plant)),
          data(phosphor), nrep)
    domcd(starsCYG, data(starsCYG), nrep)
    domcd(stack.x, data(stackloss), nrep)
    domcd(data.matrix(subset(coleman, select = -Y)), data(coleman), nrep)
    domcd(data.matrix(subset(salinity, select = -Y)), data(salinity), nrep)
    domcd(data.matrix(subset(wood, select = -y)), data(wood), nrep)
    domcd(data.matrix(subset(hbk,  select = -Y)), data(hbk), nrep)

    domcd(brain, "Animals", nrep)
    domcd(milk, data(milk), nrep)
    domcd(bushfire, data(bushfire), nrep)
    cat("========================================================\n")
    ##    domcd(x1000$X,data(x1000), nrep)
    ##    domcd(x5000$X,data(x5000), nrep)
}
