#### Psi(), Rho(), weight() etc functions for  M-Estimation and extensions

## Use an  S4 class for such function classes
## Follow a similar idea as  nlsModel() {in "stats"} which returns
##  a list of functions sharing a common {non-small!} environment

## NOTA BENE:  Experiments etc are currently in ../experi-psi-rho-funs.R
## ---------                                    ~~~~~~~~~~~~~~~~~~~~~~~~

## ---> look for 'FIXME' below !!!
##               -------

### A.  (Symmetric) Location / Regression

## A single   function(x, tuningPars)
## a. 1st argument 'x', numeric; must work vectorized on x
## b. further arguments: tuning parameters *with a default*
setClass("functionX", contains = "function",
         validity = function(object) {
             ## "function" is already because of 'contains'
             if(names(ff <- formals(object))[1] != "x")
                 return("first argument must be 'x'")
             f0 <- object(0)
             fI <- object(Inf)
             if(!identical(c(f0,fI), object(c(0,Inf))))
                 return("F(x, *) does not vectorize in 'x'")
             ## Otherwise : valid
             TRUE
         })

## A functional --- i.e. function of "tuning pars only",  such as
##  Ep(hc) = Int_{-Inf}^{+Inf} psi(x; hc)^2 dnorm(x) dx
setClass("functionXal", contains = "function",
         validity = function(object) {
             f0 <- object(0)
             fI <- object(Inf)
             if(!identical(c(f0,fI), object(c(0,Inf))))
                 return("F(x, *) does not vectorize in 'x'")
             ## Otherwise : valid
             TRUE
         })

setClass("psi_func",
         representation(rho = "functionX",
                        psi = "functionX", ## psi(x) == d/dx rho(x) = x * wgt(x)
                        wgt = "functionX", ## wgt(x) == psi(x) / x
                        Dpsi = "functionX",## psi'(x) == d/dx psi(x) = rho''(x)
                        ## tuning parameters, i.e., formals(rho)[-1]
                        tDefs = "numeric",## *named* values of tuning parameters
                        ## FIXME !! {see 4 lines below}
                        Erho =  "functionXal", # = E_X[rho(X)];   X~N(0,1);
                        Epsi2 = "functionXal", # = E_X[psi(X)^2]; X~N(0,1); "A"
                        EDpsi = "functionXal", # = E_X[psi'(X)];  X~N(0,1); "B"

                        xtras = "list" ## for flexible extensions..
                        ))
## FIXME: need other E[] than just wrt N(0,1)
## -----  e.g. for robglm(), need  E[] wrt Gamma(.)

### Constructors / "Examples" [the examples are the objects, we'll really use!]

psiFunc <- function(rho,psi,wgt, Dpsi, Erho=NULL, Epsi2=NULL, EDpsi=NULL, ...)
{
    lent <- length(dotsargs <- list(...))
    ## '...'  must contain all tuning parameters and their defaults:
    stopifnot(length(nt <- names(dotsargs)) == lent,
              all(nchar(nt)) >= 1)
    if(lent >= 1) {
        ## rho, psi,... checking: must have argument names
        argn <- c("x", nt)
        for(fnam in list("rho", "psi", "wgt", "Dpsi")) {
            f <- get(fnam, inherits = FALSE)
            ef <- environment(f)
            nf <- names(ff <- formals(f)) # "x" and "k" for Huber's
            if(!identical(nf, argn))
                stop("arguments of function '",fnam,"' are (",
                     paste(nf,  collapse=","),") but should be (",
                     paste(argn,collapse=","),").")

            formals(f)[-1] <- dotsargs
            environment(f) <- ef
            assign(fnam, f, inherits = FALSE)
        }
    }

    new("psi_func",
        rho = new("functionX", rho),
        psi = new("functionX", psi),
        wgt = new("functionX", wgt),
        Dpsi= new("functionX", Dpsi),
        ## tNams = if(lent) nt else character(0),
        tDefs = if(lent) unlist(dotsargs) else numeric(0),
        Erho= new("functionXal", Erho),
        Epsi2= new("functionXal", Epsi2),
        EDpsi= new("functionXal", EDpsi)
        )
}

chgDefaults <- function(object, ...)
    standardGeneric("chgDefaults")

setMethod("chgDefaults", signature("psi_func"),
          function(object, ...) {
              lent <- length(dotsargs <- list(...))
              ## '...'  must contain all tuning parameters and their defaults:
              stopifnot(length(nt <- names(dotsargs)) == lent,
                        all(nchar(nt)) >= 1)
              if(lent >= 1) {
                  ## rho "..." must conform to rho, etc:
                  nf <- names(ff <- formals(object@rho))
                  if(!identical(nf[-1], nt))
                     stop("invalid tuning parameter names: ",
                          paste(nt,    collapse=",")," instead of ",
                          paste(nf[-1],collapse=","),".")

                  for(fnam in list("rho", "psi", "wgt", "Dpsi")) {
                      f <- slot(object, fnam)
                      ef <- environment(f)
                      formals(f)[-1] <- dotsargs
                      environment(f) <- ef
                      ## lowlevel {faster than}: slot(..) <- new("functionX", f)
                      slot(object, fnam)@.Data <- f
                  }
                  object@tDefs <- unlist(dotsargs)
              }
              object
          })

## Huber:
huberPsi <- psiFunc(rho =
                  function(x, k) {
                      r <- u <- abs(x); I <- u < k
                      r[ I] <- u[I]^2 / 2
                      r[!I] <- k*(u[!I] - k / 2)
                      r
                  },
                  psi  = function(x, k) pmin2(k, pmax2(-k, x)),
                  wgt  = function(x, k) pmin2(1, k/abs(x)),
                  Dpsi = function(x, k) abs(x) <= k,
                  Erho = function(k) {iP <- pnorm(k, lower=FALSE)
                                      1/2 - iP + k*(dnorm(k) - k*iP)},
                  Epsi2= function(k) ifelse(k < 10,
                  1 - 2*(k*dnorm(k) + (1-k*k)*pnorm(k, lower=FALSE)), 1),
                  EDpsi= function(k) 2*pnorm(k) - 1,
                  ## the tuning pars and default:
                  k = 1.345)

## Hampel: --------- FIXME -----------
if(FALSE)
hampelPsi <-
    psiFunc(rho =
            function(x, k) {
### UNFINISHED
                u <- abs(x)
                lrg <- k[3] <= u
                I <- u < k[1]
                mid <- !I & !lrg    # contains constant and descending
                r <- numeric(length(x))
                r[ I] <- u[I]^2 / 2
                r[mid] <- k*(u[!I] - k / 2)
                r
            },
            psi = function(x, k) {
                ## this is "optimized" ==> factors faster than using ifelse()!
                u <- abs(x)
                lrg <- k[3] <= u
                mid <- k[1] < u & !lrg # constant _and_ descending
                ## x is result for |x| < k[1]
                x[lrg] <- 0
                if(any(mid))
                    x[mid] <- k[1] * sign(x[mid])*
                        pmin2(1, (u[mid] - k[3])/(k[2] - k[3]))
                x
            },
            wgt  = function(x, k) {
                x <- abs(x)
                lrg <- k[3] <= x
                I <- x < k[1]
                mid <- !I & !lrg    # contains constant and descending
                x[I] <- 1
                x[lrg] <- 0
                if(any(mid))
                    x[mid] <- k[1] / x[mid] *
                        pmin2(1, (x[mid] - k[3])/(k[2] - k[3]))
                x
            },
            Dpsi = function(x, k) {
            },
            Erho = function(k) {
            },
            Epsi2= function(k) {
            },
            EDpsi= function(k) {
            },
            ## the tuning pars and default:
            k = c(2,4,8) / 1.345)

## TODO:  Biweight :
if(FALSE)
tukeyPsi <- c() ##########


## maybe TODO: Optimal tanh() estimator for location



### B.  M-Estimators of Scale --- need chi() and slightly different functionals
### --- ----------------------
###
## one "challenge" is the  a(b)  needed in  chi(x; a,b) = [x^2 -1 -a]_b^b
## for  V-optimal  M-Estimates of scale
## --> but that's solved (!) in ./scale-chi-opt.R
##                              ~~~~~~~~~~~~~~~~~
## Then, I'd also want the optimal chi for s
