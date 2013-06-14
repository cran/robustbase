#### These are experiments, using the definitions in
#### ==================
#### ./R/psi-rho-funs.R  <<<<<<<< see there
#### ==================


## NB:  "Recycle" code in
## -- /u/maechler/R/MM/STATISTICS/robust/psi-funs.R
##      ~~~~~~~~~~  and
##    /u/maechler/R/MM/STATISTICS/robust/wt-funs.R
##      ~~~~~~~~~~

Hf0 <- function(x, c=1.35) pmax.int(-c, pmin.int(x,c))
Hf <- new("functionX", Hf0)
stopifnot(validObject(Hf)) # ok !


### psiFunc() examples
## classical {trivial, not interesting}:
F1 <- function(x) rep.int(1, length(x))
cPsi <- psiFunc(rho = function(x) x^2 / 2, psi = function(x) x,
                wgt = F1, Dpsi = F1,
                Erho = function(x) rep.int(1/2, length(x)),
                Epsi2 = F1, EDpsi = F1)


### MASS  -- ?rlm --- has
##
##-      psi.huber   (u, k = 1.345, deriv = 0)
##-      psi.hampel  (u, a = 2, b = 4, c = 8, deriv = 0)
##-      psi.bisquare(u, c = 4.685, deriv = 0)
##   where deriv = 0 :  psi(x)/x i.e.  'wgt'

## MM has more in psi-funs.R (see above)



## Reproduce Table 1, p.138 of Hampel, Rousseeuw, Ronchetti, Stahel (1986):
b <- c(seq(0, 3, by = 0.1), 4, 5, Inf)
A <- hubPsi@Epsi2(b)
B <- hubPsi@EDpsi(b)

options(width = 100, digits = 7)

cbind(b = b, A = A, B = B, V = A/B^2, e = B^2/A,
      gamma. = b/B, k. = 1 + b^2/A, lambda. = 1/B)


(hubPsi2 <- chgDefaults(hubPsi, k = 2)) # works too!

## TODO:  Hampel,  Biweight

## Experiments only:--- delete later:
PL <- list(rho =
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
                        pmin.int(1, (u[mid] - k[3])/(k[2] - k[3]))
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
                        pmin.int(1, (x[mid] - k[3])/(k[2] - k[3]))
                x
            }
           )
plot(function(x) PL$psi(x, k = c(2,4,8)), -10, 10)
## Compare psi(x) / x  with wgt(x) :
plot(function(x) PL$psi(x, k = c(2,4,8))/x, -10, 10)
curve(PL$wgt(x, k = c(2,4,8)), add = TRUE, col="pink")
abline(v = outer(c(-1,1),c(2,4,8)), col="light blue", lty = "1F")
## -> seems fine

ppsi <- PL$psi; formals(ppsi)[-1] <- list(k = c(2,4,8))
str(ppsi)
plot(ppsi, -10, 10) # all looks fine ...
pp <- new("functionX", ppsi)
## now ok

pwgt <- PL$wgt; formals(pwgt)[-1] <- list(k = c(2,4,8))
str(pwgt)
plot(pwgt, -10, 10) # all looks fine ...
pp <- new("functionX", pwgt)## now ok

prho <- PL$rho; formals(prho)[-1] <- list(k = c(2,4,8))
str(prho)
plot(prho, -10, 10) # all looks fine ...
pp <- new("functionX", prho)## now ok

###--- Compute  E[rho] empirically for Huber:
rho <- function (x, k)
{
    r <- u <- abs(x)
    I <- u < k
    r[I] <- u[I]^2/2
    r[!I] <- k * (u[!I] - k/2)
    r
}

kk <- c(0.01, 0.1, seq(0.5, 1.5, by = 0.1), 1.7, seq(2, 4, by = 0.5))
ErhoA <- lapply(kk, function(k)
                integrate(function(x)rho(x, k=k) * dnorm(x), -Inf,Inf,
                          rel.tol = 1e-12))
str(Erho <- sapply(ErhoA, "[[", "value"))
plot(Erho ~ kk, type = 'o', col = 2)

## Hampel:
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
                        pmin.int(1, (u[mid] - k[3])/(k[2] - k[3]))
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
                        pmin.int(1, (x[mid] - k[3])/(k[2] - k[3]))
                x
            },
            Dpsi = function(x, k) {
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

## to use for other types, just change 'ggw' argument
## standardized to have Dpsi(0) = 1
## to have rho(inf) = 1 use .M.chi instead (as well as deriv + 1)
## using this results in an error while preparing for lazy loading:
## (MM, MK: the error arises during the validity check)
## ** preparing package for lazy loading
## Creating a generic function from function "chgDefaults"
## Error in .M.psi(x, k, "ggw", -1) : object 'R_psifun' not found
## Error : unable to load R code in package 'robustbase'
## ERROR: lazy loading failed for package ‘robustbase’
## ('R_psifun' is the pointer to the C-function used in .M.psi)
ggwPsi <- psiFunc(rho = function(x, k) .M.psi(x, k, 'ggw', -1),
                  psi = function(x, k) .M.psi(x, k, 'ggw', 0),
                  wgt = function(x, k) .M.wgt(x, k, 'ggw'),
                  Dpsi = function(x, k) .M.psi(x, k, 'ggw', 1),
                  Erho = function(k) lmrob.E(.M.psi(r, k, 'ggw', -1),
                    use.integrate = TRUE),
                  Epsi2 = function(k) lmrob.E(.M.psi(r, k, 'ggw', 0)^2,
                    use.integrate = TRUE),
                  EDpsi = function(k) lmrob.E(.M.psi(r, k, 'ggw', 1),
                    use.integrate = TRUE),
                  k = c(-0.5, 1.5, 0.95, NA))

## maybe TODO: Optimal tanh() estimator for location


### --> scale - rho/psi functions
