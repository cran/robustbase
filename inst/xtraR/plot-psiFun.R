## Functions to plot and check psi-functions
## used in ../../tests/lmrob-psifns.R,
##	   ../../tests/psi-rho-etc.R
##     and ../../vignettes/psi_functions.Rnw  vignette


## Original Author of functions: Martin Maechler, Date: 13 Aug 2010, 10:17

p.psiFun <- function(x, psi, par, main=FALSE, ...)
{
    m.psi <- cbind(rho  = Mpsi(x, par, psi,deriv=-1),
                   psi  = Mpsi(x, par, psi,deriv= 0),
                   Dpsi = Mpsi(x, par, psi,deriv= 1),
                   wgt  = Mwgt(x, par, psi))
    robustbase:::matplotPsi(x, m.psi, psi=psi, par=par, main=main, ...) ## -> cbind(x, m.psi)
}
p.psiFun2 <- function(x, psi, par, main="short", ...)
    p.psiFun(x, psi, par, main=main, leg.loc= "bottomright", ylim = c(-2.2, 6))
## for psi_func class objects: simply use plot() method.

mids <- function(x) (x[-1]+x[-length(x)])/2

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)
identical5 <- function(a,b,c,d,e) identical(a,b) && identical4(b,c,d,e)

##' @title Check consistency of psi/chi/wgt/.. functions
##' @param m.psi matrix as from p.psiFun()
##' @param tol
##' @return concatenation of \code{\link{all.equal}} results
##' @author Martin Maechler
chkPsiDeriv <- function(m.psi, tol = 1e-4) {
    stopifnot(length(tol) > 0, tol >= 0,
              is.numeric(psi <- m.psi[,"psi"]),
              is.numeric(dx  <- diff(x <- m.psi[,"x"])))
    if(length(tol) < 2) tol[2] <- 10*tol[1]
    xn0 <- abs(x) > 1e-5
    c(all.equal(mids(psi), diff(m.psi[,"rho"])/dx, tol=tol[1]), # rho'  == psi
      all.equal(mids(m.psi[,"Dpsi"]), diff(psi)/dx, tol=tol[2]),# psi'  == psip
      all.equal(m.psi[xn0,"wgt"], (psi/x)[xn0], tol= tol[1]/10))# psi/x == wgt
}

##' @title Check consistency of psi/chi/wgt/.. functions
##' @param x  range or vector of abscissa values
##' @param psi psi() function spec., passed to  M.psi() etc
##' @param par tuning parameter,     passed to  M.psi() etc
##' @param tol tolerance for equality checking of numeric derivatives
##' @return concatenation of \code{\link{all.equal}} results
##' @author Martin Maechler
chkPsi.. <- function(x, psi, par, tol = 1e-4, quiet=FALSE) {
    stopifnot(length(tol) > 0, tol >= 0, is.numeric(x), is.finite(x))
    if(length(x) == 2) ## it is a *range* -> produce vector
        x <- seq(x[1], x[2], length = 1025L)
    dx <- diff(x)
    x0 <- sort(x)
    x <- c(-Inf, Inf, NA, NaN, x0)
    rho  <- Mpsi(x, par, psi, deriv=-1)
    psix <- Mpsi(x, par, psi, deriv= 0)
    Dpsi <- Mpsi(x, par, psi, deriv= 1)
    wgt  <- Mwgt(x, par, psi)

    chi  <- Mchi(x, par, psi)
    chi1 <- Mchi(x, par, psi, deriv=1)
    chi2 <- Mchi(x, par, psi, deriv=2)
    rho.Inf <- MrhoInf(par, psi)
    stopifnot(all.equal(rep(rho.Inf,2), rho[1:2]),
	      all.equal(chi, rho  / rho.Inf),
	      all.equal(chi1,psix / rho.Inf),
	      all.equal(chi2,Dpsi / rho.Inf)
	      )

    D2psi <- tryCatch(Mpsi(x, par, psi, deriv= 2), error=function(e)e)
    has2 <- !inherits(D2psi, "error")
    if(!quiet & !has2) message("Not checking psi''() := Mpsi(*, deriv=2)")
    stopifnot(is.numeric(psix),
              ## check NA / NaN :
              identical5(x[3:4], rho[3:4], psix[3:4], Dpsi[3:4], wgt[3:4]),
              if(has2) identical(x[3:4], D2psi[3:4]) else TRUE)

    if(length(tol) < 2) tol[2] <- 10*tol[1]
    i <- 5:length(x) # leaving away the first 4 (+-Inf, NA..)
    xn0 <- is.finite(x) & abs(x) > 1e-5
    c(all.equal(mids(psix[i]), diff(rho [i])/dx, tol=tol[1]), # rho'  == psi
      all.equal(mids(Dpsi[i]), diff(psix[i])/dx, tol=tol[2]), # psi'  == psip
      all.equal(wgt[xn0], (psix/x)[xn0], tol= tol[1]/10),     # psi/x == wgt
      if(has2)
      all.equal(mids(D2psi[i]),diff(Dpsi[i])/dx,tol=8*tol[2])# (psi')' == D2psi
      else TRUE)
}

