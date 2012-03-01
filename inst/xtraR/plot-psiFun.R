## Functions to plot and check psi-functions
## used in ../../tests/lmrob-psifns.R,
##	   ../../tests/psi-rho-etc.R
##     and ../doc/psi_functions.Rnw  vignette

psiF <- robustbase:::lmrob.psifun # deriv = -1 (rho), 0, 1
chiF <- robustbase:::lmrob.chifun # rho(.) normalized to max|.| = 1;  deriv
wgtF <- robustbase:::lmrob.wgtfun
plot.psiFun <- robustbase:::plot.psiFun

## Original Author of functions: Martin Maechler, Date: 13 Aug 2010, 10:17

p.psiFun <- function(x, psi, par, ...)
{
    m.psi <- cbind(rho  = psiF(x, par, psi,deriv=-1),
                   psi  = psiF(x, par, psi,deriv= 0),
                   Dpsi = psiF(x, par, psi,deriv= 1),
                   wgt  = wgtF(x, par, psi))
    plot.psiFun(x, m.psi, psi, par, ...)
}
p.psiFun2 <- function(x, psi, par, ...)
    p.psiFun(x, psi, par, shortMain = TRUE,
             leg.loc = "bottomright", ylim = c(-3, 5))
## for psi_func class objects: simply use plot() method.

mids <- function(x) (x[-1]+x[-length(x)])/2
chkPsiDeriv <- function(m.psi, tol = 1e-4) {
    ## m.psi: matrix as from p.psiFun()
    stopifnot(length(tol) > 0, tol >= 0,
              is.numeric(psi <- m.psi[,"psi"]),
              is.numeric(dx  <- diff(x <- m.psi[,"x"])))
    if(length(tol) < 2) tol[2] <- 10*tol[1]
    xn0 <- abs(x) > 1e-5
    c(all.equal(mids(psi), diff(m.psi[,"rho"])/dx, tol=tol[1]), # rho'  == psi
      all.equal(mids(m.psi[,"Dpsi"]), diff(psi)/dx, tol=tol[2]),# psi'  == psip
      all.equal(m.psi[xn0,"wgt"], (psi/x)[xn0], tol= tol[1]/10))# psi/x == wgt
}
