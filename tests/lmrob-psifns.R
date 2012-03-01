#### Tests psi(), chi(),... etc and  tuning.psi, tuning.chi :

library(robustbase)
source(system.file("xtraR/plot-psiFun.R", package = "robustbase", mustWork=TRUE))

### (1) Test the functions themselves --------------------------------
pdf("rob-psifns.pdf")

psiF <- robustbase:::lmrob.psifun # deriv = -1 (rho), 0, 1
chiF <- robustbase:::lmrob.chifun # rho(.) normalized to max|.| = 1;  deriv
wgtF <- robustbase:::lmrob.wgtfun

## Simple version, no error checking, no derivative, nothing:
psiGGW <- function(x, a,b,c) {
    ifelse((ax <- abs(x)) < c,
           x,
           ifelse((ea <- -((ax-c)^b)/(2*a)) < -708.4, 0, x * exp(ea)))
}
stopifnot(all.equal(psiF  (5:9, cc=c(0,a=1/8,b=2,c=1/8,NA), "GGW"),
		    psiGGW(5:9,	       a=1/8,b=2,c=1/8), tol = 1e-13))


funs <- list(psiF, chiF, wgtF)
## Check that psi(<empty>)  |->  <empty>  works
cG <- c(-.5,1,.95,NA)
d0 <- numeric()
IoI <- c(-Inf, 0, Inf)
## TODO: Do these checks for a *list* of combinations such as  (cG, "GGW"):
## ^^^^^
for(FUN in funs)
    stopifnot(identical(d0, FUN(d0, cG, "GGW")))
stopifnot(identical(c(0,0,0), psiF(IoI, cG,"GGW")),
	  identical(c(1,0,1), chiF(IoI, cG,"GGW")),
	  identical(c(0,1,0), wgtF(IoI, cG,"GGW")))

## FIXME: Check  chiF() <-> psiF(*, deriv = -1)


## Nice plots -- and check derivatives ----

head(x. <- seq(-5, 10, length=1501))
## [separate lines, for interactive "play": ]
stopifnot(chkPsiDeriv(p.psiFun(x., "LQQ", par=c(-.5,1.5,.95,NA))))
stopifnot(chkPsiDeriv(p.psiFun(x., "GGW", par= cG)))
stopifnot(chkPsiDeriv(p.psiFun(x., "optimal", par=2)))
stopifnot(chkPsiDeriv(p.psiFun(x., "Hampel",
                               par = ## Default, but rounded:
                               round(c(1.5, 3.5, 8) * 0.9016085, 1)),
                      tol = 1e-3))

stopifnot(chkPsiDeriv(p.psiFun(x., "biweight", par = 4)))
stopifnot(chkPsiDeriv(p.psiFun(x., "Welsh", par = 1.5)))

## The same 6, all in one plot:
op <- par(mfrow=c(3,2), mgp = c(1.5, .6, 0), mar = .1+c(3,3,2,.5))
p.psiFun2(x., "LQQ", par=c(-.5,1.5,.95,NA))
p.psiFun2(x., "GGW", par= cG)
p.psiFun2(x., "optimal", par=2)
p.psiFun2(x., "Hampel", par = round(c(1.5, 3.5, 8) * 0.9016085, 1))
p.psiFun2(x., "biweight", par = 4)
p.psiFun2(x., "Welsh", par = 1.5)
par(op)

### (2) Test them as  arguments of  lmrob() or  lmrob.control(): -----

data(aircraft)

set.seed(1)
summary(mp0 <- lmrob(Y ~ ., data = aircraft, psi = 'bisquare', method = 'SMDM'))

set.seed(2)
summary(mp1 <- update(mp0, psi = 'optimal'))

set.seed(3)
summary(mp2 <- update(mp0, psi = 'ggw'))

set.seed(4)
summary(mp3 <- update(mp0, psi = 'welsh'))

set.seed(5)
summary(mp4 <- update(mp0, psi = 'ggw', tuning.psi = c(-.5, 1.5, 0.85, NA),
                      tuning.chi = c(-0.5, 1.5, NA, 0.5)))

set.seed(6)
summary(mp5 <- update(mp0, psi = 'ggw', tuning.psi = c(-.5, 1.0, 0.95, NA),
                      tuning.chi = c(-0.5, 1.0, NA, 0.5)))

set.seed(7)
summary(mp6 <- update(mp0, psi = 'hampel'))

set.seed(8)
ctrl <- lmrob.control(psi = 'ggw', tuning.psi = c(-.3, 1.4, 0.95, NA),
                      tuning.chi = c(-0.3, 1.4, NA, 0.5))
ctrl$tuning.psi ## -> "constants"
ctrl$tuning.chi
summary(mp7 <-lmrob(Y ~ ., data = aircraft, control = ctrl))

set.seed(9)
summary(mp8 <- update(mp0, psi = 'lqq'))

set.seed(10)
ctrl <- lmrob.control(psi = 'lqq', tuning.psi = c(-.3, 1.4, 0.95, NA),
                      tuning.chi = c(-0.3, 1.4, NA, 0.5))
ctrl$tuning.psi
ctrl$tuning.chi
summary(mp7 <-lmrob(Y ~ ., data = aircraft, control = ctrl))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
