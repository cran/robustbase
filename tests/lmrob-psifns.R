#### Tests psi(), chi(),... etc and  tuning.psi, tuning.chi :

library(robustbase)
source(system.file("xtraR/plot-psiFun.R", package = "robustbase", mustWork=TRUE))

### (1) Test the functions themselves --------------------------------
pdf("rob-psifns.pdf")

## Simple version, no error checking, no derivative, nothing:
psiGGW <- function(x, a,b,c) {
    ifelse((ax <- abs(x)) < c,
           x,
           ifelse((ea <- -((ax-c)^b)/(2*a)) < -708.4, 0, x * exp(ea)))
}
stopifnot(all.equal(Mpsi  (5:9, cc=c(0,a=1/8,b=2,c=1/8,NA), "GGW"),
		    psiGGW(5:9,	       a=1/8,b=2,c=1/8), tol = 1e-13))


## Check that psi(<empty>)  |->  <empty>  works; ditto for +-Inf, NA,..
cG <- c(-.5, 1, .95, NA) # one of the 6 "builtin"s
d0 <- numeric()
IoI <- c(-Inf, 0, Inf)
NN <- c(NaN, NA)

cGs <- list(  c(-.4, 1.5,    0.85,  NA)
            , c(-.4, 1.5 ,   0.90,  NA)
            , c(-.4, 1.5 ,   0.95,  NA)
            , c(-.4, 1.5,    0.975, NA)
            , c(-.4, 1.5,    0.99 , NA)
            , c(-.4, 1.5,    0.995, NA)
            ##
            , c(-.4, 1.25,   0.975, NA)
            , c(-.4, 1.1,    0.975, NA)
            , c(-.4, 1.025,  0.975, NA)
            , c(-.4, 1.0125, 0.975, NA)
            ##
            ## FIXME , c(-.1, 1.25, 0.95, NA)
            ## FIXME , c(-.1, 1.25, 0.99, NA)
            )
st <- system.time(
cG.cnst <- lapply(cGs, function(cc)
                  lmrob.control(psi = "ggw", tuning.psi = cc)$tuning.psi)
)
cat('Time for constants computation of tuning.psi: ', st,'\n')
cGct <- t(sapply(cG.cnst, attr, "constants"))[,-1]
colnames(cGct) <- c("a","b","c", "rhoInf")
signif(cGct, 4)
stopifnot(all.equal(sapply(cG.cnst, function(cc) MrhoInf(cc, "ggw")),
                    cGct[,"rhoInf"]))


## Do these checks for a *list* of (c.par, psi) combinations:
c.psi.list <- list(
    list(cG, "GGW"),
    list(c(2,4,8), "Hampel"),
    list(c(1.5,3.5,8)*0.90, "Hampel"),
    list(par=c(-.5,1.5,.95,NA), "lqq"),
    list(bcs=c(1, 1, 1.25), "lqq"),
    list(1.1, "optimal"),
    list(0.1, "optimal"),
    list(2.3, "Welsh")
    )

for(c.psi in c.psi.list) {
    tPar <-  c.psi[[1]]; psi <- c.psi[[2]]
    stopifnot(is.numeric(tPar), is.character(psi))
    cat("Psi function ", psi,"; tuning par. c[]= (",
        paste(formatC(tPar, width=1), collapse=", "),")\n")
    for(FUN in list(Mpsi, Mchi, Mwgt))
	stopifnot(identical(d0, FUN(d0, tPar, psi=psi)),
                  identical(NN, FUN(NN, tPar, psi=psi)))
    stopifnot(identical(c(0,0,0), Mpsi(IoI, tPar,psi=psi)),
              identical(c(1,0,1), Mchi(IoI, tPar,psi=psi)),
              identical(c(0,1,0), Mwgt(IoI, tPar,psi=psi)))
    cat("chkPsi..(): ")
    chkPsi..(c(-5, 10), psi=psi, par=tPar)
    cat(" [Ok]\n------------------------\n\n")
}


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
p.psiFun2(x., "optimal", par=1.3)
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
ctr7 <- lmrob.control(psi = 'ggw', tuning.psi = c(-.3, 1.4, 0.95, NA),
                      tuning.chi = c(-0.3, 1.4, NA, 0.5))
ctr7$tuning.psi ## -> "constants"
ctr7$tuning.chi
summary(mp7 <-lmrob(Y ~ ., data = aircraft, control = ctr7))

set.seed(9)
summary(mp8 <- update(mp0, psi = 'lqq'))

set.seed(10)
ctr9 <- lmrob.control(psi = 'lqq', tuning.psi = c(ctr7$tuning.psi),
                      tuning.chi = c(ctr7$tuning.chi))
ctr9$tuning.psi
ctr9$tuning.chi
## Confirm these constants above (against the ones we got earlier)
## by recomputing them using higher accuracy :
(tpsi. <- robustbase:::lmrob.lqq.findc(ctr9$tuning.psi, rel.tol=1e-11, tol=1e-8))
(tchi. <- robustbase:::lmrob.lqq.findc(ctr9$tuning.chi, rel.tol=1e-11, tol=1e-8))
(tol4 <- .Machine$double.eps^.25)
stopifnot(all.equal(attr(ctr9$tuning.psi, "constants"), tpsi., tol=tol4),
	  all.equal(attr(ctr9$tuning.chi, "constants"), tchi., tol=tol4))

summary(mp9 <-lmrob(Y ~ ., data = aircraft, control = ctr9))


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
