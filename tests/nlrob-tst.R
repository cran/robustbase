commandArgs()
library(robustbase)

## for now: ----------------------------------------
if(file.exists(fil <- system.file("xtraR/test-tools.R",  package = "robustbase"))) {
    source(fil)
} else {
  identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
  identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)
  assert.EQ <- function(target, current, tol = if(showOnly) 0 else 1e-15,
                      giveRE = FALSE, showOnly = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## showOnly: if TRUE, return (and hence typically print) all.equal(...)
    T <- isTRUE(ae <- all.equal(target, current, tolerance = tol, ...))
    if(showOnly) return(ae) else if(giveRE && T) { ## don't show if stop() later:
	ae0 <- if(tol == 0) ae else all.equal(target, current, tolerance = 0, ...)
	if(!isTRUE(ae0)) writeLines(ae0)
    }
    if(!T) stop("all.equal() |-> ", paste(ae, collapse=sprintf("%-19s","\n")))
    else if(giveRE) invisible(ae0)
  }
}
if(FALSE) ## in future: ----------------------------------------
source(system.file("xtraR/test-tools.R",  package = "robustbase"))
##
## -> assert.EQ(), identical3(), ..


## From  <Rsrc>/tests/reg-tests-1d.R  (2022-08-09):
## A good guess if we have _not_ translated error/warning/.. messages:
## (should something like this be part of package tools ?)
englishMsgs <- {
    ## 1. LANGUAGE takes precedence over locale settings:
    if(nzchar(lang <- Sys.getenv("LANGUAGE")))
        lang == "en"
    else { ## 2. Query the  locale
        if(!onWindows) {
            ## sub() :
            lc.msgs <- sub("\\..*", "", print(Sys.getlocale("LC_MESSAGES")))
            lc.msgs == "C" || substr(lc.msgs, 1,2) == "en"
        } else { ## Windows
            lc.type <- sub("\\..*", "", sub("_.*", "", print(Sys.getlocale("LC_CTYPE"))))
            lc.type == "English" || lc.type == "C"
        }
    }
}
cat(sprintf("English messages: %s\n", englishMsgs))


DNase1 <- DNase[ DNase$Run == 1, ]
Y <- DNase1[,"density"] # for convenience below

## classical
fm1 <- nls(density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
	   data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), trace=TRUE)
summary(fm1)
wm1 <- update(fm1, weights = sqrt(conc)) # (weights as function of <var>)

## robust
rm1 <- nlrob(formula(fm1), data = DNase1, trace = TRUE,
	     start = list(Asym = 3, xmid = 0, scal = 1))
(sm1 <- summary(rm1))
stopifnot(all.equal(Y, fitted(fm1) + residuals(fm1), check.attributes=FALSE),
	  ## fitted(<nls>) has "label" attribute
	  identical3(c(fitted(fm1)), predict(fm1), predict(fm1, newdata=DNase1)),
	  ## robust fit :
	  identical3(fitted(rm1), predict(rm1), predict(rm1, newdata=DNase1)),
	  all.equal(Y, unname(fitted(rm1) + residuals(rm1))))
print(coef(rm1), digits=12)
## 2.35963008460 1.49945088410 1.04506391722	F19 Lx 64b
## 2.35963008460 1.49945088410 1.04506391722	Win(Serv.2003) 64b
## 2.35963008613 1.49945088600 1.04506391793	F19 Lx 32b
## 2.35963008613 1.49945088600 1.04506391793	Win(Serv.2003) 32b
assert.EQ(coef(rm1), giveRE=TRUE,
	  c(Asym=2.35963008, xmid=1.49945088, scal=1.04506392), tol = 4e-8)
assert.EQ(sqrt(diag(sm1$cov)), giveRE=TRUE,
	  ## 32b 0.08626872273,     0.0902194541,      0.03503833759
	  c(Asym=0.0862687305, xmid=0.0902194608, scal=0.0350383389),
	  tol = 7e-7)
## examples with weights:
rm. <- update(rm1, weights = NULL)# 'NULL' but not missing()
ww <- sqrt(DNase1[,"conc"])
wr1  <- update(rm1, weights = sqrt(conc), trace=FALSE)
wr1. <- update(rm1, weights = ww,         trace=FALSE)
ii <- names(rm1) != "call"
stopifnot(all.equal(rm1[ii], rm.[ii], tol = 1e-15),
          all.equal(wr1[ii],wr1.[ii], tol = 1e-15))

## When y has NA's: Add them to the end of DNase1
n <- nrow(DNase1)
DNAnase <- DNase1[c(1:n,1:3), ]
DNAnase$density[n+(1:3)] <- NA

naExp <- expression(na.fail = na.fail, # << the default in model.frame() which most use
                    na.omit = na.omit,
                    na.pass = na.pass,
                    na.exclude= na.exclude)

noAction <- function(fm) {
    if(("call" %in% names(fm)) && !is.null(fm$call))
        fm$call <- noAction(fm$call)
    if(!is.na(H <- match("na.action", names(fm))))
        fm[-H]
    else
        fm
}

rNls <- lapply(naExp, function(naAct)
    tryCatch(error = identity,
       nls(formula(fm1), data = DNAnase,
           na.action = eval(naAct),
           start = list(Asym = 3, xmid = 0, scal = 1))
       )
    )

str(sapply(rNls, class))
## $ na.fail   : chr [1:3] "simpleError" "error" "condition"
## $ na.omit   : chr "nls"
## $ na.pass   : chr [1:3] "simpleError" "error" "condition"
## $ na.exclude: chr "nls"

stopifnot(exprs = {
    inherits(rNls$na.pass, "simpleError")
    inherits(rNls$na.fail, "simpleError")
})

if(englishMsgs) stopifnot(exprs = {
    conditionMessage(rNls$na.fail) == "missing values in object"
    grepl("NA.*foreign function call", conditionMessage(rNls$na.pass))
})

sO <- summary(rNls$na.omit)
sX <- summary(rNls$na.exclude)
stopifnot(exprs = {
    identical(noAction(sO),
              noAction(sX))
    all.equal(coef(sm1), coef(sX), tol = 0.10) # 0.0857
})


robNls <- lapply(c(noCov=FALSE,doCov=TRUE), function(doCov)
    lapply(naExp, function(naAct)
        tryCatch(error = identity,
                 nlrob(formula(fm1), data = DNAnase,
                       na.action = eval(naAct),
                       doCov = doCov, ## <--
                       start = list(Asym = 3, xmid = 0, scal = 1))
                 )
        )
    )
## gives
## Warning messages:
## In old - new :
##   longer object length is not a multiple of shorter object length

str(sapply(robNls[["noCov"]], class))
## List of 4
##  $ na.fail   : chr [1:3] "simpleError" "error" "condition"
##  $ na.omit   : chr [1:3] "simpleError" "error" "condition"  <<< FIXME, should work as nls() does
##  $ na.pass   : chr [1:3] "simpleError" "error" "condition"
##  $ na.exclude: chr [1:2] "nlrob" "nls"
stopifnot(identical(sapply(robNls[["noCov"]], class),
                    sapply(robNls[["doCov"]], class)))

## same checks as for nls():
lapply(robNls, function(LL)
    stopifnot(exprs = {
        inherits(LL$na.pass, "simpleError")
        inherits(LL$na.fail, "simpleError")
    })) -> .tmp

if(englishMsgs) lapply(robNls, function(LL)
    stopifnot(exprs = {
        conditionMessage(LL$na.fail) == "missing values in object"
        ## different message than nls():
        grepl("missing.*weights not allowed", conditionMessage(LL$na.pass))
    })
    ) -> .tmp

## the only one which works currently:
robNxcl <- robNls[["noCov"]]$na.exclude

## Fix these : ===================
.t <- try( s.no <- summary(robNxcl) )
.t <- try( v.no <- vcov(robNxcl) )
## both give *same* error:
## Error in .vcov.m(object, Scale = sc, resid.sc = as.vector(object$residuals)/sc) :
##   length(resid.sc) == nobs(object) is not TRUE



## Fix these : ===================
str( s.do <- summary(robNls[["doCov"]]$na.exclude) )

try(## the "doCov" -- fails "only" when printing:
    print(s.do) #-> error
) -> .tmp
## Call:
## ....
## Residuals:
## Error in if (rdf > 5L) { : missing value where TRUE/FALSE needed

vcov(robNls[["doCov"]]$na.exclude)
##      Asym xmid scal
## Asym   NA   NA   NA
## xmid   NA   NA   NA
## scal   NA   NA   NA



try( ## debug this one
rmNAo <- nlrob(formula(fm1), data = DNAnase,
               na.action = na.omit,
               trace = TRUE,
               start = list(Asym = 3, xmid = 0, scal = 1)) ### fails
## Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  :
##   arguments imply differing number of rows: 19, 16
)


## From: "Pascal A. Niklaus" <pascal.niklaus@ieu.uzh.ch>
## To: <maechler@stat.math.ethz.ch>
## Subject: nlrob problem
## Date: Tue, 20 Dec 2011 07:04:38 +0100

## For psi-functions that can become zero (e.g. psi.bisquare), weights in
## the internal call to nls can become zero.


## Was
## psiTuk  <- robustbase:::psi.bisquare
## psiHamp <- robustbase:::psi.hampel

lmrob.control(psi="bisquare")$tuning.psi
psiTuk <- function(x, der=0) {
    ## cc:  dput( lmrob.control(psi="bisquare")$tuning.psi )
    if(der == 0)
        Mwgt(x, cc=4.685061, psi="Tukey")
    else
        Mpsi(x, cc=4.685061, psi="Tukey", deriv=1)
}

c.Ha <- lmrob.control(psi="hampel"); c.Ha$tuning.psi
psiHamp <- function(x, der=0) {
    ## cc: dput( lmrob.control(psi="hampel")$tuning.psi )
    if(der == 0)
	Mwgt(x, cc=c(1.35241275, 3.15562975, 7.212868), psi="Hampel")
    else
	Mpsi(x, cc=c(1.35241275, 3.15562975, 7.212868), psi="Hampel", deriv=1)
}

d <- data.frame(x = -6:9,
                y = 43 + c(7, 52, 21, 12, 10, -4, -5, -4, 0, -77, -8, -10, 22, 33, 38, 51))
nlr1 <- nlrob(y ~ a*(x + b*exp(-c*x)), start=list(a= 4, b= 1, c= 1.2),
              data = d,
              maxit = 50, # default 20 is *not* sufficient
              model = TRUE,
              trace=TRUE)
## These failed in robustbase version 0.8-0 and earlier
nlr2  <- update(nlr1, psi = psiTuk) # now *does* converge...
## check 'model' and dataClasses
stopifnot(is.list(mod <- nlr2$model), is.data.frame(mod),
          inherits(attr(mod, "terms"), "terms"),
          identical(dCl <- attr(attr(mod, "terms"),"dataClasses"),
                    nlr2$dataClasses),
          identical(dCl, c(y = "numeric", x = "numeric")))
## 'port' ditto:
nlr2. <- update(nlr2, algorithm= "port")
nlr3  <- update(nlr1, psi = psiHamp)   # *does* converge, too...
nlr3. <- update(nlr3, algorithm= "port")
summary(nlr2.)
summary(nlr3.)
i. <- -c(2, 15) # <- drop 'call' and 'iter' components
stopifnot(all.equal(nlr2[i.], nlr2.[i.], tolerance = 2e-5, check.environment = FALSE),
          all.equal(nlr3[i.], nlr3.[i.], tolerance = 1e-4, check.environment = FALSE),
          ## The redescending psi() give some exact 0 weights :
	  identical(which(abs(nlr2$rweights) < 1e-9), c(1L, 10 :12)),
	  identical(which(abs(nlr3$rweights) < 1e-9), c(1L, 10L,12L))
	  )

## Different example with more data:
pp <- list(a=10, b=4, c=1/4)
x <- seq(-6,9, by = 1/8)
f.x <- with(pp, a*(x + b*exp(-c*x)))
set.seed(6); y <- y0 <- f.x + 4*rnorm(x)
iO <- c(2:3,20,70:71,90); y[iO] <- y[iO] + 32*c(-1,-1,1)*(2+rlnorm(iO)); y <- round(y)
plot(x,y); lines(x, f.x, col="tomato", lty = 2)
dd <- data.frame(x,y)

nlc1 <- nls(formula(nlr1), start = coef(nlr1), data=dd, trace=TRUE)
nlR1 <- update(nlr1, data = dd)# update the model with the new data
summary(nlR1)
lines(x, predict(nlc1), col=3)
lines(x, predict(nlR1), col=4)
legend("top", c("f(x)", "least squares", "robust"),
       col=c("tomato", palette()[3:4]), lty=c(2,1,1))

## These both now *do* converge, but failed earlier
(nlbi <- update(nlR1, psi = psiTuk))
(nlFH <- update(nlR1, psi = psiHamp))
lines(x, predict(nlbi), col=5)
lines(x, predict(nlFH), col=6)

stopifnot(nlR1$status == "converged", nlbi$status == "converged",
	  nlFH$status == "converged")
assert.EQ(coef(nlR1), c(a=9.914874,    b=3.98612416,  c=0.250896252),  tol = 4e-9)
assert.EQ(coef(nlbi), c(a=9.947458207, b=3.954210623, c=0.2535835248), tol = 4e-9)
## This is suddently quite different :  ???!?!??
## assert.EQ(coef(nlFH), c(a=9.94242831, b=3.97370746, c=0.252907618))
assert.EQ(coef(nlFH),    c(a=9.952893755,b=3.949047387,c=0.2536216541), tol = 1e-7)
assert.EQ(1000*diag(vcov(nlR1)),
          c(a=16.167493, b=10.0986644, c=0.0200814189), tol = 7e-7, giveRE=TRUE)
assert.EQ(1000*local({V <- vcov(nlFH); V[lower.tri(V, diag=TRUE)]}),
          c(16.33774615, -9.704702857, 0.3149189329,
            10.03560556, -0.4079936961, 0.02039106329), tol = 7e-7)
assert.EQ(predict(nlR1), predict(nlbi), tol = 0.05, giveRE=TRUE)
assert.EQ(predict(nlR1), predict(nlFH), tol = 0.05, giveRE=TRUE)

nlFH2 <- update(nlFH, psi = .Mwgt.psi1("Hampel", c(2,4,8)))
## TODO: and compare
## TODO: same with Tukey


##----- *Vector* parameters indexed by factor  levels -------------
##----- MM: ~/R/MM/Pkg-ex/robustbase/nlrob-vectorpar.R

data(biomassTill)## see also smaller example in ../man/biomassTill.Rd

if(!dev.interactive(orNone=TRUE))  pdf("nlrob-biomT.pdf")

require(lattice)
xyplot(Biomass ~ DVS | Tillage, data = biomassTill)
xyplot(Biom.2  ~ DVS | Tillage, data = biomassTill)

## starting values
m0.st <- list(Wm = rep(200, 3),
               a = rep(  1, 3),
               b = rep(  2, 3))
##-> nls(), even with -expm1(.) fails to start properly (and hence nlrob() fails too):
try(
m.c0 <- nls(Biomass ~ Wm[Tillage] * (- expm1(-(DVS/a[Tillage])^b[Tillage])),
            data = biomassTill, start = m0.st, trace=TRUE)
)

## several other versions of the above fail similarly.  This works:
m00st <- list(Wm = rep(300,  3),
               a = rep( 1.5, 3),
               b = rep( 2.2, 3))
m.c00 <- nls(Biomass ~ Wm[Tillage] * (- expm1(-(DVS/a[Tillage])^b[Tillage])),
            data = biomassTill, start = m00st, trace=TRUE)

## These were the "true" beta for simulating in creation of {Biomass, Biom.2}:
m1.st <- list(Wm = c(219.8, 265.9, 343.4),
               a = c(1.461, 1.493, 1.294),
               b = c(2.889, 2.838, 4.046))

m.cl <- nls(Biomass ~ Wm[Tillage] * (1 - exp(-(DVS/a[Tillage])^b[Tillage])),
            data = biomassTill, start = m00st, trace=TRUE)
## this now fails to converge:
try( # "singular gradient"
m.c2 <- nls(Biom.2 ~ Wm[Tillage] * (1 - exp(-(DVS/a[Tillage])^b[Tillage])),
            data = biomassTill, start = m00st, trace=TRUE)
)

str(C1 <- nls.control(minFactor=1e-6, warnOnly=TRUE, printEval=TRUE, maxiter=500))

try(
m.c2 <- nls(Biom.2 ~ Wm[Tillage] * (1 - exp(-(DVS/a[Tillage])^b[Tillage])),
                data = biomassTill, start = m00st, trace=TRUE, control=C1)
)
## fails (!) too {numericDeriv() in iteration 129} even though we have
## 'warnOnly' ! ==> bug in nls()  !!!!!!!!!!!!!!!!!!!!!!!!!!!

## -expm1(u) is better than (1 - exp(u)) :
m.c00 <- nls(Biom.2 ~ Wm[Tillage] * (- expm1(-(DVS/a[Tillage])^b[Tillage])),
             data = biomassTill, start = m00st, trace=TRUE, control=C1)
## "fails" but returns .. very bad..
m.c00
## Use better starting values, as we have such problems:
m.c2 <- nls(Biom.2 ~ Wm[Tillage] * (- expm1(-(DVS/a[Tillage])^b[Tillage])),
            data = biomassTill, start = m1.st, trace=TRUE, control=C1)
## "fails" but returns at least: singular gradient iteration 126
m.c2

## Robust: not converging in 20 steps (only warning)
mrob <- nlrob(Biomass ~ Wm[Tillage] * (-expm1(-(DVS/a[Tillage])^b[Tillage])),
              data = biomassTill, start = m00st, trace=TRUE)
stopifnot(identical(mrob$dataClasses,
                    c(Biomass = "numeric", Tillage = "factor", DVS = "numeric")))
try(## now: singular gradient in nls
mr.2 <- nlrob(Biom.2  ~ Wm[Tillage] * (-expm1(-(DVS/a[Tillage])^b[Tillage])),
              data = biomassTill, start = m00st, trace=TRUE)
)

## Compare coeffs:
rbind(c.true = unlist(m1.st),
      cl0 = coef(m.c00),
      cl = coef(m.cl), rob = coef(mrob),
      c2 = coef(m.c2))#, r.2 = coef(mr.2))
## Compare  fit

## Now for plotting --- nice would be xyplot, but I don't easily see how:


(yl2  <- range(biomassTill[,"Biom.2"]))
(ylim <- range(biomassTill[,"Biomass"]))# --> *not* showing the two outliers!
## or even a bit more robustly:
## sfsmisc::rrange(biomassTill[,"Biom.2"]) ##-> -201.3064  394.0914

## using global data + fits from above
p.biomass.fits <- function(ylim = c(-200, 400), n = 257, f.DVS = 0.1,
                           leg.txt =
                               c(outer(c("nls()   ", "nlrob()"),
                                       c("", "[ + 2 outl.]"), paste)),
                           col = c("blue2","blue3","tomato","red3"),
                           lty = c(2,1,2,1),
                           lwd = 2)
{
    ## more and equispaced DVS values for nice plot:
    rr <- extendrange(biomassTill[,"DVS"], f=f.DVS)
    bbDVS <- seq(rr[1], rr[2], length = n)
    b.Till <- biomassTill[,"Tillage"]
    nP <- nlevels(b.Till) # == 3
    m <- length(leg.txt)
    col <- rep_len(col, m)
    lwd <- rep_len(lwd, m)
    lty <- rep_len(lty, m)
    ## Prefer xyplot() - this is ugly but works (and tests predict(*, <subset>)):
    op <- par(mfrow = c(nP,1), mar = .1 + c(3, 3, 2, 1), mgp = c(1.25, 0.6, 0))
    on.exit(par(op))
    for(lev in levels(b.Till)) {
        cat(lev,":\n--------\n")
        dsub <- subset(biomassTill, Tillage == lev)
        plot(Biom.2 ~ DVS, data = dsub, ylim=ylim, main = paste("Tillage = ", lev))
        grid()
        dd <- data.frame(Tillage = factor(rep.int(lev, n), levels=levels(b.Till)),
                         DVS = bbDVS)
        lines(predict(m.cl, dd) ~ DVS, data=dd, col=col[1], lty=lty[1], lwd=lwd[1])
        lines(predict(mrob, dd) ~ DVS, data=dd, col=col[2], lty=lty[2], lwd=lwd[2])
        lines(predict(m.c2, dd) ~ DVS, data=dd, col=col[3], lty=lty[3], lwd=lwd[3])
        ## lines(predict(mr.2, dd) ~ DVS, data=dd, col=col[4], lty=lty[4], lwd=lwd[4])
        if(lev == "CA-")
            legend("top", leg.txt, col = col, lty=lty, lwd=lwd,
                   inset=.02, bg = "gray96") #, bty="n")
    }
}

## showing all data points:
p.biomass.fits(ylim = yl2)
## more interesting:
p.biomass.fits()


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
