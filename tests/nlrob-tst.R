library(robustbase)

source(system.file("xtraR/ex-funs.R", package = "robustbase"))
## -> assert.EQ(), identical3(), ..

DNase1 <- DNase[ DNase$Run == 1, ]
Y <- DNase1[,"density"] # for convenience below

## classical
fm1 <- nls(density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
	   data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), trace=TRUE)
summary(fm1)

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
assert.EQ(coef(rm1), #  print(coef(rm1), digits=12)
	  ##	 2.35963008460	  1.49945088410	   1.04506391722  F19 Lx 64b
	  ##	 2.35963008613	  1.49945088600	   1.04506391793  F19 Lx 32b
	  c(Asym=2.35963008, xmid=1.49945088, scal=1.04506392), tol = 1e-8)
assert.EQ(sqrt(diag(sm1$cov)),
	  ## 32b 0.08626872273,     0.0902194541,      0.03503833759
	  c(Asym=0.0862687305, xmid=0.0902194608, scal=0.0350383389),
	  tol = 7e-7)

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
              trace=TRUE)
## These failed in robustbase version 0.8-0 and earlier
nlr2  <- update(nlr1, psi = psiTuk) # now *does* converge...
## 'port' ditto
nlr2. <- update(nlr2, algorithm= "port")
nlr3  <- update(nlr1, psi = psiHamp)   # *does* converge, too...
nlr3. <- update(nlr3, algorithm= "port")
summary(nlr2.)
summary(nlr3.)
i. <- -c(2, 15) # <- drop 'call' and 'iter' components
stopifnot(all.equal(nlr2[i.], nlr2.[i.], tolerance = 2e-5),
          all.equal(nlr3[i.], nlr3.[i.], tolerance = 1e-4),
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
nlR1 <- update(nlr1, data = dd)
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
assert.EQ(coef(nlR1), c(a=9.914874,    b=3.98612416,  c=0.250896252),  tol = 1e-9)
assert.EQ(coef(nlbi), c(a=9.947458207, b=3.954210623, c=0.2535835248), tol = 1e-9)
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


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
