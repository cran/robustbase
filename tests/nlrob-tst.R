library(robustbase)

identical3 <- function(x,y,z)	  identical(x,y) && identical (y,z)
identical4 <- function(a,b,c,d)   identical(a,b) && identical3(b,c,d)

DNase1 <- DNase[ DNase$Run == 1, ]
Y <- DNase1[,"density"] # for convenience below

## classical
fm1 <- nls(density ~ Asym/(1 + exp(( xmid - log(conc) )/scal ) ),
	   data = DNase1, start = list(Asym = 3, xmid = 0, scal = 1), trace=TRUE)
summary(fm1)

## robust
rm1 <- nlrob(formula(fm1), data = DNase1, trace = TRUE,
	     start = list(Asym = 3, xmid = 0, scal = 1))
summary(rm1)
stopifnot(all.equal(Y, fitted(fm1) + residuals(fm1), check.attr=FALSE),
	  ## fitted(<nls>) has "label" attribute
	  identical3(c(fitted(fm1)), predict(fm1), predict(fm1, newdata=DNase1)),
	  ## robust fit :
	  identical3(fitted(rm1), predict(rm1), predict(rm1, newdata=DNase1)),
	  all.equal(Y, unname(fitted(rm1) + residuals(rm1)))
	  )

## From: "Pascal A. Niklaus" <pascal.niklaus@ieu.uzh.ch>
## To: <maechler@stat.math.ethz.ch>
## Subject: nlrob problem
## Date: Tue, 20 Dec 2011 07:04:38 +0100

## For psi-functions that can become zero (e.g. psi.bisquare), weights in
## the internal call to nls can become zero.

d <- data.frame(x = -6:9,
                y = 43 + c(7, 52, 21, 12, 10, -4, -5, -4, 0, -77, -8, -10, 22, 33, 38, 51))
nlr1 <- nlrob(y ~ a*(x + b*exp(-c*x)), start=list(a= 4, b= 1, c= 1.2),
              data = d,
              maxit = 50, # default 20 is *not* sufficient
              trace=TRUE)
## These failed in robustbase version 0.8-0 and earlier
nlr2  <- update(nlr1, psi = robustbase:::psi.bisquare) # now *does* converge...
## 'port' ditto
nlr2. <- update(nlr2, algorithm= "port")
nlr3  <- update(nlr1, psi = robustbase:::psi.hampel)   # *does* converge, too...
nlr3. <- update(nlr3, algorithm= "port")
summary(nlr2.)
summary(nlr3.)
i. <- -c(2, 15) # <- drop 'call' and 'iter' components
stopifnot(all.equal(nlr2[i.], nlr2.[i.], tol = 2e-5),
          all.equal(nlr3[i.], nlr3.[i.], tol = 1e-4),
          ## The redescending psi() give some exact 0 weights :
	  identical(which(abs(nlr2$rweights) < 1e-9), c(1L, 10 :12)),
	  identical(which(abs(nlr3$rweights) < 1e-9), c(1L, 10L,12L))
	  )

psiFHa <- robustbase:::psi.hampel

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

## These both do *not* converge currently
update(nlR1, psi = robustbase:::psi.bisquare)# failed in earlier versions
update(nlR1, psi = robustbase:::psi.hampel)  # failed in earlier versions
