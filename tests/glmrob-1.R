library(robustbase)

source(system.file("xtraR/ex-funs.R", package = "robustbase"))
## -> assert.EQ()


###>> 1 ------------------- family = poisson ------------------------------------

### very simple model [with outliers]
set.seed(113)
y <- rpois(17, lambda = 4) ## -> target:  beta_0 = log(E[Y]) = log(4) = 1.386294

y[1:2] <- 99:100 # outliers
y
rm1 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              acc = 1e-13) # default is just 1e-4
## and check the robustness weights:
assert.EQ(c(0.0287933850640724, 0.0284930623638766,
		      0.950239140568007, 0.874115394740014),
		    local({w <- rm1$w.r; w[ w != 1 ] }), tol = 1e-14)
assert.EQ(coef(rm1), c("(Intercept)" = 1.41710946076738),tol = 1e-14)

cm1 <- glm   (y ~ 1, family = poisson, trace = TRUE)

rmMT <- glmrob(y ~ 1, family = poisson, trace = TRUE, method="MT")
(sMT <- summary(rmMT))

if(FALSE) # for manual digging:
debug(robustbase:::glmrobMqle)

allresid <- function(obj, types = c("deviance", "pearson", "working", "response"))
{
    sapply(types, residuals, object = obj)
}

okFit <- function(obj, check.attr=FALSE, ...) {
  all.equal(obj$y, obj$fitted.values + residuals(obj, "response"),
            check.attr=check.attr, ...)
}

## check validity of several methods simultaneously:
y. <- model.response(model.frame(rm1))
stopifnot(okFit(cm1), okFit(rm1), y. == y)

alr.c <- allresid(cm1)
alr.r <- allresid(rm1)

## MM --- just for now --
plot(resid(cm1),                resid(rm1), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="pearson"), resid(rm1, type="pearson"), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="working"), resid(rm1, type="working"), asp=1); abline(0,1, col=2)

## leave away the outliers --
cm0 <- glm   (y ~ 1, family = poisson, trace = TRUE, subset = -(1:2))
plot(resid(cm0),                resid(rm1)[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="pearson"), resid(rm1, type="pearson")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="working"), resid(rm1, type="working")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="response"), resid(rm1, type="response")[-(1:2)], asp=1); abline(0,1, col=2)


## Use weights (slow convergence !)
w2 <- c(rep(1,8), rep(10,9))
rm2 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              weights = w2, maxit = 500, acc = 1e-10) # default is just 1e-4
## slow convergence
stopifnot(okFit(rm2))


###>> 2 ------------------- family = binomial -----------------------------------

## Using  *factor*  y ...
x <- seq(0,5, length = 120)
summary(px <- plogis(-5 + 2*x))
set.seed(7)
(f <- factor(rbinom(length(x), 1, prob=px)))

summary(m.c0 <- glm   (f ~ x, family = binomial))
summary(m.r0 <- glmrob(f ~ x, family = binomial))

## add outliers --- in y:
f. <- f
f.[i1 <- 2:3] <- 1
f.[i0 <- 110+c(1,7)] <- 0
        m.c1 <- glm   (f. ~ x, family = binomial)
summary(m.r1 <- glmrob(f. ~ x, family = binomial)) ## hmm, not so robust?
stopifnot(m.r1$w.r[c(i0,i1)] < 1/3, # well, at least down weighted
	  ## and coefficients change less :
	  (coef(m.r1) - coef(m.c0)) / (coef(m.c1) - coef(m.c0)) < 1)
assert.EQ(c("(Intercept)" = -3.10817337603974, x = 1.31618564057790),
	  coef(m.r1), tol= 1e-14, giveRE=TRUE)

y <- as.numeric(as.character(f.))
m.r2  <- BYlogreg(x0=x, y=y, trace=TRUE, maxhalf= 10)
m.r2A <- BYlogreg(x0=x, y=y, trace= 2  , maxhalf= 15)
## different.. but not so much:
iB <- 1:5
assert.EQ(m.r2A[iB], m.r2[iB], tol = .003, giveRE=TRUE)


assert.EQ(c("Intercept" = -2.9554950286, x = 1.2574679132),
          ## 32-bit{ada-5}  -2.95549502890363   1.25746791332613
	  m.r2$coef, tol=4e-10, giveRE=TRUE)
assert.EQ( c(0.685919891749065, 0.256419206157062),
          ## 32-bit{ada-5}:
          ## 0.685919891858219, 0.256419206203016)
	  m.r2$sterror, tol=4e-10)

data(foodstamp)
str(foodstamp)
## Model with 'income' instead of log(income+1)  is "interesting"
## because BYlogreg() needs  maxhalf > 10 for convergence!
m.fs0   <- glm   (participation ~ ., family=binomial, data=foodstamp)
m.fs0QL <- glmrob(participation ~ ., family=binomial, data=foodstamp)
y.fs <- foodstamp[,"participation"]
X.fs0 <- model.matrix(m.fs0)
head(X.fs0)
## (former default) maxhalf = 10  leads to too early convergence:
m.fsWBY. <- BYlogreg(x0=X.fs0, y=y.fs,
                     addIntercept=FALSE, trace=TRUE, maxhalf=10)
m.fs.BY. <- BYlogreg(x0=X.fs0, y=y.fs, initwml=FALSE,
                     addIntercept=FALSE, trace=TRUE, maxhalf=10)
m.fsWBY <- BYlogreg(x0=X.fs0, y=y.fs,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
m.fs.BY <- BYlogreg(x0=X.fs0, y=y.fs, initwml=FALSE,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)

assert.EQ(m.fsWBY.[iB], m.fsWBY[iB], tol= 0.07)## almost 7% different
assert.EQ(m.fs.BY.[iB], m.fs.BY[iB], tol= 0.08)

foodSt <- within(foodstamp, { logInc <- log(1 + income) ; rm(income) })

m.fsML <- glm   (participation ~ ., family=binomial, data=foodSt)
m.fsQL <- glmrob(participation ~ ., family=binomial, data=foodSt)
X.fs <- model.matrix(m.fsML)
stopifnot(dim(X.fs) == c(150, 4)) # including intercept!
try(## FIXME -- Mahalanobis fails with singular matrix, here:
m.fsWBY <- BYlogreg(x0=X.fs, y=y.fs,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
)
## maxhalf=18 is too much --> no convergence
m.fs.BY <- BYlogreg(x0=X.fs, y=y.fs, initwml=FALSE,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
signif(
    rbind(ML = coef(m.fsML),   QL =coef(m.fsQL),
          WBY0=coef(m.fsWBY.), BY0=coef(m.fs.BY.),
          WBY =coef(m.fsWBY ), BY =coef(m.fs.BY)
          )
    , 4)


if(FALSE) {
## *scaling* of X (  ?? <==> ??   'sigma1' ) ------------------

## no "W" (Mahalanobis fail because of *singular* X):
m.fs.BY100 <- BYlogreg(x0=100*X.fs, initwml=FALSE,
                       y=y.fs,
                       addIntercept=FALSE, trace=TRUE, maxhalf=18)


X1c <- cbind(1, 100*X.fs[,-1])
m.fsWBY1c <- BYlogreg(x0=X1c, y=y.fs,
                      addIntercept=FALSE, trace=TRUE, maxhalf=18)

}## not yet
