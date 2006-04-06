
library(robustbase)

## EX 1
data(coleman)
## "empty model" (not really a lot of sense)
(m0 <- lmrob(Y ~ 0, data = coleman))
summary(m0)
## "Intercept" only: robust mean
(m1 <- lmrob(Y ~ 1, data = coleman))
summary(m1)


(mC <- lmrob(Y ~ ., data = coleman))
summary(mC)
## Values will change once we use R's random number generator !
stopifnot(
all.equal(unname(coef(mC)), c(29.51, -1.66, 0.0834, 0.666, 1.18, -4.01),
          tol = 1e-3) # tol =0.000146 for 64-bit, 0.000251 for 32-bit
)
str(mC)

## EX 2
gen <- function(n,p, n0, y0, x0, beta = rep(1, p))
{
    stopifnot(n >= 1, p >= 1, n0 >= 0, length(beta) == p)
    x <- matrix(rnorm(n*p),n,p) # iid x's
    y <- x %*% beta + rnorm(n)
    xc <- matrix(0,n0,p)
    xc[,1] <- x0
    xc <- xc + 0.1*matrix(rnorm(n0*p),n0,p)
    x[1:n0,] <- xc
    y[1:n0] <- y0 + .001*rnorm(n0)
    list(x=x, y=y)
}

## generate --a sample of  n  observations with  p  variables
## and 10% of outliers near (x1,y) = (10,10)
n <- 500 ; n0 <- n %/% 10
p <- 7 ## p = 20 is more impressive but too slow for "standard test"
set.seed(17)
a <- gen(n=n, p=p, n0= n0, y0=10, x0=10)
plot(a$x[,1], a$y, col = c(rep(2, n0), rep(1, n-n0)))
system.time( m1 <- lmrob(y~x, data = a) )
plot(m1) #- currently 5 plots; MM:I  don't like #3 (Response vs fitted)

## don't compute robust distances --> faster by factor of two:
system.time(m2 <- lmrob(y~x, data = a,
                        control = lmrob.control(compute.rd = FALSE)))
## ==> half of the CPU time is spent in covMcd()!
(sm2 <- summary(m2))
l1 <- lm(y~x, data = a)
cbind(robust = coef(sm2)[,1:2],
      lm = coef(summary(l1))[,1:2])

##--- Now use n > 2000 --> so we use C internal fast_s_large_n(...)
n <- 2500 ; n0 <- n %/% 10
a2 <- gen(n=n, p = 3, n0= n0, y0=10, x0=10)
plot(a2$x[,1], a2$y, col = c(rep(2, n0), rep(1, n-n0)))
system.time( m3 <- lmrob(y~x, data = a2) )
m3
system.time( m4 <- lmrob(y~x, data = a2, compute.rd = FALSE))
(sm4 <- summary(m4))

stopifnot(identical(coef(m3), coef(m4)),
          all.equal(unname(coef(m3)),
                    c(0.03802, 0.99653, 1.00555, 0.99981), tol= 4e-5),
          all.equal(unname(coef(sm4)[,"Std. Error"]),
                    c(0.0252, 0.00287, 0.0260, 0.0263), tol = 2e-3)
          )



## rm(a,m1, m2, m3, m4, sm2, l1)



cat('Time elapsed: ', proc.time(),'\n') # "stats"
