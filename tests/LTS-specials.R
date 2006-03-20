#### Test special cases  for  ltsReg()

library(robustbase)

### 1) p = 1 ----------------------------------------------------
set.seed(1)
x <- c(rnorm(50),100, 1e10)
(r1 <- ltsReg(x ~ 1)) # failed in Valentin's 1.0-3 (pre-version)
summary(r1)
(r1. <- ltsReg(y = x))
i1 <- 14:16; ii <- (1:20)[-i1]
UN <- function(lis) lapply(lis, unname)
stopifnot(all.equal(   r1[ii],     r1.[ii],  tol= 1e-15),
          all.equal(UN(r1[i1]), UN(r1.[i1]), tol= 1e-15))

## with alpha = 1
(r1.1 <- ltsReg(x ~ 1, alpha = 1))
summary(r1.1)

### 1b) p = 1, constant scale
(rc <- ltsReg(y = rep(1,12)))
str(rc)
summary(rc)
## with alpha = 1
(rc1 <- ltsReg(y = rep(1,12), alpha = 1))
summary(rc1)
stopifnot(residuals(rc) == 0,  all.equal(unname(coef(rc )), 1),
          residuals(rc1) == 0, all.equal(unname(coef(rc1)), 1))

### 2)  alpha = 1 : classical estimates --- for general cases --------


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
