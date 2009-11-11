library(robustbase)

### very simple model [with outliers]
set.seed(113)
y <- rpois(17, lambda = 4)
y[1:2] <- 99:100 # outliers
y
options(digits=10)
rm1 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              acc = 1e-13) # default is just 1e-4
stopifnot(all.equal(c(0.0287933850640724, 0.0284930623638766,
		      0.950239140568007, 0.874115394740014),
		    local({w <- rm1$w.r; w[ w != 1 ] }), tol = 1e-14),
	  all.equal(coef(rm1), c("(Intercept)" = 1.41710946076738),
		    tol = 1e-14)
	  )
if(FALSE) # for manual digging:
debug(robustbase:::glmrobMqle)

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
	  (coef(m.r1) - coef(m.c0)) / (coef(m.c1) - coef(m.c0)) < 1,
	  all.equal(coef(m.r1),
		    c("(Intercept)" = -3.10817337603974,
		      x = 1.31618564057790), tol= 1e-14)
          )

