library(robustbase)

#### Poisson examples from Eva Cantoni's paper

### Using Possum Data
### ================

data(possumDiv)

## Try to follow closely Cantoni & Ronchetti(2001), JASA
X <- possum.mat[, -1]
y <- possum.mat[, "Diversity"]

## "classical via robust: c = Inf :
Inf. <- 1e5 ## --- FIXME

## The following used to fails because glm.fit() returns NA coefficients
## now fine
glm.cr <- glmrob(y ~ X, family = "poisson", tcc = Inf.)
(scr <- summary(glm.cr))

## c = 2.0
g2 <- glmrob(y ~ X, family = "poisson", tcc = 2.0, trace=TRUE)
summary(g2)

## c = 1.6
glm.r <- glmrob(y ~ X, family = "poisson", tcc = 1.6)
coef(summary(glm.r))

str(glm.r)

###--- reduce the matrix from singularity ourselves:
X. <- possum.mat[, -c(1, match(c("E.nitens", "NW-NE"), colnames(possum.mat)))]
dim(X.)# 151 11

glm.c2 <- glmrob(y ~ X., family = "poisson", tcc = Inf.)
summary(glm.c2)

## c = 1.6,  x-weights, as in Cantoni-Ronchetti
glm.r2 <- glmrob(y ~ X., family = "poisson",
                 tcc = 1.6, weights.on.x = "hat")

## Now the same, for the direct possum data (no matrix),
## This indeed gives the same coefficients as in
## Table 3 of Cantoni+Ronchetti(2001): .. (tech.rep.):
glm.r2. <- glmrob(Diversity ~ ., family = "poisson", data=possumDiv,
                  tcc = 1.6, weights.on.x = "hat", acc = 1e-15)
## here iterate till convergence (acc = 10^(-15))
(sglm.r2 <- summary(glm.r2.))
## This is too accurate for S.E. (but we have converged to end)
cf2 <- matrix(c(-0.8978077842408, 0.26827453582,
                0.009943172430292,0.022240630032,
                -0.2514151110520, 0.28762303617,
                0.04016769926460, 0.011318137566,
                0.03999044772925, 0.014510240865,
                0.07141397382767, 0.038516287553,
                0.01777760987761, 0.010697122750,
                -0.02022762983567,0.19380237404,
                0.1269337686739,  0.27377258770,
                0.06009390334771, 0.19128320099,
                0.09491767881087, 0.19201471128,
                -0.5079313032077, 0.25052802424), 12,2, byrow=TRUE)
cfE <- unname(coef(sglm.r2)[, 1:2])
all.equal(cfE, cf2, tol=0)#-> show : ~ 1.46e-11
stopifnot(all.equal(cfE, cf2, tol = 1e-9),
          abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)




###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
