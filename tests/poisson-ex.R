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
cf2 <- matrix(c(-0.898213938628341, 0.269306882951903,
                0.00717220104127189, 0.0224349606070713,
                -0.25335520175528,  0.288588183720387,
                0.0403970350911325, 0.0113429514237665,
                0.0411096703375411, 0.0145996036305452,
                0.0730250489306713, 0.0386771060643486,
                0.0176994176433365, 0.0107414247342375,
                -0.0289935051669504,0.194215229266707,
                0.149521144883774,  0.271648514202971,
                0.0503262879663932, 0.191675979065398,
                0.0909870068741749, 0.192192515800464,
                -0.512247626309172, 0.250763990619973), 12,2, byrow=TRUE)
cfE <- unname(coef(sglm.r2)[, 1:2])
all.equal(cfE, cf2, tol=0)#-> show : ~ 1.46e-11
stopifnot(all.equal(cfE, cf2, tol = 1e-9),
          abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)




###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
