library(robustbase)

#### Poisson examples from Eva Cantoni's paper

### Using Possum Data
### ================

data(possumDiv)

## Try to follow closely Cantoni & Ronchetti(2001), JASA
X. <- possum.mat[, -1]
y  <- possum.mat[, "Diversity"]

## "classical via robust: c = Inf :
Inf. <- 1e5 ## --- FIXME

### FIXME:  The following fails because glm.fit() returns NA coefficients
### -----   We should eliminate corresponding columns, just
###         as Martin's  glm.rob(X.,y, choice= "poisson") from 'robGLM1' did!

glm.cr <- glmrob(y ~ X., family = "poisson", tcc = Inf.)
(scr <- summary(glm.cr))

## c = 2.0
g2 <- glmrob(y ~ X., family = "poisson", tcc = 2.0)
summary(g2)

## c = 1.6
glm.r <- glmrob(y ~ X., family = "poisson", tcc = 1.6)
coef(summary(glm.r))

str(glm.r)

###--- reduce the matrix from singularity ourselves:
X. <- possum.mat[, -c(1, match(c("E.nitens", "NW-NE"), colnames(possum.mat)))]
dim(X.)# 151 11

glm.c2 <- glmrob(y ~ X., family = "poisson", tcc = Inf.)
summary(glm.c2)

## c = 1.6
glm.r2 <- glmrob(y ~ X., family = "poisson", tcc = 1.6)
summary(glm.r2)


###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
