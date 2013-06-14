library(robustbase)

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...
source(system.file("xtraR/ex-funs.R", package = "robustbase"))
## -> newer assert.EQ()  {TODO: no longer needed in 2015}

#### Poisson examples from Eva Cantoni's paper

### Using Possum Data
### ================

data(possumDiv)

## Try to follow closely Cantoni & Ronchetti(2001), JASA
dim(X <- possum.mat[, -1]) # 151 13
str(y <- possum.mat[, "Diversity"])
##--- reduce the matrix from singularity ourselves:
X. <- possum.mat[, -c(1, match(c("E.nitens", "NW-NE"), colnames(possum.mat)))]
dim(X.)# 151 11

## "classical via robust: c = Inf :
Inf. <- 1e5 ## --- FIXME

## The following used to fail because glm.fit() returns NA coefficients
## now fine .. keep this as test!
glm.cr <- glmrob(y ~ X, family = "poisson", tcc = Inf.)
(scr <- summary(glm.cr))

scl <- summary(glm.cl <- glm   (Diversity ~ . , data=possumDiv, family=poisson))
sc2 <- summary(glm.c2 <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = Inf.))
assert.EQ(coef(scl), coef(sc2), tol = 6e-6, giveRE=TRUE) # 1.369e-6

## c = 2.0
summary(g2 <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = 2.0, trace=TRUE))

## c = 1.6
glm.r <- glmrob(Diversity ~ . , data=possumDiv, family=poisson, tcc = 1.6, trace=TRUE)
(s.16 <- summary(glm.r))
str(glm.r)

## Now with *smaller* X (two variablesless):
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
assert.EQ(cfE, cf2, tol = 1e-9, giveRE=TRUE)#-> show : ~ 1.46e-11
stopifnot(abs(glm.r2.$iter - 18) <= 1) # 18 iterations on 32-bit (2008)

## MT estimator -- "quick" examples

if(!robustbase:::doExtras()) {
    cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
    quit()
}
## if ( doExtras ) -----------------------------------------------------

X1 <- cbind(1, X.)

if(FALSE) ## for debugging ...
options(warn = 1, error=recover)

set.seed(57)
showSys.time(
    ## m1 <- glmrobMT(x=X1, y=y)
    m1 <- glmrob(Diversity ~ ., data=possumDiv, family=poisson, method="MT")
)
writeLines(tail(capture.output(warnings())))

stopifnot(m1$converged)
assert.EQ(m1$initial,
c(-0.851594294907422, -0.0107066895370536, -0.226958540075445, 0.0355906625338308,
  0.048010654640958, 0.0847493155436896, 0.0133604488401352, -0.024115201062159,
  0.0270535337324518, 0.146022135657894, -0.00751380783260833, -0.417638086169033)
          , tol = 1e-13, check.attr=FALSE, giveRE=TRUE)

## MM: I'm shocked how much this changes after every tweak ...

(arch <- Sys.info()[["machine"]])

dput(signif(unname(coef(m1)), 11)) ## -->
beta1 <- list(i686 =
c(-0.83703411949, 0.0085332123121, -0.16723425506, 0.040966496642,
  0.042392008333, 0.063163791628, 0.018623119547, -0.0063568893968,
  0.1143288682, 0.091609212, -0.025109234019, -0.6689999431)
, "x86_64" =
c(-0.83723213945, 0.0085385261915, -0.16697112315, 0.040985126003,
  0.042400738973, 0.063168847366, 0.01863253681, -0.0064477807228,
  0.11488937188, 0.091283185006, -0.025627390293, -0.66995658693)
)
## just FYI :
assert.EQ(beta1[[1]], beta1[[2]], tol = 0.002, check.attr=FALSE, giveRE=TRUE)

assert.EQ(coef(m1), beta1[[arch]], tol = 1e-10, check.attr=FALSE, giveRE=TRUE)

## The same, with another seed:
set.seed(64)
showSys.time(
    ## m2 <- glmrobMT(x=X1, y=y)
    m2 <- glmrob(Diversity ~ ., data=possumDiv, family=poisson, method="MT")
)
writeLines(tail(capture.output(warnings())))

stopifnot(m2$converged)
if(FALSE)
dput(signif(unname(m2$initial), 13)) ## -->
assert.EQ(m2$initial, ## so this is *not* platform (32bit/64bit) dependent:
c(-1.204304813829, 0.02776038445201, -0.3680174045842, 0.04325746912892,
  0.03895315289169, 0.04537145479989, 0.02847987541025, 0.07073207523212,
  0.355491639539, 0.1822955449528, 0.1323720331562, -0.3419939094877)
          , tol = 1e-12, check.attr=FALSE, giveRE=TRUE)

beta2 <-
 c(-0.722466858080995, 0.00853214741170862, -0.167080756780928, 0.0409708489859932,
   0.0423892852206197, 0.0631575676559838,  0.0186247550404213, -0.114368960856547,
   -0.12071086045638,  0.0913163131743179, -0.0254084789601603, -0.669355337605894)

beta2 <-
 c(-0.719639958181682, 0.00850921638739351, -0.166172980617748, 0.0410112563931713,
   0.0423638710960024, 0.0630362342605086, 0.0186359589765251, -0.116645159395719,
 -0.123130115061652, 0.0910865027225399, -0.0256449044169698, -0.67024227284216)

dput(signif(unname(coef(m2)), 11)) ## -->
beta2 <- list(i686 =
c(-0.83731166205, 0.0085694733057, -0.16776715809, 0.040959111008,
  0.042395767814, 0.063185207067, 0.018635039106, -0.0062832385739,
  0.11406869258, 0.091362008749, -0.025362102704, -0.6688654014)
, "x86_64" =
c(-0.83687097624, 0.0085341676033, -0.1674299545, 0.040968820903,
  0.042397459287, 0.063159075944, 0.018625582804, -0.0063140636571,
  0.11426134017, 0.091317308575, -0.025373078819, -0.66957444238)
)

## just FYI :
assert.EQ(beta2[[1]], beta2[[2]], tol = 0.001, check.attr=FALSE, giveRE=TRUE)

assert.EQ(coef(m2), beta2[[arch]], tol = 1e-10, check.attr=FALSE, giveRE=TRUE)
## slight changes of algorithm often change the above by ~ 4e-4 !!!

###---- Model Selection -----

## (not yet)  [ MM had this in ../../robGLM1/tests/quasi-possum.R ]

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
