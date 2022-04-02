### Have had cases where differences between large numbers lose precision, or even give Inf,
### which lead to NA

require(robustbase)

source(system.file("xtraR/styleData.R", package = "robustbase"))# -> smallD, mkMx()

stopifnot(exprs = {
    all.equal(scaleTau2(c(-4:4, 10000), consistency=FALSE),
             (scaleTau2(c(-4:4, 1e300), consistency=FALSE) -> sT), # <- gave NaN, now fine !
             tol = 1e-15) # even 0 (exact equality; Linux 64b)
    all.equal(3.41103800034854, sT, tol = 1e-14) # seen 6.5e-16
})

exL <- c(
    list( xNA  = c(NA, 1:6)
       , xMe9 = c(-4:4, 1e9)
       , xM   = c(-4:4, .Machine$double.xmax)
       , xI   = c(-4:4, Inf)
       , IxI  = c(-Inf, -4:4, Inf)
       , IxI2 = c(-Inf, -4:4, Inf,Inf))
    ,
    mkMx(c(1e6, 1e9,
           1e12, 1e14, 1e16,
           1e20, 1e40, .Machine$double.xmax, Inf))
)

madL <- vapply(exL, mad, pi)
## Initially, scaleTau2() "works" but gives  NaN  everywhere -- now fine!
sT2.L    <- vapply(exL, scaleTau2, FUN.VALUE=1, sigma0 = 1, consistency=FALSE)
sT2.i2.L <- vapply(exL, scaleTau2, FUN.VALUE=1, sigma0 = 1, consistency=FALSE, iter = 2)
sT2.i5.L <- vapply(exL, scaleTau2, FUN.VALUE=1, sigma0 = 1, consistency=FALSE, iter = 5)

cbind(madL, sT2.L)

stopifnot(exprs = {
    is.na(madL [1])
    is.na(sT2.L[1])
    2.3 < sT2.L[-1]
          sT2.L[-1] < 2.71
})


xI <- exL$xI
stopifnot(exprs = {
    mad(exL$xI, constant = 1) == 2.5
})

stopifnot(exprs = {
    all.equal(3.5471741782, scaleTau2(xI)) # gave NaN
    ## These even gave  Error in ..... : NA/NaN/Inf in foreign function call (arg 1)
    all.equal(3.5778,       Sn(xI))
    all.equal(3.1961829592, Qn(xI))
})

## From  example(mc)  {by MM} :

## Susceptibility of the current algorithm to large outliers :
dX10 <- function(X) c(1:5,7,10,15,25, X) # generate skewed size-10 with 'X'
(Xs <- c(10,20,30, 60, 10^(2:10), 1000^(4:19), 1e6^c(10:20,10*(3:5)), Inf))

mc10x <- vapply(Xs, function(x) mc(dX10(x)), 1)
## now fixed:
stopifnot(all.equal(c(4,6, rep(7,42))/12, mc10x))
plot(Xs, mc10x, type="b", main = "mc( c(1:5,7,10,15,25, X) )", xlab="X", log="x")

## so, Inf does work, indeed for mc()
mcOld <- function(x, ..., doScale=TRUE) mc(x, doScale=doScale, c.huberize=Inf, ...)
(x10I <- dX10(Inf))
set.seed(2020-12-04)# rlnorm(.)
summary(xlN <- rlnorm(100))
xII <- c(-Inf, xlN, Inf)
stopifnot(exprs = {
    all.equal(0.5,  print(mcOld(x10I)))
    all.equal(7/12, print(mc   (x10I, doScale=TRUE ))) # was 0.5 before huberization
    all.equal(7/12, print(mc   (x10I, doScale=FALSE)))
    mcOld(xII) == 0
    all.equal(0.3646680319, mc(xII))
})
