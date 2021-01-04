### Have had cases where differences between large numbers lose precision, or even give Inf,
### which lead to NA

require(robustbase)
M <- .Machine$double.xmax


stopifnot(exprs = {
    all.equal(scaleTau2(c(-4:4, 10000), consistency=FALSE),
             (scaleTau2(c(-4:4, 1e300), consistency=FALSE) -> sT), # <- gave NaN, now fine !
             tol = 1e-15) # even 0 (exact equality; Linux 64b)
    all.equal(3.41103800034854, sT, tol = 1e-14) # seen 6.5e-16
})

exL <- list( xNA  = c(NA, 1:6)
, xM   = c(-4:4, M)
, xI   = c(-4:4, Inf)
, IxI  = c(-Inf, -4:4, Inf)
, IxI2 = c(-Inf, -4:4, Inf,Inf)
  ##
, MxM1 = c(rep(-M, 3), 1:9,  rep(M, 7)) # < 50% "good"
, MxMe = c(rep(-M, 3), 1:10, rep(M, 7)) # half  "good"
, MxM2 = c(rep(-M, 3), 1:11, rep(M, 7)) # > 50% "good"
)

madL <- vapply(exL, mad, pi)
## Initially, scaleTau2() "works" but gives  NaN  everywhere !!
sT2.L <- vapply(exL, scaleTau2, FUN.VALUE=1, sigma0 = 1, consistency=FALSE)

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

## FIXME: should not give NaN :
scaleTau2(xI)

## FIXME: even give  Error in ..... : NA/NaN/Inf in foreign function call (arg 1)
try( Sn(xI) )
try( Qn(xI) )

## From  example(mc)  {by MM} :

## Susceptibility of the current algorithm to large outliers :
dX10 <- function(X) c(1:5,7,10,15,25, X) # generate skewed size-10 with 'X'
(Xs <- c(10,20,30, 60, 10^(2:10), 1000^(4:19), 1e6^c(10:20,10*(3:5)), Inf))

(mc10x <- vapply(Xs, function(x) mc(dX10(x)), 1))
plot(Xs, mc10x, type="b", main = "mc( c(1:5,7,10,15,25, X) )", xlab="X", log="x")

##--- FIXME: the above must change!

## so, Inf does work, indeed for mc()
dX10(Inf)
set.seed(2020-12-04)
stopifnot(exprs = {
    is.finite(mc(dX10(Inf))) # 0.5 currently
    mc(c(-Inf, rlnorm(100), Inf)) == 0
})

