### Have had cases where differences between large numbers lose precision, or even give Inf,
### which lead to NA

require(robustbase)

stopifnot(exprs = {
    all.equal(scaleTau2(c(-4:4, 10000), consistency=FALSE),
             (scaleTau2(c(-4:4, 1e300), consistency=FALSE) -> sT), # <- gave NaN, now fine !
             tol = 1e-15) # even 0 (exact equality; Linux 64b)
    all.equal(3.41103800034854, sT, tol = 1e-14) # seen 6.5e-16
})

mkMx <- function(M, ngood = 10, left = floor(ngood/3)) {
    stopifnot(is.numeric(ngood), ngood >= 3,
              is.numeric(M), length(M) == 1L, M >= 1000,
              is.numeric(left), 0 <= left, left <= ngood)
    right <- ngood-left
    res <- list(
        c(rep(-M, left), seq_len(ngood - 1L), rep(M, right)) # < 50% "good"
      , c(rep(-M, left), seq_len(ngood     ), rep(M, right)) # half  "good"
      , c(rep(-M, left), seq_len(ngood + 1L), rep(M, right)) # > 50% "good"
    )
    nM <- gsub("[-+]", "", formatC(M, digits=2, width=1))
    names(res) <- paste0("M", nM,"_n", c("m1", "eq", "p1"))
    res
}

exL <- c(
    list( xNA  = c(NA, 1:6)
       , xMe9 = c(-4:4, 1e9)
       , xM   = c(-4:4, .Machine$double.xmax)
       , xI   = c(-4:4, Inf)
       , IxI  = c(-Inf, -4:4, Inf)
       , IxI2 = c(-Inf, -4:4, Inf,Inf))
  ##
, mkMx(M = .Machine$double.xmax)
, mkMx(M = 1e6)
, mkMx(M = 1e9)
, mkMx(M = 1e12)
, mkMx(M = 1e14)
, mkMx(M = 1e16)
, mkMx(M = 1e20)
, mkMx(M = 1e40)
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

