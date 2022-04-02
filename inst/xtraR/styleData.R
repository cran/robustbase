### Small Test Datasets (all kinds odd/even, constant/regular/outlier):

D <- within(list(), {
    ## n = 0,1,2,3 :
    x0 <- numeric(0)
    x1 <- 3
    x1I <- Inf
    x2 <- 1:2
    x2I <- c(-Inf, 9)
    xII <- c(-Inf, Inf)
    x3 <- c(1:2,10)
    x3I <- c(-Inf, 9,11)
    x3.2I <- c(-Inf, 9, Inf)
    ## constant (0 mad) + 0--2 outliers
    xC <-  rep(1, 12)
    xC. <- rep(1, 11)
    xC1  <- c(xC,  10)
    xC1. <- c(xC., 10)
    xC2  <- c(xC1,  100)
    xC2. <- c(xC1., 100)
    ## "uniform"  + 0--2 outliers
    y  <- 1:10
    y. <- 1:11
    y1  <- c(y,  100)
    y1. <- c(y., 100)
    y2  <- c(y1,  1000)
    y2. <- c(y1., 1000)
    yI  <- c(y1,  Inf)
    yI. <- c(y1., Inf)

})
smallD <- D[order(lengths(D))]
rm(D)

## Constructor of such "stylized" small data with large ('M') values / outliers:

mk3Mx <- function(M, ngood = 10, left = floor(ngood/3)) {
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

## The one that works for a *vector* M:
mkMx <- function(M, ngood = 10, left = floor(ngood/3)) unlist(lapply(M, mk3Mx), recursive=FALSE)
