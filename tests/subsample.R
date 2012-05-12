### test subsample
### LU decomposition and singular subsamples handling
require(robustbase)

gotMatrix <- require(Matrix) ## to test lu decomposition

## R version of LU decomposition in subsample() in lmrob.c
## Modified from Golub G. H., Van Loan, C. F., Matrix Computations, 3rd edition
LU.gaxpy <- function(A, pivot=TRUE) {
    A <- as.matrix(A)
    stopifnot(ncol(A) >= nrow(A))
    m <- nrow(A)
    n <- ncol(A)
    v <- double(m) ## work matrix
    ## these matrices will contain the results
    L <- diag(m)
    U <- matrix(0, m, m)
    p <- integer(m-1) ## pivots
    idc <- 1L:n ## which columns of A are used
    idr <- 1L:m ## how to rows of A are permuted
    lenidc <- length(idc)
    for(j in 1L:m) {
        sing <- TRUE
        while(sing) {        
            if (j == 1L) {
                v[j:m] <- A[idr[j:m], idc[j]]
            } else {
                rows <- 1L:(j-1L)
                z <- forwardsolve(L[rows, rows, drop=FALSE], A[idr[rows], idc[j]])
                ## cat("Step", j, "z=", z, "\n");
                U[rows, j] <- z
                v[j:m] <- A[idr[j:m], idc[j]] - L[j:m, rows, drop=FALSE] %*% z
                ## cat("v=", v, "\n");
            }
            if (j < m) {
                mu <- j
                mu <- if (pivot) which.max(abs(v[j:m])) + j - 1L else j
                ## debug possumDiv example
                ## if (j == 37) {
                ##     print(round(abs(v[j:m]), 6)); print(max(abs(v[j:m])))
                ##     print(c(mu, v[mu]))
                ##     print(c(abs(v[47]) > abs(v[51]), abs(v[47]) < abs(v[51])))
                ## }
                if (abs(v[mu]) >= 1e-10) { ## singular: can stop here already
                    p[j] <- mu
                    if (pivot) {
                        tmp <- v[j]; v[j] <- v[mu]; v[mu] <- tmp
                        tmp <- idr[j]; idr[j] <- idr[mu]; idr[mu] <- tmp
                    }
                    L[(j+1L):m, j] <- v[(j+1L):m]/v[j]
                    if (pivot & j > 1) {
                        tmp <- L[j, rows]; L[j, rows] <- L[mu, rows]; L[mu, rows] <- tmp
                    }
                }
            }
            U[j, j] <- v[j]
            if (abs(v[j]) < 1e-10) {
                ##warning("singularity detected in step ", j, " candidate: ", idc[j])
                idc[j] <- idc[lenidc]
                lenidc <- lenidc - 1
                if (lenidc < j)
                    break
            } else sing <- FALSE
        }
    }
    list(L = L, U = U, p = p, idc = idc[1L:m], singular = sing)
}

Rsubsample <- function(x, y, mts=0) {
    x <- as.matrix(x)
    n <- length(y)
    p <- ncol(x)
    
    .C(robustbase:::R_subsample,
       x=as.double(x),
       y=as.double(y),
       n=n,
       m=p,
       beta=double(p),
       ind_space=integer(n),
       idc=integer(n),
       idr=integer(n),
       lu=double(p^2),
       v=double(p),
       pivot=integer(p-1),
       status=integer(1),
       sample=FALSE,
       mts=as.integer(mts),
       ss=as.integer(mts == 0))
}

subsample <- function(x, y=rnorm(nrow(x))) {
    x <- as.matrix(x)
    n <- length(y)
    p <- ncol(x)
    
    z <- Rsubsample(x, y)
    ## convert idc, idr and p to 1-based indexing
    idr <- z$idr + 1
    idc <- z$idc[1:p] + 1
    pivot <- z$p + 1
    ## get L and U
    L <- U <- LU <- matrix(z$lu, p, p)
    L[upper.tri(L, diag=TRUE)] <- 0
    diag(L) <- 1
    U[lower.tri(U, diag=FALSE)] <- 0
    
    cmp <- LU.gaxpy(t(x))
    if (!isTRUE(all.equal(cmp$p, pivot))) {
        print(cmp$p)
        print(pivot)
        print(which(cmp$p != pivot))
    }
    stopifnot(all.equal(cmp$L, L),
              all.equal(cmp$U, U),
              all.equal(cmp$p, pivot),
              all.equal(cmp$idc, idc))

    xsub <- x[idc, ]
    ysub <- y[idc]
    
    ## compare with Matrix result
    if (gotMatrix & !cmp$singular) {
        tmp <- lu(t(xsub))
        stopifnot(all.equal(tmp@x, z$lu))
    }

    ## test solved parameter
    if (!cmp$singular) {
      stopifnot(all.equal(z$beta, solve(xsub, ysub)))
    }

    invisible(z)
}

A <- matrix(c(0.001, 1, 1, 2), 2)
subsample(A)

A <- matrix(c(3, 2, 6, 17, 4, 18, 10, -2, 12), 3)
subsample(A)

## test some random matrix
set.seed(1002)
A <- matrix(rnorm(100), 10)
subsample(A)

## test singular matrix handling
A <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1), 4, byrow=TRUE)
subsample(A)


## test subsample with mts > 0
data <- data.frame(y = rnorm(9), expand.grid(A = letters[1:3], B = letters[1:3]))
x <- model.matrix(y ~ ., data)
y <- data$y
## this should produce a warning and return status == 2
z <- Rsubsample(x, y, mts=2)
stopifnot(z$status == 2)

## test real data example
data(possumDiv)## 151 * 9; the last two variables are factors
with(possumDiv, table(eucalyptus, aspect))

mf <- model.frame(Diversity ~ .^2, possumDiv)
X <- model.matrix(mf, possumDiv)
y <- model.response(mf)
stopifnot(qr(X)$rank == ncol(X))

## subsample(X, y) ## fails: different pivots in step 37

## test subsampling
testSubSampling <- function(X, y) {
    lX <- X[sample(nrow(X)), ]
    ## C version
    zc <- Rsubsample(lX, y)
    ## R version
    zR <- LU.gaxpy(t(lX))
    if (as.logical(zc$status)) {
        ## singularity in C detected
        if (!zR$singular)
            stop("singularity in C but not in R")
    } else {
        ## no singularity detected
        if (zR$singular)
            stop("singularity in R but not in C")
    }
    zR$singular
}

set.seed(10)
nsing <- 0
for (i in 1:200) if (testSubSampling(X, y)) nsing <- nsing + 1
stopifnot(nsing == 0)

## lmrob.S fails for this seed:
## set.seed(108)
## lmrob(Diversity ~ .^2 , data = possumDiv, cov="none", trace=3)
