## Test implementation of M-S estimator
require(robustbase)
lmrob.conv.cc  <- robustbase::: lmrob.conv.cc
lmrob.psi2ipsi <- robustbase::: lmrob.psi2ipsi
lmrob.wgtfun   <- robustbase::: lmrob.wgtfun

## dataset with factors and continuous variables:
data(education)
education <- within(education, Region <- factor(Region))
## for testing purposes:
education2 <- within(education, Group <- factor(rep(1:3, length.out=length(Region))))

## Test splitFrame (type fii is the only problematic type)
testFun <- function(formula, x1.idx) {
    obj <- lm(formula, education2)
    mf <- obj$model
    ret <- splitFrame(mf, type="fii")
    if (missing(x1.idx)) {
        print(ret$x1.idx)
        return(which(unname(ret$x1.idx)))
    }
    stopifnot(identical(x1.idx, which(unname(ret$x1.idx))))
}
testFun(Y ~ 1, integer(0))
testFun(Y ~ X1*X2*X3, integer(0))
testFun(Y ~ Region + X1 + X2 + X3, 1:4)
testFun(Y ~ 0 + Region + X1 + X2 + X3, 1:4)
testFun(Y ~ Region*X1 + X2 + X3, c(1:5, 8:10))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group, c(1:5, 8:18))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group*X2, c(1:6, 8:29))
testFun(Y ~ Region*X1 + X2 + Region*Group*X2, 1:28)
testFun(Y ~ Region*X1 + X2 + Region:Group:X2, 1:21)
testFun(Y ~ Region*X1 + X2*X3 + Region:Group:X2, c(1:6, 8:10, 12:23))
testFun(Y ~ (X1+X2+X3+Region)^2, c(1:7,10:12,14:19))
testFun(Y ~ (X1+X2+X3+Region)^3, c(1:19, 21:29))
testFun(Y ~ (X1+X2+X3+Region)^4, 1:32)
testFun(Y ~ Region:X1:X2 + X1*X2, c(1:1, 4:7))

## Test subsampling algorithm
m_s_subsample <- function(x1, x2, y, control, orthogonalize=TRUE) {
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    y <- y
    storage.mode(x1) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(y) <- "double"

    z <- .C(robustbase:::R_lmrob_M_S,
            X1=x1,
            X2=x2,
            y=y,
            res=double(length(y)),
            n=length(y),
            p1=ncol(x1),
            p2=ncol(x2),
            nResample=as.integer(control$nResample),
            scale=double(1),
            b1=double(ncol(x1)),
            b2=double(ncol(x2)),
            tuning_chi=as.double(control$tuning.chi),
            ipsi=as.integer(lmrob.psi2ipsi(control$psi)),
            bb=as.double(control$bb),
            K_m_s=as.integer(control$k.m_s),
            max_k=as.integer(control$k.max),
            rel_tol=as.double(control$rel.tol),
            converged=logical(1),
            trace_lev=as.integer(control$trace.lev),
            orthogonalize=as.logical(orthogonalize),
            subsample=TRUE,
            descent=FALSE,
            mts = 0L,
            ss = 1L)
    z[c("b1", "b2", "scale")]
}

control <- lmrob.control()
f.lm <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- splitFrame(f.lm$model)
y <- education$Y

## test orthogonalizing
x1 <- splt$x1
x2 <- splt$x2
tmp <- lmrob.lar(x1, y, control)
y.tilde <- tmp$resid
t1 <- tmp$coef
x2.tilde <- x2
T2 <- matrix(0, nrow=ncol(x1), ncol=ncol(x2))
for (i in 1:ncol(x2)) {
    tmp <- lmrob.lar(x1, x2[,i], control)
    x2.tilde[,i] <- tmp$resid
    T2[,i] <- tmp$coef
}
set.seed(10)
res1 <- m_s_subsample(x1, x2.tilde, y.tilde, control, FALSE)
res1 <- within(res1, b1 <- drop(t1 + b1 - T2 %*% b2))
set.seed(10)
res2 <- m_s_subsample(x1, x2, y, control, TRUE)
stopifnot(all.equal(res1, res2))

res <- vector("list", 100)
set.seed(0)
time <- system.time(for (i in 1:100) {
    tmp <- m_s_subsample(x1, x2.tilde, y.tilde, control, FALSE)
    res[[i]] <- unlist(within(tmp, b1 <- drop(t1 + b1 - T2 %*% b2)))
})
cat('Time elapsed in subsampling: ', time,'\n')
## show a summary of the results
summary(res1 <- do.call(rbind, res))
## compare with fast S solution
fmS <- lmrob(Y ~ Region + X1 + X2 + X3, education, init="S")
coef(fmS)
fmS$scale

## Test descent algorithm
m_s_descent <- function(x1, x2, y, control, b1, b2, scale) {
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    y <- y
    storage.mode(x1) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(y) <- "double"

    z <- .C(robustbase:::R_lmrob_M_S,
            X1=x1,
            X2=x2,
            y=y,
            res=double(length(y)),
            n=length(y),
            p1=ncol(x1),
            p2=ncol(x2),
            nResample=as.integer(control$nResample),
            scale=as.double(scale),
            b1=as.double(b1),
            b2=as.double(b2),
            tuning_chi=as.double(control$tuning.chi),
            ipsi=as.integer(lmrob.psi2ipsi(control$psi)),
            bb=as.double(control$bb),
            K_m_s=as.integer(control$k.m_s),
            max_k=as.integer(control$k.max),
            rel_tol=as.double(control$rel.tol),
            converged=logical(1),
            trace_lev=as.integer(control$trace.lev),
            orthogonalize=FALSE,
            subsample=FALSE,
            descent=TRUE,
            mts = 0L,
            ss = 1L)
    z[c("b1", "b2", "scale", "res")]
}

find_scale <- function(r, s0, n, p, control) {
    c.chi <- lmrob.conv.cc(control$psi, control$tuning.chi)

    b <- .C(robustbase:::R_lmrob_S,
            x = double(1),
            y = as.double(r),
            n = as.integer(n),
            p = as.integer(p),
            nResample = 0L,
            scale = as.double(s0),
            coefficients = double(p),
            as.double(c.chi),
            as.integer(lmrob.psi2ipsi(control$psi)),
            as.double(control$bb),
            best_r = 0L,
            groups = 0L,
            n.group = 0L,
            k.fast.s = 0L,
            k.iter = 0L,
            refine.tol = as.double(control$refine.tol),
            converged = logical(1),
            trace.lev = 0L,
            mts = 0L,
            ss = 1L
            )[c("coefficients", "scale", "k.iter", "converged")]
    b$scale
}

## what should it be:
m_s_descent_Ronly<- function(x1, x2, y, control, b1, b2, scale) {
    n <- length(y)
    p1 <- ncol(x1)
    p2 <- ncol(x2)
    p <- p1+p2
    t2 <- b2
    t1 <- b1
    rs <- drop(y - x1 %*% b1 - x2 %*% b2)
    sc <- scale
    ## do refinement steps
    ## do maximally control$k.max iterations
    ## stop if converged
    ## stop after k.fast.m_s step of no improvement
    if (control$trace.lev > 4) cat("scale:", scale, "\n")
    if (control$trace.lev > 4) cat("res:", rs, "\n")
    nnoimprovement <- nref <- 0; conv <- FALSE
    while((nref <- nref + 1) <= control$k.max && !conv &&
          nnoimprovement < control$k.m_s) {
        ## STEP 1: UPDATE B2
        y.tilde <- y - x1 %*% t1
        w <- lmrob.wgtfun(rs / sc, control$tuning.chi, control$psi)
        if (control$trace.lev > 4) cat("w:", w, "\n")
        z2 <- lm.wfit(x2, y.tilde, w)
        t2 <- z2$coef
        if (control$trace.lev > 4) cat("t2:", t2, "\n")
        rs <- y - x2 %*% t2
        ## STEP 2: OBTAIN M-ESTIMATE OF B1
        z1 <- lmrob.lar(x1, rs, control)
        t1 <- z1$coef
        if (control$trace.lev > 4) cat("t1:", t1, "\n")
        rs <- z1$resid
        ## STEP 3: COMPUTE THE SCALE ESTIMATE
        sc <- find_scale(rs, sc, n, p, control)
        if (control$trace.lev > 4) cat("sc:", sc, "\n")
        ## STEP 4: CHECK FOR CONVERGENCE
        #...
        ## STEP 5: UPDATE BEST FIT
        if (sc < scale) {
            scale <- sc
            b1 <- t1
            b2 <- t2
            nnoimprovement <- 0
        } else nnoimprovement <- nnoimprovement + 1
    }
    ## STEP 6: FINISH
    if (nref == control$k.max)
        warning("M-S estimate: maximum number of refinement steps reached.")

    list(b1=b1, b2=b2, scale=scale, res=rs)
}

control2 <- control
#control2$trace.lev <- 5
control2$k.max <- 1
stopifnot(all.equal(m_s_descent      (x1, x2, y, control2, res2$b1, res2$b2, res2$scale+10),
		    m_s_descent_Ronly(x1, x2, y, control2, res2$b1, res2$b2, res2$scale+10),
		    check.attr = FALSE))

## control$k.m_s <- 100
res3 <- vector("list", 100)
time <- system.time(for (i in 1:100) {
    res3[[i]] <- unlist(m_s_descent(x1, x2, y, control,
                                    res[[i]][1:4], res[[i]][5:7], res[[i]][8]))
})
cat('Time elapsed in descent proc: ', time,'\n')

## show a summary of the results
res4 <- do.call(rbind, res3)
summary(res4[,1:8])

plot(res1[, "scale"], res4[,"scale"])
abline(0,1, col=adjustcolor("gray", 0.5))

## Test lmrob.M.S
x <- model.matrix(fmS)
control$trace.lev <- 3
set.seed(1003)
fMS <- lmrob.M.S(x, y, control, fmS$model)
resid <- drop(y - x %*% fMS$coef)
stopifnot(all.equal(resid, fMS$resid, check.attr=FALSE))

## Test direct call to lmrob
set.seed(13)
fiMS <- lmrob(Y ~ Region + X1 + X2 + X3, education, init="M-S")
out2 <- capture.output(summary(fiMS))
writeLines(out2)

set.seed(13)
fiM.S <- lmrob(Y ~ Region + X1 + X2 + X3, education, init=lmrob.M.S)
out3 <- capture.output(summary(fiM.S))

## must be the same {apart from the "init=" in the call}:
stopifnot(identical(out2[-4], out3[-4]))
