
lmrob.control <-
    function(seed = NULL, ## '37' was hard-coded in C {when still using 'rand()'}
	     nResample = 500,
	     tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061,
	     max.it = 50, groups = 5, n.group = 400,
	     k.fast.s = 1,
	     best.r.s = 2,	# had '2'    hardwired in C
	     k.max = 200,	# had '50'   hardwired in C
	     refine.tol = 1e-7, # had '1e-7' hardwired in C {often not converged ?!}
             rel.tol = 1e-7,    # had '1e-7' hardwired in C for ABSOLUTE tol
             trace.lev = 0,
	     compute.rd = FALSE)
{
    list(seed = as.integer(seed), nResample = nResample, tuning.chi = tuning.chi,
	 bb = bb, tuning.psi = tuning.psi,
	 max.it = max.it, groups = groups, n.group = n.group,
	 best.r.s = best.r.s, k.fast.s = k.fast.s,
	 k.max = k.max, refine.tol = refine.tol, rel.tol = rel.tol,
         trace.lev = trace.lev, compute.rd = compute.rd)
}


lmrob.fit.MM <- function(x, y, control)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    ## init.S: an initial (typically 50% BP S-) estimator and the
    ## ------ associated residual scale estimator
    ## This is *not* an function argument, since vcov() inference
    ## depends on exact S estimator's constants!
    init.S <- lmrob.S(x = x, y = y, control = control)

    iCoef <- init.S$coef
    c.chi <- as.double(control$tuning.chi)
    c.psi <- as.double(control$tuning.psi)

    stopifnot(is.numeric(iCoef), length(iCoef) == p,
	      is.numeric(init.S$scale), init.S$scale >= 0,
	      length(c.chi) > 0, c.chi >= 0)
					# c.psi is checked in lmrob..M..fit

    ## find the final (95% efficient M) estimator using
    ## RWLS iterations
    MM.fit <- lmrob..M..fit(x = x, y = y, beta.initial = iCoef,
                            scale = init.S$scale,
                            c.psi = c.psi, control = control)
    ##>> comp. ("coef", "scale", "resid", "loss", "converged", "iter")
    scale <- MM.fit$scale

    if(MM.fit$converged) {
	coef  <- MM.fit$coef
        ## this inference *depends* the initial S constants !
        vv <- .vcov.MM(r0 = drop(y - x %*% iCoef),# initial scaled residuals
                       r = MM.fit$resid, #          final   scaled residuals
                       scale=scale, x=x,
                       c.psi=c.psi, c.chi=c.chi, bb=control$bb)
	cov.matrix <- vv$cov
        weights <- vv$wt
    }
    else { ## If IRWLS did not converge, use the initial (S) estimator:
	warning("IRWLS iterations (for MM) did NOT converge in ",
                control$max.it," steps; using initial estimate instead.")
	coef <- iCoef
        ## TODO?: only now compute  vcov(init.S) instead of always
	cov.matrix <- matrix(init.S$cov, p, p)

        ## still return the "unconverged weights":
        r.s <- MM.fit$resid / scale # final   scaled residuals
        weights <- tukeyPsi1(r.s, cc = c.psi) / r.s
    }

    rank <- qr(x)$rank ## << FIXME? qr(.)  should be available from earlier
    r1 <- 1:rank
    dn <- colnames(x)
    if (is.matrix(y)) { ## multivariate!
	coef[-r1, ] <- NA
	dimnames(coef) <- list(dn, colnames(y))
    }
    else {
	coef[-r1] <- NA
	names(coef) <- dn
    }
    f <- drop(x %*% coef)
    list(fitted.values = f, residuals = y - f, weights = weights,
         rank = rank, degree.freedom = n - rank,
	 coefficients = coef, initial.coefficients = iCoef,
	 scale = scale, cov = cov.matrix, control = control,
	 iter = MM.fit$iter, converged = MM.fit$converged,
         ## return initial estimate; notably its own 'converged'
         init.S = init.S)
}## lmrob.MM.fit()

.vcov.MM <- function(r0, r, scale, x, c.psi,c.chi, bb)
{
    ## scaled residuals
    n <- length(r)
    stopifnot(n == length(r0), is.matrix(x), n == nrow(x))
    p <- ncol(x)
    r.s	 <- r / scale # final   scaled residuals
    r0.s <- r0/ scale # initial scaled residuals
    w  <- tukeyPsi1( r.s, cc = c.psi, deriv = 1)
    w0 <- tukeyChi (r0.s, cc = c.chi, deriv = 1)
    ## FIXME for multivariate y :
    A <- solve(	crossprod(x, x * w) ) * (n * scale)
    a <- A %*% (crossprod(x, w * r.s) / (n * mean(w0 * r0.s)))
    w  <- tukeyPsi1( r.s, cc = c.psi)

    ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, bb)
    w0 <- tukeyChi(r0.s, cc = c.chi)
    Xww <- crossprod(x, w*w0)
    u1 <- crossprod(x, x * w^2) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*%	 crossprod(Xww, A) / n
    u3 <- A %*% tcrossprod(Xww, a) / n
    u4 <- mean(w0^2 - bb^2) * tcrossprod(a)

    list(cov = 	matrix((u1 - u2 - u3 + u4)/n, p, p),
         wt = w / r.s, a = a)
}

lmrob..M..fit <- function(x, y, beta.initial, scale, c.psi, control)
{
    stopifnot(is.matrix(x))
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(length(y) == n,
              length(c.psi) > 0, c.psi >= 0,
              scale >= 0, length(beta.initial) == p)

    .C(R_lmrob_MM,
       x = as.double(x),
       y = as.double(y),
       n = as.integer(n),
       p = as.integer(p),
       beta.initial = as.double(beta.initial),
       scale = as.double(scale),
       coef = double(p),
       resid = double(n),
       iter = as.integer(control$ max.it),
       c.psi = c.psi,
       loss = double(1),
       rel.tol= as.double(control$rel.tol),
       converged = logical(1),
       trace.lev = as.integer(control$trace.lev)
       )[c("coef", "scale", "resid", "loss", "converged", "iter")]
    ## FIXME?: Should rather warn *here* in case of non-convergence
}


lmrob.S <- function(x, y, control, trace.lev = 0)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    nResample <- as.integer(control$nResample)
    groups <- as.integer(control$groups)
    nGr <- as.integer(control$n.group)
    if(nGr <= p) stop("'control$n.group' must be larger than 'p'")
    large_n <- (n > 2000)
    if(large_n & nGr * groups > n)
	stop("'groups * n.group' must be smaller than 'n' for 'large_n' algorithm")
    if(nGr <= p + 10) ## FIXME (be smarter ..)
	warning("'control$n.group' is not much larger than 'p', probably too small")
    if(length(seed <- control$seed) > 0) {
	if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
	    seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
	    on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
	}
	assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    bb <- as.double(control$bb)
    c.chi <- as.double(control$tuning.chi)
    best.r <- as.integer(control$best.r.s)
    stopifnot(length(c.chi) > 0, c.chi >= 0, length(bb) > 0,
	      length(best.r) > 0, best.r >= 1)

    b <- .C(R_lmrob_S,
	    x = as.double(x),
	    y = as.double(y),
	    n = as.integer(n),
	    p = as.integer(p),
	    nResample = nResample,
	    scale = double(1),
	    coef = double(p),
	    c.chi,
	    bb,
	    best_r = best.r,
	    groups  = groups,
	    n.group = nGr,
	    k.fast.s= as.integer(control$k.fast.s),
	    k.max   = as.integer(control$k.max),
	    refine.tol= as.double(control$refine.tol),
	    converged = logical(1),
	    trace.lev = as.integer(trace.lev)
	    )[c("coef", "scale", "k.max", "converged")]
    scale <- b$scale
    if(scale < 0)
	stop("C function R_lmrob_S() exited prematurely")
    names(b)[names(b) == "k.max"] <- "k.iter" # maximal #{refinement iter.}
    ## FIXME: get 'res'iduals from C

    r0.s <- drop(y - x %*% b$coef) / scale
    w <- tukeyChi(r0.s, cc = c.chi, deriv = 2)
    ## FIXME for multivariate y :
    A <- solve(	crossprod(x, x * w) ) * (n * scale)
    a <- crossprod(x, w * r0.s) / n # before 'w' is re-assigned
    w  <- tukeyChi(r0.s, cc = c.chi, deriv = 1)
    a <- A %*% (a/ mean(w * r0.s))

    ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, control$bb)
    w0 <- tukeyChi(r0.s, cc = c.chi)
    Xww <- crossprod(x, w*w0)
    u1 <- crossprod(x, x * w^2) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*%	 crossprod(Xww, A) / n
    u3 <- A %*% tcrossprod(Xww, a) / n
    u4 <- mean(w0^2 - control$bb^2) * tcrossprod(a)

    c(b, list(cov = (u1 - u2 - u3 + u4)/n, control = control))
}
