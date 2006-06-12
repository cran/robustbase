
lmrob.control <-
    function(seed = NULL, ## '37' was hard-coded in C {when still using 'rand()'}
	     nResample = 500,
	     tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061,
	     max.it = 50, groups = 5, n.group = 400,
	     k.fast.s = 1,
	     best.r.s = 2,	# had '2'    hardwired in C
	     k.max = 200,	# had '50'   hardwired in C
	     refine.tol = 1e-7, # had '1e-7' hardwired in C {often not converged ?!}
	     compute.rd = FALSE)
{
    list(seed = as.integer(seed), nResample = nResample, tuning.chi = tuning.chi,
	 bb = bb, tuning.psi = tuning.psi,
	 max.it = max.it, groups = groups, n.group = n.group,
	 best.r.s = best.r.s, k.fast.s = k.fast.s,
	 k.max = k.max, refine.tol = refine.tol,
	 compute.rd = compute.rd)
}


lmrob.fit.MM <-
    function(x, y, control,
	     init.S = lmrob.S(x = x, y = y, control = control))
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    ## init.S: an initial (typically 50% BP S-) estimator and the
    ## associated residual scale estimator
    iCoef <- init.S$coef
    stopifnot(is.numeric(iCoef), length(iCoef) == p)

    ## find the final (95% efficient M) estimator using
    ## RWLS iterations
    final.MM <- lmrob.MM(x = x, y = y, beta.initial = iCoef,
			 scale = init.S$scale, control = control)

    if(final.MM$converged) {
	coef <- final.MM$coef
	cov.matrix <- matrix(final.MM$cov, p, p)
    }
    else { ## If IRWLS did not converge, use the initial (S) estimator:
	warning("IRWLS iterations did NOT converge in ", control$max.it," steps")
	coef <- iCoef
	cov.matrix <- matrix(init.S$cov, p, p)
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
    r <- y - f
    list(fitted.values = f, residuals = r, weights = final.MM$wt,
	 rank = rank, degree.freedom = n - rank,
	 coefficients = coef, initial.coefficients = iCoef,
	 scale = final.MM$scale, cov = cov.matrix, control = control,
	 iter = final.MM$iter, converged = final.MM$converged)
}


lmrob.MM <- function(x, y, beta.initial, scale, control)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    c.psi <- as.double(control$tuning.psi)
    c.chi <- as.double(control$tuning.chi)
    stopifnot(length(c.psi) > 0, length(c.chi) > 0,
	      c.psi >= 0, c.chi >= 0, scale >= 0,
	      length(beta.initial) == p)
    sigma <- scale

    b <- .C("R_lmrob_MM",
	    x = as.double(x),
	    y = as.double(y),
	    n = as.integer(n),
	    p = as.integer(p),
	    beta.initial = as.double(beta.initial),
	    scale = as.double(scale),
	    coef = double(p),
	    iter = as.integer(control$ max.it),
	    c.psi = c.psi,
	    converged = logical(1),
	    PACKAGE = "robustbase")[c("coef", "scale", "converged", "iter")]
    ## FIXME?: Should rather warn *here* in case of non-convergence
    ## FIXME: shouldn't C-code above return residuals ?!
    r.s	 <- drop(y - x %*% b$coef)	 / sigma
    r2.s <- drop(y - x %*% beta.initial) / sigma
    w	<- lmrob.Psi( r.s, cc = c.psi, deriv = 1)
    ch1 <- lmrob.Chi(r2.s, cc = c.chi, deriv = 1)
    ## FIXME for multivariate y :
    A <- solve(	crossprod(x, x * w) ) * (n * sigma)
    a <- A %*% (crossprod(x, w * r.s) / (n * mean(ch1 * r2.s)))
    w  <- lmrob.Psi( r.s, cc = c.psi)
    w2 <- lmrob.Chi(r2.s, cc = c.chi)
    Xww <- crossprod(x, w*w2)
    u1 <- crossprod(x, x * w^2) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*%	 crossprod(Xww, A) / n
    u3 <- A %*% tcrossprod(Xww, a) / n
    u4 <- mean(w2^2 - control$bb ^2) * tcrossprod(a)

    c(b, list(cov = (u1 - u2 - u3 + u4)/n,
	      wt = w / r.s, control = control))
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
	stop("'groups * n.group' must be larger than 'n' for 'large_n' algorithm")
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
	      length(best.r)> 0, best.r >= 1)

    b <- .C("R_lmrob_S",
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
	    trace.lev = as.integer(trace.lev),
	    PACKAGE = "robustbase")[c("coef", "scale", "k.max", "converged")]
    sigma <- b$scale
    if(sigma < 0)
	stop("C function R_lmrob_S() exited prematurely")
    names(b)[names(b) == "k.max"] <- "k.iter" # maximal #{refinement iter.}
    ## FIXME: get 'res'iduals from C
    r2.s <- drop(y - x %*% b$coef) / sigma
    w <- lmrob.Chi(r2.s, cc = c.chi, deriv = 2)
    A <- solve(	crossprod(x, x * w) ) * (n * sigma)
    a <- crossprod(x, w * r2.s) / n
    w  <- lmrob.Chi(r2.s, cc = c.chi, deriv = 1)
    w2 <- lmrob.Chi(r2.s, cc = c.chi)
    a <- A %*% (a/ mean(w * r2.s))
    Xww <- crossprod(x, w*w2)
    u1 <- crossprod(x, x * w^2) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*%	 crossprod(Xww, A) / n
    u3 <- A %*% tcrossprod(Xww, a) / n
    u4 <- mean(w2^2 - (control$bb)^2) * tcrossprod(a)

    c(b, list(cov = (u1 - u2 - u3 + u4)/n, control = control))
}
