
lmrob.control <-
    function(seed = 37, ## < was hard-coded in C
	     Nres = 500,
	     tuning.chi = 1.54764, bb = 0.5, tuning.psi = 4.685061,
	     max.it = 50, groups = 5, n.group = 400, k.fast.s = 1,
	     compute.rd = TRUE)
{
    list(seed = seed, Nres = Nres, tuning.chi = tuning.chi,
	 bb = bb, tuning.psi = tuning.psi, groups = groups, n.group = n.group,
	 k.fast.s = k.fast.s, max.it = max.it, compute.rd = compute.rd)
}


lmrob.fit.MM <- function(x, y, control)
{
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    ## find the initial (50% BP S) estimator and the
    ## associated residual scale estimator
    initial.S <- lmrob.S(x = x, y = y, control = control)

    ## find the final (95% efficient M) estimator using
    ## RWLS iterations
    final.MM <- lmrob.MM(x = x, y = y, beta.initial = initial.S$coef,
			 scale = initial.S$scale, control = control)

    coef.initial <- initial.S$coef

    ## If IRWLS does not converge, use the S estimator
    if( !final.MM$converged) {
	coef <- initial.S$coef
	cov.matrix <- matrix(initial.S$cov, p, p)
    } else {
	coef <- final.MM$coef
	cov.matrix <- matrix(final.MM$cov, p, p)
    }

    rank <- qr(x)$rank ## << FIXME? qr(.)  should be available from earlier
    r1 <- 1:rank
    dn <- colnames(x)
    if (is.matrix(y)) {
	coef[-r1, ] <- NA
	dimnames(coef) <- list(dn, colnames(y))
    }
    else {
	coef[-r1] <- NA
	names(coef) <- dn
    }
    f <- x %*% as.matrix(coef)
    r <- as.matrix(y) - f
    list(fitted.values = as.vector(f), residuals = as.vector(r),
	 rank = rank, degree.freedom = n-rank, coefficients = coef,
	 initial.coefficients = coef.initial,
	 scale = final.MM$scale, seed = final.MM$seed, cov = cov.matrix,
	 converged = final.MM$converged)
}


lmrob.MM <- function(x, y, beta.initial, scale, control)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    c.psi <- control$tuning.psi
    c.chi <- control$tuning.chi
    sigma <- scale

    b <- .C("R_lmrob_MM",
	    x = as.double(x),
	    y = as.double(y),
	    n = as.integer(n),
	    p = as.integer(p),
	    beta.initial = as.double(beta.initial),
	    scale = as.double(scale),
	    coef = double(p),
	    as.integer(control$ max.it),
	    tuning.psi = as.double(c.psi),
	    converged = logical(1),
	    PACKAGE = "robustbase")
    r <- as.vector(y - x %*% b$coef)
    r2 <- as.vector(y - x %*% beta.initial)
    w   <- lmrob.Psi( r/sigma, cc = c.psi, deriv = 1)
    ch1 <- lmrob.Chi(r2/sigma, cc = c.chi, deriv = 1)
    A <- solve(	 ( t(x) %*% (x * w) ) / n / sigma )
    w <- w * r / sigma
    a <- A %*% (( t(x) %*% w ) / (n * mean(ch1 * r2/sigma)))
    w  <- lmrob.Psi( r/sigma, cc = c.psi)
    w2 <- lmrob.Chi(r2/sigma, cc = c.chi)
    u1 <- ( t(x) %*% (x * (w^2) ) ) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*% t( t(x) %*% (w*w2) ) %*% A / n
    u3 <- A %*% (t(x) %*% (w*w2) ) %*% t(a) / n
    u4 <- mean(w2^2 - (control$bb)^2) * a %*% t(a)

    list(coef = b$coef, cov = (u1 - u2 - u3 + u4)/n, control = control,
	 scale = b$scale, seed = b$seed, converged = b$converged )
}


lmrob.S <- function(x, y, control)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    Nres <- control$Nres
    seed <- control$seed
    bb <- control$bb
    tuning.chi <- control$tuning.chi
    groups <- control$groups
    n.group <- control$n.group
    k.fast.s <- control$k.fast.s
    b <- .C("R_lmrob_S",
	    x = as.double(x),
	    y = as.double(y),
	    n = as.integer(n),
	    p = as.integer(p),
	    Nres = as.integer(Nres),
	    scale = as.double(0),
	    coef = double(p),
	    seed = as.integer(seed),
	    tuning.chi = as.double(tuning.chi),
	    as.double(bb),
	    groups = as.integer(groups),
	    n.group = as.integer(n.group),
	    k.fast.s = as.integer(k.fast.s),
	    PACKAGE = "robustbase"
	    )
    sigma <- b$scale
    r2 <- as.vector(y - x %*% b$coef)
    w <- lmrob.Chi(r2/sigma, cc = tuning.chi, deriv = 2)
    A <- solve(	 ( t(x) %*% (x * w) ) / n / sigma )
    w <- w * r2 / sigma
    a <- ( t(x) %*% w ) / n
    w  <- lmrob.Chi(r2/sigma, cc = tuning.chi, deriv = 1)
    w2 <- lmrob.Chi(r2/sigma, cc = tuning.chi)
    a <- a / mean(w * r2/sigma)
    a <- A %*% a
    u1 <- ( t(x) %*% (x * (w^2) ) ) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*% t( t(x) %*% (w*w2) ) %*% A / n
    u3 <- A %*% (t(x) %*% (w*w2) ) %*% t(a) / n
    u4 <- mean(w2^2 - bb^2) * a %*% t(a)

    list(coef = b$coef, cov = (u1 - u2 - u3 + u4)/n,
	 control = control, scale = b$scale, seed = b$seed)
}




