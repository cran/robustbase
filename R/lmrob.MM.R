lmrob.control <- function  (setting, seed = NULL, nResample = 500,
                            tuning.chi = NULL,  bb = 0.5,
                            tuning.psi = NULL, max.it = 50,
                            groups = 5, n.group = 400, k.fast.s = 1, best.r.s = 2,
                            k.max = 200, k.m_s = 20, refine.tol = 1e-07, rel.tol = 1e-07,
                            trace.lev = 0, mts = 1000,
                            subsampling = c("constrained", "simple"),
                            compute.rd = FALSE, method = 'MM',
                            psi = c('bisquare', 'lqq', 'welsh', 'optimal', 'hampel',
                              'ggw'),
                            numpoints = 10, cov = '.vcov.avar1',
                            split.type = c("f", "fi", "fii"),
                            ...)
{
    if (!missing(setting)) {
        if (setting == 'KS2011') {
            if (missing(method)) method <- 'SMDM'
            if (missing(psi)) psi <- 'lqq'
            if (missing(max.it)) max.it <- 500
            if (missing(k.max)) k.max <- 2000
            if (missing(cov)) cov <- '.vcov.w'
        } else {
            warning("Unknown setting '", setting, "'. Using defaults.")
            psi <- match.arg(psi)
        }
    } else {
        psi <- if (missing(psi) && grepl('D', method)) 'lqq' else match.arg(psi)
        if (missing(cov) && !method %in% c('SM', 'MM')) cov <- '.vcov.w'
    }
    subsampling <- match.arg(subsampling)
    
    if (missing(tuning.chi) || is.null(tuning.chi))
        tuning.chi <- switch(psi,
                             'bisquare' = 1.54764,
                             'welsh' = 0.5773502,
                             'ggw' = c(-0.5, 1.5, NA, .5), ## min slope, b, eff, bp
                             'lqq' = c(-0.5, 1.5, NA, .5), ## min slope, b, eff, bp
                             'optimal' = 0.4047,
                             'hampel' = c(1.5, 3.5, 8) * 0.2119163) ## a, b, c
    if (missing(tuning.psi) || is.null(tuning.psi))
        tuning.psi <- switch(psi,
                             'bisquare' = 4.685061,
                             'welsh' = 2.11,
                             'ggw' = c(-0.5, 1.5, .95, NA), ## min slope, b, eff, bp
                             'lqq' = c(-0.5, 1.5, .95, NA), ## min slope, b, eff, bp
                             'optimal' = 1.060158,
                             'hampel' = c(1.5, 3.5, 8) * 0.9016085) ## a, b, c

    ## in ggw, lqq:  tuning.psi, *.chi  are non-standard, calculate coefficients:
    if (psi %in% c('ggw', 'lqq')) {
        tuning.psi <- lmrob.const(tuning.psi, psi)
        tuning.chi <- lmrob.const(tuning.chi, psi)
    }

    c(list(seed = as.integer(seed), nResample = nResample, psi = psi,
           tuning.chi = tuning.chi, bb = bb, tuning.psi = tuning.psi,
           max.it = max.it, groups = groups, n.group = n.group,
           best.r.s = best.r.s, k.fast.s = k.fast.s,
           k.max = k.max, k.m_s = k.m_s, refine.tol = refine.tol,
           rel.tol = rel.tol, trace.lev = trace.lev, mts = mts,
           subsampling = subsampling,
           compute.rd = compute.rd, method = method, numpoints = numpoints,
           cov = cov, split.type = match.arg(split.type)),
      list(...))
}

lmrob.fit.MM <- function(x, y, control) ## deprecated
{
    .Deprecated("lmrob.fit(*, control) with control$method = 'SM'")
    control$method <- 'SM'
    lmrob.fit(x, y, control)
}## lmrob.MM.fit()

lmrob.fit <- function(x, y, control, init=NULL) {
    if(!is.matrix(x)) x <- as.matrix(x)
    ## old notation: MM -> SM
    if (control$method == "MM") control$method <- "SM"
    if (is.null(init)) {
        ## --- initial S estimator
        if (substr(control$method,1,1) != 'S') {
            warning("Initial estimator '", substr(control$method,1,1), "' not supported",
                    " using S-estimator instead")
            substr(control$method,1,1) <- 'S'
        }
        init <- lmrob.S(x,y,control=control)
        est <- 'S'
    } else {
        if (is.null(init$converged)) init$converged <- TRUE
        if (is.null(init$control)) {
            init$control <- control
            init$control$method <- ''
        }
        est <- init$control$method
    }
    stopifnot(is.numeric(init$coef), length(init$coef) == ncol(x),
              is.numeric(init$scale), init$scale >= 0)
    if (est != 'S' && control$cov == '.vcov.avar1') {
        warning("Can only use .vcov.avar1 for S as initial estimator",
                " using .vcov.w instead")
        control$cov <- ".vcov.w"
    }
    if (init$converged) {
        ## --- loop through the other estimators
        method <- sub(est, '', control$method)
        for (step in strsplit(method,'')[[1]]) {
            ## now we have either M or D steps
            est <- paste(est, step, sep = '')
            init <- switch(step,
                           ## D(AS)-Step
                           D = lmrob..D..fit(init, x),
                           ## M-Step
                           M = lmrob..M..fit(x = x, y = y, obj=init),
                           stop('only M and D are steps supported'))
            ## break if an estimator did not converge
            if (!init$converged) {
                warning(step, "-step did NOT converge. Returning unconverged ", est,
                        "-estimate.")
                break
            }
        }
    }
    ## << FIXME? qr(.)  should be available from earlier
    if (is.null(init$qr)) init$qr <- qr(x * sqrt(init$weights))
    if (is.null(init$rank)) init$rank <- init$qr$rank

    ## --- covariance estimate
    if (!init$converged || is.null(x)) {
        init$cov <- NA
        init$df <- init$degree.freedom <- NA
    } else {
        control$method <- est
        init$control <- control
        init$cov <-
            if (is.null(control$cov) || control$cov == "none") NA
            else {
                lf.cov <- if (!is.function(control$cov))
                    get(control$cov, mode='function') else control$cov
                lf.cov(init, x)
            }
        df <- NROW(y) - init$rank ## sum(init$weights)-init$rank
        init$degree.freedom <- init$df.residual <- df
    }

    init
}

if(getRversion() > "2.15.0" || as.numeric(R.Version()$`svn rev`) > 59233)
    utils::globalVariables("r", add=TRUE) ## below and in other lmrob.E() expressions

.vcov.w <- function(obj, x=obj$x, scale=obj$scale, cov.hubercorr=ctrl$cov.hubercorr,
             cov.dfcorr=ctrl$cov.dfcorr, cov.resid=ctrl$cov.resid,
             cov.corrfact=ctrl$cov.corrfact,
             cov.xwx=ctrl$cov.xwx)
{
    ## set defaults
    ctrl <- obj$control
    if (is.null(cov.hubercorr)) cov.hubercorr <- !grepl('D', ctrl$method)
    else if (!is.logical(cov.hubercorr))
        stop(':.vcov.w: cov.hubercorr has to be logical')
    if (is.null(cov.dfcorr)) {
        cov.dfcorr <- if (cov.hubercorr) 1 else -1
    } else if (!is.numeric(cov.dfcorr) || !cov.dfcorr %in% -1:3)
        stop(':.vcov.w: cov.dfcorr has to be one of -1:3')

    if (is.null(cov.corrfact)) {
        cov.corrfact <- if (cov.hubercorr) 'empirical' else 'tau'
    } else if (!cov.corrfact %in% c('tau', 'empirical', 'asympt', 'hybrid'))
	stop(":.vcov.w: cov.corrfact is not in 'tau', 'empirical', 'asympt', 'hybrid'")
    if (is.null(cov.resid)) cov.resid <- 'final'
    else if (!cov.resid %in% c('final','initial', 'trick'))
	stop(":.vcov.w: cov.corrfact is not in 'final','initial', 'trick'")
    if (is.null(cov.xwx)) cov.xwx <- TRUE
    else if (!is.logical(cov.hubercorr))
	stop(':.vcov.w: cov.xwx has to be logical')
    if (is.null(x))  x <- model.matrix(obj)
    if (cov.resid == 'initial') {
        psi <- ctrl$psi
        c.psi <- ctrl$tuning.chi
        if (is.null(psi)) stop('parameter psi is not defined')
        if (!is.numeric(c.psi)) stop('parameter tuning.psi is not numeric')
    } else {
        ## set psi and c.psi
        psi <- ctrl$psi
        c.psi <- if (ctrl$method %in% c('S', 'SD'))
            ctrl$tuning.chi else ctrl$tuning.psi
        if (is.null(psi)) stop('parameter psi is not defined')
        if (!is.numeric(c.psi)) stop('parameter tuning.psi is not numeric')
    }
    if (cov.resid == 'final' && (class(obj)[1] == 'lmrob.S'))
        warning(":.vcov.w: ignoring cov.resid == final since est != final")
    if (is.null(scale)) {
        warning(":.vcov.w: scale missing, using D scale")
        scale <- lmrob..D..fit(obj)$scale
    }
    n <- NROW(x)
    ## --- calculations: matrix part
    ## weighted xtx.inv matrix
    w <- if (cov.xwx) obj$weights else rep(1,n)
    ## use qr-decomposition from lm.wfit (this already includes the weights)
    ## update qr decomposition if it is missing or we don't want the weights
    if (!is.qr(obj$qr) || !cov.xwx) obj$qr <- qr(x * sqrt(w))
    p <- if (is.null(obj$rank)) obj$qr$rank else obj$rank
    cinv <- try( if (is.qr(obj$qr)) tcrossprod(solve(qr.R(obj$qr))) )
    if(inherits(cinv, 'try-error')) cinv <- matrix(NA,p,p)
    ## --- calculation: correction factor
    if (cov.corrfact == 'asympt') { ## asympt correction factor
        if (cov.hubercorr) warning('option hcorr is ignored for cov.corrfact = asympt')
        ## precalculated default values if applicable
        corrfact <-
            if (psi == 'ggw') {
                if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) 1.052619
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) 1.052581
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 1.176464
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 1.176479
                else lmrob.E(psi(r)^2, obj=obj) / lmrob.E(psi(r,1), obj=obj)^2
            } else if (isTRUE(all.equal(c.psi, lmrob.control(psi = psi)$tuning.psi))) {
                switch(psi,
                       bisquare = 1.052632, welsh = 1.052670, optimal = 1.052642,
                       hampel = 1.05265, lwg = 1.052628,
                       stop(':.vcov.w: unsupported psi function'))
            } else lmrob.E(psi(r)^2, obj=obj) / lmrob.E(psi(r,1), obj=obj)^2
        varcorr <- 1
    } else { ## empirical, approx or hybrid correction factor
        if (cov.resid == 'initial') {
            ## if the last estimator was a D or T estimator
            ## then use obj$init$init otherwise use obj$init
            ## that way for SMD we use the S residuals (and S scale)
            ## and for SMDM we use the M residuals (and D scale)
            lobj <-
                if (grepl('[DT]$',ctrl$method)) obj$init$init
                else obj$init
            rstand <- resid(lobj) / lobj$scale
        } else if (cov.resid == 'trick')
            ## residuals are in fact from earlier estimator, use its scale to standardize them
            rstand <- obj$init$resid / obj$init$scale
        else rstand <- obj$resid / scale
        tau <- if (cov.corrfact %in% c('tau', 'hybrid')) { ## added hybrid here
            if (!is.null(obj$tau)) obj$tau
            else if (!is.null(obj$init$tau)) obj$init$tau
            else stop(':.vcov.w: tau not found') }
        else rep(1,n)
       rstand <- rstand / tau
        r.psi <- lmrob.psifun(rstand, c.psi, psi)
        r.psipr <- lmrob.psifun(rstand, c.psi, psi, deriv = 1)
        if (any(is.na(r.psipr))) warning(":.vcov.w: Caution. Some psiprime are NA")
        mpp2 <- (mpp <- mean(r.psipr, na.rm=TRUE))^2
        ## Huber's correction
        hcorr <- if (cov.hubercorr) {
            vpp <- sum((r.psipr - mpp)^2) / n # var(r.psipr, na.rm=TRUE)
            (1+p/n*vpp/mpp2)^2
        } else 1
        ## sample size correction for var(r.psi^2)
        ## use tau if 'tau' correction factor, but only if it is available
        varcorr <- if (cov.corrfact == 'tau' && any(tau != 1))
            1 / mean(tau^2) else n / (n - p) ## changed from 1 / mean(tau)
        ## if hybrid: replace B (= mpp2) by asymptotic value
        if (cov.corrfact == 'hybrid') {
            mpp2 <- if (psi == 'ggw') {
                if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) 0.7598857
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) 0.6817983
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.85, NA)))) 0.4811596
                else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.85, NA)))) 0.411581
                else lmrob.E(psi(r, 1), obj=obj)^2
            } else if (isTRUE(all.equal(c.psi, lmrob.control(psi = psi)$tuning.psi)))
                switch(psi,
                       bisquare = 0.5742327, welsh = 0.5445068, optimal = 0.8598825,
                       hampel = 0.6775217, lqq = 0.6883393,
                       stop(':.vcov.w: unsupported psi function'))
            else lmrob.E(psi(r,1), obj=obj)^2
        }
        corrfact <- mean(r.psi^2)/mpp2 * hcorr
    }
    ## simple sample size correction
    sscorr <- if (cov.dfcorr > 0) {
        if (cov.dfcorr == 2) varcorr ## cov.dfcorr == 2
        else if (cov.dfcorr == 3) mean(w)^2 / (1 - p / sum(w)) ## cov.dfcorr == 3
        else mean(w) * varcorr ## cov.dfcorr == 1
    } else if (cov.dfcorr < 0) mean(w)
    else 1 ## cov.dfcorr == -1 and == 0

    ## scale^2 * a/b2 * Huber's correction * Cinv
    cv <- scale^2 * sscorr * corrfact * cinv
    attr(cv,"weights") <- w
    attr(cv,"scale") <- scale
    attr(cv,"scorr") <- sscorr
    cv
}

.vcov.avar1 <- function(obj, x=obj$x, posdef.meth = c("posdefify","orig"))
{ ## was .vcov.MM
    stopifnot(is.list(ctrl <- obj$control))
    ## works only for MM & SM estimates:
    if (!is.null(ctrl$method) && !ctrl$method %in% c('SM', 'MM'))
        stop('.vcov.avar1() supports only SM or MM estimates')
    ## set psi and chi constants
    psi <- chi <- ctrl$psi
    if (is.null(psi)) stop('parameter psi is not defined')
    stopifnot(is.numeric(c.chi <- ctrl$tuning.chi),
	      is.numeric(c.psi <- ctrl$tuning.psi))

    ## need (r0, r, scale, x, c.psi,c.chi, bb)
    r0 <- obj$init$resid
    r <- resid(obj)
    scale <- obj$scale
    if (is.null(x))  x <- model.matrix(obj)
    bb <- 1/2 ## this is always 1/2 for S estimates by convention
### --- start code from .vcov.MM ---
    ## scaled residuals
    n <- length(r)
    stopifnot(n == length(r0), is.matrix(x), n == nrow(x))
    p <- ncol(x)
    r.s	 <- r / scale # final   scaled residuals
    r0.s <- r0 / scale # initial scaled residuals
    w  <- lmrob.psifun(r.s, cc = c.psi, psi = psi, deriv = 1)
    w0 <- lmrob.chifun(r0.s, cc = c.chi, psi = chi, deriv = 1)
    ## FIXME for multivariate y :
    A <- solve(crossprod(x, x * w)) * (n * scale)
    a <- A %*% (crossprod(x, w * r.s) / (n * mean(w0 * r0.s)))
    w <- lmrob.psifun( r.s, cc = c.psi, psi = psi)

    ## 3) now the standard part  (w, x, r0.s,  n, A,a, c.chi, bb)
    w0 <- lmrob.chifun(r0.s, cc = c.chi, psi = chi)
    Xww <- crossprod(x, w*w0)
    u1 <- crossprod(x, x * w^2) / n
    u1 <- A %*% u1 %*% A
    u2 <- a %*% crossprod(Xww, A) / n
    u3 <- A %*% tcrossprod(Xww, a) / n
    u4 <- mean(w0^2 - bb^2) * tcrossprod(a)

    ## list(cov = matrix((u1 - u2 - u3 + u4)/n, p, p),
    ##      wt = w / r.s, a = a)
### --- end code from .vcov.MM ---
    ret <- (u1 - u2 - u3 + u4)/n

    ## this might not be a positive definite matrix
    ## check eigenvalues (symmetric: ensure non-complex)
    ev <- eigen(ret, symmetric = TRUE)
    if (any(neg.ev <- ev$values < 0)) { ## there's a problem
	posdef.meth <- match.arg(posdef.meth)
	if(ctrl$trace.lev)
	    message("fixing ", sum(neg.ev),
		    " negative eigen([",p,"])values")
	Q <- ev$vectors
	switch(posdef.meth,
	       "orig" = {
		   ## remove negative eigenvalue:
		   ## transform covariance matrix into eigenbasis
		   levinv <- solve(Q)
		   cov.eb <- levinv %*% ret %*% Q
		   ## set vectors corresponding to negative ev to zero
		   cov.eb[, neg.ev] <- 0
		   ## cov.eb[cov.eb < 1e-16] <- 0
		   ## and transform back
		   ret <- Q %*% cov.eb %*% levinv
	       },
	       "posdefify" = {
		   ## Instead of using	require("sfsmisc") and
		   ## ret <- posdefify(ret, "someEVadd",eigen.m = ev,eps.ev = 0)
		   lam <- ev$values
		   lam[neg.ev] <- 0
		   o.diag <- diag(ret)# original one - for rescaling
		   ret <- Q %*% (lam * t(Q)) ## == Q %*% diag(lam) %*% t(Q)
		   ## rescale to the original diagonal values
		   ## D <- sqrt(o.diag/diag(ret))
		   ## where they are >= 0 :
		   D <- sqrt(pmax(0, o.diag)/diag(ret))
		   ret[] <- D * ret * rep(D, each = n) ## == diag(D) %*% m %*% diag(D)
	       },
	       stop("invalid 'posdef.meth': ", posdef.meth))
    }
    attr(ret,"weights") <- w / r.s
    attr(ret,"eigen") <- ev
    ret
}## end{.vcov.avar1}

lmrob..M..fit <- function (x=obj$x, y=obj$y, beta.initial=obj$coef,
                           scale=obj$scale, control=obj$control, obj)
{
    c.psi <- lmrob.conv.cc(control$psi, control$tuning.psi)
    ipsi <- lmrob.psi2ipsi(control$psi)
    stopifnot(is.matrix(x))
    n <- nrow(x)
    p <- ncol(x)
    if (is.null(y) && !is.null(obj$model))
        y <- model.response(obj$model, "numeric")
    stopifnot(length(y) == n,
              length(c.psi) > 0, c.psi >= 0,
              scale >= 0, length(beta.initial) == p)

    ret <- .C(R_lmrob_MM,
              x = as.double(x),
              y = as.double(y),
              n = as.integer(n),
              p = as.integer(p),
              beta.initial = as.double(beta.initial),
              scale = as.double(scale),
              coefficients = double(p),
              residuals = double(n),
              iter = as.integer(control$max.it),
              c.psi = as.double(c.psi),
              ipsi = as.integer(ipsi),
              loss = double(1),
              rel.tol = as.double(control$rel.tol),
              converged = logical(1),
              trace.lev = as.integer(control$trace.lev),
              mts = as.integer(control$mts),
              ss = .convSs(control$subsampling)
              )[c("coefficients",  "scale", "residuals", "loss", "converged", "iter")]
    ## FIXME?: Should rather warn *here* in case of non-convergence
    names(ret$coefficients) <- colnames(x)
    ret$weights <- lmrob.wgtfun(ret$residuals / scale, control$tuning.psi, control$psi)
    ret$fitted.values <- drop(x %*% ret$coefficients)
    if (!grepl('M$', control$method)) {
        ## update control$method if it's not there already
        control$method <- paste(control$method, 'M', sep = '')
    }
    ret$control <- control
    if (!missing(obj)) {
        if (!is.null(obj$call)) {
            ret$call <- obj$call
            ret$call$method <- control$method
        }
        if (control$method %in% c('SM', 'MM')) {
            ret$init.S <- obj
        } else {
            ret$init <-
                obj[names(obj)[na.omit(match(c("coefficients","scale", "residuals",
                                               "loss", "converged", "iter", "weights",
                                               "fitted.values", "control", "init.S", "init",
                                               "kappa", "tau"),
                                             names(obj)))]]
            class(ret$init) <- 'lmrob'
            ret <- c(ret,
                     obj[names(obj)[na.omit(match(c("df.residual", "degree.freedom",
                                                    "xlevels", "terms", "model", "x", "y",
                                                    "na.action", "contrasts", "MD"),
                                                  names(obj)))]])
        }
        ret$qr <- qr(x * sqrt(ret$weights))
        ret$rank <- ret$qr$rank
        ## if there is a covariance matrix estimate available in obj
        ## update it, if possible, else replace it by the default .vcov.w
        if (!is.null(obj$cov)) {
            if (!control$method %in% c('SM', 'MM') &&
                ret$control$cov == '.vcov.avar1') ret$control$cov <- '.vcov.w'

            lf.cov <- if (!is.function(ret$control$cov))
                get(ret$control$cov, mode='function') else ret$control$cov
            ret$cov <- lf.cov(ret, x)
        }
    }
    class(ret) <- "lmrob"
    ret
}


lmrob.S <- function (x, y, control, trace.lev = control$trace.lev, mf = NULL)
{
    if (!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    nResample <- as.integer(control$nResample)
    groups <- as.integer(control$groups)
    nGr <- as.integer(control$n.group)
    if (nGr <= p)
        stop("'control$n.group' must be larger than 'p'")
    large_n <- (n > 2000)
    if (large_n & nGr * groups > n)
        stop("'groups * n.group' must be smaller than 'n' for 'large_n' algorithm")
    if (nGr <= p + 10) ## FIXME (be smarter ..)
        warning("'control$n.group' is not much larger than 'p', probably too small")
    if (length(seed <- control$seed) > 0) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            seed.keep <- get(".Random.seed", envir = .GlobalEnv,
                             inherits = FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
        }
        assign(".Random.seed", seed, envir = .GlobalEnv) ## why not set.seed(seed)
    }

    bb <- as.double(control$bb)
    c.chi <- lmrob.conv.cc(control$psi, control$tuning.chi)
    best.r <- as.integer(control$best.r.s)
    stopifnot(length(c.chi) > 0, c.chi >= 0, length(bb) > 0,
              length(best.r) > 0, best.r >= 1, length(y) == n, n > 0)

    b <- .C(R_lmrob_S,
            x = as.double(x),
            y = as.double(y),
            n = as.integer(n),
            p = as.integer(p),
            nResample = nResample,
            scale = double(1),
            coefficients = double(p),
            as.double(c.chi),
            as.integer(lmrob.psi2ipsi(control$psi)),
            bb,
            best_r = best.r,
            groups = groups,
            n.group = nGr,
            k.fast.s = as.integer(control$k.fast.s),
            k.iter = as.integer(control$k.max),
            refine.tol = as.double(control$refine.tol),
            converged = logical(1),
            trace.lev = as.integer(trace.lev),
            mts = as.integer(control$mts),
            ss = .convSs(control$subsampling)
            )[c("coefficients", "scale", "k.iter", "converged")]
    scale <- b$scale
    if (scale < 0)
        stop("C function R_lmrob_S() exited prematurely")
    class(b) <- 'lmrob.S'
    ## FIXME: get 'res'iduals from C

    names(b$coefficients) <- colnames(x)
    b$fitted.values <- x %*% b$coef
    b$residuals <- drop(y - b$fitted.values)
    ## robustness weights
    b$weights <- lmrob.wgtfun(b$residuals / b$scale, control$tuning.chi, control$psi)
    ## set method argument in control
    control$method <- 'S'
    b$control <- control
    b
}

lmrob..D..fit <- function(obj, x=obj$x, control = obj$control)
{
    if (is.null(control)) stop('lmrob..D..fit: control is missing')
    if (!obj$converged)
        stop('lmrob..D..fit: prior estimator did not converge, stopping')
    if (is.null(x)) x <- model.matrix(obj)
    w <- obj$weights
    if (is.null(w)) stop('lmrob..D..fit: robustness weights undefined')
    if (is.null(obj$residuals)) stop('lmrob..D..fit: residuals undefined')
    r <- obj$residuals
    psi <- control$psi
    if (is.null(psi)) stop('lmrob..D..fit: parameter psi is not defined')
    c.psi <- lmrob.conv.cc(psi, if (control$method %in% c('S', 'SD'))
                           control$tuning.chi else control$tuning.psi)
    if (!is.numeric(c.psi)) stop('lmrob..D..fit: parameter tuning.psi is not numeric')

    obj$init <- obj[names(obj)[na.omit(match(c("coefficients","scale", "residuals",
                                               "loss", "converged", "iter", "weights",
                                               "fitted.values", "control", "init.S", "init"),
                                             names(obj)))]]
    obj$init.S <- NULL

    if (is.null(obj$kappa))
        obj$kappa <- lmrob.kappa(obj, control)
    kappa <- obj$kappa
    if (is.null(obj$tau))
        obj$tau <- lmrob.tau(obj, x, control)
    tau <- obj$tau

    ## get starting value for root search (to keep breakdown point!)
    scale.1 <- sqrt(sum(w * r^2) / kappa / sum(tau^2*w))
    ret <- .C(R_find_D_scale,
              r = as.double(r),
              kappa = as.double(kappa),
              tau = as.double(tau),
              length = as.integer(length(r)),
              scale = as.double(scale.1),
              c = as.double(c.psi),
              ipsi = lmrob.psi2ipsi(psi),
              type = 3L, ## dt1 as only remaining option
              rel.tol = as.double(control$rel.tol),
              k.max = as.integer(control$k.max),
              converged = integer(1))
    obj$scale <- if(ret$converged) ret$scale else NA
    obj$converged <- ret$converged

    if (!grepl('D$', control$method)) {
        ## append "D"  to control$method if it's not there already
        method <- control$method
        if (method == 'MM') method <- 'SM'
        control$method <- paste(method, 'D', sep = '')
    }
    ## update call
    if (!is.null(obj$call)) obj$call$method <- control$method
    obj$control <- control
    class(obj) <- "lmrob"

    ## if there is a covariance matrix estimate available in obj
    ## update it, if possible, else replace it by the default
    ## .vcov.w
    if (!is.null(obj$cov)) {
        if (control$cov == '.vcov.avar1') control$cov <- '.vcov.w'

        lf.cov <- if (!is.function(control$cov))
            get(control$cov, mode='function') else control$cov
        obj$cov <- lf.cov(obj, x)
    }

    obj
}

if(getRversion() > "2.15.0" || as.numeric(R.Version()$`svn rev`) > 59233)
    utils::globalVariables(c("psi", "wgt", "r"), add=TRUE) ## <- lmrob.E( <expr> )

lmrob.kappa <- function(obj, control = obj$control)
{
    if (is.null(control)) stop('control is missing')
    if (control$method %in% c('S', 'SD')) control$tuning.psi <- control$tuning.chi

    fun.min <- function(kappa) lmrob.E(psi(r)*r - kappa*wgt(r), control = control)
    uniroot(fun.min, c(0.1, 1))$root
}

lmrob.tau <- function(obj,x=obj$x, control = obj$control, h, fast = TRUE)
{
    if (is.null(control)) stop('control is missing')
    if (missing(h)) {
        if (is.null(obj$qr))
            h <- lmrob.leverages(x, obj$weights)
        else
            h <- lmrob.leverages(x, obj$weights, wqr = obj$qr)
    }

    ## speed up: use approximation if possible
    if (fast && !control$method %in% c('S', 'SD')) {
        c.psi <- control$tuning.psi
        tfact <- tcorr <- NA
        switch(control$psi,
               optimal = if (isTRUE(all.equal(c.psi, 1.060158))) {
                   tfact <- 0.94735878
                   tcorr <- -0.09444537
               },
               bisquare = if (isTRUE(all.equal(c.psi, 4.685061))) {
                   tfact <- 0.9473684
                   tcorr <- -0.0900833
               },
               welsh = if (isTRUE(all.equal(c.psi, 2.11))) {
                   tfact <- 0.94732953
                   tcorr <- -0.07569506
               },
               ggw = if (isTRUE(all.equal(c.psi, c(-.5, 1.0, 0.95, NA)))) {
                   tfact <- 0.9473787
                   tcorr <- -0.1143846
               } else if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) {
                   tfact <- 0.94741036
                   tcorr <- -0.08424648
               },
               lqq = if (isTRUE(all.equal(c.psi, c(-.5, 1.5, 0.95, NA)))) {
                   tfact <- 0.94736359
                   tcorr <- -0.08594805
               },
               hampel = if (isTRUE(all.equal(c.psi, c(1.35241275, 3.15562975, 7.212868)))) {
                   tfact <- 0.94739770
                   tcorr <- -0.04103958
               },
           {})
        if (!is.na(tfact))
            return(sqrt(1 - tfact*h) * (tcorr*h + 1))
    }

    ## kappa
    if (is.null(obj$kappa))
        obj$kappa <- lmrob.kappa(obj, control)
    kappa <- obj$kappa
    ## local variables
    n <- length(h)
    ## set psi and c.psi
    psi <- control$psi
    if (is.null(psi)) stop('parameter psi is not defined')
    c.psi <- if (control$method %in% c('S', 'SD'))
        control$tuning.chi else control$tuning.psi
    if (!is.numeric(c.psi)) stop('parameter tuning.psi is not numeric')

    ## constant for stderr of u_{-i} part and other constants
    inta <- function(r)
        lmrob.psifun(r, c.psi, psi)^2*dnorm(r)
    ta <- (integrate(inta,-Inf,Inf))$value
    intb <- function(r)
        lmrob.psifun(r, c.psi, psi, deriv = 1)*dnorm(r)
    tb <- (integrate(intb,-Inf,Inf))$value
    intc <- function(r)
        lmrob.psifun(r, c.psi, psi)*r*dnorm(r) ## changed from psi/e to psi*e
    tE <- (integrate(intc,-Inf,Inf))$value

    ## calculate tau for unique h
    hu <- unique(h)
    nu <- length(hu)

    ## Initialize tau vector
    tau <- numeric(length=nu)

    tc <- ta/tb^2
    ## --- Gauss-Hermite integration
    gh <- ghq(control$numpoints)
    ghz <- gh$nodes
    ghw <- gh$weights
    ## Calulate each tau_i
    for (i in 1:nu) {
        ## stderr of u_{-i} part
        s <- sqrt(tc*(hu[i]-hu[i]^2))
        tc2 <- hu[i]/tb
        ## function to be integrated
        fun <- function(w, v, sigma.i) {
            t <- (v-tc2*lmrob.psifun(v,c.psi,psi)+w*s)/sigma.i
            (lmrob.psifun(t,c.psi,psi)*t - kappa*lmrob.psifun(t, c.psi, psi)/t)*
                dnorm(v)*dnorm(w)
        }
        ## integrate over w
        wint <- function(v, sigma.i) {
            ## sapply(v,function(v.j) integrate(fun,-Inf,Inf,v.j,sigma.i)$value)
            sapply(v, function(v.j) sum(fun(ghz, v.j, sigma.i)*ghw))
        }
        ## integrate over v
        vint <- function(sigma.i) {
            ## integrate(wint,-Inf,Inf,sigma.i)$value
            sum(wint(ghz, sigma.i)*ghw)
        }

        ## find tau
        tau[i] <- uniroot(vint, c(if (hu[i] < 0.9) 3/20 else 1/16, 1.1))$root
    }

    tau[match(h, hu)]
}

lmrob.tau.fast.coefs <- function(cc, psi) {
    ## function that calculates the coefficients for 'fast' mode of lmrob.tau
    ctrl <- lmrob.control(tuning.psi = cc, psi = psi)
    levs <- seq(0, 0.8, length.out = 80)
    ## calculate taus
    taus <- lmrob.tau(list(),,ctrl,h=levs,fast=FALSE)
    ## calculate asymptotic approximation of taus
    ta <- lmrob.E(psi(r)^2, ctrl, use.integrate = TRUE)
    tb <- lmrob.E(psi(r, 1), ctrl, use.integrate = TRUE)
    tE <- lmrob.E(psi(r)*r, ctrl, use.integrate = TRUE)
    tfact <- 2*tE/tb - ta/tb^2
    taus.0 <- sqrt(1 - tfact * levs)
    ## calculate correction factor
    tcorr <- coef(lmrob(taus / taus.0 - 1 ~ levs - 1))
    c(tfact = tfact, tcorr = tcorr)
}

lmrob.hatmatrix <- function(x, w = rep(1, NROW(x)), wqr = qr(sqrt(w) * x))
{
    tcrossprod(qr.Q(wqr))
}

lmrob.leverages <- function(x, w = rep(1, NROW(x)), ...)
{
    if (!is.matrix(x)) x <- as.matrix(x)

    diag(lmrob.hatmatrix(x, w, ...))
}

lmrob.psi2ipsi <- function(psi)
{
    switch(casefold(psi),
           'tukey' = , 'biweight' = , 'bisquare' = 1L,
           'welsh' = 2L,
           'optimal' = 3L,
           'hampel' = 4L,
           'ggw' = 5L,
           'lqq' = 6L,
           stop("unknown psi function ", psi))
}

lmrob.conv.cc <- function(psi, cc)
{
    if (!is.character(psi) || length(psi) != 1)
        stop("argument 'psi' must be a string (denoting a psi function)")
    if(!is.numeric(cc))
        stop("tuning constant 'cc' is not numeric")

    switch(casefold(psi),
           'ggw' = {
               ## Input: 4 parameters, (minimal slope, b, efficiency, breakdown point)
               ## Output 'k': either k in {1:6} or  k = c(0, k[2:5])
               if (isTRUE(all.equal(cc, c(-.5, 1, 0.95, NA)))) return(1)
               else if (isTRUE(all.equal(cc, c(-.5, 1, 0.85, NA)))) return(2)
               else if (isTRUE(all.equal(cc, c(-.5, 1.0, NA, 0.5)))) return(3)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA)))) return(4)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA)))) return(5)
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5)))) return(6)
               else if (length(cc) == 5 && cc[1] == 0 ||
                        (length(cc <- attr(cc, 'constants')) == 5 && cc[1] == 0))
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
                         'Use c({0, } minimal slope, b, efficiency, breakdown point)')
           },
           'lqq' = {
               ## Input: 4 parameters, (minimal slope, b, efficiency, breakdown point)
               ## Output: k[1:3]
               if (isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))))
                   return(c(1.4734061, 0.9822707, 1.5))
               else if (isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))
                   return(c(0.4015457, 0.2676971, 1.5))
               else if (length(cc) == 3 || length(cc <- attr(cc, 'constants')) == 3)
                   return(cc)
               else stop('Coefficients for ',psi,' function incorrectly specified.\n',
                         'Use c(minimal slope, b, efficiency, breakdown point) or k[1:3]')
           },
           'hampel' = {
               ## just check length of coefficients
               if (length(cc) != 3)
                   stop('Coef. for Hampel psi function not of length 3')
           }, {
               ## otherwise: should have length 1
               if (length(cc) != 1)
                   stop('Coef. for psi function ', psi,' not of length 1')
           })

    return(cc)
}

lmrob.ggw.mx <- function(a, b, c) ## find x with minimal slope
    optimize(lmrob.psifun, c(c, max(a+b+2*c, 0.5)),
             cc=c(0, a, b, c, 1), psi = 'ggw', deriv = 1)[[1]]

lmrob.ggw.ms <- function(a, b, c) ## find minimal slope
    lmrob.psifun(lmrob.ggw.mx(a, b, c), c(0, a, b, c, 1), 'ggw', 1)

lmrob.ggw.finda <- function(ms, b, c) ## find a constant
{
    val <- uniroot(function(a) lmrob.ggw.ms(1/a, b, c) - ms,
                   c(200, if (b > 1.4) 1/400 else if (b > 1.3) 1/50 else 1/20))
    1/val$root
}

lmrob.ggw.ac <- function(a, b, c) ## calculate asymptotic efficiency
{
    lmrob.E(lmrob.psifun(r, c(0, a, b, c, 1), 'ggw', 1), use.integrate = TRUE)^2 /
        lmrob.E(lmrob.psifun(r, c(0, a, b, c, 1), 'ggw')^2, use.integrate = TRUE)
}

lmrob.ggw.bp <- function(a, b, c) { ## calculate kappa
    nc <- integrate(function(x) lmrob.psifun(x, c(0, a, b, c, 1), 'ggw'), 0, Inf)$value
    lmrob.E(lmrob.chifun(r, c(0, a, b, c, nc), 'ggw'), use.integrate = TRUE)
}

lmrob.ggw.findc <- function(ms, b, eff = NA, bp = NA) {
    ## find c by eff for bp
    c <- if (!is.na(eff)) {
        if (!is.na(bp))
            warning('tuning constants for ggw psi: both eff and bp specified, ignoring bp')
        ## find c by eff
        uniroot(function(x) lmrob.ggw.ac(lmrob.ggw.finda(ms, b, x), b, x) - eff,
                c(0.15, if (b > 1.61) 1.4 else 1.9))$root
    } else {
        if (is.na(bp))
            stop('Error: neither breakdown point nor efficiency specified')
        ## find c by bp
        uniroot(function(x) lmrob.ggw.bp(lmrob.ggw.finda(ms, b, x), b, x) - bp,
                c(0.08, if (ms < -0.4) 0.6 else 0.4))$root
    }

    a <- lmrob.ggw.finda(ms, b, c)
    nc <- integrate(function(x) lmrob.psifun(x, c(0, a, b, c, 1), 'ggw'), 0, Inf)$value
    return(c(0, a, b, c, nc))
}

lmrob.efficiency <-  function(psi, cc) {
  integrate(function(x) lmrob.psifun(x, cc, psi, 1)*dnorm(x), -Inf, Inf)$value^2 /
    integrate(function(x) lmrob.psifun(x, cc, psi)^2*dnorm(x), -Inf, Inf)$value
}

lmrob.bp <- function(psi, cc)
  integrate(function(x) lmrob.chifun(x, cc, psi)*dnorm(x), -Inf, Inf)$value

lmrob.lqq.findc <- function(cc) {
    ## cc = c(min slope, b, eff, bp)
    ## constants for c function: c(b*c, c, s = 1 - min slope)
    t.fun <- if (!is.na(cc[3])) {
        if (!is.na(cc[4]))
            warning('tuning constants for lqq psi: both eff and bp specified, ignoring bp')
        ## find c by b, s and eff
        function(c)
            lmrob.efficiency('lqq', c(cc[2]*c, c, 1-cc[1])) - cc[3]
    } else {
        if (is.na(cc[4]))
            stop('Error: neither breakdown point nor efficiency specified')
        function(c)
            lmrob.bp('lqq', c(cc[2]*c, c, 1-cc[1])) - cc[4]
    }
    c <- try(uniroot(t.fun, c(0.1, 4))$root, silent = TRUE)

    if (inherits(c, 'try-error'))
        stop('lmrob.lqq.findc: unable to find constants for psi function')
    else return(c(cc[2]*c, c, 1-cc[1]))
}

lmrob.const <- function(cc, psi)
{
    switch(psi,
           "ggw" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1,   0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1,   NA,  0.5))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, 0.85, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
                   attr(cc, 'constants') <- lmrob.ggw.findc(cc[1],cc[2],cc[3],cc[4])
               }
           },
           "lqq" = { ## only calculate for non-standard coefficients
               if (!(isTRUE(all.equal(cc, c(-.5, 1.5, 0.95, NA))) ||
                     isTRUE(all.equal(cc, c(-.5, 1.5, NA, 0.5))))) {
                   attr(cc, 'constants') <- lmrob.lqq.findc(cc)
               }
           },
           stop("method for psi function ",psi, " not implemented"))

    return(cc)
}

lmrob.psifun <- function(x, cc, psi, deriv=0)
{
    cc <- lmrob.conv.cc(psi, cc)

    ## catch NAs
    idx <- !is.na(x)

    if (any(idx))
        x[idx] <- .C(R_psifun, x = as.double(x[idx]), cc = as.double(cc),
                     ipsi = as.integer(lmrob.psi2ipsi(psi)), NAOK= TRUE, # for +- Inf
                     deriv = as.integer(deriv), length = as.integer(length(x[idx])))$x
    x
}

lmrob.chifun <- function(x, cc, psi, deriv=0)
{
    cc <- lmrob.conv.cc(psi, cc)

    ## catch NAs
    idx <- !is.na(x)

    if (any(idx))
        x[idx] <- .C(R_chifun, x = as.double(x[idx]), cc = as.double(cc),
                     ipsi = as.integer(lmrob.psi2ipsi(psi)), NAOK= TRUE, # for +- Inf
                     deriv = as.integer(deriv), length = as.integer(length(x[idx])))$x
    x
}

lmrob.wgtfun <- function(x, cc, psi)
{
    cc <- lmrob.conv.cc(psi, cc)

    ## catch NAs
    idx <- !is.na(x)

    if (any(idx))
        x[idx] <- .C(R_wgtfun, x = as.double(x[idx]), cc = as.double(cc),
                     ipsi = as.integer(lmrob.psi2ipsi(psi)), NAOK= TRUE, # for +- Inf
                     length = as.integer(length(x[idx])))$x
    x
}

residuals.lmrob.S <- function(obj)
    obj$residuals

lmrob.E <- function(expr, control, dfun = dnorm, use.integrate = FALSE, obj)
{
    expr <- substitute(expr)
    if (missing(control) && !missing(obj))
        control <- obj$control

    lenvir <-
      if (!missing(control)) {
        psi <- control$psi
        if (is.null(psi)) stop('parameter psi is not defined')
	c.psi <- control[[if (control$method %in% c('S', 'SD'))
			  "tuning.chi" else "tuning.psi"]]
        if (!is.numeric(c.psi)) stop('tuning parameter (chi/psi) is not numeric')

        list(psi = function(r, deriv = 0) lmrob.psifun(r, c.psi, psi, deriv),
             chi = function(r, deriv = 0) lmrob.chifun(r, c.psi, psi, deriv), ## change?
             wgt = function(r) lmrob.wgtfun(r, c.psi, psi)) ## change?

      } else list()

    if (use.integrate) {
        pf <- parent.frame()
        integrate(function(r)
                  eval(expr, envir = c(list(r = r), lenvir),
                       enclos = pf)*dfun(r),-Inf,Inf)$value
    } else {
        ## initialize Gauss-Hermite Integration
        gh <- ghq(if (is.null(control$numpoints)) 13 else control$numpoints)
        ghz <- gh$nodes
        ## integrate
	sum(gh$weights * eval(expr, envir = c(list(r = ghz), lenvir),
			      enclos = parent.frame()) * dfun(ghz))
    }
}

ghq <- function(n = 1, modify = TRUE) {
    ## Adapted from gauss.quad in statmod package
    ## which itself has been adapted from Netlib routine gaussq.f
    ## Gordon Smyth, Walter and Eliza Hall Institute

    n <- as.integer(n)
    if(n<0) stop("need non-negative number of nodes")
    if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
    ## i <- seq_len(n) # 1 .. n
    i1 <- seq_len(n-1L)

    muzero <- sqrt(pi)
    ## a <- numeric(n)
    b <- sqrt(i1/2)

    A <- numeric(n*n)
    ## A[(n+1)*(i-1)+1] <- a # already 0
    A[(n+1)*(i1-1)+2] <- b
    A[(n+1)*i1] <- b
    dim(A) <- c(n,n)
    vd <- eigen(A,symmetric=TRUE)
    n..1 <- n:1L
    w <- vd$vectors[1, n..1]
    w <- muzero * w^2
    x <- vd$values[n..1] # = rev(..)
    list(nodes=x, weights= if (modify) w*exp(x^2) else w)
}

.convSs <- function(ss)
    switch(ss,
           simple=0L,
           constrained=1L,
           stop("unknown setting for parameter ss"))
           
