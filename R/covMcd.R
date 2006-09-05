#### This is originally from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## hidden in namespace:
quan.f <- function(alpha, n, rk) {
    ## Compute size of subsample, given alpha
    ## Same function for covMcd() and ltsReg()
    n2 <- (n+rk+1) %/% 2
    floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

covMcd <- function(x,
		   cor = FALSE,
		   alpha = 1/2,
		   nsamp = 500,
		   seed = NULL,
		   trace = FALSE,
		   use.correction = TRUE,
		   control)
{
    ##	 Analyze and validate the input parameters ...

    ## if a control object was supplied, take the option parameters
    ## from it, but if single parameters were passed (not defaults)
    ## they will override the control object.
    if(!missing(control)) {
	defCtrl <- rrcov.control()	# default control
	if(alpha == defCtrl$alpha)	 alpha <- control$alpha
	if(nsamp == defCtrl$nsamp)	 nsamp <- control$nsamp
	if(identical(seed, defCtrl$seed)) seed <- control$seed
	if(trace == defCtrl$trace) trace <- control$trace
	if(use.correction == defCtrl$use.correction)
	    use.correction <- control$use.correction
    }

    if(length(seed) > 0) {
	if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
	    seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
	    on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
	}
	assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    ##	 vt::03.02.2006 - added options "best" and "exact" for nsamp
    ##	 nsamp will be further analized in the wrapper .fastmcd()
    if(!missing(nsamp) && is.numeric(nsamp) && nsamp <= 0)
	stop("Invalid number of trials nsamp = ",nsamp, "!")


    ## vt:: tolerance to be used for computing the mahalanobis distances (default = 1e-7)
    tol <- 1e-10

    if(is.data.frame(x))
	x <- data.matrix(x)
    else if (!is.matrix(x))
	x <- matrix(x, length(x), 1,
		    dimnames = list(names(x), deparse(substitute(x))))

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
	stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]
    if(n < 2 * p)
	stop("Need at least 2*(number of variables) observations ")
    jmin <- (n + p + 1) %/% 2
    if(alpha < 1/2) ## FIXME? shouldn't we rather test	'alpha < jmin/n' ?
	stop("The MCD must cover at least", jmin, "observations")
    else if(alpha > 1)
	stop("alpha is out of range")

    quantiel <- qchisq(0.975, p)
    quan <- quan.f(alpha, n, p)

    ## vt::03.02.2006 - raw.cnp2 and cnp2 are vectors of size 2 and  will
    ##	 contain the correction factors (concistency and finite sample)
    ##	 for the raw and reweighted estimates respectively. Set them
    ##	 initially to 1.  If use.correction is set to FALSE
    ##	 (default=TRUE), the finite sample correction factor will not
    ##	 be used (neither for the raw estimates nor for the reweighted)
    raw.cnp2 <- cnp2 <- c(1,1)

    ans <- list(method = "Minimum Covariance Determinant Estimator.",
		call = match.call())

    if(alpha == 1) { ## alpha=1: Just compute the classical estimates --------
	mcd <- cov.wt(x)$cov
	loc <- as.vector(colMeans(x))
	obj <- determinant(mcd, log = TRUE)$modulus[1]
	if ( -obj/p > 50 ) {
	    ans$cov <- mcd
	    dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
	    if (cor)
		ans$cor <- cov2cor(ans$cov)
	    ans$center <- loc
	    if(length(dimn[[2]]))
		names(ans$center) <- dimn[[2]]
	    ans$n.obs <- n
	    msg <- "The classical covariance matrix is singular."
	    ans$method <- paste(ans$method, msg, sep="\n")
	    if(trace)
		cat(msg,"\n")

	    weights <- 1
	}
	else {
	    mah <- mahalanobis(x, loc, mcd, tol = tol)
	    ## VT:: 01.09.2004 - bug in alpha=1
	    ##	    (tol instead of tol.inv as parameter name)
	    weights <- as.numeric(mah < quantiel) # 0/1
	    sum.w <- sum(weights)
	    ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
	    ans$cov <- sum.w/(sum.w - 1) * ans$cov

	    ## Consistency factor for reweighted MCD
	    if(sum.w != n) {
		qdelta.rew <- qchisq(sum.w/n, p)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1) / (sum.w/n)
		cnp2[1] <- 1/cdeltainvers.rew
		ans$cov <- ans$cov * cnp2[1]
	    }
	    if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
		msg <- "The reweighted MCD scatter matrix is singular."
		ans$method <- paste(ans$method, msg, sep="\n")
		if(trace)
		    cat(msg,"\n")
	    }
	    else {
		mah <- mahalanobis(x, ans$center, ans$cov, tol = tol)
		weights <- as.numeric(mah < quantiel) # 0/1
	    }
	}

	ans$alpha <- alpha
	ans$quan <- quan
	ans$raw.cov <- mcd
	dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
	ans$raw.center <- loc
	if(length(dimn[[2]]))
	    names(ans$raw.center) <- dimn[[2]]
	ans$crit <- exp(obj)
	ans$method <- paste(ans$method,
			    "\nThe minimum covariance determinant estimates based on",
			    n, "observations \nare equal to the classical estimates.")
	ans$mcd.wt <- rep(NA, length(ok))
	ans$mcd.wt[ok] <- weights
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	ans$wt <- NULL
	ans$X <- x
	if(length(dimn[[1]]))
	    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
	else
	    dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
	if(trace)
	    cat(ans$method, "\n")
	ans$raw.cnp2 <- raw.cnp2
	ans$cnp2 <- cnp2
	class(ans) <- "mcd"
	return(ans)
    } ## end {alpha=1} --

    mcd <- .fastmcd(x, quan, nsamp)

    ## Compute the consistency correction factor for the raw MCD
    ##	(see calfa in Croux and Haesbroeck)
    qalpha <- qchisq(quan/n, p)
    calphainvers <- pgamma(qalpha/2, p/2 + 1)/(quan/n)
    raw.cnp2[1] <- calpha <- 1/calphainvers
    raw.cnp2[2] <- correct <- MCDcnp2(p, n, alpha)
    if(!use.correction)	  # do not use finite sample correction factor
	raw.cnp2[2] <- correct <- 1.0

    if(p == 1) {
	## ==> Compute univariate location and scale estimates
	ans$method <- "Univariate location and scale estimation."

	scale <- sqrt(calpha * correct) * as.double(mcd$initcovariance)
	center <- as.double(mcd$initmean)
	if(abs(scale - 0) < 1e-07) {
	    ## VT:: 22.12.04 - ans$cov and ans$raw.cov must be a matrices
	    ans$cov <- matrix(0)
	    names(ans$cov) <- dimn[[2]][1]
	    ans$center <- center
	    names(ans$center) <- dimn[[2]][1]
	    ans$n.obs <- n
	    ans$method <- paste(ans$method,"\nMore than", quan,
				"of the observations are identical.")
	    ans$alpha <- alpha
	    ans$quan <- quan
	    ans$raw.cov <- matrix(0)
	    names(ans$raw.cov) <- dimn[[2]][1]
	    ans$raw.center <- center
	    names(ans$raw.center) <- dimn[[2]][1]
	    ans$crit <- 0
	    weights <- as.numeric(abs(x - center) < 1e-07) # 0 / 1
	} ## end { scale ~= 0 }
	else {
	    ## Compute the weights for the raw MCD in case p=1
	    weights <- as.numeric(((x - center)/scale)^2 < quantiel) # 0/1
	    sum.w <- sum(weights)
	    ans <- c(ans, cov.wt(x, wt = weights, cor = cor))
	    ans$cov <- sum.w/(sum.w - 1) * ans$cov

	    ## Apply the correction factor for the reweighted cov
	    if(sum.w == n) {
		cdelta.rew <- 1
		correct.rew <- 1
	    }
	    else {
		qdelta.rew <- qchisq(sum.w/n, p)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum.w/n)
		cnp2[1] <- cdelta.rew <- 1/cdeltainvers.rew
		cnp2[2] <- correct.rew <- MCDcnp2.rew(p, n, alpha)
		if(!use.correction) # do not use finite sample correction factor
		    cnp2[2] <- correct.rew <- 1.0
	    }
	    ans$cov <- ans$cov * cdelta.rew * correct.rew
	    ans$alpha <- alpha
	    ans$quan <- quan
	    ans$raw.cov <- scale^2
	    names(ans$raw.cov) <- dimn[[2]][1]
	    ans$raw.center <- as.vector(center)
	    names(ans$raw.center) <- dimn[[2]][1]
	    ans$crit <- 1/(quan - 1) *
                sum(sort((x - as.double(mcd$initmean))^2, partial = quan)[1:quan])
	    center <- ans$center
	    scale <- as.vector(sqrt(ans$cov))
	    weights <- as.numeric(((x - center)/scale)^2 < quantiel)
	} ## end{ scale > 0 }
    } ## end p=1

    else { ## p >= 2 : ----------------------------------------------------------

	## Apply correction factor to the raw estimates and use them to compute weights
	mcd$initcovariance <- calpha * mcd$initcovariance * correct
	dim(mcd$initcovariance) <- c(p, p)

	## If not all observations are in general position, i.e. more than
	## h observations lie on a hyperplane, the program still yields
	## the MCD location and scatter matrix, the latter being singular
	## (as it should be), as well as the equation of the hyperplane.
	if(mcd$exactfit != 0) {
	    dim(mcd$coeff) <- c(5, p)
	    ans$cov <- mcd$initcovariance
	    dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
	    ans$center <- as.vector(mcd$initmean)
	    if(length(dimn[[2]]))
		names(ans$center) <- dimn[[2]]
	    ans$n.obs <- n
## no longer relevant:
##	    if(mcd$exactfit == -1)
##		stop("The program allows for at most ", mcd$kount, " observations.")
##	    if(mcd$exactfit == -2)
##		stop("The program allows for at most ", mcd$kount, " variables.")
	    if(mcd$exactfit == 1) {
		msg <- "The covariance matrix of the data is singular."
		ans$method <- paste(ans$method, msg, sep = "\n")
		if(trace)
		    cat(msg, "\n")
	    }
	    if(mcd$exactfit == 2) {
		msg <- paste("The covariance matrix has become singular during",
			     "the iterations of the MCD algorithm.",
			     collapse = "\n")
		ans$method <- paste(ans$method, msg, sep = "\n")
		if(trace)
		    cat(msg, "\n")
	    }
	    if(p == 2) {
		msg <- paste("There are", mcd$kount,
			     "observations in the entire dataset of\n", n,
			     "observations that lie on the line with equation\n",
			     signif(mcd$coeff[1,1], digits= 5), "(x_i1-m_1) +",
			     signif(mcd$coeff[1,2], digits= 5), "(x_i2-m_2)=0\n",
			     "with (m_1,m_2) the mean of these observations.")
		ans$method <- paste(ans$method, msg, sep = "\n")
		if(trace)
		    cat(msg, "\n")
	    }
	    else if(p == 3) {
		msg <- paste("There are", mcd$kount,
			     "observations in the entire dataset of\n", n,
			     "observations that lie on the plane with equation\n",
			     signif(mcd$coeff[1,1], digits= 5), "(x_i1-m_1) +",
			     signif(mcd$coeff[1,2], digits= 5), "(x_i2-m_2) +",
			     signif(mcd$coeff[1,3], digits= 5), "(x_i3-m_3)=0\n",
			     "with (m_1,m_2) the mean of these observations."
			     )
		ans$method <- paste(ans$method, msg, sep = "\n")
		if(trace)
		    cat(msg, "\n")
	    }
	    else { ##  p > 3 -----------
		msg1 <- paste("There are", mcd$kount,
			     "observations in the entire dataset of\n", n,
			     "observations that lie on the hyperplane with equation\n",
			     "a_1*(x_i1-m_1)+...+a_p*(x_ip-m_p)=0 \n",
			     "with (m_1,...,m_p) the mean\n",
			     "of these observations and coefficients a_i equal to: \n")
                msg <- paste(msg1,
                             paste(formatC(mcd$coeff[1, ], digits= 5), collapse=","))
		if(trace) {
		    cat(msg1)
		    print(signif(mcd$coeff[1, ], digits= 5))
		}
                ans$method <- paste(ans$method, msg, sep = "\n")
	    } ## end {p > 3}

	    ans$alpha <- alpha
	    ans$quan <- quan
	    ans$raw.cov <- mcd$initcovariance
	    dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
	    ans$raw.center <- as.vector(mcd$initmean)
	    if(length(dimn[[2]]))
		names(ans$raw.center) <- dimn[[2]]
	    ans$crit <- 0
	    weights <- mcd$weights
	} ## end exact fit <==>	 (mcd$exactfit != 0)

	else { ## exactfit == 0 : have general position ------------------------

	    mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tol)
	    mcd$weights <- weights <- as.numeric(mah < quantiel)
	    sum.w <- sum(weights)

	    ## Compute and apply the consistency correction factor for
	    ## the reweighted cov
	    if(sum.w == n) {
		cdelta.rew <- 1
		correct.rew <- 1
	    }
	    else {
		qdelta.rew <- qchisq(sum.w/n, p)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum.w/n)
		cnp2[1] <- cdelta.rew <- 1/cdeltainvers.rew
		cnp2[2] <- correct.rew <- MCDcnp2.rew(p, n, alpha)
		if(!use.correction) # do not use finite sample correction factor
		    cnp2[2] <- correct.rew <- 1.0
	    }

	    ans <- c(ans, cov.wt(x, wt = weights, cor))
	    ans$cov <- sum.w/(sum.w - 1) * ans$cov
	    ans$cov <- ans$cov * cdelta.rew * correct.rew

	    ##vt:: add also the best found subsample to the result list
	    ans$best <- sort(as.vector(mcd$best))

	    ans$alpha <- alpha
	    ans$quan <- quan
	    ans$raw.cov <- mcd$initcovariance
	    dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
	    ans$raw.center <- as.vector(mcd$initmean)
	    if(length(dimn[[2]]))
		names(ans$raw.center) <- dimn[[2]]
	    ans$raw.weights <- weights
	    ans$crit <- mcd$mcdestimate
	    ans$raw.mah <- mahalanobis(x, ans$raw.center, ans$raw.cov, tol = tol)

	    ## Check if the reweighted scatter matrix is singular.
	    if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {

		msg <- "The reweighted MCD scatter matrix is singular."
		ans$method <- paste(ans$method, msg, sep="\n")
		if(trace)
		    cat(msg,"\n")
		ans$mah <- ans$raw.mah
	    }
	    else {
		mah <- mahalanobis(x, ans$center, ans$cov, tol = tol)
		ans$mah <- mah
		weights <- as.numeric(mah < quantiel)
	    }
	} ## end{ not exact fit }

    } ## end{ p >= 2 }

    ans$mcd.wt <- rep(NA, length(ok))
    ans$mcd.wt[ok] <- weights
    if(length(dimn[[1]]))
	names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
	dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else
	dimnames(ans$X) <- list(seq(along = ok)[ok], NULL)
    if(trace)
	cat(ans$method, "\n")
    ans$raw.cnp2 <- raw.cnp2
    ans$cnp2 <- cnp2
    class(ans) <- "mcd"
    return(ans)
}

print.mcd <- function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    cat("Minimum Covariance Determinant (MCD) estimator.\n")
    if(!is.null(cl <- x$call)) {
	cat("Call:\n")
	dput(cl)
    }
    cat("-> Method: ", x$method, "\n")
    cat("\nLog(Det.): ", format(log(x$crit), digits = digits) ,"\n")
    cat("Robust Estimate of Location:\n")
    print(x$center, digits = digits, print.gap = print.gap, ...)
    cat("Robust Estimate of Covariance:\n")
    print(x$cov, digits = digits, print.gap = print.gap, ...)
    invisible(x)
}

summary.mcd <- function(object, ...)
{
    class(object) <- c("summary.mcd", class(object))
    object
}

print.summary.mcd <-
    function(x, digits = max(3, getOption("digits") - 3), print.gap = 2, ...)
{
    print.mcd(x, digits = digits, print.gap = print.gap, ...) # see above

    ## hmm, maybe not *such* a good idea :
    if(!is.null(x$cor)) {
      cat("\nRobust Estimate of Correlation: \n")
      dimnames(x$cor) <- dimnames(x$cov)
      print(x$cor, digits = digits, print.gap = print.gap, ...)
    }

    cat("\nEigenvalues:\n")
    print(eigen(x$cov, only.values = TRUE)$values, digits = digits, ...)

    if(!is.null(x$mah)) {
	cat("\nRobust Distances: \n")
	print(as.vector(x$mah), digits = digits, ...)
    }
    invisible(x)
}

## NOTE:  plot.mcd() is in ./covPlot.R !
## ----                    ~~~~~~~~~~~

### --- Namespace hidden (but parsed once and for all) : -------------

MCDcnp2 <- function(p, n, alpha)
{
    if(p > 2) {
	##				"alfaq"		   "betaq"    "qwaarden"
	coeffqpkwad875 <- matrix(c(-0.455179464070565, 1.11192541278794, 2,
				   -0.294241208320834, 1.09649329149811, 3), ncol = 2)
	coeffqpkwad500 <- matrix(c(-1.42764571687802,  1.26263336932151, 2,
				   -1.06141115981725,  1.28907991440387, 3), ncol = 2)

	y.500 <- log( - coeffqpkwad500[1, ] / p^coeffqpkwad500[2, ] )
	y.875 <- log( - coeffqpkwad875[1, ] / p^coeffqpkwad875[2, ] )

	A.500 <- cbind(1, - log(coeffqpkwad500[3, ] * p^2))
	A.875 <- cbind(1, - log(coeffqpkwad875[3, ] * p^2))
	coeffic.500 <- solve(A.500, y.500)
	coeffic.875 <- solve(A.875, y.875)
	fp.500.n <- 1 - exp(coeffic.500[1]) / n^coeffic.500[2]
	fp.875.n <- 1 - exp(coeffic.875[1]) / n^coeffic.875[2]
    }
    else { ## p <= 2
	if(p == 2) {
	    fp.500.n <- 1 - exp( 0.673292623522027) / n^0.691365864961895
	    fp.875.n <- 1 - exp( 0.446537815635445) / n^1.06690782995919
	}
	if(p == 1) {
	    fp.500.n <- 1 - exp( 0.262024211897096) / n^0.604756680630497
	    fp.875.n <- 1 - exp(-0.351584646688712) / n^1.01646567502486
	}
    }

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
	fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
	fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)

    return(1/fp.alpha.n)
} ## end{ MCDcnp2 }

MCDcnp2.rew <- function(p, n, alpha)
{
    if(p > 2) {
	##				"alfaq"		   "betaq"    "qwaarden"
	coeffrewqpkwad875 <- matrix(c(-0.544482443573914, 1.25994483222292, 2,
				      -0.343791072183285, 1.25159004257133, 3), ncol = 2)
	coeffrewqpkwad500 <- matrix(c(-1.02842572724793,  1.67659883081926, 2,
				      -0.26800273450853,  1.35968562893582, 3), ncol = 2)

	y.500 <- log( - coeffrewqpkwad500[1, ] / p^ coeffrewqpkwad500[2, ] )
	y.875 <- log( - coeffrewqpkwad875[1, ] / p^ coeffrewqpkwad875[2, ] )

	A.500 <- cbind(1, - log(coeffrewqpkwad500[3, ] * p^2))
	coeffic.500 <- solve(A.500, y.500)
	A.875 <- cbind(1, - log(coeffrewqpkwad875[3, ] * p^2))
	coeffic.875 <- solve(A.875, y.875)
	fp.500.n <- 1 - exp(coeffic.500[1]) / n^ coeffic.500[2]
	fp.875.n <- 1 - exp(coeffic.875[1]) / n^ coeffic.875[2]
    }
    else {
	if(p == 2) {
	    fp.500.n <- 1 - exp( 3.11101712909049 ) / n^ 1.91401056721863
	    fp.875.n <- 1 - exp( 0.79473550581058 ) / n^ 1.10081930350091
	}
	if(p == 1) {
	    fp.500.n <- 1 - exp( 1.11098143415027 ) / n^ 1.5182890270453
	    fp.875.n <- 1 - exp( -0.66046776772861) / n^ 0.88939595831888
	}
    }

    stopifnot(0.5 <= alpha, alpha <= 1)
    if(alpha <= 0.875)
	fp.alpha.n <- fp.500.n + (fp.875.n - fp.500.n)/0.375 * (alpha - 0.5)
    else ##  0.875 < alpha <= 1
	fp.alpha.n <- fp.875.n + (1 - fp.875.n)/0.125 * (alpha - 0.875)
    return(1/fp.alpha.n)
} ## end{ MCDcnp2.rew }


.fastmcd <- function(x, quan, nsamp)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    ##	 parameters for partitioning
    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini

    ##	 vt::03.02.2006 - added options "best" and "exact" for nsamp
    if(!missing(nsamp)) {
	if(is.numeric(nsamp) && (nsamp < 0 || nsamp == 0 && p > 1)) {
	    warning("Invalid number of trials nsamp= ",nsamp," !Using default.\n")
	    nsamp <- -1
	} else if(nsamp == "exact" || nsamp == "best") {
	    myk <- p
	    if(n > 2*nmini-1) {
		warning("Options 'best' and 'exact' not allowed for n greater than ",
			2*nmini-1,". \nUsing nsamp= ",nsamp,"\n")
		nsamp <- -1
	    } else {
		nall <- choose(n, myk)
		if(nall > 5000 && nsamp == "best") {
		    nsamp <- 5000
		    warning("'nsamp = \"best\"' allows maximally 5000 subsets;\n",
			    "computing these subsets of size ",
                            myk," out of ",n,"\n")
		} else {
		    nsamp <- 0 ## all subsamples
		    if(nall > 5000)
			warning("Computing all ", nall, " subsets of size ",myk,
				" out of ",n,
				"\n This may take a very long time!\n",
				immediate. = TRUE)
		}
	    }
	}

	if(!is.numeric(nsamp) || nsamp == -1) { # still not defined - set it to the default
	    defCtrl <- rrcov.control() # default control
	    if(!is.numeric(nsamp))
		warning("Invalid number of trials nsamp= ",nsamp,
			" ! Using default nsamp= ",defCtrl$nsamp,"\n")
	    nsamp <- defCtrl$nsamp	# take the default nsamp
	}
    }

    storage.mode(x) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(quan) <- "integer"
    storage.mode(nsamp) <- "integer"

    ##	 Allocate temporary storage for the Fortran implementation,
    ##	 directly in the .Fortran() call.
    ##	  (if we used C, we'd rather allocate there, and be quite faster!)

    .Fortran("rffastmcd",
	     x,
	     n,
	     p,
	     nhalff = quan,
	     nsamp,
	     initcovariance = double(p * p),
	     initmean	    = double(p),
	     best	    = rep.int(as.integer(10000), quan),
	     mcdestimate = double(1),
	     weights   = integer(n),
	     exactfit  = integer(1), # output indicator: 0: ok; 1: ..., 2: ..
	     coeff     = matrix(double(5 * p), nrow = 5, ncol = p), ## plane
	     kount     = integer(1),
	     adjustcov = double(p * p),
	     integer(1),## << 'seed' no longer used -- FIXME
	     temp   = integer(n),
	     index1 = integer(n),
	     index2 = integer(n),
	     nmahad = double(n),
	     ndist  = double(n),
	     am	    = double(n),
	     am2    = double(n),
	     slutn  = double(n),

	     med   = double(p),
	     mad   = double(p),
	     sd	   = double(p),
	     means = double(p),
	     bmeans= double(p),
	     w	   = double(p),
	     fv1   = double(p),
	     fv2   = double(p),

	     rec   = double(p+1),
	     sscp1 = double((p+1)*(p+1)),
	     cova1 = double(p * p),
	     corr1 = double(p * p),
	     cinv1 = double(p * p),
	     cova2 = double(p * p),
	     cinv2 = double(p * p),
	     z	   = double(p * p),

	     cstock = double(10 * p * p),# (10,nvmax2)
	     mstock = double(10 * p),	 # (10,nvmax)
	     c1stock = double(km10 * p * p), # (km10,nvmax2)
	     m1stock = double(km10 * p), # (km10,nvmax)
	     dath = double(nmaxi * p),	 # (nmaxi,nvmax)

	     cutoff = qchisq(0.975, p),
	     chimed = qchisq(0.5,   p),

	     PACKAGE = "robustbase")[ ## keep the following ones:
	     c("initcovariance", "initmean", "best", "mcdestimate",
	       "weights", "exactfit", "coeff", "kount", "adjustcov") ]
}
