#### This is from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov


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

##  I would like to thank Peter Rousseeuw and Katrien van Driessen for
##  providing the initial code of this function.


## hidden in namespace:
quan.f <- function(alpha, n, rk) {
    ## Compute size of subsample -- Same function for covMcd() and ltsReg()
    n2 <- floor((n+rk+1)/2)
    floor(2 * n2 - n + 2 * (n - n2) * alpha)
}

covMcd <- function(x,
		   cor = FALSE,
		   alpha = 1/2,
		   nsamp = 500,
		   seed = 0,
		   print.it = FALSE,
		   control)
{
    ##	 Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    if(!missing(control)) {
	defcontrol <- rrcov.control()	# default control
	if(alpha == defcontrol$alpha)	alpha <- control$alpha
	if(nsamp == defcontrol$nsamp)	nsamp <- control$nsamp
	if(seed	 == defcontrol$seed)	 seed <- control$seed
	if(print.it == defcontrol$print.it) print.it <- control$print.it
    }

    ## vt:: tolerance to be used for computing the mahalanobis distances (default = 1e-7)
    tol = 1e-10

    if(is.vector(x) || (is.matrix(x) && !is.data.frame(x))) {
	if(!is.numeric(x))
	    stop("x is not a numeric dataframe or matrix.")
    }
    if((!is.vector(x) && !is.matrix(x)) || is.data.frame(x)) {
	if((!is.data.frame(x) && !is.numeric(x)) || (!all(sapply(x,data.class) == "numeric")))
	    stop("x is not a numeric dataframe or matrix.")
    }

    ##vt:: if the data is supplied as a data.frame, the following expressions results in an error
    ## as workaround convert the data.frame to a matrix
    if(is.data.frame(x))
	x <- as.matrix(x)
    else if(!is.matrix(x)) {
	x <- array(x, c(length(x), 1),
		   list(names(x), deparse(substitute(data))))
	x <- as.matrix(x)
    }
    dimn <- dimnames(x)
    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
	stop("All observations have missing values!")
    n <- dx[1]
    p <- dx[2]
    if(n < 2 * p)
	stop("Need at least 2*(number of variables) observations ")
    jmin <- floor((n + p + 1)/2)
    if(alpha < 1/2) {
	stop("The MCD must cover at least", jmin, "observations")
    }
    else if(alpha > 1)
	stop("alpha is out of range")
    quan <- quan.f(alpha, n, p)

    ## alpha = 1 : Compute the classical estimates
    if(alpha == 1) {
	mcd <- cov.wt(x)$cov
	loc <- as.vector(colMeans(x))
	obj <- determinant(mcd, log = TRUE)$modulus[1]
	if(( - obj/p) > 50) {
	    ans <- list()
	    ans$cov <- mcd
	    dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
	    ans$center <- loc
	    if(length(dimn[[2]]))
		names(ans$center) <- dimn[[2]]
	    ans$n.obs <- n
	    ans$call <- match.call()
	    ans$method <- paste("Minimum Covariance Determinant Estimator.")
	    ans$method <- paste(ans$method,
				"\nThe classical covariance matrix is singular.")
	    if(!print.it) {
		cat("The classical covariance matrix is singular.\n")
	    }
	    ans$alpha <- alpha
	    ans$quan <- quan
	    ans$raw.cov <- mcd
	    dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
	    ans$raw.center <- loc
	    if(length(dimn[[2]]))
		names(ans$raw.center) <- dimn[[2]]
	    ans$crit <- exp(obj)
	    ans$mcd.wt <- rep(NA, length(na.x))
	    ans$mcd.wt[ok] <- rep(1, sum(ok == TRUE))
	}
	else {
	    mah <- mahalanobis(x, loc, mcd, tol = tol) # VT:: 01.09.2004 - bug in alpha=1
	    ## (tol instead of tol.inv as parameter name)
	    weights <- ifelse(mah < qchisq(0.975, p), 1, 0)
	    ans <- cov.wt(x, wt = weights, cor)
	    ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov

	    ## Consistency factor for reweighted MCD
	    if(sum(weights) == n)
		cdelta.rew <- 1
	    else {
		qdelta.rew <- qchisq(sum(weights)/n, p)
		cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 +
					   1)/(sum(weights)/n)
		cdelta.rew <- 1/cdeltainvers.rew
	    }
	    ans$cov <- ans$cov * cdelta.rew
	    ans$call <- match.call()
	    ans$method <- paste("Minimum Covariance Determinant Estimator.")
	    if( - (determinant(ans$cov, log = TRUE)$modulus[1] - 0)/p > 50) {
		ans$method <- paste(ans$method, "\nThe reweighted MCD scatter matrix is singular.")
		if(!print.it) {
		    cat("The reweighted MCD scatter matrix is singular.\n")
		}
		ans$alpha <- alpha
		ans$quan <- quan
		ans$raw.cov <- mcd
		dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
		ans$raw.center <- loc
		if(length(dimn[[2]]))
		    names(ans$raw.center) <- dimn[[2]]
		ans$crit <- exp(obj)
		ans$mcd.wt <- rep(NA, length(na.x))
		ans$mcd.wt[ok] <- weights
		if(length(dimn[[1]]))
		    names(ans$mcd.wt) <- dimn[[1]]
		ans$wt <- NULL
		ans$X <- x
		if(length(dimn[[1]]))
		    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
		else {
		    xx <- seq(1, length(na.x))
		    dimnames(ans$X) <- list(NULL, NULL)
		    dimnames(ans$X)[[1]] <- xx[ok]
		}
		ans$method <- paste(ans$method, "\nThe minimum covariance determinant estimates based on", n,
				    "observations \nare equal to the classical estimates.")
		if(print.it) {
		    cat(ans$method, "\n")
		}
		class(ans) <- "mcd"
		## have '$call' already: attr(ans, "call") <- sys.call()
		return(ans)
	    }
	    else {
		ans$alpha <- alpha
		ans$quan <- quan
		ans$raw.cov <- mcd
		dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
		ans$raw.center <- loc
		if(length(dimn[[2]]))
		    names(ans$raw.center) <- dimn[[2]]
		ans$crit <- exp(obj)
		mah <- mahalanobis(x, ans$center, ans$cov, tol = tol)
	    }
	    ans$mcd.wt <- rep(NA, length(na.x))
	    ans$mcd.wt[ok] <- ifelse(mah < qchisq(0.975, p), 1, 0)
	}
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	ans$wt <- NULL
	ans$X <- x
	if(length(dimn[[1]]))
	    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
	else {
	    xx <- seq(1, length(na.x))
	    dimnames(ans$X) <- list(NULL, NULL)
	    dimnames(ans$X)[[1]] <- xx[ok]
	}
	ans$method <- paste(ans$method,
			    "\nThe minimum covariance determinant estimates based on",
			    n, "observations \nare equal to the classical estimates."
			    )
	if(print.it) {
	    cat(ans$method, "\n")
	}
	class(ans) <- "mcd"
	## have '$call' already: attr(ans, "call") <- sys.call()
	return(ans)
    } ## end {alpha=1} --

    mcd <- .fastmcd(x, quan, nsamp, seed)

    ## Compute the consistency correction factor for the raw MCD
    ##	(see calfa in Croux and Haesbroeck)
    qalpha <- qchisq(quan/n, p)
    calphainvers <- pgamma(qalpha/2, p/2 + 1)/(quan/n)
    calpha <- 1/calphainvers
    correct <- MCDcorrfactor.s(p, n, alpha)

    if(p == 1) {
	## The number of variables is 1 - compute univariate location and scale estimates
	scale <- sqrt(calpha) * as.double(mcd$initcovariance) * sqrt(correct)

	center <- as.double(mcd$initmean)
	if(abs(scale - 0) < 1e-07) {
	    ans <- list()
	    ## VT:: 22.12.04 - ans$cov must be a matrix. For example a subsequent
	    ##	 call to determinant() will raise an error if it is a double.
	    ##	 The same for ans$raw.cov - see below

	    ## ans$cov <- 0
	    ans$cov <- matrix(0)
	    names(ans$cov) <- dimn[[2]][1]
	    ans$center <- center
	    names(ans$center) <- dimn[[2]][1]
	    ans$n.obs <- n
	    ans$call <- match.call()
	    ans$method <- paste("Univariate location and scale estimation.\nMore than",
				quan, "of the observations are identical.")
	    ans$alpha <- alpha
	    ans$quan <- quan
	    ## ans$raw.cov <- 0
	    ans$raw.cov <- matrix(0)
	    names(ans$raw.cov) <- dimn[[2]][1]
	    ans$raw.center <- center
	    names(ans$raw.center) <- dimn[[2]][1]
	    ans$crit <- 0
	    ans$mcd.wt <- rep(NA, length(na.x))
	    ans$mcd.wt[ok] <- as.vector(ifelse(abs(x - center) < 1e-07, 1, 0))
	    if(length(dimn[[1]]))
		names(ans$mcd.wt) <- dimn[[1]]
	    if(print.it) {
		cat(ans$method, "\n")
	    }
	    ans$X <- x
	    if(length(dimn[[1]]))
		dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
	    else {
		xx <- seq(1, length(na.x))
		dimnames(ans$X) <- list(NULL, NULL)
		dimnames(ans$X)[[1]] <- xx[ok]
	    }

	    class(ans) <- "mcd"
	    ## have '$call' already: attr(ans, "call") <- sys.call()
	    return(ans)
	}

	## Compute the weights for the raw MCD in case p=1
	quantiel <- qchisq(0.975,p)
	weights <- ifelse(((x - center)/scale)^2  < quantiel, 1, 0)
	ans <- cov.wt(x, wt = weights, cor = cor)
	ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov

	## Apply the correction factor for the reweighted cov
	if(sum(weights) == n) {
	    cdelta.rew <- 1
	    correct.rew <- 1
	}
	else {
	    qdelta.rew <- qchisq(sum(weights)/n, p)
	    cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(weights)/n)
	    cdelta.rew <- 1/cdeltainvers.rew
	    correct.rew <- MCDcorrfactor.rew.s(p, n, alpha)
	}
	ans$cov <- ans$cov * cdelta.rew * correct.rew
	ans$call <- match.call()
	ans$method <- paste("Univariate location and scale estimation.")
	ans$alpha <- alpha
	ans$quan <- quan
	ans$raw.cov <- scale^2
	names(ans$raw.cov) <- dimn[[2]][1]
	ans$raw.center <- as.vector(center)
	names(ans$raw.center) <- dimn[[2]][1]
	ans$crit <- 1/(quan - 1) *
	    sum(sort((x - as.double(mcd$initmean))^2, quan)[1:quan])
	center <- ans$center
	scale <- as.vector(sqrt(ans$cov))
	ans$mcd.wt <- rep(NA, length(na.x))
	weights <- ifelse(((x - center)/scale)^2 < qchisq(0.975, p), 1, 0)

	ans$mcd.wt[ok] <- weights
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	ans$wt <- NULL
	if(print.it) {
	    cat(ans$method, "\n")
	}
	ans$X <- x
	if(length(dimn[[1]]))
	    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
	else {
	    xx <- seq(1, length(na.x))
	    dimnames(ans$X) <- list(NULL, NULL)
	    dimnames(ans$X)[[1]] <- xx[ok]
	}
	class(ans) <- "mcd"
	## have '$call' already: attr(ans, "call") <- sys.call()
	return(ans)
    } ## end p=1

    ## else  p >= 2 : ----------------------------------------------------------

    ## Apply correction factor to the raw estimates and use them to compute weights
    mcd$initcovariance <- calpha * mcd$initcovariance * correct
    dim(mcd$initcovariance) <- c(p, p)

    msg <- paste("Minimum Covariance Determinant Estimator")

## If not all observations are in general position, i.e. more than h observations lie on
## a hyperplane, the program still yields the MCD location and scatter matrix,
## the latter being singular (as it should be), as well as the equation of the hyperplane.
    if(mcd$exactfit != 0) {
	dim(mcd$coeff) <- c(5, p)
	ans <- list()
	ans$cov <- mcd$initcovariance
	dimnames(ans$cov) <- list(dimn[[2]], dimn[[2]])
	ans$center <- as.vector(mcd$initmean)
	if(length(dimn[[2]]))
	    names(ans$center) <- dimn[[2]]
	ans$n.obs <- n
	ans$call <- match.call()
	ans$method <- msg
	if(mcd$exactfit == -1) {
	    stop("The program allows for at most ", mcd$kount,
		 " observations.")
	}
	if(mcd$exactfit == -2) {
	    stop("The program allows for at most ", mcd$kount,
		 " variables.")
	}
	if(mcd$exactfit == 1) {
	    ans$method <- paste(ans$method,
				"\nThe covariance matrix of the data is singular.")
	    if(!print.it) {
		cat("The covariance matrix of the data is singular.\n")
	    }
	}
	if(mcd$exactfit == 2) {
	    ans$method <- paste(ans$method,
				"\nThe covariance matrix has become singular during\nthe iterations of the MCD algorithm."
				)
	    if(!print.it) {
		cat("The covariance matrix has become singular during\nthe iterations of the MCD algorithm.\n"
		    )
	    }
	}
	if(p == 2) {
	    ans$method <- paste(ans$method, "\nThere are", mcd$
				kount,
				"observations in the entire dataset of\n", n,
				"observations that lie on the line with equation\n",
				round(mcd$coeff[1, 1], digits = 4),
				"(x_i1-m_1)+", round(mcd$coeff[1, 2], digits =
						     4),
				"(x_i2-m_2)=0 \nwith (m_1,m_2) the mean of these observations."
				)
	    if(!print.it) {
		cat("There are", mcd$kount,
		    "observations in the entire dataset of\n", n,
		    "observations that lie on the line with equation\n",
		    round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+",
		    round(mcd$coeff[1, 2], digits = 4),
		    "(x_i2-m_2)=0 \nwith (m_1,m_2) the mean of these observations.\n"
		    )
	    }
	}
	else if(p == 3) {
	    ans$method <- paste(ans$method, "\nThere are", mcd$
				kount,
				"observations in the entire dataset of\n", n,
				"observations that lie on the plane with equation \n",
				round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+",
				round(mcd$coeff[1, 2], digits = 4), "(x_i2-m_2)+",
				round(mcd$coeff[1, 3], digits = 4),
				"(x_i3-m_3)=0 \nwith (m_1,m_2) the mean of these observations."
				)
	    if(!print.it) {
		cat("There are", mcd$kount,
		    "observations in the entire dataset of\n",
		    n, "observations that lie on the plane with equation\n",
		    round(mcd$coeff[1, 1], digits = 4), "(x_i1-m_1)+",
		    round(mcd$coeff[1, 2], digits = 4), "(x_i2-m_2)+",
		    round(mcd$coeff[1, 3], digits = 4), "(x_i3-m_3)=0 \n",
		    "with (m_1,m_2) the mean of these observations.\n"
		    )
	    }
	}
	else { ##  p > 3 -----------
	    ans$method <-
		paste(ans$method, "\nThere are",
		      mcd$kount, " observations in the entire dataset of\n",
		      n, "observations that lie on the hyperplane with equation\n",
		      "a_1*(x_i1-m_1)+...+a_p*(x_ip-m_p)=0 \n",
		      "with (m_1,...,m_p) the mean\n",
		      "of these observations and coefficients a_i equal to:\n")
	    if(!print.it) {
		cat("There are", mcd$kount,
		    " observations in the entire dataset of\n", n,
		    "observations that lie on the hyperplane with equation \na_1*(x_i1-m_1)+...+a_p*(x_ip-m_p)=0 \nwith (m_1,...,m_p) the mean\nof these observations and coefficients a_i equal to: \n"
		    )
	    }

	    for(i in 1:p) {
		ans$method <-
		    paste(ans$method, round(mcd$coeff[ 1, i], digits = 4))
	    }
	    if(!print.it)
		print(round(mcd$coeff[1,  ], digits = 4))
	} # end {p > 3}
	ans$alpha <- alpha
	ans$quan <- quan
	ans$raw.cov <- mcd$initcovariance
	dimnames(ans$raw.cov) <- list(dimn[[2]], dimn[[2]])
	ans$raw.center <- as.vector(mcd$initmean)
	if(length(dimn[[2]]))
	    names(ans$raw.center) <- dimn[[2]]
	ans$crit <- 0
	ans$mcd.wt <- rep(NA, length(na.x))
	ans$mcd.wt[ok] <- mcd$weights
	if(length(dimn[[1]]))
	    names(ans$mcd.wt) <- dimn[[1]]
	ans$wt <- NULL
	ans$X <- x
	if(length(dimn[[1]]))
	    dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
	else {
	    xx <- seq(1, length(na.x))
	    dimnames(ans$X) <- list(xx[ok], NULL)
	}

	if(print.it) {
	    cat(ans$method, "\n")
	}
	class(ans) <- "mcd"
	## have '$call' already: attr(ans, "call") <- sys.call()
	return(ans)
    } #end exact fit <==>  (mcd$exactfit != 0)

    ## else ------ exactfit == 0 ----------------------------------------------

    mah <- mahalanobis(x, mcd$initmean, mcd$initcovariance, tol = tol)
    mcd$weights <- ifelse(mah < qchisq(0.975, p), 1, 0)

    weights <- mcd$weights
    weights <- as.vector(weights)

## Compute and apply the consistency correction factor for the reweighted cov
    if(sum(weights) == n) {
	cdelta.rew <- 1
	correct.rew <- 1
    }
    else {
	qdelta.rew <- qchisq(sum(weights)/n, p)
	cdeltainvers.rew <- pgamma(qdelta.rew/2, p/2 + 1)/(sum(weights)/n)
	cdelta.rew <- 1/cdeltainvers.rew
	correct.rew <- MCDcorrfactor.rew.s(p, n, alpha)
    }

    ans <- cov.wt(x, wt = weights, cor)
    ans$cov <- sum(weights)/(sum(weights) - 1) * ans$cov
    ans$cov <- ans$cov * cdelta.rew * correct.rew
    ans$call <- match.call()
    ans$method <- msg

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
	ans$method <- paste(ans$method, "\nThe reweighted MCD scatter matrix is singular.")
	if(!print.it) {
	    cat("The reweighted MCD scatter matrix is singular.\n")
	}
	ans$mah <- ans$raw.mah
    }
    else {
	mah <- mahalanobis(x, ans$center, ans$cov, tol = tol)
	ans$mah <- mah
	weights <- ifelse(mah < qchisq(0.975, p), 1, 0)
    }
    ans$mcd.wt <- rep(NA, length(na.x))
    ans$mcd.wt[ok] <- weights
    if(length(dimn[[1]]))
	names(ans$mcd.wt) <- dimn[[1]]
    ans$wt <- NULL
    ans$X <- x
    if(length(dimn[[1]]))
	dimnames(ans$X)[[1]] <- names(ans$mcd.wt)[ok]
    else {
	xx <- seq(1, length(na.x))
	dimnames(ans$X) <- list(xx[ok], NULL)
    }
    if(print.it)
	cat(ans$method, "\n")
    class(ans) <- "mcd"
    ## have '$call' already: attr(ans, "call") <- sys.call()
    return(ans)
}

### --- Namespace hidden (but parsed once and for all) : -------------

MCDcorrfactor.s <- function(p, n, alpha)
{
    if(p > 2) {
        ##                              "alfaq"            "betaq"    "qwaarden"
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
} ## end{ MCDcorrfactor.s }

MCDcorrfactor.rew.s <- function(p, n, alpha)
{
    if(p > 2) {
        ##                              "alfaq"            "betaq"    "qwaarden"
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
} ## end{ MCDcorrfactor.rew.s }

.fastmcd <- function(x, quan, nsamp, seed)
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    deter <-  fit <- kount <- 0
    cutoff <- qchisq(0.975,p)
    chimed <- qchisq(0.5, p)

    storage.mode(x) <- "double"
    storage.mode(n) <- "integer"
    storage.mode(p) <- "integer"
    storage.mode(quan) <- "integer"
    storage.mode(nsamp) <- "integer"
    storage.mode(seed) <- "integer"

    storage.mode(deter) <- "double"
    storage.mode(fit) <- "integer"
    storage.mode(kount) <- "integer"
    storage.mode(cutoff) <- "double"
    storage.mode(chimed) <- "double"


    initcov <- matrix(0, nrow = p * p, ncol = 1)
    adcov <- matrix(0, nrow = p * p, ncol = 1)
    initmean <- matrix(0, nrow = p, ncol = 1)
    inbest <- matrix(10000, nrow = quan, ncol = 1)
    plane <- matrix(0, nrow = 5, ncol = p)
    weights <- matrix(0, nrow = n, ncol = 1)

    storage.mode(initcov) <- "double"
    storage.mode(adcov) <- "double"
    storage.mode(initmean) <- "double"
    storage.mode(inbest) <- "integer"
    storage.mode(plane) <- "double"
    storage.mode(weights) <- "integer"

    ##	 Allocate temporary storage for the fortran implementation
    temp <- matrix(0, nrow = n, ncol = 1)
    index1 <- matrix(0, nrow = n, ncol = 1)
    index2 <- matrix(0, nrow = n, ncol = 1)
    nmahad <- matrix(0, nrow = n, ncol = 1)
    ndist <- matrix(0, nrow = n, ncol = 1)
    am <- matrix(0, nrow = n, ncol = 1)
    am2 <- matrix(0, nrow = n, ncol = 1)
    slutn <- matrix(0, nrow = n, ncol = 1)

    med <- matrix(0, nrow = p, ncol = 1)
    mad <- matrix(0, nrow = p, ncol = 1)
    sd <- matrix(0, nrow = p, ncol = 1)
    means <- matrix(0, nrow = p, ncol = 1)
    bmeans <- matrix(0, nrow = p, ncol = 1)
    w <- matrix(0, nrow = p, ncol = 1)
    fv1 <- matrix(0, nrow = p, ncol = 1)
    fv2 <- matrix(0, nrow = p, ncol = 1)

    rec <- matrix(0, nrow = p+1, ncol = 1)
    sscp1 <- matrix(0, nrow = (p+1)*(p+1), ncol = 1)
    cova1 <- matrix(0, nrow = p*p, ncol = 1)
    corr1 <- matrix(0, nrow = p*p, ncol = 1)
    cinv1 <- matrix(0, nrow = p*p, ncol = 1)
    cova2 <- matrix(0, nrow = p*p, ncol = 1)
    cinv2 <- matrix(0, nrow = p*p, ncol = 1)
    z <- matrix(0, nrow = p*p, ncol = 1)

    kmini <- 5
    nmini <- 300
    km10 <- 10*kmini
    nmaxi <- nmini*kmini
    cstock <- matrix(0, nrow = 10*p*p, ncol = 1)    #(10,nvmax2)
    mstock <- matrix(0, nrow = 10*p, ncol = 1)	    #(10,nvmax)
    c1stock <- matrix(0, nrow = km10*p*p, ncol = 1) #(km10,nvmax2)
    m1stock <- matrix(0, nrow = km10*p, ncol = 1)   #(km10,nvmax)
    dath <- matrix(0, nrow = nmaxi*p, ncol = 1)	    #(nmaxi,nvmax)

    storage.mode(temp) <- "integer"
    storage.mode(index1) <- "integer"
    storage.mode(index2) <- "integer"
    storage.mode(nmahad) <- "double"
    storage.mode(ndist) <- "double"
    storage.mode(am) <- "double"
    storage.mode(am2) <- "double"
    storage.mode(slutn) <- "double"

    storage.mode(med) <- "double"
    storage.mode(mad) <- "double"
    storage.mode(sd) <- "double"
    storage.mode(means) <- "double"
    storage.mode(bmeans) <- "double"
    storage.mode(w) <- "double"
    storage.mode(fv1) <- "double"
    storage.mode(fv2) <- "double"

    storage.mode(rec) <- "double"
    storage.mode(sscp1) <- "double"
    storage.mode(cova1) <- "double"
    storage.mode(corr1) <- "double"
    storage.mode(cinv1) <- "double"
    storage.mode(cova2) <- "double"
    storage.mode(cinv2) <- "double"
    storage.mode(z) <- "double"

    storage.mode(cstock) <- "double"
    storage.mode(mstock) <- "double"
    storage.mode(c1stock) <- "double"
    storage.mode(m1stock) <- "double"
    storage.mode(dath) <- "double"

    mcd <- .Fortran("rffastmcd",
		    x,
		    n,
		    p,
		    quan,
		    nsamp,
		    initcovariance = initcov,
		    initmean = initmean,
		    best = inbest,
		    mcdestimate = deter,
		    weights = weights,
		    exactfit = fit,
		    coeff = plane,
		    kount = kount,
		    adjustcov = adcov,
		    seed,
		    temp,
		    index1,
		    index2,
		    nmahad,
		    ndist,
		    am,
		    am2,
		    slutn,
		    med,
		    mad,
		    sd,
		    means,
		    bmeans,
		    w,
		    fv1,
		    fv2,
		    rec,
		    sscp1,
		    cova1,
		    corr1,
		    cinv1,
		    cova2,
		    cinv2,
		    z,
		    cstock,
		    mstock,
		    c1stock,
		    m1stock,
		    dath,
		    cutoff,
		    chimed,
		    PACKAGE = "robustbase")
    ## FIXME? -- do not return *everything*
    mcd
}
