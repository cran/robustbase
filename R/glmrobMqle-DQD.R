#### Quasi-Deviance Differences	 --- used for Model Selection
#### -------------------------------------------------------- ./modsel.glmrob.R

## MM: These function names are really too long
##     but then, they are hidden in the name space ...

## (Maybe it would be nice to do this as one function with "family" .. )

glmrobMqleDiffQuasiDevB <- function(mu, mu0, y, ni, w.x, tcc)
{
    ##
    f.cnui <- function(u, y, ni, tcc)
    {
	## First part: nui
	pr <- u/ni
	Vmu <- pr * (1 - pr) ## = binomial()$variance
	residP <- (y-pr)*sqrt(ni/Vmu)
	nui <- pmax(-tcc, pmin(residP,tcc))

	## Second part: Enui
	H <- floor(u - tcc*sqrt(ni*Vmu))
	K <- floor(u + tcc*sqrt(ni*Vmu))
	## Actually, floor is not needed because pbinom() can cope
	## with noninteger values in the argument q!
	## what follows is similar to glmrob.Mqle.EpsiB except a
	## different vectorisation
	h1 <- (if(ni == 1) as.numeric(- (H < 0) + (K >= 1) ) * sqrt(Vmu)
	       else
	       (pbinom(K-1,1,pr) - pbinom(H-1,ni-1,pr)
		- pbinom(K,ni,pr) + pbinom(H,ni,pr)) * pr * sqrt(ni/Vmu))
	## pmax was needed to get numeric returns from pbinom

	Enui <- (tcc*(1 - pbinom(K,ni,pr) - pbinom(H,ni,pr)) + h1)

	return((nui - Enui) / sqrt(ni*Vmu))
    } ## f.cnui()

    nobs <- length(mu)
    stopifnot(nobs > 0)
    QMi <- numeric(nobs)
    ## Numerical integrations
    for(i in 1:nobs)
	QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], tcc = tcc,
			    subdivisions = 200,
			    lower = mu[i]*ni[i], upper = mu0[i]*ni[i])$value
    ## robust quasi-deviance
    ## -2*(sum(QMi1)-sum(QMi2))	   ## Andreas' interpretation of (4) and (5)
    ## -2*(sum(QMi1)-sum(QMi2)/nobs)  ## Eva's interpretation of (4) and (5)
    ## According to Andreas' interpretation
    -2*sum(QMi*w.x)
}

glmrobMqleDiffQuasiDevPois <- function(mu, mu0, y, ni, w.x, tcc)
{
    ##
    f.cnui <- function(u, y, ni, tcc)
    {
	## First part: nui
	Vmu <- u ## = poisson()$variance
	residP <- (y-u)/sqrt(Vmu)
	nui <- pmax(-tcc,pmin(residP,tcc))

	## Second part: Enui
	H <- floor(u - tcc*sqrt(Vmu))
	K <- floor(u + tcc*sqrt(Vmu))
	## what follows is similar to glmrob.Mqle.Epsipois except a
	## different vectorisation
	h1 <- u/sqrt(Vmu)*(dpois(H,u)- dpois(K,u))
	Enui <- tcc*(1 - ppois(K,u) - ppois(H,u)) + h1

	return((nui - Enui) / sqrt(Vmu))
    }

    nobs <- length(mu)
    stopifnot(nobs > 0)
    QMi <- numeric(nobs)
    ## Numerical integrations
    for(i in 1:nobs)
	QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], tcc = tcc,
			    lower = mu[i], upper = mu0[i])$value

    ## robust quasi-deviance
    ## -2*(sum(QMi1)-sum(QMi2))	  ## Andreas' interpretation of (4) and (5)
    ## -2*(sum(QMi1)-sum(QMi2)/nobs) ## Eva's interpretation of (4) and (5)
    ## According to Andreas' interpretation
    -2*sum(QMi*w.x)
}

