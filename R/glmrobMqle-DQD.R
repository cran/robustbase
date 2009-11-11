#### Quasi-Deviance Differences	 --- for Model Selection
#### --------------------------------------------------- -> ./anova-glmrob.R

## MM: These function names are really too long
##     but then, they are hidden in the name space ...

## (Maybe it would be nice to do this as one function with "family" .. )

glmrobMqleDiffQuasiDevB <- function(mu, mu0, y, ni, w.x, phi, tcc)
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

glmrobMqleDiffQuasiDevPois <- function(mu, mu0, y, ni, w.x, phi, tcc)
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
	## what follows is similar to Epsipois except a
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

glmrobMqleDiffQuasiDevGamma <- function(mu, mu0, y, ni, w.x, phi, tcc)
{
    ## Notation similar to the discrete case (Cantoni & Ronchetti, 2001)
    f.cnui <- function(u, y, ni, phi, tcc)
    {
	## First part: nui
	sV <- sqrt(phi) * u ## = sqrt(dispersion * Gamma()$variance)
	residP <- (y-u)/sV
	nui <- pmax(-tcc,pmin(residP,tcc))

	## Second part: Enui
        ## what follows is similar to glmrob.Mqle.Epsipois except a
	## different vectorisation
        nu <- 1/phi      ## form parameter nu
        snu <- 1/sqrt(phi) ## sqrt (nu)

        pPtmc <- pgamma(snu - tcc, shape=nu, rate=snu)
        pPtpc <- pgamma(snu + tcc, shape=nu, rate=snu)
        #GLtcc <- Gmn(-tcc,nu)
        #GUtcc <- Gmn( tcc,nu)

        Enui <- tcc*(1-pPtpc-pPtmc) + Gmn(-tcc,nu) - Gmn( tcc,nu)

	return((nui/ sV - Enui/u*sqrt(phi) ))
    }
    f.cnui1 <- function(u, y, ni, phi, tcc)
    {
	## First part: nui
	sV <- sqrt(phi) * u ## = sqrt(dispersion * Gamma()$variance)
	residP <- (y-u)/sV
	nui <- pmax(-tcc,pmin(residP,tcc))

       	return(nui  / sV)
    }
    f.cnui2 <- function(u, y, ni, phi, tcc)
    {
	## First part: nui
	sV <- sqrt(phi) * u ## = sqrt(dispersion * Gamma()$variance)
        snu <- 1/sqrt(phi) ## sqrt (nu)

	## Second part: Enui
        ## what follows is similar to EpsiGamma except a
	## different vectorisation
        nu <- 1/phi      ## form parameter nu

        pPtmc <- pgamma(snu - tcc, shape=nu, rate=snu)
        pPtpc <- pgamma(snu + tcc, shape=nu, rate=snu)

        Enui <- tcc*(1-pPtpc-pPtmc) + Gmn(-tcc,nu) - Gmn( tcc,nu)
	return(Enui/u * sqrt(phi))
    }
    f.cnui22 <- function(u, y=y[1], ni=1, phi=phi, tcc=1.5)
    {
	## First part: nui
	sV <- sqrt(phi) * u ## = sqrt(dispersion * Gamma()$variance)

	## Second part: Enui
        ## what follows is similar to glmrob.Mqle.Epsipois except a
	## different vectorisation
        nu <- 1/phi      ## form parameter nu
        snu <- 1/sqrt(phi)

        pPtmc <- pgamma(snu - tcc, shape=nu, rate=snu)
        pPtpc <- pgamma(snu + tcc, shape=nu, rate=snu)

        Enui <- tcc*(1-pPtpc-pPtmc) + Gmn(-tcc,nu) - Gmn( tcc,nu)
	return((Enui) )
    }



    nobs <- length(mu)
    stopifnot(nobs > 0)
    QMi <- numeric(nobs)
    ## Numerical integrations
    for(i in 1:nobs)
	QMi[i] <- integrate(f.cnui, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
			    lower = mu[i], upper = mu0[i])$value

    ## robust quasi-deviance
    ## -2*(sum(QMi1)-sum(QMi2))	  ## Andreas' interpretation of (4) and (5)
    ## -2*(sum(QMi1)-sum(QMi2)/nobs) ## Eva's interpretation of (4) and (5)
    ## According to Andreas' interpretation
    QM1i <- QM2i <- numeric(nobs)
for(i in 1:nobs)
    QM1i[i] <- integrate(f.cnui1, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
			    lower = mu[i], upper = mu0[i])$value

for(i in 1:nobs)
    QM2i[i] <- integrate(f.cnui2, y = y[i], ni = ni[i], phi=phi, tcc = tcc,
			    lower = mu[i], upper = mu0[i])$value

## browser()

    -2*sum(QMi*w.x)
}

