#### These are extension of the MASS functions of the same name
#### (and the same "flaw": returning  psi(x)/x for deriv = 0 !)
#### with the additional option	 rho = FALSE / TRUE
#### needed for	 glmrob() model selection (deviance - difference computations)

### by Andreas Ruckstuhl

### FIXME: replace by using the objects in  ./psi-rho-funs.R
###        and                              ./biweight-funs.R

psi.bisquare <- function (u, c=4.685, deriv=0, rho=FALSE)
{
    t <- (u/c)^2
    if (deriv)
	return(ifelse(t < 1, (1 - t) * (1 - 5 * t), 0))
    if(rho)
	return((1 - (1-pmin(1, t))^3)/6)
    (1 - pmin(1, abs(u/c))^2)^2
}


psi.hampel <- function (u, a=2, b=4, c=8, deriv=0, rho=FALSE)
{
    U <- pmin(abs(u) + 1e-50, c)
    if (deriv)
	return(ifelse(abs(u) <= c, ifelse(U <= a, 1,
			 ifelse(U <= b, 0, -a/(c - b))), 0))
    if(rho)
	return(NULL)
    ifelse(U <= a, U, ifelse(U <= b, a, a * (c - U)/(c - b)))/U
}


psi.huber <- function (u, k=1.345, deriv=0, rho=FALSE)
{
    if (deriv)
	return(abs(u) <= k)
    if(rho)
	return(ifelse(abs(u)<=k, u^2/2, k*(abs(u) - k/2)))
    pmin(1, k/abs(u))
}

