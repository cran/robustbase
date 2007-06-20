#### These Chi() and Psi() functions are currently used by lmrob() functions
#### FIXME: integrate these with the psi-rho objects --> ./psi-rho-funs.R

## Chi'() is just a scaled version of psi() :
## with current scale (new for psi()):
##	 i)  Chi'(x, c) == (6/c^2) Psi(x,c)
## ==>	 ii) Chi''(x,c) == (6/c^2) Psi'(x,c)
## and       Chi (x, c) == (6/c^2) Rho(x,c), where Psi(.) = Rho'(.)


tukeyChi <- function(x, cc, deriv = 0)
{
    x <- x / cc
    x2 <- x*x
    out <- x2 > 1
    switch(deriv + 1,
       {  ## deriv = 0
	   r <- x2*(3 + x2*(-3 + x2))
	   r[out] <- 1
       },
       {  ## deriv = 1
	   r <- 6/cc * x * (1-x2)^2
	   r[out] <- 0
       },
       {  ## deriv = 2
	   r <- 6/(cc^2) * (1 - x2) * (1 - 5*x2)
	   r[out] <- 0
       },
       stop("deriv must be in {0,1,2}"))
    r
}

## we call this  '*Psi1'  such as to not be confounded with
## the (future!) S4 object tukeyPsi1() !
tukeyPsi1 <- function(x, cc, deriv = 0)
{
    ## This version of psi() is scaled such that psi'(0) = 1
    x2 <- (x / cc)^2
    out <- x2 > 1
    switch(deriv + 2,
       {  ## deriv = -1
           c. <- cc^2/6
	   r <- c.*(1 - (1- x2)^3)
	   r[out] <- c.
       },
       {  ## deriv = 0
	   r <- x * (1-x2)^2
	   r[out] <- 0
       },
       {  ## deriv = 1
	   r <- (1 - x2) * (1 - 5*x2)
	   r[out] <- 0
       },
       stop("deriv must be in {-1,0,1}"))
    r
}

if(FALSE)
tukeyPsi1Ex <- function (x, cc, deriv = 0)
## tukeyPsi1Ex <- function (x, cc = 4.685, deriv = 0)
##                               ^^^^^^^^^
{
  ## This version of psi() is scaled such that psi'(0) = 1
  u <- pmin((x/cc)^2, 1)
  if(deriv < 0)
    return((1 - (1-u)^3)*cc^2/6)
  if(deriv == 0)
    return(x * (1 - u)^2)
  return((1 - u) * (1 - 5 * u))
}
