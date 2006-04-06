c#### These Chi() and Psi() functions are currently used by lmrob() functions
#### FIXME: integrate these with the psi-rho objects --> ./psi-rho-funs.R

lmrob.Chi <- function(x, cc, deriv = 0)
{

    x <- x / cc
    out <- abs(x) > 1
    switch(deriv + 1,
       {  ## deriv = 0
           x <- x*x
           r <- x*(3 + x*(-3 + x))
           r[out] <- 1
       },
       {  ## deriv = 1
           r <- 6/cc *x*(1-x^2)^2
           r[out] <- 0
       },
       {  ## deriv = 2
           x <- x*x # x^2
           r <- 6/(cc^2) * (1 - x) * (1 - 5*x)
           r[out] <- 0
       },
       stop("deriv must be in {0,1,2}"))
    r
}

## MM (FIXME?): why use  Psi() additional to Chi() ?
## --  Chi'() is just a scaled version of psi(); i.e. chi() = rho() ??

lmrob.Psi <- function(x, cc, deriv = 0)
{
    x <- x / cc
    r <- if(deriv == 0) x*(1 - x^2)^2 else (1 - x^2)*(1 - 5* x^2) / cc
    r[ abs(x) > 1 ] <- 0
    r
}
