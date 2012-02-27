require(robustbase)
## see also ./lmrob-psifns.R <<<<<<<<

## Demonstrate that  one of  tukeyChi() / tukeyPsi1() is superfluous

EQ <- function(x,y) all.equal(x,y, tol = 1e-13)

x <- seq(-4,4, length=201)
c. <- pi

stopifnot(EQ(tukeyChi(x, c.),
             6/c.^2* tukeyPsi1(x, c., deriv=-1)),
          EQ(tukeyChi(x, c., deriv= 1),
             6/c.^2* tukeyPsi1(x, c., deriv= 0)),
          EQ(tukeyChi(x, c., deriv= 2),
             6/c.^2* tukeyPsi1(x, c., deriv= 1)))

## Test if default arguments are used
tPsi <- chgDefaults(huberPsi, k = 2)

x <- 1:10
stopifnot(tPsi@ rho(x, k=2) == tPsi@ rho(x),
          tPsi@ psi(x, k=2) == tPsi@ psi(x),
          tPsi@Dpsi(x, k=2) == tPsi@Dpsi(x),
          tPsi@ wgt(x, k=2) == tPsi@ wgt(x))

## Test default arguments for E... slots
stopifnot(EQ(tPsi@Erho (), 0.49423127328548),
          EQ(tPsi@Epsi2(), 0.920536925636323),
          EQ(tPsi@EDpsi(), 0.954499736103642))

stopifnot(EQ(1, huberPsi@psi(1, k = 1e16)),
          huberPsi@wgt(0.1591319494080224, 0.5 + 1/13) <= 1)
## both used to fail because of numeric instability in pmin2/pmax2

