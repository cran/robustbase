require(robustbase)

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
