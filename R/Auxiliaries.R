### Fast versions of pmin() and pmax() for 2 arguments only:

### FIXME: should rather add these to R
pmin2 <- function(k,x) (x+k - abs(x-k))/2
pmax2 <- function(k,x) (x+k + abs(x-k))/2

