## trimmed mad := trimmed *[M]ean* of [A]bsolute [D]eviations (from the median; more generally 'center')
## =       ===    {no default for trim on purpose}
tmad <- function(x, center = median(x), trim, na.rm = FALSE)
{
    stopifnot(is.numeric(trim), length(trim) == 1L, 0 <= trim, trim <= 0.5)
    if(na.rm)
        x <- x[!is.na(x)]
    ## TODO: consistency correction (for non-large) 'n' as a function of trim
    ## ----  not needed for huberize() though
    ## n <- length(x)
    mean(abs(x - center), trim=trim)
}

## Estimates mu (optionally, via huberM(.)) and sigma of x
##   x: without NA: na.rm=TRUE must have happened
## sets boundaries at  M +/- c*sigma
## sets outliers to be equal to lower/upper boundaries
huberize <- function(x, M = huberM(x, k=k)$mu, c = k,
                     trim = (5:1)/16, # Lukas Graz' MSc thesis had c(0.25, 0.15, 0.075)
                     k = 1.5,
                     warn0 = getOption("verbose"), saveTrim = TRUE)
{
    stopifnot(is.numeric(M), length(M) == 1,
              length(trim) >= 1, trim >= 0, diff(trim) < 0) # trim must be strictly decreasing
    qn. <- Qn(x)
    j <- 0L
    while(!is.na(qn.) && qn. == 0 && j < length(trim))
        qn. <- tmad(x, center = M, trim = trim[j <- j+1L])
        ##     ~~~~
    if(qn. == 0 && warn0)
        warning(sprintf("Qn(x) == 0 and tmad(x, trim=%g) == 0", trim[j]))
    upper <- M + qn.*c   # qnorm(c,lower.tail = F)
    lower <- M - qn.*c
    x[x > upper] <- upper
    x[x < lower] <- lower
    ## store the final 'trim' used (if there was one) as attribute:
    if(j && saveTrim) structure(x, trim = trim[j]) else x
}
