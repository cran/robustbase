require(robustbase)

## to check:
## - for the empty model
summary(lmrob(Y ~ 0, coleman))
## - with and without an intercept in the  model
summary(lmrob(Y ~ 1, coleman))
summary(lmrob(Y ~ ., coleman))
summary(lmrob(Y ~ ., coleman, model=FALSE, x=FALSE, y=FALSE))
summary(lmrob(Y ~ . - 1, coleman, model=FALSE, x=FALSE, y=FALSE))
## - when prior-weights are included
wts <- c(rep(0.05, 10), rep(2, 10))
summary(lmrob(Y ~ . - 1, coleman, model=FALSE, x=FALSE, y=FALSE,
              weights = wts))
## - should work for object with NA in the coefficients (done NAcoef.R)
## - should work for object is NA in the observations (done NAcoef.R)

## check equality with lm() for classical model
## FIXME: SE column of coefficients does not match (.vcov.w)
##        because corrfact does not reduce to 1.
test <- function(formula, data,
                 items=c(## "coefficients", "residuals", "df", "scale",
                         "r.squared", "adj.r.squared"),
                 tolerance = 1e-4, ...) {
    ## FIXME: also include weights
    sc <- summary(lm(formula, data))
    sr <- summary(lmrob(formula, data,
                        control=lmrob.control(psi = "hampel",
                            tuning.psi = c(1000, 2000, 3000),
                            method="SMDM", ...)))
    names(sc)[names(sc) == "sigma"] <- "scale"
    ret <- all.equal(sc[items], sr[items], tolerance=tolerance)
    if (!isTRUE(ret)) {
        print(sr)
        for (i in seq_along(items)) {
            print(sc[items[i]])
            print(sr[items[i]])
        }
        print(ret)
        stop("all.equal(sc[items], sr[items], tolerance = tolerance) are not all TRUE")
    }
    ret
}

test(Y ~ 0, coleman, c(## "residuals", "df", "coefficients",
                       "r.squared", "adj.r.squared"))
test(Y ~ 1, coleman)
test(Y ~ ., coleman)
test(Y ~ . - 1, coleman)
