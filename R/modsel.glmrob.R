
modsel.glmrob <-
    function(full.mfit, reduced.mfit, quad.form = FALSE,
             test = c("Wald", "Quasi-Deviance"))
{
    if (!inherits(full.mfit, "glmrob"))
        stop("currently only works for 'glmrob' objects")
    test <- match.arg(test)
    X <- model.matrix(full.mfit)
    asgn <- attr(X, "assign")
    tt <- terms(full.mfit)
    tt0 <- terms(reduced.mfit)
    tl <- attr(tt, "term.labels")
    tl0 <- attr(tt0, "term.labels")
    numtl0 <- match(tl0 , tl, nomatch = -1)
    if(attr(tt0, "intercept") == 1) numtl0 <- c(0, numtl0)
    if(any(is.na(match(numtl0, unique(asgn)))))
        stop("Models are not nested!")
    mod0 <- seq(along = asgn)[!is.na(match(asgn, numtl0))]
    if (length(asgn) == length(mod0))
        stop("Models are not strictly nested")
    H0ind <- setdiff(seq(along = asgn), mod0)
    H0coef <- coef(full.mfit)[H0ind]
    df <- length(H0coef)
    ## null.value {same names, etc}:
    c0 <- H0coef ; c0[] <- 0

    if(test == "Wald") {
        t.cov <- full.mfit$cov
        t.chisq <- sum(H0coef * solve(t.cov[H0ind, H0ind], H0coef))
        statistic <- c(chisq = t.chisq)
        method <- "robust Wald test"
    }
    else if(test == "Quasi-Deviance") {

        matM <- full.mfit$matM
        if(quad.form) {
            ## Difference of robust quasi-deviances
            ## via the asymptotically equivalent quadratic form
            matM11 <- matM[mod0,mod0]
            matM12 <- matM[mod0,H0ind]
            matM22 <- matM[H0ind, H0ind]
            matM22.1 <- matM22-t(matM12) %*% solve(matM11) %*% matM12
            Dquasi.dev <- nrow(X) * c(H0coef %*% matM22.1 %*% H0coef)
        }
        else {
            quasiDev <- switch(full.mfit$family$family,
                               poisson = glmrobMqleDiffQuasiDevPois,
                               binomial =  glmrobMqleDiffQuasiDevB,
                               stop("This family is not implemented"))

            ## note that qdev and qdev0 do depend on an incorrectly specified
            ## lower limits in the integration. But this does't matter in
            ## the following difference, because the difference does not
            ## deepend on it! (Hence I could use the centered nui
            ## (cnui= nui - Enui) in quasiDev as the function to be integrated.
            Dquasi.dev <- quasiDev(mu = full.mfit$fitted.values,
                                   mu0 = reduced.mfit$fitted.values,
                                   y = full.mfit$y, ni = full.mfit$ni,
                                   w.x = full.mfit$w.x,
                                   tcc = full.mfit$tcc)
        }
        ## Asymptotic distribution: variance and weights of the sum of chi2
        matQ <- full.mfit$matQ
        matM11inv <- solve(matM[mod0,mod0])
        pp <- df + length(mod0)
        Mplus <- matrix(0, ncol = pp, nrow = pp)
        Mplus[mod0, mod0] <- matM11inv

        d.ev <- Re(eigen(matQ %*% (solve(matM)-Mplus))$values)
        d.ev <- d.ev[1:df] ## just the q (=df) lagest eigenvalues are needed
        if(any(d.ev < 0)) warning("not all eigenvalues are non-negative")

        ## p-value: exact computation for q=1, approximated for q>1 (q=df)
        statistic <- c(quasi.dev = Dquasi.dev/mean(d.ev))
        method <- "Robust Quasi-Deviance Test"

    } else stop("non-implemented test method:", test)

    ## return
    structure(list(statistic = statistic, df = df,
                   data.name = paste("from", deparse(full.mfit$call)),
                   method = method,
                   alternative = "two.sided", null.value = c0,
                   p.value =
                   pchisq(as.vector(statistic), df = df, lower.tail = FALSE)),
              class = "htest")
}

