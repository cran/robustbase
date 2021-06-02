#### -*- mode: R; kept-new-versions: 30; kept-old-versions: 20 -*-

#### MC-adjusted Outlyingness
#### ------------------------
###
### Original code from  the web site from the Antwerpen Statistics Group :
###  http://www.agoras.ua.ac.be/Robustn.htm
### which has a "MC" section and for the software links to
### ftp://ftp.win.ua.ac.be/pub/software/agoras/newfiles/mc.tar.gz
### and that contains  mcrsoft/adjoutlyingness.R

##_NEW_ (2014): moved from Antwerpen to Leuwen,
## ===> http://wis.kuleuven.be/stat/robust/software
##      has several links to 'robustbase', and S-plus code
## http://wis.kuleuven.be/stat/robust/Programs/adjusted-boxplot/adjusted-boxplot.ssc
## (copy in ../misc/Adjusted-boxplot.ssc

## MM [ FIXME ]:
## -----------

## 1)   Use  *transposed*  B[] and A[] (now called 'E') matrices   -- DONE
## 2)   use  IQR() instead of   quantile(., .75) - quantile(., .25)

##-->  but only *after* testing original code
##     ^^^^^^^^^^^^^^^^^^^^^^^^


adjOutlyingness <- function(x, ndir = 250, p.samp = p, clower=4, cupper=3,
                            IQRtype = 7,
                            alpha.cutoff = 0.75, coef = 1.5,
                            qr.tol = 1e-12, keep.tol = 1e-12,
                            only.outlyingness = FALSE, maxit.mult = max(100, p),
                            trace.lev = 0,
                            ## these are all passed to mc() {when applied to the projected data}
                            mcReflect = n <= 100, mcScale = TRUE, mcMaxit = 2*maxit.mult,
                            mcEps1 = 1e-12, mcEps2 = 1e-15,
                            mcTrace = max(0, trace.lev-1))
## Skewness-Adjusted Outlyingness
{
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n >= 1, p >= 1, p.samp >= p, is.numeric(x))
    if (p <= n) {
        B <- matrix(0, p, ndir)
        E <- matrix(1, p.samp, 1)
        x. <- unname(x) # for speed in subsequent subsetting and solve
        maxit <- as.integer(maxit.mult * ndir)
        ## ^^ original code had 'Inf', i.e. no iter.count check;
        ## often, maxit == ndir would suffice
        if(trace.lev >= 2) p10 <- 10 ^ max(0, min(6 - trace.lev, floor(log10(maxit))))
	i <- 1L
	it <- 0L
	while (i <= ndir  &&  (it <- it+1L) < maxit) {
            ## we sort to get *identical* projections instead of "almost"
            P <- x.[sort(sample.int(n, p.samp)), , drop=FALSE]
            if ((qrP <- qr(P, tol = qr.tol))$rank == p) {
                B[,i] <- solve(qrP, E, tol = qr.tol)
                ## if(trace.lev >= 2) cat(" it=",it,"; found direction # ", i,"\n", sep="")
		i <- i+1L
            } else if(trace.lev >= 2) {
                if(it %% p10 == 0)
                    cat(" it=",it,": rank(qr(P ..)) = ", qrP$rank, " < p = ",p,"\n", sep="")
            }
        }
	if(it == maxit) {
	    rnk.x <- qr(x, tol = qr.tol)$rank
	    if(rnk.x < p)
		stop("Matrix 'x' is not of full rank: rankM(x) = ",rnk.x," < p = ",p,
                     "\n Use fullRank(x) instead")
	    ## else
	    stop("** direction sampling iterations were not sufficient. Maybe try increasing 'maxit.mult'")
	}
        Bnorm <- sqrt(colSums(B^2))
        Nx <- mean(abs(x.)) ## so the comparison is scale-equivariant:
        keep <- Bnorm*Nx > keep.tol
        if(!all(keep)) {
            if(trace.lev)
                cat("keeping ", sum(keep), "(out of", length(keep),") normalized directions\n")
            Bnorm <- Bnorm[  keep ]
            B     <-     B[, keep , drop=FALSE]
        } else if(trace.lev) cat("keeping *all* ",length(keep)," normalized directions\n")

        B <- B / rep(Bnorm, each = nrow(B)) # normalized B {called 'A' in orig.code}
    }
    else {
        stop('More dimensions than observations: not yet implemented')
        ## MM: In LIBRA(matlab) they have it implemented:
        ##    seed=0;
        ##    nrich1=n*(n-1)/2;
        ##    ndirect=min(250,nrich1);
        ##    true = (ndirect == nrich1);
        ##    B=extradir(x,ndir,seed,true); % n*ri
        ##      ======== % Calculates ndirect directions through
        ##               % two random choosen data points from data
        ##    for i=1:size(B,1)
        ##        Bnorm(i)=norm(B(i,:),2);
        ##    end
        ##    Bnorm=Bnorm(Bnorm > 1.e-12); %ndirect*1
        ##    B=B(Bnorm > 1.e-12,:);       %ndirect*n
        ##    A=diag(1./Bnorm)*B;         %ndirect*n

    }
    ## NB:  colSums( B^2 ) == 1
    Y <- x %*% B # (n x p) %*% (p, nd') == (n x nd');
    ##  nd' = ndir.final := ndir - {those not in 'keep'}

    ## Compute and sweep out the median
    med <- colMedians(Y)
    if(trace.lev) {
        cat("med <- colMedians(Y): ", length(med), " values; summary():\n")
        print(summary(med))
    }
    Y <- Y - rep(med, each=n)
    ## central :<==> non-adjusted  <==> "classical" outlyingness
    central <- clower == 0 && cupper == 0
    if(!central) {
        ## MM: mc() could be made faster if we could tell it that med(..) = 0
        ##                  vv
        tmc <- apply(Y, 2L, mc, doReflect=mcReflect, doScale=mcScale,
                     maxit=mcMaxit, eps1=mcEps1, eps2=mcEps2, trace.lev = mcTrace)
        ## original Antwerpen *wrongly*: tmc <- mc(Y)
        if(trace.lev) {
            cat("Columnwise mc() got ", length(tmc), " values; summary():\n")
            print(summary(tmc))
        }
    }
    Q13 <- apply(Y, 2, quantile, c(.25, .75), names=FALSE, type = IQRtype)
    Q1 <- Q13[1L,]; Q3 <- Q13[2L,]
    IQR <- Q3 - Q1
    ## NOTA BENE(MM): simplified definition of tup/tlo here and below
    ## 2014-10-18: "flipped" sign (which Pieter Setaert (c/o Mia H) proposed, Jul.30,2014:
    tup <- Q3 + coef*
	(if(central) IQR else IQR*exp( cupper*tmc*(tmc >= 0) + clower*tmc*(tmc < 0)))
    tlo <- Q1 - coef*
	(if(central) IQR else IQR*exp(-clower*tmc*(tmc >= 0) - cupper*tmc*(tmc < 0)))
    ## Note: all(tlo < med & med < tup) # where med = 0
    if(trace.lev >= 3) {
        print( cbind(tlo, Q1, Q3, tup) )
    }

    ## Instead of the loop:
    ##  for (i in 1:ndir) {
    ##      tup[i] <-  max(Y[Y[,i] < tup[i], i])
    ##      tlo[i] <- -min(Y[Y[,i] > tlo[i], i])
    ##      ## MM:  FIXED typo-bug :       ^^^ this was missing!
    ##      ## But after the fix, the function stops "working" for longley..
    ##      ## because tlo[] becomes 0 too often, YZ[.,.] = c / 0 = Inf !
    ##  }
    Yup <- Ylo <- Y
    Yup[Y >= rep(tup, each=n)] <- -Inf
    Ylo[Y <= rep(tlo, each=n)] <-  Inf
    y.up <-  apply(Yup, 2, max) # =  max{ Y[i,] ; Y[i,] < tup[i] }
    y.lo <- -apply(Ylo, 2, min) # = -min{ Y[i,] ; Y[i,] > tlo[i] }
    if(trace.lev) {
        cat(length(y.up), "lower & upper Y (:= X - med(.)) values:\n")
        print(summary(y.lo))
        print(summary(y.up))
    }
    tY <- t(Y)
    ## Note: column-wise medians are all 0 : "x_i > m" <==> y > 0
    ## Note: this loop is pretty fast
    for (j in 1:n) { # when y = (X-med) = 0  ==> adjout = 0 rather than
	## 0 / 0 --> NaN; e.g, in  set.seed(3); adjOutlyingness(longley)
	non0 <- 0 != (y <- tY[,j]); y <- y[non0]; I <- (y > 0)
	D <- I*y.up[non0] + (1 - I)*y.lo[non0]
        if(trace.lev >= 3) {
            cat(sprintf("j=%2d: #{non0}= %2d; quantile(D)=\n", j, sum(non0)))
            print(quantile(D), digits=3)
        }
	tY[non0, j] <- abs(y) / D
    }
    ## We get +Inf above for "small n"; e.g. set.seed(11); adjOutlyingness(longley)
    if(trace.lev) {
        cat("outlyingnesses for all directions (of which max(.) will be chosen:\n")
        print(quantile(tY, digits=3))
    }
    adjout <- apply(tY, 2, function(x) max(x[is.finite(x)]))
    ##----                             ---
    if(abs(trace.lev %% 1  - 0.7) < 1e-3) { ## really not for the end user ..
        cat("Plotting outlyingnesses vs. observation i:\n")
        matplot(t(tY), log="y");           axis(2, at=1); abline(h=1, lty=3)
        Sys.sleep(2)
        ## somewhat revealing: 3  groups:  very large |  medium | very small (incl. 0 which are *not* plotted)
        matplot(t(tY), log="y", type="b"); axis(2, at=1); abline(h=1, lty=3)
        browser()  ## <<<<<<<<<<<<<<<<<<<
    }

    if(only.outlyingness)
	adjout
    else {
	Qadj <- quantile(adjout, probs = c(1 - alpha.cutoff, alpha.cutoff))
	mcadjout <- if(cupper != 0) mc(adjout, doScale=mcScale, eps1=mcEps1, eps2=mcEps2) else 0
	##			    ===
	cutoff <- Qadj[2] + coef * (Qadj[2] - Qadj[1]) *
	    (if(mcadjout > 0) exp(cupper*mcadjout) else 1)

	list(adjout = adjout, iter = it, ndir.final = sum(keep),
	     MCadjout = mcadjout, Qalph.adjout = Qadj, cutoff = cutoff,
	     nonOut = (adjout <= cutoff))
    }
}

##' Compute a "full rank" version of matrix x,
##' by removing columns (or rows when nrow(x) < ncol(x)), using qr() and it's pivots
fullRank <- function(x, tol = 1e-7, qrx = qr(x, tol=tol)) {
    d <- dim(x)
    n <- d[[1L]]; p <- d[[2L]]
    if(n < p)
        return( t(fullRank(t(x), tol=tol)) )
    ## else n >= p >= rank(.)
    rnk <- qrx$rank
    if(rnk == p)
        x
    else
        x[, qrx$pivot[seq_len(rnk)], drop=FALSE]
}

