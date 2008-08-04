#### -*- mode: R; kept-new-versions: 30; kept-old-versions: 20 -*-

#### MC-adjusted Outlyingness
#### ------------------------
###
### Original code from  the web site from the Antwerpen Statistics Group :
###  http://www.agoras.ua.ac.be/Robustn.htm
### which has a "MC" section and for the software links to
### ftp://ftp.win.ua.ac.be/pub/software/agoras/newfiles/mc.tar.gz
### and that contains  mcrsoft/adjoutlyingness.R

## MM [ FIXME ]:
## -----------

## 1)   Use  *transposed*  B[] and A[] matrices   -- done

## 2)   use  IQR() instead of   quantile(., .75) - quantile(., .25)

##-->  but only *after* testing original code
##     ^^^^^^^^^^^^^^^^^^^^^^^^

adjOutlyingness <- function(x, ndir=250, clower=3, cupper=4,
                            alpha.cutoff = 0.75, coef = 1.5, qr.tol = 1e-12)
## Skewness-Adjusted Outlyingness
{
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    stopifnot(n >= 1, p >= 1, is.numeric(x))
    if (p < n) {
        i <- 1
        it <- 0
        B <- matrix(0, p, ndir)
        E <- matrix(1, p, 1)
        x. <- unname(x) # for speed in subsequent subsetting and solve
        maxit <- 100 * ndir
        ## ^^ original code had 'Inf', i.e. no iter.count check;
        ## often, maxit == ndir would suffice
        while (i <= ndir  &&  (it <- it+1) < maxit) {
            P <- x.[sample(n,p), , drop=FALSE]
            if ((qrP <- qr(P, tol = qr.tol))$rank == p) {
                B[,i] <- solve(qrP, E, tol =  qr.tol)
                i <- i+1
            }
        }
        if(it == maxit)
            stop("**** sampling iterations were not sufficient. Please report")

        Bnorm <- sqrt(colSums(B^2))
        ## MM: FIXME '1e-12' is not ok scale-equivariant!
        Bnormr <- Bnorm[  Bnorm > 1e-12]
        B      <-     B[, Bnorm > 1e-12 , drop=FALSE]

        A <- B / rep(Bnormr, each = nrow(B))
    }
    else {
        stop('More dimensions than observations: not yet implemented')
        ## MM: In LIBRA(matlab) they have it implemented:
        ##    seed=0;
        ##    nrich1=n*(n-1)/2;
        ##    ndirect=min(250,nrich1);
        ##    true = (ndirect == nrich1);
        ##    B=extradir(x,ndir,seed,true); %n*ri
        ##      ======== % Calculates ndirect directions through
        ##               % two random choosen data points from data
        ##    for i=1:size(B,1)
        ##        Bnorm(i)=norm(B(i,:),2);
        ##    end
        ##    Bnormr=Bnorm(Bnorm > 1.e-12); %ndirect*1
        ##    B=B(Bnorm > 1.e-12,:);       %ndirect*n
        ##    A=diag(1./Bnormr)*B;         %ndirect*n

    }
    Y <- x %*% A # (n x p) %*% (p, ndir) == (n x ndir)

    ## Compute and sweep out the median
    med <- apply(Y, MARGIN = 2, median)
    Y <- Y - rep(med, each=n)
    ## MM: mc() could be made faster if we could tell it that med(..) = 0
    tmc <- apply(Y, MARGIN = 2, mc)
    ##                          ==
    Q3 <-  apply(Y, MARGIN = 2, quantile, 0.75)
    Q1 <-  apply(Y, MARGIN = 2, quantile, 0.25)
    IQR <- Q3-Q1
    ## NOTA BENE(MM): simplified definition of tup/tlo here and below
    tup <- Q3 + coef*IQR*exp( cupper*tmc*(tmc >= 0) - clower*tmc*(tmc < 0))
    tlo <- Q1 - coef*IQR*exp(-clower*tmc*(tmc >= 0) + cupper*tmc*(tmc < 0))
    ## Note: all(tlo < med & med < tup)

    ## Instead of the loop:
    ##  for (i in 1:ndir) {
    ##      tup[i] <-  max(Y[Y[,i] < tup[i], i])
    ##      tlo[i] <- -min(Y[Y[,i] > tlo[i], i])
    ##      ## MM:  FIXED typo-bug :       ^^^ this was missing!
    ##      ## But after the fix, the function stops "working" for longley..
    ##      ## because tlo[] becomes 0 too often, YZ[.,.] = c / 0 = Inf !
    ##  }
    Yup <- Ylo <- Y
    Yup[!(Y < rep(tup, each=n))] <- -Inf
    Ylo[!(Y > rep(tlo, each=n))] <-  Inf
    tup <-  apply(Yup, 2, max) # =  max{ Y[i,] ; Y[i,] < tup[i] }
    tlo <- -apply(Ylo, 2, min) # = -min{ Y[i,] ; Y[i,] > tlo[i] }

    tY <- t(Y)
    Ypos <- (tY >= 0) ## note that all column-wise medians are 0
    ## Note: this loop is pretty fast
    for (j in 1:n)
        tY[, j] <- abs(tY[,j]) / (Ypos[,j]*tup + (1 - Ypos[,j])*tlo)
    ## FIXME -- what if denominator is 0 ? happens often in small samples
    ##	 e.g in  set.seed(3); adjOutlyingness(longley)
    ## even  have  0/0 -> NaN there  --> is.finite(.) below.. hmm, FIXME!

    adjout <- apply(tY, 2, function(x) max(x[is.finite(x)]))
    Qadj <- quantile(adjout, probs = c(1 - alpha.cutoff, alpha.cutoff))
    mcadjout <- mc(adjout)
    ##          ===
    cutoff <- Qadj[2] + coef* (Qadj[2] - Qadj[1])*
        (if(mcadjout > 0) exp(cupper*mcadjout) else 1)

    list(adjout = adjout,
	 MCadjout = mcadjout, Qalph.adjout = Qadj, cutoff = cutoff,
	 nonOut = (adjout <= cutoff))
}
