#### This is from the R package
####
####  rrcov : Scalable Robust Estimators with High Breakdown Point
####
#### by Valentin Todorov

### This program is free software; you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation; either version 2 of the License, or
### (at your option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this program; if not, write to the Free Software
### Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

##   I would like to thank Peter Filzmoser for providing the initial code of
##   this function.

plot.mcd <-
    function(x,
             which=c("all", "dd","distance","qqchi2","tolEllipsePlot","screeplot"),
             classic= FALSE,
             ask = (which=="all" && dev.interactive()),
             cutoff,
             id.n,
             tol = 1e-7, ...)
{  ## VT:: 16.04.2005 - change for 2.1.0 - use 'tol' instead of 'tol.inv'
    if (!inherits(x, "mcd"))
        stop("Use only with 'mcd' objects")

    covPlot(x$X, which=which, classic=classic, ask=ask,
            cutoff=cutoff, mcd=x, id.n=id.n, tol=tol, ...)
}

covPlot <- function(x,
                   which=c("all", "dd","distance","qqchi2","tolEllipsePlot","screeplot"),
                   classic=FALSE,
                   ask=FALSE,
                   mcd,
                   cutoff,
                   id.n,
                   tol = 1e-7, ...)
{

    ##@bdescr
    ##  Make plots based on the covariance structure of a data set:
    ##  dd       -  distance-distance plot: Robust distances versus
    ##              Mahalanobis distances
    ##  distance -  a plot of the robust distances
    ##  qqchi2   -  a qq-plot of the robust distances versus the
    ##              quantiles of the chi-squared distribution
    ##  tolEllipsePlot- a tolerance ellipse
    ##  screeplot-  a screeplot of the eigenvalues ov the covariance matrix
    ##
    ## Distance Plot:
    ## Draw a Distance-Distance Plot: Plots the robust distances
    ## versus the classical Mahalanobis distances as introduced by
    ## Rousseeuw, P. J., and van Zomeren, B. C. (1990). Unmasking
    ## Multivariate Outliers and Leverage Points. Journal of the American
    ## Statistical Association, 85, 633-639.
    ##
    ## The dashed line is the set of points where the robust distance is
    ## equal to the classical distance.
    ## The horizontal and vertical dotted lines are drawn at values equal cutoff
    ## which defaults to square root of the 97.5% quantile of a chi-squared
    ## distribution with p degrees of freedom. Points beyond these lines can
    ## be considered outliers.
    ##
    ##@edescr
    ##
    ##@in  x                 : [matrix] A data.frame or matrix, n > 2*p
    ##@in  which          : [character] A plot option, one of:
    ##                            classic: index plot of the classical mahalanobis distances
    ##                            robust:  index plot of the robust mahalanobis distances
    ##                            dd:      distance-distance plot
    ##                            index:   parallel index plot of classical and robust distances
    ##                            all:     all three plots
    ##                          default is "all"
    ##@in  classic           : [logical] If true the classical plot will be displayed too
    ##                                   default is classic=FALSE
    ##@in  mcd               : [mcd object] An object of type mcd - its attributes
    ##                                      center and cov will be used
    ##@in  cutoff            : [number] The cutoff value for the distances
    ##@in  id.n              : [number] number of observations to be identified with a label.
    ##                                  Defaults to the number of observations with distance
    ##                                  larger than cutoff
    ##@in  tol               : [number] tolerance to be used for computing the inverse - see 'solve'.
    ##                                  defaults to 1e-7

    ## NOTE: The default tolerance 1e-7, will not work for some example
    ##       data sets, like milk or aircraft

    myscreeplot <- function(x, mcd) {

        if(missing(mcd))
            mcd <- covMcd(x)
        erob <- eigen(mcd$cov,symmetric=TRUE,only.values=TRUE)$values
        eclass <- eigen(var(x),symmetric=TRUE,only.values=TRUE)$values

        leg.txt <- c("Robust", "Classical")
        leg.col <- c("green", "red")
        leg.pch = c(1,24)
        leg.lty = c("solid", "dotted")

        eall <- c(erob,eclass)
        ylim <- c( min(eall), max(eall))

        plot(erob, ylim=ylim, ylab="Eigenvalues", xlab="Index", type="n")
        legend("topright", leg.txt, pch = leg.pch, lty = leg.lty, col = leg.col)

        lines(erob, type="o", pch = leg.pch[1], lty = leg.lty[1], col=leg.col[1])
        lines(eclass, type="o", pch = leg.pch[2], lty = leg.lty[2], col=leg.col[2])

        title(main="Scree plot")
    }

    mydistplot <- function(x, cutoff, classic = FALSE, id.n) {
        ##  Index Plot:
        ##  Plot the vector x (robust or mahalanobis distances) against
        ##  the observation indexes. Identify by a label the id.n
        ##  observations with largest value of x. If id.n is not supplied,
        ##  calculate it as the number of observations larger than cutoff.
        ##  Use cutoff to draw a horisontal line.
        ##  Use classic=FALSE/TRUE to choose the label of the vertical axes

        n <- length(x)
        if(missing(id.n))
            id.n <- length(which(x>cutoff))
        if(classic)
            ylab="Square Root of Mahalanobis distance"
        else
            ylab="Square Root of Robust distance"
        plot(x, ylab=ylab, xlab="Index", type="p")
        label(1:n, x, id.n)
        abline(h=cutoff)

        title(main="Distance Plot")
    }

    myddplot <- function(md, rd, cutoff, id.n) {
        ##  Distance-Distance Plot:
        ##  Plot the vector y=rd (robust distances) against
        ##  x=md (mahalanobis distances). Identify by a label the id.n
        ##  observations with largest rd. If id.n is not supplied, calculate
        ##  it as the number of observations larger than cutoff. Use cutoff
        ##  to draw a horisontal and a vertical line. Draw also a dotted line
        ##  with a slope 1.
        n <- length(md)
        if(missing(id.n))
            id.n <- length(which(rd>cutoff))
        xlab <- "Mahalanobis distance"
        ylab <- "Robust distance"
        plot(md, rd, xlab=xlab, ylab=ylab, type="p")
        label(md,rd,id.n)
        abline(0, 1, lty=2)
        abline(v=cutoff)
        abline(h=cutoff)

        title(main="Distance-Distance Plot")
    }

    qqplot <- function(x, p, cutoff, classic=FALSE, id.n) {
        ##  Chisquare QQ-Plot:
        ##  Plot the vector x (robust or mahalanobis distances) against
        ##  the square root of the quantiles of the chi-squared distribution
        ##  with p degrees of freedom.
        ##  Identify by a label the id.n observations with largest value of x.
        ##  If id.n is not supplied, calculate it as the number of observations
        ##  larger than cutoff.
        ##  Use classic=FALSE/TRUE to choose the label of the vertical axes


        ##  parameters and preconditions

        n <- length(x)

        if(missing(cutoff))
            cutoff <- sqrt(qchisq(0.975, p))

        if(missing(id.n))
            id.n <- length(which(x>cutoff))

        qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

        x <- sort(x, index.return=TRUE)
        ix <- x$ix
        x <- x$x

        if(classic)
            ylab="Mahalanobis distance"
        else
            ylab="Robust distance"

        plot(qq, x, xlab="Square root of the quantiles of the chi-squared distribution", ylab=ylab, type="p")
        if(id.n > 0) {
            ind <- (n-id.n+1):n
            xrange <- par("usr")
            xrange <- xrange[2] - xrange[1]
            text(qq[ind] + xrange/50, x[ind], ix[ind])
        }
        abline(0, 1, lty=2)
        title(main="Chisquare QQ-Plot")
    }

    label <- function(x, y, id.n=3) {
        if(id.n > 0) {
            xrange <- par("usr")
            xrange <- xrange[2] - xrange[1]
            n <- length(y)
            ind <- sort(y, index.return=TRUE)$ix
            ind <- ind[(n-id.n+1):n]
            text(x[ind] + xrange/50, y[ind], ind)
        }
    }

    ##  parameters and preconditions

    if(is.vector(x) || is.matrix(x)) {
        if(!is.numeric(x))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(x)) {
        if(!all(sapply(x,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(x)[1]
    p <- dim(x)[2]

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    if(missing(mcd))
        mcd <- covMcd(x)

    if(length(mcd$center)  == 0 ||  length(mcd$cov) == 0)
        stop( "Invalid mcd object: attributes center and cov missing!")

    if(length(mcd$center)  != p)
        stop( "Data set and provided center have different dimensions!")

    if(mcd$crit == 0)
        stop( "The covariance matrix is singular!")

    md <- mahalanobis(x, colMeans(x), var(x), tol=tol)
    md <- sqrt(md)

    rd <- mahalanobis(x, mcd$center, mcd$cov, tol=tol)
    rd <- sqrt(rd)

    which <- match.arg(which)

    ## *Never* here : par(mfrow=c(1,1), pty="m")
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    if(which == "all" || which == "distance") {
        if(classic) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        mydistplot(rd, cutoff, id.n=id.n) # index plot of mahalanobis distances
        if(classic) {
            mydistplot(md, cutoff, classic=TRUE, id.n=id.n) # index plot of robust distances
            par(opr)
        }
    }

    if(which == "all" || which == "dd") {
        myddplot(md, rd, cutoff=cutoff, id.n=id.n) # distance-distance plot
    }

    if(which == "all" || which == "qqchi2") {
        if(classic) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        qqplot(rd, p, cutoff=cutoff, id.n=id.n) # qq-plot of the robust distances versus the
                                        # quantiles of the chi-squared distribution
        if(classic) {
            qqplot(md, p, cutoff=cutoff, classic=TRUE, id.n=id.n)
                                        # qq-plot of the mahalanobis distances
            par(opr)
        }
    }

    if(which == "all" || which == "tolEllipsePlot") {
        if(length(dim(x)) >= 2 && dim(x)[2] == 2)
            tolEllipsePlot(x, mcd=mcd, cutoff=cutoff,
                           id.n=id.n, classic=classic, tol=tol)
        else
            warning("Warning: For tolerance ellipses the dimension must be 2!")
    }

    if(which == "all" || which == "screeplot") {
        myscreeplot(x, mcd=mcd)
    }

} ## end { covPlot }

## ddplot <- function(x,...) {
##     covPlot(x, which="dd", ...)
## }

## distplot <- function(x,...) {
##     covPlot(x, which="distance", ...)
## }

## chi2qqplot <- function(x,...) {
##     covPlot(x, which="qqchi2", ...)
## }

## ellipse() exists in other packages
## ellipse <- function(x,...) {
##     covPlot(x, which="tolEllipsePlot", ...)
## }
