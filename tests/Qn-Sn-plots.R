library(robustbase)


### Back-compatibility check of consistency & finite sample correction for Qn() :

set.seed(153)
x <- sort(c(rnorm(80), rt(20, df = 1)))

ix <- c(27, 57, 15, 1, 26, 13, 23, 70, 9, 54, 6, 12, 8, 80, 11, 69,
        41, 53, 10, 37, 50, 21, 48, 51, 71, 17, 30, 16, 45, 28, 55, 5,
        59, 77, 61, 40, 63, 42, 72, 66)
QnA. <- c(0, 0.72307295, 1.2926498, 1.596857, 1.0979815, 0.84209457,
          1.0719335, 0.88620416, 1.0905118, 0.99056842, 1.2229216, 1.0626517,
          1.1738174, 1.1433873, 1.2071829, 1.1562513, 1.2182886, 1.1587793,
          1.1585524, 1.1555462, 1.1376428, 1.0532134, 1.0447343, 1.0200998,
          1.0495224, 1.0120569, 1.0094172, 0.9749928, 0.9530458, 0.92767184,
          0.90922667, 0.98987601, 0.98223857, 1.0053697, 0.98792848, 0.951908,
          0.92226488, 0.92312857, 0.92313406, 0.92733413)
QnAx <- sapply(seq_along(ix), function(n) Qn(x[ix[1:n]]))
stopifnot( all.equal(QnA., QnAx) )


### -------------- Plots -----------------------------------------------

if(!dev.interactive(orNone=TRUE)) pdf("Qn-Sn-plots.pdf")

plot(QnA., type="o", main="Qn(<random(n)>)")
abline(h = 1, lty=2)

n <- 1:50
(qnn <- sapply(n, function(n)Qn(1:n, const=1)))
plot(n, qnn, type = 'b', col = 2,
     ylab = "Qn", main = "Qn(1:n) [unscaled]")

(snn <- sapply(n, function(n)Sn(1:n, const=1)))
plot(n, snn, type = 'b', col = 2,
     ylab = "Sn", main = "Sn(1:n) [unscaled]")

matplot(n, cbind(qnn, snn), type = 'b',
        ylab = "Qn & Sn", main = "Qn(1:n) & Sn(1:n) [unscaled]")
legend("topleft", c("Qn", "Sn"), col=1:2, lty=1:2, bty="n", pch=paste(1:2))

(sdn <- c(1, sapply(n[-1], function(n)sd(1:n)/n)))
## sd(1) => NA

for(Sample in c(function(n) ppoints(n),
                function(n) qnorm(ppoints(n))))
{
    ##mult.fig(2) :
    op <- par(mfrow=c(2,1), mgp = c(1.5,.6,0), mar = .1 + c(4,4,2,1))
    for(N in c(50, 200)) {
        n <- 1:N
        sdn <- c(1, sapply(n[-1], function(m)sd(Sample(m))))
        r <- cbind(Qn = sapply(n, function(m)Qn(Sample(m))),
                   Sn = sapply(n, function(m)Sn(Sample(m)))) / sdn
        matplot(n, r, type = 'b', col = 2:3, lty = 1, ylab = "Qn & Sn",
                main = "Qn(Sample(n)) & Sn(..) [consistently scaled]")
        legend(.85*N, 0.4, c("Qn()", "Sn()"), col = 2:3, lty = 1,
               pch = paste(1:2))
        abline(h=1, col = "gray", lty = 2)
    }
    par(op)

    ## Hmm,  the above does not look 100% consistent to *my* eyes...
    ## Investigate:

    matplot(n, r, ylim = c(0.9, 1.1), type = 'b', col = 2:3, lty = 1)
    abline(h=1, col = "gray", lty = 2)

    matplot(n, r^2, ylim = c(0.7, 1.3), type = 'b', col = 2:3, lty = 1)
    abline(h=1, col = "gray", lty = 2)
}
rownames(r) <- paste(n)
r
