\name{fullRank}
\alias{fullRank}
\title{Remove Columns (or Rows) From a Matrix to Make It Full Rank}
\description{
  From the QR decomposition with pivoting, (\code{\link{qr}(x, tol)} if
  \eqn{n \ge p}), if
  the matrix is not of full rank, the corresponding columns (\eqn{n \ge
    p}{n >= p}) or rows (\eqn{n < p}) are omitted to form a full rank matrix.
}
\usage{% -> ../R/adjoutlyingness.R
fullRank(x, tol = 1e-7, qrx = qr(x, tol=tol))
}
\arguments{
  \item{x}{a numeric matrix of dimension \eqn{n \times p}{n * p}, or a
    similar object for which \code{\link{qr}()} works.}
  \item{tol}{tolerance for determining rank (deficiency).  Currently is
    simply passed to \code{\link{qr}}.}
  \item{qrx}{optionally may be used to pass a \code{\link{qr}(x, ..)};
    only used when \code{p <= n}.}
}
\value{
  a version of the matrix \code{x}, with less columns or rows if
  \code{x}'s rank was smaller than \code{min(n,p)}.

  If \code{x} is of full rank, it is returned unchanged.
}
\author{Martin Maechler}
\note{
  This is useful for robustness algorithms that rely on \eqn{X} matrices
  of full rank, e.g., \code{\link{adjOutlyingness}}.

  This also works for numeric data frames and whenever \code{qr()} works correctly.
}
\seealso{
  \code{\link{qr}}; for more sophisticated rank determination,
  \code{\link[Matrix]{rankMatrix}} from package \CRANpkg{Matrix}.
}
\examples{
stopifnot(identical(fullRank(wood), wood))

## More sophisticated and delicate
dim(T <- tcrossprod(data.matrix(toxicity))) # 38 x 38
dim(T. <- fullRank(T)) # 38 x 10
if(requireNamespace("Matrix")) {
  rMmeths <- eval(formals(Matrix::rankMatrix)$method)
  rT. <- sapply(rMmeths, function(.m.) Matrix::rankMatrix(T., method = .m.))
  print(rT.) # "qr" (= "qrLinpack"): 13,  others rather 10
}
dim(T.2 <- fullRank(T, tol = 1e-15))# 38 x 18
dim(T.3 <- fullRank(T, tol = 1e-12))# 38 x 13
dim(T.3 <- fullRank(T, tol = 1e-10))# 38 x 13
dim(T.3 <- fullRank(T, tol = 1e-8 ))# 38 x 12
dim(T.) # default from above          38 x 10
dim(T.3 <- fullRank(T, tol = 1e-5 ))# 38 x 10 -- still

plot(svd(T, 0,0)$d, log="y", main = "singular values of T", yaxt="n")
axis(2, at=10^(-14:5), las=1)
## pretty clearly indicates that  rank 10  is "correct" here.
}
\keyword{algebra}
\keyword{array}
