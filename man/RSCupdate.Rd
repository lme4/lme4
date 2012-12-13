\name{RSCupdate}
\alias{RSCupdate}
\title{Update for the penalized least squares problem}
\usage{
  RSCupdate(rv, xv, theta, lower, resid, m, ZtXty)
}
\arguments{
  \item{rv}{the matrix of row indices for the regular
  sparse column representation of ZtXt}

  \item{xv}{the non-zero values in ZtXt}

  \item{theta}{the covariace parameter vector}

  \item{lower}{lower bounds on the covariace parameter
  vector}

  \item{resid}{current residual}

  \item{m}{the sparse matrix to be updated.  Must have the
  correct pattern.}

  \item{ZtXty}{the product of the updated LambdatZt and Xt
  with resid.}
}
\description{
  Update for the penalized least squares problem
}
\examples{
fm1 <- lmer(Yield ~ 1|Batch, Dyestuff, REML=FALSE)
rv  <- rbind(with(Dyestuff, as.integer(Batch) - 1L), 6L)
xv  <- matrix(1, nrow=nrow(rv), ncol=ncol(rv))
A   <- tcrossprod(sparseMatrix(i = as.vector(rv),
                              j = as.vector(col(rv) - 1L),
                              x = as.vector(xv), index1=FALSE))
ubeta <- numeric(ncol(A))
with(Dyestuff, RSCupdate(rv[-nrow(rv),,drop=FALSE], xv, getME(fm1,"theta"), 0., Yield, A, ubeta))
L <- Cholesky(A, perm=FALSE, LDL=FALSE)
vv <- solve(L, ZtXty, system="A")
all.equal(vv[1:6], getME(fm1, "u"))
all.equal(vv[7], getME(fm1, "beta"))
LL <- as(L, "Matrix")
## the transpose below is because a diagonal dtCMatrix is declared upper triangular
all.equal(t(LL[1:6,1:6]), as(getME(fm1, "L"), "Matrix"))
all.equal(LL[7,7], as.vector(getME(fm1, "RX")))
all.equal(LL[7,1:6], as.vector(getME(fm1, "RZX")))
}

