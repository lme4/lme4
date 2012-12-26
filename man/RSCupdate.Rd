\name{RSCupdate}
\alias{RSCupdate}
\title{Update for the penalized least squares problem}
\usage{
  RSCupdate(pred, resid)
}
\arguments{
  \item{pred}{an RSC object}

  \item{resid}{current residual}
}
\value{
  the solution of the PLS problem as matrix of 1 column
}
\description{
  Update for the penalized least squares problem
}
\details{
  This function has a side-effect of updating both the
  sparse symmetric system matrix in the \code{A} slot and
  its Cholesky factors in its \code{factors} slot.  Note
  that it updates in place, thus violating the functional
  language requirement of not changing arguments.  This is
  to conserve memory.  Usually the sparse Cholesky factor
  is the largest object in the model representation and it
  doesn't make sense to keep an out-of-date \code{A} slot
  and its factors in an updated object
}
\examples{
pred <- createRSC(fl=Dyestuff$Batch, x=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
str(pred)
pred@theta[] <- 0.782
with(Dyestuff, RSCupdate(pred, Yield))
str(pred)
pred@A
as(pred@A@factors[[1]], "sparseMatrix")
}

