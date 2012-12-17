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
  the 'ldL2' value for the updated A
}
\description{
  Update for the penalized least squares problem
}
\examples{
pred <- createRSC(rv=as.integer(Dyestuff$Batch),
                  xv=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
str(pred)
with(Dyestuff, RSCupdate(pred, Yield))
str(pred)
}

