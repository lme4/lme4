\name{createRSC}
\alias{createRSC}
\title{Create an RSC object from the rv and xv components}
\usage{
  createRSC(fl, xv, theta = rep.int(1, k),
    lower = numeric(k))
}
\arguments{
  \item{fl}{the matrix of row indices for the regular
  sparse column representation of Zt}

  \item{xv}{the non-zero values in ZtXt}

  \item{theta}{an optional variance-component parameter
  vector}

  \item{lower}{optional lower bounds on the
  variance-component parameter vector}
}
\value{
  an RSC object
}
\description{
  Create an RSC object from the rv and xv components
}
\examples{
pred <- createRSC(fl=Dyestuff$Batch, x=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
str(pred)
}

