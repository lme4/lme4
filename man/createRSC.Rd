\name{createRSC}
\alias{createRSC}
\title{Create an RSC object from the formula and the model.frame}
\usage{
  createRSC(bb, X, fr)
}
\arguments{
  \item{bb}{list of random-effects terms (ie. result of
  findbars)}

  \item{X}{the fixed-effects model matrix}

  \item{fr}{the model frame for the mixed-effects model}
}
\value{
  an RSC object
}
\description{
  Create an RSC object from the formula and the model.frame
}
\examples{
pred <- createRSC(fl=Dyestuff$Batch, x=matrix(1, nrow=2L, ncol=nrow(Dyestuff)))
str(pred)
}

