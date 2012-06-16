\name{refit}
\alias{refit}
\alias{refit.merMod}
\title{Refit a model by maximum likelihood criterion}
\usage{
  refit(object, newresp, ...)

  \method{refit}{merMod} (object, newresp = NULL, ...)
}
\arguments{
  \item{object}{a fitted model, usually of class
  \code{"\linkS4class{lmerMod}"}, to be refit with a new
  response}

  \item{newresp}{a numeric vector providing the new
  response. Must be of the same length as the original
  response.}

  \item{...}{optional additional parameters.  None are used
  at present.}
}
\value{
  an object like \code{x} but fit by maximum likelihood
}
\description{
  Refit a model with a different response vector
}
\details{
  Refit a model after modifying the response vector.  This
  could be done using an \code{\link{update}} method but
  this approach should be faster because it bypasses the
  creation of the model representation and goes directly to
  the optimization step.
}
\examples{
## using refit() to fit each column in a matrix of responses
set.seed(101)
Y <- matrix(rnorm(1000),ncol=10)
res <- list()
d <- data.frame(y=Y[,1],x=rnorm(100),f=rep(1:10,10))
fit1 <- lmer(y~x+(1|f),data=d)
res <- c(fit1,lapply(as.data.frame(Y[,-1]),
        refit,object=fit1))
}

