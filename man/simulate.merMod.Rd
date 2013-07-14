\name{simulate.merMod}
\alias{simulate.merMod}
\title{Simulate responses from a \code{\linkS4class{merMod}} object}
\usage{
  \method{simulate}{merMod} (object, nsim = 1, seed = NULL,
    use.u = FALSE, ...)
}
\arguments{
  \item{object}{a fitted model object}

  \item{nsim}{positive integer scalar - the number of
  responses to simulate}

  \item{seed}{an optional seed to be used in
  \code{set.seed} immediately before the simulation so as
  to generate a reproducible sample.}

  \item{use.u}{(logical) if \code{TRUE}, generate a
  simulation conditional on the current random-effects
  estimates; if \code{FALSE} generate new Normally
  distributed random-effects values}

  \item{...}{optional additional arguments, none are used
  at present}
}
\description{
  Simulate responses from the model represented by a fitted
  model object
}
\examples{
## test whether fitted models are consistent with the
##  observed number of zeros in CBPP data set:
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
gg <- simulate(gm1,1000)
zeros <- sapply(gg,function(x) sum(x[,"incidence"]==0))
plot(table(zeros))
abline(v=sum(cbpp$incidence==0),col=2)
}

