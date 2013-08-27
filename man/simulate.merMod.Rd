\name{simulate.merMod}
\alias{simulate.merMod}
\title{Simulate responses from a \code{\linkS4class{merMod}} object}
\usage{
  \method{simulate}{merMod} (object, nsim = 1, seed = NULL,
    use.u = FALSE, ReForm=NA, newdata=NULL, newparams=NULL,
                            allow.new.levels=FALSE, \dots)
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
  distributed random-effects values. (Redundant with \code{ReForm},
  which is preferred: \code{TRUE} corresponds to \code{ReForm=NULL}
  (condition on all random effects), while \code{TRUE} corresponds to
  \code{ReForm=~0} (condition on none of the random effects).)}
\item{ReForm}{formula for random effects to condition on.  If \code{NULL},
  include all random effects; if \code{NA} or \code{~0},
  include no random effects.}
  \item{newdata}{data frame for which to evaluate
    predictions}
  \item{newparams}{new parameters to use in evaluating predictions,
    specified as in the \code{start} parameter for \code{\link{lmer}}
    or \code{\link{glmer}} -- a list with components \code{theta} and/or
    (for GLMMs) \code{beta}}
    \item{allow.new.levels}{(logical) if FALSE (default),
    then any new levels (or NA values) detected in
    \code{newdata} will trigger an error; if TRUE, then the
    prediction will use the unconditional (population-level)
    values for data with previously unobserved levels (or
    NAs)}
  \item{\dots}{optional additional arguments: none are used at present}
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
##
## simulate from a non-fitted model; in this case we are just
## replicating the previous model, but
params <- list(theta=0.5,beta=c(2,-1,-2,-3))
simdat <- with(cbpp,expand.grid(herd=levels(herd),period=factor(1:4)))
simdat$size <- 15
simdat$incidence <- sample(0:1,size=nrow(simdat),replace=TRUE)
lme4:::simulate.merMod(formula=formula(gm1),newdata=simdat,family=binomial,
    newparams=params)
}

