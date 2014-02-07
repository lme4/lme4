\name{simulate.merMod}
\alias{simulate.merMod}
\alias{simulate.formula}
\alias{.simulateFun}
\title{Simulate responses from a \code{\linkS4class{merMod}} object}
\usage{
  \method{simulate}{merMod} (object, nsim = 1, seed = NULL,
    use.u = FALSE, re.form=NA, ReForm, REForm, REform,
    newdata=NULL, newparams=NULL,
    family=NULL,
    allow.new.levels=FALSE, na.action=na.pass, \dots)
  \method{simulate}{formula} (object, nsim = 1 , seed = NULL, family, weights=NULL, offset=NULL, \dots)
.simulateFun(object, nsim = 1, seed = NULL, use.u = FALSE,
                            re.form=NA, ReForm, REForm, REform,
                            newdata=NULL, newparams=NULL,
                            formula=NULL,family=NULL,
                            weights=NULL, offset=NULL,
                            allow.new.levels=FALSE, na.action=na.pass, \dots)
}
\arguments{
  \item{object}{(for \code{simulate.merMod}) a fitted model object or
    (for \code{simulate.formula}) a (one-sided) mixed model
    formula, as described for \code{\link{lmer}}}
  \item{formula}{a (one-sided) mixed model
    formula, as described for \code{\link{lmer}}}
  \item{nsim}{positive integer scalar - the number of
    responses to simulate}
  \item{seed}{an optional seed to be used in
    \code{set.seed} immediately before the simulation so as
    to generate a reproducible sample.}
  \item{use.u}{(logical) if \code{TRUE}, generate a
    simulation conditional on the current random-effects
    estimates; if \code{FALSE} generate new Normally
    distributed random-effects values. (Redundant with \code{re.form},
    which is preferred: \code{TRUE} corresponds to \code{re.form=NULL}
    (condition on all random effects), while \code{FALSE} corresponds to
    \code{re.form=~0} (condition on none of the random effects).)}
  \item{re.form}{formula for random effects to condition on.  If \code{NULL},
    include all random effects; if \code{NA} or \code{~0},
    include no random effects.}
  \item{ReForm}{allowed for backward compatibility: \code{re.form} is
      now the preferred argument name}
  \item{REForm}{allowed for backward compatibility: \code{re.form} is
    now the preferred argument name}
  \item{REform}{allowed for backward compatibility: \code{re.form} is
      now the preferred argument name}
  \item{newdata}{data frame for which to evaluate
    predictions}
  \item{newparams}{new parameters to use in evaluating predictions,
    specified as in the \code{start} parameter for \code{\link{lmer}}
    or \code{\link{glmer}} -- a list with components \code{theta} and
    \code{beta} and (for LMMs or GLMMs that estimate a scale parameter)
    \code{sigma}}
  \item{family}{a GLM family, as in \code{\link{glmer}}}
  \item{weights}{prior weights, as in \code{\link{lmer}}
    or \code{\link{glmer}}}
  \item{offset}{offset, as in \code{\link{glmer}}}
  \item{allow.new.levels}{(logical) if FALSE (default),
    then any new levels (or NA values) detected in
    \code{newdata} will trigger an error; if TRUE, then the
    prediction will use the unconditional (population-level)
    values for data with previously unobserved levels (or
    NAs)}
  \item{na.action}{what to do with \code{NA} values in new data:
    see \code{\link{na.fail}}}
  \item{\dots}{optional additional arguments: none are used at present}
}
\seealso{
  \code{\link{bootMer}} for \dQuote{simulestimate}, i.e., where each
  simulation is followed by refitting the model.
}
\description{
  Simulate responses from the model represented by a \code{"merMod"} fitted
  model object.
}
\details{
  \itemize{
    \item{ordinarily \code{simulate} is used to generate new
      values from an existing, fitted model (\code{merMod} object):
      however, if \code{formula}, \code{newdata}, and \code{newparams} are
      specified, \code{simulate} generates the appropriate model
      structure to simulate from}
  }
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
form <- formula(gm1)[-2]
simulate(form,newdata=simdat,family=binomial,
    newparams=params)
}
