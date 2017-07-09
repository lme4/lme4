\name{simulate.merMod}
\title{Simulate Responses From \code{\linkS4class{merMod}} Object}
\alias{simulate.merMod}
\alias{simulate.formula}
\alias{.simulateFun}
\description{
  Simulate responses from a \code{"merMod"} fitted model object, i.e.,
  from the model represented by it.
}
\usage{
\method{simulate}{merMod}(object, nsim = 1, seed = NULL,
	 use.u = FALSE, re.form = NA, ReForm, REForm, REform,
	 newdata=NULL, newparams=NULL, family=NULL,
	 allow.new.levels = FALSE, na.action = na.pass, \dots)

\method{simulate}{formula}(object, nsim = 1 , seed = NULL,
         family, weights=NULL, offset=NULL, \dots)

.simulateFun(object, nsim = 1, seed = NULL, use.u = FALSE,
             re.form = NA, ReForm, REForm, REform,
             newdata=NULL, newparams=NULL,
             formula=NULL, family=NULL, weights=NULL, offset=NULL,
             allow.new.levels = FALSE, na.action = na.pass,
             cond.sim = TRUE, \dots)
}
\arguments{
  \item{object}{(for \code{simulate.merMod}) a fitted model object or
    (for \code{simulate.formula}) a (one-sided) mixed model formula, as
    described for \code{\link{lmer}}.}
  \item{nsim}{positive integer scalar - the number of responses to simulate.}
  \item{seed}{an optional seed to be used in \code{\link{set.seed}}
    immediately before the simulation so as to generate a reproducible sample.}
  \item{use.u}{(logical) if \code{TRUE}, generate a simulation
    conditional on the current random-effects estimates; if \code{FALSE}
    generate new Normally distributed random-effects values. (Redundant
    with \code{re.form}, which is preferred: \code{TRUE} corresponds to
    \code{re.form = NULL} (condition on all random effects), while
    \code{FALSE} corresponds to \code{re.form = ~0} (condition on none
    of the random effects).)}
  \item{re.form}{formula for random effects to condition on.  If
    \code{NULL}, condition on all random effects; if \code{NA} or \code{~0},
    condition on no random effects.  See Details.}
  \item{ReForm, REForm, REform}{deprecated: \code{re.form} is
    now the preferred argument name.}
  \item{newdata}{data frame for which to evaluate predictions.}
  \item{newparams}{new parameters to use in evaluating predictions,
    specified as in the \code{start} parameter for \code{\link{lmer}} or
    \code{\link{glmer}} -- a list with components \code{theta} and
    \code{beta} and (for LMMs or GLMMs that estimate a scale parameter)
    \code{sigma}}
  \item{formula}{a (one-sided) mixed model formula, as described for
    \code{\link{lmer}}.}
  \item{family}{a GLM family, as in \code{\link{glmer}}.}
  \item{weights}{prior weights, as in \code{\link{lmer}} or
    \code{\link{glmer}}.}
  \item{offset}{offset, as in \code{\link{glmer}}.}
  \item{allow.new.levels}{(logical) if FALSE (default), then any new
    levels (or \code{NA} values) detected in \code{newdata} will trigger an
    error; if TRUE, then the prediction will use the unconditional
    (population-level) values for data with previously unobserved levels
    (or \code{NA}s).}
  \item{na.action}{what to do with \code{NA} values in new data: see
    \code{\link{na.fail}}}
  \item{cond.sim}{(experimental) simulate the conditional
    distribution?  if \code{FALSE}, simulate only random effects; do not
    simulate from the conditional distribution, rather return the
    predicted group-level values}
  \item{\dots}{optional additional arguments: none are used at present.}
}
\seealso{
  \code{\link{bootMer}} for \dQuote{simulestimate}, i.e., where each
  simulation is followed by refitting the model.
}
\details{
  \itemize{
    \item{ordinarily \code{simulate} is used to generate new
      values from an existing, fitted model (\code{merMod} object):
      however, if \code{formula}, \code{newdata}, and \code{newparams} are
      specified, \code{simulate} generates the appropriate model
      structure to simulate from.}

    \item{The \code{re.form} argument allows the user to specify how the
      random effects are incorporated in the simulation.  All of the
      random effects terms included in \code{re.form} will be
      \emph{conditioned on} - that is, the conditional modes of those
      random effects will be included in the deterministic part of the
      simulation. (If new levels are used (and \code{allow.new.levels}
      is \code{TRUE}), the conditional modes for these levels will be
      set to the population mode, i.e. values of zero will be used for
      the random effects.)  Conversely, the random effect terms that are
      \emph{not} included in \code{re.form} will be \emph{simulated
      from} - that is, new values will be chosen for each group based on
      the estimated random-effects variances.

      The default behaviour (using \code{re.form=NA}) is to condition on
      none of the random effects, simulating new values for all of the
      random effects.

    }
    \item{For Gaussian fits, \code{sigma} specifies the residual
      standard deviation; for Gamma fits, it specifies the shape
      parameter (the rate parameter for each observation i
      is calculated as shape/mean(i)).  For negative binomial fits,
      the overdispersion parameter is specified via the family,
      e.g. \code{simulate(..., family=negative.binomial(theta=1.5))}.
    }
    \item{For binomial models, \code{simulate.formula} looks for the
      binomial size first in the \code{weights} argument (if it's supplied),
      second from the left-hand side of the formula (if the formula has been
      specified in success/failure form), and defaults to 1 if neither of
      those have been supplied.
      Simulated responses will be given as proportions, unless the supplied
      formula has a matrix-valued left-hand side, in which case they will be
      given in matrix form. If a left-hand side is given, variables in that
      expression must be available in \code{newdata}.
    }
    \item{For negative binomial models, use the \code{negative.binomial}
      family (from the \code{MASS} package)
      and specify the overdispersion parameter via the
      \code{theta} (sic) parameter of the family function, e.g.
      \code{simulate(...,family=negative.binomial(theta=1))} to simulate
      from a geometric distribution (negative binomial with
      overdispersion parameter 1).
    }
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
## replicating the previous model, but starting from scratch
params <- list(theta=0.5,beta=c(2,-1,-2,-3))
simdat <- with(cbpp,expand.grid(herd=levels(herd),period=factor(1:4)))
simdat$size <- 15
simdat$incidence <- sample(0:1,size=nrow(simdat),replace=TRUE)
form <- formula(gm1)[-2]  ## RHS of equation only
simulate(form,newdata=simdat,family=binomial,
    newparams=params)
## simulate from negative binomial distribution instead
simulate(form,newdata=simdat,family=negative.binomial(theta=2.5),
    newparams=params)
}
