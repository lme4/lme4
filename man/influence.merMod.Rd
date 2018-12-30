\name{influence.merMod}
\alias{influence.merMod}
\alias{dfbeta.influence.merMod}
\alias{dfbetas.influence.merMod}
\alias{cooks.distance.influence.merMod}
\alias{cooks.distance.merMod}
\title{
  Influence Diagnostics for Mixed-Effects Models
}
\description{
  These functions compute deletion influence diagnostics for linear
  (fit by \code{\link{lmer}} and generalized linear mixed-effects models
  (fit by \code{\link[lme4]{glmer}}). The main functions are methods for the \code{\link{influence}} generic function. Other functions are provided for
computing \code{\link{dfbeta}}, \code{\link{dfbetas}}, \code{\link{cooks.distance}}, and influence on variance-covariance components based
on the objects computed by \code{influence.merMod}
}


\usage{
\method{influence}{merMod}(model, groups, data, maxfun=1000, do.coef =
TRUE, ...)
\method{cooks.distance}{influence.merMod}(model, ...)
\method{dfbeta}{influence.merMod}(model, which = c("fixed", "var.cov"), ...)
\method{dfbetas}{influence.merMod}(model, ...)
}

\arguments{
  \item{model}{in the case of \code{influence.merMod}, a model of class \code{"merMod"};
  in the case of \code{cooks.distance}, \code{dfbeta}, or \code{dfbetas}, an object returned by
  \code{influence.merMod}}
  \item{groups}{a character vector containing the name of a grouping factor or names of grouping factors; if more than one name
  is supplied, then groups are defined by all combinations of levels of the grouping factors that appear in
  the data. If omitted, then each individual row of the data matrix is treated as a "group" to be deleted in turn.}
  \item{data}{an optional data frame with the data to which \code{model}
  was fit; \code{influence.merMod} can usually retrieve the data used to
  fit the model, unless it can't be found in the current environment, so it's usually unnecessary to supply this argument.}
  \item{maxfun}{The maximum number of function evaluations (for \code{influence.merMod})
  to perform after deleting each group; the defaults are large enough so that the iterations will typically continue to convergence.
  Setting to \code{maxfun=20} for an \code{lmer} model or \code{100} for a \code{glmer} model will typically produce a faster reasonable approximation. 
  An even smaller value can be used if interest is only in influence on the fixed effects.}
  \item{which}{if \code{"fixed.effects"} (the default), return influence
  on the fixed effects; if \code{"var.cov"}, return influence on the
  variance-covariance components.}
\item{do.coef}{if \code{FALSE}, skip potentially time-consuming
  computations; return just a list containing hat values}
  \item{\dots}{ignored.}
}

\details{
\code{influence.merMod} start with the estimated variance-covariance components from \code{model} and then refit
the model omitting each group in turn, not necessarily iterating to completion. For example, \code{maxfun=20} takes up to 20 function evaluations
step away from the ML or REML solution for the full data, which usually provides decent approximations to the fully iterated estimates.

The other functions are methods for the \code{\link{dfbeta}}, \code{\link{dfbetas}}, and \code{\link{cooks.distance}} generics, to be applied to the
\code{"influence.merMod"} object produced by the \code{influence} function; the \code{dfbeta} methods can also return
influence on the variance-covariance components.
}

\value{\code{influence.merMod} returns objects of class
  \code{"influence.merMod"}, which contain the following elements:
\describe{
    \item{\code{"fixed.effects"}}{the estimated fixed effects for the model.}
    \item{\code{"fixed.effects[-groups]"}}{a matrix with columns corresponding to the fixed-effects coefficients and rows corresponding to groups, giving the
    estimated fixed effects with each group deleted in turn; \emph{groups} is formed from the name(s) of the grouping factor(s).}
    \item{\code{"var.cov.comps"}}{the estimated variance-covariance parameters for the model.}
    \item{\code{"var.cov.comps[-groups]"}}{a matrix with the estimated covariance parameters (in columns) with each group deleted in turn.}
    \item{\code{"vcov"}}{The estimated covariance matrix of the fixed-effects coefficients.}
    \item{\code{"vcov[-groups]"}}{a list each of whose elements is the estimated covariance matrix of the fixed-effects coefficients with one group deleted.}
    \item{\code{"groups"}}{a character vector giving the names of the grouping factors.}
    \item{\code{"deleted"}}{the possibly composite grouping factor, each of whose elements is deleted in turn.}
    \item{\code{"converged"}}{for \code{influence.merMod}, a logical vector indicating whether the computation converged for each group.}
    \item{\code{"function.evals"}}{for \code{influence.merMod}, a vector of the number of function evaluations performed for each group.}
}
For plotting \code{"influence.merMod"} objects, see \code{\link[car]{infIndexPlot}}.
}
\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{
J. Fox \email{jfox@mcmaster.ca}
}

\seealso{
\code{\link[car]{infIndexPlot}}, \code{\link{influence.measures}}
}

\examples{
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
inf_fm1 <- influence(fm1, "Subject")
if (require("car")) {
    infIndexPlot(inf_fm1)
}
dfbeta(inf_fm1)
dfbetas(inf_fm1)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
inf_gm1 <- influence(gm1, "herd", maxfun=100)
gm1.11 <- update(gm1, subset = herd != 11) # check deleting herd 11
if (require("car")) {
    infIndexPlot(inf_gm1)
    compareCoefs(gm1, gm1.11)
}
## FIXME: dfbetas will fail at present *if* car is available, because
##  influence.merMod and dfbetas.influence.merMod methods are being
##  gotten from car, which has a sqrt(diag(x)) applied to a dpoMatrix ...
## dfbeta(inf_gm1)
## dfbetas(inf_gm1)
}
\keyword{models}
