\name{simulate-mer}
\alias{simulate-mer}
\alias{simulate,mer-method}
\alias{refit,mer-method}
\title{Simulate based on mer fits}
\description{
  These generic functions (1) generate simulations based on the
  estimated fitted models (conditional on the estimated values of both
  the random and fixed effects) and (2) efficiently refit a
  specified model based on a new vector (or matrix) of responses
}
\usage{
\S4method{simulate}{mer}(object, nsim = 1 , seed = NULL, \dots)
\S4method{refit}{mer}(object, newresp, \dots)
}
\section{Methods}{
  \describe{
    \item{refit}{\code{signature(object = "mer", newresp = "numeric")}:
      Update the response vector only and refit the model. See
      \code{\link{refit}}.}
    \item{simulate}{\code{signature(object = "mer")}: simulate
      \code{nsim} (defaults to 1) responses from the theoretical
      distribution corresponding to the fitted model.  The \code{refit}
      method is particularly useful in combination with this method.
      See also \code{\link{simulate}}.}
  }
}

\arguments{
  \item{object}{An object of a suitable class - usually an
    \code{"\linkS4class{mer}"} object.}
  \item{nsim}{(integer) number of simulations to generate}
  \item{seed}{(integer or NULL) an object specifying if and how the
    random number generator
    should be initialized ('seeded').
    If an integer, \code{seed} is
    used in a call to \code{set.seed} before simulating the response
    vectors.  If set, the value is saved as the 'seed'
    attribute of the returned value.  The default, 'NULL', will
    not change the random generator state, and return
    '.Random.seed' as the 'seed' attribute}
  \item{newresp}{(numeric) a vector (for most models) or matrix (for
  binomial models where the response was specified as a matrix)
  containing new values for the response variable}
  \item{\dots}{Some methods may take additional, optional arguments.}
}
\value{
  Numeric vector, matrix, or list representing one or more
  simulations from the fitted model, conditional on both
  the fixed and random effect estimates.
  
  For binomial models where the response is specified as
  a two-column matrix, either a two-column numeric matrix
  with the same number of rows as the original data set
  (if \code{nsim==1}) or a list of such matrices (if \code{nsim>1}).

  For all other models, either a numeric vector matching the
  length of the original response vector (if \code{nsim==1}) or
  a matrix with \code{nsim} columns and the appropriate number
  of rows.
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

## generic parametric bootstrapping function; return a single simulated deviance
##  difference between full (m1) and reduced (m0) models under the
##  null hypothesis that the reduced model is the true model
pboot <- function(m0,m1) {
  s <- simulate(m0)
  L0 <- logLik(refit(m0,s))
  L1 <- logLik(refit(m1,s))
  2*(L1-L0)
}

obsdev <- c(2*(logLik(fm1)-logLik(fm2)))
\dontrun{
## parametric bootstrap test of significance of correlation between
##   random effects of `(Intercept)` and Days
## Timing: approx. 70 secs on a 2.66 GHz Intel Core Duo laptop
set.seed(1001)
sleepstudy_PB <- replicate(500,pboot(fm2,fm1))
}
library(lattice)
## null value for correlation parameter is *not* on the boundary
## of its feasible space, so we would expect chisq(1)
qqmath(sleepstudy_PB,distribution=function(p) qchisq(p,df=1),
     type="l",
            prepanel = prepanel.qqmathline,
            panel = function(x, ...) {
               panel.qqmathline(x, ...)
               panel.qqmath(x, ...)
            })
## classical test
pchisq(obsdev,df=1,lower.tail=FALSE)
## parametric bootstrap-based test
mean(sleepstudy_PB>obsdev)
## pretty close in this case

## cbpp data
## PB test of significance of main effect of period
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              family = binomial, data = cbpp)
gm0 <- update(gm1, . ~. -period)

\dontrun{
cbpp_PB <- replicate(500,pboot(gm0,gm1))
}

obsdev <- c(2*(logLik(gm1)-logLik(gm0)))
qqmath(cbpp_PB,distribution=function(p) qchisq(p,df=3),
     type="l",
            prepanel = prepanel.qqmathline,
            panel = function(x, ...) {
               panel.qqmathline(x, ...)
               panel.qqmath(x, ...)
            })
## classical test
pchisq(obsdev,df=3,lower.tail=FALSE)
## parametric bootstrap-based test
mean(cbpp_PB>obsdev)

## alternative plot
nsim <- length(cbpp_PB)
pdat <- data.frame(PB=(1:nsim)/(nsim+1),
       nominal=pchisq(sort(cbpp_PB,decreasing=TRUE),3,lower.tail=FALSE))
xyplot(nominal~PB,data=pdat,type="l",grid=TRUE,
       scales=list(x=list(log=TRUE),y=list(log=TRUE)),
       panel=function(...) {
          panel.abline(a=0,b=1,col="darkgray")
          panel.xyplot(...)
       })
}
\keyword{methods}
\keyword{datagen}
