\name{sleepstudy_PB}
\docType{data}
\alias{sleepstudy_PB}
\title{
  Stored parametric bootstrap samples for the \code{sleepstudy} data
}
\description{
  Parametric bootstrap distribution
  computed from models fitted to the \code{\link{sleepstudy}}
  data set
}
\usage{sleepstudy_PB}
\format{
  Contains a vector of 500 computed deviances
  under the (simulated) null hypothesis
}
\seealso{\code{\link{simulate-mer}} for examples of
  what to do with the parametric boostrap results
}
\examples{
fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy)
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)

## generic parametric bootstrapping function; return a single simulated deviance
##  difference between full (`m1') and reduced (`m0') models under the
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
}
\keyword{datasets}
