\name{cbpp_PB}
\docType{data}
\alias{cbpp_PB}
\title{
  Stored parametric bootstrap samples for \code{\link{cbpp}} data
}
\description{
  Examples of parametric bootstrap distributions
  computed from models fitted to the
  \code{\link{cbpp}} data set
}
\usage{cbpp_PB}
\format{
  Contains a vector of 500 computed deviances
  under the (simulated) null hypothesis
}
\seealso{\code{\link{simulate-mer}} for examples of
  what to do with the parametric boostrap results
}
\examples{
## PB test of significance of main effect of period
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              family = binomial, data = cbpp)
gm0 <- update(gm1, . ~. -period)

## generic parametric bootstrapping function; return a single simulated deviance
##  difference between full (`m1') and reduced (`m0') models under the
##  null hypothesis that the reduced model is the true model
pboot <- function(m0,m1) {
  s <- simulate(m0)
  L0 <- logLik(refit(m0,s))
  L1 <- logLik(refit(m1,s))
  2*(L1-L0)
}
obsdev <- c(2*(logLik(gm1)-logLik(gm0)))
\dontrun{
## parametric bootstrap test of significance of correlation between
##   random effects of `(Intercept)` and Days
## Timing approx. 240 secs on a 2.66 GHz Intel Core Duo laptop
set.seed(1001)
cbpp_PB <- replicate(500,pboot(gm0,gm1))
}
}
\keyword{datasets}
