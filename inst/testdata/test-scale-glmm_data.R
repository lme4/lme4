## Script to generate reference values from glmmTMB for scale-parameter GLMM tests.
## These values are used by tests/testthat/test-glmFamily.R to check that lme4's
## corrected sigma estimation (issue #936) agrees with glmmTMB.
##
## Run this script from the package root directory to regenerate the .rds:
##   Rscript inst/testdata/test-scale-glmm_data.R

library(lme4)
library(glmmTMB)
library(nlme)   # for Rail

## --- gaussian(link="log") with Rail data (original issue #936) ---
data("Rail", package = "nlme")
m_glmmTMB_gauss <- glmmTMB(travel ~ 1 + (1|Rail), data = Rail,
                             family = gaussian(link = "log"))

## --- Gamma with simulated data (same simfun_gam as test-glmFamily.R) ---
simfun_gam <- function(ngrp = 50, nrep = 50, shape_gam = 2, intercept = 1,
                       use_simulate = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dd <- expand.grid(group = 1:ngrp, rep = 1:nrep)
  if (use_simulate) {
    dd$y <- simulate(~ 1 + (1 | group), newdata = dd,
                     family = Gamma(link = "log"),
                     newparams = list(
                       theta = 1, beta = 1, sigma = 1/sqrt(shape_gam)))[[1]]
    return(dd)
  }
  b <- rnorm(ngrp)
  eta <- intercept + b
  mu <- exp(eta)
  y <- rgamma(nrow(dd), shape = shape_gam, scale = mu/shape_gam)
  data.frame(dd, y)
}
dd2 <- simfun_gam(seed = 101, use_simulate = TRUE)
m_glmmTMB_gam <- glmmTMB(y ~ 1 + (1|group), family = Gamma(link = "log"),
                          data = dd2)

scale_glmm_test <- list(
  ## gaussian(link="log") Rail reference values
  gauss_sigma  = sigma(m_glmmTMB_gauss),
  gauss_fixef  = unname(glmmTMB::fixef(m_glmmTMB_gauss)$cond),
  ## Gamma reference sigma
  gamma_sigma  = sigma(m_glmmTMB_gam)
)

saveRDS(scale_glmm_test,
        "inst/testdata/test-scale-glmm_data.rds")
