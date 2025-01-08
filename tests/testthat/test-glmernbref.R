library(testthat)
## DON'T load lme4; test is to see if glmer.nb works when
## lme4 is not loaded
## this does *not* work properly in a devtools::test environment
##  (lme4 is not really detached)
## see tests/test-glmernbref.R for working test ...

test_that("glmer.nb ref to glmer", {
  set.seed(101)
  dd <- data.frame(x=runif(200), f= rep(1:20, each=10))
  b <- rnorm(20)
  dd <- transform(dd, y = rnbinom(200, mu  = exp(1 + 2*x + b[f]), size = 2))
  ## lme4 may not be attached if running tests via devtools::test()
  if ("package:lme4" %in% search()) detach("package:lme4")
  g <- lme4::glmer.nb(y~x + (1|f), data = dd)
  expect_is(g, "glmerMod")
})
