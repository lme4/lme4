library(lme4)
library(testthat)
load(system.file("testdata","survdat_reduced.Rda",package="lme4"))

test_that('Step-halving works properly', {
  # this example is known to require step-halving (or at least has in the past
  # required step-halving)
  form <- survprop~(1|nobs)
  m <- glmer(form,weights=eggs,data=survdat_reduced,family=binomial,nAGQ=1L)
  expect_that(m, is_a("glmerMod"))
})
