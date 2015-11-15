library(lme4)
library(testthat)
load(system.file("testdata","crabs_randdata00.Rda",package="lme4"))

test_that('RZX is being calculated properly', {
  # this is a test for an old problem, documented here:
  # http://stevencarlislewalker.github.io/notebook/RZX_problems.html
  fr <- cbind(final.snail.density, snails.lost) ~ crab.speciesS + crab.sizeS + 
    crab.speciesS:crab.sizeS + (snail.size | plot)
  m <- glmer(fr, data = randdata00, family = binomial)
  expect_that(m, is_a("glmerMod"))
})
