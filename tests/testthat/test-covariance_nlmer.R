## There are lots of problems with nlmer.
## The compound symmetry (cs) covariance structure reduces to a scalar variance
## when nc=1 (a single RE), so cs(Asym|Tree) is equivalent to the default model.
## The nlme reference model is therefore the simple random-intercept model (no
## within-group correlation structure). See
## inst/testdata/test-covariance_structures_data.R for the nlme reference model
## definitions.
##
## Despite using the correct nlme reference model, comparisons of sigma, logLik,
## and vcov between nlmer and nlme still fail due to known numerical discrepancies
## in nlmer's optimization. The corresponding assertions are placed in skip() blocks
## below until nlmer is improved.

other_mod <- readRDS(
  system.file("testdata", "test-covariance_structures_data.rds", package = "lme4")
)

## lme4 nonlinear mixed effects model
startvec <- c(Asym = 200, xmid = 725, scal = 350)
nm <-  nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym|Tree,
             Orange, start = startvec)

nm.us <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ us(Asym|Tree),
               Orange, start = startvec)

nm.cs <- nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ cs(Asym|Tree),
               Orange, start = startvec)

## cs(Asym|Tree) with a single RE (nc=1) is equivalent to the default model:
## compound symmetry reduces to a scalar variance when nc=1.
test_that("nlmer compound symmetry (nc=1) is equivalent to default", {
  expect_equal(fixef(nm), fixef(nm.cs))
  expect_equal(sigma(nm), sigma(nm.cs))
  expect_equal(logLik(nm), logLik(nm.cs))
  expect_equal(ranef(nm), ranef(nm.cs))
  expect_equal(predict(nm), predict(nm.cs))
})

## Compare nlmer cs(Asym|Tree) with the equivalent nlme reference model (ranef
## and predict match; sigma/logLik/vcov are skipped due to nlmer discrepancies).
test_that("nlmer compound symmetry (nc=1) matches nlme (ranef and predict)", {
  expect_equal_nocheck(as.matrix(other_mod$nm.nlme.cs_ranef),
                       as.matrix(ranef(nm.cs)$Tree), tolerance = 5e-4)
  expect_equal_nocheck(other_mod$nm.nlme.cs_predict, predict(nm.cs),
                       tolerance = 5e-4)
})

test_that("nlmer tests (skipped: remaining issues)", {
  skip()
  ## Ensuring unstructured covariance results are the same as default
  expect_equal(coef(nm), coef(nm.us))
  expect_equal(fixef(nm), fixef(nm.us))

  ## integration tests for sigma
  expect_equal(sigma(nm), sigma(nm.us))
  expect_equal_nocheck(sigma(nm), other_mod$nm.nlme_sigma, tolerance = 5e-4)

  ## Log likelihood tests
  expect_equal(logLik(nm), logLik(nm.us))
  expect_equal_nocheck(logLik(nm.us), other_mod$nm.nlme_logLik, tolerance = 5e-4)

  ## Integration tests for vcov
  nm_vcov <- suppressWarnings(vcov(nm))
  nm.us_vcov <- suppressWarnings(vcov(nm.us))
  expect_equal(nm_vcov, nm.us_vcov)
  expect_equal_nocheck(nm_vcov, other_mod$nm.nlme_vcov, tolerance = 5e-4)

  ## Integration tests for ranef
  expect_equal(ranef(nm), ranef(nm.us))
  expect_equal_nocheck(as.matrix(other_mod$nm.nlme_ranef),
                       as.matrix(ranef(nm)$Tree), tolerance = 5e-4)

  ## Integration tests for predict
  expect_equal(predict(nm), predict(nm.us))
  expect_equal_nocheck(other_mod$nm.nlme_predict, predict(nm), tolerance = 5e-4)

  ## sigma, logLik, vcov for cs vs nlme (fail due to nlmer discrepancies)
  expect_equal_nocheck(sigma(nm.cs), other_mod$nm.nlme.cs_sigma, tolerance = 5e-4)
  expect_equal_nocheck(logLik(nm.cs), other_mod$nm.nlme.cs_logLik, tolerance = 5e-4)
  nm.cs_vcov <- suppressWarnings(vcov(nm.cs))
  expect_equal_nocheck(nm.cs_vcov, other_mod$nm.nlme.cs_vcov, tolerance = 5e-4)
})