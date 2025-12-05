## There are lots of problems with nlmer.
## If it eventually gets fixed, comparing the covariance structures between
## lme4::nlmer and the nlme::nlme may be worthwhile. 
## The following tests are skipped since (most) of them fail.

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

test_that("nlmer tests", {
  skip()
  ## Ensuring unstructured covariance results are the same as default
  expect_equal(coef(nm), coef(nm.us))
  expect_equal(fixef(nm), fixef(nm.us))
  ## Ensuring unstructured covariance results are the same as default
  expect_equal(fixef(nm), fixef(nm.us))
  
  # integration tests for sigma
  expect_equal(sigma(nm), sigma(nm.us))
  expect_equal_nocheck(sigma(nm), other_mod$nm.nlme_sigma, tolerance = 5e-4)
  expect_equal_nocheck(sigma(nm.cs), other_mod$nm.nlme.cs_sigma)

  ## Log likelihood tests 
  expect_equal(logLik(nm), logLik(nm.us))
  expect_equal_nocheck(logLik(nm.us), other_mod$nm.nlme_logLik)
  expect_equal_nocheck(logLik(nm.cs), other_mod$nm.nlme.cs_logLik)
  
  ## Integration tests for vcov
  nm_vcov <- suppressWarnings(vcov(nm))
  nm.us_vcov <- suppressWarnings(vcov(nm.us))
  expect_equal(nm_vcov, nm.us_vcov)
  expect_equal_nocheck(nm_vcov, other_mod$nm.nlme_vcov)
  nm.cs_vcov <- suppressWarnings(vcov(nm.cs))
  expect_equal_nocheck(nm.cs_vcov, other_mod$nm.nlme.cs_vcov)
  
  ## Integration tests for ranef
  expect_equal(ranef(nm), ranef(nm.us))
  expect_equal_nocheck(as.matrix(other_mod$nm.nlme_ranef), 
                         as.matrix(ranef(nm)$Tree))
  expect_equal_nocheck(as.matrix(other_mod$nm.nlme.cs_ranef), 
                         as.matrix(ranef(nm.cs)$Tree))
  
  ## Integration tests for predict
  expect_equal(predict(nm), predict(nm.us))
  expect_equal_nocheck(other_mod$nm.nlme_predict, predict(nm))
  expect_equal_nocheck(other_mod$nm.nlme.cs_predict, predict(nm.cs))
})