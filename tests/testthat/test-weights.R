# unit tests for weights 

test_that("weights for lmer is always 1", {
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  ## the default should be just the prior weights,
  ## tried both anyways in case it changes.
  wgt_default <- weights(fm1)
  wgt_prior <- weights(fm1, "prior")
  wgt_working <- weights(fm1, "working")
  
  expect_equal(length(wgt_default), nrow(sleepstudy))
  expect_equal(length(wgt_prior), nrow(sleepstudy))
  expect_equal(length(wgt_working), nrow(sleepstudy))
  ## checking it is indeed always a vector of 1
  expect_equal(wgt_default, rep(1, nrow(sleepstudy)))
  expect_equal(wgt_prior, rep(1, nrow(sleepstudy)))
  expect_equal(wgt_working, rep(1, nrow(sleepstudy)))
})

test_that("testing na.action", {
  sstud <- sleepstudy
  sstud$Reaction[1] <- NA
  fm1 <- lmer(Reaction ~ Days + (Days | Subject), sstud)
  fm1.exclude <- lmer(Reaction ~ Days + (Days | Subject), sstud, 
              na.action = na.exclude)
  fm1.omit <- lmer(Reaction ~ Days + (Days | Subject), sstud, 
                      na.action = na.omit)

  expect_equal(length(weights(fm1)), nrow(sleepstudy) - 1)
  expect_equal(length(weights(fm1.exclude)), nrow(sleepstudy))
  expect_equal(length(weights(fm1.omit)), nrow(sleepstudy) - 1)
  ## this is not necessarily for weights, but to ensure na.exclude and na.omit
  ## are acting as expected/
  expect_equal(length(residuals(fm1.exclude)), nrow(sleepstudy))
  expect_equal(length(residuals(fm1.omit)), nrow(sleepstudy) - 1)
  
  cbpp3 <- cbpp
  cbpp3$incidence[1] <- NA
  cbpp3$size[1] <- NA
  
  gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                    data = cbpp3, family = binomial)
  gm1.exclude <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                       data = cbpp3, family = binomial, na.action = na.exclude)
  gm1.omit <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                    data = cbpp3, family = binomial, na.action = na.omit)
  expect_equal(length(weights(gm1)), nrow(cbpp) - 1)
  expect_equal(length(weights(gm1.exclude)), nrow(cbpp))
  expect_equal(length(weights(gm1.omit)), nrow(cbpp) - 1)
})
