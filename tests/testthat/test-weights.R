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

test_that("a zero prior weight in lmer is the same as dropping the row (GH #880)", {
  set.seed(1)
  d <- sleepstudy
  d$w <- runif(nrow(d), 0.5, 2)
  drop <- c(3, 50, 120)
  d$w0 <- replace(d$w, drop, 0)
  for (reml in c(TRUE, FALSE)) {
    mz <- lmer(Reaction ~ Days + (Days | Subject), d, weights = w0, REML = reml)
    md <- lmer(Reaction ~ Days + (Days | Subject), d[-drop, ], weights = w, REML = reml)
    expect_equal(fixef(mz), fixef(md), tolerance = 1e-6)
    expect_equal(sigma(mz), sigma(md), tolerance = 1e-6)
    expect_equal(getME(mz, "theta"), getME(md, "theta"), tolerance = 1e-5)
    expect_equal(c(logLik(mz)), c(logLik(md)), tolerance = 1e-8)
    expect_equal(deviance(mz, REML = !reml), deviance(md, REML = !reml),
                 tolerance = 1e-8)
    expect_equal(nobs(mz), nobs(md))
    expect_equal(df.residual(mz), df.residual(md))
    expect_equal(BIC(mz), BIC(md), tolerance = 1e-8)
  }
  mz <- lmer(Reaction ~ Days + (Days | Subject), d, weights = w0)
  md <- lmer(Reaction ~ Days + (Days | Subject), d[-drop, ], weights = w)
  expect_equal(c(logLik(refitML(mz))), c(logLik(refitML(md))), tolerance = 1e-8)
  expect_equal(sigma(refitML(mz)), sigma(refitML(md)), tolerance = 1e-6)
})

test_that("deviance on the other REML scale is finite under na.exclude", {
  d <- sleepstudy
  d$Reaction[5] <- NA
  m_excl <- lmer(Reaction ~ Days + (Days | Subject), d, na.action = na.exclude)
  m_omit <- lmer(Reaction ~ Days + (Days | Subject), d, na.action = na.omit)
  expect_equal(deviance(m_excl, REML = FALSE), deviance(m_omit, REML = FALSE))
})

test_that("nobs() of glmer leaves out zero prior weights, as nobs.glm() does", {
  cb <- cbpp
  cb$p <- cb$incidence / cb$size
  cb$size0 <- replace(cb$size, c(2, 10), 0)
  gz <- glmer(p ~ period + (1 | herd), cb, binomial, weights = size0)
  gd <- glmer(p ~ period + (1 | herd), cb[-c(2, 10), ], binomial, weights = size)
  expect_equal(nobs(gz), nobs(gd))
  expect_equal(BIC(gz), BIC(gd), tolerance = 1e-6)
})
