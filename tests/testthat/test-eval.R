## examples for eval lookup
testthat::skip_on_cran()

if (require(car, quietly = TRUE)) {
  test_that("infIndexPlot env lookup OK", {
    fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
    ## silly test; the point is to see if this errors out with
    ## Error in as.list.environment(X[[i]], ...) :
    ## promise already under evaluation: recursive default argument reference or earlier problems?
    ##   Calls: infIndexPlot -> influence -> influence.merMod -> lapply -> FUN
    expect_equal(car::infIndexPlot(influence(fm1, "Subject")), NULL)
  })
}

if (require(rr2, quietly = TRUE)) {
  test_that("rr2 env lookup OK", {
    ## Error under alternate eval lookup
    ##     Error: bad 'data': object 'd' not found

    set.seed(123456)
    p1 <- 10; nsample <- 20; n <- p1 * nsample
    d <- data.frame(x1 = rnorm(n = n),
                    x2 = rnorm(n = n),
                    u1 = rep(1:p1, each = nsample),
                    u2 = rep(1:p1, times = nsample))
    d$u1 <- as.factor(d$u1); d$u2 <- as.factor(d$u2)

    ## LMM: y with random intercept
    b1 <- 1; b2 <- -1; sd1 <- 1.5
    d$y_re_intercept <- b1 * d$x1 + b2 * d$x2 +
      rep(rnorm(n = p1, sd = sd1), each = nsample) +  # random intercept u1
      rep(rnorm(n = p1, sd = sd1), times = nsample) + # random intercept u2
      rnorm(n = n)

    z.f2 <- lme4::lmer(y_re_intercept ~ x1 + x2 + (1 | u1) + (1 | u2), data = d, REML = T)
    ## NOTE, fails to produce warnings on second run of devtools::test()
    ## (possible interference from lmerTest methods being loaded ...?)
    expect_warning(R2(z.f2), "mod updated with REML = F")
})
}

## semEff::VIF() is here being applied to a previously fitted lmer
## model (shipley.growth[[3]])
## previous messing around with env evaluation had messed this up
if (suppressWarnings(require(semEff))) {
  ## suppress warning about 'cov2cor' import replacement
  test_that("semEff env lookup OK", {
    ##   Error in as.list.environment(X[[i]], ...) :
    ##     promise already under evaluation: recursive default argument reference or earlier problems?
    ##   Calls: VIF ... update.merMod -> do.call -> lapply -> FUN -> as.list.environment
  m <- shipley.growth[[3]]
  expect_equal(VIF(m),
               c(Date = 6.06283840168881,
                 DD = 6.07741017455859, lat = 1.01215136160858)
               )
  })
}
