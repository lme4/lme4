(testLevel <- if (nzchar(s <- Sys.getenv("LME4_TEST_LEVEL"))) as.numeric(s) else 1)

if (testLevel > 1) {

  ## read system file 
  glmmTMB_glmer_ar1 <- readRDS(
    system.file("testdata", "glmmTMB_glmer_ar1_comparison.rds", package = "lme4")
  )
  Contraception <- readRDS(
    system.file("testdata", "Contraception.rds", package = "lme4")
  )

  test_that("glmer + ar1 does not return zero values for fixed effects coefficients", {
    ## Test 1; make sure the results match as expected
    form <- out_bin_1 ~ ftime + trt + ar1(ftime + 0 | id_cluster) +
      (1 | id_individual)
    dd <- glmmTMB_glmer_ar1$glmer_ar1_dd
    lme4_1 <- glmer(form, data = dd, family = binomial)
    
    thetas <- c(getReCovs(lme4_1)[[1]]@par, getReCovs(lme4_1)[[2]]@par)
    betas <- fixef(lme4_1)
    expect_equal_nocheck(c(thetas, betas), lme4_1@optinfo$val)
    
    ## Test 2; make sure we still have nonzero results
    form2 <- update(form, out_gauss ~ .)
    lme4_2 <- lmer(form2, data = dd, REML = FALSE)
    expect_true(!all(fixef(lme4_2) == 0.0))
    
    ## Test 3; make sure ncol(X) matches beta
    gm.ar1_0 <- glmer(use ~ age + urban + ar1(0 + urban | district),
                      data = Contraception,
                      family = binomial)
    expect_equal(length(colnames(getME(gm.ar1_0, "X"))), 
                 length(getME(gm.ar1_0, "beta")))
  })
}

test_that("npar test works for models with 0 FE parameters", {
  expect_is(glmer(round(Reaction) ~ 0 + (1|Subject), family = poisson,
                  data = sleepstudy),
            "merMod")
})
