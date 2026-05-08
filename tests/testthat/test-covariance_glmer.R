## read system file 
glmmTMB_glmer_ar1 <- readRDS(
  system.file("testdata", "glmmTMB_glmer_ar1_comparison.rds", package = "lme4")
)
Contraception <- readRDS(
  system.file("testdata", "Contraception.rds", package = "lme4")
)

test_that("glmer + ar1 does not return zero values for fixed effects coefficients", {
  ## Test 1; just making sure the results match as expected
  form <- out_bin_1 ~ ftime + trt + ar1(ftime + 0 | id_cluster) +
    (1 | id_individual)
  dd <- glmmTMB_glmer_ar1$glmer_ar1_dd
  lme4_1 <- glmer(form, data = dd, family = binomial)
  
  thetas <- c(getReCovs(lme4_1)[[1]]@par, getReCovs(lme4_1)[[2]]@par)
  betas <- fixef(lme4_1)
  lme4_1@optinfo$val
  
  expect_true(all.equal(c(thetas, betas), lme4_1@optinfo$val, 
            check.attributes = FALSE, check.class = FALSE))
  
  ## Test 2; just making sure we still have nonzero results
  form2 <- update(form, out_gauss ~ .)
  lme4_2 <- lmer(form2, data = dd, REML = FALSE)
  expect_true(!0 %in% fixef(lme4_2))
  
  ## Test 3; just making sure the length of getME matches beta
  gm.ar1_0 <- glmer(use ~ age + urban + ar1(0 + urban | district),
                    data = Contraception,
                    family = binomial)
  expect_equal(length(colnames(getME(gm.ar1_0, "X"))), 
               length(getME(gm.ar1_0, "beta")))
})





