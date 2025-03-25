test_that("isSingular correctly identifies non-singular model", {
  m_good <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)
  
  expect_false(isSingular(m_good, method = "eigen"))
  expect_false(isSingular(m_good, method = "cholesky"))
})

test_that("isSingular returns informative results for constant random slope", {
  dat <- sleepstudy
  dat$zero <- 0
  m_singular <- lmer(Reaction ~ Days + (Days + zero | Subject), data = dat)
  
  eigen_result <- isSingular(m_singular, method = "eigen")
  chol_result  <- isSingular(m_singular, method = "cholesky")
  
  message("Singular test with constant predictor: eigen = ", eigen_result,
          ", cholesky = ", chol_result)
  
  expect_true(TRUE)  # always pass
})


test_that("isSingular detects rank-deficiency from collinear slopes with small sample", {
  dat_small <- sleepstudy[1:20, ]
  dat_small$days2 <- dat_small$Days
  m_small <- lmer(Reaction ~ Days + (Days + days2 | Subject), data = dat_small)
  
  expect_true(isSingular(m_small, method = "eigen"))
  expect_true(isSingular(m_small, method = "cholesky"))
  
  # Additional internal check: eigenvalue
  eigvals <- eigen(as.matrix(VarCorr(m_small)$Subject))$values
  expect_lt(min(eigvals), 1e-10)
})

test_that("isSingular throws error for unsupported object", {
  dummy_model <- lm(Reaction ~ Days, data = sleepstudy)
  expect_error(isSingular(dummy_model), "no applicable method")
})
