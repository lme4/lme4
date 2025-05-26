library("testthat")
library("lme4")

context("vcconv - testing safe_chol()")

test_that("safe_chol returns upper triangular matrix for a simple PD matrix", {
  M <- diag(c(3, 2, 1))
  L <- lme4:::safe_chol(M)
  expect_true(all(L[lower.tri(L)] == 0),
              info = "Matrix L should be upper triangular")
  expect_equal(crossprod(L), M, tolerance = 1e-8,
               info = "Reconstructed matrix should match the original PD matrix")
})

test_that("safe_chol works with a generic PD matrix", {
  M <- matrix(c(4, 1, 1,
                1, 3, 2,
                1, 2, 5), 3, 3)
  M <- (M + t(M))/2  ## ensure symmetry
  L <- lme4:::safe_chol(M)
  expect_true(all(L[lower.tri(L)] == 0),
              info = "Matrix L should be upper triangular")
  expect_equal(crossprod(L), M, tolerance = 1e-8,
               info = "Reconstructed matrix should equal M")
})

test_that("safe_chol handles a semi-definite (PSD) matrix", {
  V <- matrix(c(1600, 1160, 1760,
                1160, 3442, 2806,
                1760, 2806, 2836), 3, 3)
  V <- (V + t(V))/2  ## ensure symmetry
  L <- lme4:::safe_chol(V)
  expect_true(all(L[lower.tri(L)] == 0),
              info = "Matrix L should be upper triangular")
  expect_equal(crossprod(L), V, tolerance = 1e-8,
               info = "Reconstructed matrix should equal the PSD matrix")
})

test_that("safe_chol handles a rank-deficient matrix", {
  V <- matrix(c(1, 0, 0,
                0, 0, 0,
                0, 0, 0), 3, 3)
  L <- lme4:::safe_chol(V)
  expect_true(all(L[lower.tri(L)] == 0),
              info = "Matrix L should be upper triangular")
  expect_equal(crossprod(L), V, tolerance = 1e-8,
               info = "Reconstructed matrix should equal the rank-deficient matrix")
})

test_that("safe_chol handles a nearly singular matrix", {
  eps <- 1e-10
  M <- matrix(c(1, 1-eps,
                1-eps, 1), 2, 2)
  M <- (M + t(M))/2
  L <- lme4:::safe_chol(M)
  expect_true(all(L[lower.tri(L)] == 0),
              info = "Matrix L should be upper triangular")
  expect_equal(crossprod(L), M, tolerance = 1e-7,
               info = "Reconstructed matrix should match the nearly singular matrix within tolerance")
})

test_that("safe_chol throws an error for an indefinite matrix", {
  Bad <- matrix(c(-1, 0,
                  0, 2), 2, 2)
  expect_error(lme4:::safe_chol(Bad),
               regexp = "chol",
               info = "safe_chol() should fail for an indefinite matrix")
})
