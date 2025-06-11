library(lme4) 
library(testthat)

# UNIT TESTS FOR DiagonalCov

context("Testing DiagonalCov Class and Methods")

test_that("DiagonalCov object creation and initilialization works", {
	  # Test creation with a specified dimension
	  dim <- 3
	  diag_cov_obj <- new("DiagonalCov", dimension = dim)

	  expect_S4_class(diag_cov_obj, "DiagonalCov")
	  expect_equal(diag_cov_obj@dimension, dim)

	  # Test that start values are set correctly
	  start_vals <- get_start_values(diag_cov_obj)
	  expect_equal(length(start_vals), dim)
	  expect_equal(start_vals, rep(log(1), dim))
})




